#' Process an .fcs files to remove debris and doublets
#'
#' \code{process_fcs} uses mixture models to cluster bacteria from background
#' debris and fits a linear model to SSC-H vs SSC-A to remove doublets.
#'
#' @param fcs_file path of an .fcs files to be processed.
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to "BL1-H".
#' @param do_plot a Boolean flag to determine whether to produce plots showing
#' the trimming of each \code{flowFrame} Defaults to \code{FALSE}.
#' @param pre_cleaned have you pre removed background debris.
#'
#' @return nothing is returned. A processed .fcs file is produced and a plot if
#' \code{do_plot = TRUE}.
#'
#' @export
process_fcs <- function(fcs_file, flu_channels=c("BL1-H"), do_plot = F,
                        pre_cleaned = F) {

  flow_frame <- flowCore::read.FCS(fcs_file, emptyValue = F)

  ## Trim 0 values and log10 transform
  prepped_flow_frame <- prep_flow_frame(flow_frame, flu_channels)

  ## Try to remove background debris by clustering
  bacteria_flow_frame <- get_bacteria(prepped_flow_frame, pre_cleaned)
  try(if (bacteria_flow_frame == 0) {
    stop()
  }, silent = T) # if we haven't found bacteria, stop

  ## Try to remove doublets
  singlet_flow_frame <- get_singlets(bacteria_flow_frame)

  ## Save processed flowFrames to a new folder
  out_name <- paste(substr(x = fcs_file, start = 1, stop = nchar(fcs_file) - 4),
                    "_processed.fcs", sep = "")

  flowCore::write.FCS(singlet_flow_frame, out_name)

  ##  Plot a grid of graphs showing the stages of trimming
  if (do_plot) {
    plot_trimming(flow_frame, bacteria_flow_frame, singlet_flow_frame,
                  NA, flu_channels, out_name, F)
  }
}

#' Process .fcs files within a directory to remove debris and doublets, and
#' optionally calibrate
#'
#' \code{process_fcs_dir} uses mixture models to cluster bacteria from
#' background debris and fits a linear model to SSC-H vs SSC-A to remove doublets.
#'
#' @param dir_path a directory path containing the .fcs files  to be parsed or
#' folders to be recursed through.
#' @param pattern a regex pattern to match particular .fcs files. Default is
#' \code{"*.fcs"} matching all .fcs files.
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to "BL1-H".
#' @param calibrate a Boolean flag to determine whether to convert fluorescence
#' to MEF values. Requires an .fcs file with named \code{"*beads*.fcs"}.
#' Defaults to \code{FALSE}.
#' @param do_plot a Boolean flag to determine whether to produce plots showing
#' the trimming of each \code{flowFrame}. Defaults to \code{FALSE}.
#' @param pre_cleaned have you pre removed background debris
#' @param neg_fcs filename of .fcs file for a negative control used for
#' fluorescence normalisation.
#' @param mef_peaks a list of lists in the form \code{list(list(channel="BL1-H",
#'  peaks=c(0, 200, ...)} of MEF fluorescence values for the calibration beads.
#'   Default values for BL1-H and YL2-H.
#' @param bead_dens_bw the bandwidth for the kernel density of the bead peak
#' data. Default = 0.025. Increase if erroneous peaks are being found, decrease
#'  if not enough peaks are found.
#' @param manual_peaks if bead peaks are not being found by the EM algorithm,
#' one can manually specify the positions of the bead peaks (on a log10 scale).
#'  A list of lists in the form
#'  \code{list(list(channel="BL1-H", peaks=c(2.1, 2.8, ...)} of log10
#'  fluorescence values for the calibration beads.
#'
#' @return nothing is returned. A new folder is created with the processed .fcs
#' files and plot if \code{do_plot = TRUE}.
#' @export
process_fcs_dir <- function(dir_path, pattern = "*.fcs", flu_channels=c("BL1-H"),
                            do_plot = F, pre_cleaned = F, neg_fcs = NA, calibrate = F,
                            mef_peaks = NA, bead_dens_bw = 0.025, manual_peaks = NA) {
  ## Create directory for processed flowFrames
  if (!dir.exists(paste(dir_path, "processed", sep = "_"))) {
    dir.create(paste(dir_path, "processed", sep = "_"), recursive = T)
  }

  if (calibrate) {
    ## First step is to get calibration standard curves
    bead_files <- list.files(path = dir_path,
                             pattern = utils::glob2rx("*beads*.fcs"),
                             full.names = T, recursive = T, include.dirs = T)

    if(length(bead_files) > 1){
      warning("More than one beads file found.")
      bead_file <- bead_files[1]

    } else if(length(bead_files) == 1){
      bead_file <- bead_files[1]
    } else {
      stop("No beads file found.")
    }

    calibration_parameters <- get_calibration(bead_file, flu_channels,
                                              mef_peaks, manual_peaks,
                                              bead_dens_bw)
  }

  if (!is.na(neg_fcs)) {
    print("Getting autofluorescence normalisation parameter")

    flow_frame <- flowCore::read.FCS(paste(dir_path, neg_fcs, sep = "/"), emptyValue = F)

    ## Trim 0 values and log10 transform
    prepped_flow_frame <- prep_flow_frame(flow_frame, flu_channels)

    ## Try to remove background debris by clustering
    bacteria_flow_frame <- get_bacteria(prepped_flow_frame, pre_cleaned)
    try(if (bacteria_flow_frame == 0) {
      next
    }, silent = T) # if we haven't found bacteria move on to the next flow frame

    ## Try to remove doublets
    negative_flow_frame <- get_singlets(bacteria_flow_frame)
  }

  all_files <- list.files(path = dir_path, pattern = utils::glob2rx(pattern),
                          full.names = T, recursive = T, include.dirs = T)
  print(paste("Trimming ", length(all_files), " .fcs files.", sep = ""))

  for (next_fcs in all_files) {
    flow_frame <- flowCore::read.FCS(next_fcs, emptyValue = F)

    ## Trim 0 values and log10 transform
    prepped_flow_frame <- prep_flow_frame(flow_frame, flu_channels)

    ## Try to remove background debris by clustering
    bacteria_flow_frame <- get_bacteria(prepped_flow_frame, pre_cleaned)
    try(if (bacteria_flow_frame == 0) {
      next
    }, silent = T) # if we haven't found bacteria move on to the next flow frame

    ## Try to remove doublets
    singlet_flow_frame <- get_singlets(bacteria_flow_frame)

    ## normalise autofluorescence
    normalised_flow_frame <- singlet_flow_frame
    if (!is.na(neg_fcs)) {
      for(flu in flu_channels) {

        ## get geometric mean of negative control for given fluorescence channel
        neg_mean <- exp(mean(log(negative_flow_frame@exprs[,flu]), trim = 0.05, na.rm = T))

        ## hack to create a new channel
        normalised_flow_frame <- flowCore::transform(normalised_flow_frame,
                                                     normFlu = 1)

        ## add normalised fluorescence to the channel
        normalised_flow_frame@exprs[,"normFlu"] <- normalised_flow_frame@exprs[,flu] - neg_mean

        ## rename the channel
        normalised_flu = paste("normalised_", flu, sep = "")
        old_names <- flowCore::colnames(normalised_flow_frame)
        old_names[length(old_names)] <- normalised_flu
        flowCore::colnames(normalised_flow_frame) <- old_names
      }
    }

    ## Convert to MEF
    calibrated_flow_frame <- normalised_flow_frame
    if (calibrate) {
      calibrated_flow_frame <- to_mef(calibrated_flow_frame, flu_channels,
                                      calibration_parameters, !is.na(neg_fcs))
    }

    out_flow_frame <- calibrated_flow_frame

    ## Save processed flowFrames to a new folder
    out_name <- paste(dirname(next_fcs), "_processed/", basename(next_fcs),
                      sep = "")
    flowCore::write.FCS(out_flow_frame, out_name)

    ##  Plot a grid of graphs showing the stages of trimming
    if (do_plot) {
      plot_trimming(flow_frame, bacteria_flow_frame, singlet_flow_frame,
                    normalised_flow_frame, out_flow_frame, flu_channels,
                    out_name, !is.na(neg_fcs), calibrate)
    }
  }
}


#' Get fluorescence calibration curve parameters
#'
#' @param bead_file path to beads .fcs file.
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to \code{"BL1-H"}.
#' @param mef_peaks a list of lists in the form \code{list(list(channel="BL1-H",
#'  peaks=c(0, 200, ...)} of MEF fluorescence values for the calibration beads.
#'   Default values for BL1-H and YL2-H.
#' @param manual_peaks if bead peaks are not being found by the EM algorithm,
#' one can manually specify the positions of the bead peaks (on a log10 scale).
#'  A list of lists in the form
#'  \code{list(list(channel="BL1-H", peaks=c(2.1, 2.8, ...)} of log10
#'  fluorescence values for the calibration beads.
#' @param bead_dens_bw the bandwidth for the kernel density of the bead peak
#' data. Default = 0.025. Increase if erroneous peaks are being found, decrease
#'  if not enough peaks are found.
#'
#' @return
#' @importFrom rlang .data
#'
get_calibration <- function(bead_file, flu_channels, mef_peaks, manual_peaks,
                            bead_dens_bw) {
  bead_frame <- flowCore::read.FCS(bead_file, emptyValue = F)

  out_name <- paste(dirname(bead_file),
                    "_processed/",
                    basename(unlist(strsplit(bead_file, split = "[.]"))[1]),
                    sep = "")

  ## default peak postions if none provided
  if (is.na(mef_peaks)) {
    mef_peaks <- list(list(channel = "BL1-H",
                           peaks = c(0, 822, 2114, 5911, 17013, 41837, 145365, 287558)),
                      list(channel = "YL2-H",
                           peaks = c(0, 218, 581, 1963, 6236, 15267, 68766, 181945)))
  }

  ## make a vector of channels for which to produce calibration parameters
  peak_channels <- c()
  for (i in seq_len(length(mef_peaks))) {
    if (mef_peaks[[i]]$channel %in% flu_channels) {
      peak_channels <- c(peak_channels, mef_peaks[[i]]$channel)
    }
  }

  ## Remove events on the boundary that would cause -Inf when log transformed
  bf <- flowCore::boundaryFilter(x = flu_channels, side = "lower")
  bounded_bead_filter <- flowCore::filter(bead_frame, bf)
  bounded_bead_frame <- flowCore::Subset(bead_frame, bounded_bead_filter)

  ## Log10 transform bead data
  log10_bead_frame <- flowCore::transform(bounded_bead_frame,
                                          flowCore::transformList(from = flowCore::colnames(bounded_bead_frame),
                                                                  tfun = log10))

  ## gate to remove debris and aggregates
  singlet_cluster <- flowClust::flowClust(log10_bead_frame,
                                          varNames = c("FSC-H", "SSC-H"),
                                          K = 1,
                                          level = 0.8);
  singlet_beads <- log10_bead_frame[singlet_cluster, ]

  singlet_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = as.data.frame(log10_bead_frame@exprs),
                        ggplot2::aes(`FSC-H`, `SSC-H`),
                        size = 0.2,
                        alpha = 0.01) +
    ggplot2::geom_density_2d(data = as.data.frame(singlet_beads@exprs),
                             ggplot2::aes(`FSC-H`, `SSC-H`), size = 0.1) +
    ggplot2::scale_x_continuous(name = "log10 FSC-H") +
    ggplot2::scale_y_continuous(name = "log10 SSC-H") +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(plot = singlet_plot, filename = paste(out_name,
                                                        "_singlet_beads.pdf",
                                                        sep = ""),
                  width = 60, height = 60, units = "mm")

  calibration_parameters <- c()
  for (i in seq_len(length(peak_channels))) {
    hist_plt <- ggplot2::ggplot() +
      ggplot2::geom_density(data = as.data.frame(log10_bead_frame@exprs),
                            ggplot2::aes(log10_bead_frame[, peak_channels[i]]@exprs),
                            fill = "black", alpha = 0.25, bw = bead_dens_bw) +
      ggplot2::geom_density(data = as.data.frame(singlet_beads@exprs),
                            ggplot2::aes(singlet_beads[, peak_channels[i]]@exprs),
                            fill = "green", alpha = 0.75, bw = bead_dens_bw) +
      ggplot2::scale_x_continuous(paste("log10", peak_channels[i], "(a.u.)")) +
      ggplot2::theme_bw(base_size = 8)

    # find peaks --------------------------------------------------------------
    if (is.na(manual_peaks)) {
      ## find peaks based on density estimate
      dens_d <- stats::density(singlet_beads@exprs[, peak_channels[i]],
                               bw = bead_dens_bw)
      peak_table <- data.frame(dens_d[c("x", "y")])[c(F, diff(diff(dens_d$y) >= 0) < 0), ] ## get peaks and heights
      peaks <- dplyr::top_n(peak_table, n = length(mef_peaks[[i]]$peaks),
                            wt = .data$y)$x ## select only as many peaks as we have defined

      hist_plt <- hist_plt +
        ggplot2::geom_vline(xintercept = peaks,
                            linetype = 2, size = 0.2)

      peaks <- peaks[-1]
    } else {  ## OR peaks are being manually identified
      peaks <- NA
      for (j in seq_len(length(manual_peaks))) {
        if (manual_peaks[[j]]$channel == peak_channels[i]) {
          peaks <- (manual_peaks[[j]]$peaks)[-1]
          break
        }
      }
      if (is.na(peaks)) { ## if we don't have calibration data, skip this channel
        next
      }

      hist_plt <- hist_plt +
        ggplot2::geom_vline(xintercept = peaks, linetype = 2, size = 0.2)
    }

    cal_peaks <- NA
    for (j in seq_len(length(mef_peaks))) {
      if (mef_peaks[[j]]$channel == peak_channels[i]) {
        cal_peaks <- log10(mef_peaks[[j]]$peaks)[-1]
        break
      }
    }
    if (is.na(cal_peaks)) { ## if we don't have calibration data, skip this channel
      next
    }

    model <- stats::nls(cal_peaks ~ m * peaks + b, start = list(m = 1, b = 0))

    calibration_parameters <- rbind(calibration_parameters,
                                    data.frame(flu = peak_channels[i],
                                               m = stats::coef(model)["m"],
                                               b = stats::coef(model)["b"]))

    new_data <- data.frame(peaks = seq(0, 6, len = 100))

    cal_plt <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = peaks, y = cal_peaks)) +
      ggplot2::geom_line(ggplot2::aes(new_data$peaks,
                                      stats::predict(model,
                                                     newdata = new_data)),
                         size = 0.2) +
      ggplot2::scale_x_continuous(name = paste("log10", peak_channels[i],
                                               "(a.u.)"), limits = c(0, 6.5)) +
      ggplot2::scale_y_continuous(name = paste("log10", peak_channels[i],
                                               "(MEF)"), limits = c(0, 6.5)) +
      ggplot2::theme_bw(base_size = 8)

    plt <- gridExtra::arrangeGrob(grobs = list(hist_plt, cal_plt), ncol = 2)
    ggplot2::ggsave(plot = plt, filename = paste(out_name, "_",
                                                 peak_channels[i],
                                                 "_calibration.pdf", sep = ""),
                    width = 130, height = 60, units = "mm")
  }

  return(calibration_parameters)
}

#' Log10 transform and remove -Inf values from \code{flowFrame}
#'
#' @param flow_frame a \code{flowFrame} for processing
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to \code{"BL1-H"}.
#'
#' @return
#'
prep_flow_frame <- function(flow_frame, flu_channels) {
  # log10 transform the values
  # log10_flow_frame <-
  #   flowCore::transform(flow_frame,
  #                       flowCore::transformList(from = flowCore::colnames(flow_frame),
  #                                               tfun = log10))

  log10_flow_frame <-
    flowCore::transform(flow_frame,
                        flowCore::transformList(from = c("FSC-A", "SSC-A", "FSC-H", "SSC-H"),
                                                tfun = log10))

  # Remove -Inf values produced when log10 is applied
  processed_flow_frame <- log10_flow_frame
  for (marker in c("FSC-H", "SSC-H", "SSC-A", flu_channels)) {
    processed_flow_frame <-
      flowCore::Subset(processed_flow_frame,
                       is.finite(flowCore::exprs(processed_flow_frame[, marker]))[, 1])
  }

  return(processed_flow_frame)
}

#' Remove background debris events
#'
#' @param flow_frame a \code{flowFrame} for processing
#' @param pre_cleaned a flag to identify if debris has already been manually
#' gated
#'
#' @return
#'
get_bacteria <- function(flow_frame, pre_cleaned) {
  if (pre_cleaned) {
    best_clusters <- flowClust::flowClust(flow_frame,
                                          varNames = c("FSC-H", "SSC-H"),
                                          K = 1,
                                          level = 0.90);

    bact_flow_frame <- flow_frame[best_clusters, ]

    return(bact_flow_frame)
  } else {

    ## calculate clusters for K=1 and K=2
    all_clusters <- flowClust::flowClust(flow_frame,
                                         varNames = c("FSC-H", "SSC-H"),
                                         K = 1:2,
                                         criterion = "ICL",
                                         level = 0.90);

    ## get the results for the K with the best ICL
    best_clusters <- all_clusters[[all_clusters@index]]
    print(paste(best_clusters@K, "clusters found"))

    if (best_clusters@K == 2) { ## if there are two cluster
      bact_indx <- which(flowClust::getEstimates(best_clusters, flow_frame)$locationsC[, 2]
                         == max(flowClust::getEstimates(best_clusters, flow_frame)$locationsC[, 2]))

      clst_fsc <- flowClust::getEstimates(best_clusters, flow_frame)$locationsC[[bact_indx, 1]]
      clst_ssc <- flowClust::getEstimates(best_clusters, flow_frame)$locationsC[[bact_indx, 2]]

      if ((clst_fsc < 4) & (clst_ssc < 3)) { # TODO: remove hardcoding of thresholds
        print("only debris found")
        bact_flow_frame <- 0
      } else {
        ks <- seq(1:2)
        debris_clusts <- setdiff(ks, bact_indx)
        split_flow_frame <- flowClust::split(flow_frame, best_clusters,
                                             population = list(debris = debris_clusts,
                                                               non_debris = bact_indx))
        bact_flow_frame <- split_flow_frame$non_debris
      }

    }
    else {  ## if there is only one cluster
      clst_fsc <- flowClust::getEstimates(best_clusters, flow_frame)$locationsC[[1]]
      clst_ssc <- flowClust::getEstimates(best_clusters, flow_frame)$locationsC[[2]]
      if ((clst_fsc < 4) & (clst_ssc < 3)) { # TODO: remove hardcoding of thresholds
        print("only debris found")
        bact_flow_frame <- 0
      } else {
        bact_flow_frame <- flow_frame[best_clusters, ]
      }
    }

    return(bact_flow_frame)
  }
}

#' Remove doublet events
#'
#' @param flow_frame a \code{flowFrame} for processing
#'
#' @return
#'
get_singlets <- function(flow_frame) {

  # ## method using function from openCYto package
  # sg <- flowStats::singletGate(flow_frame, area = "SSC-A", height = "SSC-H",
  #                              prediction_level = 0.95)
  # fres <- flowCore::filter(flow_frame, sg)
  # singlet_flow_frame <- flowCore::Subset(flow_frame, fres)

  singlet_flow_frame <- flowCore::Subset(flow_frame,
                                         ((flow_frame[, "SSC-H"]@exprs - flow_frame[, "SSC-A"]@exprs)[,1])^2 < 0.01)

  return(singlet_flow_frame)
}

#' Calibrate fluorescence measurements
#'
#' @param flow_frame a \code{flowFrame} for processing
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to "BL1-H".
#' @param calibration_parameters parameters of a model returned by
#' get_calibration(...)
#'
#' @return
#'
to_mef <- function(flow_frame, flu_channels, calibration_parameters, normalise) {
  out_flow_frame <- flow_frame
  for (fl in flu_channels) {
    fl_idx <- which(calibration_parameters$flu == fl)
    if (length(fl_idx) == 1) {
      linear_trans <- flowCore::linearTransform(transformationId = "Linear-transformation",
                                                a = calibration_parameters[fl_idx, ]$m,
                                                b = calibration_parameters[fl_idx, ]$b)

      out_flow_frame[,paste("calibrated_", fl, sep = "")] <- out_flow_frame[, fl]
      out_flow_frame <- flowCore::transform(out_flow_frame,
                                            flowCore::transformList(paste("calibrated_", fl, sep = ""),
                                                                    linear_trans))

      if(normalise){
        out_flow_frame[, paste("calibrated_normalised_", fl, sep = "")] <- out_flow_frame[, paste("normalised_", fl, sep = "")]
        out_flow_frame <- flowCore::transform(out_flow_frame,
                                              flowCore::transformList(paste("calibrated_normalised_", fl, sep = ""),
                                                                      linear_trans))
      }
    }
  }

  return(out_flow_frame)
}

#' Plot trimming process
#'
#' @param flow_frame the original \code{flowFrame}
#' @param bacteria_flow_frame the \code{flowFrame} with debris removed
#' @param singlet_flow_frame the \code{flowFrame} with debris and doublets removed
#' @param normalised_flow_frame the normalised \code{flowFrame}
#' @param out_flow_frame the calibrated \code{flowFrame} with debris and doublets
#' removed
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to "BL1-H".
#' @param out_name the filename for the outputted \code{flowFrame}
#' @param normalised a Boolean flag to determine whether to normalise
#' @param calibrate a Boolean flag to determine whether to convert fluorescence
#' to MEF values. Requires an .fcs file with named \code{"*beads*.fcs"}.
#' Defaults to \code{FALSE}.
#'
plot_trimming <- function(flow_frame, bacteria_flow_frame, singlet_flow_frame,
                          normalised_flow_frame, out_flow_frame, flu_channels,
                          out_name, normalised, calibrate) {
  ##  This function allows us to take a legend from a plot
  get_legend <- function(myggplot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  plts <- list()

  # plotting theme ----------------------------------------------------------

  apatheme <- ggplot2::theme_bw() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 8),
                   strip.background = ggplot2::element_rect(colour = "white"),
                   axis.text = ggplot2::element_text(size = 8),
                   axis.text.x = ggplot2::element_text(angle = -40,
                                                       vjust = 0.5),
                   axis.title = ggplot2::element_text(size = 8),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.title = ggplot2::element_blank())



  # FSC-H vs SSC-H ----------------------------------------------------------

  plt_main <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(flow_frame[, c("FSC-H", "SSC-H")]@exprs),
                                               size = min(2000, flowCore::nrow(flow_frame))),
                        ggplot2::aes(x = log10(`FSC-H`), y = log10(`SSC-H`),
                                     color = "all_data"),
                        alpha = 0.1) +
    ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(bacteria_flow_frame[, c("FSC-H", "SSC-H")]@exprs),
                                               size = min(2000, flowCore::nrow(bacteria_flow_frame))),
                        ggplot2::aes(x = `FSC-H`, y = `SSC-H`,
                                     color = "bacteria"),
                        alpha = 0.1) +
    ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(singlet_flow_frame[, c("FSC-H", "SSC-H")]@exprs),
                                               size = min(2000, flowCore::nrow(singlet_flow_frame))),
                        ggplot2::aes(x = `FSC-H`, y = `SSC-H`,
                                     color = "single_bacteria"),
                        alpha = 0.1) +
    ggplot2::xlab("log10(FSC-H)") +
    ggplot2::ylab("log10(SSC-H)") +
    ggplot2::xlim(1, 6) +
    ggplot2::ylim(1, 6) + apatheme

  ## Grab the legend to use seperately
  bacteria_legend <- get_legend(plt_main)
  plt_main <- plt_main + ggplot2::theme(legend.position = "none")

  plts[[1]] <- plt_main


  # SSC-H vs SSC-A ----------------------------------------------------------

  plt_single <- ggplot2::ggplot() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(flow_frame[, c("SSC-H", "SSC-A")]@exprs),
                                               size = min(2000, flowCore::nrow(flow_frame))),
                        ggplot2::aes(x = log10(`SSC-H`), y = log10(`SSC-A`),
                                     color = "all_data"),
                        alpha = 0.1) +
    ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(bacteria_flow_frame[, c("SSC-H", "SSC-A")]@exprs),
                                               size = min(2000, flowCore::nrow(bacteria_flow_frame))),
                        ggplot2::aes(x = `SSC-H`, y = (`SSC-A`),
                                     color = "bacteria"),
                        alpha = 0.1) +
    ggplot2::geom_point(data = dplyr::sample_n(as.data.frame(singlet_flow_frame[, c("SSC-H", "SSC-A")]@exprs),
                                               size = min(2000, flowCore::nrow(singlet_flow_frame))),
                        ggplot2::aes(x = `SSC-H`, y = (`SSC-A`),
                                     color = "single_bacteria"),
                        alpha = 0.1) +
    ggplot2::xlab("log10(SSC-H)") +
    ggplot2::ylab("log10(SSC-A)") +
    ggplot2::xlim(1, 6) +
    ggplot2::ylim(1, 6)  + apatheme +
    ggplot2::theme(legend.position = "none")

  plts[[2]] <- plt_single

  plts[[3]] <- bacteria_legend


  # Fluorescence channels ---------------------------------------------------

  ## NOTE: the local is necessary here as R is not good at variable scope.
  ## Without the local, each plot gets overwritten by the final plot in the
  ## loop. Also note the "<<-" to assign the plot outside of the local scope
  flu_legend <- c()
  for (f_count in seq_len(length(flu_channels)))
    local({
      f_count <- f_count
      filt <- flu_channels[f_count]
      plts[[f_count + 3]] <<- ggplot2::ggplot() +
        ggplot2::geom_density(data = as.data.frame(flow_frame[, filt]@exprs),
                              ggplot2::aes(x = log10(flow_frame[, filt]@exprs),
                                           y = ..count.., fill = "all_data"),
                              alpha = 0.5, bw = 0.01) +
        ggplot2::geom_density(data = as.data.frame(bacteria_flow_frame[, filt]@exprs),
                              ggplot2::aes(x = bacteria_flow_frame[, filt]@exprs,
                                           y = ..count.., fill = "bacteria"),
                              alpha = 0.5, bw = 0.01) +
        ggplot2::geom_density(data = as.data.frame(singlet_flow_frame[, filt]@exprs),
                              ggplot2::aes(x = singlet_flow_frame[, filt]@exprs,
                                           y = ..count.., fill = "single_bacteria"),
                              alpha = 0.5, bw = 0.01) +
        ggplot2::xlab(paste("log10(", filt, ")", sep = "")) +
        ggplot2::xlim(0, 6)  + apatheme +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())

      if (normalised) {
        plts[[f_count + 3]] <<- plts[[f_count + 3]] +
          ggplot2::geom_density(data = as.data.frame(normalised_flow_frame[, paste("normalised_", filt, sep="")]@exprs),
                                ggplot2::aes(x = normalised_flow_frame[, paste("normalised_", filt, sep="")]@exprs,
                                             y = ..count.., fill = "normalised"),
                                alpha = 0.5, bw = 0.01)
        if (calibrate) {
          plts[[f_count + 3]] <<- plts[[f_count + 3]] +
            ggplot2::geom_density(data = as.data.frame(out_flow_frame[, paste("calibrated_normalised_", filt, sep="")]@exprs),
                                  ggplot2::aes(x = out_flow_frame[, paste("calibrated_normalised_", filt, sep="")]@exprs,
                                               y = ..count.., fill = "calibrated_normalised"),
                                  alpha = 0.5, bw = 0.01)
        }
      }
      if (calibrate) {
        plts[[f_count + 3]] <<- plts[[f_count + 3]] +
          ggplot2::geom_density(data = as.data.frame(out_flow_frame[, paste("calibrated_", filt, sep="")]@exprs),
                                ggplot2::aes(x = out_flow_frame[, paste("calibrated_", filt, sep="")]@exprs,
                                             y = ..count.., fill = "calibrated"),
                                alpha = 0.5, bw = 0.01)
      }

      ## Grab the legend to use seperately
      flu_legend <<- get_legend(plts[[f_count + 3]])
      plts[[f_count + 3]] <<- plts[[f_count + 3]] +
        ggplot2::theme(legend.position = "none")
    })

  plts[[4 + length(flu_channels)]] <- flu_legend


  # Assemble plots ----------------------------------------------------------

  plt <- gridExtra::arrangeGrob(grobs = plts, ncol = 3)
  title <- grid::textGrob(paste("Trimming of flow data to remove background and doublets:\n",
                                flowCore::identifier(flow_frame)),
                          gp = grid::gpar(fontsize = 10))
  padding <- grid::unit(5, "mm")
  plt <- gtable::gtable_add_rows(plt,
                                 heights = grid::grobHeight(title) + padding,
                                 pos = 0)
  plt <- gtable::gtable_add_grob(plt, title, 1, 1, 1, ncol(plt))

  print(paste("Plotting processed flowFrame ",
              flowCore::identifier(flow_frame)))

  ggplot2::ggsave(filename = paste(dirname(out_name),
                                   gsub(".fcs",
                                        "_processed.pdf",
                                        basename(out_name)),
                                   sep = "/"),
                  plot = plt,
                  width = 150,
                  height = 50 * ceiling((length(plts) / 3)) + 20,
                  units = "mm")
}
