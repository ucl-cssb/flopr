### VISUALISATION FOR process_fcs()/process_fcs_dir()'S TRIMMING STAGES
### split out from flowCytometryFunctions.R (code review backlog 013)

#' Plot trimming process
#'
#' @param flow_frame the original \code{flowFrame}
#' @param bacteria_flow_frame the \code{flowFrame} with debris removed
#' @param singlet_flow_frame the \code{flowFrame} with debris and doublets removed
#' @param normalised_flow_frame the normalised \code{flowFrame}
#' @param calibrated_flow_frame the calibrated \code{flowFrame} with debris and doublets
#' removed
#' @param flu_channels a vector of strings of the fluorescence channels to keep
#' in the processed data and plotting. Defaults to "BL1-H".
#' @param out_name the filename for the outputted \code{flowFrame}
#' @param normalised a Boolean flag to determine whether to normalise
#' @param calibrate a Boolean flag to determine whether to convert fluorescence
#' to MEF values. Requires an .fcs file with named \code{"*beads*.fcs"}.
#' Defaults to \code{FALSE}.
#'
#' @noRd
plot_trimming <-
  function(flow_frame,
           bacteria_flow_frame,
           singlet_flow_frame,
           normalised_flow_frame,
           calibrated_flow_frame,
           flu_channels,
           out_name,
           normalised,
           calibrate) {
    ##  This function allows us to take a legend from a plot
    get_legend <- function(myggplot) {
      tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
      leg <-
        which(sapply(tmp$grobs, function(x)
          x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }

    plts <- list()

    # plotting theme ----------------------------------------------------------

    apatheme <- ggplot2::theme_bw() +
      ggplot2::theme(
        strip.text.x = ggplot2::element_text(size = 8),
        strip.background = ggplot2::element_rect(colour = "white"),
        axis.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(angle = -40,
                                            vjust = 0.5),
        axis.title = ggplot2::element_text(size = 8),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(),
        legend.title = ggplot2::element_blank()
      )


    # FSC-H vs SSC-H ----------------------------------------------------------

    plt_main <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = dplyr::sample_n(
          as.data.frame(flow_frame[, c("FSC-H", "SSC-H")]@exprs),
          size = min(2000, flowCore::nrow(flow_frame))
        ),
        ggplot2::aes(
          x = log10(`FSC-H`),
          y = log10(`SSC-H`),
          color = "all_data"
        ),
        alpha = 0.1
      ) +
      ggplot2::geom_point(
        data = dplyr::sample_n(
          as.data.frame(bacteria_flow_frame[, c("FSC-H", "SSC-H")]@exprs),
          size = min(2000, flowCore::nrow(bacteria_flow_frame))
        ),
        ggplot2::aes(x = `FSC-H`, y = `SSC-H`,
                     color = "bacteria"),
        alpha = 0.1
      ) +
      ggplot2::geom_point(
        data = dplyr::sample_n(
          as.data.frame(singlet_flow_frame[, c("FSC-H", "SSC-H")]@exprs),
          size = min(2000, flowCore::nrow(singlet_flow_frame))
        ),
        ggplot2::aes(x = `FSC-H`, y = `SSC-H`,
                     color = "single_bacteria"),
        alpha = 0.1
      ) +
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
      ggplot2::geom_point(
        data = dplyr::sample_n(
          as.data.frame(flow_frame[, c("SSC-H", "SSC-A")]@exprs),
          size = min(2000, flowCore::nrow(flow_frame))
        ),
        ggplot2::aes(
          x = log10(`SSC-H`),
          y = log10(`SSC-A`),
          color = "all_data"
        ),
        alpha = 0.1
      ) +
      ggplot2::geom_point(
        data = dplyr::sample_n(
          as.data.frame(bacteria_flow_frame[, c("SSC-H", "SSC-A")]@exprs),
          size = min(2000, flowCore::nrow(bacteria_flow_frame))
        ),
        ggplot2::aes(
          x = `SSC-H`,
          y = (`SSC-A`),
          color = "bacteria"
        ),
        alpha = 0.1
      ) +
      ggplot2::geom_point(
        data = dplyr::sample_n(
          as.data.frame(singlet_flow_frame[, c("SSC-H", "SSC-A")]@exprs),
          size = min(2000, flowCore::nrow(singlet_flow_frame))
        ),
        ggplot2::aes(
          x = `SSC-H`,
          y = (`SSC-A`),
          color = "single_bacteria"
        ),
        alpha = 0.1
      ) +
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
        flu_plt <- ggplot2::ggplot() +
          ggplot2::geom_density(
            data = as.data.frame(flow_frame[, filt]@exprs),
            ggplot2::aes(
              x = flow_frame[, filt]@exprs,
              y = ggplot2::after_stat(count),
              fill = "all_data",
              colour = "all_data"
            ),
            bw = 0.05,
            alpha = 0.4
          ) +
          ggplot2::geom_density(
            data = as.data.frame(bacteria_flow_frame[, filt]@exprs),
            ggplot2::aes(
              x = bacteria_flow_frame[, filt]@exprs,
              y = ggplot2::after_stat(count),
              fill = "bacteria",
              colour = "bacteria"
            ),
            bw = 0.05,
            alpha = 0.4
          ) +
          ggplot2::geom_density(
            data = as.data.frame(singlet_flow_frame[, filt]@exprs),
            ggplot2::aes(
              x = singlet_flow_frame[, filt]@exprs,
              y = ggplot2::after_stat(count),
              fill = "single_bacteria",
              colour = "single_bacteria"
            ),
            bw = 0.05,
            alpha = 0.4
          ) +
          ggplot2::scale_colour_discrete(guide = F) +
          apatheme +
          ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )

        if (normalised) {
          flu_plt <- flu_plt +
            ggplot2::geom_density(
              data = as.data.frame(normalised_flow_frame[, paste("normalised_", filt, sep =
                                                                   "")]@exprs),
              ggplot2::aes(
                x = normalised_flow_frame[, paste("normalised_", filt, sep = "")]@exprs,
                y = ggplot2::after_stat(count),
                fill = "normalised",
                colour = "normalised"
              ),
              bw = 0.05,
              alpha = 0.4
            )
          if (calibrate) {
            try(flu_plt <- flu_plt +
                  ggplot2::geom_density(
                    data = as.data.frame(calibrated_flow_frame[, paste("calibrated_normalised_",
                                                                       filt,
                                                                       sep = "")]@exprs),
                    ggplot2::aes(
                      x = calibrated_flow_frame[, paste("calibrated_normalised_",
                                                        filt,
                                                        sep = "")]@exprs,
                      y = ggplot2::after_stat(count),
                      fill = "calibrated_normalised",
                      colour = "calibrated_normalised"
                    ),
                    bw = 0.05,
                    alpha = 0.4
                  ))
          }
        }
        if (calibrate) {
          try(flu_plt <- flu_plt +
                ggplot2::geom_density(
                  data = as.data.frame(calibrated_flow_frame[, paste("calibrated_",
                                                                     filt,
                                                                     sep = "")]@exprs),
                  ggplot2::aes(
                    x = calibrated_flow_frame[, paste("calibrated_",
                                                      filt,
                                                      sep = "")]@exprs,
                    y = ggplot2::after_stat(count),
                    fill = "calibrated",
                    colour = "calibrated"
                  ),
                  bw = 0.05,
                  alpha = 0.4
                ))
        }

        ## Grab the legend to use separately
        flu_legend <<- get_legend(flu_plt)
        flu_plt_main <- flu_plt +
          ggplot2::scale_x_continuous(filt,
                                      trans = "log10",
                                      limits = c(1e0, 1e6)) +
          ggplot2::theme(legend.position = "none")

        plts[[f_count + 3]] <<- flu_plt_main
      })

    plts[[length(flu_channels) + 4]] <- flu_legend


    # Assemble plots ----------------------------------------------------------

    plt <- gridExtra::arrangeGrob(grobs = plts, ncol = 3)
    title <-
      grid::textGrob(
        paste(
          "Trimming of flow data to remove background and doublets:\n",
          flowCore::identifier(flow_frame)
        ),
        gp = grid::gpar(fontsize = 10)
      )
    padding <- grid::unit(5, "mm")
    plt <- gtable::gtable_add_rows(plt,
                                   heights = grid::grobHeight(title) + padding,
                                   pos = 0)
    plt <- gtable::gtable_add_grob(plt, title, 1, 1, 1, ncol(plt))

    print(paste(
      "Plotting processed flowFrame ",
      flowCore::identifier(flow_frame)
    ))

    ggplot2::ggsave(
      filename = paste(
        dirname(out_name),
        gsub(".fcs",
             "_processed.pdf",
             basename(out_name)),
        sep = "/"
      ),
      plot = plt,
      width = 150,
      height = 50 * ceiling((length(plts) / 3)) + 20,
      units = "mm"
    )
  }
