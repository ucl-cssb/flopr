#' Plate reader normalisation and calibration
#'
#' @param data_csv path to a .csv file containing parsed plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param od_name the column name for the optical density data
#' @param flu_names the column names for the fluorescence data
#' @param af_model model used to fit negative control autofluorescence.
#' For now these include "polynomial", "invers_poly", "exponential", "spline" and "loess".
#' @param to_MEFL a Boolean to determine whether to attempt to convert OD and
#' GFP reading to calibrated units
#' @param flu_gains if to_MEFL=T, the gain values at which the fluorophores
#' specified in flu_names was recorded. If there isn't calibration data for a
#' fluorophore, do not speficy a gain value
#' @param conversion_factors_csv if to_MEFL=T, path of the csv file containing
#' conversion factors from plate reader calibration
#'
#' @return a data.frame with columns for raw plate reader data, normalised data
#' and, if to_MEFL = T, calibrated OD and GFP data
#' @export
#' @importFrom rlang .data
#'
#' @examples
process_plate <- function(data_csv, blank_well = "A1", neg_well = "A2",
                          od_name = "OD", flu_names = c("GFP"),
                          af_model = "spline", to_MEFL = F,
                          flu_gains, conversion_factors_csv) {

  pr_data <- utils::read.csv(data_csv, check.names = F)

  od_norm_pr_data <- od_norm(pr_data, blank_well, od_name)

  plt_od <- ggplot2::ggplot(od_norm_pr_data) +
    ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data[[od_name]],
                                    colour = "raw"), size = 0.2) +
    ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$normalised_OD,
                                    colour = "normalised"), size = 0.2) +
    ggplot2::scale_x_continuous("time") +
    ggplot2::scale_colour_discrete("") +
    ggplot2::facet_grid(row~column) +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(filename = gsub(".csv", "_OD.pdf", data_csv),
                  plot = plt_od, height = 160,
                  width = 240, units = "mm")

  flu_norm_pr_data <- od_norm_pr_data
  if(all(!is.na(flu_names))){
    if(length(flu_names) >= 1){
      for (flu_idx in seq_len(length(flu_names))) {
        flu_norm_pr_data <- flu_norm(flu_norm_pr_data, neg_well, blank_well,
                                     flu_names[flu_idx], af_model, data_csv)

        plt_flu <- ggplot2::ggplot(flu_norm_pr_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data[[flu_names[flu_idx]]],
                                          colour = "raw"), size = 0.2) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                          y = .data[[paste("normalised_",
                                                           flu_names[flu_idx],
                                                           sep = "")]],
                                          colour = "normalised"),
                             size = 0.2) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_colour_discrete("") +
          ggplot2::facet_grid(row~column) +
          ggplot2::theme_bw(base_size = 8)
        ggplot2::ggsave(filename = gsub(".csv",
                                        paste("_", flu_names[flu_idx], ".pdf", sep = ""),
                                        data_csv),
                        plot = plt_flu, height = 160,
                        width = 240, units = "mm")
      }
    }
  }

  out_data <- flu_norm_pr_data

  if (to_MEFL) {
    out_data <- calibrate_od(out_data, od_name,
                             conversion_factors_csv)

    if(all(!is.na(flu_names))){
      for (flu_idx in seq_len(length(flu_names))) {
        if(length(flu_gains) >= flu_idx){
          out_data <- calibrate_flu(out_data, flu_names[flu_idx],
                                    flu_gains[flu_idx], od_name,
                                    conversion_factors_csv)
        }
        else {break}
      }
    }
  }

  utils::write.csv(x = out_data,
                   file = gsub(".csv", "_processed.csv", data_csv),
                   row.names = FALSE)
  return(out_data)
}


#' Normalisation absorbance against blank well
#'
#' @param pr_data a long data.frame containing you plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param od_name the column name for the optical density data
#'
#' @return an updated data.frame with an additional column "normalised_OD"
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
od_norm <- function(pr_data, blank_well, od_name) {
  pr_data$normalised_OD <- pr_data[, od_name]

  # remove background absorbance signal -------------------------------------

  pr_data <- pr_data %>%
    dplyr::group_by(.data$time) %>%
    dplyr::mutate(normalised_OD = .data$normalised_OD -
                    mean(.data$normalised_OD[.data$well %in% blank_well]))

  return(as.data.frame(pr_data))
}


#' Normalise fluorescence against negative well
#'
#' @param pr_data a long data.frame containing you plate reader data with OD
#' normalised
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param blank_well the well coordinates of a media blank
#' @param flu_name the column name of the fluorescence chanel to normalise
#' @param af_model model used to fit negative control autofluorescence.
#' For now these include "polynomial", "invers_poly", "exponential", "spline" or "loess".
#' @param data_csv path to the original data. Used for saving normalisation curve plots.
#'
#' @return
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
flu_norm <- function(pr_data, neg_well, blank_well, flu_name, af_model, data_csv) {
  pr_data$v1 <- pr_data[, flu_name]

  # fit autofluorescence model to negative control --------------------------

  negative_data <- pr_data %>% dplyr::filter(.data$well %in% neg_well)

  if (af_model == "polynomial") {
    model <- stats::nls(v1 ~ (a * normalised_OD + b * normalised_OD ^ 2 + c),
                        start = c(a = 1, b = 1, c = 1), data = negative_data)

  } else if (af_model == "inverse_poly") {
    model <- stats::nls(normalised_OD ~ (a * v1 + b * v1 ^ 2 + c),
                        start = c(a = 1, b = 1, c = 1), data = negative_data)

  } else if (af_model == "exponential") {
    ## ae^(bx) + c
    ## intial parameter estimation
    model_0 <- stats::lm(log(v1) ~ normalised_OD, data = negative_data)
    start <- list(a = exp(stats::coef(model_0)[1]),
                  b = stats::coef(model_0)[2],
                  c = -1)

    model <- stats::nls(v1 ~ (a * exp(b * normalised_OD) + c),
                        start = start, data = negative_data)

  } else if (af_model == "bi_exponential") {
    ## a exp(bx) + c exp(dx) + e

    model_0 <- stats::lm(log(v1) ~ normalised_OD, data = negative_data)
    start <- list(a = exp(stats::coef(model_0)[1]*0.2),
                  b = stats::coef(model_0)[2]*0.2,
                  c = exp(stats::coef(model_0)[1])*0.8,
                  d = stats::coef(model_0)[2]*0.8,
                  e = 1)

    model <- stats::nls(v1 ~ (a * exp(b * normalised_OD) +
                                c * exp(d * normalised_OD) + e),
                        start = start,
                        data = negative_data)

  } else if (af_model == "linear_exponential") {
    ## ax + be^cx + d

    model_01 <- stats::lm(v1 ~ normalised_OD, data = negative_data)
    model_02 <- stats::lm(log(v1) ~ normalised_OD, data = negative_data)
    start <- list(a = stats::coef(model_01)[2],
                  b = exp(stats::coef(model_02)[1]),
                  c = stats::coef(model_02)[2],
                  d = 1)

    model <- stats::nls(v1 ~ (a * normalised_OD +
                                b * exp(c * normalised_OD) + d),
                        start = start,
                        data = negative_data)

  } else if (af_model == "power") {
    ## ax^b + c
    model_0 <- stats::lm(log(v1) ~ log(normalised_OD), data = negative_data)
    start <- list(a = exp(stats::coef(model_0)[1]),
                  b = stats::coef(model_0)[2],
                  c = 1)

    model <- stats::nls(v1 ~ (a * normalised_OD ^ b + c),
                        start = start, data = negative_data)

  } else if (af_model == "linear_power") {
    ## ax + bx^c + d
    model_01 <- stats::lm(v1 ~ normalised_OD, data = negative_data)
    model_02 <- stats::lm(log(v1) ~ log(normalised_OD), data = negative_data)
    start <- list(a = stats::coef(model_01)[2],
                  b = exp(stats::coef(model_02)[1]),
                  c = stats::coef(model_02)[2],
                  d = 1)

    model <- stats::nls(v1 ~ (a * normalised_OD + b * normalised_OD ^ c + d),
                        start = start,
                        data = negative_data)
  } else if (af_model == "loess") {
    model <- stats::loess(v1 ~ normalised_OD,
                          data = negative_data,
                          span = 0.5)
  } else if (af_model == "spline") {
    model <- mgcv::gam(v1 ~ s(normalised_OD), data = negative_data)
  }

  # plot model fit curves ---------------------------------------------------

  if (af_model == "polynomial" | af_model == "power" |
      af_model == "exponential" | af_model == "bi_exponential" |
      af_model == "linear_exponential" | af_model == "linear_power" |
      af_model == "loess") {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                      y = stats::predict(model,
                                                         negative_data))) +
      ggplot2::geom_point(ggplot2::aes(x = negative_data$normalised_OD,
                                       y = negative_data$v1)) +
      ggplot2::scale_x_continuous("normalised_OD") +
      ggplot2::scale_y_continuous(flu_name) +
      ggplot2::theme_bw()
  } else if (af_model == "spline") {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                      y = mgcv::predict.gam(model, negative_data))) +
      ggplot2::geom_point(ggplot2::aes(x = negative_data$normalised_OD,
                                       y = negative_data$v1)) +
      ggplot2::scale_x_continuous("normalised_OD") +
      ggplot2::scale_y_continuous(flu_name) +
      ggplot2::theme_bw()
  }
  else if (af_model == "inverse_poly") {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                      y = ((- (stats::coef(model)[1]) +
                                              sqrt((stats::coef(model)[1]) ^ 2 -
                                                     4 *
                                                     (stats::coef(model)[2]) *
                                                     (stats::coef(model)[3]) +
                                                     4 *
                                                     (stats::coef(model)[2]) *
                                                     negative_data$normalised_OD)) /
                                             (2 * (stats::coef(model)[2]))))) +
      ggplot2::geom_point(ggplot2::aes(y = negative_data$v1,
                                       x = negative_data$normalised_OD)) +
      ggplot2::scale_x_continuous("normalised_OD") +
      ggplot2::scale_y_continuous(flu_name) +
      ggplot2::theme_bw()
  }
  ggplot2::ggsave(filename = gsub(".csv",
                                  paste("_norm-curve_", flu_name, ".pdf", sep = ""),
                                  data_csv),
                  plot = plt)

  # normalise fluorescence data ---------------------------------------------

  if (af_model == "polynomial" | af_model == "power" |
      af_model == "exponential" | af_model == "bi_exponential" |
      af_model == "linear_exponential" | af_model == "linear_power" |
      af_model == "loess") {
    pr_data$v1 <- pr_data$v1 - stats::predict(model, pr_data)
  } else if (af_model == "spline") {
    pr_data$v1 <- pr_data$v1 - mgcv::predict.gam(model, pr_data)
  }
  else if (af_model == "inverse_poly") {
    pr_data$v1 <- pr_data$v1 - ((- (stats::coef(model)[1]) +
                                   sqrt((stats::coef(model)[1]) ^ 2 - 4 *
                                          (stats::coef(model)[2]) *
                                          (stats::coef(model)[3]) +
                                          4 * (stats::coef(model)[2]) *
                                          pr_data$normalised_OD)) /
                                  (2 * (stats::coef(model)[2])))
  }

  # rename normalised fluorescence column
  names(pr_data)[ncol(pr_data)] <- paste("normalised_", flu_name, sep = "")

  return(pr_data)
}

#' Convert arbitrary absorbance units to calibrated units
#'
#' @param pr_data a data.frame of parsed plate reader data
#' @param od_name the column name for the optical density data
#' @param conversion_factors_csv if to_MEFL=T, path of the csv file containing
#' conversion factors from plate reader calibration
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return
calibrate_od <- function(pr_data, od_name, conversion_factors_csv) {
  conversion_factors <- utils::read.csv(conversion_factors_csv)

  # Get conversion factor for OD --------------------------------------------

  od_cf <- unlist(conversion_factors %>%
                    dplyr::filter(.data$measure == od_name) %>%
                    dplyr::select(.data$cf))

  pr_data$calibrated_OD <- pr_data$normalised_OD / od_cf

  return(pr_data)
}


#' Convert arbitrary fluorescence units to calibrated units
#'
#' @param pr_data a data.frame of parsed plate reader data
#' @param flu_name name of fluorophore to be calibrated
#' @param flu_gain gain at which the fluorophore was measured
#' @param od_name the column name for the optical density data
#' @param conversion_factors_csv if to_MEFL=T, path of the csv file containing
#' conversion factors from plate reader calibration
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return
calibrate_flu <- function(pr_data, flu_name, flu_gain, od_name, conversion_factors_csv) {
  conversion_factors <- utils::read.csv(conversion_factors_csv)


  # Get conversion factor for fluorophore ------------------------------------
  flu_cfs <- conversion_factors %>%
    dplyr::filter(.data$fluorophore == flu_name)

  # if there is no calibration for the fluorophore
  if(nrow(flu_cfs) == 0) {
    return(pr_data)
  }

  tryCatch(this_cf <- flu_cfs[which(flu_cfs$measure == paste(flu_name, flu_gain)),]$cf,
           finally = this_cf <- NA)


  # if a conversion factor doesn't exist at the measured gain, try a --------
  if(is.na(this_cf)){
    flu_cfs$gain <- as.numeric(gsub(paste(flu_name, " ", sep=""), "", flu_cfs$measure))

    # Fit cf to Gain relation to get cf for specific gain ---------------------
    model <- stats::lm(log10(cf) ~ poly(gain, 2), data = flu_cfs)
    this_cf <- 10 ^ stats::predict(model, data.frame(gain = flu_gain))
    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = flu_cfs$gain,
                                      y = 10^stats::predict(model, flu_cfs))) +
      ggplot2::geom_point(ggplot2::aes(x = flu_cfs$gain,
                                       y = flu_cfs$cf)) +
      ggplot2::geom_vline(xintercept = flu_gain, linetype = 2) +
      ggplot2::geom_hline(yintercept = 10 ^ stats::predict(model, data.frame(gain = flu_gain)),
                          linetype = 2) +
      ggplot2::geom_point(ggplot2::aes(x = flu_gain,
                                       y = 10 ^ stats::predict(model, data.frame(gain = flu_gain))),
                          colour = "red", shape = 1, size = 2) +
      ggplot2::scale_x_continuous("Gain") +
      ggplot2::scale_y_continuous("Conversion factor (MEFL/a.u.)",
                                  trans = "log10") +
      ggplot2::theme_bw()
  }

  pr_data[,paste("calibrated_", flu_name, sep="")] <-
    (pr_data[,paste("normalised_", flu_name, sep="")] / this_cf)

  return(pr_data)
}


#' Generate Conversion Factors
#'
#' @param calibration_csv path of a .csv file of your calibration data
#'
#' @return Saves a CSV file with columns: cf (conversion factor), beta (pipetting error parameter), calibrant, fluorophore, and measure. Each row corresponds to a fitted conversion factor for a calibrant/fluorophore/measure combination.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
generate_cfs <- function(calibration_csv) {
  calibration_data <- utils::read.csv(calibration_csv, header = T, check.names = F)
  fit_values <- c()

  # get types of measure ----------------------------------------------------
  well_idx <- which(names(calibration_data) == "well")
  row_idx <- which(names(calibration_data) == "row")
  measures <- names(calibration_data)[(well_idx+1):(row_idx-1)]

  long_calibration_data <- tidyr::pivot_longer(calibration_data,
                                               cols = tidyselect::all_of(measures),
                                               names_to = 'measure')

  for(calib in unique(long_calibration_data$calibrant)){
    if(calib == "" | is.na(calib)){
      next
    }

    calib_output <- c()
    processed_data <- c()

    # get values only for this calibrant
    temp_calib_values <- long_calibration_data %>%
      dplyr::filter(.data$calibrant == calib)

    for(meas in unique(temp_calib_values$measure)){
      # get values only for this measure
      temp_calib_meas_values <- temp_calib_values %>%
        dplyr::filter(.data$measure == meas) %>%
        dplyr::arrange(dplyr::desc(.data$concentration)) %>%
        dplyr::group_by(.data$calibrant) %>%
        dplyr::mutate(norm_value = .data$value - mean(.data$value[.data$concentration == 0])) %>%
        dplyr::ungroup()

      # remove saturated values -------------------------------------------------
      non_sat_values <- c()

      # get concentrations at which we have measurements
      concentrations <- sort(unique(temp_calib_meas_values$concentration), decreasing = T)

      # calculate fold dilution used
      fold_dilution <- concentrations[2] / concentrations[3]

      # work out dilution order from concentrations
      temp_calib_meas_values$dilution_ratio <- 1 / fold_dilution
      temp_calib_meas_values$max_concentration <- max(concentrations)
      temp_calib_meas_values$dilution_idx <- as.integer(round(- log(temp_calib_meas_values$max_concentration / temp_calib_meas_values$concentration) / log(temp_calib_meas_values$dilution_ratio)))

      blank_mean <- mean(temp_calib_meas_values[temp_calib_meas_values$concentration == 0,]$value, na.rm = T)
      blank_sd <- stats::sd(temp_calib_meas_values[temp_calib_meas_values$concentration == 0,]$value, na.rm = T)

      saturated_list <- c()
      for(rplct in unique(temp_calib_meas_values$replicate)){
        values <- temp_calib_meas_values[temp_calib_meas_values$replicate == rplct,]$value

        # set up array to track saturated values
        rep_saturated <- rep('FALSE', length(values))

        top_down <- TRUE
        for(i in 1:(length(rep_saturated)-1)){
          if(is.na(values[i])){
            rep_saturated[i] <- 'NA'
          }
          else{
            if(top_down){
              if((i < length(rep_saturated)) & (values[i] <= (values[i+1] * (fold_dilution * 0.75)))){
                rep_saturated[i] <- 'TOP'
              } else {
                top_down <- FALSE
              }
            }
            else {
              if((i > 1) & (values[i] >= (values[i-1] / (fold_dilution * 0.75)))){
                rep_saturated[i] <- 'BOTTOM'
              }
            }

            if(values[i] <= (blank_mean + 2*blank_sd)){
              rep_saturated[i] <- 'BLANK'
            }
          }
        }

        df_rep_saturated <- data.frame(saturated=rep_saturated, concentration=concentrations, replicate=rplct)
        saturated_list <- dplyr::bind_rows(saturated_list, df_rep_saturated)
      }
      temp_calib_meas_values <- dplyr::full_join(temp_calib_meas_values, saturated_list)
      non_sat_values <- rbind(non_sat_values, temp_calib_meas_values %>% filter(.data$saturated == 'FALSE'))
      processed_data <- rbind(processed_data, temp_calib_meas_values)

      # calculate mean of replicates -------------

      summ_values <- non_sat_values %>%
        dplyr::group_by(.data$calibrant, .data$fluorophore, .data$media, .data$measure,
                        .data$concentration, .data$dilution_ratio,
                        .data$max_concentration, .data$dilution_idx, .drop = F) %>%
        dplyr::summarise(mean_value = mean(.data$value, na.rm = TRUE)) %>%
        dplyr::filter(!is.na(.data$concentration))

      # normalise data ----------------------------------------------------------

      norm_values <- summ_values %>%
        dplyr::group_by(.data$calibrant) %>%
        dplyr::mutate(norm_value = .data$mean_value - .data$mean_value[.data$concentration == 0]) %>%
        dplyr::filter(.data$concentration != 0) %>%
        dplyr::ungroup()

      if(nrow(norm_values) < 3){
        next
      }

      # fit pipetting error model for conversion factors ------------------------
      # error model from Beal et al. 2019 bioRxiv

      error_func <- function(x){

        data <- norm_values
        cf <- x[1]
        beta <- x[2]
        error <- 0

        for(i in data$dilution_idx){
          data_i <- data[data$dilution_idx == i,]

          b_i <- data_i$max_concentration * (1 - data_i$dilution_ratio - beta) *
            (data_i$dilution_ratio + beta) ^ (data_i$dilution_idx - 1)

          e_i <- abs(log10(cf * b_i / data_i$norm_value))^2
          error <- error + e_i
        }

        return(error)
      }

      for(init_cf in 10^seq(-8, -16)){
        res <- stats::optim(c(init_cf, 0), error_func)

        if(res$convergence == 0){
          new_fit <- data.frame(cf = res$par[1], beta = res$par[2],
                                calibrant = calib,
                                fluorophore = norm_values$fluorophore[1],
                                measure = meas)
          fit_values <- rbind(fit_values, new_fit)
          break
        }
      }

      calib_output <- rbind(calib_output, dplyr::full_join(norm_values, new_fit))
    }

    # plot the mean normalized values -----------------------------------------
    plt <- ggplot2::ggplot() +
      ggplot2::geom_point(data = processed_data %>%
                            dplyr::filter(.data$calibrant == calib),
                          ggplot2::aes(x = .data$dilution_idx,
                                       y = .data$norm_value),
                          colour = 'red') +
      ggplot2::geom_point(data = calib_output %>%
                            dplyr::filter(.data$calibrant == calib),
                          ggplot2::aes(x = .data$dilution_idx,
                                       y = .data$norm_value)) +
      ggplot2::geom_line(data = calib_output %>%
                           dplyr::filter(.data$calibrant == calib),
                         ggplot2::aes(x = .data$dilution_idx,
                                      y = .data$cf * .data$max_concentration *
                                        (1 - .data$dilution_ratio - .data$beta) *
                                        (.data$dilution_ratio + .data$beta) ^ (.data$dilution_idx - 1))) +
      ggplot2::scale_y_continuous("Normalised measurement", trans = "log10") +
      ggplot2::scale_x_continuous("Dilution index") +
      ggplot2::facet_wrap(~measure) +
      ggplot2::theme_bw(base_size = 8)

    # calculate plot width and height by how many facets are created
    n_facets <- length(unique(processed_data$measure))
    facets_per_row <- ceiling(sqrt(n_facets))
    num_rows <- ceiling(n_facets / facets_per_row)
    plot_width <- 30 * facets_per_row
    plot_height <- 30 * num_rows
    ggplot2::ggsave(gsub(".csv", paste("_", calib, "_cfs.pdf", sep = ""), calibration_csv), plot = plt,
                    width = plot_width, height = plot_height, units = "mm")
  }

  utils::write.csv(fit_values, gsub(".csv", "_cfs.csv", calibration_csv), row.names = FALSE)
  return(fit_values)
}
