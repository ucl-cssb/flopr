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
                          af_model = "polynomial", to_MEFL = F,
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

  out_data <- flu_norm_pr_data

  if (to_MEFL) {
    out_data <- calibrate_od(out_data, od_name,
                             conversion_factors_csv)

    for (flu_idx in seq_len(length(flu_names))) {
      if(length(flu_gains) >= flu_idx){
        out_data <- calibrate_flu(out_data, flu_names[flu_idx],
                                  flu_gains[flu_idx], od_name,
                                  conversion_factors_csv)
      }
      else {break}
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
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
generate_cfs <- function(calibration_csv) {
  calibration_data <- utils::read.csv(calibration_csv, header = T, check.names = F)


  # get types of measure ----------------------------------------------------
  well_idx <- which(names(calibration_data) == "well")
  row_idx <- which(names(calibration_data) == "row")
  measures <- names(calibration_data)[(well_idx+1):(row_idx-1)]


  # remove saturated values -------------------------------------------------
  # using similar approach to Beal et al. 2019 bioRxiv to assess validitiy of measurements

  non_sat_values <- c()
  for(calib in unique(calibration_data$calibrant)){
    temp_calib_values <- calibration_data %>%
      dplyr::filter(.data$calibrant == calib)

    concentrations <- sort(unique(temp_calib_values$concentration))
    fold_dilution <- concentrations[3] / concentrations[2]

    high_saturation_threshold <- fold_dilution * 0.75

    temp_calib_values$dilution_ratio <- 1 / fold_dilution
    temp_calib_values$max_concentration <- max(concentrations)
    temp_calib_values$dilution_idx <- - log(temp_calib_values$max_concentration / temp_calib_values$concentration) / log(temp_calib_values$dilution_ratio)

    for(meas in measures){
      blank_mean <- mean(temp_calib_values[temp_calib_values$concentration == 0,][[meas]], na.rm = T)
      blank_sd <- stats::sd(temp_calib_values[temp_calib_values$concentration == 0,][[meas]], na.rm = T)

      for(rplct in unique(temp_calib_values$replicate)){
        prev_value <- 0
        for(conc in concentrations){

          this_value <- temp_calib_values[temp_calib_values$concentration == conc & temp_calib_values$replicate == rplct,][[meas]]
          if(is.na(this_value)){ next }

          if(conc != 0){
            ## check high saturation
            if(this_value <= prev_value * high_saturation_threshold){
              temp_calib_values[temp_calib_values$concentration == conc & temp_calib_values$replicate == rplct, meas] <- NA
            }
            ## check low saturation
            if(this_value <= blank_mean + 2 * blank_sd){
              temp_calib_values[temp_calib_values$concentration == conc & temp_calib_values$replicate == rplct, meas] <- NA
            }
          }
          prev_value <- this_value
        }
      }
    }
    non_sat_values <- rbind(non_sat_values, temp_calib_values)
  }


  # calculate mean of 4 replicates -------------
  #
  summ_values <- non_sat_values %>%
    dplyr::group_by(.data$calibrant, .data$fluorophore, .data$media,
                    .data$concentration, .data$dilution_ratio,
                    .data$max_concentration, .data$dilution_idx, .drop = F) %>%
    dplyr::summarise_at(measures, mean, na.rm = TRUE) %>%
    dplyr::filter(!is.na(.data$concentration))


  # normalise data ----------------------------------------------------------

  norm_values <- summ_values
  for(meas in measures){
    norm_values <- norm_values %>%
      dplyr::group_by(.data$calibrant) %>%
      dplyr::mutate({{meas}} := .data[[meas]] -
                      .data[[meas]][.data$concentration == 0])
  }
  norm_values <- norm_values %>% dplyr::filter(.data$concentration != 0)


  # fit pipetting error model for conversion factors ------------------------
  # error model from Beal et al. 2019 bioRxiv

  long_values <- stats::na.omit(norm_values %>%
                                  tidyr::pivot_longer(tidyselect::all_of(measures),
                                                      names_to = "measure",
                                                      values_to = "normalised_value"))

  fit_values <- c()
  for(calib in unique(long_values$calibrant)){
    temp_calib_values <- long_values %>% dplyr::filter(.data$calibrant == calib)
    for(meas in unique(temp_calib_values$measure)){
      temp_meas_calib_values <- temp_calib_values %>%
        dplyr::filter(.data$measure == meas)

      if(nrow(temp_meas_calib_values) < 3){
        next
      }

      model <- 0

      error_func <- function(x){
        data <- temp_meas_calib_values

        cf <- x[1]
        beta <- x[2]
        error <- 0

        for(i in data$dilution_idx){
          data_i <- data[data$dilution_idx == i,]

          b_i <- data_i$max_concentration * (1 - data_i$dilution_ratio - beta) *
            (data_i$dilution_ratio + beta) ^ (data_i$dilution_idx - 1)

          e_i <- abs(log10(cf * b_i / data_i$normalised_value))^2

          error <- error + e_i
        }

        return(error)
      }

      if(calib == "microspheres"){ # n.b. initial cf value for microspheres needs to be much lower than for fluorescein to acheive a fit
        res <- stats::optim(c(1e-8,0), error_func)
      } else if(calib == "fluorescein"){
        res <- stats::optim(c(1,0), error_func)
      }

      if(res$convergence == 0){
        new_fit <- data.frame(cf = res$par[1], beta = res$par[2],
                              calibrant = calib,
                              fluorophore = temp_meas_calib_values$fluorophore[1],
                              measure = meas)
        fit_values <- rbind(fit_values, new_fit)
      }
    }
  }

  long_values <- dplyr::full_join(long_values, fit_values)


  # plot the mean normalized values -----------------------------------------

  abs_plt <-
    ggplot2::ggplot(data = long_values %>%
                      dplyr::filter(.data$calibrant == "microspheres")) +
    ggplot2::geom_point(ggplot2::aes(x = dilution_idx,
                                     y = normalised_value)) +
    ggplot2::geom_line(ggplot2::aes(x = dilution_idx,
                                    y = cf * max_concentration *
                                      (1 - dilution_ratio - beta) *
                                      (dilution_ratio + beta) ^ (dilution_idx - 1))) +
    ggplot2::scale_y_continuous("Normalised Absorbance", trans = "log10") +
    ggplot2::scale_x_continuous("Microspheres dilution") +
    ggplot2::facet_wrap(~measure) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(gsub(".csv", "_absorbance_cfs.pdf", calibration_csv),
                  plot = abs_plt)

  flu_plt <-
    ggplot2::ggplot(data = long_values %>%
                      dplyr::filter(.data$calibrant == "fluorescein")) +
    ggplot2::geom_point(ggplot2::aes(x = dilution_idx,
                                     y = normalised_value)) +
    ggplot2::geom_line(ggplot2::aes(x = dilution_idx,
                                    y = cf * max_concentration *
                                      (1 - dilution_ratio - beta) *
                                      (dilution_ratio + beta) ^ (dilution_idx - 1))) +
    ggplot2::scale_y_continuous("Normalised Fluorescence", trans = "log10") +
    ggplot2::scale_x_continuous("Fluorescein dilution") +
    ggplot2::facet_wrap(~measure) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(gsub(".csv", "_fluorescence_cfs.pdf", calibration_csv),
                  plot = flu_plt)


  # save conversion factors to a csv ----------------------------------------

  utils::write.csv(fit_values, gsub(".csv", "_cfs.csv", calibration_csv), row.names = FALSE)
}
