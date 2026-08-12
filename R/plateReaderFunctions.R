#' Plate reader normalisation and calibration
#'
#' @param data_csv path to a .csv file containing parsed plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param od_name the column name for the optical density data
#' @param flu_names the column names for the fluorescence data
#' @param af_model model used to fit negative control autofluorescence. One of
#' "polynomial", "inverse_poly", "exponential", "bi_exponential",
#' "linear_exponential", "power", "linear_power", "loess" or "spline".
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
#' \dontrun{
#' process_plate(data_csv = "examples/plate_reader/tecan_spark/200228_example_data_parsed.csv",
#'               blank_well = c("C12", "D12"),
#'               neg_well = c("C6", "D6", "E6"),
#'               od_name = "OD700",
#'               flu_names = c("GFP", "mCherry"),
#'               af_model = "spline",
#'               to_MEFL = TRUE,
#'               flu_gains = 135,
#'               conversion_factors_csv =
#'                 "examples/plate_reader/tecan_spark/191219_calibration_membrane_parsed_cfs.csv")
#' }
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
#' @noRd
od_norm <- function(pr_data, blank_well, od_name) {
  pr_data$normalised_OD <- pr_data[, od_name]

  # remove background absorbance signal -------------------------------------

  pr_data <- pr_data %>%
    dplyr::group_by(.data$time) %>%
    dplyr::mutate(normalised_OD = .data$normalised_OD -
                    mean(.data$normalised_OD[.data$well %in% blank_well]))

  return(as.data.frame(pr_data))
}


#' Autofluorescence model definitions used by \code{flu_norm}
#'
#' Each entry provides a \code{fit} function (negative-control data -> a
#' fitted model) and a \code{predict} function (model, newdata -> predicted
#' fluorescence), so \code{flu_norm} only has one place that dispatches on
#' \code{af_model} instead of three parallel chains.
#' @importFrom mgcv gam predict.gam
#' @noRd
af_models <- list(
  polynomial = list(
    fit = function(data) {
      stats::nls(v1 ~ (a * normalised_OD + b * normalised_OD ^ 2 + c),
                start = c(a = 1, b = 1, c = 1), data = data)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  inverse_poly = list(
    fit = function(data) {
      stats::nls(normalised_OD ~ (a * v1 + b * v1 ^ 2 + c),
                start = c(a = 1, b = 1, c = 1), data = data)
    },
    predict = function(model, newdata) {
      ((- (stats::coef(model)[1]) +
          sqrt((stats::coef(model)[1]) ^ 2 -
                 4 * (stats::coef(model)[2]) * (stats::coef(model)[3]) +
                 4 * (stats::coef(model)[2]) * newdata$normalised_OD)) /
         (2 * (stats::coef(model)[2])))
    }
  ),
  exponential = list(
    ## ae^(bx) + c
    fit = function(data) {
      ## intial parameter estimation
      model_0 <- stats::lm(log(v1) ~ normalised_OD, data = data)
      start <- list(a = exp(stats::coef(model_0)[1]),
                    b = stats::coef(model_0)[2],
                    c = -1)
      stats::nls(v1 ~ (a * exp(b * normalised_OD) + c),
                start = start, data = data)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  bi_exponential = list(
    ## a exp(bx) + c exp(dx) + e
    fit = function(data) {
      model_0 <- stats::lm(log(v1) ~ normalised_OD, data = data)
      start <- list(a = exp(stats::coef(model_0)[1]*0.2),
                    b = stats::coef(model_0)[2]*0.2,
                    c = exp(stats::coef(model_0)[1])*0.8,
                    d = stats::coef(model_0)[2]*0.8,
                    e = 1)
      stats::nls(v1 ~ (a * exp(b * normalised_OD) +
                          c * exp(d * normalised_OD) + e),
                start = start, data = data)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  linear_exponential = list(
    ## ax + be^cx + d
    fit = function(data) {
      model_01 <- stats::lm(v1 ~ normalised_OD, data = data)
      model_02 <- stats::lm(log(v1) ~ normalised_OD, data = data)
      start <- list(a = stats::coef(model_01)[2],
                    b = exp(stats::coef(model_02)[1]),
                    c = stats::coef(model_02)[2],
                    d = 1)
      stats::nls(v1 ~ (a * normalised_OD +
                          b * exp(c * normalised_OD) + d),
                start = start, data = data)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  power = list(
    ## ax^b + c
    fit = function(data) {
      model_0 <- stats::lm(log(v1) ~ log(normalised_OD), data = data)
      start <- list(a = exp(stats::coef(model_0)[1]),
                    b = stats::coef(model_0)[2],
                    c = 1)
      stats::nls(v1 ~ (a * normalised_OD ^ b + c),
                start = start, data = data)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  linear_power = list(
    ## ax + bx^c + d
    fit = function(data) {
      model_01 <- stats::lm(v1 ~ normalised_OD, data = data)
      model_02 <- stats::lm(log(v1) ~ log(normalised_OD), data = data)
      start <- list(a = stats::coef(model_01)[2],
                    b = exp(stats::coef(model_02)[1]),
                    c = stats::coef(model_02)[2],
                    d = 1)
      stats::nls(v1 ~ (a * normalised_OD + b * normalised_OD ^ c + d),
                start = start, data = data)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  loess = list(
    fit = function(data) {
      stats::loess(v1 ~ normalised_OD, data = data, span = 0.5)
    },
    predict = function(model, newdata) stats::predict(model, newdata)
  ),
  spline = list(
    fit = function(data) {
      mgcv::gam(v1 ~ s(normalised_OD), data = data)
    },
    predict = function(model, newdata) mgcv::predict.gam(model, newdata)
  )
)

#' Normalise fluorescence against negative well
#'
#' @param pr_data a long data.frame containing you plate reader data with OD
#' normalised
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param blank_well the well coordinates of a media blank
#' @param flu_name the column name of the fluorescence chanel to normalise
#' @param af_model model used to fit negative control autofluorescence. One of
#' "polynomial", "inverse_poly", "exponential", "bi_exponential",
#' "linear_exponential", "power", "linear_power", "loess" or "spline".
#' @param data_csv path to the original data. Used for saving normalisation curve plots.
#'
#' @return an updated data.frame with an additional \code{normalised_<flu_name>}
#' column
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @noRd
flu_norm <- function(pr_data, neg_well, blank_well, flu_name, af_model, data_csv) {
  pr_data$v1 <- pr_data[, flu_name]

  # fit autofluorescence model to negative control --------------------------

  negative_data <- pr_data %>% dplyr::filter(.data$well %in% neg_well)

  model_spec <- af_models[[af_model]]
  if (is.null(model_spec)) {
    stop("Unknown af_model \"", af_model, "\". Must be one of: ",
         paste(names(af_models), collapse = ", "))
  }
  model <- model_spec$fit(negative_data)

  # plot model fit curve ------------------------------------------------------

  plt <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                    y = model_spec$predict(model, negative_data))) +
    ggplot2::geom_point(ggplot2::aes(x = negative_data$normalised_OD,
                                     y = negative_data$v1)) +
    ggplot2::scale_x_continuous("normalised_OD") +
    ggplot2::scale_y_continuous(flu_name) +
    ggplot2::theme_bw()
  ggplot2::ggsave(filename = gsub(".csv",
                                  paste("_norm-curve_", flu_name, ".pdf", sep = ""),
                                  data_csv),
                  plot = plt)

  # normalise fluorescence data ---------------------------------------------

  pr_data$v1 <- pr_data$v1 - model_spec$predict(model, pr_data)

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
#' @return an updated data.frame with an additional \code{calibrated_OD} column
#' @noRd
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
#' @return an updated data.frame with an additional \code{calibrated_<flu_name>}
#' column
#' @noRd
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
