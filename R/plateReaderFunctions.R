#' Plate reader normalisation and calibration
#'
#' @param pr_data a long data.frame containing you plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param od_name the column name for the optical density data
#' @param flu_names the column names for the fluorscence data
#' @param af_model model used to fit negative control autofluorescence.
#' For now these include "polynomial", "invers_poly", "exponential", or "loess".
#' @param to_MEFL a Boolean to determine whether to attempt to convert OD and
#' GFP reading to calibrated units
#' @param GFP_gain if to_MEFL=T, the gain value at which GFP was recorded
#' @param lid_type if to_MEFL=T, the sealing mechanism used on the
#' microtitre-plate. See the conversion_factors_csv for available options
#' @param conversion_factors_csv if to_MEFL=T, path of the csv file containing
#' conversion factors from plate reader calibration
#'
#' @return a data.frame with columns for raw plate reader data, normalised data
#' and, if to_MEFL = T, calibrated OD and GFP data
#' @export
#' @importFrom rlang .data
#'
#' @examples
process_plate <- function(pr_data, blank_well = "A1", neg_well = "A2",
                          od_name = "OD", flu_names = c("GFP"),
                          af_model = "polynomial", to_MEFL = F, GFP_gain = 0,
                          lid_type = "nolid", conversion_factors_csv = NA) {
  od_norm_pr_data <- od_norm(pr_data, blank_well, od_name)

  plt_od <- ggplot2::ggplot(od_norm_pr_data) +
    ggplot2::geom_line(ggplot2::aes_string(x = "time", y = od_name,
                                           colour = '"raw"'), size = 0.2) +
    ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$normalised_OD,
                                    colour = "normalised"), size = 0.2) +
    ggplot2::scale_x_continuous("time") +
    ggplot2::scale_colour_discrete("") +
    ggplot2::facet_grid(row~column) +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(filename = "plot-plate-OD.pdf", plot = plt_od, height = 160,
                  width = 240, units = "mm")

  flu_norm_pr_data <- od_norm_pr_data
  for (flu_idx in seq_len(length(flu_names))) {
    flu_norm_pr_data <- flu_norm(flu_norm_pr_data, neg_well, blank_well,
                                 flu_names[flu_idx], af_model)

    plt_flu <- ggplot2::ggplot(flu_norm_pr_data) +
      ggplot2::geom_line(ggplot2::aes_string(x = "time", y = flu_names[flu_idx],
                                             colour = '"raw"'), size = 0.2) +
      ggplot2::geom_line(ggplot2::aes_string(x = "time",
                                             y = paste("normalised_",
                                                       flu_names[flu_idx],
                                                       sep = ""),
                                             colour = '"normalised"'),
                         size = 0.2) +
      ggplot2::scale_x_continuous("time") +
      ggplot2::scale_colour_discrete("") +
      ggplot2::facet_grid(row~column) +
      ggplot2::theme_bw(base_size = 8)
    ggplot2::ggsave(filename = paste("plot-plate-", flu_names[flu_idx], ".pdf",
                                     sep = ""), plot = plt_flu, height = 160,
                    width = 240, units = "mm")
  }

  out_data <- flu_norm_pr_data

  if (to_MEFL) {
    out_data <- calibrate_plate(out_data, lid_type, GFP_gain, od_name,
                                conversion_factors_csv)
  }

  return(out_data)
}

#' Plate reader normalisation and calibration
#'
#' Deprecated. Please use \code{process_plate()}.
#'
#' @param pr_data a long data.frame containing you plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param OD_name the column name for the optical density data
#' @param flu_names the column names for the fluorscence data
#' @param af_model model used to fit negative control autofluorescence.
#' For now these include "polynomial", "invers_poly", "exponential", or "loess".
#' @param to_MEFL a Boolean to determine whether to attempt to convert OD and
#' GFP reading to calibrated units
#' @param GFP_gain if to_MEFL=T, the gain value at which GFP was recorded
#' @param lid_type if to_MEFL=T, the sealing mechanism used on the
#' microtitre-plate. See the conversion_factors_csv for available options
#' @param conversion_factors_csv if to_MEFL=T, path of the csv file containing
#' conversion factors from plate reader calibration
#'
#' @return a data.frame with columns for raw plate reader data, normalised data
#' and, if to_MEFL = T, calibrated OD and GFP data
#' @export
#'
#'
#' @examples
prNorm <- function(pr_data, blank_well = "A1", neg_well = c("A2"),
                   OD_name = "OD", flu_names = c("GFP"),
                   af_model = "polynomial", to_MEFL = F, GFP_gain = 0,
                   lid_type = "nolid", conversion_factors_csv = NA) {

  return(process_plate(pr_data, blank_well, neg_well, OD_name, flu_names, af_model,
                       to_MEFL, GFP_gain, lid_type, conversion_factors_csv))
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
#' For now these include "polynomial", "invers_poly", "exponential", or "loess".
#'
#' @return
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
flu_norm <- function(pr_data, neg_well, blank_well, flu_name, af_model) {
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
  } else if (af_model == "inverse_poly") {
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
  ggplot2::ggsave(filename = paste("plot-norm-curve-", flu_name, ".pdf", sep = ""),
                  plot = plt)

  # normalise fluorescence data ---------------------------------------------

  if (af_model == "polynomial" | af_model == "power" |
      af_model == "exponential" | af_model == "bi_exponential" |
      af_model == "linear_exponential" | af_model == "linear_power" |
      af_model == "loess") {
    pr_data$v1 <- pr_data$v1 - stats::predict(model, pr_data)
  } else if (af_model == "inverse_poly") {
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

#' Convert arbitrary fluorescence units to calibrated units
#'
#' @param pr_data a data.frame of parsed plate reader data
#' @param gfp_gain if to_MEFL=T, the gain value at which GFP was recorded
#' @param od_name the column name for the optical density data
#' @param lid if to_MEFL=T, the sealing mechanism used on the microtitre-plate.
#' See the conversion_factors_csv for available options
#' @param conversion_factors_csv if to_MEFL=T, path of the csv file containing
#' conversion factors from plate reader calibration
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return
calibrate_plate <- function(pr_data, lid, gfp_gain, od_name, conversion_factors_csv) {
  conversion_factors <- utils::read.csv(conversion_factors_csv)

  # Get conversion factor for OD --------------------------------------------

  od_cf <- unlist(conversion_factors %>%
                    dplyr::filter(.data$measure == od_name) %>%
                    dplyr::filter(.data$lid_type == lid) %>%
                    dplyr::select(.data$slope))


  # Get conversion factor for GFP -------------------------------------------

  gfp_cfs <- conversion_factors %>%
    dplyr::filter(.data$Calibrant == "fluorescein") %>%
    dplyr::filter(.data$lid_type == lid)

  gfp_cfs$measure <- as.numeric(gsub("Gain ", "", gfp_cfs$measure))

  if (length(which(gfp_cfs$measure == gfp_gain)) > 0) {
    gfp_cf <- gfp_cfs[which(gfp_cfs$measure == gfp_gain), ]$slope
  } else {
    model <- stats::lm(log10(slope) ~ poly(measure, 2), data = gfp_cfs)

    gfp_cf <- 10 ^ stats::predict(model, data.frame(measure = gfp_gain))

    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = gfp_cfs$measure,
                                      y = 10^stats::predict(model, gfp_cfs))) +
      ggplot2::geom_point(ggplot2::aes(x = gfp_cfs$measure,
                                       y = gfp_cfs$slope)) +
      ggplot2::geom_vline(xintercept = gfp_gain, linetype = 2) +
      ggplot2::geom_hline(yintercept = 10 ^ stats::predict(model, data.frame(measure = gfp_gain)),
                          linetype = 2) +
      ggplot2::geom_point(ggplot2::aes(x = gfp_gain,
                                       y = 10 ^ stats::predict(model, data.frame(measure = gfp_gain))),
                          colour = "red", shape = 1, size = 2) +
      ggplot2::scale_x_continuous("Gain") +
      ggplot2::scale_y_continuous("Conversion factor (MEFL/a.u.)",
                                  trans = "log10") +
      ggplot2::theme_bw()
  }

  pr_data$calibrated_OD <- pr_data$normalised_OD / od_cf

  MFL_per_uM <- 6.02e13 #MEFL/uM fluorescein (value from Jacob Beal's iGEM conversion excel spreadsheet)
  pr_data$calibrated_GFP <- (pr_data$normalised_GFP / gfp_cf) * MFL_per_uM

  return(pr_data)
}

#' Generate Conversion Factors
#'
#' @param calibration_dir path of the directory in which your calibration data
#' is held
#' @param date input date as YYYYMMDD on which calibration was carried out.
#' n.b. this should also be the name of the folder holding the csv files
#'
#' @return saves a csv data frame with the conversion factors and two png plots
#' of the conversion factors (for Absorbance and fluorescence)
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data
generate_cfs <- function(calibration_dir, date) {
  csv_files <- c("film.csv", "nolid.csv")
  plate_layout <- utils::read.csv(paste(calibration_dir,
                                        "calibration_plate_layout.csv",
                                        sep = "/"))

  # parse data --------------------------------------------------------------
  all_values <- c()
  for (csv_file in csv_files) {
    csv_data <- utils::read.table(paste(calibration_dir, date, csv_file,
                                        sep = "/"),
                                  sep = ",", blank.lines.skip = T, header = F,
                                  stringsAsFactors = F) #read the csv file

    start_time_idx <- which(csv_data[, 1] == "Start Time")  # get start and end time ids
    end_idx <- which(csv_data[, 1] == "End Time")
    names_idx <- which(csv_data[, 1] == "Name")
    names_idx <- names_idx[2:length(names_idx)]  # remove the first start time entry which just details plate type

    lid_type <- gsub(pattern = ".csv", replacement = "", x = csv_file)
    for (i in seq_len(length(start_time_idx))) {
      block_name <- csv_data[names_idx[i], 2]  # record name of what is being measured

      block_start <- start_time_idx[i] + 4  # find start and end of measurement block
      block_end_idx <- end_idx[i] - 3

      new_block <- csv_data[(block_start):(block_end_idx), 1:2]  # grab and name the data
      names(new_block)[1] <- "well"
      names(new_block)[2] <- "value"

      joined_block <- dplyr::full_join(plate_layout, new_block)  # join to plate layout csv, add measurement category
      joined_block$measure <- block_name
      joined_block$lid_type <- lid_type

      all_values <- rbind(all_values, joined_block)  # add to all data
    }
  }

  # remove saturating values and calculate mean of 4 replicates -------------

  all_values$value <- as.numeric(all_values$value)
  summ_values <- all_values %>%
    dplyr::filter(!is.na(.data$value)) %>%
    dplyr::group_by(.data$measure, .data$concentration, .data$Calibrant,
                    .data$lid_type) %>%
    dplyr::summarise(mean_value = mean(.data$value))


  # remove unnecessary observations -----------------------------------------

  summ_values <- summ_values %>%
    dplyr::filter(((.data$Calibrant == "microspheres") &
                     ((.data$measure == "OD600") |
                        (.data$measure == "OD700"))) |
                    ((.data$Calibrant == "fluorescein") &
                       ((.data$measure != "OD600") &
                          (.data$measure != "OD700"))))


  # normalise data ----------------------------------------------------------

  norm_values <- summ_values %>%
    dplyr::group_by(.data$measure, .data$Calibrant, .data$lid_type) %>%
    dplyr::mutate(normalised_value = .data$mean_value -
                    .data$mean_value[.data$concentration == 0]) %>%
    dplyr::filter(((.data$Calibrant == "fluorescein") &
                     (.data$concentration != 0)) |             # exclude fluorescein concentration of 0 to avoid dividing by 0
                    ((.data$Calibrant == "microspheres") &
                       (.data$concentration > 8.139678e6)))  # only take the top 8 data points to avoid low concentration variability

  # fit linear models to the normalised data  --------

  fit_values <- norm_values %>%
    dplyr::group_by(.data$lid_type, .data$measure, .data$Calibrant) %>%
    dplyr::do(broom::tidy(stats::lm(normalised_value~concentration,
                                    data = .data))) %>%  # fit linear model (lm) and extract coefficients (broom::tidy)
    dplyr::select(1:5) %>%                                                   # only keep useful columns
    tidyr::spread(.data$term, .data$estimate) %>%                                        # spread the data so we have a column for intercept coeffs and one for slope coeffs
    dplyr::rename("intercept" = .data$`(Intercept)`,
                  "slope" = .data$concentration)


  # plot the mean normalized values -----------------------------------------

  abs_plt <-
    ggplot2::ggplot(data = norm_values %>%
                      dplyr::filter(.data$Calibrant == "microspheres"),
                    ggplot2::aes(x = .data$concentration,
                                 y = .data$normalised_value)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = y~x) + #fit a line
    ggplot2::scale_y_continuous("Normalised Absorbance", trans = "log10") +
    ggplot2::scale_x_continuous("Number of Microspheres", trans = "log10") +
    ggplot2::facet_grid(measure ~ lid_type) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(paste(calibration_dir, date,
                        "Absorbance_conversion_factors.pdf", sep = "/"),
                  plot = abs_plt)

  flu_plt <-
    ggplot2::ggplot(data = norm_values %>%
                      dplyr::filter(.data$Calibrant == "fluorescein"),
                    ggplot2::aes(x = .data$concentration,
                                 y = .data$normalised_value)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = y~x) +
    ggplot2::scale_y_continuous("Normalized Fluorescence (a.u.)",
                                trans = "log10") +
    ggplot2::scale_x_continuous("Fluorescein Concentration (uM)",
                                trans = "log10") +
    ggplot2::facet_grid(measure ~ lid_type) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(paste(calibration_dir, date,
                        "Fluorescence_conversion_factors.pdf", sep = "/"),
                  plot = flu_plt, width = 12, height = 25, units = "cm")


  # save conversion factors to a csv ----------------------------------------

  utils::write.csv(fit_values, paste(calibration_dir, date, "cfs_generated.csv",
                                     sep = "/"), row.names = FALSE)

  return()
}

#' Generate conversion factors
#'
#' @param calibration_dir path of the directory in which your calibration data
#' is held
#' @param date input date as YYYYMMDD on which calibration was carried out.
#' n.b. this should also be the name of the folder holding the csv files
#'
#' @return saves a csv data frame with the conversion factors and two png plots
#' of the conversion factors (for Absorbance and fluorescence)
#'
#' @export
generateCfs <- function(calibration_dir, date) {
  generate_cfs(calibration_dir, date)
}
