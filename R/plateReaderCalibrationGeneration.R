### GENERATE CALIBRATION CONVERSION FACTORS FROM RAW CALIBRATION-PLATE DATA
### split out from plateReaderFunctions.R (code review backlog 013)

#' Generate Conversion Factors
#'
#' @param calibration_csv path of a .csv file of your calibration data
#'
#' @return Saves a CSV file with columns: cf (conversion factor), beta (pipetting error parameter), calibrant, fluorophore, and measure. Each row corresponds to a fitted conversion factor for a calibrant/fluorophore/measure combination.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#'
#' @examples
#' \dontrun{
#' generate_cfs(calibration_csv =
#'   "examples/plate_reader/tecan_spark/191219_calibration_membrane_parsed.csv")
#' }
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
      non_sat_values <- rbind(non_sat_values, temp_calib_meas_values %>% dplyr::filter(.data$saturated == 'FALSE'))
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
    ggplot2::ggsave(gsub(".csv", paste("_", safe_filename_fragment(calib), "_cfs.pdf", sep = ""),
                        calibration_csv), plot = plt,
                    width = plot_width, height = plot_height, units = "mm")
  }

  utils::write.csv(fit_values, gsub(".csv", "_cfs.csv", calibration_csv), row.names = FALSE)
  return(fit_values)
}
