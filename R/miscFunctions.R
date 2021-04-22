#' Parser for Tecan Spark plate reader data
#'
#' @param data_csv path to csv file from Tecan Spark plate reader
#' @param layout_csv path to csv file containing plate layout information
#' @param timeseries Boolean flag indicating whether the data is a timeseries or
#' single recording. The Tecan software outputs the two scenarios differently.
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#' @importFrom rlang .data
#'
spark_parse <- function(data_csv, layout_csv, timeseries=F) {
  data <- utils::read.table(data_csv, sep = ",", blank.lines.skip = T,
                            header = F, stringsAsFactors = F)

  plate_layout <- utils::read.csv(layout_csv)

  if(timeseries == TRUE){
    start_time_idx <- which(data[, 1] == "Start Time")
    if (length(start_time_idx) > 1) {
      start_time_idx <- start_time_idx[length(start_time_idx)]
    }

    next_block_start_idx <- start_time_idx + 2

    end_of_file <- F
    all_data <- c()
    while (!end_of_file) {
      # find what is being measured
      block_name <- data[next_block_start_idx, 1]

      # check if we are at the end of the file
      if (block_name == "End Time") {
        end_of_file <- T
        break
      }

      # find where the end of the current measurement block is
      block_end_idx <- next_block_start_idx
      while (data[block_end_idx, 1] != "") {
        block_end_idx <- block_end_idx + 1
      }

      # grab the data only for that measurement
      new_block <- data[(next_block_start_idx + 1):(block_end_idx - 1), ]

      # trim unecessary readings i.e. temp and cycle number
      # and rename columns
      times <- new_block[2, ]
      new_block <- new_block[c(-1, -2, -3), ]
      names(new_block) <- times
      names(new_block)[1] <- "well"

      # add info for each well
      joined_block <- dplyr::full_join(plate_layout, new_block)
      joined_block$measure <- block_name

      #
      all_data <- rbind(all_data, joined_block)

      #
      next_block_start_idx <- block_end_idx + 1
    }

    # rearrange data ----------------------------------------------------------
    well_idx <- which(names(all_data) == "well")
    gathered_data <- tidyr::gather(all_data, key = "time", value = p,
                                   -c(1:well_idx, ncol(all_data)))
    gathered_data$time <- as.numeric(gathered_data$time)
    gathered_data$p <- as.numeric(gathered_data$p)
    spread_data <- tidyr::spread(gathered_data, key = .data$measure,
                                 value = .data$p)

    spread_data$row <- substr(x = spread_data$well, start = 1, stop = 1)
    spread_data$column <- as.numeric(substr(x = spread_data$well, start = 2,
                                            stop = nchar(spread_data$well)))
    spread_data <- dplyr::arrange_at(spread_data, dplyr::vars(.data$time,
                                                              .data$row,
                                                              .data$column))



    # write parsed data to csv ------------------------------------------------
    out_name <- gsub(".csv", "_parsed.csv", data_csv)
    utils::write.csv(x = spread_data, file = out_name, row.names = FALSE)

    return(spread_data)
  }
  else if (timeseries == FALSE){
    start_time_idx <- which(data[, 1] == "Start Time")  # get start and end time ids
    end_idx <- which(data[, 1] == "End Time")
    names_idx <- which(data[, 1] == "Name")
    names_idx <- names_idx[2:length(names_idx)]  # remove the first start time entry which just details plate type

    all_data <- c()
    for (i in seq_len(length(start_time_idx))) {
      block_name <- data[names_idx[i], 2]  # record name of what is being measured

      block_start <- start_time_idx[i] + 4  # find start and end of measurement block
      block_end_idx <- end_idx[i] - 3

      new_block <- data[(block_start):(block_end_idx), 1:2]  # grab and name the data
      names(new_block)[1] <- "well"
      names(new_block)[2] <- "value"

      new_block$value <- as.numeric(new_block$value)

      joined_block <- dplyr::full_join(plate_layout, new_block)  # join to plate layout csv, add measurement category
      joined_block$measure <- block_name

      all_data <- rbind(all_data, joined_block)  # add to all data
    }

    # rearrange data ----------------------------------------------------------
    spread_data <- tidyr::pivot_wider(all_data, names_from = .data$measure,
                                      values_from = .data$value)
    spread_data$row <- substr(x = spread_data$well, start = 1, stop = 1)
    spread_data$column <- as.numeric(substr(x = spread_data$well, start = 2,
                                            stop = nchar(spread_data$well)))
    spread_data <- dplyr::arrange_at(spread_data, dplyr::vars(.data$row,
                                                              .data$column))

    # write parsed data to csv ------------------------------------------------
    out_name <- gsub(".csv", "_parsed.csv", data_csv)
    utils::write.csv(x = spread_data, file = out_name, row.names = FALSE)

    return(spread_data)
  }
}

#' Parser for Tecan Spark plate reader data
#'
#' @param data_csv path to csv file from Tecan Spark plate reader
#' @param layout_csv path to csv file containing plate layout information
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#'
sparkParse <- function(data_csv, layout_csv) {
  return(spark_parse(data_csv, layout_csv))
}

#' Parser for Biotek Neo 2 plate reader data
#'
#' @param data_xl path to xls or xlsx file from Biotek Neo 2 plate reader
#' @param layout_csv path to csv file containing plate layout information
#' @param timeseries Boolean flag indicating whether the data is a timeseries or
#' single recording. Currently only works for timeseries=T.
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#'
biotek_parse <- function(data_xl, layout_csv, timeseries=T) {
  data <- readxl::read_excel(path = data_xl,
                             col_names = F,
                             trim_ws = T)

  plate_layout <- utils::read.csv(layout_csv)

  if(timeseries == TRUE){
    start_time_idx <- which(data[, 1] == "End Kinetic")

    next_block_start_idx <- start_time_idx + 2

    end_of_file <- F
    all_data <- c()
    measure_times <- c()
    while (!end_of_file) {
      # find what is being measured
      block_name <- unlist(data[next_block_start_idx, 1])

      # check if we are at the end of the file
      if (block_name == "Results") {
        end_of_file <- T
        break
      }

      # find where the end of the current measurement block is
      block_start_idx <- next_block_start_idx + 2
      block_end_idx <- block_start_idx + 1
      while (!is.na(data[block_end_idx, 2])) {
        block_end_idx <- block_end_idx + 1
      }

      # grab the data only for that measurement
      new_block <- data[(block_start_idx):(block_end_idx - 1), -1]

      # trim unnecessary readings i.e. temp
      # and rename columns
      wells <- new_block[1, -2]
      trimmed_new_block <- new_block[-1, -2]
      names(trimmed_new_block) <- wells
      names(trimmed_new_block)[1] <- "time"

      # hack to unify times of readings across different measurements
      if(length(measure_times) == 0){
        measure_times <- trimmed_new_block[1]
      }
      trimmed_new_block[1] <- measure_times

      long_trimmed_new_block <- tidyr::pivot_longer(data = trimmed_new_block,
                                                    cols = 2:ncol(trimmed_new_block),
                                                    names_to = "well")

      # add metadata for each well
      joined_block <- dplyr::full_join(plate_layout, long_trimmed_new_block, by="well")
      joined_block$measure <- block_name

      #
      all_data <- rbind(all_data, joined_block)

      #
      next_block_start_idx <- block_end_idx + 1
    }

    all_data$value <- as.numeric(all_data$value)
    all_data$time <- as.numeric(all_data$time) * 86400  # conversion to seconds - coefficient calculated empirically

    # rearrange data ----------------------------------------------------------
    # includes a bit of a hack to get around duplicate rows (there seems to be
    # lots of empty rows with time = 0)
    out_data <- all_data %>%
      dplyr::group_by(.data$measure) %>%
      dplyr::mutate(id = dplyr::row_number()) %>%
      tidyr::pivot_wider(names_from = .data$measure,
                         values_from = .data$value) %>%
      dplyr::select(-.data$id)

    # get split row and column names for convenience
    out_data$row <- substr(x = out_data$well, start = 1, stop = 1)
    out_data$column <- as.numeric(substr(x = out_data$well, start = 2,
                                         stop = nchar(out_data$well)))

    # write parsed data to csv ------------------------------------------------
    out_name <- paste(tools::file_path_sans_ext(data_xl), "_parsed.csv", sep="")
    utils::write.csv(x = out_data, file = out_name, row.names = FALSE)

    return(out_data)
  }
  else if (timeseries == FALSE){
    stop("We can currently only parse timeseries from Biotek Neo 2 plate
         readers. Please send me some single timepoint data so I can update the
         parser.")
  }
}
