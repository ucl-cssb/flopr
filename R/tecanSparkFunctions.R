#' Parser for Tecan Spark plate reader data
#'
#' @param data_csv path to .csv, .xls or .xlsx file from Tecan Spark plate reader
#' @param layout_csv path to csv file containing plate layout information
#' @param timeseries Boolean flag indicating whether the data is a timeseries or
#' single recording. The Tecan software outputs the two scenarios differently.
#' @param wells_as_columns Boolean flag indicating whether blocks of data are
#' oriented with wells as columns or rows
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' spark_parse(data_csv = "examples/plate_reader/tecan_spark/191219_calibration_membrane.csv",
#'             layout_csv = "examples/plate_reader/tecan_spark/calibration_plate_layout.csv",
#'             timeseries = FALSE)
#' }
spark_parse <- function(data_csv, layout_csv, timeseries=T, wells_as_columns=F) {

  if(stringr::str_ends(data_csv, ".xlsx") | stringr::str_ends(data_csv, ".xls")){
    data <- as.data.frame(readxl::read_excel(path = data_csv,
                               col_names = F,
                               col_types = "text"))
  } else if(stringr::str_ends(data_csv, ".csv")){
    data <- utils::read.table(data_csv, sep = ",",
                              na.strings = c(""),
                              header = F,
                              stringsAsFactors = F)
  } else {
    stop("data_csv is must be a .csv, .xls or .xlsx file.")
  }

  plate_layout <- utils::read.csv(layout_csv)

  if(timeseries == TRUE){
    start_time_idx <- which(data[, 1] == "Start Time")
    if (length(start_time_idx) > 1) {
      start_time_idx <- start_time_idx[length(start_time_idx)]
    }

    # find where the next block starts
    next_block_start_idx <- next_filled_cell(start_time_idx + 1, data, col=1)

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
      block_end_idx <- next_blank_cell(next_block_start_idx, data, col=1)

      # grab the data only for that measurement
      new_block <- data[(next_block_start_idx + 1):(block_end_idx - 1), ]
      # new_block <- new_block[-c(1,3), ]  # remove cycle no. and temp.

      if(!wells_as_columns){
        # trim unnecessary readings i.e. temp and cycle number
        # and rename columns
        times <- as.character(new_block[2,])
        new_block <- new_block[-c(1:3), ]
        names(new_block) <- times
        names(new_block)[1] <- "well"
        new_block <- new_block %>%
          tidyr::pivot_longer(cols = 2:ncol(new_block),
                              names_to = "time",
                              values_to = "value",
                              values_transform = list(value = as.numeric)) %>%
          dplyr::mutate(time = as.numeric(.data$time))
      } else if(wells_as_columns){
        wells <- new_block[1, -c(1,3)]
        new_block <- new_block[-1,-c(1,3)]
        names(new_block) <- wells
        names(new_block)[1] <- "time"
        new_block <- new_block %>%
          tidyr::pivot_longer(cols = 2:ncol(new_block),
                              names_to = "well",
                              values_to = "value",
                              values_transform = list(value = as.numeric)) %>%
          dplyr::mutate(time = as.numeric(.data$time))
      }

      # add info for each well
      joined_block <- dplyr::full_join(plate_layout, new_block)
      joined_block$measure <- block_name

      #
      all_data <- rbind(all_data, joined_block)

      #
      next_block_start_idx <- next_filled_cell(block_end_idx + 1, data, col=1)
    }

    # rearrange data ----------------------------------------------------------
    out_data <- all_data %>%
      tidyr::pivot_wider(names_from = measure, values_from = value) %>%  # reshape so we have a column for each measurement type
      add_row_column() %>%  # make "row" and "column" columns from the "well" column
      dplyr::arrange(time, row, column)  # order the rows

    # write parsed data to csv ------------------------------------------------
    out_name <- parsed_out_name(data_csv)
    utils::write.csv(x = out_data, file = out_name, row.names = FALSE)

    return(out_data)
  }
  else if (timeseries == FALSE){
    # We use 'Name' in the first column to identify the start of a measurement block
    names_idx <- which(data[, 1] == "Name")
    if (length(names_idx) <= 1) {
      # only the plate-type entry (or nothing) was found, no real measurement blocks
      stop("No measurement blocks found in the data.")
    }
    names_idx <- names_idx[-1]  # remove the first start time entry which just details plate type

    all_data <- c()
    for (i in seq_len(length(names_idx))) {
      block_name <- data[names_idx[i], 2]  # record name of what is being measured

      # the next '<>' following 'Name' is used as the start point for the data in the measurement block
      block_start <- next_char_row(start_idx = names_idx[i] + 1, string = "<>", data = data) + 1

      # the next '' following '<>' is used as the end point for the data in the measurement block
      block_end <- next_char_row(start_idx = block_start + 1, string = NA, data = data) - 1

      new_block <- data[(block_start):(block_end), 1:2]  # grab and name the data
      names(new_block)[1] <- "well"
      names(new_block)[2] <- "value"

      new_block$value <- as.numeric(new_block$value)

      joined_block <- dplyr::full_join(plate_layout, new_block)  # join to plate layout csv, add measurement category
      joined_block$measure <- block_name

      all_data <- rbind(all_data, joined_block)  # add to all data
    }

    # rearrange data ----------------------------------------------------------
    spread_data <- tidyr::pivot_wider(all_data, names_from = measure,
                                      values_from = value)
    spread_data <- add_row_column(spread_data)
    spread_data <- dplyr::arrange(spread_data, row, column)

    # write parsed data to csv ------------------------------------------------
    out_name <- parsed_out_name(data_csv)
    utils::write.csv(x = spread_data, file = out_name, row.names = FALSE)

    return(spread_data)
  }
}

#' Parser for Tecan Spark plate reader data
#'
#' @description
#' Deprecated alias for \code{\link{spark_parse}}. Only forwards \code{data_csv}
#' and \code{layout_csv} - use \code{spark_parse()} directly to access its
#' \code{timeseries} and \code{wells_as_columns} arguments.
#'
#' @param data_csv path to csv file from Tecan Spark plate reader
#' @param layout_csv path to csv file containing plate layout information
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#'
#' @examples
#' \dontrun{
#' sparkParse(data_csv = "examples/plate_reader/tecan_spark/191219_calibration_membrane.csv",
#'            layout_csv = "examples/plate_reader/tecan_spark/calibration_plate_layout.csv")
#' }
sparkParse <- function(data_csv, layout_csv) {
  .Deprecated("spark_parse")
  return(spark_parse(data_csv, layout_csv))
}
