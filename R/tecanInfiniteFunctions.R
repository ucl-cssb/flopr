#' Parser for Tecan Infinite plate reader data
#'
#' @param data_csv path to .csv, .xls or .xlsx file from Tecan Infinit plate reader
#' @param layout_csv path to csv file containing plate layout information
#' @param timeseries Boolean flag indicating whether the data is a timeseries or
#' single recording. The Tecan software outputs the two scenarios differently.
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#' @importFrom rlang .data
#'
infinite_parse <- function(data_csv, layout_csv, timeseries=F) {

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
    start_time_idx <- which(data[, 1] == "Start Time:")
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
      if (block_name == "End Time:") {
        end_of_file <- T
        break
      }

      # find where the end of the current measurement block is
      block_end_idx <- next_blank_cell(next_block_start_idx, data, col=1)

      # grab the data only for that measurement
      new_block <- data[(next_block_start_idx + 1):(block_end_idx - 1), ]
      # new_block <- new_block[-c(1,3), ]  # remove cycle no. and temp.

      # trim unnecessary readings i.e. temp and cycle number
      # and rename columns
      times <- as.character(new_block[2,])
      new_block <- new_block[-c(1:3), ]
      names(new_block) <- times
      names(new_block)[1] <- "well"
      new_block <- new_block[,names(new_block) != "NA"]  # remove erroneous NA columns

      # reshape time into long data
      new_block <- new_block %>%
        tidyr::pivot_longer(cols = 2:ncol(new_block),
                            names_to = "time", names_transform = as.numeric,
                            values_to = "value", values_transform = as.numeric)

      # add info for each well
      joined_block <- dplyr::full_join(plate_layout, new_block)
      joined_block$measure <- block_name

      #
      all_data <- rbind(all_data, joined_block)

      #
      next_block_start_idx <- next_filled_cell(block_end_idx + 1, data, col=1)
    }

    # rearrange data ----------------------------------------------------------
    layout_cols <- ncol(plate_layout)
    out_data <- all_data %>%
      tidyr::pivot_wider(names_from = .data$measure, values_from = .data$value) %>%  # reshape so we have a column for each measurement type
      add_row_column() %>%  # make "row" and "column" columns from the "well" column
      dplyr::arrange_at(dplyr::vars(.data$time,  # order the rows
                                    .data$row,
                                    .data$column))

    # write parsed data to csv ------------------------------------------------
    out_name <- parsed_out_name(data_csv)
    utils::write.csv(x = out_data, file = out_name, row.names = FALSE)

    return(out_data)
  }
  else if (timeseries == FALSE){
    stop("We can currently only parse timeseries from Tecan Infinite plate
         readers. Please send me some single timepoint data so I can update the
         parser.")
  }
}
