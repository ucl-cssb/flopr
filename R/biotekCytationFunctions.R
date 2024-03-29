#' Parser for Biotek Cytation plate reader data
#'
#' @param data_xl path to xls or xlsx file from Biotek Cytation plate reader
#' @param layout_csv path to csv file containing plate layout information
#' @param timeseries Boolean flag indicating whether the data is a timeseries or
#' single recording. Currently only works for timeseries=T.
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#'
cytation_parse <- function(data_csv, layout_csv, timeseries=T) {
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

    # find start of data blocks
    end_kinetic_idx <- which(data[, 1] == "End Kinetic")

    # find where the next block starts
    next_blank_idx <- next_blank_cell(end_kinetic_idx, data, col=1)
    next_block_start_idx <- next_filled_cell(next_blank_idx, data, col=1)

    end_of_file <- F
    all_data <- c()
    while (!end_of_file) {
      # find what is being measured
      block_name <- data[next_block_start_idx, 1]

      # start of data
      next_block_start_idx <- next_filled_cell(next_block_start_idx, data, col=2)

      # find where the end of the current measurement block is
      block_end_idx <- next_blank_row(next_block_start_idx, data)
      if(is.na(block_end_idx)){  # if we're on the last block, there is no blank row at the end, so the last row is the end of the block
        block_end_idx <- nrow(data)
        end_of_file <- T
      }

      # grab the data only for that measurement
      new_block <- data[(next_block_start_idx):(block_end_idx - 1), ]

      # manipulate the data
      wells <- new_block[1, -c(1,3)]
      new_block <- new_block[-1,-c(1,3)]
      names(new_block) <- wells
      names(new_block)[1] <- "time"
      new_block <- new_block %>%
        tidyr::pivot_longer(cols = 2:ncol(new_block),
                            names_to = "well",
                            values_to = "value",
                            values_transform = list(value = as.numeric)) %>%
        dplyr::mutate(time = as.numeric(time),
                      time = time*24*60*60,
                      time = round(time/600)*600) # round to nearest 10 minute

      # add info for each well
      joined_block <- dplyr::full_join(plate_layout, new_block)
      joined_block$measure <- block_name

      #
      all_data <- rbind(all_data, joined_block)

      #
      next_block_start_idx <- next_filled_cell(block_end_idx + 1, data, col=1)

      # if we're at the end
      if(is.na(next_block_start_idx)){
        end_of_file <- T
      }
    }

    # rearrange data ----------------------------------------------------------
    out_data <- all_data %>%
      tidyr::pivot_wider(names_from = .data$measure, values_from = .data$value) %>%  # reshape so we have a column for each measurement type
      dplyr::mutate(row = substr(x = .data$well, start = 1, stop = 1)) %>%  # make a "row" column from the "well" column
      dplyr::mutate(column = as.numeric(substr(x = .data$well, start = 2,  # and make a "column" column
                                               stop = nchar(.data$well)))) %>%
      dplyr::arrange_at(dplyr::vars(.data$time,  # order the rows
                                    .data$row,
                                    .data$column))

    # write parsed data to csv ------------------------------------------------
    if(stringr::str_ends(data_csv, ".xlsx")){
      out_name <- gsub(".xlsx", "_parsed.csv", data_csv)
    } else if(stringr::str_ends(data_csv, ".xls")){
      out_name <- gsub(".xls", "_parsed.csv", data_csv)
    } else if(stringr::str_ends(data_csv, ".csv")){
      out_name <- gsub(".xls", "_parsed.csv", data_csv)
    }
    utils::write.csv(x = out_data, file = out_name, row.names = FALSE)

    return(out_data)
  }
  else if (timeseries == FALSE){
    # get start and end block idxs
    start_block_idx <- which(data[, 2] == "Well")
    end_block_idx <- next_blank_row(start_idx = start_block_idx, data = data)
    if(is.na(end_block_idx)){
      end_block_idx <- nrow(data)
    }

    # grab the data
    all_data <- data[start_block_idx:end_block_idx, 2:ncol(data)]

    # simplify names
    all_data <- all_data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(...2 = unlist(strsplit(.data$...2, split = ':'))[1])  %>% # take first section of name
      dplyr::ungroup()

    # transpose the data
    all_data <- as.data.frame(t(as.matrix(all_data)))

    # 1st row contains column names
    names(all_data) <- all_data[1,]
    all_data <- all_data[-1,]
    all_data <- all_data[, -length(all_data)]

    # convert to numeric
    all_data <- all_data %>%
      dplyr::mutate(dplyr::across(!Well, as.numeric)) %>%
      dplyr::rename(well = Well)

    # join to plate layout csv
    joined_block <- dplyr::full_join(plate_layout, all_data, by = 'well')

    # split row and column from well
    joined_block$row <- substr(x = joined_block$well, start = 1, stop = 1)
    joined_block$column <- as.numeric(substr(x = joined_block$well, start = 2,
                                             stop = nchar(joined_block$well)))
    joined_block <- dplyr::arrange_at(joined_block, dplyr::vars(.data$row,
                                                                .data$column))

    # write parsed data to csv ------------------------------------------------

    if(stringr::str_ends(data_csv, ".xlsx")){
      out_name <- gsub(".xlsx", "_parsed.csv", data_csv)
    } else if(stringr::str_ends(data_csv, ".xls")){
      out_name <- gsub(".xls", "_parsed.csv", data_csv)
    } else if(stringr::str_ends(data_csv, ".csv")){
      out_name <- gsub(".xls", "_parsed.csv", data_csv)
    }
    utils::write.csv(x = joined_block, file = out_name, row.names = FALSE)

    return(joined_block)
  }
}
