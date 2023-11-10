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
    stop("We can currently only parse endpoint from Biotek Cytation plate
         readers.")
  }
  else if (timeseries == FALSE){
    # get start and end block idxs
    start_block_idx <- which(data[, 2] == "Well")
    end_block_idx <- next_blank(start_idx = start_block_idx, data = data)

    # grab the data
    all_data <- data[start_block_idx:end_block_idx, 2:ncol(data)]

    # simplify names
    all_data <- all_data %>%
      rowwise() %>%
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
  }
}
