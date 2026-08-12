### SHARED HELPERS FOR INSTRUMENT PARSERS
### Consolidated from the Cytation/Infinite/Neo2/Spark parser files, which
### each re-implemented these independently (see code review backlog 010,
### 011, 012 for the history).

#' Split a well identifier into row and column
#'
#' @param df a data.frame with a `well` column (e.g. "A1", "H12")
#'
#' @return df with `row` and `column` columns added
#' @noRd
add_row_column <- function(df) {
  df$row <- substr(df$well, 1, 1)
  df$column <- as.numeric(substr(df$well, 2, nchar(df$well)))
  df
}

#' Derive the "_parsed.csv" output filename for a parser's input file
#'
#' @param path path to the raw input file (.csv, .xls, or .xlsx)
#'
#' @return the corresponding "<input>_parsed.csv" path
#' @noRd
parsed_out_name <- function(path) {
  paste0(tools::file_path_sans_ext(path), "_parsed.csv")
}

#' Find next line containing a given string
#'
#' @param start_idx Row index at which to start search
#' @param string String containing a regular expression to search for. Empty cell might be stored as NA in the dataframe
#' @param data dataframe to search within
#' @param col Optional column to search within, searches whole row by default
#'
#' @return row index
#' @noRd
next_char_row <- function(start_idx, string=NA, data, col=NA) {
  if (is.na(col)) {
    col <- seq_len(ncol(data))
  }
  for (i in start_idx:nrow(data)) {
    if (is.na(string)) {
      if(all(is.na(data[i, col]))) {
        return(i)
      }
    } else {
      if (any(grepl(string, data[i, col], fixed = TRUE))) {
        return(i)
      }
    }
  }
  return(NA)
}

#' Find next blank line
#'
#' @param start_idx row index to start searching from
#' @param data dataframe to search within
#'
#' @return row index of next blank line
#' @noRd
next_blank_row <- function(start_idx, data){
  next_start_idx <- start_idx
  while (any(!is.na(data[next_start_idx, ]))) {
    next_start_idx <- next_start_idx + 1
    if(next_start_idx >= nrow(data)){
      return(NA)
    }
  }
  return(next_start_idx)
}

#' Find next blank cell in a given column
#'
#' @param start_idx row index to start searching from
#' @param data dataframe to search within
#' @param col column to search within
#'
#' @return row index of next blank cell
#' @noRd
next_blank_cell <- function(start_idx, data, col){
  next_start_idx <- start_idx
  while (!is.na(data[next_start_idx, col])) {
    next_start_idx <- next_start_idx + 1
    if(next_start_idx >= nrow(data)){
      return(NA)
    }
  }
  return(next_start_idx)
}

#' Find next non-blank cell in a given column
#'
#' @param start_idx row index to start searching from
#' @param data dataframe to search within
#' @param col column to search within
#'
#' @return row index of next non-blank cell
#' @noRd
next_filled_cell <- function(start_idx, data, col){
  next_start_idx <- start_idx
  while (is.na(data[next_start_idx, col])) {
    next_start_idx <- next_start_idx + 1
    if(next_start_idx >= nrow(data)){
      return(NA)
    }
  }
  return(next_start_idx)
}
