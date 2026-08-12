### DECLARE NAMES THAT ARE NOT REAL GLOBAL VARIABLES
### These are all either bare tidyselect column references (dplyr::arrange(),
### tidyr::pivot_wider(names_from=, values_from=), dplyr::select(),
### dplyr::group_by() - see code review backlog 036/049, which moved these off
### the .data$x pattern since it's separately deprecated in tidyselect
### contexts), ggplot2::aes() column names that require backticks because they
### contain a hyphen (FSC-H, SSC-H, SSC-A), or a ggplot2::after_stat() binding
### (count). R CMD check's static analysis can't tell these apart from real
### undefined globals.
utils::globalVariables(c(
  "measure", "value", "id", "time", "column", "Well",
  "FSC-H", "SSC-H", "SSC-A", "count"
))
