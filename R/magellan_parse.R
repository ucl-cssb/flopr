#' Parser for Spark/Magellan data
#'
#' Parser for data exported from a Tecan Spark plate reader using SparkControl Magellan software.
#'
#' @param data_csv path to csv file from Tecan Spark plate reader
#' @param layout_csv path to csv file containing plate layout information
#' @param timeseries Boolean flag indicating whether the data is a timeseries or
#'   single recording. The Tecan software outputs the two scenarios differently.
#' @param timestart string indicating the timepoint specified in column1 of the
#'   export file corresponding to the first row of data. "0s" by default.
#' @param interval time interval in minutes between readings in a kinetic loop.
#'   Default is 10.
#' @param mode mode has two options: "read_first" and "incubate_first". Setting
#'   "read_first" mode tells the script that the plate reader method started
#'   with reading the relevant channels, followed by incubation for the interval
#'   time, meaning the first timepoint was at 0 min. Setting "incubate_first"
#'   mode tells the script the opposite was true, making the first timepoint
#'   equal to the length of the set interval, eg 10 min.
#' @param metadata numeric value corresponding to the number of types of
#'   metadata requested during creation of the Excel export file in Magellan.
#'   These can include Well positions, Layout, Replicate info, etc.
#' @param custom Boolean flag indicating whether script should deviate from the
#'   default of collecting data from columns 2:97. If TRUE, script looks at
#'   arguments `insert_wells_above`, `insert_wells_below`, `startcol`, `endcol`.
#'   The total number of columns needs to add up to 96 if the layout_csv file
#'   includes rows A1-H12 and they must be in the same order as the layout file
#'   because this script joins positionally, not by recorded well value (wells
#'   aren't exported by default).
#' @param startcol numeric value corresponding to first column of data_csv
#'   corresponding to data
#' @param endcol numeric value corresponding to last column of data_csv
#'   corresponding to data
#' @param insert_wells_above numeric value corresponding to number of empty
#'   entries to insert before data in custom mode. This can be useful if a
#'   subplate was read instead of an entire plate, meaning the number of rows
#'   created by startcol:endcol does not add up to 96. For example, if the data
#'   starts at B1, but the plate layout starts at A1, can set this to 12 to add
#'   12 empty rows above the data which allows correct joining of data and
#'   layout tables.
#' @param insert_wells_below numeric value corresponding to number of empty
#'   entries to insert after data in custom mode.
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples parsed_calib_plate <- magellan_parse(data_csv = "calibrations/20210104_calibration_data.csv", layout_csv = "calibrations/20210104_calibration_plate_layout.csv", timeseries = FALSE)
#' parsed_data <- magellan_parse(data_csv = "data/20210104_data.csv", layout_csv = "data/20210104_data_layout.csv", timeseries = TRUE, timestart = "0s", interval = 30, mode = "read_first")

magellan_parse <- function(data_csv, layout_csv, timeseries = FALSE,
                           timestart = "0s",
                           interval = 10, # minutes.
                           mode = "read_first", # mode can only be "read_first" or "incubate_first"
                           metadata = 0,
                           custom = FALSE, startcol = 2, endcol = 97, insert_wells_above = 0, insert_wells_below = 0
                           ) {

  data <- utils::read.table(data_csv, sep = ",", blank.lines.skip = TRUE,
                            header = FALSE, stringsAsFactors = FALSE)

  plate_layout <- utils::read.csv(layout_csv)

  # TimeSeries TRUE ----------------------------------------------------------

  if(timeseries == TRUE){

    ## Work out where data begins and the number of channels (eg. OD600, GFP) used
    data_start <- which(data[, 1] == timestart)
    num_channels <- data_start-1-(2*metadata)

    ## Work out the number of rows (timepoints) per channel
    j <- data_start # row number of "0s" or equivalent
    while (!is.na(as.numeric(unlist(strsplit(data[j, 1], split = "s", fixed = TRUE)[1])))) {
      j <- j+1
    }
    # j # is one higher than the row number of the last data row
    # j-1 # is the row number of the last data row
    total_rows <- (j-1) - data_start +1 # total number of rows of data is last row-first row, +1 to be inclusive of first row
    timepoints <- total_rows/num_channels

    ## Check that number of timepoints is an integer, flag if not
    if(!isTRUE(all.equal(timepoints, as.integer(timepoints), check.attributes = FALSE))){
      message("Error: number of data rows (timepoints) per channel is not equal.")
      return()
    }

    all_data <- c()

    for (i in seq_len(num_channels)) {
      ## For each channel..

      ## Find channel name
      block_name <- data[i+metadata, 1]

      ## Find rows of the channel
      # Data is in blocks of rows(all timepoints in channel) * cols(wells A1-H12)
      ## starting row = starting timepoint = first row of all data (data_start) + ith channel-1(i-1)*number timepoints (timepoints)
      block_startrow <- (data_start) + (i-1)*(timepoints)
      ## ending row = ending timepoint = starting timepoint + (timepoints-1)
      block_endrow <- block_startrow + (timepoints-1)

      if(custom == FALSE){
        ## Default

        ## Get data
        new_block <- data[block_startrow:block_endrow, 2:97] # 96 columns = 96 wells
        new_block <- as.data.frame(t(new_block)) # 96 rows

        ## Add timepoints as column names

        ## Cannot use column1 times ("0s" etc) as colnames since you end up with slightly different timepoints
        ## for each channel, which cannot then be combined into a single table.
        ## Instead, need to work out timepoints in minutes from the specified interval time, and the total # of timepoints.

        ## Calculate lag time before first timepoint
        if(mode == "read_first"){ lag <- 0 }
        if(mode == "incubate_first"){ lag <- interval }

        ## Last timepoint in minutes:
        last_time <- (interval*timepoints)-interval+lag
        ## "read_first" mode: OD600, GFP etc readings taken first, then cells are incubated for interval time
        ## first timepoint is 0min
        ## "incubate_first" mode: cells are incubated for interval time first, then OD600, GFP etc readings taken
        ## first timepoint will be however long the interval is, eg. 30min

        ## Times to use as colnames:
        times <- seq(from = 0, to = last_time, by = interval)
        # length(times) # should equal timepoints
        names(new_block) <- times

      } else if(custom == TRUE){
        ## Custom column numbering: use cols startcol:endcol (wells) and add wells above/below

        ## Get data
        preblock <- data.frame(matrix(NA, nrow = timepoints, ncol = insert_wells_above))
        ## add specified number of rows before the data
        datablock <- data[block_startrow:block_endrow, startcol:endcol]
        postblock <- data.frame(matrix(NA, nrow = timepoints, ncol = insert_wells_below))
        ## add specified number of rows after the data
        new_block <- cbind(preblock, datablock, postblock) ## nb. timeseries false rbinds, timeseries true cbinds
        names(new_block) <- seq(from=1, to=96, by=1)

        new_block <- as.data.frame(t(new_block)) ## transform

        ## Add timepoints as column names
        if(mode == "read_first"){ lag <- 0 }
        if(mode == "incubate_first"){ lag <- interval }
        last_time <- (interval*timepoints)-interval+lag
        times <- seq(from = 0, to = last_time, by = interval)
        names(new_block) <- times

      }

      ## Bind data with plate layout:
      combined_block <- cbind(plate_layout, new_block)

      ## Add block_name (eg. OD600) as new 'measure' column
      combined_block$measure <- block_name
      ## this combined_block is its own 96-rows (corresponding to wells) * however many columns

      ## Add block to all data
      all_data <- rbind(all_data, combined_block)  # add to all data
      ## this creates blocks of 96-row-tibbles joined below one another

    }

    # rearrange data ----------------------------------------------------------
    well_idx <- which(names(all_data) == "well")
    gathered_data <- tidyr::gather(all_data, key = "time", value = p,
                                   -c(1:all_of(well_idx), ncol(all_data)))
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

  # TimeSeries FALSE ----------------------------------------------------------

  # else
  if (timeseries == FALSE){

    ## Work out where data begins and the number of channels (eg. OD600, GFP) used
    data_start <- which(data[, 1] == timestart)
    num_channels <- data_start-1-(2*metadata)

    all_data <- c()

    for (i in seq_len(num_channels)) {
      ## For each channel..

      ## Find channel name
      block_name <- data[i+metadata, 1]

      if(custom == FALSE){
        ## Default

        ## Get data
        new_block <- data[(data_start-1 +i), 2:97]
        new_block <- as.data.frame(t(new_block)) ## transform
        names(new_block)[1] <- "value"

      } else if(custom == TRUE){
        ## Custom column numbering

        ## Get data
        preblock <- tibble::tibble(value = rep(NA, insert_wells_above)) ## add specified number of rows above the data
        datablock <- data[(data_start-1 +i), startcol:endcol]
        datablock <- as.data.frame(t(datablock)) ## transform
        names(datablock)[1] <- "value"
        postblock <- tibble::tibble(value = rep(NA, insert_wells_below)) ## add specified number of rows below the data

        new_block <- rbind(preblock, datablock, postblock) ## rowbind these to make 96 rows, one for each well of the data_layout

      }

      new_block$value <- as.numeric(new_block$value) ## required for fluorescence calibrations that include "Overflow" wells

      ## Bind data with plate layout:
      combined_block <- cbind(plate_layout, new_block)

      ## Add block_name (eg OD600) as new 'measure' column
      combined_block$measure <- block_name
      ## this combined_block is its own 96-rows (corresponding to wells) * few columns

      ## Add block to all data
      all_data <- rbind(all_data, combined_block)
      ## this creates blocks of 96-row-tibbles joined below one another
    }

    # rearrange data ----------------------------------------------------------
    spread_data <- tidyr::pivot_wider(all_data, names_from = .data$measure,
                                      values_from = .data$value)
    spread_data
    spread_data$row <- substr(x = spread_data$well, start = 1, stop = 1)
    spread_data$column <- as.numeric(substr(x = spread_data$well, start = 2,
                                            stop = nchar(spread_data$well)))
    spread_data <- dplyr::arrange_at(spread_data, dplyr::vars(.data$row,
                                                              .data$column))

    spread_data

    # write parsed data to csv ------------------------------------------------
    out_name <- gsub(".csv", "_parsed.csv", data_csv)
    utils::write.csv(x = spread_data, file = out_name, row.names = FALSE)

    return(spread_data)
  }

}
