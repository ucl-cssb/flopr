# Example: parsing and calibrating Biotek Cytation plate reader data
#
# Run this with your working directory set to this folder
# (examples/plate_reader/biotek_cytation), e.g. via
# setwd("examples/plate_reader/biotek_cytation") from the repo root.

library(flopr)

# Parse the raw Cytation export into a standard long-format .csv, merging in
# the plate layout metadata (calibrant, fluorophore, concentration, etc. -
# see plate_layout.csv). This is a single-timepoint calibration read, not a
# timeseries, so timeseries = FALSE.
flopr::cytation_parse('231031_pr_calibration_membrane.xlsx',
                       layout_csv = 'plate_layout.csv',
                       timeseries = FALSE)

# Fit calibration curves (fluorescein -> GFP, microspheres -> OD) to the
# parsed data and write out the conversion factors ("_cfs.csv"), along with
# .pdf plots of each fitted curve for a sanity check.
flopr::generate_cfs('231031_pr_calibration_membrane_parsed.csv')
