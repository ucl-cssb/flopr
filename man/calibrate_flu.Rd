% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plateReaderFunctions.R
\name{calibrate_flu}
\alias{calibrate_flu}
\title{Convert arbitrary fluorescence units to calibrated units}
\usage{
calibrate_flu(pr_data, flu_name, flu_gain, od_name, conversion_factors_csv)
}
\arguments{
\item{pr_data}{a data.frame of parsed plate reader data}

\item{flu_name}{name of fluorophore to be calibrated}

\item{flu_gain}{gain at which the fluorophore was measured}

\item{od_name}{the column name for the optical density data}

\item{conversion_factors_csv}{if to_MEFL=T, path of the csv file containing
conversion factors from plate reader calibration}
}
\description{
Convert arbitrary fluorescence units to calibrated units
}
