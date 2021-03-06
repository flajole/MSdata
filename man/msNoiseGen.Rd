% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MissingValues.R
\docType{methods}
\name{msNoiseGen}
\alias{msNoiseGen}
\title{Noise Generator}
\usage{
\S4method{msNoiseGen}{MSdata}(msdata)
}
\arguments{
\item{msdata}{\code{\link{MSdata-class}}}
}
\description{
Generates background signal in replicate groups completely filled with missing values (NAs).
To estimate global noise, the median [M] and the 95% quantile [max]
are computed for replication groups with less than a half of measurements (minority rule).
The standard error [SE] is estimated as 95% quantile for all combinations with
more than a half of measurements (majority rule) revealing a groupwise maximum peak height of
less than [max].
The resulting global noise is  considered as [M ± SE] (as mean and SE) and used to random
generate values from a normal distribution constrained with the global noise parameter.
}
\seealso{
\code{\link{msFillNA}}
}

