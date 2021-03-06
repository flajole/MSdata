% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetRepGroup.R
\docType{methods}
\name{SetRepGroup}
\alias{SetRepGroup}
\title{Set replication groups numbers}
\usage{
\S4method{SetRepGroup}{MSdata}(msdata, repFac = NULL, impFac = NULL)
}
\arguments{
\item{msdata}{\code{\link{MSdata-class}} object}

\item{repFac}{If you already have replication group numbers in sample metadata, set \code{repFac} - the name of this factoring variable in table.}

\item{impFac}{The vector of the names of factors which are taken into account 
during automatic replication group labeling. If \code{NULL}, all grouping factors are used.}
}
\value{
\code{\link{MSdata-class}} object with \code{$ReplicationGroup} in sample data table
}
\description{
If there are data about replication groups, corresponding column is renamed 
in a standard way to "ReplicationGroup".\cr
Otherwise function automatically adds a \code{ReplicationGroup} column into \code{sampleData}
table according to group factors combinations.
}

