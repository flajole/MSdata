## Normalization by the list of standards written in file
## as table: name, mz, retention time

setGeneric("StandNorm", 
           function(msdata, ...) standardGeneric("StandNorm"))


#' Normalisation by standards
#'
#' Normalisation by the list of external or internal standards. 
#' The list of the standards is provided as a file with three columns: compound, m/z, retention time. 
#' By these MZ and RT values corresponding peaks in dataset are determined 
#' and their sum intensity is used for normalisation.\cr
#' Afterwards, standards are excluded from feature list.
#' @param msdata \code{\link{MSdata-class}} object
#' @param standards.list The link to the file with the table of standards looking like: compound, m/z, retention time
#' @param meanInt Mean sum intensity of all the standards (usually, got from previous experiments).
#' @param mzwindow Range (in ppm) for searching the peak corresponding to a standard from the list.
#' @param rtwindow Range (in seconds) for searching the peak corresponding to a standard from the list.
#' @param recalculateMean If \code{TRUE} then \code{meanInt} is not used, but is recalculated from these particular data.
#' 
#' @return \code{\link{MSdata-class}} object with normalised intensity matrix
#' @name StandNorm
#' @export 		   		   

setMethod("StandNorm", "MSdata",
          function(msdata, 
				   standards.list,
				   mzwindow = 0.01,
				   rtwindow = 10,
				   meanInt = 2000,
				   recalculateMean = FALSE) {
				   
              .intMatrix <- intMatrix(msdata)
              .peakData <- peakData(msdata)
              peak.rt <- .peakData[, "RT"]
              peak.mz <- .peakData[, "MZ"]
              #load standards list, convert into seconds if necessary
              standards.table <- read.table(standards.list, col.names = c("compound", "mz","rt"))
              if (all(standards.table[,"rt"] < 60)) standards.table[,"rt"] <- standards.table[,"rt"]*60
              
              standards <- c()
              
              #find standards in data matrix
              for (i in 1:dim(standards.table)[1]) {
                  found.standard <- which(peak.mz > standards.table[i, "mz"] - 0.01 & 
                                          peak.mz < standards.table[i, "mz"] + 0.01 & 
                                          peak.rt > standards.table[i, "rt"] - 10   & 
                                          peak.rt < standards.table[i, "rt"] + 10)
                  if (length(found.standard) > 1) warning("Found more than one peak corresponding to standard ", i, " - ", standards.table[i, "chemical"]
                  ) else if (length(found.standard) < 1) warning("Peak corresponding to standard ", i, " - ", standards.table[i, "chemical"], " is not found"
                  ) else standards <- c(standards, found.standard)
              }
              
              
              if (is.null(standards)) {
                  stop ("No reliable standards found, data are left native")
              }
              
              #recalculate mean for normalization
              if (recalculateMean) meanInt <- mean(colSums(msdata[standards,], na.rm=TRUE))
              #recalculate intensities
              standards.sums <- matrix(colSums(.intMatrix[standards,], na.rm=TRUE), 
                                       nrow=nrow(.intMatrix), 
                                       ncol=ncol(.intMatrix), byrow=TRUE)
              .intMatrix <- .intMatrix * meanInt / standards.sums 
              intMatrix(msdata) <- .intMatrix
              #exclude standard from data matrix
              msdata <- msdata[-standards, ]
              
              msdata@processLog <- paste0(msdata@processLog, "/n/nData are normalised by internal standards; the following list of standards is used: /n", standards.list)
              return(msdata)
          })