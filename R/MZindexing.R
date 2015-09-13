setGeneric("MZindexing", 
           function(object, ...) standardGeneric("MZindexing"))

#' Mass-charge indexing
#'
#' Perform mass-charge correction in \code{\link{MSdata-class}} object.
#' Function takes a given list of internal standards, corresponding peaks are found 
#' in data set and basing on their M/Z shift the shifts for the rest of peaks
#' are calculated. The resulting data set is realigned.\cr\cr
#' It's better to perform MZ indexing after \code{\link{RTindexing}}.
#'
#' @param object \code{\link{MSdata-class}} object. 
#' For correct processing there should be "RT" and "MZ" columns in peak data (\code{peakData(object)})
#' @param targets.list File path to simple table of the internal standards. 
#' Table should have three columns: compoundID, mass-charge and retention time (with or without header).
#' @param mz.window Peak corresponding to standard compound is searched in this range 
#' of MZs around standard MZ value.
#' @param rt.window Peak corresponding to standard compound is searched in this range 
#' of RTs around standard RT value.
#'
#' @return \code{\link{MSdata-class}} object with recalculated MZ values in peak data
#' @seealso \code{\link{RTindexing}}
#' @name MZindexing
#' @export
setMethod("MZindexing", "MSdata",
          function(object, 
                   targets.list, 
                   mz.window   = 0.02, 
                   rt.window   = 3) {
				   
			  .peakData <- peakData(object)
			  if (!("RT" %in% names(.peakData)) | !("MZ" %in% names(.peakData))) 
				  stop("There have to be \"RT\" and \"MZ\" columns in peak data table.")
				  
              peak.rt <- .peakData[, "RT"]
              peak.mz <- .peakData[, "MZ"]
              targets <- read.table(targets.list, col.names = c("chemical", "mz", "rt"))
              if (all(targets[, "rt"] < 60)) targets[, "rt"] <- targets[, "rt"] * 60
              ntrg <- dim(targets)[1]
              for (i in 1:ntrg) {
                  targeted.peak <- peak.mz[peak.mz > targets[i, "mz"] - mz.window &
                                           peak.mz < targets[i, "mz"] + mz.window & 
                                           peak.rt > targets[i, "rt"] - rt.window & 
                                           peak.rt < targets[i, "rt"] + rt.window]
                  targeted.peaks <- c(targeted.peaks, targeted.peak)
              }
              coef.a <- (targets[, "mz"] - c(0, targets[-ntrg, "mz"]))/(targeted.peaks - c(0, targeted.peaks[-ntrg]))
              coef.b <- targets[, "mz"] - targeted.peaks * coef.a
              coef.a <- c(coef.a, coef.a[length(coef.a)])
              coef.b <- c(coef.b, coef.b[length(coef.b)])
              c.a <- sapply(peak.mz, function(x) coef.a[sum(x > targeted.peaks) + 1])
              c.b <- sapply(peak.mz, function(x) coef.b[sum(x > targeted.peaks) + 1])
              .peakData["MZ"] <- peak.mz * c.a + c.b
			  msg <- paste0("Mass-charge is indexed by internal standards list ", targets.list,
						    "\n", paste(readLines(targets.list, warn = FALSE), collapse="\n"))
              
			  peakData(msdata) <- .peakData
			  processLog(msdata) <- msg
			  cat(msg)
			  return(msdata)
          })