setGeneric("PeakFilter", 
           function(msdata, ...) standardGeneric("PeakFilter"))

setMethod("PeakFilter", "MSdata",
          function(msdata, 
                   blanks = NULL, 
                   fold.blank = NULL,
                   above.blank = 2,
                   min.int = 1000,
                   min.number.repgroup = 3,
                   ...){
              require(abind)
              peaks <- msdata@intMatrix
              reps  <- msdata@sampleData$ReplicationGroup
              filt  <- array(dim = c(length(levels(reps)), dim(peaks)[1], 0))
              
              se <- function(x) sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))
              
              repapply <- function(foo) {
                  apply(peaks, 1, function(peak) tapply(peak, reps, foo, peak))
              }
              
#               rep.apply <- function(foo) {
#                   sapply(levels(reps), function(repGroup)
#                          apply(peaks[ , reps == repGroup], 1, foo))
#               }
              
              if (is.null(blank) & (!is.null(fold.blank) | !is.null(above.blank))) {
                  warning("Blank samples are not selected, filtering by blanks has not been conducted")
              } else if (length(setdiff(blanks, colnames(peaks)) != 0) &
                        (length(setdiff(blanks, 1:ncol(peaks))   != 0) {
                  warning("Blank samples set does not match sample names or numbers, filtering by blanks has not been conducted")
              } else {       
                  if (!is.null(above.blank) & !is.null(blanks)) {
                      filt <- abind(along = 3, filt, repapply(function(peakgr, peak) 
                          mean(peakgr, na.rm = TRUE) - above.blank * se(peakgr) >
                          mean(peak[blanks], na.rm = TRUE) + above.blank * se(peak[blanks])))
                  }
                  if (!is.null(fold.blank) & !is.null(blanks)) {
                      filt <- abind(along = 3, filt, repapply(function(repGroup, peak)
                          TRUE))
                  }
              }
              
              if (!is.null(min.int)) {
                  filt <- abind(along = 3, filt, repapply(function(peakgr, ...) mean(peakgr, na.rm = T) >= min.int))
              }
              
              if (!is.null(min.number.repgroup)) {
                  filt <- abind(along = 3, filt, repapply(function(peakgr, ...) sum(!is.na(peakgr)) >= min.number.repgroup))
              }
              filtered.peaks <- apply(apply(filt, 2, apply, 1, all, na.rm = FALSE), 2, any, na.rm = FALSE)
              msdata <- msdata[filtered.peaks, ]
              msdata@progressLog <- c(msdata@progressLog,
                                      "\n\n", sum(filtered.peaks), " are filtered")
              return(msdata)
          })