setGeneric("RTrename", 
           function(msdata, ...) standardGeneric("RTrename"))

setMethod("RTrename", "MSdata",
          function(msdata, RT = NULL) {
              .peakData <- peakData(msdata)
              peak.vars <- tolower(names(peakData(msdata)))
              msg <- c()
              if (!is.null(RT)) {
                  col.rt <- names(.peakData) == RT
                  if (any(col.rt)) {
                      names(.peakData)[col.rt] <- "RT"
                      peakData(msdata) <- .peakData
                      msg <- paste0(msg, "Column ", col.rt, " in peak data is renamed to RT.")
                      processLog(msdata) <- msg
					  cat(msg)
                      return(msdata)
                  } else {
                      msg <- c(msg, "InputRT column name was not found./n")
                  }
              } 
              
			  matches <- sapply(c("reten", "rt", "time"), pmatch, peak.vars)
              if (length(matches) == 0) {
                  stop("No RT column is found.")
              } else {
                  col.rt <- na.omit(matches)[1]
                  names(.peakData)[col.rt] <- "RT"
                  msg <- paste0(msg, "Column ", col.rt, " in peak data is renamed to RT.")
                  processLog(msdata) <- msg
				  cat(msg)
                  return(msdata)
              }
          })
