setGeneric("MZindexing", 
           function(object, ...) standardGeneric("MZindexing"))

setMethod("MZindexing", "MSdata",
          function(object, 
                   targets.list, 
                   model       = "steplinear", 
                   mz.window   = 0.02, 
                   rt.window   = 3) {
              peak.rt <- object@peakData[, "rt"]
              peak.mz <- object@peakData[, "mz"]
              targets <- read.table(targets.list, sep = "\t", col.names = c("chemical", "mz", "rt"))
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
              object@peakData["mz"] <- peak.mz * c.a + c.b
              object@progressLog <- c(object@progressLog,
                                      "\n\nMass-charge is indexed by internal standards list ",
                                      targets.list, "\n",
                                      paste(readLines(targets.list, warn = FALSE), collapse="\n"))
              return(object)
          })