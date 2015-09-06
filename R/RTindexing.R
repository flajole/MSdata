setGeneric("RTindexing", 
           function(object, ...) standardGeneric("RTindexing")
		   
#' Retention time indexing
#'
#' Perform retention time correction in \code{\link{MSdata-class}} object.
#' Function takes a given list of internal standards, corresponding peaks are found 
#' in data set and basing on their RT shift the shifts for the rest of peaks
#' are calculated. The resulting data set is realigned.
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
#' @return \code{\link{MSdata-class}} object with recalculated RT values in peak data 
#' @seealso \code{\link{MZindexing}}
#' @export
setMethod("RTindexing", "MSdata",
          function(object, 
                   targets.list, 
                   mz.window   = 0.02, 
                   rt.window   = 30) {
              peak.rt <- object@peakData[, "RT"]
              peak.mz <- object@peakData[, "MZ"]
              targets <- read.table(targets.list, sep = "\t", col.names = c("compound", "mz", "rt"))
              if (all(targets[, "rt"] < 60)) targets[, "rt"] <- targets[, "rt"] * 60
              ntrg <- dim(targets)[1]
              for (i in 1:ntrg) {
                  targeted.peak <- peak.rt[peak.mz > targets[i, "mz"] - mz.window & 
                                           peak.mz < targets[i, "mz"] + mz.window & 
                                           peak.rt > targets[i, "rt"] - rt.window & 
                                           peak.rt < targets[i, "rt"] + rt.window]
                  targeted.peaks <- c(targeted.peaks, targeted.peak)
              }
              coef.a <- (targets[, "rt"] - c(0, targets[-ntrg, "rt"]))/(targeted.peaks - c(0, targeted.peaks[-ntrg]))
              coef.b <- targets[, "rt"] - targeted.peaks * coef.a
              coef.a <- c(coef.a, coef.a[length(coef.a)])
              coef.b <- c(coef.b, coef.b[length(coef.b)])
              c.a <- sapply(peak.rt, function(x) coef.a[sum(x > targeted.peaks) + 1])
              c.b <- sapply(peak.rt, function(x) coef.b[sum(x > targeted.peaks) + 1])
              object@peakData["rt"] <- peak.rt * c.a + c.b
              object@progressLog <- c(object@progressLog,
                                      "\n\nRetention time is indexed by internal standards list ",
                                      targets.list, "\n",
                                      paste(readLines(targets.list, warn = FALSE), collapse="\n"))
              return(object)
          })
		  
## Retention time indexing in xcmsSet data, with possibility to rewrite raw files

# setMethod("RTindexing", "xcmsSet",
          # function(object, 
                   # targets.list, 
                   # model       = "steplinear", 
                   # mz.window   = 0.02, 
                   # rt.window   = 30, 
                   # write.mzXML = FALSE) {
              # xset <- object
              # targets <- read.table(targets.list, sep = "\t", col.names = c("mz", "rt"))
              # npks <- dim(xset@peaks)[1]
              # ntrg <- dim(targets)[1]
              # nsmp <- length(xset@filepaths)
              # if (all(targets[, "rt"] < 60)) targets[, "rt"] <- targets[, "rt"] * 60
              # targeted.peaks <- array(, dim = c(nsmp, dim(xset@peaks)[2], 0), 
                                      # dimnames = dimnames(xset@peaks))
              # for (i in 1:ntrg) {
                  # targeted.peak <- xset@peaks[xset@peaks[, "mz"] > targets[i, "mz"] - mz.window & 
                                              # xset@peaks[, "mz"] < targets[i, "mz"] + mz.window & 
                                              # xset@peaks[, "rt"] > targets[i, "rt"] - rt.window & 
                                              # xset@peaks[, "rt"] < targets[i, "rt"] + rt.window, ]
                  # if (identical(as.integer(targeted.peak[, "sample"]), 1:nsmp)) 
                  # targeted.peaks <- abind(targeted.peaks, targeted.peak, along = 3)
              # }
              # if (model == "steplinear") {
                  # coef.a <- matrix(targets[, "rt"] - c(0, targets[-ntrg, "rt"]), 
                                   # ncol = ntrg, 
                                   # nrow = nsmp, 
                                   # byrow = TRUE)/(targeted.peaks[, "rt", ] - cbind(rep(0, nsmp), 
                                                                                   # targeted.peaks[, "rt", -ntrg]))
                  # coef.b <- matrix(targets[, "rt"], ncol = ntrg, nrow = nsmp, 
                                   # byrow = TRUE) - targeted.peaks[, "rt", ] * coef.a
                  # coef.a <- cbind(coef.a, coef.a[, dim(coef.a)[2]])
                  # coef.b <- cbind(coef.a, coef.a[, dim(coef.a)[2]])
                  # xset.index <- xset
                  # for (j in 1:npks) {
                      # peak <- xset.deshift@peaks[j, ]
                      # smpl <- peak["sample"]
                      # c.a <- coef.a[smpl, sum(peak["rt"] > targeted.peaks[smpl, "rt", ]) + 1]
                      # c.b <- coef.b[smpl, sum(peak["rt"] > targeted.peaks[smpl, "rt", ]) + 1]
                      # xset.index@peaks[j, c("rt", "rtmin", "rtmax")] <- peak[c("rt", "rtmin", "rtmax")] * c.a + c.b
                  # }
                  # for (j in 1:smpl) {
                      # smpl.rt <- xset.index@rt$raw[[j]]
                      # c.a <- sapply(smpl.rt, function(x) coef.a[j, sum(x > targeted.peaks[j, "rt", ]) + 1])
                      # c.b <- sapply(smpl.rt, function(x) coef.b[j, sum(x > targeted.peaks[j, "rt", ]) + 1])
                      # xset.index@rt$corrected[[j]] <- smpl.rt * c.a + c.b
                  # }
                  
              # }
              
              # if (write.mzXML) for (i in 1:nsmp) {
                  # old.file <- file(xset@filepaths[i], open="r")
                  # ind.file <- file(paste0(
                      # substr(xset@filepaths[i], 1, nchar(xset@filepaths[i]) - 6),
                      # "_rtindex_", model, ".mzXML"), open="w")
                  # while (length(body <- readLines(old.file, n = 19))>0) {
                      # if (any(grepl("retentionTime", body))) {
                          # old.rt <- grep(pattern = "retentionTime=", x = body, value = TRUE)
                          # old.rt <- as.numeric(strsplit(old.rt, split = "\\D{2,}")[[1]][2])
                          # #old.rt <- as.numeric(lapply(strsplit(grep(pattern="retentionTime=",x=body,value=TRUE), split="\\D{2,}"), "[", 2))
                          # c.a <- coef.a[i, sum(oldrt > targeted.peaks[i, "rt", ]) + 1]
                          # c.b <- coef.b[i, sum(oldrt > targeted.peaks[i, "rt", ]) + 1]
                          # ind.rt <- old.rt * c.a + c.b
                          # ind.rt <- ind.rt*(sign(ind.rt) + 1)/2
                          # body[grep("retentionTime", body)] <- paste("          retentionTime=\"PT", as.character(format(ind.rt, digits=6)), "S\"", sep="")		  
                      # }
                      # writeLines(body,ind.file)
                  # }
                  # close(old.file);close(ind.file);
              # }
              # return(xset.index)
          # })


