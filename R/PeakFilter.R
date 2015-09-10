setGeneric("PeakFilter", 
           function(msdata, ...) standardGeneric("PeakFilter"))

#' Peak filtering
#' 
#' Remove peaks fitting some criteria from data set. Criteria are stated by function arguments.
#' If an argument is equal \code{NULL}, corresponding filtering criterion is not applied.\cr
#' Usually each peak in each replicate group considered separately in terms of each criterion. 
#' After that this peak/replicate group pair can be marked as "filterable" according to this criterion.
#' If for at least one selecting criterion all replicate groups are marked as "filterable", 
#' this peak is removed from data set.
#' @param msdata \code{\link{MSdata-class}} object to be filtered
#' @param blanks The vector of blank samples. Either vector of sample numbers or sample names. 
#' (check names by command: \code{sampleData(msdata)}).
#' These sample will be removed from dataset after filtering.
#' @param above.blank Filter peaks with intensities close to blanks. 
#' \code{mean - above.blank*SE} for certain peak's intensities in replicate group have to be more than 
#' \code{mean + above.blank*SE} for this peak's intensities in blank group.
#' (where \code{SE} - standard error of the mean.)\cr
#' Note: if \code{above.blank = 0} then just mean values are compared.
#' 
#' @param min.int Filter peaks by total intensity. Mean peak intensity in replicate group 
#' have to be higher than \code{min.int} value.
#' @param min.nonNAnum.repgroup Filter peaks with too many missing values. \code{min.number.repgroup} 
#' is minimal number of non-NA values in replicate group.
#' @param min.nonNApercent Filter peaks with too many missing values. \code{min.allgroups} 
#' is minimal allowed quotient of non-NA values for each peak.
#' @return \code{\link{MSdata-class}} object without filtered peaks and blank samples
#' @name PeakFilter
#' @export

setMethod("PeakFilter", "MSdata",
          function(msdata, 
                   blanks = NULL, 
                   above.blank = 2,
                   min.int = 1000,
                   min.nonNAnum.repgroup = 3,
				   min.nonNApercent = 0.4){
			  fold.blank<-NULL
              .intMatrix <- intMatrix(msdata)
              se <- function(x) sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))
              repapply <- function(foo) {
                  apply(.intMatrix, 1, function(peak) tapply(peak, reps, foo, peak))
              }
			  
#               rep.apply <- function(foo) {
#                   sapply(levels(reps), function(repGroup)
#                          apply(.intMatrix[ , reps == repGroup], 1, foo))
#               }

			  if (!is.null(blanks)) {
					blank.intMatrix <- .intMatrix[ , blank]
					.intMatrix <- .intMatrix[ , -blank]
					reps  <- msdata@sampleData[-blanks, ]$ReplicationGroup
			  } else {
					reps  <- msdata@sampleData$ReplicationGroup
			  }
			  filt  <- array(dim = c(length(levels(reps)), nrow(.intMatrix), 0))
			  

			  
              if (is.null(blanks) & (!is.null(fold.blank) | !is.null(above.blank))) {
                  warning("Blank samples are not selected, filtering by blanks has not been conducted")
              
			  } else if ((length(setdiff(blanks, colnames(.intMatrix))) != 0) &
                         (length(setdiff(blanks, 1:ncol(.intMatrix)))   != 0)) {
                  warning("Blank samples set does not match sample names or numbers. Filtering by blanks will not be conducted")
              
			  } else {       
                  if (!is.null(above.blank) & !is.null(blanks)) {
                      filt <- abind::abind(along = 3, filt, repapply(function(peakgr, peak) 
                          mean(peakgr, na.rm = TRUE) - above.blank * se(peakgr) >
                          mean(blank.intMatrix, na.rm = TRUE) + above.blank * se(blank.intMatrix)))
                  }
                  # if (!is.null(fold.blank) & !is.null(blanks)) {
                      # filt <- abind::abind(along = 3, filt, repapply(function(repGroup, peak)
                          # TRUE))
                  # }
              }
              
              if (!is.null(min.int)) {
                  filt <- abind::abind(along = 3, filt, repapply(function(peakgr, ...) mean(peakgr, na.rm = T) >= min.int))
              }
              
              if (!is.null(min.nonNAnum.repgroup)) {
                  filt <- abind::abind(along = 3, filt, repapply(function(peakgr, ...) sum(!is.na(peakgr)) >= min.nonNAnum.repgroup))
              }
			  
			  if (!is.null(min.nonNApercent)) {
				  comparison <- apply(.intMatrix, 1, function(peak) sum(!is.na(peak))/length(peak) >= min.nonNApercent)
                  comparison <- matrix(comparison,
                                    nrow = length(levels(reps)), 
                                    ncol = nrow(.intMatrix), 
                                    byrow = TRUE)
                  filt <- abind::abind(along = 3, filt, comparison)
              }
			  
              filtered.peaks <- apply(apply(filt, 2, apply, 1, all, na.rm = FALSE), 2, any, na.rm = FALSE)
			  msdata <- msdata[filtered.peaks, ]
			  
			  if (!is.null(blanks)) msdata <- msdata[ , -blanks]
              
			  msdata@processLog <- c(msdata@processLog,
                                     "\n\n", sum(!filtered.peaks), " are filtered")
              return(msdata)
          })
		  

		  
# the final variable should be less than 5000 for effective computing
FilterVariable <- function(filter){
     int.mat <- as.matrix(dataSet$proc);
     feat.num <- ncol(int.mat);
     nm <- NULL;
    if(filter == "none" && feat.num <= 2000){ # only allow for less than 2000
        remain <- rep(TRUE, feat.num);
        #dataSet$proc <<- as.data.frame(int.mat);
        msg <- "No data filtering was applied";
    }else{
        if (filter == "rsd" ){
            sds <- apply(int.mat, 2, sd, na.rm=T);
            mns <- apply(int.mat, 2, mean, na.rm=T);
            filter.val <- abs(sds/mns);
            nm <- "Relative standard deviation";
        }else if (filter == "nrsd" ){
            mads <- apply(int.mat, 2, mad, na.rm=T);
            meds <- apply(int.mat, 2, median, na.rm=T);
            filter.val <- abs(mads/meds);
            nm <- "Non-paramatric relative standard deviation";
        }else if (filter == "mean"){
            filter.val <- apply(int.mat, 2, mean, na.rm=T);
            nm <- "mean";
        }else if (filter == "sd"){
            filter.val <- apply(int.mat, 2, sd, na.rm=T);
            nm <- "standard deviation";
        }else if (filter == "mad"){
            filter.val <- apply(int.mat, 2, mad, na.rm=T);
            nm <- "Median absolute deviation";
        }else if (filter == "median"){
            filter.val <- apply(int.mat, 2, median, na.rm=T);
            nm <- "median";
        }else{ # iqr
            filter.val <- apply(int.mat, 2, IQR, na.rm=T);
            nm <- "Interquantile Range";
        }

        # get the rank of the
        rk <- rank(-filter.val, ties.method='random');

        var.num <- ncol(int.mat);
        if(var.num < 250){ # reduce 5%
            remain <- rk < var.num*0.95;
        #    dataSet$proc <<- as.data.frame(int.mat[,rk < var.num*0.95]);
            msg <- paste("Reduce 5\\% features (", sum(!(rk < var.num*0.95)), ") based on", nm);
        }else if(ncol(int.mat) < 500){ # reduce 10%
            remain <- rk < var.num*0.9;
         #   dataSet$proc <<- as.data.frame(int.mat[,rk < var.num*0.9]);
            msg <- paste("Reduce 10\\% features (", sum(!(rk < var.num*0.9)), ") based on", nm);
        }else if(ncol(int.mat) < 1000){ # reduce 25%
            remain <- rk < var.num*0.75;
         #   dataSet$proc <<- as.data.frame(int.mat[,rk < var.num*0.75]);
            msg <- paste("Reduce 25\\% features (", sum(!(rk < var.num*0.75)), ") based on", nm);
        }else{ # reduce 40%, if still over 5000, then only use top 5000
            remain <- rk < var.num*0.6;
            msg <- paste("Reduce 40\\% features (", sum(!remain), ") based on", nm);
            if(sum(remain) > 5000){
                remain <-rk < 5000;
                msg <- paste("Reduced to 5000 features based on", nm);
            }
           # dataSet$proc <<- as.data.frame(int.mat[,remain]);
        }
    }

    dataSet$remain <<- remain; 
    dataSet$filter.msg <<- msg;
    print(msg);
}
