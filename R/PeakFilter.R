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
#' @param min.nonNAnum.repgroup Filter peaks with too many missing values. \code{nonNAnum.repgroup} 
#' is minimal number of non-NA values in replicate group.
#' @param min.nonNApercent Filter peaks with too many missing values. \code{min.nonNApercent} 
#' is minimal allowed quotient of non-NA values for each peak.
#' @return \code{\link{MSdata-class}} object without filtered peaks and blank samples
#' @name PeakFilter
#' @export

setMethod("PeakFilter", "MSdata",
          function(msdata, 
                   blanks = NULL, 
                   above.blank = NULL,
                   min.int = NULL,
                   min.nonNAnum.repgroup = NULL,
                   min.nonNApercent = 0.4){
              fold.blank<-NULL #temporary
              msg <- c()
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
              
              
              if (dim(filt)[3] == 0) {
                  msg <- paste0(msg, "\n\nNo data filtering was applied");
              } else {
                  filtered.peaks <- apply(apply(filt, 2, apply, 1, all, na.rm = FALSE), 2, any, na.rm = FALSE)
                  msdata <- msdata[filtered.peaks, ]
                  msg <- paste0(msg, "\n\nFiltering: ", sum(!filtered.peaks), " are filtered")
              }
              
              if (!is.null(blanks)) {
                  msdata <- msdata[ , -blanks]
                  msg <- paste0(msg, "\nBlank samples has been removed from dataset.")
              }
              
              msdata@processLog <- paste0(msdata@processLog, msg)
              return(msdata)
          })




setGeneric("BasicFilter", 
           function(msdata, ...) standardGeneric("BasicFilter"))
		   
#' Basic peak filtering
#'
#' The purpose of the data filtering is to identify and remove variables that are unlikely 
#' to be of use when modeling the data. No phenotype information are used in the filtering process, 
#' so the result can be used with any downstream analysis. 
#' This step is strongly recommended for untargeted metabolomics datasets 
#' (i.e. spectral binning data, peak lists) with large number of variables, 
#' many of them are from baseline noises. Filtering can usually improve the results.\cr
#'
#' Non-informative variables can be characterized in two groups: \itemize{
#' \item variables of very small values - can be detected using mean or median; 
#' \item variables that are near-constant throughout the experiment conditions - 
#' can be detected using different variance measures.
#' }
#'
#' The following empirical rules are applied during data filtering:\itemize{
#' \item less than 250 variables: 5\% will be filtered;\cr
#' \item between 250 - 500 variables: 10\% will be filtered;\cr
#' \item between 500 - 1000 variables: 25\% will be filtered;\cr
#' \item over 1000 variables: 40\% will be filtered.
#' }
#'
#' Please note, that \code{"none"} option is only for less than 2000 features. 
#' Over that, if you choose \code{"none"}, the  IQR filter will still be applied.
#' 
#' The maximum allowed number of variables is 5000. 
#' If over 5000 variables were left after filtering, only the top 5000 will be used in the 
#' subsequent analysis.
#' @param msdata \code{\link{MSdata-class}} object to be filtered
#' @param method Method of filtering one of:\cr
#' \code{"none"} - no filtering applied\cr
#' \code{"iqr"} - interquantile range (IQR)\cr
#' \code{"sd"} - standard deviation (SD)\cr
#' \code{"mad"} - median absolute deviation (MAD)\cr
#' \code{"rsd"} - relative standard deviation (RSD = SD/mean)\cr
#' \code{"nprsd"} - non-parametric relative standard deviation (MAD/median)\cr
#' \code{"mean"} - mean intensity value\cr
#' \code{"median"} - median intensity value
#' none (less than 2000 features)
#' the final variable should be less than 5000 for effective computing
#' @return \code{\link{MSdata-class}} object without filtered peaks
#' @name BasicFilter
#' @export

setMethod("BasicFilter", "MSdata",
          function(msdata,
                   method = "none"){
              
              match.arg(method, c("none", "iqr", "sd", "mad", "rsd", "nprsd", "mean", "median"))
              .intMatrix <- intMatrix(msdata)
              npks <- nrow(.intMatrix);
              
              if (method == "none") {
                  if (npks <= 2000) {
                      msg <- "No data filtering was applied"
                      msdata@processLog <- paste0(msdata@processLog, "\n\n", msg)
                      return(msdata)
                  } else 
                      method <- "iqr"
              }
              
              if (method == "rsd" ){
                  sds <- apply(.intMatrix, 1, sd, na.rm=T);
                  mns <- apply(.intMatrix, 1, mean, na.rm=T);
                  filter.val <- abs(sds/mns);
                  nm <- "Relative standard deviation";
              } else if (method == "nprsd" ){
                  mads <- apply(.intMatrix, 1, mad, na.rm=T);
                  meds <- apply(.intMatrix, 1, median, na.rm=T);
                  filter.val <- abs(mads/meds);
                  nm <- "Non-parametric relative standard deviation";
              } else if (method == "mean"){
                  filter.val <- apply(.intMatrix, 1, mean, na.rm=T);
                  nm <- "mean";
              } else if (method == "sd"){
                  filter.val <- apply(.intMatrix, 1, sd, na.rm=T);
                  nm <- "Standard deviation";
              } else if (method == "mad"){
                  filter.val <- apply(.intMatrix, 1, mad, na.rm=T);
                  nm <- "Median absolute deviation";
              } else if (method == "median"){
                  filter.val <- apply(.intMatrix, 1, median, na.rm=T);
                  nm <- "median";
              } else if (method == "iqr"){
                  filter.val <- apply(.intMatrix, 1, IQR, na.rm=T);
                  nm <- "Interquantile Range";
              }
              
              # get the rank of the
              rk <- rank(-filter.val, ties.method='random');
              
              if(npks < 250){ # reduce 5%
                  remain <- rk < npks*0.95;
                  msg <- paste("Reduce 5% features (", sum(!remain), ") based on", nm);
              } else if(npks < 500){ # reduce 10%
                  remain <- rk < npks*0.9;
                  msg <- paste("Reduce 10% features (", sum(!remain), ") based on", nm);
              } else if(npks < 1000){ # reduce 25%
                  remain <- rk < npks*0.75;
                  msg <- paste("Reduce 25% features (", sum(!remain), ") based on", nm);
              } else{ # reduce 40%, if still over 5000, then only use top 5000
                  remain <- rk < npks*0.6;
                  msg <- paste("Reduce 40% features (", sum(!remain), ") based on", nm);
                  if(sum(remain) > 5000){
                      remain <- rk < 5000;
                      msg <- paste("Reduced to 5000 features based on", nm);
                  }
              }
              
              msdata <- msdata[remain, ]
              msdata@processLog <- paste0(msdata@processLog, "\n\n", msg)
              return(msdata)
          })