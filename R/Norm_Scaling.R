setGeneric("DataScaling", 
           function(msdata, ...) standardGeneric("DataScaling"))

#' Data scaling
#'
#' This function realizes two steps: mean-centring and scaling.\cr
#' Mean-centring means that for each feature (peak/compound) all samples intensities are considered 
#' as differences from the mean intensity of this feature.\cr
#' Scaling means that these differences are standardized for all the features by dividing by 
#' standard deviation, or range, or another measure of variance.\cr
#' Therefore, after scaling all the features have mean intensity 0 and corresponding variance measure 1.
#'
#' @param msdata \code{\link{MSdata-class}} object
#' @param method The method of scaling, one of:\cr
#' \code{"auto"} - Autoscaling, mean-centring and dividing by the standard deviation of each variable;\cr
#' \code{"pareto"} - Pareto scaling, mean-centring and dividing by the square root of standard deviation of each variable; \cr
#' \code{"range"} - Range scaling, mean-centring and dividing by the range of each variable 
#' @return \code{\link{MSdata-class}} object with normalised intensity matrix
#' @name DataScaling
#' @export 

setMethod("DataScaling", "MSdata",
          function(msdata, method = "pareto") {
              match.arg(method, c("pareto", "auto", "range"))
              .intMatrix <- intMatrix(msdata)
              if (method == "pareto") {
                  .intMatrix <- apply(.intMatrix, 1, ParetoNorm)
                  
              } else if (method == "auto") {
                  .intMatrix <- apply(.intMatrix, 1, AutoNorm)
                  
              } else if (method == "range") {
                  .intMatrix <- apply(.intMatrix, 1, RangeNorm)
              }
			  .intMatrix <- t(.intMatrix)
              #rownames(.intMatrix) <- rownames(intMatrix(msdata))
              #colnames(.intMatrix) <- colnames(intMatrix(msdata))
              msg <- paste0("Data are rescaled by ", method, " scaling")
			  intMatrix(msdata) <- .intMatrix
              processLog(msdata) <- paste0("Normalisation:\n", msg)
			  cat(msg)
              return(msdata)
          })



# normalize to zero mean and unit variance
AutoNorm <- function(x){
    (x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but varaince/SE
ParetoNorm <- function(x){
    (x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but varaince/SE
RangeNorm <- function(x){
    if(max(x) == min(x)){
        x;
    }else{
        (x - mean(x))/(max(x)-min(x));
    }
}