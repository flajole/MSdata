setGeneric("DataTransform", 
           function(msdata, ...) standardGeneric("DataTransform"))
		   
#' Data transformation
#'
#' Logarithmical or cube root data transformation.
#' @param msdata \code{\link{MSdata-class}} object
#' @param method The method of transformation, one of:\cr
#' \code{"log10"} - log10-transformation; \cr
#' \code{"log2"} - log2-transformation;\cr
#' \code{"glog10"} - generalised log10, tolerant to zeros (zeros are replaced with 1/10 of the minimal value);\cr
#' \code{"glog2"} - generalised log2, tolerant to zeros (zeros are replaced with 1/10 of the minimal value); \cr
#' \code{"cuberoot"} - cube root transformation \cr
#' 
#' @return \code{\link{MSdata-class}} object with transformed intensity matrix
#' @name DataTransform
#' @export 		   
setMethod("DataTransform", "MSdata",
          function(msdata, method = "glog10") {
              match.arg(method, c("glog10", "glog2", "log10", "log2", "cuberoot"))
              .intMatrix <- intMatrix(msdata)
              min.val <- min(abs(.intMatrix[.intMatrix!=0]))/10;
              if (method == "log10"){
                  .intMatrix <- apply(.intMatrix, 1, log10)
                  msg <- "Common Logarithm Transformation"
              } else if (method == "log2") {
                  .intMatrix <- apply(.intMatrix, 1, log2)
                  msg <- "Binary Logarithm Transformation"
              } else if (method == "glog10") {
                  .intMatrix <- apply(.intMatrix, 1, glog10, min.val)
                  msg <- "Generalised Common Logarithm Transformation"
              } else if (method == "glog2") {
                  .intMatrix <- apply(.intMatrix, 1, glog2, min.val)
                  msg <- "Generalised Binary Logarithm Transformation"
              } else if (method == "cuberroot") {
                  .intMatrix <- .intMatrix^1/3
                  msg <- "Cube Root Transformation"
              }
			  .intMatrix <- t(.intMatrix)
              .processLog <- paste0(msdata@processLog, "Data are transformed by ", msg, "\n")
              MSdata(intMatrix  = .intMatrix,
                     peakData   = peakData(msdata),
                     sampleData = sampleData(msdata),
                     processLog = .processLog)
          })

glog2 <- function(x, min.val){
    log10((x + sqrt(x^2 + min.val^2))/2)
}

glog10 <- function(x, min.val){
    log2((x + sqrt(x^2 + min.val^2))/2)
}