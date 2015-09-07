setGeneric("EvalMissVal", 
           function(msdata, ...) standardGeneric("EvalMissVal"))
		   
		   
#' Evaluate Missing Values
#' 
#' Fill NA positions in data set.
#' @param msdata \code{\link{MSdata-class}}
#' @param method Method of evaluation, one of:
#' "Min", "FeatureMin", "FeatureMean", "FeatureMedian", "PCA", "GausSim"
#' @param ... Other optional parameters: 
#' \itemize{
#' \item  \code{pcamethods} - character vector containing the name of applied PCA method, one of:
#' 
#' }
#' @include MSdata_class.R
#' @name EvalMissVal
#' @export
setMethod("EvalMissVal", "MSdata",
          function(msdata, 
				   method = "Min", 
				   ...) {
              .intMatrix <- msdata@intMatrix
              npks <- nrows(.intMatrix)
              match.arg(method, c("Min", "FeatureMin", "FeatureMean", "FeatureMedian", "PCA", "GausSim"))
              
              if (method == "GausSim") {
                  randomset <- as.integer(rnorm(n = sum(is.na(.intMatrix)), mean = Mean, sd = Sd))
                  .intMatrix[is.na(.intMatrix)] <- randomset
                  msg <- "Missing variables were replaced with a Gaussian simulated value."
                  
              } else if (method == "Min") {
                  .intMatrix[is.na(.intMatrix)] <- percent * min(.intMatrix)
                  msg <- "Missing variables were replaced with a small value."
                  
              } else if (method == "FeatureMin") {
                  for (i in 1:npks) {
                      pkmin <- min(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the half of minimum values for each feature column."
                  
              } else if (method == "FeatureMean") {
                  for (i in 1:npks) {
                      pkmean <- mean(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the mean concentration of the corresponding compound."
                  
              } else if (method == "FeatureMedian") {
                  for (i in 1:npks) {
                      pkmedian <- median(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the median concentration of the corresponding compound."
                  
              }	else if (method == "PCA") {      
                  # log2 and mean-centering
                  pc <- log2(.intMatrix)
                  pc <- lapply(pcamethod, function(x) pcaMethods::completeObs(pca(pc, nPcs=3, center = TRUE, method=x)))
                  # back transformation
                  pc <- 2 ^ (rowMeans(sapply(pc, '[', is.na(.intMatrix))))
                  .intMatrix[is.na(.intMatrix)] <- rowMeans(sapply(allpc, '[', is.na(.intMatrix)))
				  msg <- "Missing variables were replaced with PCA-predicted values."
              }
              .processLog <- paste0(processLog(msdata), "\n\n EvalMissVal: ", msg)
              MSdata(intMatrix  = .intMatrix,
                     peakData   = peakData(msdata),
                     sampleData = sampleData(msdata),
                     processLog = .processLog)
          })