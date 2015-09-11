setGeneric("EvalMissVal", 
           function(msdata, ...) standardGeneric("EvalMissVal"))
		   
#' Evaluate Missing Values
#' 
#' Fill NA positions in data set. See more: \code{\link[pcaMethods]{pca}}
#' @param msdata \code{\link{MSdata-class}}
#' @param method Method of evaluation, one of:
#' \code{"Min"},\cr
#' \code{"FeatureMin"}, \code{"FeatureMean"}, \code{"FeatureMedian"}\cr
#' \code{"GausSim"}\cr
#' \code{"knn"}\cr
#' \code{"bpca"}, \code{"ppca"} or \code{"svdImpute"}
#' 
#' @param percent For \code{method = "Min"}. Quotient of minimal value used to replace missing values.
#' @include MSdata_class.R
#' @name EvalMissVal
#' @export
setMethod("EvalMissVal", "MSdata",
          function(msdata, 
				   method = "Min", 
				   percent = 0.5) {
              .intMatrix <- msdata@intMatrix
              npks <- nrow(.intMatrix)
              match.arg(method, c("Min", "FeatureMin", "FeatureMean", "FeatureMedian", 
								  "PCA", "GausSim", "knn", "bpca", "ppca", "svdImpute"))
              
              if (method == "GausSim") {
				  Mean <- mean(.intMatrix, na.rm = T)
				  Sd <- sd(.intMatrix, na.rm = T)
                  randomset <- as.integer(rnorm(n = sum(is.na(.intMatrix)), mean = Mean, sd = Sd))
                  .intMatrix[is.na(.intMatrix)] <- randomset
                  msg <- "Missing variables were replaced with a Gaussian simulated value."
                  
              } else if (method == "Min") {
                  .intMatrix[is.na(.intMatrix)] <- percent * min(.intMatrix)
                  msg <- "Missing variables were replaced with a small value."
                  
              } else if (method == "FeatureMin") {
                  for (i in 1:npks) {
                      pkmin <- min(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmin
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
                  
              }	else if (method == "knn"){
				 .intMatrix <- impute::impute.knn(.intMatrix)$data
				 msg <- "Missing variables were replaced with KNN-predicted values."
				 
			  }	else if (method %in% c("ppca", "bpca", "svdImpute")) {
				   
				  .intMatrix <- pcaMethods::completeObs(pcaMethods::pca(.intMatrix, nPcs = 5, center = TRUE, method = method))
                  # log2 and mean-centering
                  # pc <- log2(.intMatrix)
                  # pc <- lapply(pcamethod, function(x) pcaMethods::completeObs(pca(pc, nPcs=3, center = TRUE, method=x)))
                  # # back transformation
                  # pc <- 2 ^ (rowMeans(sapply(pc, '[', is.na(.intMatrix))))
                  # .intMatrix[is.na(.intMatrix)] <- rowMeans(sapply(allpc, '[', is.na(.intMatrix)))
				  
				  msg <- paste0("Missing variables were replaced with PCA-predicted values (", method, " method)")
              }
              .processLog <- paste0(msdata@processLog, "\n\nEvalMissVal:\n", msg)
              MSdata(intMatrix  = .intMatrix,
                     peakData   = peakData(msdata),
                     sampleData = sampleData(msdata),
                     processLog = .processLog)
          })