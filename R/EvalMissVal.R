setGeneric("EvalMissVal", 
           function(msdata, ...) standardGeneric("EvalMissVal"))
		   
#' Evaluate Missing Values
#' 
#' Fill NA positions in data set. See more: \code{\link[pcaMethods]{pca}}
#' @param msdata \code{\link{MSdata-class}}
#' @param method Method of evaluation, one of:
#' \itemize{
#' \item \code{"min"}, - certain \code{percent} of minimal intensity throughout all dataset;\cr
#' \item \code{"featureMin"}, \code{"featureMean"}, \code{"featureMedian"} - NAs in feature are replaced by corresponding statistics of this feature intensity values;\cr
#' \item \code{"GausSim"} - random generation of normally distributed values;\cr
#' \item \code{"knn"} - nearest neighbour averaging. \code{\link[impute]{impute.knn}} function is used;\cr
#' \item \code{"bpca"}, \code{"ppca"}, \code{"svdImpute"} - imputation based on PCA analysis. \code{\link[pcaMethods]{pca}} function is used.
#' }
#' @param percent For \code{method = "min"}. Quotient of minimal value used to replace missing values.
#' @include MSdata_class.R
#' @name EvalMissVal
#' @export
setMethod("EvalMissVal", "MSdata",
          function(msdata, 
				   method = "min", 
				   percent = 0.5) {
              .intMatrix <- msdata@intMatrix
              npks <- nrow(.intMatrix)
              match.arg(method, c("min", "featureMin", "featureMean", "featureMedian", 
								  "PCA", "GausSim", "knn", "bpca", "ppca", "svdImpute"))
              
              if (method == "GausSim") {
				  Mean <- mean(.intMatrix, na.rm = T)
				  Sd <- sd(.intMatrix, na.rm = T)
                  randomset <- as.integer(rnorm(n = sum(is.na(.intMatrix)), mean = Mean, sd = Sd))
                  .intMatrix[is.na(.intMatrix)] <- randomset
                  msg <- "Missing variables were replaced with a Gaussian simulated value."
                  
              } else if (method == "min") {
                  .intMatrix[is.na(.intMatrix)] <- percent * min(.intMatrix)
                  msg <- "Missing variables were replaced with a small value."
                  
              } else if (method == "featureMin") {
                  for (i in 1:npks) {
                      pkmin <- min(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmin
                  }
                  msg <- "Missing variables were replaced with the half of minimum values for each feature column."
                  
              } else if (method == "featureMean") {
                  for (i in 1:npks) {
                      pkmean <- mean(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the mean concentration of the corresponding compound."
                  
              } else if (method == "featureMedian") {
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
			  intMatrix(msdata) <- .intMatrix
              processLog(msdata) <- paste0("EvalMissVal:\n", msg)
			  cat(msg)
              return(msdata)
          })