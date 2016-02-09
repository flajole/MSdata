setGeneric("msFillNA", 
           function(msdata, ...) standardGeneric("msFillNA"))
		   
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
#' @name msFillNA
#' @seealso \code{\link{msNoiseGen}} 
#' @export
setMethod("msFillNA", "MSdata",
          function(msdata, 
				   method = "min", 
				   percent = 0.5) {
              .intMatrix <- msdata@intMatrix
              npks <- nrow(.intMatrix)
              match.arg(method, c("min", "featureMin", "featureMean", "featureMedian", "GausSim", "knn", "bpca", "ppca", "svdImpute"), several.ok = T)
              
              pcas <- c("ppca", "bpca", "svdImpute")
              
              if (any(method %in% pcas) & !all(pcas %in% method)) {
                  stop("Only ppca, bpca and svdImpute methods can be chosen together")
                  
              } else if (any(method %in% c("ppca", "bpca", "svdImpute"))) {
                  logint <- log2(.intMatrix)
                  pc <- lapply(method,
                               function(meth) t(pcaMethods::completeObs(pcaMethods::pca(t(logint), center = TRUE, method=meth))))
                  #knn <- impute::impute.knn(t(centint))$data
                  imp <- sapply(pc, '[', is.na(.intMatrix))
                  # back transformation and mean
                  imp <- apply(imp^2, 1, tukey)
                  .intMatrix[is.na(.intMatrix)] <- imp
                  
                  msg <- paste0("Missing variables were replaced with PCA-predicted values (", paste(method, collapse = "/"), " method)")
              } else if (method == "GausSim") {
				  Mean <- mean(.intMatrix, na.rm = T)
				  Sd <- sd(.intMatrix, na.rm = T)
                  randomset <- as.integer(rnorm(n = sum(is.na(.intMatrix)), mean = Mean, sd = Sd))
                  .intMatrix[is.na(.intMatrix)] <- randomset
                  msg <- "Missing variables were replaced with a Gaussian simulated value."
                  
              } else if (method == "min") {
                  .intMatrix[is.na(.intMatrix)] <- percent * min(abs(.intMatrix), na.rm = TRUE)
                  msg <- "Missing variables were replaced with a small value."
                  
              } else if (method == "featureMin") {
                  for (i in 1:npks) {
                      pkmin <- min(abs(.intMatrix[i, ]), na.rm = T);
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
				 
			  }
			  intMatrix(msdata) <- .intMatrix
              processLog(msdata) <- paste0("EvalMissVal:\n", msg)
			  cat(msg)
              return(msdata)
          })



setGeneric("msNoiseGen", 
           function(msdata, ...) standardGeneric("msNoiseGen"))

#' Noise Generator
#' 
#' Generates background signal in replicate groups completely filled with missing values (NAs).
#' To estimate global noise, the median [M] and the 95% quantile [max]
#' are computed for replication groups with less than a half of measurements (minority rule).
#' The standard error [SE] is estimated as 95% quantile for all combinations with
#' more than a half of measurements (majority rule) revealing a groupwise maximum peak height of
#' less than [max].
#' The resulting global noise is  considered as [M ± SE] (as mean and SE) and used to random
#' generate values from a normal distribution constrained with the global noise parameter.
#' 
#' @param msdata \code{\link{MSdata-class}}
#' @include MSdata_class.R
#' @name msNoiseGen
#' @seealso \code{\link{msFillNA}}
#' @export
setMethod("msNoiseGen", "MSdata",
          function(msdata) {
              .intMatrix <- intMatrix(msdata)
              if (is.null(sampleData(msdata)$ReplicationGroup)) {
                  stop("Set replication group factor first (for details see ?SetRepGroup)")
              }

              reps <- sampleData(msdata)$ReplicationGroup
              
              notna <- repapply(msdata, function(group) sum(!is.na(group)))
              suppressWarnings(maxint <- repapply(msdata, max, na.rm=T))
              repsize <- repapply(msdata, length)
              noise.mean <- median(.intMatrix[notna <= repsize%/%2], na.rm=T)
              noise.amp <- quantile(.intMatrix[notna <= repsize%/%2], 0.95, na.rm=T)
              noise.se <- quantile(.intMatrix[(notna > repsize - repsize%/%2) & (maxint < noise.amp)], 0.95, na.rm=T)
              
              noise <- rnorm(sum(notna == 0), noise.mean, noise.se)
              while(any(noise < 0)) {
                  noise[noise <= 0] <- rnorm(sum(noise <= 0),  noise.mean, noise.se)
              }
              
              .intMatrix[notna == 0] <- noise
              intMatrix(msdata) <- .intMatrix
              return(msdata)
          })