setGeneric("EvalMissVal", 
           function(msdata, ...) standardGeneric("EvalMissVal"))
#'

setMethod("EvalMissVal", "MSdata",
          function(msdata, EvalMethod = "Min", ...){
              .intMatrix <- msdata@intMatrix
              npks <- nrows(.intMatrix)
              match.arg <- (EvalMethod, c("Min", "FeatureMin", "FeatureMean", "FeatureMedian", "PCA", "GausSim"))
              
              if (EvalMethod == "GausSim") {
                  randomset <- as.integer(rnorm(n = sum(is.na(.intMatrix)), mean = Mean, sd = Sd))
                  .intMatrix[is.na(.intMatrix)] <- randomset
                  msg <- "Missing variables were replaced with a Gaussian simulated value."
                  
              } else if (EvalMethod == "Min") {
                  .intMatrix[is.na(.intMatrix)] <- percent * min(.intMatrix)
                  msg <- "Missing variables were replaced with a small value."
                  
              } else if (EvalMethod == "FeatureMin") {
                  for (i in 1:npks) {
                      pkmin <- min(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the half of minimum values for each feature column."
                  
              } else if (EvalMethod == "FeatureMean") {
                  for (i in 1:npks) {
                      pkmean <- mean(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the mean concentration of the corresponding compound."
                  
              } else if (EvalMethod == "FeatureMedian") {
                  for (i in 1:npks) {
                      pkmedian <- median(.intMatrix[i, ], na.rm = T);
                      .intMatrix[i, ][is.na(.intMatrix[i, ])] <- pkmean
                  }
                  msg <- "Missing variables were replaced with the median concentration of the corresponding compound."
                  
              }	else if (EvalMethod == "PCA") {
                  require(pcaMethods)
                  
                  # log2 and mean-centering
                  pc <- log2(.intMatrix)
                  pc <- lapply(pcamethod, function(x) completeObs(pca(pc, nPcs=3, center = TRUE, method=x)))
                  # back transformation
                  pc <- 2 ^ (rowMeans(sapply(pc, '[', is.na(.intMatrix))))
                  .intMatrix[is.na(.intMatrix)] <- rowMeans(sapply(allpc, '[', is.na(.intMatrix)))
              }
              .processLog <- paste0(processLog(msdata), "\n\n EvalMissVal: ", msg)
              MSdata(intMatrix  = .intMatrix,
                     peakData   = peakData(msdata),
                     sampleData = sampleData(msdata),
                     processLog = .processLog)
          })