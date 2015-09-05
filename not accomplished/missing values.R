setGeneric("MissVal", 
           function(msdata, ...) standardGeneric("MissVal"))

setMethod("MissVal", "MSdata",
          function(msdata, MissReplNum = NULL) {
              if (MissRepl == NULL) {
                  return(which(is.null(msdata@intMatrix)))
              } else if (MissRepl){
                  msdata@sampleData$RepGroup
              }
              
          })

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

ImputeVar<-function(int.mat=dataSet$preproc, method="min"){

    new.mat<-NULL;
    msg<-dataSet$replace.msg;

    if(method=="exclude"){
        good.inx<-apply(is.na(int.mat), 2, sum)==0
        new.mat<-int.mat[,good.inx];
        msg <- c(msg,"Variables with missing values were excluded.");
    }else if(method=="min"){
        minConc<-dataSet$minConc;
        int.mat[int.mat==0 | is.na(int.mat)] <- minConc;
        new.mat <- int.mat;
        msg <- c(msg,"Variables with missing values were replaced with a small value.");
    }else if (method=="colmin"){
            new.mat<-apply(int.mat, 2, function(x){
                if(sum(is.na(x))>0){
                    x[is.na(x)]<-min(x,na.rm=T)/2;
                }
                x;
            });
        msg <- c(msg,"Missing variables were replaced with the half of minimum values for each feature column.");
    }else if (method=="mean"){
            new.mat<-apply(int.mat, 2, function(x){
                if(sum(is.na(x))>0){
                    x[is.na(x)]<-mean(x,na.rm=T);
                }
                x;
            });
        msg <- c(msg,"Missing variables were replaced with mean.");
    }else if (method == "median"){
            new.mat<-apply(int.mat, 2, function(x){
                if(sum(is.na(x))>0){
                    x[is.na(x)]<-median(x,na.rm=T);
                }
                x;
            });
        msg <- c(msg,"Missing variables were replaced with median.");
    }else {
        if(method == "knn"){
            suppressMessages(require(impute));
            #print("loading for KNN...");
            new.mat<-t(impute.knn(t(int.mat))$data);
        }else{
            suppressMessages(require(pcaMethods));
            if(method == "bpca"){
                new.mat<-pca(int.mat, nPcs =5, method="bpca", center=T)@completeObs;
            }else if(method == "ppca"){
                new.mat<-pca(int.mat, nPcs =5, method="ppca", center=T)@completeObs;
            }else if(method == "svdImpute"){
                new.mat<-pca(int.mat, nPcs =5, method="svdImpute", center=T)@completeObs;
            }
        }
        msg <- c(msg, paste("Missing variables were imputated using", toupper(method)));
    }

