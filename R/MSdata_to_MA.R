setGeneric("MSdata_to_MA", 
           function(msdata, ...) standardGeneric("MSdata_to_MA"))
#' Convert MSdata object to MetaboAnalyst object
#' 
#' Create an object for storing data for processing in MetaboAnalysis (\code{MApckg}).
#'
#' @param msdata An object of \code{\link[MSdata]{MSdata}} class.
#' @param designType \code{"time"} if the data is time-series data; \code{"regular"} otherwise. 
#' In the first case one of the arguments \code{facA} or \code{facB} have to be equal \code{"Time"}
#' @param facA Grouping factor, one of the names of sample data columns in \code{msdata@@sampleData} 
#' @param facB Optional second factor.
#' @name MSdata_to_MA
#' @export
setMethod("MSdata_to_MA", "MSdata",
          function (msdata, 
                    designType = "regular", 
                    facA = names(sampleData(msdata))[1], 
                    facB = NULL){
              match.arg(designType, c("time", "regular"))
              dataSet <- list();
              dataSet$type <- "conc";
              dataSet$design.type <- designType;
              dataSet$cls.type <- "disc"; # default
              dataSet$format <- "colu";
              dataSet$paired <- FALSE
              dataSet$combined.method <- TRUE
              
              facA.match <- charmatch(facA, names(msdata@sampleData))
              if (is.na(facA.match))
                  stop("Factor ", facA, " is not found. Please, check sample data.");
              if (facA.match == 0)
                  stop("Multiple factors matching ", facA, " are found. Please, check sample data.");
              facA <- as.factor(msdata@sampleData[[facA.match]])
              facA <- factor(facA, facA[sort(match(levels(facA), facA))])
              dataSet$facA.lbl <- names(msdata@sampleData)[facA.match];
              dataSet$facA <- facA;
              cls.lbl <- dataSet$facA;
              dataSet$cls <- as.factor(cls.lbl)
              
              dataSet <- SetColor(dataSet, msdata = msdata, fac = dataSet$facA.lbl)
              dataSet <- SetShape(dataSet, msdata = msdata, fac = dataSet$facA.lbl)
              
              if (!is.null(facB)) {
                  dataSet$format <- "colts";
                  facB.match <- charmatch(facB, names(msdata@sampleData))
                  if (is.na(facB.match))
                      stop("Factor ", facB, " is not found. Please, check sample data.");
                  if (facB.match == 0)
                      stop("Multiple factors matching ", facB, " are found. Please, check sample data.");
                  facB <- as.factor(msdata@sampleData[[facB.match]])
                  facB <- factor(facB, facB[sort(match(levels(facB), facB))])
                  dataSet$facB.lbl <- names(msdata@sampleData)[facB.match];
                  dataSet$facB <- facB;
              }
              
              if(dataSet$design.type == "time"){
                  dataSet$format <- "colts";
                  # determine time factor
                  if(tolower(facA) == "time"){
                      dataSet$time.lbl <- "facA";
                  }else if(tolower(facB) == "time"){
                      dataSet$time.lbl <- "facB";
                  }else{
                      warning("No time points found in your data");
                      warning("The time points group must be labeled as \"Time\"");
					  dataSet$design.type <- "regular"
                  }
              }
              
              conc <- msdata@intMatrix        
              var.nms <- rownames(conc)
              smpl.nms <- colnames(conc)
              
              if(sum(is.na(iconv(smpl.nms)))>0)
                  stop("No special letters (i.e. Latin, Greek) are allowed in sample names!");
              if(sum(is.na(iconv(var.nms)))>0)
                  stop("No special letters (i.e. Latin, Greek) are allowed in feature names!");
              if(length(unique(smpl.nms))!=length(smpl.nms))
                  stop("Duplicate sample names are detected!");
              if(length(unique(var.nms))!=length(var.nms))
                  stop("Duplicate feature names are detected!");  
              
              empty.inx <- is.na(smpl.nms) | smpl.nms == "";
              if(sum(empty.inx)>0) {
                  warning("Empty sample names are detected. These samples will be removed");
                  conc <- conc[!empty.inx, ]
              }
              
              empty.inx <- is.na(var.nms) | var.nms == "";
              if(sum(empty.inx)>0) {
                  warning("Empty feature names are detected. These features will be removed");
                  conc <- conc[ , !empty.inx];
                  conc <- conc[!empty.inx, ]
              }
              
              dataSet$cmpd <- var.nms;
              dataSet$norm <- t(conc);
              dataSet$cls <- as.factor(cls.lbl)
              dataSet$cls.num <- length(levels(dataSet$cls))
			  cat("MetaboAnalyst dataSet object is created")
              return(dataSet);
          })