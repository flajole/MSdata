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
#' 
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
              
              if (is.na(match(facA, names(msdata@sampleData))))
                  stop("Factor ", facA, " is not found. Please, check sample data.");
              dataSet$facA.lbl <- facA;
              dataSet$facA <- as.factor(msdata@sampleData[[facA]]);
              cls.lbl <- dataSet$facA;
              dataSet$cls <- as.factor(cls.lbl)
              
              if (!is.null(facB)) {
                  dataSet$format <- "colts";
                  if (is.na(match(facB, names(msdata@sampleData))))
                      stop("Factor ", facB, " is not found. Please, check sample data.");
                  dataSet$facB.lbl <- facB;
                  dataSet$facB <- as.factor(msdata@sampleData[[facB]]);
              }
              
              if(dataSet$design.type == "time"){
                  dataSet$format <- "colts";
                  # determine time factor
                  if(tolower(facA) == "time"){
                      dataSet$time.lbl <- "facA";
                  }else if(tolower(facB) == "time"){
                      dataSet$time.lbl <- "facB";
                  }else{
                      Warning("No time points found in your data");
                      Warning("The time points group must be labeled as \"Time\"");
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
              return(dataSet);
          })