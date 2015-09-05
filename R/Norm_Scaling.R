setGeneric("DataScaling", 
           function(msdata, ...) standardGeneric("DataScaling"))

#' Data scaling
#'
#' @param msdata MSdata-class object
#' @param method The method of scaling, one of:\cr
#' \code{"auto"} - Autoscaling, mean-centering and dividing by the standard deviation of each variable;\cr
#' \code{"pareto"} - Pareto scaling, mean-centering and dividing by the square root of standard deviation of each variable; \cr
#' \code{"range"} - Range scaling, mean-centering and dividing by the range of each variable 
#' @return MSdata-class object with normalised intensity matrix
#' @export 

setMethod("DataScaling", "MSdata",
          function(msdata, method = "pareto", )) {
    match.arg(method, c("pareto", "auto", "range"))
    .int.Matrix <- intMatrix(msdata)
    if (method == "pareto") {
        .intMatrix <- apply(.int.Matrix, 1, ParetoNorm)
            
    } else if (method == "auto") {
        .intMatrix <- apply(.int.Matrix, 1, AutoNorm)
            
    } else if (method == "range") {
        .intMatrix <- apply(.int.Matrix, 1, RangeNorm)
    }
    
    
    .processLog <- paste0(.processLog, "Data are rescaled by ", method, " scaling\n")
    MSdata(intMatrix  = .intMatrix,
           peakData   = peakData(msdata),
           sampleData = sampleData(msdata),
           processLog = .processLog)
}



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

# normalize by a sum of each sample, assume constant sum (1000)
# return: normalized data
SumNorm<-function(x){
    1000*x/sum(x, na.rm=T);
}

# normalize by median
MedianNorm<-function(x){
    x/median(x, na.rm=T);
}

# normalize by a reference sample (probability quotient normalization)
# ref should be the name of the reference sample
ProbNorm<-function(x, ref.smpl){
    x/median(as.numeric(x/ref.smpl), na.rm=T)
}

# normalize by a reference reference (i.e. creatinine)
# ref should be the name of the cmpd
CompNorm <- function(x, ref){
    1000*x/x[ref];
}