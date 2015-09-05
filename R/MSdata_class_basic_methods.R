#' Class of MS data table with meta-daat of peaks and samples
#'
#' @slot intMatrix 
#' @slot peakData metadata 
#' @slot sampleData
#' @slot processLog
#' 
#' @return \code{\link{MSdata}} class object
#' 
#' @examples
#' None
#'
#' @export

MSdata <- setClass("MSdata",
                     slots = c(intMatrix  = "matrix",
                               peakData   = "data.frame",
                               sampleData = "data.frame",
                               processLog = "character"))

setMethod("show",
          signature = "MSdata",
          definition = function(object) {
              cat("An object of class ", class(object), "\n", sep = "")
              cat(" ", nrow(object@intMatrix), " peaks by ",
                  ncol(object@intMatrix), " samples.\n", sep = "")
              invisible(NULL)
              })

setMethod("[", "MSdata",
            function(x, i, j, drop = "missing") {
                .intMatrix  <- x@intMatrix[i, j]
                .peakData   <- x@peakData[i, ]
                .sampleData <- x@sampleData[j, ]
                .processLog <- x@processLog
                MSdata(intMatrix  = .intMatrix,
                       peakData   = .peakData,
                       sampleData = .sampleData,
                       processLog = .processLog)
                })

setValidity("MSdata", function(object) {
    msg <- NULL
    valid <- TRUE
    if (nrow(object@intMatrix) != nrow(object@peakData)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of data and peak meta-data rows must be identical.")
    }
    if (ncol(object@intMatrix) != nrow(object@sampleData)) {
        valid <- FALSE
        msg <- c(msg,
                 "Number of data rows and sample meta-data columns must be identical.")
    }
    if (!identical(rownames(object@intMatrix), rownames(object@peakData))) {
        valid <- FALSE
        msg <- c(msg,
                 "Data and peak meta-data row names must be identical.")
    }
    if (!identical(colnames(object@intMatrix), rownames(object@sampleData))) {
        valid <- FALSE
        msg <- c(msg,
                 "Data row names and sample meta-data columns names must be identical.")
    }
    if (valid) TRUE else msg
})



##==============================================================================

## Peak meta-data accesors and replacement methods
setGeneric("peakData", function(object, ...) standardGeneric("peakData"))

setMethod("peakData", "MSdata", 
          function(object) object@peakData)

setGeneric("peakData<-",
           function(object, value) standardGeneric("peakData<-"))

setMethod("peakData<-", "MSdata",
          function(object, value) {
              object@peakData <- value
              if (validObject(object))
                  return(object)
          })

## Sample meta-data accesors and replacement methods
setGeneric("sampleData", function(object, ...) standardGeneric("sampleData"))

setMethod("sampleData", "MSdata",
          function(object) object@sampleData)

setGeneric("sampleData<-",
           function(object, value) standardGeneric("sampleData<-"))

setMethod("sampleData<-", "MSdata",
          function(object, value) {
              object@sampleData <- value
              if (validObject(object))
                  return(object)
          })

## Intensity matrix accesors and replacement methods
setGeneric("intMatrix", function(object, ...) standardGeneric("intMatrix"))

setMethod("intMatrix", "MSdata",
          function(object) object@intMatrix)

setGeneric("intMatrix<-",
           function(object, value) standardGeneric("intMatrix<-"))

setMethod("intMatrix<-", "MSdata",
          function(object, value) {
              object@intMatrix <- value
              if (validObject(object))
                  return(object)
          })

##==============================================================================

