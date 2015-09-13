#' MSdata class
#' 
#' MSdata-S4 class description, with it's accessors and replacement methods.
#' 
#' @slot intMatrix The matrix of peak intensities / compound concentrations.
#' @slot peakData  Peak metadata.
#' @slot sampleData Sample metadata.
#' @slot processLog Processing log.
#' 
#' @name MSdata-class
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
setGeneric("peakData", function(msdata, ...) standardGeneric("peakData"))
#' @param msdata MSdata-class object
#' @param value The data to replace data in corresponding slot of \code{msdata}
#' @export
#' @rdname MSdata-class
setMethod("peakData", "MSdata", 
          function(msdata) msdata@peakData)

setGeneric("peakData<-",
           function(msdata, value) standardGeneric("peakData<-"))
#' @export
#' @rdname MSdata-class
setMethod("peakData<-", "MSdata",
          function(msdata, value) {
              msdata@peakData <- value
              if (validObject(msdata))
                  return(msdata)
          })

## Sample meta-data accesors and replacement methods
setGeneric("sampleData", function(msdata, ...) standardGeneric("sampleData"))

#' @export
#' @rdname MSdata-class
setMethod("sampleData", "MSdata",
          function(msdata) msdata@sampleData)

setGeneric("sampleData<-",
           function(msdata, value) standardGeneric("sampleData<-"))
#' @export
#' @rdname MSdata-class
setMethod("sampleData<-", "MSdata",
          function(msdata, value) {
              msdata@sampleData <- value
              if (validObject(msdata))
                  return(msdata)
          })

## Intensity matrix accesors and replacement methods
setGeneric("intMatrix", function(msdata, ...) standardGeneric("intMatrix"))
#' @export
#' @rdname MSdata-class
setMethod("intMatrix", "MSdata",
          function(msdata) msdata@intMatrix)

setGeneric("intMatrix<-",
           function(msdata, value) standardGeneric("intMatrix<-"))

#' @export
#' @rdname MSdata-class
setMethod("intMatrix<-", "MSdata",
          function(msdata, value) {
              msdata@intMatrix <- value
              if (validObject(msdata))
                  return(msdata)
          })

## Process Log accessors and addition methods
setGeneric("processLog", function(msdata, ...) standardGeneric("processLog"))
#' @export
#' @rdname MSdata-class
setMethod("processLog", "MSdata", 
          function(msdata) cat(msdata@processLog))

setGeneric("processLog<-",
           function(msdata, value) standardGeneric("processLog<-"))
#' @export
#' @rdname MSdata-class
setMethod("processLog<-", "MSdata",
          function(msdata, value) {
              msdata@processLog <- paste0(msdata@processLog, "\n\n", value)
              if (validObject(msdata))
                  return(msdata)
          })
##==============================================================================

setGeneric("peakNames", function(msdata, ...) standardGeneric("peakNames"))
#' @export
#' @rdname MSdata-class
setMethod("peakNames", "MSdata", 
          function(msdata) rownames(msdata@peakData))
		  
setGeneric("sampleNames", function(msdata, ...) standardGeneric("sampleNames"))
#' @export
#' @rdname MSdata-class
setMethod("sampleNames", "MSdata", 
          function(msdata) rownames(msdata@sampleData))