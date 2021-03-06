#' Apply Functions Over Replicate Groups
#' 
#' Return a matrix of the same format as intensity matrix of \code{\link[MSdata]{MSdata}}
#' object with values 
#' @param msdata An object of \code{\link[MSdata]{MSdata}} class.
#' @param fun The function to be applied.
#' @export
setGeneric("repapply", 
           function(msdata, ...) standardGeneric("repapply"))

setMethod("repapply", "MSdata",
          function(msdata, fun, ...) {
              .intMatrix <- intMatrix(msdata)
              if (is.null(sampleData(msdata)$ReplicationGroup)) {
                  stop("Set replication group factor first (for details see ?SetRepGroup)")
              }
              reps <- sampleData(msdata)$ReplicationGroup
              
              out <- matrix(ncol = ncol(.intMatrix), nrow=nrow(.intMatrix))
              colnames(out) <- sampleNames(msdata)
              rownames(out) <- peakNames(msdata)
              for (i in levels(reps)) {
                  .out <- apply(.intMatrix[ , reps == i], 1, fun, ...)
                  if (class(out) != class(.out)) { class(out) <- class(.out) }
                  out[, reps==i] <- .out
              }
              
              return(out)
          })