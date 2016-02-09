setGeneric("msNorm", 
           function(msdata, ...) standardGeneric("msNorm"))
#' MSdata normalisation
#'
#' Sample-wise normalisation methods.
#' @param msdata \code{\link{MSdata-class}} object
#' @param method Method of normalisation. In each sample data are normalised by one of:\cr
#' \code{"sum"} -  sum concentration of all compounds\cr
#' \code{"median"} - median concentration of all compounds\cr
#' \code{"refCompound"} - concentration of reference compound\cr
#' \code{"robustPeaks"} - normalization by robust peaks (see details in Meissner, Steinhauser, Dittmann 2015)
#' @param ref.cmpd For \code{"refCompound"} method, the name or the order number of reference compound.
#' 
#' @return \code{\link{MSdata-class}} object with normalised intensity matrix
#' @name msNorm
#' @seealso \code{\link{msNormStandard}}, \code{\link{msNormBiomass}}, \code{\link{msScaling}}, \code{\link{msTransform}}
#' @export     	   	
setMethod("msNorm", "MSdata",
          function(msdata,
                   method,
                   ref.cmpd = NULL) {
              match.arg(method, c("sum", "median", "refCompound", "robustPeaks"))
              .intMatrix <- msdata@intMatrix
                  
              if (method == "median") {
                  medians <- matrix(apply(.intMatrix, 2, median, na.rm = TRUE),
                                    nrow=nrow(.intMatrix), 
                                    ncol=ncol(.intMatrix), 
                                    byrow=TRUE)
                  meanInt <- mean(medians)
                  .intMatrix <- .intMatrix * meanInt / medians
                  msg <- "Data are normalised by sample median" 
                  
              } else if (method == "sum") {
                  sums <- matrix(colSums(.intMatrix, na.rm=TRUE),
                                 nrow=nrow(.intMatrix), 
                                 ncol=ncol(.intMatrix), 
                                 byrow=TRUE)
                  meanInt <- mean(sums)
                  .intMatrix <- .intMatrix * meanInt / sums
                  msg <- "Data are normalised by sample sum" 
                  
              } else if (method == "robustPeaks") {
                  tukeys <- repapply(msdata, tukey)
                  logratios <- log2(.intMatrix/tukeys)
                  robust <- apply(logratios, 2, tukey.full)
                  robust <- apply(robust, 1, function(row) !any(is.na(row)))
                  robust <- logratios[robust, ]
                  robust <- matrix(2^apply(robust, 2, mean),
                                   nrow=nrow(.intMatrix), 
                                   ncol=ncol(.intMatrix), 
                                   byrow=TRUE)
                  .intMatrix <- .intMatrix/robust
                  msg <- paste0("Normalization by ",  , " robust peaks")
                  
              } else if (method == "refCompound") {
                  if (is.null(ref.cmpd)) 
                      stop("Method is set to normalisation by feature, 
                           but reference feature is not selected!")
                  
                  if (!(ref.cmpd %in% rownames(.intMatrix)) & !(ref.cmpd %in% 1:nrow(.intMatrix))) {
                      stop("Reference feature is not found!")
                  } else if (ref.cmpd %in% rownames(.intMatrix)) { 
                      ref.num <- match("Protein_TSS", peakNames(msdata))
                      ref.name <- ref.cmpd
                  } else {
                      ref.name <- rownames(.intMatrix)[ref.cmpd]
                      ref.num <- ref.cmpd
                  }
                  
                  refints <- matrix(.intMatrix[ref.num, ],
                                    nrow=nrow(.intMatrix), 
                                    ncol=ncol(.intMatrix), 
                                    byrow=TRUE)
                  meanInt <- mean(.intMatrix[ref.num, ])
                  .intMatrix <- .intMatrix * meanInt / refints
                  msg <- paste0("Data are normalised by reference compound ", ref.name)
              } 
              
			  intMatrix(msdata) <- .intMatrix
              processLog(msdata) <- paste0("Normalisation:\n", msg)
              if (method == "refCompound") msdata <- msdata[-ref.num, ]
			  cat(msg)
              return(msdata)
          })