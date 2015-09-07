setGeneric("Norm", 
           function(msdata, ...) standardGeneric("Norm"))
#' MSdata normalisation
#'
#' Sample-wise normalisations.
#' @param msdata \code{\link{MSdata-class}} object
#' @param method Method of normalisation. In each sample data are normalised by one of:\cr
#' \code{"Sum"} -  sum concentration of all compounds\cr
#' \code{"Median"} - median concentration of all compounds\cr
#' \code{"RefCompound"} - concentration of reference compound\cr
#' \code{"Biomass"} - sample biomass/volume/etc. \cr
#' @param biomass.list For \code{"Biomass"} method, the link to the file containing the simple table with two columns: sampleID and biomasses/volumes/etc.
#' @param ref.cmpd For \code{"RefCompound"} method, the name or number of reference compound.
#' 
#' @return \code{\link{MSdata-class}} object with normalised intensity matrix
#' @name Norm
#' @export 		   	
setMethod("Norm", "MSdata",
          function(msdata,
                   method,
                   biomass.list = NULL,
                   ref.cmpd = NULL) {
              match.arg(method, c("Biomass", "Sum", "Median", "RefCompound"))
			  .intMatrix <- msdata@intMatrix
              
              if (method == "Biomass") {
                  if (is.null(biomass.list)) 
                      stop("Method is set to normalisation by biomass, 
                           but biomasses list is not selected!")
                  biomass.table <- suppressWarnings(read.table(biomass.list, row.names = 1, col.names = c("Biomass")))
				  
				  miss.samples <- setdiff(colnames(.intMatrix), rownames(biomass.table))  
				  if (length(miss.samples) > 0)
					  stop("Not all sample names are found in biomass table!")
					  
				  biomass.table <- biomass.table[colnames(.intMatrix), ]  
				  
                  masses <- matrix(biomass.table,
                                   nrow=nrow(.intMatrix), 
                                   ncol=ncol(.intMatrix), byrow=TRUE)
                  meanMass <- mean(biomass.table)
                  .intMatrix <- .intMatrix * meanMass / masses 
                  msg <- paste0("Data are normalised by biomass; the file used as list: ", biomass.list)
                  
              } else if (method == "Median") {
                  medians <- matrix(apply(.intMatrix, 2, median, na.rm = TRUE),
                                    nrow=nrow(.intMatrix), 
                                    ncol=ncol(.intMatrix), 
                                    byrow=TRUE)
                  meanInt <- mean(medians)
                  .intMatrix <- .intMatrix * meanInt / medians
                  msg <- "Data are normalised by sample median" 
                  
              } else if (method == "Sum") {
                  sums <- matrix(colSums(.intMatrix, na.rm=TRUE),
                                 nrow=nrow(.intMatrix), 
                                 ncol=ncol(.intMatrix), 
                                 byrow=TRUE)
                  meanInt <- mean(sums)
                  .intMatrix <- .intMatrix * meanInt / sums
                  msg <- "Data are normalised by sample sum" 
                  
              } else if (method == "RefCompound") {
                  if (is.null(ref.cmpd)) 
                      stop("Method is set to normalisation by feature, 
                           but reference feature is not selected!")
                  if (!(ref.cmpd %in% rownames(.intMatrix)) & !(ref.cmpd %in% 1:nrow(.intMatrix))) 
                      stop("Reference feature is not found!")
                  refints <- matrix(.intMatrix[ref.cmpd, ],
                                    nrow=nrow(.intMatrix), 
                                    ncol=ncol(.intMatrix), 
                                    byrow=TRUE)
                  meanInt <- mean(.intMatrix[ref.cmpd, ])
                  .intMatrix <- .intMatrix * meanInt / refints
                  msg <- paste0("Data are normalised by reference compound ", ref.cmpd)
              } 
              
              .processLog <- paste0(processLog(msdata), "\n\n Normalisation: ", msg)
              
              MSdata(intMatrix  = .intMatrix,
                     peakData   = peakData(msdata),
                     sampleData = sampleData(msdata),
                     processLog = .processLog)
               
          })