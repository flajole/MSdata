setGeneric("BiomassNorm", 
           function(msdata, ...) standardGeneric("BiomassNorm"))
#' Normalisation by biomass
#'
#' Normalisation of intensities by the list of sample masses/volumes/etc.
#' @param msdata \code{\link{MSdata-class}} object
#' @param biomass.list One of: 
#' \enumerate{
#' \item a path to the file containing the simple table with two columns: sample ID and biomass/volum/etc.
#' \item a name of the corresponding column in sample data table \code{sampleData(msdata)} (if biomass were uploaded there)
#' }
#' @return \code{\link{MSdata-class}} object with normalised intensity matrix
#' @name BiomassNorm
#' @export     	   
setMethod("BiomassNorm", "MSdata",
          function(msdata,
                   biomass.list) {
              .intMatrix <- msdata@intMatrix
              
			  if (biomass.list %in% names(sampleData(msdata))) {
					biomass.table <- sampleData(msdata)[biomass.list]
					msg <- ("Data are normalised by biomass")
			  } else if (file.exists(biomass.list)) {
					biomass.table <- suppressWarnings(read.table(biomass.list, row.names = 1, col.names = c("Biomass")))
					miss.samples <- setdiff(colnames(.intMatrix), rownames(biomass.table))  
					if (length(miss.samples) > 0)
						stop("Not all sample names are found in biomass table!")
					biomass.table <- biomass.table[colnames(.intMatrix), ]  
					msg <- paste0("Data are normalised by biomass; the following file is used as a list:\n", biomass.list)
			  } else {
				stop("Biomass list is neither a factor name in sample data nor existing file path.")
			  }
              
              masses <- matrix(biomass.table[[1]],
                               nrow=nrow(.intMatrix), 
                               ncol=ncol(.intMatrix), byrow=TRUE)
              meanMass <- mean(biomass.table[[1]])
              .intMatrix <- .intMatrix * meanMass / masses 
			  
			  intMatrix(msdata) <- .intMatrix
              processLog(msdata) <- paste0("Normalisation:\n", msg)
			  cat(msg)
              return(msdata)
          })