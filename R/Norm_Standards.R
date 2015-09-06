## Normalization by the list of standards written in file
## as table: name, mz, retention time

setGeneric("StandNorm", 
           function(msdata, ...) standardGeneric("StandNorm"))
		  
		   
#' Normalisation by the list internal standards
#'
#' Normalisation by the list of external or internal standards. 
#' In each sample all the intensities are divided by sum of standards intensities.
#' @param msdata \code{\link{MSdata-class}} object
#' @param standards.list The link to the file with the table of standards looking like: name, mz, retention time
#' @param meanInt Mean sum intensity of all the standards (ususally, from previous experiments).
#' @param recalculateMean If \code{TRUE} \{meanInt} is not used, but is recalculated from this particular data.
#' 
#' @return \code{\link{MSdata-class}} object with normalised intensity matrix
#' @export 		   		   

setMethod("StandNorm", "MSdata",
          function(msdata, standards.list, meanInt = 2000, recalculateMean = FALSE) {
              #load standards list, convert into seconds if necessary
              standards.table <- read.table(standards.list, col.names = c("chemical", "mz","rt"))
              if (all(standards.table[,"rt"] < 60)) standards.table[,"rt"] <- standards.table[,"rt"]*60
              
              standards <- c()
			  
              #find standards in data matrix
              for (i in 1:dim(standards.table)[1]) {
                  found.standard <- which(msdata@peakData[,"mz"] > standards.table[i, "mz"] - 0.01 & 
                                          msdata@peakData[,"mz"] < standards.table[i, "mz"] + 0.01 & 
                                          msdata@peakData[,"rt"] > standards.table[i, "rt"] - 10   & 
                                          msdata@peakData[,"rt"] < standards.table[i, "rt"] + 10)
                  if (length(found.standard) > 1) warning("Found more than one peak corresponding to standard ", i, " - ", standards.table[i, "chemical"]
                  ) else if (length(found.standard) < 1) warning("Peak corresponding to standard ", i, " - ", standards.table[i, "chemical"], " is not found"
                  ) else standards <- c(standards, found.standard)
              }
              
              
              if (is.null(standards)) {
                  warning ("No reliable standards found, data are left native")
              } else {
                  #recalculate mean for normalization
                  if (recalculateMean) meanInt <- mean(colSums(msdata[standards,], na.rm=TRUE))
                  #recalculate intensities
                  standards.sums <- matrix(colSums(msdata@intMatrix[standards,], na.rm=TRUE), 
                                           nrow=nrow(msdata@intMatrix), 
                                           ncol=ncol(msdata@intMatrix), byrow=TRUE)
                  msdata@intMatrix <- msdata@intMatrix * meanInt / standards.sums 
                  #exclude standard from data matrix
                  msdata <- msdata[-standards,]
              }
              
              return(msdata)
          })

		  
