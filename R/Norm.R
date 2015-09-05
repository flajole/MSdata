setGeneric("Norm", 
           function(msdata, ...) standardGeneric("Norm"))

setMethod("Norm", "MSdata",
          function(msdata,
                   method,
                   biomass.list = NULL,
                   ref.cmpd = NULL) {
              .intMatrix <- msdata@intMatrix
              match.arg(method, c("Biomass", "Sum", "Median", "RefCompound"))
              
              if (method == "Biomass") {
                  if (is.null(biomass.list)) 
                      stop("Method is set to normalisation by biomass, 
                           but biomasses list is not selected!")
                  biomass.table <- read.table(biomass.list, sep="\t", col.names = c("Sample", "Biomass"))
                  if (all(standards.table[,"rt"] < 60)) standards.table[,"rt"] <- standards.table[,"rt"]*60
                  
                  masses <- matrix(biomass.table[["Biomass"]],
                                   nrow=nrow(.intMatrix), 
                                   ncol=ncol(.intMatrix), byrow=TRUE)
                  meanMass <- mean(biomass.table["Biomass"])
                  .intMatrix <- .intMatrix * meanMass / masses 
                  msg <- paste0("Data are normalised by biomass, file used as list: ", biomass.list)
                  
              } else if (method == "Median") {
                  medians <- matrix(apply(.intMatrix, 2, median, na.rm = TRUE)
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