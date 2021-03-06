#' Upload MSdata 
#' 
#' Create a \code{\link{MSdata-class}} object from external data tables.
#' Automatically adds a \code{ReplicationGroup} column in \code{sampleData} table according to group factors combinations.
#' #' Remember that in file paths you should use either "/" or double "\", not just "\".
#' @param object One of: 
#' \enumerate{
#' \item a character vector of length 1 - just a file path name of one .csv or .txt data frame.
#' \item a list of paths to three files: matrix of intensities, sample metadata and peak metadata.
#' }

#' @param orientation Orientation of table, one of \code{"SamplesInCol"} or \code{"SamplesInRow"}
#' @param zeros.as.NA If \code{TRUE} then zero intensities are treated as NAs.
#' @return \code{\link{MSdata-class}} object
#' @name MSupload

setGeneric("MSupload", 
           function(object, ...) standardGeneric("MSupload"))

#' @param sampleDataLines The number of lines containing data about samples.
#' @param peakDataLines The number of lines containing data about peaks/compounds.
#' @param sampleNames If \code{TRUE}, sample names are taken from the first column/row of table. \cr
#' If \code{FALSE}, standard names like \code{"sample1"} are created.
#' @param peakNames If \code{TRUE}, sample names are taken from the first row/column of table. \cr
#' If \code{FALSE}, standard names like \code{"peak1"} are created.
#' @describeIn MSupload Upload intensity matrix with several lines of peak and sample metadata.
#' @export
setMethod("MSupload", "character",
          function(object, 
                   orientation = "SamplesInCol",
				   zeros.as.NA = TRUE,
                   sampleDataLines = 1, 
                   peakDataLines   = 1,
                   sampleNames = TRUE,
                   peakNames = TRUE) {
              
              match.arg(orientation, c("SamplesInCol", "SamplesInRow"))
              
              format.check <- readLines(object, n = 10)
              
              if (all(grepl(";", format.check))) {
                  sep <- ";"
              } else if (all(grepl(",", format.check))) {
                  sep <- ","
              } else {
                  sep <- ""
              }
              
              tab <- read.table(object,
                                sep = sep,
                                header = FALSE, 
                                check.names = FALSE, 
                                stringsAsFactors = FALSE,
                                colClasses = "character",
                                comment.char = "");
              
              tab <- as.matrix(tab)
              
              if (orientation == "SamplesInRow") tab <- t(tab)
              if (sampleDataLines == 1) {
                  .sampleData <- data.frame(tab[ (1:sampleDataLines), -(1:peakDataLines)], stringsAsFactors = FALSE)
              } else {
                  .sampleData <- data.frame(t(tab[ (1:sampleDataLines), -(1:peakDataLines)]), stringsAsFactors = FALSE)
              }
              
              .peakData   <- data.frame(tab[-(1:sampleDataLines),  (1:peakDataLines)], stringsAsFactors = FALSE)
              
              names(.sampleData) <- tab[1:sampleDataLines, peakDataLines]
              names(.peakData)   <- tab[sampleDataLines, 1:peakDataLines]
              
              if (sampleNames) {
                  rownames(.sampleData) <- .sampleData[ , 1]
                  .sampleData <- .sampleData[-1]
              } else {
                  rownames(.sampleData) <- paste0("sample", 1:nrow(.sampleData));
              }
              
              if (peakNames) {
                  rownames(.peakData) <- .peakData[ , 1]         
                  .peakData <- .peakData[-1]
              } else {
                  rownames(.peakData) <- paste0("peak",   1:nrow(.peakData))        
              }
              
              .intMatrix <- tab[-(1:sampleDataLines), -(1:peakDataLines)]
              .intMatrix <- apply(.intMatrix, 2, gsub, pattern = ",", replacement = ".")
              class(.intMatrix) <- "numeric"
              
              rownames(.intMatrix) <- rownames(.peakData)
              colnames(.intMatrix) <- rownames(.sampleData)
			  
              if (zeros.as.NA) .intMatrix[.intMatrix==0] <- NA
              #.sampleData$ReplicateGroup <- RepGroup(.sampleData)
              .processLog <- paste0("Data uploaded from the file:\n", object)
			  cat(.processLog)
              MSdata(intMatrix  = .intMatrix,
                     peakData   = .peakData,
                     sampleData = .sampleData,
                     processLog = .processLog)
          })


#' @describeIn MSupload Upload intensity matrix, sample data, peak data separately
#' @export
setMethod("MSupload", "list",
          function(object = list(intFile = "",
                                 sampleFile = "",
                                 peakFile = ""), 
                   orientation = "SamplesInCol",
				   zeros.as.NA = TRUE) {
              
              match.arg(orientation, c("SamplesInCol", "SamplesInRow"))
              
              if (length(object$intFile) == 0) {
                  stop("Please, enter the file path for the matrix of intensities/concentrations.")
              }
              
              fileFormat <- substr(object$intFile, nchar(object)-2, nchar(object))
              if (fileFormat == "txt"){
                  tab <- read.table(object, header=FALSE, check.names=F, colClasses = "character");
              } else if (fileFormat == "csv") { # note, read.csv is more than read.table with sep=","
                  tab <- read.csv(object, header=FALSE, check.names=F, colClasses = "character");
              } else{
                  stop("File format is not .txt or .csv")
              }
              
              tab <- as.matrix(tab)
              if (orientation == "SamplesInRow") tab <- t(tab)
              
              sampleNames <- tab[1, -1]
              peakNames <- tab[-1, 1]
              
              .intMatrix <- tab[-1,-1]
              class(.intMatrix) <- "numeric"
              rownames(.intMatrix) <- peakNames
              rownames(.intMatrix) <- sampleNames
              
              if (length(object$sampleFile) == 0) {
                  .sampleData <- as.data.frame(matrix(nrow = length(sampleNames), ncol = 0), row.names = sampleNames)
              } else {
                  fileFormat <- substr(object$sampleFile, nchar(object)-2, nchar(object))
                  if (fileFormat == "txt"){
                      .sampleData <- read.table(object, header=TRUE, check.names=F, colClasses = "character", row.names = sampleNames);
                  } else if (fileFormat == "csv") { # note, read.csv is more than read.table with sep=","
                      .SampleData <- read.csv(object, header=TRUE, check.names=F, colClasses = "character", row.names = sampleNames);
                  } else{
                      stop("File format is not .txt or .csv")
                  }
                  if (sampleNames == as.character(.sampleData[[1]])) .sampleData <- .sampleData[-1]
              }
              
              
              if (length(object$peakFile) == 0) {
                  .peakData <- as.data.frame(matrix(nrow = length(peakNames), ncol = 0), row.names = peakNames)
              } else {
                  fileFormat <- substr(object$peakFile, nchar(object)-2, nchar(object))
                  if (fileFormat == "txt"){
                      .peakData <- read.table(object, header=TRUE, check.names=F, colClasses = "character", row.names = peakNames);
                  } else if (fileFormat == "csv") { # note, read.csv is more than read.table with sep=","
                      .peakData <- read.csv(object, header=TRUE, check.names=F, colClasses = "character", row.names = peakNames);
                  } else{
                      stop("File format is not .txt or .csv")
                  }
                  if (peakNames == as.character(.peakData[[1]])) .peakData <- .peakData[-1]
                  
              }
              
			  if (zeros.as.NA) .intMatrix[.intMatrix==0] <- NA
              #.sampleData$ReplicateGroup <- RepGroup(.sampleData)
              .processLog <- paste0("Data uploaded from the files:\n", 
                                    object$intFile, " - intensity matrix\n",
                                    object$sampleFile, " - sample metadata\n",
                                    object$peakFile, " - peak metadata\n\n")
			  cat(.processLog)
              MSdata(intMatrix  = .intMatrix,
                     peakData   = .peakData,
                     sampleData = .sampleData,
                     processLog = .processLog)
          })




	  

# #' Transforms an object of xcmsSet class into object of MSdata class
# #'
# #' @param object method to order colors 
# #' 
# #' @param intensity the measure of intensity:
# #' \code{"into"} - integrated area of original (raw) peak
# #' \code{"intf"} - integrated area of filtered peak
# #' \code{"maxo"} - maximum intensity of original (raw) peak
# #' \code{"maxf"} - maximum intensity of filtered peak
# #' 
# #'
# #' @return \code{\link{MSdata-class}} object
# #' 
# #' @seealso \code{\link{MSdata-class}}
# #'
# #' @export

# into:  integrated area of original (raw) peak
# intf:  integrated area of filtered peak
# maxo:  maximum intensity of original (raw) peak
# maxf:  maximum intensity of filtered peak

# setMethod("MSupload", "xcmsSet",
          # function(object, 
                   # intensity = c("into", "intf", "maxo", "maxf")){
              # xset <- object
              # intensity <- match.arg(intensity)
              # .sampleData <- xset@phenoData
              # .peakData <- data.frame(xset@groups[,1:7])
              # rownames(.sampleData) = rownames(xset@phenoData)
              # if (is.NULL(rownames(xset@peaks))) {
                  # rownames(.peakData) <- paste0("peak",   1:nrow(.peakData))
              # } else {
                  # rownames(.peakData) <- rownames(xset@peaks)
              # }
              
              # .intMatrix <- matrix(NA, nrow = nrow(.peakData), ncol = nrow(.sampleData))
              # rownames(.intMatrix) <- rownames(.peakData)
              # colnames(.intMatrix) <- rownames(.sampleData)
              # for (i in 1:length(xset@groupidx)) {
                  # smpl <- xset@peaks[xset@groupidx[[i]], c("sample", intensity)]
                  # .intMatrix[i, smpl[ , "sample"]] <- smpl[ , intensity]
              # }
              # .processLog <- paste0("Data uploaded from xcmsSet\n\nList of original files:\n", paste(xset@filepaths, collapse="\n"))
              # MSdata(intMatrix  = .intMatrix,
                     # peakData   = .peakData,
                     # sampleData = .sampleData,
                     # processLog = .processLog)
          # })

