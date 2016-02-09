setGeneric("MSoutput", 
           function(msdata, ...) standardGeneric("MSoutput"))
#' Output MSdata in files
#' 
#' Writes MSdata to four files, containing intensity matrix, sample metadata,
#' peak metadata and processing log file.
#'
#' @param msdata An object of \code{\link[MSdata]{MSdata}} class.
#' @param dir The output directory path.
#' @param file File name prefix.
#' @name MSoutput
#' @export
setMethod("MSoutput", "MSdata",
          function (msdata, 
                    dir = "",
                    file = "MSdata"){
              file <- paste0(file, 
                             c("_matrix.csv", 
                               "_samples.csv",
                               "_peaks.csv", 
                               "_log.txt"))
              file <- file.path(dir, file)
              
              write.table(msdata@intMatrix,
                          file = file[1],
                          col.names = FALSE, row.names = FALSE, sep = ",")
              write.table(msdata@sampleData,
                          file = file[2],
                          col.names = FALSE, row.names = FALSE, sep = ",")
              write.table(msdata@peakData,
                          file = file[3],
                          col.names = FALSE, row.names = FALSE, sep = ",")
              cat(msdata@processLog, file = file[4])
              cat("MSdata object is saved in directory", dir)
              return();
          })