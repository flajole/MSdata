setGeneric("SetRepGroup", 
           function(msdata, ...) standardGeneric("SetRepGroup"))

#' Set replication groups numbers 
#' 
#' If there are data about replication groups, corresponding column is renamed 
#' in a standard way to "ReplicationGroup".\cr
#' Otherwise function automatically adds a \code{ReplicationGroup} column into \code{sampleData}
#' table according to group factors combinations.
#' @param msdata \code{\link{MSdata-class}} object
#' @param repFac If you already have replication group numbers in sample metadata, set \code{repFac} - the name of this factoring variable in table.
#' @param impFac The vector of the names of factors which are taken into account 
#' during automatic replication group labeling. If \code{NULL}, all grouping factors are used.
#' @return \code{\link{MSdata-class}} object with \code{$ReplicationGroup} in sample data table
#' @name SetRepGroup
#' @export

setMethod("SetRepGroup", "MSdata",
          function(msdata,
                   repFac = NULL,
                   impFac = NULL) {
              .sampleData <- sampleData(msdata)
              if (is.null(impFac)) {
				impFac <- names(.sampleData)
				impFac <- impFac[impFac != "ReplicationGroup"]
			  }
			  
              if (!is.null(repFac)) {
                  if (repFac %in% names(.sampleData)) {
                     names(.sampleData)[names(.sampleData)==repFac] <- "ReplicationGroup"
					 msg <- paste0(repFac, " column in sampleData is set as ReplicationGroup column")
                  } else {
                      stop("Stated replication factor is not found is sample data.")
                  }
              } else {
                  .sampleData$ReplicationGroup <- RepGroup(.sampleData, impFac)
				  msg <- "ReplicationGroup column is created in sampleData"
              }
              sampleData(msdata) <- .sampleData
			  msdata@processLog <- paste0(msgata@processLog, "\n\n", msg)
			  cat(msg)
              return(msdata)
          })
		  
RepGroup <- function(sampleData, impFac = names(sampleData)) {
    missFac <- setdiff(impFac, names(sampleData))
    if (length(missFac) != 0) {
        impFac <- impFac[!(impFac %in% missFac)]
        warning ("Not all important factors are found in sample meta-data.
                 Missing factors are excluded.")
    }
    combFac <- expand.grid(sapply(impFac, 
                                  function(x) levels(as.factor(sampleData[[x]]))))
    repGroup <- apply(sampleData[,impFac], 1,
                      function(x) which(apply(combFac, 1, 
                                              function(y) all(x == y))))
    repGroup <- as.factor(repGroup)
    levels(repGroup)<- 1:length(levels(repGroup))
    return(repGroup)
    }