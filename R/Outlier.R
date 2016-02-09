setGeneric("msOutlierUni", 
           function(msdata, ...) standardGeneric("msOutlierUni"))

#' Detect Univariate Outliers
#' 
#' The function detects univariate outliers among peak intensity values
#' within each replication group. Outlying values are replaced with NAs.
#' 
#' @param msdata \code{\link{MSdata-class}} object
#' @param method Method of outlier detection, one of:\cr
#' \code{"tukey"} - Tukey's biweight.
#' \code{"box"} - boxplot statistic; values higher than 
#' (75% quantile + 1.5*IQR) or lower than (25% quantile - 1.5*IQR) 
#' are considered as outliers. If the number of values is less than 6,
#' only the most outlying one will be removed.
#' \code{"both"} - only values detected by both algorithms are marked as outliers.
#' @seealso \code{\link{msOutlierMulti}}, \code{\link{msOutlierPeakCV}}
#' @export

setMethod("msOutlierUni", "MSdata",
          function(msdata, method = "both") {
              .intMatrix <- intMatrix(msdata)
              
              match.arg(method, c("tukey", "box", "both"))
              out.box <- repapply(msdata, outlier.box)
              out.tukey <- repapply(msdata, outlier.tukey)
              
              if (method == "tukey") {
                  outlier <- out.tukey
              } else if (method == "box") {
                  outlier <- out.box
              } else {
                  outlier <- out.tukey & out.box
              }
              
              .intMatrix[outlier] <- NA
              intMatrix(msdata) <- .intMatrix
              msg <- paste0(sum(outlier), " univariate outliers are replaced with NA")
              processLog(msdata) <- msg
              cat(msg)
              return(msdata)
          })

setGeneric("msOutlierPeakCV", 
           function(msdata, ...) standardGeneric("msOutlierPeakCV"))

#' Detect Outlying Peaks
#' 
#' The function detects reliably quantified peaks using their groupwise 
#' coefficient of variation (CV). For this, the mean average, the Tukey's biweight,
#' and the median over all groupwise CV values are estimated and outliers 
#' are removed by iterative boxplot statistics. 
#' Peaks revealing outlying values in the majority of cases (2 out of 3) are 
#' removed from further analysis. 
#' 
#' @param msdata \code{\link{MSdata-class}} object
#' @param method Method of outlier detection, one of:\cr
#' \code{"tukey"} - Tukey's biweight.
#' \code{"box"} - boxplot statistic; values more than 
#' (75% quantile + 1.5*IQR) or less than (25% quantile - 1.5*IQR) 
#' are considered as outliers. If the number of values is less than 6,
#' only the most outlying one will be removed.
#' \code{"both"} - only values detected by both algorithms are marked as outliers.
#' @seealso \code{\link{msOutlierUni}}, \code{\link{msOutlierPeakCV}}
#' @export

setMethod("msOutlierPeakCV", "MSdata",
          function(msdata) {
              .intMatrix <- intMatrix(msdata)
              if (is.null(sampleData(msdata)$ReplicationGroup)) {
                  stop("Set replication group factor first (for details see ?SetRepGroup)")
              }
              reps <- sampleData(msdata)$ReplicationGroup
              cv <- .intMatrix
              cv[ , ] <- 0
              for (i in levels(reps)) {
                  cv[ , reps == i] <- apply(.intMatrix[ , reps == i], 1, function(r) sd(r, na.rm = T)/mean(r, na.rm = T))
              }
              cv.mean <- apply(cv, 1, mean, na.rm = T)
              cv.tukey <- apply(cv, 1, tukey, na.rm = T)
              cv.median <- apply(cv, 1, median, na.rm = T)
              outlier <- outlier.box(cv.mean)+outlier.box(cv.tukey)+outlier.box(cv.median) >=2
              msdata <- msdata[!outlier, ]
              msg <- paste0(sum(outlier), " peaks are removed as outliers based on coefficient of variation")
              processLog(msdata) <- msg
              cat(msg)
              return(msdata)
          })

setGeneric("msOutlierMulti", 
           function(msdata, ...) standardGeneric("msOutlierMulti"))

#' Detect Multivariate outliers
#' 
#' The function detects multivariate outliers among samples.
#' 
#' @param msdata \code{\link{MSdata-class}} object
#' @param method Method of outlier detection, one of:\cr
#' \code{"mahdist"} - based on the Mahalanodis distance (\code{\link[rrcovHD]{OutlierMahdist}} is used)\cr
#' \code{"pcdist"} - algorithm proposed by Shieh and Hung (2009) (\code{\link[rrcovHD]{OutlierPCDist}} is used)\cr
#' \code{"pcout"} - algorithm proposed by Filzmoser, Maronna and Werner (2008) (\code{\link[rrcovHD]{OutlierPCOut}} is used)\cr
#' \code{"sign1"} - using spatial sings, computation is based on the Mahalanodis distance (\code{\link[rrcovHD]{OutlierSign1}} is used)\cr
#' \code{"sign2"} - using spatial sings, computation is based on principal components (\code{\link[rrcovHD]{OutlierSign2}} is used)
#' @param ... - other arguments passed to the function of the respective method
#' @seealso \code{\link{msOutlierUni}}, \code{\link{msOutlierMulti}}
#' @export

setMethod("msOutlierMulti", "MSdata",
          function(msdata, method = "pcdist", ...) {
            require("rrcovHD")
              .intMatrix <- intMatrix(msdata)
              match.arg(method, c("mahdist", "pcdist", "pcout", "sign1", "sign2"))
              if (method == "mahdist") {
                  outlier <- rrcovHD::OutlierMahdist(t(.intMatrix), ...)
              } else if (method == "pcdist") {
                  outlier <- rrcovHD::OutlierPCDist(t(.intMatrix), ...)
              } else if (method == "pcout") {
                  outlier <- rrcovHD::OutlierPCOut(t(.intMatrix), ...)
              } else if (method == "sign1") {
                  outlier <- rrcovHD::OutlierSign1(t(.intMatrix), ...)
              } else if (method == "sign2") {
                  outlier <- rrcovHD::OutlierSign2(t(.intMatrix), ...)
              }
              outlier <- !as.logical(getFlag(outlier))
              msg <- paste0("By using ", method, "method ", sum(outlier),
                            " samples are considered as multivariate outliers and removed: \n",
                            paste(sampleNames(msdata)[outlier], collapse = " "))
              msdata <- msdata[, !outlier]
              processLog(msdata) <- msg
              cat(msg)
              return(msdata)
          })


tukey.w <- function(x, c = 5, e = 0.0001, na.rm = FALSE) {
    M <- median(x, na.rm=na.rm)
    S <- median(abs(x - M), na.rm=na.rm)
    u <- (x - M)/(c*S + e)
    nonout <- abs(u) <= 1
    nonout[is.na(nonout)] <- FALSE
    w <- rep(0, length(x))
    w[nonout] <- (1 - u[nonout]^2)^2
    w[is.na(u)] <- NA
    return(w)
}
    
tukey <- function(x, c = 5, e = 0.0001, na.rm = FALSE) {
    w <- tukey.w(x, c, e, na.rm)
    return(sum(w*x, na.rm=na.rm)/sum(w, na.rm=na.rm))
}

outlier.box <- function(x) {
    upper <- quantile(x, 0.75, na.rm=T) + 1.5*IQR(x, na.rm=T)
    lower <- quantile(x, 0.25, na.rm=T) - 1.5*IQR(x, na.rm=T)
    outlier <- (x > upper) | (x < lower)
    if ((length(na.omit(x)) > 6) | (sum(outlier, na.rm=T) <= 1)) {
        return(outlier)
    } else {
        outlier <- rep(FALSE, length(x))
        outlier[which.max(abs(x - median(x, na.rm=T)))] <- TRUE
        return(outlier)
    }
}

outlier.tukey <- function(x) {
    outlier <- tukey.w(x, na.rm = T) == 0
    outlier[is.na(outlier)] <- FALSE
    return(outlier)
}

tukey.full <- function(x) {
    prev <- NA
    filt <- x
    while (!identical(filt, prev)){
        prev <- filt
        filt[outlier.tukey(filt)] <- NA
    }
    return(filt)
}