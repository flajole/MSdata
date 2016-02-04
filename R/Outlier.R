setGeneric("OutlierUni", 
           function(msdata, ...) standardGeneric("OutlierUni"))

# 
#
#
# Set replication group factor
#
# Univariate outliers
# Each rep group each peak boxplot & Tukey biweight
# 

setMethod("OutlierUni", "MSdata",
          function(msdata) {
              .intMatrix <- intMatrix(msdata)
              if (is.null(sampleData(msdata)$ReplicationGroup)) {
                  stop("Set replication group factor first (for details see ?SetRepGroup)")
              }
              
              reps <- sampleData(msdata)$ReplicationGroup
              outlier <- .intMatrix
              outlier[ , ] <- NA
              class(outlier) <- "logical"
              for (i in levels(reps)) {
                  out.box <- apply(.intMatrix[ , reps == i], 1, outlier.box)
                  out.tukey <- apply(.intMatrix[ , reps == i], 1, outlier.tukey)
                  outlier[, reps==i] <- t(out.box & out.tukey)
              }
              
              .intMatrix[outlier] <- NA
              intMatrix(msdata) <- .intMatrix
              msg <- paste0(sum(outlier), " univariate outliers are replaced with NA")
              processLog(msdata) <- msg
              cat(msg)
              return(msdata)
          })

setGeneric("OutlierCV", 
           function(msdata, ...) standardGeneric("OutlierCV"))

setMethod("OutlierCV", "MSdata",
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

setMethod("Outlier.mult", "MSdata",
          function(msdata) {
              .intMatrix <- intMatrix(msdata)
              if (is.null(sampleData(msdata)$ReplicationGroup)) {
                  stop("Set replication group factor first (for details see ?SetRepGroup)")
              }
              
              # HAVE TO sort by replication group
              
              reps <- sampleData(msdata)$ReplicationGroup
              cv <- .intMatrix
              cv[ , ] <- 0
              for (i in levels(reps)) {
                  cv[ , reps == i] <- apply(.intMatrix[ , reps == i], 1, function(r) sd(r)/mean(r)))
              }
              cv.mean <- apply(cv, 1, mean)
              cv.tukey <- apply(cv, 1, tukey)
              cv.median <- apply(cv, 1, median)
              outlier <- 
                  outlier.box(cv.mean)+outlier.box(cv.tukey)+outlier.box(cv.median) >=2
              return(msdata[!outlier, ])
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