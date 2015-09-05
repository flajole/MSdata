
setGeneric("QuantNorm", 
           function(object, ...) standardGeneric("StandNorm"))
setMethod("QuantNorm", "MSdata",
          require('preprocessCore');






if(!is.null(ref)){
    if(grpRef == 'T'){
        grp.inx <- cls == ref;
        ref.smpl <- apply(proc[grp.inx, ], 2, mean);
    }else{
        ref.smpl <- proc[ref,];
    }
    data<-t(apply(data, 1, ProbNorm, ref.smpl));
    rownm<-"Probabilistic Quotient Normalization";
    
    
}
QuantileNormalize <- function(){
    data<-dataSet$proc;
    cls <- dataSet$proc.cls;
    cls.lvl <- levels(cls);
    
    # first log normalize
    data <- glog(data);
    
    require('preprocessCore');
    
    # normalize within replicates
    #for (lv in cls.lvl){
    #    sub.inx <- dataSet$proc.cls == lv;
    #    data[sub.inx, ] <- t(normalize.quantiles(t(data[sub.inx, ]), copy=FALSE));
    #}
    data <- t(normalize.quantiles(t(data), copy=FALSE));
    
    dataSet$norm <<- as.data.frame(data);
    dataSet$cls <<- cls;
    dataSet$rownorm.method<<-NULL;
    dataSet$colnorm.method<<-NULL;
    dataSet$combined.method<<-TRUE;
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n<-45
m<-501
rownames(data.mat) <- paste0("peak", 1:m)
colnames(data.mat) <- paste0("sample", 1:n)

smeta <- read.table("clipboard", header=T)
rownames(smeta) <- paste0("sample", 1:n)
pmeta <- read.table("clipboard", sep="\t", header=F)
rownames(pmeta) <- paste0("peak", 1:m)
names(pmeta) <- c("peakID", "m/z", "RI") 
ex <- MSdata(intMatrix = data.mat, peakData = pmeta, sampleData = smeta, processLog = "")
ex@smeta

> n <- 10
> m <- 6
> marray <- matrix(rnorm(n * m, 10, 5), ncol = m)
> pmeta <- data.frame(sampleId = 1:m,
                      + condition = rep(c("WT", "MUT"), each = 3))
> rownames(pmeta) <- colnames(marray) <- LETTERS[1:m]
> fmeta <- data.frame(geneId = 1:n,
                      + pathway = sample(LETTERS, n, replace = TRUE))
> rownames(fmeta) <-
    
###############################################################################
#MEMORY FUNCTION FOR SIMULATING
remember <- function() {
    memory <- list()
    f <- function(...) {
        # This is inefficient!
        memory <<- append(memory, list(...))
        invisible()
    }
    
    structure(f, class = "remember")
}  
    # Simulate XCMS N times with different parameters (N=200 by default)
The randomized ???200 parameter
combinations were uniformly distributed inside the parameter's
hyper-cube in the ranges as follows:

(1) Parameters of the "xcmsSet" method: 
fwhm = (5:50)
step = (0.02:0.2)
steps = (1:5)
mzdiff = (-0.1:0.2)
snthresh = (2:20)
profmethod = ("bin","binlin")
max = (500:5000)

(2) Parameters of the "group" method: 
bw = (1:50)
mzwid = (0.05:0.3)
minfrac =(0.1:0.5)
minsamp = (1:3)
max = (2:10)

(3) Parameters of the "retcor" method: 
missing = (1:29)
extra = (1:29)
method = ("loess", "linear")
span = (0.1:2.0)
family = ("symmetric", "gaussian")

# Write parameters and calculated Zcorr to the table
# If all Zcorr<150 repeat simulation till Zcorr>150 or K more times

data.mat<-data.matrix(read.table("clipboard"))
median.zero.control <- apply(data.mat[,1:5], MARGIN=1, median, na.rm=TRUE)
median.all.controls <- apply(data.mat[,1:15], MARGIN=1, median, na.rm=TRUE)

 
data.norm<-data.mat

median0 <- apply(data[,ref], MARGIN=1, median, na.rm=TRUE)

median1 <-apply(data.mat[,6:10], MARGIN=1, median, na.rm=TRUE)
median3 <-apply(data.mat[,11:15], MARGIN=1, median, na.rm=TRUE)

data.norm[,16:20]<-data.norm[,16:20]/matrix(median0, ncol=5, nrow=501, byrow=F)
data.norm[,31:35]<-data.norm[,31:35]/matrix(median0, ncol=5, nrow=501, byrow=F)
data.norm[,21:25]<-data.norm[,21:25]/matrix(median1, ncol=5, nrow=501, byrow=F)
data.norm[,36:40]<-data.norm[,36:40]/matrix(median1, ncol=5, nrow=501, byrow=F)
data.norm[,26:30]<-data.norm[,26:30]/matrix(median3, ncol=5, nrow=501, byrow=F)
data.norm[,41:45]<-data.norm[,41:45]/matrix(median3, ncol=5, nrow=501, byrow=F)
data.norm[,16:45] <- data.mat[,16:45]/matrix(apply(data.norm[,16:45], MARGIN=2, median, na.rm=TRUE), ncol=30, nrow=501)

normd<-data.norm/matrix(median0, ncol=45, nrow=501, byrow=F)
data.norm <- data.norm/matrix(apply(normd, MARGIN=2, median, na.rm=TRUE), ncol=45, nrow=501)


QuotNorm

function (X, method) 
{
    if (nargs() < 2) {
        stop("incorrect number of input parameters")
    }
    obs <- dim(X)[1]
    dimm <- dim(X)[2]
    factors <- rep(NaN, obs)
    for (i in 1:obs) {
        factors[i] <- sum(X[i, ])
        X[i, ] <- X[i, ]/factors[i]
    }
    switch(method, total = return(), prob = {
        X[0 == X] <- 1e-08
        if (nargs() < 3) {
            normRef <- X[1, ]
        }
        F <- X/(matrix(rep(normRef, each = obs), ncol = length(normRef)))
        for (i in 1:obs) {
            X[i, ] <- 10000 * X[i, ]/median(F[i, ])
            factors[i] <- (factors[i] * median(F[i, ]))/10000
        }
    })
    return(list(Sp = X, factors = factors, NormSp = normRef))
}