library(xcms)
mzfiles <- list.files("C:/DATA/new plant lipids", pattern="*pos.mzXML$", full.names = TRUE)
targets.list <- ("C:/DATA/new plant lipids")
xset <- xcmsSet(mzfiles)



#generate random lists of parameters
pargen <- function(n, 
                   goal = c("findPeaks", "group", "retcor"),
                   
                   # Parameters of the "xcmsSet" method: 
                   fwhm = (5:50),
                   step = (0.02:0.2),
                   steps = (1:5),
                   mzdiff = (-0.1:0.2),
                   snthresh = (2:20),
                   profmethod = c("bin","binlin"),
                   max = (500:5000),
                   
                   # Parameters of the "group" method:
                   bw = (1:50),
                   mzwid = (0.05:0.3),
                   minfrac =(0.1:0.5),
                   minsamp = (1:3),
                   max = (2:10),
                   
                   #Parameters of the "retcor" method:
                   missing = (1:29),
                   extra = (1:29),
                   method = c("loess", "linear"),
                   span = (0.1:2.0),
                   family = c("symmetric", "gaussian")) {
    
    arglist <- as.list(environment(), all = TRUE)
    match.arg(goal)
    prunif <- function(n, x) {
        switch (class(x),
                numeric = runif(n, min(x), max(x)),
                integer = runif(n, min(x), max(x)),
                character = sample(x, n, replace = TRUE))
    }
    pnames <- switch(goal,
                     findPeaks = c("fwhm", "step", "steps", "mzdiff", "snthersh", "profmethod", "max"),
                     group = c("bw", "mzwid", "minfrac", "minsamp", "max"),
                     retcor = c("missing", "extra", "method", "span", "family"))
    pars <- as.data.frame.list(lapply(arglist[pnames], prunif, n = n))
}

pars.findPeaks <- pargen(n, "findPeaks")
pars.retcor <- pargen(n, "retcor")
pars.group <- pargen(n, "group")

#pars <- lapply(c("findPeaks", "group", "retcor"), pargen, n = n)
#names(par) c("findPeaks", "group", "retcor")
    
for (i in 1:n) {
    xset <- do.call (xcmsSet, c(list(files = mzfiles), par.findPeaks[i, ]))
    xset <- do.call (group, c(list(object = xset), par.group[i, ]))
    xset <- do.call (retcor, c(list(object = xset), par.retcor[i, ]))
    diffreport(xset, filebase = paste0("xcms_Re", as.character(i)))
}


(shell("C:\\DATA\\MetaboQC\\ac101216e_si_003\\runMetaboQCb.bat"))
QCout <- read.table("C:\\DATA\\MetaboQC\\ac101216e_si_003\\Arab4TissueMFmulti_RTmz_perSampleQC.txt", header=TRUE)
for (i in levels(QCout$InputFile)) {
    mean(subset(QCout, InputFile=i, select=Zcorr))
}

fwhm =  
step =  runif(0.02:0.2),
steps = runif(1:5),
mzdiff = runif(-0.1:0.2),
snthresh = runif(2:20),
profmethod = ("bin","binlin"),
max = runif(500:5000),

(2) Parameters of the "group" method: 
par.group
names() <- c("missing", "extra", "method", "span", "family)

bw = runif(1:50),
mzwid = runif(0.05:0.3),
minfrac =runif(0.1:0.5),
minsamp = runif(1:3),
max = runif(2:10),

(3) Parameters of the "retcor" method: 
par.retcor <- list()
missing = runif(1:29),
extra = runif(1:29),
method = runif("loess"),
span = runif(0.1:2.0),
family = runif("symmetric")