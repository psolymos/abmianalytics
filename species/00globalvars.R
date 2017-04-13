DATE <- as.Date(Sys.time(), tz=Sys.timezone(location = TRUE))
## load mefa4 package
stopifnot(require(mefa4))
## load package RODBC
#stopifnot(require(RODBC)) # working from csv files
## definition of global variables
do.prof <- FALSE  # do profiling or not
do.image <- TRUE  # save workspace image or not
use.date <- TRUE # use date as part of file names or not
check.completeness <- TRUE # use site list to check for completeness
combine.tables <- TRUE # combine site label info and spp matrix
DBVERSION <- "" # "y:/Oracle_access/DatabaseCombined_Clean_2014-04-28.accdb"
#D <- file.path(ROOT, "out") # directory to save files into
D <- OUTDIR
if (use.date) {
    #d <- paste(strsplit(date(), " ")[[1]][c(2,3,5)], collapse="-")
    d <- paste("_", DATE, sep="")
} else {
    d <- ""
}

LabelFun <- function(res) {
    ## ALPAC-SK issue
    levels(res$SITE_LABEL) <- gsub("ALPAC-SK", "ALPACSK", levels(res$SITE_LABEL))

    rr <- data.frame(t(sapply(levels(res$SITE_LABEL), function(z) strsplit(z, "_", fixed=TRUE)[[1]])))
    ## mos lichen PL-MH issue
    if (ncol(rr)>8) {
        tmp <- unname(apply(rr[,8:ncol(rr)], 1, paste, collapse="_"))
        rr <- rr[,1:8]
        rr[,8] <- as.factor(tmp)
    }
    colnames(rr) <- c("Protocol","OnOffGrid","DataProvider","SiteLabel","ABMIYear","Visit","SubType","SubTypeID")
    for (i in 1:ncol(rr))
        attr(rr[,i], "names") <- NULL
    tmp <- lapply(levels(rr$SiteLabel), function(z) strsplit(z, "-", fixed=TRUE)[[1]])
    rr$ClosestABMISite <- as.factor(sapply(tmp, function(zz)
        if (length(zz) == 1) zz[1] else zz[3])[as.integer(rr$SiteLabel)])
    rr$OGSeqenceID <- as.factor(sapply(tmp, function(zz)
        if (length(zz) == 1) NA else zz[4])[as.integer(rr$SiteLabel)])
    rr$Label <- as.factor(with(rr,
        paste(Protocol, OnOffGrid, DataProvider, SiteLabel, ABMIYear, Visit, SubType, SubTypeID, sep="_")))
    rr$Label2 <- as.factor(with(rr,
        paste(Protocol, OnOffGrid, DataProvider, SiteLabel, ABMIYear, Visit, sep="_")))

    ## ALPAC-SK issue
    ii <- grep("ALPACSK", rr$SiteLabel)
    rr$ClosestABMISite[ii] <- NA
    rr$ClosestABMISite <- droplevels(rr$ClosestABMISite)
    ## restore labels to avoide cascading problems
    levels(rr$SiteLabel) <- gsub("ALPACSK", "ALPAC-SK", levels(rr$SiteLabel))
    levels(rr$Label) <- gsub("ALPACSK", "ALPAC-SK", levels(rr$Label))
    levels(rr$Label2) <- gsub("ALPACSK", "ALPAC-SK", levels(rr$Label2))

    rr
}
