library(parallel)
library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"
#ROOT <- "~/Dropbox/Public"

level <- 0.9

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
#    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

en <- new.env()
load(file.path(ROOT, "data", "data-north.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY
off <- en$OFF
bb <- en$BB

## names etc
e <- new.env()
load(file.path(ROOT, "data", "data-wrsi.Rdata"), envir=e)
TAX <- droplevels(e$TAX)
TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")

rm(e, en)

## model for species
## subset: good models on spp website
splt <- read.csv("~/repos/abmispecies/_data/birds.csv")
SPP <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATTW", "BAOR", "BARS",
    "BAWW", "BBMA", "BBWA", "BBWO", "BCCH", "BHCO", "BHVI", "BLJA",
    "BLPW", "BOCH", "BRBL", "BRCR", "BTNW", "CAWA", "CCSP", "CEDW",
    "CHSP", "CMWA", "COGR", "CONW", "CORA", "COYE", "DEJU", "DOWO",
    "EAKI", "EAPH", "EUST", "EVGR", "FOSP", "GCKI", "GRAJ", "GRCA",
    "GRYE", "HAWO", "HETH", "HOWR", "KILL", "LCSP", "LEFL", "LEYE",
    "LISP", "MAWA", "MODO", "MOWA", "NOFL", "NOWA", "OCWA", "OSFL",
    "OVEN", "PAWA", "PHVI", "PIGR", "PISI", "PIWO", "PUFI", "RBGR",
    "RBNU", "RCKI", "RECR", "REVI", "RUBL", "RUGR", "RWBL", "SAVS",
    "SOSA", "SOSP", "SPSA", "SWSP", "SWTH", "TEWA", "TOSO", "TRES",
    "VATH", "VEER", "VESP", "WAVI", "WBNU", "WCSP", "WETA", "WEWP",
    "WISN", "WIWA", "WIWR", "WTSP", "WWCR", "YBFL", "YBSA", "YEWA",
    "YRWA")

tax <- droplevels(TAX[SPP, c("Spp","English_Name","Scientific_Name","Family_Sci")])
compare_sets(tax$Spp, as.character(splt$sppid))
compare_sets(tax$Spp, as.character(splt$sppid[splt$veghf.north]))
SPPkeep <- sort(intersect(tax$Spp, as.character(splt$sppid[splt$veghf.north])))
tax <- droplevels(tax[tax$Spp %in% SPPkeep, ])
SPP <- rownames(tax)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))

## spp specific output

#spp <- "BTNW"

all_yr <- list()
for (spp in SPP) {
cat(spp, "\n");flush.console()

resn <- loadSPP(file.path(ROOT, "results", "north", paste0("birds_abmi-north_", spp, ".Rdata")))
estYr <- getEst(resn, stage=which(names(modsn)=="Year"), na.out=FALSE, Xnn)

## Boreal year effect estimates

## 0.1* because it is annual and not decadal
apc <- 100 * (exp(0.1*estYr[,"YR"]) - 1)
all_yr[[spp]] <- apc
apcstat <- round(c(fstat(apc), summary(apc)), 3)

png(file.path(ROOT, "josm2", "yr", paste0(spp, ".png")))
d <- density(apc)
hist(apc, col="grey", xlab="% annual population change",
    main=tax[spp,"English_Name"], freq=FALSE, border=NA, ylim=c(0, max(d$y)))
lines(d)
rug(apc)
i <- which.min(abs(d$x - apcstat[2]))
lines(c(d$x[c(i,i)]), c(d$y[i], -0), col=2, lwd=2)
i <- which.min(abs(d$x - apcstat[3]))
lines(c(d$x[c(i,i)]), c(d$y[i], -0), col=2, lwd=1)
i <- which.min(abs(d$x - apcstat[4]))
lines(c(d$x[c(i,i)]), c(d$y[i], -0), col=2, lwd=1)
dev.off()
}
save(all_yr, tax, file=file.path(ROOT, "josm2", "josm-yreffects.Rdata"))

## Residual trend estimates

yr_fun <- function(i, subset=NULL, part=c("all", "bbs", "bam"), colD="Dhf") {
    part <- match.arg(part)
    if (is.null(subset))
        subset <- rep(TRUE, nrow(DAT))
    dat <- DAT
    dat$SUBSET <- subset
    dat$D <- dat[[colD]]
    dat <- dat[bb[,i],]
    dat <- dat[dat$SUBSET,,drop=FALSE]
    if (part=="bbs") # BBS only
        dat <- dat[dat$isBBS,,drop=FALSE]
    if (part=="bam") # non-BBS excluding roadside surveys
        dat <- dat[!dat$isBBS & dat$ROAD01==0,,drop=FALSE]
    if (part=="all") # non-BBS excluding roadside surveys
        dat <- dat[dat$isBBS | (!dat$isBBS & dat$ROAD01==0),,drop=FALSE]
    if (nrow(dat) < 1)
        return(NA)
    dat$logDoff <- log(dat$D) + dat$off
    mod <- glm(Y ~ YR, data=dat, offset=dat$logDoff, family=poisson)
    out <- 100 * (exp(0.1*coef(mod)[2]) - 1)
    #attr(out, "n") <- nrow(dat)
    out
}

vals <- expand.grid(part=c("both", "BBS", "offroad", "CL", "both-noCL", "offroad-noCL"),
    dens=c("D0", "Dhb", "Dcl", "Dhf"))
rownames(vals) <- paste(vals$part, vals$dens, sep="_")
Bmax <- 100

all_res <- list()
for (spp in SPP) {
cat("\n------------", spp, "------------\n");flush.console()

resn <- loadSPP(file.path(ROOT, "results", "north", paste0("birds_abmi-north_", spp, ".Rdata")))
est0 <- sapply(resn, "[[", "null")
estHb <- getEst(resn, stage=which(names(modsn)=="ARU"), na.out=FALSE, Xnn)
estCl <- getEst(resn, stage=which(names(modsn)=="Space"), na.out=FALSE, Xnn)
estHF <- getEst(resn, stage=which(names(modsn)=="HF"), na.out=FALSE, Xnn)

pr0 <- exp(est0)
col_keep <- colSums(abs(estHb) > 0) != 0
prHb <- exp(sapply(1:nrow(estHb), function(j)
    Xnn[,colnames(estHb[,col_keep,drop=FALSE]),drop=FALSE] %*%
    estHb[j,col_keep]))
col_keep <- colSums(abs(estCl) > 0) != 0
prCl <- exp(sapply(1:nrow(estCl), function(j)
    Xnn[,colnames(estCl[,col_keep,drop=FALSE]),drop=FALSE] %*%
    estCl[j,col_keep]))
col_keep <- colSums(abs(estHF) > 0) != 0
prHF <- exp(sapply(1:nrow(estHF), function(j)
    Xnn[,colnames(estHF[,col_keep,drop=FALSE]),drop=FALSE] %*%
    estHF[j,col_keep]))

DAT <- droplevels(xnn[,c("PKEY","SS","PCODE","YEAR","YR","ROAD01")])
rownames(DAT) <- rownames(xnn)
DAT$Y <- yyn[,spp]
DAT$Y1 <- ifelse(yyn[,spp]>0, 1, 0)
DAT$off <- off[,spp]
DAT$isBBS <- DAT$PCODE == "BBSAB"
with(DAT, table(isBBS, ROAD01))
DAT$D0 <- mean(pr0)
DAT$Dhb <- rowMeans(prHb)
DAT$Dcl <- rowMeans(prCl)
DAT$Dhf <- rowMeans(prHF)
## North only
#DAT$sset <- xnn$NRNAME != "Grassland" & xnn$POINT_Y > 50


cl <- makeCluster(4)
tmp <- clusterExport(cl, c("DAT", "bb"))
res <- list()
for (j in 1:nrow(vals)) {
    jj <- rownames(vals)[j]
    cat(jj, "\n");flush.console()
    PART <- switch(as.character(vals$part[j]),
        "both"="all",
        "BBS"="bbs",
        "offroad"="bam",
        "CL"="bam",
        "both-noCL"="all",
        "offroad-noCL"="bam")
    SUBSET <- NULL
    if (as.character(vals$part[j]) == "CL")
        SUBSET <- DAT$PCODE == "CL"
    if (as.character(vals$part[j]) %in% c("both-noCL", "offroad-noCL"))
        SUBSET <- DAT$PCODE != "CL"
    res[[jj]] <- pbsapply(1:Bmax, yr_fun, cl=cl,
        subset=SUBSET, part=PART, colD=as.character(vals$dens[j]))
    #res[[j]] <- yr_fun(1, subset=SUBSET, part=PART, colD=as.character(vals$dens[j]))
}
stopCluster(cl)

#vals$est <- unlist(res)
#vals$n <- sapply(res, attr, "n")
#vals

all_res[[spp]] <- res
}
save(all_res, tax, file=file.path(ROOT, "josm2", "josm-reseffects.Rdata"))
