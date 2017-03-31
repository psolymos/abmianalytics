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
save(all_res, tax, vals, file=file.path(ROOT, "josm2", "josm-reseffects.Rdata"))

## determining associations for species

library(ResourceSelection)

#spp <- "AMRO"

s_fun <- function(spp) {
    DAT <- xnn[,c("ROAD01"),drop=FALSE]
    rownames(DAT) <- rownames(xnn)
    DAT$Y <- yyn[,spp]
    DAT$ST <- 0
    #DAT <- DAT[bb[,j],]
    DAT1 <- DAT[rep(seq_len(nrow(DAT)), DAT$Y),]
    DAT1$ST <- 1
    DAT <- rbind(DAT1, DAT)
    m <- rsf(ST ~ ROAD01, DAT, m=0, B=0)
    #exp(coef(m))
    #mm <- glm(ST ~ ROAD01, DAT, family="binomial")
    #exp(coef(mm))

    #fU <- sum(DAT$ROAD01[DAT$ST==1]) / sum(DAT$ST)
    #fA <- sum(DAT$ROAD01[DAT$ST==0]) / sum(1-DAT$ST)
    #s <- fU / fA
    #c(s=s, fU=fU, fA=fA)
    unname(coef(m))
}

sroad <- pbsapply(SPP, s_fun)
pdet <- colSums(yyn[,SPP])/nrow(yyn)
## sroad=log(pdet1/pdet0)
ii1 <- xnn$PCODE == "BBSAB"
ii0 <- !ii1 & xnn$ROAD01==0
pdet1 <- colSums(yyn[ii1,SPP])/sum(ii1)
pdet0 <- colSums(yyn[ii0,SPP])/sum(ii0)

save(sroad, pdet, pdet0, pdet1, file=file.path(ROOT, "josm2", "josm-sroad.Rdata"))


## combining all estimates

library(plotrix)
library(mefa4)
setwd("e:/peter/AB_data_v2016/out/birds/josm2")

## year effects
load("josm-yreffects.Rdata")

LEVEL <- 0.9
fstat <- function(x, level=0.95) {
    c(Mean=mean(x), Median=median(x), quantile(x, c((1-level)/2, 1 - (1-level)/2)))
}
tyr <- t(sapply(all_yr, fstat, level=LEVEL))

## species list
SPP <- rownames(tax)
tyr <- tyr[SPP,]
load("josm-sroad.Rdata")
tax$sroad <- sroad[rownames(tax)]
tax$pdet <- pdet[rownames(tax)]
TAX <- tax

## official BBS trend in BCR6/AB
tbbs <- read.csv("ron_bbs_t20170330.csv")
compare_sets(SPP, tbbs$SpeciesID)
setdiff(SPP, tbbs$SpeciesID)
bbs_lt <- droplevels(tbbs[tbbs$timeFrame=="Long-term",])
bbs_st <- droplevels(tbbs[tbbs$timeFrame=="Short-term",])
bbs_lt <- bbs_lt[match(SPP, bbs_lt$SpeciesID),c("annualTrend", "lowerLimit", "upperLimit")]
bbs_st <- bbs_st[match(SPP, bbs_st$SpeciesID),c("annualTrend", "lowerLimit", "upperLimit")]
rownames(bbs_lt) <- rownames(bbs_st) <- SPP

## residual trend
load("josm-reseffects.Rdata")

all_res2 <- list()
for (i in rownames(vals))
    all_res2[[i]] <- t(sapply(all_res, function(z) fstat(z[[i]], level=LEVEL)))

D0 <- sapply(all_res2[1:6], function(z) z[,2])
Dhb <- sapply(all_res2[7:12], function(z) z[,2])
Dcl <- sapply(all_res2[13:18], function(z) z[,2])
Dhf <- sapply(all_res2[19:24], function(z) z[,2])
D0[D0 > 100] <- 100
Dhb[Dhb > 100] <- 100
Dcl[Dcl > 100] <- 100
Dhf[Dhf > 100] <- 100
colnames(D0) <- colnames(Dhb) <- colnames(Dcl) <- colnames(Dhf) <- levels(vals$part)

## changing averages across the board
op <- par(mfrow=c(2,2), las=1, mar=c(5,8,2,2))
boxplot(D0[,6:1], horizontal=TRUE, col="#ABD9E9", main="D0", xlab="% annual change")
abline(v=0, col="#D7191C", lwd=2)
boxplot(Dhb[,6:1], horizontal=TRUE, col="#ABD9E9", main="Dhb", xlab="% annual change")
abline(v=0, col="#D7191C", lwd=2)
boxplot(Dcl[,6:1], horizontal=TRUE, col="#ABD9E9", main="Dcl", xlab="% annual change")
abline(v=0, col="#D7191C", lwd=2)
boxplot(Dhf[,6:1], horizontal=TRUE, col="#ABD9E9", main="Dhf", xlab="% annual change")
abline(v=0, col="#D7191C", lwd=2)
par(op)

## D0-Dhf effects
tmp <- cbind(D0=D0[,"offroad-noCL"], Dhb=Dhb[,"offroad-noCL"],
    Dcl=Dcl[,"offroad-noCL"], Dhf=Dhf[,"offroad-noCL"])
ladderplot(tmp, pch=NA, col="#2C7BB6", ylab="% annual change", main="offroad-noCL")
abline(h=0, col="#D7191C", lwd=2)

## compare with BBS

rn <- intersect(rownames(Dhf), tbbs$SpeciesID)
tmp <- cbind(BBS_lt=bbs_lt[rn,1], BBS_st=bbs_st[rn,1],
    Dhf[rn,c("both-noCL", "BBS", "offroad-noCL")])
ladderplot(tmp, pch=NA, col="#2C7BB6", ylab="% annual change", main="")
abline(h=0, col="#D7191C", lwd=2)

tmp2 <- cbind(BBS_lt=bbs_lt[rn,1], BBS_st=bbs_st[rn,1],
    Dhf[rn,c("both", "BBS", "offroad", "CL")])

e <- new.env()
load("e:/peter/AB_data_v2016/out/3x7/veg-hf_3x7_fix-fire_fix-age0.Rdata", envir=e)
slt <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
nrn <- as.character(slt$SITE_ID[slt$NATURAL_REGIONS %in% c("Boreal", "Parkland", "Foothills")])
hfn <- c("RoadHardSurface", "RoadTrailVegetated", "RoadVegetatedVerge")
allhf <- c("BorrowpitsDugoutsSumps",
    "Canals", "CultivationCropPastureBareground", "HighDensityLivestockOperation",
    "IndustrialSiteRural", "MineSite", "MunicipalWaterSewage", "OtherDisturbedVegetation",
    "PeatMine", "Pipeline", "RailHardSurface", "RailVegetatedVerge",
    "Reservoirs", "RoadHardSurface", "RoadTrailVegetated", "RoadVegetatedVerge",
    "RuralResidentialIndustrial", "SeismicLine", "TransmissionLine",
    "Urban", "WellSite", "WindGenerationFacility", "CCDecid0", "CCDecidR",
    "CCDecid1", "CCDecid2", "CCDecid3", "CCDecid4", "CCMixwood0",
    "CCMixwoodR", "CCMixwood1", "CCMixwood2", "CCMixwood3", "CCMixwood4",
    "CCConif0", "CCConifR", "CCConif1", "CCConif2", "CCConif3", "CCConif4",
    "CCPine0", "CCPineR", "CCPine1", "CCPine2", "CCPine3", "CCPine4")
rd1999 <- sum(e$yearly_vhf[["1999"]][["veg_current"]][nrn, hfn])
rd2014 <- sum(e$yearly_vhf[["2014"]][["veg_current"]][nrn, hfn])
hf1999 <- sum(e$yearly_vhf[["1999"]][["veg_current"]][nrn, allhf])
hf2014 <- sum(e$yearly_vhf[["2014"]][["veg_current"]][nrn, allhf])
all1999 <- sum(e$yearly_vhf[["1999"]][["veg_current"]][nrn, ])
all2014 <- sum(e$yearly_vhf[["2014"]][["veg_current"]][nrn, ])

## annual rate of change
100*((rd2014/rd1999)^(1/(2014-1999))-1) # 0.8023384
(rd2014/rd1999)^(1/(2014-1999)) # 1.008023

(hf2014/hf1999)^(1/(2014-1999)) # 1.01053

d <- (1+tmp2[,"BBS"]/100) / (1+tmp2[,"offroad"]/100)
dc <- d / Deltap
plot(d, sroad[names(d)])
cor.test(d, sroad[names(d)]) # no correlation

## subsets make a difference, that is consistent across how residual is defined
## residual definition impacts mostly extreme estimates
## CL is influential, correlates best with off-road data
## with and without CL (Dhf), explain correlations
## lack of pattern re road associations
## HF calculations: no huge change in road strata to drive changes
## emphasize that BBS and bbs is correlated (but bbs < BBS)
## strata specific trend: math indicates that
## but is off-road trend reliable given temporal sparsity?
## based on CL vs. offroad there is big scatter but agreement on average
## reasons why we see strata specific trend?
## - geographic shift relative to road network (climate change)
## - habitat related rearrangements: suboptimal habitat density declining more
## but why is offroad trend positive??? -- not very reliable (CL average is neutral)
## is it a data artifact or is it real?
## - Hard to decide, but there is strong relationship with pdet
##   indicating that data is driving the extreme trends
##   or that rare species decline more: no because of funnel shape
## can it be disturbance other than roads? not much bigger changes there either
## correcting for pdet indicates that roadside trend might be -5%
##   offroad trend might be around 0%, definitely pointing towards concentration
##   in better habitats, need to come up with ways of testing it
##   e.g. annual variation in good/bad habitats, and trends over time in these

plot(BBS ~ BBS_st, tmp2) # bbs is smaller than BBS
abline(0,1,lty=2)
abline(lm(BBS ~ BBS_st, data.frame(tmp2)))

plot(CL ~ offroad, tmp2) # 1:1 but huge scatter
abline(0,1,lty=2)
abline(lm(CL ~ offroad, data.frame(tmp2)))

#plot(abs(tmp[,"offroad-noCL"])~pdet[rownames(tmp)]) # strong pattern
plot(tmp[,"offroad-noCL"]~pdet[rownames(tmp)]) # strong pattern
abline(h=0, lty=2)
abline(lm(tmp[,"offroad-noCL"]~pdet[rownames(tmp)]))

plot(tmp[,"BBS"]~pdet[rownames(tmp)]) # strong pattern
abline(h=0, lty=2)
abline(lm(tmp[,"BBS"]~pdet[rownames(tmp)]))
