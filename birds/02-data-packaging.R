#' ---
#' title: "Bird data packaging"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Nov 28, 2018"
#' output: pdf_document
#' ---
#'
#' # Preamble
#'
#' In this script we filter, transform, mutate, and get some offsets calculated.
#'
library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")
knitr::opts_chunk$set(eval=FALSE)
load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData")


## dealing with sine NAs
## checking veg+soil+hf summaries
dd$vshf <- ifelse(rowSums(is.na(vc1)) > 0, NA, 0)
(aa <- data.frame(n_of_NAs=colSums(is.na(dd))))
## ROAD: needs HF info
## DATE/DATI: use constant sra
## XY: should probably be dropped too
dd <- droplevels(dd[!is.na(dd$X) & !is.na(dd$vshf),])
dd$MAXDUR[is.na(dd$MAXDUR)] <- 1 # ABMISM bits
## climate: probably outside of AB bound: can use nearest if distance is small
dd <- dd[!is.na(dd$pAspen) & !is.na(dd$NRNAME),]

dd <- dd[dd$YEAR %[]% c(1993, 2017),]
dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
dd$start <- as.POSIXlt(dd$DATI)$hour + as.POSIXlt(dd$DATI)$min / 60
keep <- is.na(dd$DATI)
keep[dd$JULIAN %[]% c(125, 200)] <- TRUE
keep[dd$start %[]% c(3, 12)] <- TRUE

dd <- droplevels(dd[keep,])

dd$YR <- dd$YEAR - min(dd$YEAR)

## shrink yy
yy <- yy[rownames(dd), colSums(yy > 0) >= 20]

dd <- data.frame(dd, vs0[rownames(dd),])

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]

dd$pWater <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_water"]])
dd$pRoad <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_road"]])
dd$pWell <- row_std(vc1[rownames(dd),])[,"WellSite"]
dd$pRoadVeg <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_road_veg"]])

lcc4 <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"LCC4"]))
tmp <- find_max(lcc4)
dd$LCC4 <- tmp$index
dd$LCC2 <- dd$LCC4
levels(dd$LCC2) <- c("OpenWet", "Forest", "Forest", "OpenWet")
dd$pClosed <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_closed"]])
dd$pHarest <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_harvest"]])
dd$pHF <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_HF"]])

dd$pWater_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_water"]])
dd$pWet_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_wet"]])
dd$pWetWater_KM <-dd$pWater_KM + dd$pWet_KM

dd$vegpt <- as.factor(tv$UseInAnalysisFine[match(dd$VEGHFAGEclass, rownames(tv))])
dd$soilpt <- as.factor(ts$UseInAnalysisCoarse[match(dd$SOILHFclass, rownames(ts))])

vc1r <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"UseInAnalysisFine"]))
vr1r <- row_std(groupSums(vr1[rownames(dd),], 2, tv[colnames(vr1),"UseInAnalysisFine"]))
sc1r <- row_std(groupSums(sc1[rownames(dd),], 2, ts[colnames(sc1),"UseInAnalysisCoarse"]))
sr1r <- row_std(groupSums(sr1[rownames(dd),], 2, ts[colnames(sr1),"UseInAnalysisCoarse"]))

## road
dd$ROAD[is.na(dd$ROAD) & dd$pRoad > 0.04] <- 1
dd$ROAD[is.na(dd$ROAD)] <- 0
table(dd$ROAD, dd$pRoad > 0.04, useNA="a")

## dominant veg type
tmp <- find_max(vc1r[,colnames(vc1r) %ni% c("Water","HWater", "SnowIce", "Bare",
    "Seismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
dd$vegc <- factor(as.character(tmp$index), levels(dd$vegpt))
#dd$vegc <- tmp$index
dd$vegv <- tmp$value
dd$vegw <- pmax(0, pmin(1, 2*dd$vegv-0.5))
## how often in the center
a <- table(pt=dd$vegpt,bf=droplevels(dd$vegc))
a <- a[colnames(a),]
sum(diag(a))/sum(a)

## modifier variables
## north only
dd$mEnSft <- vc1r[,"EnSoftLin"]
dd$mTrSft <- vc1r[,"TrSoftLin"]
dd$mSeism <- vc1r[,"Seismic"]
## north & south
dd$mWell <- vc1r[,"Well"]
dd$mHard <- vc1r[,"HardLin"] # optional, use ROAD instead
## south only
dd$mSoft <- dd$mSeism + dd$mEnSft + dd$mTrSft

## use in north
dd$useNorth <- TRUE
dd$useNorth[dd$NRNAME == "Grassland"] <- FALSE
dd$useNorth[dd$useNorth & (dd$pWater > 0.5 | dd$vegpt == "Water")] <- FALSE
dd$useNorth[dd$vegw == 0] <- FALSE

## dominant soil type
tmp <- find_max(sc1r[,colnames(sc1r) %ni% c("SoilWater", "SoilUnknown", "HWater",
    "SoftLin", "HardLin", "HFor", "Well")])
dd$soilc <- factor(as.character(tmp$index), levels(dd$soilpt))
#dd$soilc <- tmp$index
dd$soilv <- tmp$value
dd$soilw <- pmax(0, pmin(1, 2*dd$soilv-0.5))

a <- table(pt=dd$soilpt,bf=droplevels(dd$soilc))
a <- a[colnames(a),]
sum(diag(a))/sum(a)


dd$useSouth <- FALSE
dd$useSouth[dd$NRNAME %in% c("Grassland", "Parkland")] <- TRUE
dd$useSouth[dd$NSRNAME %in% c("Dry Mixedwood")] <- TRUE
dd$useSouth[dd$useSouth & dd$Y > 56.7] <- FALSE # ~80 points
dd$useSouth[dd$useSouth & (dd$pWater > 0.5 | dd$vegpt == "Water")] <- FALSE
dd$useSouth[dd$useSouth & sr1r[,"SoilUnknown"] > 0] <- FALSE
dd$useSouth[dd$useSouth & sc1r[,"HFor"] > 0] <- FALSE
dd$useSouth[dd$soilw == 0] <- FALSE

## explore
if (FALSE) {
ys <- groupSums(yy[dd$useSouth,], 1, dd$SS[dd$useSouth])
ys[ys > 0] <- 1
ys <- ys[,colSums(ys) >= 20]
ds <- nonDuplicated(dd, SS, TRUE)[rownames(ys),]

table(ds$soilc, ds$WELL)
hist(ds$soilv)

library(opticut)
library(parallel)
cl <- makeCluster(8)
o <- opticut(as.matrix(ys) ~ ROAD, data=ds, strata=ds$soilc, dist="binomial:cloglog",
    weights=ds$soilv, cl=cl)
stopCluster(cl)
plot(o,cex.axis=0.5)

yn <- groupSums(yy[dd$useNorth,], 1, dd$SS[dd$useNorth])
yn[yn > 0] <- 1
yn <- yn[,colSums(yn) >= 100]
dn <- nonDuplicated(dd, SS, TRUE)[rownames(yn),]

data.frame(table(dn$vegc))
hist(dn$vegv)

library(opticut)
library(parallel)
cl <- makeCluster(8)
o <- opticut(as.matrix(yn) ~ ROAD, data=dn, strata=dn$vegc, dist="binomial:cloglog",
    weights=dn$vegw, cl=cl)
stopCluster(cl)
plot(o,cex.axis=0.5)
}

## offsets
## JDAY
dd$JDAY <- dd$JULIAN / 365
## TSSR
library(maptools)
Coor <- as.matrix(dd[,c("X", "Y")])
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dd$srise <- NA
dd$srise[subset] <- sr
dd$TSSR <- (dd$start - dd$srise) / 24

dd$JDAY2 <- dd$JDAY^2
dd$TSSR2 <- dd$TSSR^2
dd$MAXDIS <- dd$MAXDIS / 100
table(dd$MAXDIS)

Xp <- cbind("(Intercept)"=1, as.matrix(dd[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=dd$pClosed,
    LCC2OpenWet=ifelse(dd$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(dd$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(dd$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(dd$LCC4=="Wet", 1, 0))
summary(Xp)
summary(Xq)

library(QPAD)
load_BAM_QPAD(version=3)
sppp <- intersect(colnames(yy), getBAMspecieslist())

OFF <- matrix(NA, nrow(dd), length(sppp))
rownames(OFF) <- rownames(dd)
colnames(OFF) <- sppp

#spp <- "OVEN"
for (spp in sppp) {
    cat(spp, "\n");flush.console()
    p <- rep(NA, nrow(dd))
    A <- q <- p

    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model
    mi <- bestmodelBAMspecies(spp, type="BIC",
        model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    #vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0

    p[!OKp] <- sra_fun(dd$MAXDUR[!OKp], cf0[1])
    unlim <- ifelse(dd$MAXDIS[!OKq] == Inf, TRUE, FALSE)
    A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dd$MAXDIS[!OKq]^2)
    q[!OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[!OKq], cf0[2]))

    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    p[OKp] <- sra_fun(dd$MAXDUR[OKp], phi1)
    unlim <- ifelse(dd$MAXDIS[OKq] == Inf, TRUE, FALSE)
    A[OKq] <- ifelse(unlim, pi * tau1^2, pi * dd$MAXDIS[OKq]^2)
    q[OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[OKq], tau1))

    ii <- which(p == 0)
    p[ii] <- sra_fun(dd$MAXDUR[ii], cf0[1])

    OFF[,spp] <- log(p) + log(A) + log(q)

}

(Ra <- apply(OFF, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))

## make BB

## make xlat/xlong/xClimate

## make wtAge and stand type indicators and fCC2

## wet/water at( KM scale in veg models?)
