library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

load("d:/abmi/AB_data_v2019/data/misc/bg/veghf-summaries.RData")
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]

## --------------- point level variables

## 150m point
vc1 <- ddp17$veg_current
vc1r <- row_std(groupSums(vc1, 2, tv[colnames(vc1),"UseInAnalysisFine"]))

tmp <- find_max(vc1r[,colnames(vc1r) %ni% c("Water","HWater", "SnowIce", "Bare",
    "Seismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
dd <- data.frame(vegc = factor(as.character(tmp$index), colnames(vc1r)))

cc <- paste0(ifelse(tv[colnames(vc1), "is_harvest"], "CC", ""),
    as.character(tv[colnames(vc1), "UseInAnalysisFine"]))
tmp <- row_std(groupSums(vc1, 2, cc))
tmp <- find_max(tmp[,colnames(tmp) %ni% c("Water","HWater", "SnowIce", "Bare",
    "Seismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
dd$vegccc <- factor(as.character(tmp$index), c(levels(dd$vegc),
    paste0("CC", c(c("Spruce","Decid","Mixedwood","Pine")))))
table(cc=dd$vegccc, not=dd$vegc)

#levels(dd$vegc) <- levels(dd$vegccc)
#' Indicator variable for harvest area
dd$isCC <- startsWith(as.character(dd$vegccc), "CC")
tmp <- dd$vegccc
levels(tmp) <- gsub("CC", "", levels(dd$vegccc))
dd$vegc[dd$isCC] <- tmp[dd$isCC]

dd$vegc <- droplevels(dd$vegc)
dd$vegc <- relevel(dd$vegc, "Decid")
dd$vegccc <- droplevels(dd$vegccc)
dd$vegccc <- relevel(dd$vegccc, "Decid")

dd$isMix <- ifelse(dd$vegc == "Mixedwood", 1L, 0L)
dd$isWSpruce <- ifelse(dd$vegc == "Spruce", 1L, 0L)
dd$isPine <- ifelse(dd$vegc == "Pine", 1L, 0L)
dd$isBSpruce <- ifelse(dd$vegc == "BSpr", 1L, 0L)
dd$isLarch <- ifelse(dd$vegc == "Larch", 1L, 0L)
dd$isBSLarch <- ifelse(dd$vegc %in% c("BSpr", "Larch"), 1L, 0L)
dd$isUpCon <- ifelse(dd$vegc %in% c("Spruce", "Pine"), 1L, 0L)
dd$isCon <- ifelse(dd$vegc %in% c("BSpr", "Larch",
    "Spruce", "Pine"), 1L, 0L)

ac <- as.character(tv[colnames(vc1), "AGE"])
ac[is.na(ac)] <- ""
vc1age <- row_std(groupSums(vc1, 2, ac))
## exclude unknown (0) and non-forest (blank)
AgePtCr <- t(vc1age[,c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")])
AgeMin <- structure(c(0,10,20,40,60,80,100,120,140,160)/200,
    names=c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
dd$wtAge <- colSums(AgePtCr * AgeMin) / colSums(AgePtCr)
dd$wtAge[is.na(dd$wtAge)] <- 0
dd$isFor <- dd$vegc %in% c("Spruce","Decid","Mixedwood","Pine","BSpr", "Larch")
dd$wtAge[!dd$isFor] <- 0
dd$wtAge2 <- dd$wtAge^2
dd$wtAge05 <- sqrt(dd$wtAge)

MAXFOR <- 50/200
dd$fCC1 <- 0
dd$fCC1[dd$isCC==1] <- pmax(0, 1 - (dd$isCC * dd$wtAge/MAXFOR)[dd$isCC==1])
#plot(fCC1 ~ wtAge, dd[dd$isCC==1,])
## fCC2: Dave Huggard's recovery trajectories
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3,
    82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85,
    88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100
#data.frame(Age=age*200, wtAge=age, fCC2_decid=decid, fCC2_decid=decid)
dd$fCC2 <- 0
tmp1 <- approxfun(age, decid)(dd$wtAge)
tmp1[is.na(tmp1)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- dd$isFor & !dd$isCon & dd$isCC==1
dd$fCC2[ii] <- tmp1[ii]
tmp2 <- approxfun(age, conif)(dd$wtAge)
tmp2[is.na(tmp2)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- dd$isCon & dd$isCC==1
dd$fCC2[ii] <- tmp2[ii]
#plot(dd$wtAge, dd$fCC1, col=2, pch=".")
#points(dd$wtAge, dd$fCC2, col=4, pch=".")
sum(is.na(dd$fCC2))
by(dd$wtAge*200, list(veg=interaction(dd$isCC,dd$vegc,drop=TRUE)), fstat, level=1)

dd$mEnSft <- vc1r[,"EnSoftLin"]
dd$mTrSft <- vc1r[,"TrSoftLin"]
dd$mSeism <- vc1r[,"Seismic"]
#' Modifiers used in the north and the south
dd$mWell <- vc1r[,"Well"]
dd$mHard <- vc1r[,"HardLin"] # optional, use ROAD instead
dd$mSoft <- dd$mSeism + dd$mEnSft + dd$mTrSft
dd$ROAD <- ifelse(dd$mHard >= 0.01, 1, 0)
dd$CMETHOD <- factor("SM", c("HS", "RF", "SM"))
rownames(dd) <- rownames(vc1)

dd$YR <- 2017 - 1993

## ----------- climate stuff

gg <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/bg-from-peter-GRIDS-2019-11-14.csv")
vv <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/bg-from-peter-IDS-2019-11-14.csv")

library(cure4insect)
load_common_data()
rt <- .read_raster_template()
XY <- get_id_locations()

xy <- vv[,c("x_gis", "y_gis")]
colnames(xy) <- c("X", "Y")
coordinates(xy) <- ~ X + Y
proj4string(xy) <- proj4string(rt)
xy2 <- spTransform(xy, proj4string(XY))
xyc <- data.frame(coordinates(xy2))
rownames(xyc) <- vv$id_final

for (i in c("AHM", "FFP", "Eref", "MAP", "MAT", "MCMT", "MWMT")) {
    r <- raster(paste0("d:/spatial/ab-climate/", i, ".asc"))
    xyc[[i]] <- extract(r, xy)
}
colnames(xyc)[colnames(xyc)=="Eref"] <- "PET"

xyt <- transform_clim(xyc)
xyt <- xyt[rownames(dd),]
dd <- data.frame(dd, xyt)

## ----------- landscape level variables: SSH

## 600x600m square
vc2 <- dd17$veg_current
ao <- paste0(as.character(tv[colnames(vc2), "UseInAnalysisFine"]),
    ifelse(tv[colnames(vc2), "MatureOld"], "O", ""))
SSH <- row_std(groupSums(vc2[rownames(dd),], 2, ao))
SSH <- SSH[,colnames(SSH) %ni% c("Water","HWater", "SnowIce", "Bare",
    "Seismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")]
SSH <- SSH[rownames(dd),]

dd$SSH_KM <- 0
dd$SSH05_KM <- sqrt(dd$SSH_KM)
dd$pWater_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_water"]])
dd$pWater2_KM <- dd$pWater_KM^2

## ----------- landscape level variables: KM_HF

#' Surrounding footprint at the 1 km$^2$ level
dd$THF_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_HF"]])
dd$Lin_KM <- rowSums(row_std(vc2[rownames(dd),])[,
    tv[colnames(vc2), "is_HF"] & tv[colnames(vc2), "is_linear"]])
## note: no abandoned or rough pasture here
dd$Cult_KM <- rowSums(row_std(vc2[rownames(dd),])[,c("CultivationCrop",
    "CultivationTamePasture", "HighDensityLivestockOperation")])
dd$Alien_KM <- rowSums(row_std(vc2[rownames(dd),])[,
    tv[colnames(vc2), "is_HF"] & tv[colnames(vc2), "is_alien"]])
dd$Nonlin_KM <- dd$THF_KM - dd$Lin_KM
dd$Noncult_KM <- dd$THF_KM - dd$Cult_KM
dd$Succ_KM <- dd$THF_KM - dd$Alien_KM
dd$THF2_KM <- dd$THF_KM^2
dd$Succ2_KM <- dd$Succ_KM^2
dd$Alien2_KM <- dd$Alien_KM^2
dd$Noncult2_KM <- dd$Noncult_KM^2
dd$Nonlin2_KM <- dd$Nonlin_KM^2


## ----------- counts

yy <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/allgrids-bird-counts.csv")

Date <- as.POSIXct(yy$lubridated)

tmp <- as.character(yy$Site)
tmp <- gsub("BG-", "", tmp)
Site <- as.integer(tmp)
tmp <- sapply(strsplit(as.character(yy$StationKey), "-"), "[[", 3)
Station <- as.integer(tmp)
Key <- paste0("BG_", Site, "_", Station)

pp <- data.frame(Key=Key, Site=Site, Station=Station, Date=Date)
compare_sets(rownames(dd), Key)
ii <- pp$Key %in% rownames(dd)
SPP <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATTW", "BAOR", "BARS",
    "BAWW", "BBMA", "BBWA", "BBWO", "BCCH", "BHCO", "BHVI", "BLJA",
    "BLPW", "BOCH", "BRBL", "BRCR", "BTNW", "CAWA", "CCSP", "CEDW",
    "CHSP", "CMWA", "CONW", "CORA", "COYE", "DEJU", "DOWO", "EAKI",
    "EVGR", "FOSP", "GCKI", "GRAJ", "GRCA", "HAWO", "HETH", "HOLA",
    "HOWR", "LCSP", "LEFL", "LISP", "MAWA", "MODO", "MOWA", "NOFL",
    "NOWA", "OCWA", "OSFL", "OVEN", "PAWA", "PHVI", "PISI", "PIWO",
    "PUFI", "RBGR", "RBNU", "RCKI", "REVI", "RUBL", "RWBL", "SAVS",
    "SOSP", "SWSP", "SWTH", "TEWA", "TRES", "VATH", "VEER", "VESP",
    "WAVI", "WBNU", "WCSP", "WETA", "WEWP", "WIWA", "WIWR", "WTSP",
    "WWCR", "YBFL", "YBSA", "YEWA", "YRWA")
yy <- yy[ii,intersect(colnames(yy), SPP)]
pp <- pp[ii,]
rownames(yy) <- NULL
rownames(pp) <- rownames(yy)

## ----------- offsets


if (!requireNamespace("QPAD")) {
  if (!requireNamespace("remotes"))
    install.packages("remotes")
  remotes::install_github("psolymos/QPAD")
}
if (!requireNamespace("sp"))
  install.packages("sp")
if (!requireNamespace("maptools"))
  install.packages("maptools")
if (!requireNamespace("raster"))
  install.packages("raster")
if (!requireNamespace("intrval"))
  install.packages("intrval")

library(QPAD)
library(maptools)
library(intrval)
library(raster)

load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")

od <- setwd("~/repos/recurring/offset")

rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

dtm <- as.POSIXlt(pp$Date)
dur <- 3
dis <- Inf
day <- as.integer(dtm$yday)
hour <- as.numeric(round(dtm$hour + dtm$min/60, 2))

## intersect here
xy3 <- xy2[match(pp$Key, vv$id_final),]
rownames(xy3@coords) <- rownames(yy)

## LCC4 and LCC2
vlcc <- extract(rlcc, xy3)
vlcc[vlcc==18] <- 6

lcclevs <- c("0"="", "1"="Conif", "2"="Conif", "3"="", "4"="",
  "5"="DecidMixed", "6"="DecidMixed", "7"="", "8"="Open", "9"="",
  "10"="Open", "11"="Open", "12"="Open", "13"="Open", "14"="Wet",
  "15"="Open", "16"="Open", "17"="Open", "18"="", "19"="")
lcc4 <- factor(lcclevs[vlcc+1], c("DecidMixed", "Conif", "Open", "Wet"))
lcc2 <- lcc4
levels(lcc2) <- c("Forest", "Forest", "OpenWet", "OpenWet")

## TREE
vtree <- extract(rtree, xy3)
TREE <- vtree / 100
TREE[TREE %)(% c(0, 1)] <- 0

## extract seedgrow value (this is rounded)
d1 <- extract(rd1, xy3)
## UTC offset + 7 makes Alberta 0 (MDT offset)
tz <- extract(rtz, xy3) + 7

## transform the rest
JDAY <- day / 365
TREE <- vtree / 100
MAXDIS <- rep(dis / 100, nrow(yy))
MAXDUR <- rep(dur, nrow(yy))

## sunrise time adjusted by offset
sr <- sunriset(coordinates(xy3),
  as.POSIXct(dtm, tz="America/Edmonton"),
  direction="sunrise", POSIXct.out=FALSE) * 24
TSSR <- round(unname((hour - sr + tz) / 24), 4)

## days since local spring
DSLS <- (day - d1) / 365

off <- matrix(0, nrow(yy), ncol(yy))
dimnames(off) <- dimnames(yy)


for (spp in colnames(off)) {
    cat(spp, "\n")
    flush.console()
    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model (includes DSLS)
    #mi <- bestmodelBAMspecies(spp, type="BIC",
    #    model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
    mi <- bestmodelBAMspecies(spp, type="BIC")
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)

    ## make Xp and Xq
    #' Design matrices for singing rates (`Xp`) and for EDR (`Xq`)
    Xp <- cbind(
      "(Intercept)"=1,
      "TSSR"=TSSR,
      "JDAY"=JDAY,
      "TSSR2"=TSSR^2,
      "JDAY2"=JDAY^2)
    Xq <- cbind("(Intercept)"=1,
      "TREE"=TREE,
      "LCC2OpenWet"=ifelse(lcc4 %in% c("Open", "Wet"), 1, 0),
      "LCC4Conif"=ifelse(lcc4=="Conif", 1, 0),
      "LCC4Open"=ifelse(lcc4=="Open", 1, 0),
      "LCC4Wet"=ifelse(lcc4=="Wet", 1, 0))

    p <- as.numeric(rep(NA, nrow(off)))
    A <- q <- p
    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model
    mi <- bestmodelBAMspecies(spp, type="BIC",
        model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    ## design matrices matching the coefs
    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0
    ## calculate p, q, and A based on constant phi and tau for the respective NAs
    p[!OKp] <- sra_fun(MAXDUR[!OKp], cf0[1])
    unlim <- ifelse(MAXDIS[!OKq] == Inf, TRUE, FALSE)
    A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * MAXDIS[!OKq]^2)
    q[!OKq] <- ifelse(unlim, 1, edr_fun(MAXDIS[!OKq], cf0[2]))
    ## calculate time/lcc varying phi and tau for non-NA cases
    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    p[OKp] <- sra_fun(MAXDUR[OKp], phi1)
    unlim <- ifelse(MAXDIS[OKq] == Inf, TRUE, FALSE)
    A[OKq] <- ifelse(unlim, pi * tau1^2, pi * MAXDIS[OKq]^2)
    q[OKq] <- ifelse(unlim, 1, edr_fun(MAXDIS[OKq], tau1))
    ## log(0) is not a good thing, apply constant instead
    ii <- which(p == 0)
    p[ii] <- sra_fun(MAXDUR[ii], cf0[1])

    off[,spp] <- log(p) + log(A) + log(q)

}

sum(is.na(off))
range(off)

## ----------- checks/saving

en <- new.env()
load("d:/abmi/AB_data_v2018/data/analysis/birds/data/ab-birds-north-2018-12-07.RData", envir=en)
Xn <- get_model_matrix(dd, en$mods)

compare_sets(colnames(SSH), colnames(en$SSH))

names(en)

YY <- yy
OFF <- off
DAT <- data.frame(pp, dd[match(pp$Key, rownames(dd)),])

all(rownames(YY)==rownames(OFF))
all(rownames(YY)==rownames(DAT))

save(DAT, OFF, YY, gg, vv, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-data-package.RData")
