##% Processing data for intactness, all-in-one, website things
##% P Solymos
##% Aug 18, 2015

library(RODBC)
library(mefa4)
library(raster)
library(sp)
library(rgdal)
library(maptools)

## use preprocessed national data set for BAM & BBS

ROOT <- "c:/bam/May2015"
load(file.path(ROOT, "out", "data_package_2015-08-14.Rdata"))
load(file.path(ROOT, "out", "offsets_allspp_BAMBBS_2015-07-24.Rdata"))
colnames(OFF)[colnames(OFF) == "YWAR"] <- "YEWA"
rownames(TAX) <- TAX$Species_ID
DAT <- data.frame(PKEY, SS[match(PKEY$SS, SS$SS),])
DAT$SS.1 <- NULL
DAT$PCODE.1 <- NULL
#rm(PCTBL, PKEY, SS)
table(duplicated(DAT$PKEY))
DAT <- DAT[!is.na(DAT$JURS) & DAT$JURS == "AB",]
table(duplicated(DAT$PKEY))
rownames(DAT) <- DAT$PKEY

## ABMI data
load(file=file.path(ROOT, "out",
    paste0("abmi_data_package_2015-08-18.Rdata")))

YY <- Xtab(ABUND ~ PKEY + SPECIES, PCTBL)
ii <- sort(intersect(rownames(DAT), rownames(YY)))
DAT <- DAT[ii,]
YY <- YY[ii,]
#YY <- YY[,colSums(YY) > 0]

YY2 <- Xtab(ABUND ~ PKEY + SPECIES, pc2)
ii <- sort(intersect(rownames(dat2), rownames(YY2)))
DAT2 <- dat2[ii, ]
YY2 <- YY2[ii,]

ii <- sort(intersect(colnames(YY), colnames(YY2)))
YY <- YY[,ii]
YY2 <- YY2[,ii]

rn <- paste0("T_", rownames(DAT2))
rn <- gsub("_PC_", "_PT_", rn)
rownames(DAT2) <- rn
rownames(YY2) <- rn

## habitat info

## Function to compare sets (factors are left untouched)
compare.sets <- function(x, y) {
    x <- as.factor(x)
    y <- as.factor(y)
    xl <- levels(x)
    yl <- levels(y)
    xa <- levels(droplevels(x))
    ya <- levels(droplevels(y))
    lab <- c(xlength=length(xl), ylength=length(yl),
        intersect=length(intersect(xl, yl)),
        union=length(union(xl, yl)),
        xbutnoty=length(setdiff(xl, yl)),
        ybutnotx=length(setdiff(yl, xl)))
    act <- c(xlength=length(xa), ylength=length(ya),
        intersect=length(intersect(xa, ya)),
        union=length(union(xa, ya)),
        xbutnoty=length(setdiff(xa, ya)),
        ybutnotx=length(setdiff(ya, xa)))
    rbind(labels=lab, unique=act)
}

ROOT2 <- "c:/p/AB_data_v2015"

load(file.path(ROOT2, "out", "abmi_onoff", 
    "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0.Rdata"))
rm(dd1ha)
load(file.path(ROOT2, "out", "bambbs", 
    "veg-hf_bambbs_fix-fire_fix-age0.Rdata"))

compare.sets(rownames(YY), rownames(climPoint_bambbs))
compare.sets(rownames(YY2), rownames(climPoint))

grep("-11_", rownames(climPoint))
rn <- gsub("-11_", "-2_", rownames(climPoint))
compare.sets(rownames(YY2), rownames(climPoint))
compare.sets(rownames(YY2), rn)
compare.sets(rownames(dd150m[[1]]), rownames(climPoint))
compare.sets(rownames(dd1km[[1]]), rownames(climPoint))
rownames(climPoint) <- rn
rownames(dd150m[[1]]) <- rownames(dd150m[[2]]) <- rn
rownames(dd150m[[3]]) <- rownames(dd150m[[4]]) <- rn
rownames(dd1km[[1]]) <- rownames(dd1km[[2]]) <- rn
rownames(dd1km[[3]]) <- rownames(dd1km[[4]]) <- rn

rn <- sort(intersect(rownames(YY), rownames(climPoint_bambbs)))
DAT <- DAT[rn,]
YY <- YY[rn,]
climPoint_bambbs <- climPoint_bambbs[rn,]
for (i in 1:4) {
    dd150m_bambbs[[i]] <- dd150m_bambbs[[i]][rn,]
    dd1km_bambbs[[i]] <- dd1km_bambbs[[i]][rn,]
}

rn <- sort(intersect(rownames(YY2), rownames(climPoint)))
DAT2 <- DAT2[rn,]
YY2 <- YY2[rn,]
climPoint <- climPoint[rn,]
for (i in 1:4) {
    dd150m[[i]] <- dd150m[[i]][rn,]
    dd1km[[i]] <- dd1km[[i]][rn,]
}

tmp <- strsplit(rownames(DAT2), "_")
DAT2$SS <- as.factor(sapply(tmp, function(z) paste("ABMI", z[4], z[8], sep="_")))
DAT2$SITE <- DAT2$SS

## join
YY <- rbind(YY, YY2)
for (i in 1:4) {
    dd150m[[i]] <- rbind(dd150m_bambbs[[i]], dd150m[[i]])
    dd1km[[i]] <- rbind(dd1km_bambbs[[i]], dd1km[[i]])
}
climPoint$YEAR <- climPoint$Year
climPoint$PKEY <- rownames(climPoint)
climPoint_bambbs$POINT_X <- DAT$X[match(rownames(climPoint_bambbs), rownames(DAT))]
climPoint_bambbs$POINT_Y <- DAT$Y[match(rownames(climPoint_bambbs), rownames(DAT))]
cn <- intersect(colnames(climPoint_bambbs),colnames(climPoint))
climPoint <- rbind(climPoint_bambbs[,cn], climPoint[,cn])
cn <- intersect(colnames(DAT),colnames(DAT2))
DAT <- rbind(DAT[,cn], DAT2[,cn])

all(rownames(DAT) == rownames(YY))
all(rownames(DAT) == rownames(climPoint))
all(rownames(DAT) == rownames(dd150m[[1]]))
all(rownames(DAT) == rownames(dd1km[[1]]))

DAT <- data.frame(DAT, climPoint)
DAT$PKEY.1 <- NULL
YY <- YY[,colSums(YY) > 0]
TAX <- droplevels(TAX[colnames(YY),])
## now we have a clean version of: YY, DAT, dd150m, dd1km
DAT$HAB_NALC1 <- NULL
DAT$HAB_NALC2 <- NULL
DAT$TREE3 <- NULL
## imputing
DAT$JDAY[is.na(DAT$JDAY)] <- mean(DAT$JDAY, na.rm=TRUE)
DAT$TSSR[is.na(DAT$TSSR)] <- mean(DAT$TSSR, na.rm=TRUE)
data.frame(x=colSums(is.na(DAT)))

DAT <- droplevels(DAT)

## lookup tables & reclassed tables

ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")
ts$UseInAnalysis <- as.character(ts$UseInAnalysis)
ts$UseInAnalysis[is.na(ts$UseInAnalysis)] <- as.character(ts$Levels4)[is.na(ts$UseInAnalysis)]
ts$Levels1 <- as.character(ts$Levels1)
ts$Levels1[is.na(ts$Levels1)] <- ts$UseInAnalysis[is.na(ts$Levels1)]

SoilKmCr <- groupSums(dd1km$soil_current, 2, ts[colnames(dd1km$soil_current), "UseInAnalysis"])
SoilKmRf <- groupSums(dd1km$soil_reference, 2, ts[colnames(dd1km$soil_reference), "UseInAnalysis"])
SoilPcCr <- groupSums(dd150m$soil_current, 2, ts[colnames(dd150m$soil_current), "UseInAnalysis"])
SoilPcRf <- groupSums(dd150m$soil_reference, 2, ts[colnames(dd150m$soil_reference), "UseInAnalysis"])

psoilhf <- groupSums(dd150m$soil_current, 2, ts[colnames(dd150m$soil_current), "Levels1"])
psoilhf <- psoilhf[,!(colnames(psoilhf) %in% c("SoilUnknown", "SoilWater",
    "SoilWetland", "HWater", "HFor", "SoftLin", "HardLin"))]
psoilhf <- as.matrix(psoilhf / rowSums(psoilhf))

#SoilKmCr <- SoilKmCr[,colnames(SoilKmCr) != "HFor"] # exclude forestry
SoilKmCr <- as.matrix(SoilKmCr / rowSums(SoilKmCr))
SoilKmRf <- as.matrix(SoilKmRf / rowSums(SoilKmRf))
#SoilPcCr <- SoilPcCr[,colnames(SoilPcCr) != "HFor"] # exclude forestry
SoilPcCr <- as.matrix(SoilPcCr / rowSums(SoilPcCr))
SoilPcRf <- as.matrix(SoilPcRf / rowSums(SoilPcRf))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv$UseInAnalysis <- as.character(tv$UseInAnalysis)
ii <- is.na(tv$UseInAnalysis) | tv$UseInAnalysis == "HFor"
tv$UseInAnalysis[ii] <- as.character(tv$VEGAGE_use)[ii]
tv$UseInAnalysis[tv$UseInAnalysis == "WetBare"] <- "NonVeg"
tv$UseInAnalysis[tv$UseInAnalysis %in% c("WetGrassHerb", "WetShrub")] <- "Wetland"

VegKmCr <- groupSums(dd1km$veg_current, 2, tv[colnames(dd1km$veg_current), "UseInAnalysis"])
VegKmRf <- groupSums(dd1km$veg_reference, 2, tv[colnames(dd1km$veg_reference), "UseInAnalysis"])
VegPcCr <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "UseInAnalysis"])
VegPcRf <- groupSums(dd150m$veg_reference, 2, tv[colnames(dd150m$veg_reference), "UseInAnalysis"])

VegKmCr <- as.matrix(VegKmCr / rowSums(VegKmCr))
VegKmRf <- as.matrix(VegKmRf / rowSums(VegKmRf))
VegPcCr <- as.matrix(VegPcCr / rowSums(VegPcCr))
VegPcRf <- as.matrix(VegPcRf / rowSums(VegPcRf))

find_max <- function(x) {
    tmp <- apply(x, 1, which.max)
    tmp <- factor(tmp, levels=seq_len(ncol(x)))
    levels(tmp) <- colnames(x)
    tmp
}

## hab & CC in north

tv$UseInAnalysis2 <- as.character(tv$Type)
ii <- is.na(tv$UseInAnalysis2)
tv$UseInAnalysis2[ii] <- as.character(tv$UseInAnalysis)[ii]
tv$UseInAnalysis2[tv$UseInAnalysis2 == "WetBare"] <- "NonVeg"
tv$UseInAnalysis2[tv$UseInAnalysis2 %in% c("WetGrassHerb", "WetShrub")] <- "Wetland"
tv$UseInAnalysis3 <- tv$UseInAnalysis2
ii <- !is.na(tv$HF) & tv$HF == "CutBlocks"
tv$UseInAnalysis3[ii] <- "CC" #paste0("CC", tv$UseInAnalysis3[ii])
tv$UseAge <- as.character(tv$AGE)
tv$UseAge[is.na(tv$UseAge)] <- ""
tv$EC_AGE <- as.character(tv$EC_AGE)
tv$EC_AGE[is.na(tv$EC_AGE)] <- ""
tv$LCC5 <- tv$EC_AGE
tv$LCC5[tv$EC_AGE == ""] <- "5"
tv$LCC5[tv$EC_AGE %in% c("A","B","R") & tv$UseInAnalysis2 %in% 
    c("Decid","Mixwood")] <- "4"
tv$LCC5[tv$EC_AGE %in% c("A","B","R") & tv$UseInAnalysis2 %in% 
    c("Conif","Pine","BSpr", "Larch")] <- "3"
tv$LCC5[tv$EC_AGE %in% c("C","D") & tv$UseInAnalysis2 %in% 
    c("Decid","Mixwood")] <- "2"
tv$LCC5[tv$EC_AGE %in% c("C","D") & tv$UseInAnalysis2 %in% 
    c("Conif","Pine","BSpr", "Larch")] <- "1"


tmp <- groupSums(dd150m$veg_reference, 2, tv[colnames(dd150m$veg_reference), "UseInAnalysis2"])
tmp <- tmp[,!(colnames(tmp) %in% c("NonVeg", "Water"))]
tmp <- as.matrix(tmp / rowSums(tmp))
DAT$hab0 <- find_max(tmp)

tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "UseInAnalysis2"])
pveghf <- tmp[,!(colnames(tmp) %in% c("NonVeg", "Water","HWater"))]
pveghf <- as.matrix(pveghf / rowSums(pveghf))

tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "UseInAnalysis3"])
tmp <- tmp[,!(colnames(tmp) %in% c("NonVeg", "Water","HWater","SoftLin","HardLin"))]
tmp <- as.matrix(tmp / rowSums(tmp))
DAT$hab1cc <- find_max(tmp)

DAT$hab1 <- as.character(DAT$hab1cc)
ii <- !is.na(DAT$hab1cc) & !is.na(DAT$hab0) &
    DAT$hab1cc == "CC" & DAT$hab0 %in% c("Conif","Decid","Mixwood","Pine")
DAT$hab1[ii] <- as.character(DAT$hab0[ii])
DAT$isCC <- ifelse(ii, 1L, 0L)
ii <- !is.na(DAT$hab1cc) & !is.na(DAT$hab0) & DAT$hab1 == "CC"
DAT$hab1[ii] <- as.character(DAT$hab0[ii])

DAT$hab0 <- relevel(DAT$hab0, "Decid")
DAT$hab1 <- relevel(as.factor(DAT$hab1), "Decid")
DAT$hab1cc <- relevel(DAT$hab1cc, "Decid")

table(DAT$hab1, DAT$hab1cc, useNA="a")
table(DAT$hab1,DAT$isCC)

## age 

AgePtCr <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "UseAge"])
AgePtCr <- as.matrix(AgePtCr / rowSums(AgePtCr))
## exclude unknown (0) and non-forest (blank)
AgePtCr <- t(AgePtCr[,c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")])
AgeMin <- structure(c(0,10,20,40,60,80,100,120,140,160)/200,
    names=c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
DAT$wtAge <- colSums(AgePtCr * AgeMin) / colSums(AgePtCr)
DAT$wtAge[is.na(DAT$wtAge)] <- 0
FORclasses <- c("Conif","Decid","Mixwood","Pine","BSpr", "Larch")
DAT$wtAge[!(DAT$hab1 %in% FORclasses)] <- 0
## recent fire (non CC)
DAT$isRR <- ifelse(DAT$wtAge < 0.1 & !DAT$isCC, 1L, 0L)
DAT$isFor <- DAT$hab1 %in% FORclasses
DAT$isConif <- DAT$hab1 %in% c("Conif","Pine","BSpr", "Larch")
DAT$isRR[!DAT$isFor] <- 0
DAT$hab1b <- DAT$hab1
levels(DAT$hab1b) <- c(levels(DAT$hab1b), "Burn")
DAT$hab1b[DAT$isRR == 1] <- "Burn"

table(burn=DAT$isRR, cc=DAT$isCC)
table(DAT$hab1b,DAT$hab1)

## forestry convergence
## fCC1: linear
MAXFOR <- 50/200
DAT$fCC1 <- 0
DAT$fCC1[DAT$isCC==1] <- pmax(0, 1 - (DAT$isCC * DAT$wtAge/MAXFOR)[DAT$isCC==1])
plot(fCC1 ~ wtAge, DAT[DAT$isCC==1,])
## Dave`s recovery trajectories
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3, 
    82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85, 
    88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100
DAT$fCC2 <- 0
tmp1 <- approxfun(age, decid)(DAT$wtAge)
tmp1[is.na(tmp1)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- DAT$isFor & !DAT$isConif & DAT$isCC==1
DAT$fCC2[ii] <- tmp1[ii]
tmp2 <- approxfun(age, conif)(DAT$wtAge)
tmp2[is.na(tmp2)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- DAT$isConif & DAT$isCC==1
DAT$fCC2[ii] <- tmp2[ii]
plot(DAT$wtAge, DAT$fCC1, col=2, pch=".")
points(DAT$wtAge, DAT$fCC2, col=4, pch=".")
sum(is.na(DAT$fCC2))

table(DAT$hab1,round(DAT$wtAge,1),useNA="a")
table(DAT$hab1,DAT$isCC,useNA="a")
table(DAT$hab1,DAT$isRR,useNA="a")
table(DAT$isRR,DAT$isCC,useNA="a")

tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "LCC5"])
tmp <- as.matrix(tmp / rowSums(tmp))
DAT$hab_lcc <- find_max(tmp)
DAT$hab_lcc <- as.integer(as.character(DAT$hab_lcc))
DAT$hab_lcc <- factor(DAT$hab_lcc, 1:5)

DAT$ECage <- factor("", levels=c("", "R","A", "B", "C", "D"))
DAT$ECage[DAT$isFor & DAT$wtAge >= 0] <- "R"
DAT$ECage[DAT$wtAge >= 10/200] <- "A"
DAT$ECage[DAT$wtAge >= 20/200] <- "B"
DAT$ECage[DAT$wtAge >= 40/200] <- "C"
DAT$ECage[DAT$isFor & !DAT$isConif & DAT$wtAge >= 60/200] <- "D" # decid
DAT$ECage[DAT$isConif & DAT$wtAge >= 80/200] <- "D" # conif
DAT$hab1ec <- interaction(DAT$hab1, DAT$ECage, drop=TRUE, sep="")
DAT$hab1ec <- relevel(DAT$hab1ec, "DecidD")

## LCC_combo & TREE: use hab info

tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "EC_AGE"])
tmp <- as.matrix(tmp / rowSums(tmp))
DAT$ClosedCanopy <- tmp[,"C"] + tmp[,"D"]

DAT$TREE[is.na(DAT$TREE)] <- DAT$ClosedCanopy[is.na(DAT$TREE)]
DAT$TREE[is.na(DAT$TREE)] <- DAT$ClosedCanopy[is.na(DAT$TREE)]

DAT$LCC_combo[is.na(DAT$LCC_combo)] <- DAT$hab_lcc[is.na(DAT$LCC_combo)] 

## soil in south

round(100*colMeans(SoilPcRf[SoilPcRf[,"SoilUnknown"]==0,]),2)
round(100*colMeans(SoilPcCr[SoilPcRf[,"SoilUnknown"]==0,]),2)

tmp <- SoilPcRf[,!(colnames(SoilPcRf) %in% c("SoilWater","SoilWetland"))]
DAT$soil0 <- find_max(tmp / rowSums(tmp))
DAT$soil0[DAT$soil0 == "SoilUnknown"] <- NA
DAT$soil0 <- droplevels(DAT$soil0)
DAT$soil0 <- relevel(DAT$soil0, "Productive")
table(DAT$soil0, useNA="a")

tmp <- SoilPcCr[,!(colnames(SoilPcCr) %in% c("SoilWater","SoilWetland",
    "HWater","SoftLin", "HardLin", "HFor"))]
DAT$soil1 <- find_max(tmp / rowSums(tmp))
DAT$soil1[DAT$soil1 == "SoilUnknown"] <- NA
DAT$soil1 <- droplevels(DAT$soil1)
DAT$soil1 <- relevel(DAT$soil1, "Productive")
table(DAT$soil1, useNA="a")

table(DAT$soil1, DAT$soil0, useNA="a")

DAT$soil1v <- DAT$soil1
levels(DAT$soil1v)[levels(DAT$soil1) %in% c("RapidDrain","SalineAndClay")] <- "NonProductive"
table(DAT$soil1, DAT$soil1v, useNA="a")

## Water

DAT$pWater <- rowSums(dd150m$veg_current[,c("Water",
    "BorrowpitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals")]) /
    rowSums(dd150m$veg_current)

## ROAD and sof lin PC

DAT$Road_PC <- rowSums(dd150m$veg_current[,c("RailHardSurface",
    "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / rowSums(dd150m$veg_current)
DAT$ROAD01 <- ifelse(DAT$Road_PC >= 0.1, 1L, 0L)
DAT$ROAD01[DAT$PCODE == "BBSAB"] <- 1L
DAT$SoftLin_PC <- rowSums(dd150m$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge")]) /
    rowSums(dd150m$veg_current)

## ARU

DAT$ARU <- ifelse(DAT$PCODE %in% c("EMCLA","EMCLA2014"), 1L, 0L)

## transformations

## mapping projection
XYlatlon <- DAT[,c("POINT_X", "POINT_Y")]
coordinates(XYlatlon) <- ~ POINT_X + POINT_Y
proj4string(XYlatlon) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
XY <- as.data.frame(spTransform(XYlatlon, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))
DAT$X <- XY[,"POINT_X"]
DAT$Y <- XY[,"POINT_Y"]

transform_CLIM <- function(x, ID="PKEY") {
    z <- x[,ID,drop=FALSE]
    z$xlong <- (x$POINT_X - (-113.7)) / 2.15
    z$xlat <- (x$POINT_Y - 53.8) / 2.28
    z$xAHM <- (x$AHM - 0) / 50
    z$xPET <- (x$PET - 0) / 800
    z$xFFP <- (x$FFP - 0) / 130
    z$xMAP <- (x$MAP - 0) / 2200
    z$xMAT <- (x$MAT - 0) / 6
    z$xMCMT <- (x$MCMT - 0) / 25
    z$xMWMT <- (x$MWMT - 0) / 20
    z
}
DAT$MAP <- gsub(",", "", DAT$MAP)
DAT$MAP <- as.numeric(DAT$MAP)
DAT <- data.frame(DAT, transform_CLIM(DAT))
DAT$PKEY.1 <- NULL


## surrounding HF at 1km scale

DAT$THF_KM <- rowSums(VegKmCr[,setdiff(colnames(VegKmCr), colnames(VegKmRf))])
DAT$Lin_KM <- rowSums(dd1km$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailHardSurface", "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / rowSums(dd1km$veg_current)
DAT$Nonlin_KM <- DAT$THF_KM - DAT$Lin_KM
DAT$Cult_KM <- rowSums(dd1km$veg_current[,c("CultivationCropPastureBareground",
    "HighDensityLivestockOperation")]) / rowSums(dd1km$veg_current)
DAT$Noncult_KM <- DAT$THF_KM - DAT$Cult_KM
CClabs <- colnames(dd1km$veg_current)[grep("CC", colnames(dd1km$veg_current))]
DAT$Succ_KM <- rowSums(dd1km$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge", 
    CClabs)]) / rowSums(dd1km$veg_current)
DAT$Alien_KM <- DAT$THF_KM - DAT$Succ_KM

DAT$THF2_KM <- DAT$THF_KM^2
DAT$Succ2_KM <- DAT$Succ_KM^2
DAT$Alien2_KM <- DAT$Alien_KM^2
DAT$Noncult2_KM <- DAT$Noncult_KM^2
DAT$Nonlin2_KM <- DAT$Nonlin_KM^2

## restricted data use

PCODE_useOK <- c("ABCAWAWEST", "MGLE", "FTL", "THIN", 
    "LMWELL", "EMB-ASP", "EMB-BS", "EMB-NOISE", 
    "ECJOSM", "ECJOSM_JRB", 
    "EMCLA", "EMCLA2014", "ABMI", "BBSAB")
DAT$useOK <- DAT$PCODE %in% PCODE_useOK
DAT$useOK[DAT$YEAR > 2007 & DAT$PCODE == "CL"] <- TRUE # Calling Lake


## YR

DAT$YR <- (DAT$YEAR - 1997) / 10
DAT$YR5 <- cut(DAT$YEAR, c(1996, 2001, 2006, 2010, 2015))
table(DAT$YEAR, DAT$YR5)
table(DAT$YR5)

## Bootstrap blocks

## NR and lat/long
DAT$bootid <- as.character(DAT$NRNAME)
tmp <- as.character(interaction(ifelse(DAT$POINT_Y < (56.5), "S", "N"),
    ifelse(DAT$POINT_X < (-115.5), "W", "E"), sep=""))
DAT$bootid[DAT$bootid %in% c("Boreal","Canadian Shield")] <-
    tmp[DAT$bootid %in% c("Boreal","Canadian Shield")]
#with(DAT, plot(X, Y, col=as.factor(bootid)))
table(DAT$bootid, DAT$YR5)
DAT$bootid <- interaction(DAT$bootid, as.integer(DAT$YR5))

## use in south


DAT$useSouth <- FALSE
DAT$useSouth[DAT$NRNAME %in% c("Grassland", "Parkland")] <- TRUE
DAT$useSouth[DAT$NSRNAME %in% c("Dry Mixedwood")] <- TRUE
DAT$useSouth[DAT$useSouth & DAT$POINT_Y > 56.7] <- FALSE
DAT$useSouth[DAT$useSouth & DAT$pWater > 0.5] <- FALSE
DAT$useSouth[DAT$useSouth & SoilPcRf[,"SoilUnknown"] > 0] <- FALSE
DAT$useSouth[DAT$useSouth & is.na(DAT$soil1)] <- FALSE

## use in north
DAT$useNorth <- DAT$NRNAME != "Grassland"

## within year visits

DAT <- DAT[sample.int(nrow(DAT), nrow(DAT)),]
DAT$SS_YR <- interaction(DAT$SS, DAT$YEAR, drop=TRUE)
table(table(DAT$SS_YR))
dup <- unique(DAT$SS_YR[duplicated(DAT$SS_YR)])
tmp1 <- DAT[DAT$SS_YR %in% dup, c("PKEY","SS_YR")]
rownames(tmp1) <- tmp1$PKEY
dim(tmp1)
tmp2 <- nonDuplicated(tmp1[sample.int(nrow(tmp1)),], SS_YR)
rownames(tmp2) <- tmp2$PKEY
dim(tmp2)
dim(tmp1) - dim(tmp2)
tmp1 <- tmp1[setdiff(rownames(tmp1), rownames(tmp2)),]
dim(tmp1)
dim(tmp1) + dim(tmp2)
DAT$Revisit <- DAT$PKEY %in% rownames(tmp1)
table(DAT$Revisit)
table(DAT$PCODE,DAT$Revisit)
DAT$SS_YR <- NULL

DAT$SITE <- as.character(DAT$SITE)
DAT$SITE[is.na(DAT$SITE)] <- as.character(DAT$SS)[is.na(DAT$SITE)]
DAT$SITE <- as.factor(DAT$SITE)

## derived variables

DAT$xlat2 <- DAT$xlat^2
DAT$xlong2 <- DAT$xlong^2
DAT$wtAge2 <- DAT$wtAge^2
DAT$wtAge05 <- DAT$wtAge^0.5

DAT$hab_lcc3 <- DAT$hab_lcc
levels(DAT$hab_lcc3) <- c("1", "1", "2", "2", "3")
DAT$hab_lcc2 <- DAT$hab_lcc
levels(DAT$hab_lcc2) <- c("1", "1", "1", "1", "2")

## simply treat Mixed as intercept with Decid as reference
DAT$isMix <- ifelse(DAT$hab1 == "Mixwood", 1L, 0L)
DAT$isWSpruce <- ifelse(DAT$hab1 == "Conif", 1L, 0L)
DAT$isPine <- ifelse(DAT$hab1 == "Pine", 1L, 0L)
DAT$isBSpruce <- ifelse(DAT$hab1 == "BSpr", 1L, 0L)
DAT$isLarch <- ifelse(DAT$hab1 == "Larch", 1L, 0L)
DAT$isBSLarch <- ifelse(DAT$hab1 %in% c("BSpr", "Larch"), 1L, 0L)
DAT$isUpCon <- ifelse(DAT$hab1 %in% c("Conif", "Pine"), 1L, 0L)
DAT$isCon <- ifelse(DAT$hab1 %in% c("BSpr", "Larch",
    "Conif", "Pine"), 1L, 0L)

source("~/repos/abmianalytics/R/analysis_models.R")
source("~/repos/bragging/R/glm_skeleton.R")
compare.sets(getTerms(modsSoil, "list"), colnames(DAT))
setdiff(getTerms(modsSoil, "list"), colnames(DAT))
compare.sets(getTerms(modsVeg, "list"), colnames(DAT))
setdiff(getTerms(modsVeg, "list"), colnames(DAT))




## offsets

offdat <- DAT[,c("JDAY","TSSR","TREE","LCC_combo","MAXDUR","MAXDIS")]

library(detect)
load_BAM_QPAD(version=1)
BAMspp <- getBAMspecieslist()
load("~/Dropbox/abmi/intactness/dataproc/BAMCOEFS25.Rdata")
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

(sppp <- union(BAMspp, BAMCOEFS25$spp))

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- rownames(offdat)
colnames(OFF) <- sppp
for (i in sppp) {
    cat(i, date(), "\n");flush.console()
    tmp <- try(offset_fun(j=1, i, offdat))
    if (!inherits(tmp, "try-error"))
        OFF[,i] <- tmp
}
## 99-100 percentile can be crazy high (~10^5), thus reset
for (i in sppp) {
    q <- quantile(OFF[,i], 0.99, na.rm=TRUE)
    OFF[!is.na(OFF[,i]) & OFF[,i] > q, i] <- q
}
colSums(is.na(OFF))/nrow(OFF)
apply(exp(OFF), 2, range, na.rm=TRUE)

OFFmean <- log(rowMeans(exp(OFF)))

compare.sets(rownames(OFF),rownames(DAT))

## subsets

keep <- DAT$YEAR >= 1997 & !is.na(DAT$hab1) & !DAT$Revisit

DAT$keep <- keep
YY <- YY[rownames(DAT),]
pveghf <- pveghf[rownames(DAT),]
psoilhf <- psoilhf[rownames(DAT),]
save(DAT, YY, OFF, OFFmean, TAX, pveghf, psoilhf,
    file=file.path(ROOT2, "out", "birds", "data", "data-full-withrevisit.Rdata"))

DAT <- droplevels(DAT[keep,])
YY <- YY[rownames(DAT),]

plot(DAT$X, DAT$Y, col=ifelse(DAT$useOK, 1, 2), pch=19, cex=0.2)

OFF <- OFF[rownames(DAT),]
OFFmean <- OFFmean[rownames(DAT)]

compare.sets(rownames(OFF),rownames(DAT))
pveghf <- pveghf[rownames(DAT),]
#save(DAT, YY, OFF, OFFmean, TAX, pveghf,
#    file=file.path(ROOT2, "out", "birds", "data", "data-full.Rdata"))


DATSfull <- DAT[DAT$useSouth,]
DATNfull <- DAT[DAT$useNorth,]

aaS <- colSums(is.na(DATSfull))
aaN <- colSums(is.na(DATNfull))
aaS[aaS>0]
aaN[aaN>0]

DATS <- DATSfull[DATSfull$useOK,]
DATN <- DATNfull[DATNfull$useOK,]

sapply(list(DATSfull,DATNfull,DATS,DATN), nrow)

YYS <- YY[rownames(DATS),]
YYN <- YY[rownames(DATN),]
YYSfull <- YY[rownames(DATSfull),]
YYNfull <- YY[rownames(DATNfull),]

## bootids

#DAT1 <- DAT
source("~/repos/detect/R/hbootindex.R")
B <- 239

bbfun <- function(DAT1, B) {
    set.seed(1234)
    ## make sure that time intervals are considered as blocks
    ## keep out 10% of the data for validation
    id2 <- list()
    for (l in levels(DAT1$bootid)) {
        sset <- which(DAT1$bootid == l)
        id2[[l]] <- sample(sset, floor(length(sset) * 0.9), FALSE)
    }
    KEEP_ID <- unname(unlist(id2))
    HOLDOUT_ID <- setdiff(seq_len(nrow(DAT1)), KEEP_ID)

    DAT1k <- DAT1[KEEP_ID,]
    DAT1 <- DAT1[c(KEEP_ID, HOLDOUT_ID),]
    DAT1k$SITE <- droplevels(DAT1k$SITE)
    BB1 <- hbootindex(DAT1k$SITE, DAT1k$bootid, B=B)
    BB1
}
BBS <- bbfun(DATS, B)
BBN <- bbfun(DATN, B)
#BBSfull <- bbfun(DATSfull, B)
#BBNfull <- bbfun(DATNfull, B)

## figure out sets of species to analyze

## pa maps and habitat suitability: DAT[useOK,]
## S/N: nmin=25
## look at taxonomy???


OFF0 <- OFF
OFFmean0 <- OFFmean
nmin <- 25

DAT <- DATS
YY <- YYS
YY <- YY[,colSums(YY>0) >= nmin]
BB <- BBS
OFF <- OFF0[rownames(DAT),colnames(OFF0) %in% colnames(YY)]
OFFmean <- OFFmean0[rownames(DAT)]
mods <- modsSoil
save(DAT, YY, OFF, OFFmean, mods, BB,
    file=file.path(ROOT2, "out", "birds", "data", "data-useok-south.Rdata"))

DAT <- DATN
YY <- YYN[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= nmin]
BB <- BBN
OFF <- OFF0[rownames(DAT),colnames(OFF0) %in% colnames(YY)]
OFFmean <- OFFmean0[rownames(DAT)]
mods <- modsVeg
save(DAT, YY, OFF, OFFmean, mods, BB,
    file=file.path(ROOT2, "out", "birds", "data", "data-useok-north.Rdata"))

