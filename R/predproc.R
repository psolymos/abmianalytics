library(mefa4)
#library(raster)
#library(sp)
#library(rgdal)
#library(maptools)

ROOT <- "c:/p/AB_data_v2015"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))

## climate
transform_CLIM <- function(x, ID="Row_Col") {
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
kgrid <- data.frame(kgrid, transform_CLIM(kgrid))
kgrid$xlong2 <- kgrid$xlong^2
kgrid$xlat2 <- kgrid$xlat^2

kgrid$Row_Col.1 <- NULL
kgrid$OBJECTID <- NULL
kgrid$Row <- NULL
kgrid$Col <- NULL
kgrid$AHM <- NULL
kgrid$PET <- NULL
kgrid$FFP <- NULL
kgrid$MAP <- NULL
kgrid$MAT <- NULL
kgrid$MCMT <- NULL
kgrid$MWMT <- NULL
kgrid$Row10 <- NULL
kgrid$Col10 <- NULL
all(rownames(kgrid) == rownames(dd1km_pred$veg_current))

## surrounding HF at 1km scale
kgrid$THF_KM <- rowSums(dd1km_pred$veg_current[,setdiff(colnames(dd1km_pred$veg_current), 
    colnames(dd1km_pred$veg_reference))])
kgrid$Lin_KM <- rowSums(dd1km_pred$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailHardSurface", "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / rowSums(dd1km_pred$veg_current)
kgrid$Nonlin_KM <- kgrid$THF_KM - kgrid$Lin_KM
kgrid$Cult_KM <- rowSums(dd1km_pred$veg_current[,c("CultivationCropPastureBareground",
    "HighDensityLivestockOperation")]) / rowSums(dd1km_pred$veg_current)
kgrid$Noncult_KM <- kgrid$THF_KM - kgrid$Cult_KM
CClabs <- colnames(dd1km_pred$veg_current)[grep("CC", colnames(dd1km_pred$veg_current))]
kgrid$Succ_KM <- rowSums(dd1km_pred$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge", 
    CClabs)]) / rowSums(dd1km_pred$veg_current)
kgrid$Alien_KM <- kgrid$THF_KM - kgrid$Succ_KM

kgrid$THF2_KM <- kgrid$THF_KM^2
kgrid$Succ2_KM <- kgrid$Succ_KM^2
kgrid$Alien2_KM <- kgrid$Alien_KM^2
kgrid$Noncult2_KM <- kgrid$Noncult_KM^2
kgrid$Nonlin2_KM <- kgrid$Nonlin_KM^2
rm(dd1km_pred)

## design matrices

## 0 out in North
cnn0 <- c("ROAD01", "SoftLin_PC", "ARU", "YR", 
    "hab_lcc22:ROAD01", "hab_lcc32:ROAD01", "hab_lcc33:ROAD01", 
    "hab_lcc2:ROAD01", "hab_lcc3:ROAD01", "hab_lcc4:ROAD01", "hab_lcc5:ROAD01")
## habitat in North
cnnHab <- c("(Intercept)", "hab1BSpr", "hab1Conif", "hab1Cult", "hab1GrassHerb", 
    "hab1Larch", "hab1Mixwood", "hab1Pine", "hab1Shrub", "hab1UrbInd", 
    "hab1Wetland", "hab1bBSpr", "hab1bConif", "hab1bCult", "hab1bGrassHerb", 
    "hab1bLarch", "hab1bMixwood", "hab1bPine", "hab1bShrub", "hab1bUrbInd", 
    "hab1bWetland", "hab1bBurn", "wtAge", "wtAge2", "wtAge05", "fCC2", 
    "isCon:wtAge", "isCon:wtAge2", "isUpCon:wtAge", 
    "isBSLarch:wtAge", "isUpCon:wtAge2", "isBSLarch:wtAge2", "isMix:wtAge", 
    "isPine:wtAge", "isWSpruce:wtAge", "isMix:wtAge2", "isPine:wtAge2", 
    "isWSpruce:wtAge2", "isCon:wtAge05", "isUpCon:wtAge05", "isBSLarch:wtAge05", 
    "isMix:wtAge05", "isPine:wtAge05", "isWSpruce:wtAge05")
## 0 out in South
cns0 <- c("pAspen", "ROAD01", "SoftLin_PC", "YR", "hab_lcc22:ROAD01")
## habitat in South
cnsHab <- c("(Intercept)", "soil1RapidDrain", "soil1SalineAndClay", "soil1Cult", "soil1UrbInd", 
    "soil1vNonProductive", "soil1vCult", "soil1vUrbInd")
## climate (North & South)
cnClim <- c("xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT", "xMCMT", 
    "xlat", "xlong", "xlat2", "xlong2", 
    "THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", 
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", "Noncult2_KM", 
    "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xlat:xlong")

## model matrix for Clim & SurroundingHF
fclim <- as.formula(paste("~", paste(cnClim, collapse=" + ")))
#Xclim <- model.matrix(fclim, kgrid)


## example for structure (trSoil, trVeg)
load(file.path(ROOT, "out", "transitions", "LowerPeace_LowerBorealHighlands.Rdata"))







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

if (FALSE) {
## Dave`s recovery trajectories
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3, 
    82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85, 
    88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100
Age <- seq(0, 100, by=0.5) / 200
plot(Age*200, approxfun(age, decid)(Age), type="l", col=2, lwd=2, 
    xlab="Years since disturbance", ylab="Transformed value")
lines(Age*200, approxfun(age, conif)(Age), col=4, lwd=2)
lines(Age*200, pmax(0, 1 - (Age * 200 / 50)), col=3, lwd=2)
lines(Age*200, 1/((Age*200+1)^(1/3)), col=1, lwd=2)
legend("topright", col=1:4, lty=1, lwd=2, 
    legend=c("cubic-root","expert decid","linear","expert conif"))
}

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

