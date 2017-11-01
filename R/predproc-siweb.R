library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"
OUTDIR <- "e:/peter/AB_data_v2016/out/birds/web"

do1 <- FALSE#TRUE # do only 1st run (100%)
doB <- FALSE # do bootstrap (10%)
do10 <- TRUE#FALSE # aggregated info for bootstrap (1%)
## seismic lines can be treated as early seral when no surrounding effects considered
## or apply the backfilled value when surrounding hf is considered
SEISMIC_AS_EARLY_SERAL <- TRUE

PROP <- 100
BMAX <- 100
if (do10 && doB)
    doB <- FALSE

## set here if shf is wanted
load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata")) # kgrid
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

if (FALSE) {
## define 10km level predictions for alternative bootstrap
library(raster)
rt <- raster(file.path("e:/peter/AB_data_v2016", "data", "kgrid", "AHM1k.asc"))
projection(rt) <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
mat0 <- as.matrix(rt)
mat1 <- as.matrix(Xtab(1:nrow(kgrid) ~ Row + Col, kgrid))
mat1[is.na(mat0)] <- NA
r1 <- raster(x=mat1, template=rt)
r10 <- aggregate(r1, fact=10, fun=min)
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- c("POINT_X", "POINT_Y")
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, projection(rt))
matx <- as.matrix(r10)
matx <- t(matx)
matx[!is.na(matx)] <- seq_len(sum(!is.na(matx)))
matx <- t(matx)
rr <- raster(x=matx, template=r10)
kgrid$id10 <- extract(rr, xy)
na <- is.na(kgrid$id10)
#tb <- table(kgrid$id10)
#tb100 <- as.integer(names(tb))[tb == 100]
for (i in which(na)) {
    xy0 <- coordinates(xy)[i,]
    d <- sqrt((xy0[1]-coordinates(xy)[,1])^2+(xy0[2]-coordinates(xy)[,2])^2)
    d[na] <- Inf
#    d[kgrid$id10 %in% tb100] <- Inf
    kgrid$id10[i] <- kgrid$id10[which.min(d)]
}
table(table(kgrid$id10))
}

## climate
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

    z$xASP <- x$ASP
    z$xSLP <- log(x$SLP + 1)
    z$xTRI <- log(x$TRI / 5)
    z$xCTI <- log((x$CTI + 1) / 10)
    z
}
kgrid <- data.frame(kgrid, transform_CLIM(kgrid, "Row_Col"))
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
kgrid$ASP <- NULL
kgrid$TRI <- NULL
kgrid$SLP <- NULL
kgrid$CTI <- NULL
all(rownames(kgrid) == rownames(dd1km_pred$veg_current))

## surrounding HF at 1km scale
kgrid$THF_KM <- rowSums(dd1km_pred$veg_current[,setdiff(colnames(dd1km_pred$veg_current),
    colnames(dd1km_pred$veg_reference))]) / rowSums(dd1km_pred$veg_current)
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

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
kgrid$WetKM <- rowSums(dd1km_pred$veg_current[,tv[colnames(dd1km_pred$veg_current), "WET"]==1]) / rowSums(dd1km_pred$veg_current)
#kgrid$WaterKM <- rowSums(dd1km_pred$veg_current[,tv[colnames(dd1km_pred$veg_current), "WATER"]==1]) / rowSums(dd1km_pred$veg_current)
kgrid$WetWaterKM <- rowSums(dd1km_pred$veg_current[,tv[colnames(dd1km_pred$veg_current), "WETWATER"]==1]) / rowSums(dd1km_pred$veg_current)
kgrid$WetKM0 <- rowSums(dd1km_pred$veg_reference[,tv[colnames(dd1km_pred$veg_reference), "WET"]==1]) / rowSums(dd1km_pred$veg_reference)
kgrid$WetWaterKM0 <- rowSums(dd1km_pred$veg_reference[,tv[colnames(dd1km_pred$veg_reference), "WETWATER"]==1]) / rowSums(dd1km_pred$veg_reference)

## 10x10km units
kgrid0 <- NULL
regs <- levels(kgrid$LUFxNSR)
if (do10) {
    kgrid0 <- kgrid
    kgrid <- as.matrix(kgrid[,sapply(kgrid, is.numeric)])
    #kgrid <- data.frame(groupMeans(kgrid, 1, kgrid0$id10))
    kgrid <- data.frame(groupMeans(kgrid, 1, kgrid0$Row10_Col10))
}

rm(dd1km_pred)

## design matrices

## 0 out in North
cnn0 <- c("ROAD01", "SoftLin_PC", "ARU2ARU", "ARU3SM", "ARU3RF", "YR", "habCl:ROAD01")
## habitat in North
cnnHab <- c("(Intercept)", "hab1BSpr", "hab1Conif", "hab1Cult", "hab1GrassHerb",
    "hab1Larch", "hab1Mixwood", "hab1Pine", "hab1Shrub", "hab1Swamp",
    "hab1UrbInd", "hab1WetGrass", "hab1WetShrub", "wtAge", "wtAge2",
    "wtAge05", "fCC2",
    "isCon:wtAge", "isCon:wtAge2",
    "isUpCon:wtAge", "isBSLarch:wtAge", "isUpCon:wtAge2", "isBSLarch:wtAge2",
    "isMix:wtAge", "isPine:wtAge", "isWSpruce:wtAge", "isMix:wtAge2",
    "isPine:wtAge2", "isWSpruce:wtAge2", "isCon:wtAge05", "isUpCon:wtAge05",
    "isBSLarch:wtAge05", "isMix:wtAge05", "isPine:wtAge05", "isWSpruce:wtAge05")
## 0 out in South
cns0 <- c("pAspen", "ROAD01", "SoftLin_PC", "ARU2ARU", "YR", "habCl:ROAD01")
## habitat in South
cnsHab <- c("(Intercept)", "soil1RapidDrain", "soil1Saline", "soil1Clay",
    "soil1Cult", "soil1UrbInd", "soil1vRapidDrain", "soil1vSalineAndClay",
    "soil1vCult", "soil1vUrbInd")
## climate (North & South)
## ASP and CTI included here
#cnClim <- c("xASP", "xCTI", "xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
#    "xMCMT", "xlat", "xlong", "xlat2", "xlong2", "xASP:xCTI", "xFFP:xMAP",
#    "xMAP:xPET", "xAHM:xMAT", "xlat:xlong", "WetKM", "WetWaterKM")
## ASP and CTI *not* included here
cnClim <- c("xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
    "xMCMT", "xlat", "xlong", "xlat2", "xlong2", "xFFP:xMAP",
    "xMAP:xPET", "xAHM:xMAT", "xlat:xlong", "WetKM", "WetWaterKM")
cnHF <- c("THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM",
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM",
    "Noncult2_KM")
cnClimHF <- c(cnClim, cnHF)

## model matrix for Clim & SurroundingHF
fclim <- as.formula(paste("~ - 1 +", paste(cnClimHF, collapse=" + ")))


en <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-north.Rdata"), envir=en)
xnn <- en$DAT[1:500,]
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-south.Rdata"), envir=es)
xns <- es$DAT[1:500,]
modss <- es$mods
yys <- es$YY
rm(en, es)

## model for species
fln <- list.files(file.path(ROOT, "out", "birds", "results", "north"))
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- list.files(file.path(ROOT, "out", "birds", "results", "south"))
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

## example for structure (trSoil, trVeg)
load(file.path(ROOT, "out", "transitions", "LowerPeace_LowerBorealHighlands.Rdata"))
tv$Sector2 <- factor(ifelse(is.na(tv$Sector), "NATIVE", as.character(tv$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))

ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")
ts$All <- as.character(ts$HF)
ts$All[is.na(ts$All)] <- as.character(ts$Levels1)[is.na(ts$All)]
ts$Sector2 <- factor(ifelse(is.na(ts$Sector), "NATIVE", as.character(ts$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))

ch2veg <- t(sapply(strsplit(colnames(trVeg), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2veg <- data.frame(ch2veg)
colnames(ch2veg) <- c("rf","cr")
rownames(ch2veg) <- colnames(trVeg)
## strata where abundance is assumed to be 0
ch2veg$rf_zero <- ch2veg$rf %in% c("NonVeg","Water")
ch2veg$cr_zero <- ch2veg$cr %in% c("NonVeg","Water",
    "BorrowpitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals",
    "RailHardSurface","RoadHardSurface",
    "MineSite", "PeatMine")
## strata that do not count into mean density calculation
ch2veg$HardLin <- ch2veg$cr %in% c("RailHardSurface","RailVegetatedVerge",
    "RoadHardSurface","RoadTrailVegetated","RoadVegetatedVerge")
#ch2veg$SoftLin <- ch2veg$cr %in% c("SeismicLine","TransmissionLine","Pipeline")
ch2veg$VegetatedLinear <- ch2veg$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "SeismicLine","TransmissionLine","Pipeline")
EARLY_SERAL <- c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "TransmissionLine","Pipeline")
if (SEISMIC_AS_EARLY_SERAL)
    EARLY_SERAL <- c(EARLY_SERAL, "SeismicLine")
ch2veg$EarlySeralLinear <- ch2veg$cr %in% EARLY_SERAL
ch2veg$CutLine <- ch2veg$cr == "SeismicLine"
## nonveg is terrestrial stratum, do not exclude here
## water is not terrestrial, age0 is all redistributed (so =0)
ch2veg$exclude <- ch2veg$cr %in% c("Water",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0",
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")
ch2veg$exclude[ch2veg$rf %in% c("Water",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0",
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")] <- TRUE
ch2veg$isHF <- ch2veg$cr %in% c("BorrowpitsDugoutsSumps",
    "Canals", "CCConif0", "CCConif1", "CCConif2", "CCConif3", "CCConif4",
    "CCConifR", "CCDecid0", "CCDecid1", "CCDecid2", "CCDecid3", "CCDecid4",
    "CCDecidR", "CCMixwood0", "CCMixwood1", "CCMixwood2", "CCMixwood3",
    "CCMixwood4", "CCMixwoodR", "CCPine0", "CCPine1", "CCPine2",
    "CCPine3", "CCPine4", "CCPineR",
    "CultivationCropPastureBareground", "HighDensityLivestockOperation",
    "IndustrialSiteRural",
    "MineSite",
    "MunicipalWaterSewage", "OtherDisturbedVegetation",
    "PeatMine", "Pipeline", "RailHardSurface",
    "RailVegetatedVerge", "Reservoirs", "RoadHardSurface", "RoadTrailVegetated",
    "RoadVegetatedVerge", "RuralResidentialIndustrial", "SeismicLine",
    "TransmissionLine", "Urban", "WellSite",
    "WindGenerationFacility")
ch2veg$Sector <- tv$Sector2[match(ch2veg$cr, tv$Combined)]
ch2veg$Sector[is.na(ch2veg$Sector)] <- "NATIVE" # strage Swamp thing
table(ch2veg$Sector,useNA="a")

ch2soil <- t(sapply(strsplit(colnames(trSoil), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2soil <- data.frame(ch2soil)
colnames(ch2soil) <- c("rf","cr")
rownames(ch2soil) <- colnames(trSoil)
## strata where abundance is assumed to be 0
ch2soil$rf_zero <- ch2soil$rf %in% c("NonVeg","Water",
    "SoilWater","SoilWetland","SoilUnknown")
ch2soil$cr_zero <- ch2soil$cr %in% c("NonVeg","Water",
    "SoilWater","SoilWetland","SoilUnknown","CutBlocks",
    "BorrowpitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals",
    "RailHardSurface","RoadHardSurface",
    "MineSite", "PeatMine")
## strata that do not count into mean density calculation
ch2soil$HardLin <- ch2soil$cr %in% c("RailHardSurface","RailVegetatedVerge",
    "RoadHardSurface","RoadTrailVegetated","RoadVegetatedVerge")
#ch2soil$SoftLin <- ch2soil$cr %in% c("SeismicLine","TransmissionLine","Pipeline")
ch2soil$VegetatedLinear <- ch2soil$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "SeismicLine","TransmissionLine","Pipeline")
ch2soil$EarlySeralLinear <- ch2soil$cr %in% EARLY_SERAL
ch2soil$CutLine <- ch2soil$cr == "SeismicLine"
ch2soil$exclude <- ch2soil$cr %in% c("SoilWater","SoilUnknown",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0",
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")
ch2soil$exclude[ch2soil$rf %in% c("SoilWater","SoilUnknown",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0",
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")] <- TRUE
ch2soil$isHF <- ch2soil$cr %in% c("BorrowpitsDugoutsSumps", "Canals",
    "CultivationCropPastureBareground",
    "CutBlocks", "HighDensityLivestockOperation", "IndustrialSiteRural",
    "MineSite", "MunicipalWaterSewage", "OtherDisturbedVegetation",
    "PeatMine", "Pipeline", "RailHardSurface", "RailVegetatedVerge",
    "Reservoirs", "RoadHardSurface", "RoadTrailVegetated",
    "RoadVegetatedVerge", "RuralResidentialIndustrial",
    "SeismicLine", "TransmissionLine",
    "Urban", "WellSite", "WindGenerationFacility")
ch2soil$Sector <- ts$Sector2[match(ch2soil$cr, ts$All)]
table(ch2soil$Sector,useNA="a")

XNhab <- as.matrix(read.csv("~/repos/abmianalytics/lookup/xn-veg-noburn.csv"))
colnames(XNhab) <- gsub("\\.", ":", colnames(XNhab))
colnames(XNhab)[1] <- "(Intercept)"
setdiff(rownames(XNhab), ch2veg$cr)
setdiff(ch2veg$cr, rownames(XNhab))
setdiff(cnnHab, colnames(XNhab))

## vegetated stuff in linear features is early seral based on backfilled
## and early seral means that age is 0
XNhab_es <- XNhab
XNhab_es[,"fCC2"] <- 0
XNhab_es[,grepl("wtAge", colnames(XNhab_es))] <- 0

XShab <- as.matrix(read.csv("~/repos/abmianalytics/lookup/xn-soil.csv"))
colnames(XShab)[1] <- "(Intercept)"
setdiff(rownames(XShab), ch2soil$cr)
setdiff(ch2soil$cr, rownames(XShab))
setdiff(cnsHab, colnames(XShab))


SPP0 <- union(fln, fls)
TAX <- read.csv("~/repos/abmispecies/_data/birds.csv")
SPP <- sort(as.character(TAX$AOU)[TAX$map.pred])
compare_sets(SPP,SPP0)
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL","RWBL")
#SPP <- SPP[!(SPP %in% c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL","RWBL"))]

STAGE <- list(
    veg =which(names(modsn) == "Space"),
    soil=which(names(modss) == "Space"))

for (spp in SPP) { # species START

cat("\n\n---", spp, which(spp==SPP), "/", length(SPP), "---\n")

fn <- file.path(ROOT, "out", "birds", "results", "north",
    paste0("birds_abmi-north_", spp, ".Rdata"))
resn <- loadSPP(fn)
fs <- file.path(ROOT, "out", "birds", "results", "south",
    paste0("birds_abmi-south_", spp, ".Rdata"))
ress <- loadSPP(fs)
estn <- suppressWarnings(getEst(resn, stage=STAGE$veg, na.out=FALSE, Xnn))
ests <- suppressWarnings(getEst(ress, stage=STAGE$soil, na.out=FALSE, Xns))

NSest <- c(north=!is.null(resn), south=!is.null(ress))

#regi <- "LowerAthabasca_CentralMixedwood"
#date()
for (regi in regs) { # regions START

t0 <- proc.time()
cat(spp, regi, which(regs==regi), "/", length(regs), "\n")
flush.console()
gc()

load(file.path(ROOT, "out", "transitions", paste0(regi,".Rdata")))

if (do10) {
    #trVeg <- groupMeans(trVeg, 1, kgrid0[rownames(trVeg), "id10"])
    #trSoil <- groupMeans(trSoil, 1, kgrid0[rownames(trSoil), "id10"])
    trVeg <- groupMeans(trVeg, 1, kgrid0[rownames(trVeg), "Row10_Col10"])
    trSoil <- groupMeans(trSoil, 1, kgrid0[rownames(trSoil), "Row10_Col10"])
}

ii <- rownames(kgrid) %in% rownames(trVeg)

iib <- if (do10)
    ii else ii & kgrid$Rnd10 <= PROP

Xclim <- model.matrix(fclim, kgrid[ii,,drop=FALSE])
colnames(Xclim) <- fixNames(colnames(Xclim))
XclimB <- model.matrix(fclim, kgrid[iib,,drop=FALSE])
colnames(XclimB) <- fixNames(colnames(XclimB))
## reference has 0 surrounding HF
## but not surrounding Wet etc, which needs to reflect backfilled
Xclim0 <- Xclim
Xclim0[,cnHF] <- 0
Xclim0[,"WetKM"] <- kgrid$WetKM0[ii]
Xclim0[,"WetWaterKM"] <- kgrid$WetWaterKM0[ii]
XclimB0 <- XclimB
XclimB0[,cnHF] <- 0
XclimB0[,"WetKM"] <- kgrid$WetKM0[iib]
XclimB0[,"WetWaterKM"] <- kgrid$WetWaterKM0[iib]
estnClim <- estn[,colnames(Xclim),drop=FALSE]
estsClim <- ests[,colnames(Xclim),drop=FALSE]

## north - current
logPNclim1 <- Xclim %*% estnClim[1,]
#logPNclimB <- apply(estnClim, 1, function(z) XclimB %*% z)
logPNclimB <- do.call(cbind, lapply(1:min(BMAX, nrow(estnClim)),
    function(i) XclimB %*% estnClim[i,]))

## south - current
logPSclim1 <- Xclim %*% estsClim[1,]
#logPSclimB <- apply(estsClim, 1, function(z) XclimB %*% z)
logPSclimB <- do.call(cbind, lapply(1:min(BMAX, nrow(estsClim)),
    function(i) XclimB %*% estsClim[i,]))

## north - reference
logPNclim01 <- Xclim0 %*% estnClim[1,]
#logPNclim0B <- apply(estnClim, 1, function(z) XclimB0 %*% z)
logPNclim0B <- do.call(cbind, lapply(1:min(BMAX, nrow(estnClim)),
    function(i) XclimB0 %*% estnClim[i,]))

## south - reference
logPSclim01 <- Xclim0 %*% estsClim[1,]
#logPSclim0B <- apply(estsClim, 1, function(z) XclimB0 %*% z)
logPSclim0B <- do.call(cbind, lapply(1:min(BMAX, nrow(estsClim)),
    function(i) XclimB0 %*% estsClim[i,]))

estnHab <- estn[,colnames(XNhab),drop=FALSE]
estsHab <- ests[,colnames(XShab),drop=FALSE]
## north
logPNhab1 <- XNhab %*% estnHab[1,]
logPNhabB <- apply(estnHab, 1, function(z) XNhab %*% z)
rownames(logPNhabB) <- rownames(logPNhab1)
logPNhab_es1 <- XNhab_es %*% estnHab[1,]
logPNhab_esB <- apply(estnHab, 1, function(z) XNhab_es %*% z)
rownames(logPNhab_esB) <- rownames(logPNhab_es1)
## south
logPShab1 <- XShab %*% estsHab[1,]
logPShabB <- apply(estsHab, 1, function(z) XShab %*% z)
rownames(logPShabB) <- rownames(logPShab1)

Aveg1 <- trVeg[rownames(kgrid)[ii],,drop=FALSE]
Aveg1[,ch2veg$exclude] <- 0
rs <- rowSums(Aveg1)
rs[rs <= 0] <- 1
Aveg1 <- Aveg1 / rs
AvegB <- trVeg[rownames(kgrid)[iib],,drop=FALSE]
AvegB[,ch2veg$exclude] <- 0
rs <- rowSums(AvegB)
rs[rs <= 0] <- 1
AvegB <- AvegB / rs

Asoil1 <- trSoil[rownames(kgrid)[ii],,drop=FALSE]
pSoil1 <- 1-Asoil1[,"SoilUnknown"]
Asoil1[,ch2soil$exclude] <- 0
rs <- rowSums(Asoil1)
rs[rs <= 0] <- 1
Asoil1 <- Asoil1 / rs
AsoilB <- trSoil[rownames(kgrid)[iib],,drop=FALSE]
pSoilB <- 1-AsoilB[,"SoilUnknown"]
AsoilB[,ch2soil$exclude] <- 0
rs <- rowSums(AsoilB)
rs[rs <= 0] <- 1
AsoilB <- AsoilB / rs

## pixel level results, North
pxNcr1 <- matrix(0, nrow(logPNclim1), 1)
pxNrf1 <- matrix(0, nrow(logPNclim1), 1)
pxNcrB <- matrix(0, nrow(logPNclimB), BMAX)
pxNrfB <- matrix(0, nrow(logPNclimB), BMAX)
## pixel level results, South
pxScr1 <- matrix(0, nrow(logPSclim1), 1)
pxSrf1 <- matrix(0, nrow(logPSclim1), 1)
pxScrB <- matrix(0, nrow(logPSclimB), BMAX)
pxSrfB <- matrix(0, nrow(logPSclimB), BMAX)
## habitat level results, North
hbNcr1 <- matrix(0, nrow(ch2veg), 1)
hbNrf1 <- matrix(0, nrow(ch2veg), 1)
hbNcrB <- matrix(0, nrow(ch2veg), BMAX)
hbNrfB <- matrix(0, nrow(ch2veg), BMAX)
## habitat level results, South
hbScr1 <- matrix(0, nrow(ch2soil), 1)
hbSrf1 <- matrix(0, nrow(ch2soil), 1)
hbScrB <- matrix(0, nrow(ch2soil), BMAX)
hbSrfB <- matrix(0, nrow(ch2soil), BMAX)

estAsp <- ests[,"pAspen"]
## pAspen and pSoil is needed for combo
pAspen1 <- kgrid$pAspen[ii]
pAspenB <- kgrid$pAspen[iib]
## keeping track of cells
Cells <- ifelse(iib, 1L, 0L)[ii]
names(Cells) <- rownames(kgrid)[ii]

## 1st run
if (do1) {
    j <- 1
    ## North
    D_hab_cr <- exp(logPNhab1[match(ch2veg$cr, rownames(logPNhab1)),j])
    D_hab_rf <- exp(logPNhab1[match(ch2veg$rf, rownames(logPNhab1)),j])
    ## vegetated linear treated as early seral
    if (any(ch2veg$EarlySeralLinear))
        D_hab_cr[ch2veg$EarlySeralLinear] <- exp(logPNhab_es1[match(ch2veg$rf,
            rownames(logPNhab_es1)),j][ch2veg$EarlySeralLinear])
    ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption) --- ?????
    if (!SEISMIC_AS_EARLY_SERAL && any(ch2veg$CutLine))
        D_hab_cr[ch2veg$CutLine] <- D_hab_rf[ch2veg$CutLine]
    ## 0 density where either cr or rf hab is water or hard linear surface
    if (any(ch2veg$cr_zero))
        D_hab_cr[ch2veg$cr_zero] <- 0
    if (any(ch2veg$rf_zero))
        D_hab_cr[ch2veg$rf_zero] <- 0 # things like Water->Road
    if (any(ch2veg$rf_zero))
        D_hab_rf[ch2veg$rf_zero] <- 0

    AD_cr <- t(D_hab_cr * t(Aveg1)) * exp(logPNclim1[,j])
    AD_rf <- t(D_hab_rf * t(Aveg1)) * exp(logPNclim01[,j])
    pxNcr1[,j] <- rowSums(AD_cr)
    pxNrf1[,j] <- rowSums(AD_rf)
    hbNcr1[,j] <- colSums(AD_cr) / colSums(Aveg1)
    hbNrf1[,j] <- colSums(AD_rf) / colSums(Aveg1)
    pxNcrS <- groupSums(AD_cr, 2, ch2veg$Sector)
    pxNrfS <- groupSums(AD_rf, 2, ch2veg$Sector)

    ## South
    D_hab_cr <- exp(logPShab1[match(ch2soil$cr, rownames(logPShab1)),j])
    D_hab_rf <- exp(logPShab1[match(ch2soil$rf, rownames(logPShab1)),j])
    ## vegetated linear treated as early seral
    if (any(ch2soil$EarlySeralLinear))
        D_hab_cr[ch2soil$EarlySeralLinear] <- exp(logPShab1[match(ch2soil$rf,
            rownames(logPShab1)),j][ch2soil$EarlySeralLinear])
    ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption) --- ?????
    if (!SEISMIC_AS_EARLY_SERAL && any(ch2soil$CutLine))
        D_hab_cr[ch2soil$CutLine] <- D_hab_rf[ch2soil$CutLine]
    ## 0 density where either cr or rf hab is water or hard linear surface
    if (any(ch2soil$cr_zero))
        D_hab_cr[ch2soil$cr_zero] <- 0
    if (any(ch2soil$rf_zero))
        D_hab_cr[ch2soil$rf_zero] <- 0 # things like Water->Road
    if (any(ch2soil$rf_zero))
        D_hab_rf[ch2soil$rf_zero] <- 0

    AD_cr <- t(D_hab_cr * t(Asoil1)) * exp(logPSclim1[,j])
    AD_rf <- t(D_hab_rf * t(Asoil1)) * exp(logPSclim01[,j])
    ## pAspen based weighted average
    Asp <- exp(estAsp[j] * pAspen1)
    AD_cr <- pAspen1 * (Asp * AD_cr) + (1-pAspen1) * AD_cr
    AD_rf <- pAspen1 * (Asp * AD_rf) + (1-pAspen1) * AD_rf
    pxScr1[,j] <- rowSums(AD_cr)
    pxSrf1[,j] <- rowSums(AD_rf)
    hbScr1[,j] <- colSums(AD_cr) / colSums(Asoil1)
    hbSrf1[,j] <- colSums(AD_rf) / colSums(Asoil1)
    pxScrS <- groupSums(AD_cr, 2, ch2soil$Sector)
    pxSrfS <- groupSums(AD_rf, 2, ch2soil$Sector)
}

## BMAX runs
if (doB || do10) {
    for (j in 1:BMAX) {
        if (NSest["north"]) {
            ## North
            D_hab_cr <- exp(logPNhabB[match(ch2veg$cr, rownames(logPNhabB)),j])
            D_hab_rf <- exp(logPNhabB[match(ch2veg$rf, rownames(logPNhabB)),j])
            ## vegetated linear treated as early seral
            if (any(ch2veg$EarlySeralLinear))
                D_hab_cr[ch2veg$EarlySeralLinear] <- exp(logPNhab_esB[match(ch2veg$rf,
                    rownames(logPNhab_esB)),j][ch2veg$EarlySeralLinear])
            ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption)
            if (!SEISMIC_AS_EARLY_SERAL && any(ch2veg$CutLine))
                D_hab_cr[ch2veg$CutLine] <- D_hab_rf[ch2veg$CutLine]
            ## 0 density where either cr or rf hab is water or hard linear surface
            if (any(ch2veg$cr_zero))
                D_hab_cr[ch2veg$cr_zero] <- 0
            if (any(ch2veg$rf_zero))
                D_hab_cr[ch2veg$rf_zero] <- 0 # things like Water->Road
            if (any(ch2veg$rf_zero))
                D_hab_rf[ch2veg$rf_zero] <- 0

            AD_cr <- t(D_hab_cr * t(AvegB)) * exp(logPNclimB[,j])
            AD_rf <- t(D_hab_rf * t(AvegB)) * exp(logPNclim0B[,j])
            pxNcrB[,j] <- rowSums(AD_cr)
            pxNrfB[,j] <- rowSums(AD_rf)
            hbNcrB[,j] <- colSums(AD_cr) / colSums(AvegB)
            hbNrfB[,j] <- colSums(AD_rf) / colSums(AvegB)
        } else {
            pxNcrB[,j] <- 0
            pxNrfB[,j] <- 0
            hbNcrB[,j] <- 0
            hbNrfB[,j] <- 0
        }

        ## South
        if (NSest["south"]) {
            D_hab_cr <- exp(logPShabB[match(ch2soil$cr, rownames(logPShabB)),j])
            D_hab_rf <- exp(logPShabB[match(ch2soil$rf, rownames(logPShabB)),j])
            ## vegetated linear treated as early seral
            if (any(ch2soil$EarlySeralLinear))
                D_hab_cr[ch2soil$EarlySeralLinear] <- exp(logPShabB[match(ch2soil$rf,
                    rownames(logPShabB)),j][ch2soil$EarlySeralLinear])
            ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption)
            if (!SEISMIC_AS_EARLY_SERAL && any(ch2soil$CutLine))
                D_hab_cr[ch2soil$CutLine] <- D_hab_rf[ch2soil$CutLine]
            ## 0 density where either cr or rf hab is water or hard linear surface
            if (any(ch2soil$cr_zero))
                D_hab_cr[ch2soil$cr_zero] <- 0
            if (any(ch2soil$rf_zero))
                D_hab_cr[ch2soil$rf_zero] <- 0 # things like Water->Road
            if (any(ch2soil$rf_zero))
                D_hab_rf[ch2soil$rf_zero] <- 0

            AD_cr <- t(D_hab_cr * t(AsoilB)) * exp(logPSclimB[,j])
            AD_rf <- t(D_hab_rf * t(AsoilB)) * exp(logPSclim0B[,j])
            ## pAspen based weighted average
            Asp <- exp(estAsp[j] * pAspenB)
            AD_cr <- pAspenB * (Asp * AD_cr) + (1-pAspenB) * AD_cr
            AD_rf <- pAspenB * (Asp * AD_rf) + (1-pAspenB) * AD_rf
            pxScrB[,j] <- rowSums(AD_cr)
            pxSrfB[,j] <- rowSums(AD_rf)
            hbScrB[,j] <- colSums(AD_cr) / colSums(AsoilB)
            hbSrfB[,j] <- colSums(AD_rf) / colSums(AsoilB)
        } else {
            pxScrB[,j] <- 0
            pxSrfB[,j] <- 0
            hbScrB[,j] <- 0
            hbSrfB[,j] <- 0
        }
    }
}

TIME <- proc.time() - t0
if (do1) {
    if (!dir.exists(file.path(OUTDIR, "do1", spp)))
        dir.create(file.path(OUTDIR, "do1", spp))
    save(TIME, NSest,
        pxNcr1,pxNrf1,
        pxScr1,pxSrf1,
        hbNcr1,hbNrf1,
        hbScr1,hbSrf1,
        pxNcrS,pxNrfS,
        pxScrS,pxSrfS,
        pAspen1,pSoil1,Cells,
        file=file.path(OUTDIR, "do1", spp, paste0(regi, ".Rdata")))
}
if (doB) {
    if (!dir.exists(file.path(OUTDIR, "doB", spp)))
        dir.create(file.path(OUTDIR, "doB", spp))
    save(TIME, NSest,
        pxNcrB,pxNrfB,
        pxScrB,pxSrfB,
        hbNcrB,hbNrfB,
        hbScrB,hbSrfB,
        pAspenB,pSoilB,Cells,
        file=file.path(OUTDIR, "doB", spp, paste0(regi, ".Rdata")))
}
if (do10) {
    if (!dir.exists(file.path(OUTDIR, "do10", spp)))
        dir.create(file.path(OUTDIR, "do10", spp))
    save(TIME, NSest,
        pxNcrB,pxNrfB,
        pxScrB,pxSrfB,
        hbNcrB,hbNrfB,
        hbScrB,hbSrfB,
        pAspenB,pSoilB,Cells,
        file=file.path(OUTDIR, "do10", spp, paste0(regi, ".Rdata")))
}

} # regions END
} # species END
