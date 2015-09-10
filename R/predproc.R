library(mefa4)
#library(raster)
#library(sp)
#library(rgdal)
#library(maptools)

ROOT <- "c:/p/AB_data_v2015"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

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
fclim <- as.formula(paste("~ - 1 +", paste(cnClim, collapse=" + ")))


en <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-useok-north.Rdata"), envir=en)
xnn <- en$DAT[1:500,]
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-useok-south.Rdata"), envir=es)
xns <- es$DAT[1:500,]
modss <- es$mods
yys <- es$YY
rm(en, es)

## model for species
fl <- list.files(file.path(ROOT, "out", "birds", "results"))
fln <- fl[grep("-north_", fl)]
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- fl[grep("-south_", fl)]
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

regs <- levels(kgrid$LUFxNSR)

## example for structure (trSoil, trVeg)
load(file.path(ROOT, "out", "transitions", "LowerPeace_LowerBorealHighlands.Rdata"))

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
ch2veg$SoftLin <- ch2veg$cr %in% c("SeismicLine","TransmissionLine","Pipeline")
ch2veg$VegetatedLinear <- ch2veg$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "SeismicLine","TransmissionLine","Pipeline")

ch2soil <- t(sapply(strsplit(colnames(trSoil), "->"), 
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2soil <- data.frame(ch2soil)
colnames(ch2soil) <- c("rf","cr")
rownames(ch2soil) <- colnames(trSoil)
## strata where abundance is assumed to be 0
ch2soil$rf_zero <- ch2soil$rf %in% c("NonVeg","Water")
ch2soil$cr_zero <- ch2soil$cr %in% c("NonVeg","Water",
    "BorrowpitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals",
    "RailHardSurface","RoadHardSurface",
    "MineSite", "PeatMine")
## strata that do not count into mean density calculation
ch2soil$HardLin <- ch2soil$cr %in% c("RailHardSurface","RailVegetatedVerge",
    "RoadHardSurface","RoadTrailVegetated","RoadVegetatedVerge")
ch2soil$SoftLin <- ch2soil$cr %in% c("SeismicLine","TransmissionLine","Pipeline")
ch2soil$VegetatedLinear <- ch2soil$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "SeismicLine","TransmissionLine","Pipeline")

XNhab <- as.matrix(read.csv("~/repos/abmianalytics/lookup/xn-veg.csv"))
colnames(XNhab) <- gsub("\\.", ":", colnames(XNhab))
colnames(XNhab)[1] <- "(Intercept)"
setdiff(rownames(XNhab), ch2veg$cr)
setdiff(ch2veg$cr, rownames(XNhab))
setdiff(cnnHab, colnames(XNhab))

## vegetated stuff in linear features is early seral based on backfilled
XNhab_es <- XNhab
XNhab_es[,"fCC2"] <- 0
XNhab_es[,grepl("wtAge", colnames(XNhab_es))] <- 0

XShab <- as.matrix(read.csv("~/repos/abmianalytics/lookup/xn-soil.csv"))
colnames(XShab)[1] <- "(Intercept)"
setdiff(rownames(XShab), ch2soil$cr)
setdiff(ch2soil$cr, rownames(XShab))
setdiff(cnsHab, colnames(XShab))


PROP <- 10
BMAX <- 100

regi <- "LowerPeace_LowerBorealHighlands"
spp <- "CAWA"

resn <- loadSPP(file.path(ROOT, "out", "birds", "results", 
    paste0("birds_abmi-north_", spp, ".Rdata")))
estn <- suppressWarnings(getEst(resn, stage=length(modsn)-2, na.out=FALSE, Xnn))
ress <- loadSPP(file.path(ROOT, "out", "birds", "results", 
    paste0("birds_abmi-south_", spp, ".Rdata")))
ests <- suppressWarnings(getEst(ress, stage=length(modss)-2, na.out=FALSE, Xns))

ii <- kgrid$LUFxNSR == regi
iib <- ii & kgrid$Rnd10 <= PROP
Xclim <- model.matrix(fclim, kgrid[ii,,drop=FALSE])
colnames(Xclim) <- fixNames(colnames(Xclim))
XclimB <- model.matrix(fclim, kgrid[iib,,drop=FALSE])
colnames(XclimB) <- fixNames(colnames(XclimB))
estnClim <- estn[1:BMAX,colnames(Xclim)]
estsClim <- ests[1:BMAX,colnames(Xclim)]
## north
logPNclim1 <- Xclim %*% estnClim[1,]
logPNclimB <- apply(estnClim, 1, function(z) XclimB %*% z)
## south
logPSclim1 <- Xclim %*% estsClim[1,]
logPSclimB <- apply(estsClim, 1, function(z) XclimB %*% z)

estnHab <- estn[1:BMAX,colnames(XNhab)]
estsHab <- ests[1:BMAX,colnames(XShab)]
## north
logPNhab1 <- XNhab %*% estnHab[1,]
logPNhabB <- apply(estnHab, 1, function(z) XNhab %*% z)
logPNhab_es1 <- XNhab_es %*% estnHab[1,]
logPNhab_esB <- apply(estnHab, 1, function(z) XNhab_es %*% z)
## south
logPShab1 <- XShab %*% estsHab[1,]
logPShabB <- apply(estsHab, 1, function(z) XShab %*% z)

Aveg1 <- trVeg[rownames(kgrid)[ii],,drop=FALSE]
Aveg1 <- Aveg1 / rowSums(Aveg1)
AvegB <- trVeg[rownames(kgrid)[iib],,drop=FALSE]
AvegB <- AvegB / rowSums(AvegB)

Asoil1 <- trSoil[rownames(kgrid)[ii],,drop=FALSE]
Asoil1 <- Asoil1 / rowSums(Asoil1)
AsoilB <- trSoil[rownames(kgrid)[iib],,drop=FALSE]
AsoilB <- AsoilB / rowSums(AsoilB)

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

## pAspen and pSoil is needed for combo
pAspen1 <- kgrid$pAspen[ii]
pAspenB <- kgrid$pAspen[iib]
pSoil1 <- 1-Asoil1[,"SoilUnknown"]
pSoilB <- 1-AsoilB[,"SoilUnknown"]
## keeping track of cells
Cells <- ifelse(iib, 1L, 0L)[ii]
names(Cells) <- rownames(kgrid)[ii]

## 1st run
j <- 1
## North
D_hab_cr <- exp(logPNhab1[match(ch2veg$cr, rownames(logPNhab1)),j])
if (any(ch2veg$VegetatedLinear))
    D_hab_cr[ch2veg$VegetatedLinear] <- exp(logPNhab_es1[match(ch2veg$rf, 
        rownames(logPNhab_es1)),j][ch2veg$VegetatedLinear])
if (any(ch2veg$cr_zero))
    D_hab_cr[ch2veg$cr_zero] <- 0
if (any(ch2veg$rf_zero))
    D_hab_cr[ch2veg$rf_zero] <- 0 # things like Water->Road
D_hab_rf <- exp(logPNhab1[match(ch2veg$rf, rownames(logPNhab1)),j])
if (any(ch2veg$rf_zero))
    D_hab_rf[ch2veg$rf_zero] <- 0
AD_cr <- t(D_hab_cr * t(Aveg1)) * exp(logPNclim1[,j])
AD_rf <- t(D_hab_rf * t(Aveg1)) * exp(logPNclim1[,j])
pxNcr1[,j] <- rowSums(AD_cr)
pxNrf1[,j] <- rowSums(AD_rf)
hbNcr1[,j] <- colSums(AD_cr) / colSums(Aveg1)
hbNrf1[,j] <- colSums(AD_rf) / colSums(Aveg1)
## South
D_hab_cr <- exp(logPShab1[match(ch2soil$cr, rownames(logPShab1)),j])
if (any(ch2soil$VegetatedLinear))
    D_hab_cr[ch2soil$VegetatedLinear] <- exp(logPShab1[match(ch2soil$rf, 
        rownames(logPShab1)),j][ch2soil$VegetatedLinear])
if (any(ch2soil$cr_zero))
    D_hab_cr[ch2soil$cr_zero] <- 0
if (any(ch2soil$rf_zero))
    D_hab_cr[ch2soil$rf_zero] <- 0 # things like Water->Road
D_hab_rf <- exp(logPShab1[match(ch2soil$rf, rownames(logPShab1)),j])
if (any(ch2soil$rf_zero))
    D_hab_rf[ch2soil$rf_zero] <- 0
AD_cr <- t(D_hab_cr * t(Asoil1)) * exp(logPSclim1[,j])
AD_rf <- t(D_hab_rf * t(Asoil1)) * exp(logPSclim1[,j])
pxScr1[,j] <- rowSums(AD_cr)
pxSrf1[,j] <- rowSums(AD_rf)
hbScr1[,j] <- colSums(AD_cr) / colSums(Asoil1)
hbSrf1[,j] <- colSums(AD_rf) / colSums(Asoil1)

## BMAX runs
for (j in 1:BMAX) {
    ## North
    D_hab_cr <- exp(logPNhabB[match(ch2veg$cr, rownames(logPNhabB)),j])
    if (any(ch2veg$VegetatedLinear))
        D_hab_cr[ch2veg$VegetatedLinear] <- exp(logPNhab_esB[match(ch2veg$rf, 
            rownames(logPNhab_esB)),j][ch2veg$VegetatedLinear])
    if (any(ch2veg$cr_zero))
        D_hab_cr[ch2veg$cr_zero] <- 0
    if (any(ch2veg$rf_zero))
        D_hab_cr[ch2veg$rf_zero] <- 0 # things like Water->Road
    D_hab_rf <- exp(logPNhabB[match(ch2veg$rf, rownames(logPNhabB)),j])
    if (any(ch2veg$rf_zero))
        D_hab_rf[ch2veg$rf_zero] <- 0
    AD_cr <- t(D_hab_cr * t(AvegB)) * exp(logPNclimB[,j])
    AD_rf <- t(D_hab_rf * t(AvegB)) * exp(logPNclimB[,j])
    pxNcrB[,j] <- rowSums(AD_cr)
    pxNrfB[,j] <- rowSums(AD_rf)
    hbNcrB[,j] <- colSums(AD_cr) / colSums(AvegB)
    hbNrfB[,j] <- colSums(AD_rf) / colSums(AvegB)
    ## South
    D_hab_cr <- exp(logPShabB[match(ch2soil$cr, rownames(logPShabB)),j])
    if (any(ch2soil$VegetatedLinear))
        D_hab_cr[ch2soil$VegetatedLinear] <- exp(logPShabB[match(ch2soil$rf, 
            rownames(logPShabB)),j][ch2soil$VegetatedLinear])
    if (any(ch2soil$cr_zero))
        D_hab_cr[ch2soil$cr_zero] <- 0
    if (any(ch2soil$rf_zero))
        D_hab_cr[ch2soil$rf_zero] <- 0 # things like Water->Road
    D_hab_rf <- exp(logPShabB[match(ch2soil$rf, rownames(logPShabB)),j])
    if (any(ch2soil$rf_zero))
        D_hab_rf[ch2soil$rf_zero] <- 0
    AD_cr <- t(D_hab_cr * t(AsoilB)) * exp(logPSclimB[,j])
    AD_rf <- t(D_hab_rf * t(AsoilB)) * exp(logPSclimB[,j])
    pxScrB[,j] <- rowSums(AD_cr)
    pxSrfB[,j] <- rowSums(AD_rf)
    hbScrB[,j] <- colSums(AD_cr) / colSums(AsoilB)
    hbSrfB[,j] <- colSums(AD_rf) / colSums(AsoilB)
}

