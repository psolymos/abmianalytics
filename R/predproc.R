library(mefa4)

shf <- FALSE
doB <- TRUE

PROP <- 10
BMAX <- 100
if (!doB)
    BMAX <- 1
BMAX

ROOT <- "c:/p/AB_data_v2015"

## surrounding hf
if (shf) {
    OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-shf-1/"
    OUTDIRB <- "e:/peter/sppweb2015/birds-pred-shf-B/"
} else {
    OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-1/"
    OUTDIRB <- "e:/peter/sppweb2015/birds-pred-B/"
}

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
    "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xlat:xlong")
cnHF <- c("THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", 
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", "Noncult2_KM")
cnClimHF <- c(cnClim, cnHF)

## model matrix for Clim & SurroundingHF
fclim <- as.formula(paste("~ - 1 +", paste(cnClimHF, collapse=" + ")))


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
ch2veg$EarlySeralLinear <- ch2veg$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "TransmissionLine","Pipeline")
ch2veg$CutLine <- ch2veg$cr == "SeismicLine"
ch2veg$exclude <- ch2veg$cr %in% c("Water",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0", 
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")
ch2veg$exclude[ch2veg$rf %in% c("Water",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0", 
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")] <- TRUE

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
ch2soil$SoftLin <- ch2soil$cr %in% c("SeismicLine","TransmissionLine","Pipeline")
ch2soil$VegetatedLinear <- ch2soil$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "SeismicLine","TransmissionLine","Pipeline")
ch2soil$EarlySeralLinear <- ch2soil$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "TransmissionLine","Pipeline")
ch2soil$CutLine <- ch2soil$cr == "SeismicLine"
ch2soil$exclude <- ch2soil$cr %in% c("SoilWater","SoilUnknown",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0", 
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")
ch2soil$exclude[ch2soil$rf %in% c("SoilWater","SoilUnknown",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0", 
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")] <- TRUE

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


#spp <- "BTNW"
SPP <- union(fln, fls)
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL")
for (spp in SPP) { # species START

cat("\n\n---", spp, which(spp==SPP), "/", length(SPP), "---\n")

STAGE <- list(
    veg =length(modsn) - ifelse(shf, 1, 2), 
    soil=length(modss) - ifelse(shf, 1, 2))

resn <- loadSPP(file.path(ROOT, "out", "birds", "results", 
    paste0("birds_abmi-north_", spp, ".Rdata")))
ress <- loadSPP(file.path(ROOT, "out", "birds", "results", 
    paste0("birds_abmi-south_", spp, ".Rdata")))
estn <- suppressWarnings(getEst(resn, stage=STAGE$veg, na.out=FALSE, Xnn))
ests <- suppressWarnings(getEst(ress, stage=STAGE$soil, na.out=FALSE, Xns))

NSest <- c(north=!is.null(resn), south=!is.null(ress))


#if (spp == "PUMA") {
#    estn <- if (is.null(resn))
#        estn[rep(1, BMAX),,drop=FALSE] else estn[1:BMAX+1,,drop=FALSE]
#    ests <- if (is.null(ress))
#        ests[rep(1, BMAX),,drop=FALSE] else ests[1:BMAX+1,,drop=FALSE]
#} else {
    estn <- if (is.null(resn))
        estn[rep(1, BMAX),,drop=FALSE] else estn[1:BMAX,,drop=FALSE]
    ests <- if (is.null(ress))
        ests[rep(1, BMAX),,drop=FALSE] else ests[1:BMAX,,drop=FALSE]
#}


#regi <- "LowerAthabasca_CentralMixedwood"
#date()
for (regi in regs) { # regions START

t0 <- proc.time()
cat(spp, regi, which(regs==regi), "/", length(regs), "\n")
flush.console()
gc()

load(file.path(ROOT, "out", "transitions", paste0(regi,".Rdata")))


#ks <- kgrid[rownames(trVeg),]
#with(ks, plot(X,Y))

ii <- kgrid$LUFxNSR == regi
iib <- ii & kgrid$Rnd10 <= PROP
Xclim <- model.matrix(fclim, kgrid[ii,,drop=FALSE])
colnames(Xclim) <- fixNames(colnames(Xclim))
XclimB <- model.matrix(fclim, kgrid[iib,,drop=FALSE])
colnames(XclimB) <- fixNames(colnames(XclimB))
## reference has 0 surrounding HF
Xclim0 <- Xclim
Xclim0[,cnHF] <- 0
XclimB0 <- XclimB
XclimB0[,cnHF] <- 0
estnClim <- estn[,colnames(Xclim),drop=FALSE]
estsClim <- ests[,colnames(Xclim),drop=FALSE]


## north - current
logPNclim1 <- Xclim %*% estnClim[1,]
logPNclimB <- apply(estnClim, 1, function(z) XclimB %*% z)
## south - current
logPSclim1 <- Xclim %*% estsClim[1,]
logPSclimB <- apply(estsClim, 1, function(z) XclimB %*% z)
## north - reference
logPNclim01 <- Xclim0 %*% estnClim[1,]
logPNclim0B <- apply(estnClim, 1, function(z) XclimB0 %*% z)
## south - reference
logPSclim01 <- Xclim0 %*% estsClim[1,]
logPSclim0B <- apply(estsClim, 1, function(z) XclimB0 %*% z)

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
Asoil1[,ch2soil$exclude] <- 0
rs <- rowSums(Asoil1)
rs[rs <= 0] <- 1
Asoil1 <- Asoil1 / rs
AsoilB <- trSoil[rownames(kgrid)[iib],,drop=FALSE]
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
pSoil1 <- 1-Asoil1[,"SoilUnknown"]
pSoilB <- 1-AsoilB[,"SoilUnknown"]
## keeping track of cells
Cells <- ifelse(iib, 1L, 0L)[ii]
names(Cells) <- rownames(kgrid)[ii]

## 1st run
if (TRUE) {
    j <- 1
    ## North
    D_hab_cr <- exp(logPNhab1[match(ch2veg$cr, rownames(logPNhab1)),j])
    ## vegetated linear (not cutline) treated as early seral
    if (any(ch2veg$EarlySeralLinear))
        D_hab_cr[ch2veg$EarlySeralLinear] <- exp(logPNhab_es1[match(ch2veg$rf, 
            rownames(logPNhab_es1)),j][ch2veg$EarlySeralLinear])
    ## 0 density where either cr or rf hab is water or hard linear surface
    if (any(ch2veg$cr_zero))
        D_hab_cr[ch2veg$cr_zero] <- 0
    if (any(ch2veg$rf_zero))
        D_hab_cr[ch2veg$rf_zero] <- 0 # things like Water->Road
    D_hab_rf <- exp(logPNhab1[match(ch2veg$rf, rownames(logPNhab1)),j])
    if (any(ch2veg$rf_zero))
        D_hab_rf[ch2veg$rf_zero] <- 0
    ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption) --- ?????
    if (any(ch2veg$CutLine))
        D_hab_cr[ch2veg$CutLine] <- D_hab_rf[ch2veg$CutLine]

 #   if (any(ch2veg$exclude)) {
 #       D_hab_cr[ch2veg$exclude] <- 0
 #       D_hab_rf[ch2veg$exclude] <- 0
 #   }

    AD_cr <- t(D_hab_cr * t(Aveg1)) * exp(logPNclim1[,j])
    AD_rf <- t(D_hab_rf * t(Aveg1)) * exp(logPNclim01[,j])
    pxNcr1[,j] <- rowSums(AD_cr)
    pxNrf1[,j] <- rowSums(AD_rf)
    hbNcr1[,j] <- colSums(AD_cr) / colSums(Aveg1)
    hbNrf1[,j] <- colSums(AD_rf) / colSums(Aveg1)
    ## South
    D_hab_cr <- exp(logPShab1[match(ch2soil$cr, rownames(logPShab1)),j])
    ## vegetated linear (not cutline) treated as early seral
    if (any(ch2soil$EarlySeralLinear))
        D_hab_cr[ch2soil$EarlySeralLinear] <- exp(logPShab1[match(ch2soil$rf, 
            rownames(logPShab1)),j][ch2soil$EarlySeralLinear])
    if (any(ch2soil$cr_zero))
        D_hab_cr[ch2soil$cr_zero] <- 0
    if (any(ch2soil$rf_zero))
        D_hab_cr[ch2soil$rf_zero] <- 0 # things like Water->Road
    D_hab_rf <- exp(logPShab1[match(ch2soil$rf, rownames(logPShab1)),j])
    if (any(ch2soil$rf_zero))
        D_hab_rf[ch2soil$rf_zero] <- 0
    ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption) --- ?????
    if (any(ch2soil$CutLine))
        D_hab_cr[ch2soil$CutLine] <- D_hab_rf[ch2soil$CutLine]
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
}

## BMAX runs
if (doB) {
    for (j in 1:BMAX) {
        ## North
        D_hab_cr <- exp(logPNhabB[match(ch2veg$cr, rownames(logPNhabB)),j])
        ## vegetated linear (not cutline) treated as early seral
        if (any(ch2soil$EarlySeralLinear))
            D_hab_cr[ch2veg$EarlySeralLinear] <- exp(logPNhab_esB[match(ch2veg$rf, 
                rownames(logPNhab_esB)),j][ch2veg$EarlySeralLinear])
        if (any(ch2veg$cr_zero))
            D_hab_cr[ch2veg$cr_zero] <- 0
        if (any(ch2veg$rf_zero))
            D_hab_cr[ch2veg$rf_zero] <- 0 # things like Water->Road
        D_hab_rf <- exp(logPNhabB[match(ch2veg$rf, rownames(logPNhabB)),j])
        if (any(ch2veg$rf_zero))
            D_hab_rf[ch2veg$rf_zero] <- 0
        ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption)
        if (any(ch2veg$CutLine))
            D_hab_cr[ch2veg$CutLine] <- D_hab_rf[ch2veg$CutLine]
        AD_cr <- t(D_hab_cr * t(AvegB)) * exp(logPNclimB[,j])
        AD_rf <- t(D_hab_rf * t(AvegB)) * exp(logPNclim0B[,j])
        pxNcrB[,j] <- rowSums(AD_cr)
        pxNrfB[,j] <- rowSums(AD_rf)
        hbNcrB[,j] <- colSums(AD_cr) / colSums(AvegB)
        hbNrfB[,j] <- colSums(AD_rf) / colSums(AvegB)
        ## South
        D_hab_cr <- exp(logPShabB[match(ch2soil$cr, rownames(logPShabB)),j])
        ## vegetated linear (not cutline) treated as early seral
        if (any(ch2soil$EarlySeralLinear))
            D_hab_cr[ch2soil$EarlySeralLinear] <- exp(logPShabB[match(ch2soil$rf, 
                rownames(logPShabB)),j][ch2soil$EarlySeralLinear])
        if (any(ch2soil$cr_zero))
            D_hab_cr[ch2soil$cr_zero] <- 0
        if (any(ch2soil$rf_zero))
            D_hab_cr[ch2soil$rf_zero] <- 0 # things like Water->Road
        D_hab_rf <- exp(logPShabB[match(ch2soil$rf, rownames(logPShabB)),j])
        if (any(ch2soil$rf_zero))
            D_hab_rf[ch2soil$rf_zero] <- 0
        ## cutlines are backfilled (surrounding HF effect applies, + behavioural assumption)
        if (any(ch2soil$CutLine))
            D_hab_cr[ch2soil$CutLine] <- D_hab_rf[ch2soil$CutLine]
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
    }
}

TIME <- proc.time() - t0
if (TRUE) {
    if (!dir.exists(file.path(OUTDIR1, spp)))
        dir.create(file.path(OUTDIR1, spp))
    save(TIME, NSest,
        pxNcr1,pxNrf1,
        pxScr1,pxSrf1,
        hbNcr1,hbNrf1,
        hbScr1,hbSrf1,
        pAspen1,pSoil1,Cells,
        file=file.path(OUTDIR1, spp, paste0(regi, ".Rdata")))
}
if (doB) {
    if (!dir.exists(file.path(OUTDIRB, spp)))
        dir.create(file.path(OUTDIRB, spp))
    save(TIME, NSest,
        pxNcrB,pxNrfB,
        pxScrB,pxSrfB,
        hbNcrB,hbNrfB,
        hbScrB,hbSrfB,
        pAspenB,pSoilB,Cells,
        file=file.path(OUTDIRB, spp, paste0(regi, ".Rdata")))
}

} # regions END
} # species END

