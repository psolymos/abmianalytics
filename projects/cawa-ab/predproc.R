library(mefa4)

shf <- TRUE
doB <- TRUE

PROP <- 100
BMAX <- 240
if (!doB)
    BMAX <- 1
BMAX

ROOT <- "e:/peter/AB_data_v2016"
#ROOT2 <- "~/Dropbox/josm/2016/wewp"

OUTDIR1 <- "e:/peter/AB_data_v2016/out/birds/results/cawa/pred1"
OUTDIRB <- "e:/peter/AB_data_v2016/out/birds/results/cawa/predB"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

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
kgrid$WetPT <- kgrid$WetKM
kgrid$WetPT0 <- kgrid$WetKM0

## Dec/Mix
## 1km Dec80 for CAWA
cn_d <- c("DecidR", paste0("Decid", 1:9), "CCDecidR", paste0("CCDecid", 1:4))
cn_dm <- c(cn_d, "MixwoodR", paste0("Mixwood", 1:9), "CCMixwoodR", paste0("CCMixwood", 1:4))
cn_d80 <- paste0("Decid", 5:9)
cn_dm80 <- c(cn_d80, paste0("Mixwood", 5:9))
kgrid$DecKM <- rowSums(dd1km_pred$veg_current[,cn_d]) / rowSums(dd1km_pred$veg_current)
kgrid$DecMixKM <- rowSums(dd1km_pred$veg_current[,cn_dm]) / rowSums(dd1km_pred$veg_current)
kgrid$Dec80KM <- rowSums(dd1km_pred$veg_current[,cn_d80]) / rowSums(dd1km_pred$veg_current)
kgrid$DecMix80KM <- rowSums(dd1km_pred$veg_current[,cn_dm80]) / rowSums(dd1km_pred$veg_current)
cn_d <- cn_d[cn_d %in% colnames(dd1km_pred$veg_reference)]
cn_dm <- cn_dm[cn_dm %in% colnames(dd1km_pred$veg_reference)]
kgrid$DecKM0 <- rowSums(dd1km_pred$veg_reference[,cn_d]) / rowSums(dd1km_pred$veg_reference)
kgrid$DecMixKM0 <- rowSums(dd1km_pred$veg_reference[,cn_dm]) / rowSums(dd1km_pred$veg_reference)
kgrid$Dec80KM0 <- rowSums(dd1km_pred$veg_reference[,cn_d80]) / rowSums(dd1km_pred$veg_reference)
kgrid$DecMix80KM0 <- rowSums(dd1km_pred$veg_reference[,cn_dm80]) / rowSums(dd1km_pred$veg_reference)

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
    "xMAP:xPET", "xAHM:xMAT", "xlat:xlong")
cnHF <- c("THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", 
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", 
    "Noncult2_KM")
cnClimHF <- c(cnClim, cnHF, "xCTI", "WetPT", "WetPT:xCTI",
    "DecKM", "DecMixKM", "Dec80KM", "DecMix80KM")

## model matrix for Clim & SurroundingHF
fclim <- as.formula(paste("~ - 1 +", paste(cnClimHF, collapse=" + ")))


en <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-cawa.Rdata"), envir=en)
xnn <- en$DAT[1:500,]
modsn <- en$mods
yyn <- en$YY

## terms and design matrices
nTerms <- getTerms(modsn, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))

#regs <- levels(kgrid$LUFxNSR)
regs <- levels(droplevels(kgrid$LUFxNSR[kgrid$NRNAME != "Grassland"]))

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
#ch2veg$SoftLin <- ch2veg$cr %in% c("SeismicLine","TransmissionLine","Pipeline")
ch2veg$VegetatedLinear <- ch2veg$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "SeismicLine","TransmissionLine","Pipeline")
ch2veg$EarlySeralLinear <- ch2veg$cr %in% c("RailVegetatedVerge",
    "RoadTrailVegetated","RoadVegetatedVerge",
    "TransmissionLine","Pipeline")
ch2veg$CutLine <- ch2veg$cr == "SeismicLine"
## nonveg is terrestrial stratum, do not exclude here
## water is not terrestrial, age0 is all redistributed (so =0)
ch2veg$exclude <- ch2veg$cr %in% c("Water",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0", 
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")
ch2veg$exclude[ch2veg$rf %in% c("Water",
    "Conif0", "Decid0", "Mixwood0", "Pine0", "BSpr0", "Larch0", 
    "CCConif0", "CCDecid0", "CCMixwood0", "CCPine0")] <- TRUE

XNhab <- as.matrix(read.csv("~/repos/abmianalytics/lookup/xn-veg-noburn.csv"))
colnames(XNhab) <- gsub("\\.", ":", colnames(XNhab))
colnames(XNhab)[1] <- "(Intercept)"
setdiff(rownames(XNhab), ch2veg$cr)
setdiff(ch2veg$cr, rownames(XNhab))
setdiff(cnnHab, colnames(XNhab))

## vegetated stuff in linear features is early seral based on backfilled
XNhab_es <- XNhab
XNhab_es[,"fCC2"] <- 0
XNhab_es[,grepl("wtAge", colnames(XNhab_es))] <- 0

spp <- "CAWA"

## extra things to consider:
## CTI and WetPT
modsn$Wet
## CTI is taken from point intersections
## IDEALLY WetPT is either 0 (upland) or 1 (lowland), as part of ch2veg and the estimate
## BUT that makes habitat and spatial terms tangled up, so do the following:
## use WetKM instead of WetPT (they are highly correlated, and effect is linear
## so it might be justified)

STAGE <- list(veg =length(modsn) - ifelse(shf, 1, 2))

fn <- file.path(ROOT, "out", "birds", "results", "cawa", 
    paste0("birds_abmi-cawa_", spp, ".Rdata"))
resn <- loadSPP(fn)
estn <- suppressWarnings(getEst(resn, stage=STAGE$veg, na.out=FALSE, Xnn))

#regi <- "LowerAthabasca_LowerBorealHighlands"
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
## but not surrounding Wet etc, which needs to reflect backfilled
Xclim0 <- Xclim
Xclim0[,cnHF] <- 0
Xclim0[,"WetPT"] <- kgrid$WetKM0[ii]
Xclim0[,"WetPT:xCTI"] <- kgrid$WetKM0[ii] * kgrid$xCTI[ii]
Xclim0[,"DecKM"] <- kgrid$DecKM0[ii]
Xclim0[,"DecMixKM"] <- kgrid$DecMixKM0[ii]
Xclim0[,"Dec80KM"] <- kgrid$Dec80KM0[ii]
Xclim0[,"DecMix80KM"] <- kgrid$DecMix80KM0[ii]
XclimB0 <- XclimB
XclimB0[,cnHF] <- 0
XclimB0[,"WetPT"] <- kgrid$WetKM0[iib]
XclimB0[,"WetPT:xCTI"] <- kgrid$WetKM0[iib] * kgrid$xCTI[iib]
XclimB0[,"DecKM"] <- kgrid$DecKM0[iib]
XclimB0[,"DecMixKM"] <- kgrid$DecMixKM0[iib]
XclimB0[,"Dec80KM"] <- kgrid$Dec80KM0[iib]
XclimB0[,"DecMix80KM"] <- kgrid$DecMix80KM0[iib]
estnClim <- estn[,colnames(Xclim),drop=FALSE]

## north - current
logPNclim1 <- Xclim %*% estnClim[1,]
logPNclimB <- apply(estnClim, 1, function(z) XclimB %*% z)
## north - reference
logPNclim01 <- Xclim0 %*% estnClim[1,]
logPNclim0B <- apply(estnClim, 1, function(z) XclimB0 %*% z)

estnHab <- estn[,colnames(XNhab),drop=FALSE]
## north
logPNhab1 <- XNhab %*% estnHab[1,]
logPNhabB <- apply(estnHab, 1, function(z) XNhab %*% z)
rownames(logPNhabB) <- rownames(logPNhab1)
logPNhab_es1 <- XNhab_es %*% estnHab[1,]
logPNhab_esB <- apply(estnHab, 1, function(z) XNhab_es %*% z)
rownames(logPNhab_esB) <- rownames(logPNhab_es1)

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

## pixel level results, North
pxNcr1 <- matrix(0, nrow(logPNclim1), 1)
pxNrf1 <- matrix(0, nrow(logPNclim1), 1)
pxNcrB <- matrix(0, nrow(logPNclimB), BMAX)
pxNrfB <- matrix(0, nrow(logPNclimB), BMAX)
## habitat level results, North
hbNcr1 <- matrix(0, nrow(ch2veg), 1)
hbNrf1 <- matrix(0, nrow(ch2veg), 1)
hbNcrB <- matrix(0, nrow(ch2veg), BMAX)
hbNrfB <- matrix(0, nrow(ch2veg), BMAX)

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
}

## BMAX runs
if (doB) {
    for (j in 1:BMAX) {
        ## North
        D_hab_cr <- exp(logPNhabB[match(ch2veg$cr, rownames(logPNhabB)),j])
        ## vegetated linear (not cutline) treated as early seral
        if (any(ch2veg$EarlySeralLinear))
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
    }
}

TIME <- proc.time() - t0
if (TRUE) {
    if (!dir.exists(file.path(OUTDIR1)))
        dir.create(file.path(OUTDIR1))
    save(TIME, #NSest,
        pxNcr1,pxNrf1,
        #pxScr1,pxSrf1,
        hbNcr1,hbNrf1,
        #hbScr1,hbSrf1,
        #pAspen1,pSoil1,
        Cells,
        file=file.path(OUTDIR1, paste0(regi, ".Rdata")))
}
if (doB) {
    if (!dir.exists(file.path(OUTDIRB)))
        dir.create(file.path(OUTDIRB))
    save(TIME, #NSest,
        pxNcrB,pxNrfB,
        #pxScrB,pxSrfB,
        hbNcrB,hbNrfB,
        #hbScrB,hbSrfB,
        #pAspenB,pSoilB,
        Cells,
        file=file.path(OUTDIRB, paste0(regi, ".Rdata")))
}

} # regions END
