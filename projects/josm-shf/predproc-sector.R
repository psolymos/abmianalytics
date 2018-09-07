library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"
OUTDIR <- paste0("e:/peter/josm/2018/hshfix")

STAGE <- list(veg = 7) # hab=5, hab+clim=6, hab+clim+shf=7
shf <- TRUE # surrounding HF
do1 <- TRUE # do only 1st run
doB <- FALSE # do bootstrap
## seismic lines can be treated as early seral when no surrounding effects considered
## or apply the backfilled value (age set to 0) when surrounding hf is considered
SEISMIC_AS_EARLY_SERAL <- TRUE

PROP <- 100
BMAX <- 100
if (!doB)
    BMAX <- 1
BMAX

#SPP <- "CAWA"

## surrounding HF only includes `sectors` only
#sect <- "All"
#sect <- "Agriculture"
#sect <- "EnergyLin"
#sect <- "EnergyMW"
#sect <- "Forestry"
#sect <- "Misc"
#sect <- "RuralUrban"
#sect <- "Transportation"
sectors_all <- c("Agriculture", "EnergyLin", "EnergyMW", "Forestry", "Misc",
    "RuralUrban","Transportation")

for (sect in c("All", sectors_all)) { # sect start

stopifnot(length(sect)==1)
sectors <- if (sect != "All")
    sect else sectors_all

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata")) # kgrid
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

en <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-josmshf.Rdata"), envir=en)
xnn <- en$DAT[1:500,]
modsn <- en$mods
yyn <- en$YY

## model for species
fln <- list.files(file.path(ROOT, "out", "birds", "results", "josmshf"))
fln <- sub("birds_abmi-josmshf_", "", fln)
fln <- sub(".Rdata", "", fln)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))

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
stopifnot(all(rownames(kgrid) == rownames(dd1km_pred$veg_current)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv$SectorForSeMs <- factor(ifelse(is.na(tv$SectorForSeMs), "NATIVE", as.character(tv$SectorForSeMs)),
    c("NATIVE", "Agriculture", "EnergyLin", "EnergyMW",
    "Forestry", "Misc", "RuralUrban","Transportation"))
tv$hab1ec <- paste0(tv$Type, tv$EC_AGE)

VEGALL <- dd1km_pred$veg_current
VEGALL0 <- dd1km_pred$veg_reference
RS <- rowSums(VEGALL)
rm(dd1km_pred)

## partial backfill here: surrounding HF should only include focals
stopifnot(all(colnames(VEGALL) == rownames(tv)))
hfcn <- tv$SectorForSeMs %in% sectors_all & !(tv$SectorForSeMs %in% sectors)
addmargins(table(tv$SectorForSeMs,hfcn))
if (any(hfcn)) {
    V1 <- VEGALL[,!hfcn]
    V2 <- VEGALL[,hfcn]
    V2[] <- 0
    V3 <- cbind(V1, V2)
    VEGALL <- V3[,colnames(VEGALL)]
    rm(V1, V2, V3)
}

## surrounding HF at 1km scale
kgrid$THF_KM <- rowSums(VEGALL[,setdiff(colnames(VEGALL),
    colnames(VEGALL0))]) / RS
kgrid$Lin_KM <- rowSums(VEGALL[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailHardSurface", "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / RS
kgrid$Nonlin_KM <- kgrid$THF_KM - kgrid$Lin_KM
kgrid$Cult_KM <- rowSums(VEGALL[,c("CultivationCropPastureBareground",
    "HighDensityLivestockOperation")]) / RS
kgrid$Noncult_KM <- kgrid$THF_KM - kgrid$Cult_KM
CClabs <- colnames(VEGALL)[grep("CC", colnames(VEGALL))]
kgrid$Succ_KM <- rowSums(VEGALL[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge",
    CClabs)]) / RS
kgrid$Alien_KM <- kgrid$THF_KM - kgrid$Succ_KM

kgrid$THF2_KM <- kgrid$THF_KM^2
kgrid$Succ2_KM <- kgrid$Succ_KM^2
kgrid$Alien2_KM <- kgrid$Alien_KM^2
kgrid$Noncult2_KM <- kgrid$Noncult_KM^2
kgrid$Nonlin2_KM <- kgrid$Nonlin_KM^2

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
## ASP and CTI *not* included here
cnClim <- c("xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
    "xMCMT", "xlat", "xlong", "xlat2", "xlong2", "xFFP:xMAP",
    "xMAP:xPET", "xAHM:xMAT", "xlat:xlong")
cnHF <- c("THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM",
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM",
    "Noncult2_KM")
cnClimHF <- c(cnClim, cnHF)

## model matrix for Clim & SurroundingHF
fclim <- as.formula(paste("~ - 1 +", paste(cnClimHF, collapse=" + ")))

regtmp <- nonDuplicated(kgrid[,c("LUFxNSR", "NRNAME", "NSRNAME")], kgrid$LUFxNSR, TRUE)
#regs <- rownames(regtmp)[regtmp$NRNAME != "Grassland"]
regs <- rownames(regtmp)

## example for structure (trSoil, trVeg)
load(file.path(ROOT, "out", "transitions", "LowerPeace_LowerBorealHighlands.Rdata"))

ch2veg <- t(sapply(strsplit(colnames(trVeg), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2veg <- data.frame(ch2veg)
colnames(ch2veg) <- c("rf","cr")
rownames(ch2veg) <- colnames(trVeg)

ch2veg$Sector <- tv$SectorForSeMs[match(ch2veg$cr, tv$Combined)]
ch2veg$Sector[is.na(ch2veg$Sector)] <- "NATIVE" # strage Swamp thing
table(ch2veg$Sector,useNA="a")

## partial backfill here: need to backfill everything that is not part of focus
hfcn <- ch2veg$Sector %in% sectors_all & !(ch2veg$Sector %in% sectors)
addmargins(table(ch2veg$Sector,hfcn))
ch2veg$cr[hfcn] <- ch2veg$rf[hfcn]

## hsh definitions
ch2veg$hsh_cr <- tv$hab1ec[match(ch2veg$cr, tv$Combined)]
ch2veg$hsh_rf <- tv$hab1ec[match(ch2veg$rf, tv$Combined)]

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

#SPP <- "CAWA"
do_hsh <- FALSE
#do_veg <- TRUE
#spp <- SPP

SPP <- fln

#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL","RWBL")
#SPP <- SPP[!(SPP %in% c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL","RWBL"))]

for (spp in SPP) { # species START

cat("\n\n---", spp, which(spp==SPP), "/", length(SPP), "---\n")

fn <- file.path(ROOT, "out", "birds", "results", "josmshf",
    paste0("birds_abmi-josmshf_", spp, ".Rdata"))
resn <- loadSPP(fn)
estn <- suppressWarnings(getEst(resn, stage=STAGE$veg, na.out=FALSE, Xnn))
estn <- if (is.null(resn))
    estn[rep(1, BMAX),,drop=FALSE] else estn[1:BMAX,,drop=FALSE]

est1 <- c(estn[1,], "HSH_KM"=0, "HSH2_KM"=0, "HSH05_KM"=0)
de <- new.env()
load(paste0("e:/peter/josm/2018/hsh-estimates4x/", spp, ".Rdata"), envir=de)
hsh_lab <- NULL
if (de$out$mid > 0) {
    # 1: HSH, 2: HSH05, 3: HSH+HSH2, 4: HSH05+HSH
    est1[] <- 0
    est1[names(de$out$hsh_coef[[de$out$mid]])] <- de$out$hsh_coef[[de$out$mid]]
    hsh_lab <- de$out$hsh_labels
}

    #regi <- "LowerAthabasca_CentralMixedwood"
    #date()
    for (regi in regs) { # regions START

        t0 <- proc.time()
        cat("Predict -", sect, spp, regi, which(regs==regi), "/", length(regs), "\n")
        flush.console()
        gc()

        load(file.path(ROOT, "out", "transitions", paste0(regi,".Rdata")))
        ii <- kgrid$LUFxNSR == regi

        Aveg1all <- trVeg[rownames(kgrid)[ii],,drop=FALSE]
        Aveg1    <- trVeg[rownames(kgrid)[ii],,drop=FALSE]
        Aveg1[,ch2veg$exclude] <- 0
        rs <- rowSums(Aveg1)
        rs[rs <= 0] <- 1
        Aveg1 <- Aveg1 / rs

        ## HSH definition
        PhshCr <- rowSums(Aveg1[, ch2veg$hsh_cr %in% hsh_lab])
        PhshRf <- rowSums(Aveg1[, ch2veg$hsh_rf %in% hsh_lab])

        Xclim <- model.matrix(fclim, kgrid[ii,,drop=FALSE])
        Xclim <- cbind(Xclim, HSH_KM=PhshCr, HSH2_KM=PhshCr^2, HSH05_KM=sqrt(PhshCr))
        colnames(Xclim) <- fixNames(colnames(Xclim))

        ## reference has 0 surrounding HF
        ## but not surrounding Wet etc, which needs to reflect backfilled
        Xclim0 <- Xclim
        Xclim0[,cnHF] <- 0
        Xclim0[,"HSH_KM"] <- PhshRf
        Xclim0[,"HSH2_KM"] <- PhshRf^2
        Xclim0[,"HSH05_KM"] <- sqrt(PhshRf)
        estnClim <- est1[colnames(Xclim)]

        ## north - current
        logPNclim1 <- Xclim %*% estnClim
        ## north - reference
        logPNclim01 <- Xclim0 %*% estnClim

        estnHab <- est1[colnames(XNhab)]
        logPNhab1 <- XNhab %*% estnHab
        logPNhab_es1 <- XNhab_es %*% estnHab

        ## keeping track of cells
        Cells <- rep(1L, sum(ii))
        names(Cells) <- rownames(kgrid)[ii]

        ## 1st run
        if (do1) {
            j <- 1
            ## North
            D_hab_cr <- exp(logPNhab1[match(ch2veg$cr, rownames(logPNhab1)),j])
            D_hab_rf <- exp(logPNhab1[match(ch2veg$rf, rownames(logPNhab1)),j])
            ## vegetated linear (not cutline) treated as early seral
            if (any(ch2veg$EarlySeralLinear))
                D_hab_cr[ch2veg$EarlySeralLinear] <- exp(logPNhab_es1[match(ch2veg$rf,
                    rownames(logPNhab_es1)),j][ch2veg$EarlySeralLinear])
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

            AD_cr <- t(D_hab_cr * t(Aveg1)) * exp(logPNclim1[,j])
            AD_rf <- t(D_hab_rf * t(Aveg1)) * exp(logPNclim01[,j])
            pxNcrS <- groupSums(AD_cr, 2, ch2veg$Sector)
            pxNrfS <- groupSums(AD_rf, 2, ch2veg$Sector)
            if (sect == "All") {
                abund <- list(
                    Ntr=colSums(100 * AD_cr),
                    Ncr=colSums(groupSums(100 * AD_cr, 2, ch2veg$cr)),
                    Nrf=colSums(groupSums(100 * AD_rf, 2, ch2veg$rf)),
                    Atr=colSums(Aveg1all/10^6),
                    Acr=colSums(groupSums(Aveg1all/10^6, 2, ch2veg$cr)),
                    Arf=colSums(groupSums(Aveg1all/10^6, 2, ch2veg$rf)))
            }

        }

        TIME <- proc.time() - t0
        if (do1) {
            if (!dir.exists(file.path(OUTDIR, sect)))
                dir.create(file.path(OUTDIR, sect))
            if (!dir.exists(file.path(OUTDIR, sect, spp)))
                dir.create(file.path(OUTDIR, sect, spp))
            if (!dir.exists(file.path(OUTDIR, sect, spp)))
                dir.create(file.path(OUTDIR, sect, spp))
            toSave <- c("TIME", "pxNcrS", "pxNrfS", "Cells")
            if (sect == "All")
                toSave <- c(toSave, "abund")
            save(list=toSave,
                file=file.path(OUTDIR, sect, spp, paste0(regi, ".Rdata")))
        }


    } # regions END

    ## putting things together (pixel x sector)
    fl <- list.files(file.path(OUTDIR, sect, spp))
    cn <- c("NATIVE", sectors_all)
    OUTcr <- matrix(0, nrow(kgrid), length(cn))
    rownames(OUTcr) <- rownames(kgrid)
    colnames(OUTcr) <- cn
    OUTrf <- OUTcr
    for (i in 1:length(fl)) {
        cat("Combine -", sect, spp, i, "/", length(fl), "\n");flush.console()
        e <- new.env()
        load(file.path(OUTDIR, sect, spp, fl[i]), envir=e)
        Cells <- names(e$Cells)
        OUTcr[Cells,] <- as.matrix(e$pxNcrS[,cn])
        OUTrf[Cells,] <- as.matrix(e$pxNrfS[,cn])
    }
#    for (i in 1:length(cn)) {
#        q <- quantile(c(OUTcr[,i], OUTrf[,i]), 0.99)
#        OUTcr[OUTcr[,i] > q,i] <- q
#        OUTrf[OUTrf[,i] > q,i] <- q
#    }
    colnames(OUTcr)[colnames(OUTcr) == "NATIVE"] <- "Native"
    colnames(OUTrf)[colnames(OUTrf) == "NATIVE"] <- "Native"
    SA.Curr <- OUTcr
    SA.Ref <- OUTrf
    save(SA.Curr, SA.Ref,
        file=file.path(OUTDIR, sect, paste0(spp, ".RData")))

    if (sect == "All")
        unlink(file.path(OUTDIR, sect, spp), recursive=TRUE, force=TRUE)

} # species END

} # end sect

library(sendmailR)
sendmail(sprintf("<%s>", "bot@abmi.ca"),
            sprintf("<%s>", "psolymos@gmail.com"),
            "prediction done",
            list("done"),
            control=list(smtpServer="ASPMX.L.GOOGLE.COM"))

if (sect == "All") {
    NAD <- list()
    load(file.path(OUTDIR, "All", "OVEN", "LowerAthabasca_AthabascaPlain.Rdata"))
    for (spp in SPP) {
        fl <- list.files(file.path(OUTDIR, sect, spp))
        NOUTtr <- matrix(0, length(fl), length(abund$Ntr))
        NOUTcr <- matrix(0, length(fl), length(abund$Ncr))
        NOUTrf <- matrix(0, length(fl), length(abund$Nrf))
        rownames(NOUTtr) <- rownames(NOUTcr) <- rownames(NOUTrf) <- gsub(".Rdata", "", fl)
        colnames(NOUTtr) <- names(abund$Ntr)
        colnames(NOUTcr) <- names(abund$Ncr)
        colnames(NOUTrf) <- names(abund$Nrf)
        AOUTtr <- NOUTtr
        AOUTcr <- NOUTcr
        AOUTrf <- NOUTrf
        for (i in 1:length(fl)) {
            cat("Combine -", sect, spp, i, "/", length(fl), "\n");flush.console()
            e <- new.env()
            load(file.path(OUTDIR, sect, spp, fl[i]), envir=e)
            j <- gsub(".Rdata", "", fl[i])
            NOUTtr[j,] <- e$abund$Ntr
            NOUTcr[j,] <- e$abund$Ncr
            NOUTrf[j,] <- e$abund$Nrf
            AOUTtr[j,] <- e$abund$Atr
            AOUTcr[j,] <- e$abund$Acr
            AOUTrf[j,] <- e$abund$Arf
        }
        ## N in males, A in km^2, D in males/ha
        NOUTtr <- colSums(NOUTtr)
        NOUTcr <- colSums(NOUTcr)
        NOUTrf <- colSums(NOUTrf)
        AOUTtr <- colSums(AOUTtr)
        AOUTcr <- colSums(AOUTcr)
        AOUTrf <- colSums(AOUTrf)
        DOUTtr <- NOUTtr / (100 * AOUTtr)
        DOUTcr <- NOUTcr / (100 * AOUTcr)
        DOUTrf <- NOUTrf / (100 * AOUTrf)

        NAD[[spp]] <- list(
            tr=cbind(N=NOUTtr, A=AOUTtr, D=DOUTtr),
            cr=cbind(N=NOUTcr, A=AOUTcr, D=DOUTcr),
            rf=cbind(N=NOUTrf, A=AOUTrf, D=DOUTrf))
    }
}

save(NAD, file=file.path(OUTDIR, "NAD.RData"))

Dcr <- sapply(NAD, function(z) z$cr[,"D"])
Dcr <- Dcr[!grepl("0", rownames(Dcr)),]
Ncr <- sapply(NAD, function(z) z$cr[,"N"])
Ncr <- Ncr[!grepl("0", rownames(Ncr)),]
Acr <- NAD[[1]]$cr[,"A"] * 100 # ha
Acr <- Acr[!grepl("0", names(Acr))]
write.csv(data.frame(LC=rownames(Dcr), Dcr),
    row.names=FALSE, file=file.path(OUTDIR, "pop-density-by-landcover.csv"))
write.csv(data.frame(LC=rownames(Ncr), Area_ha=Acr, Ncr),
    row.names=FALSE, file=file.path(OUTDIR, "pop-abundance-by-landcover.csv"))


SPP <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATTW", "BANS", "BAOR",
    "BARS", "BAWW", "BBMA", "BBWA", "BBWO", "BCCH", "BHCO", "BHVI",
    "BLJA", "BLPW", "BOCH", "BRBL", "BRCR", "BTNW", "CAWA", "CCSP",
    "CEDW", "CHSP", "CLSW", "CMWA", "COGR", "CONW", "CORA", "COYE",
    "DEJU", "DOWO", "DUFL", "EAKI", "EAPH", "EUST", "EVGR", "FOSP",
    "GCKI", "GRAJ", "GRCA", "GRYE", "HAWO", "HETH", "HOSP", "HOWR",
    "KILL", "LCSP", "LEFL", "LEYE", "LISP", "MAWA", "MOBL", "MODO",
    "MOWA", "NESP", "NOFL", "NOWA", "OCWA", "OSFL", "OVEN", "PAWA",
    "PHVI", "PIGR", "PISI", "PIWO", "PUFI", "RBGR", "RBNU", "RCKI",
    "RECR", "REVI", "ROPI", "RUBL", "RUGR", "RWBL", "SAVS", "SOSA",
    "SOSP", "SPSA", "SWSP", "SWTH", "TEWA", "TOSO", "TRES", "VATH",
    "VEER", "VESP", "WAVI", "WBNU", "WCSP", "WETA", "WEWP", "WISN",
    "WIWA", "WIWR", "WTSP", "WWCR", "YBFL", "YBSA", "YEWA", "YHBL",
    "YRWA")

distr <- list()

for (spp in SPP) {
   fl <- list.files(file.path(OUTDIR, "All", spp))
   d <- numeric(0)
   for (i in 1:length(fl)) {
        cat("Combine -", sect, spp, i, "/", length(fl), "\n");flush.console()
        e <- new.env()
        load(file.path(OUTDIR, sect, spp, fl[i]), envir=e)
        d <- c(d, rowSums(e$pxNcrS))
    }
    distr[[spp]] <- d
}
save(distr, file=file.path(OUTDIR, "pixel-level-values.RData"))
