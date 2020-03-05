## Sector effects

rm(list=ls())

PROJ <- "north"

#STAGE <- "Space"
STAGE <- "HF"

## SEC describe what HF was not backfilled
## can be: All, None, Energy
#SEC <- "All"
#SEC <- "None"
#SEC <- "Agr"
#SEC <- "Transp"
#SEC <- "Energy"
#SEC <- "EnS"
#SEC <- "EnH"
#SEC <- "Urban"
#SEC <- "For"

sectors_list <- list(
    Agr=c("Crop", "RoughP", "TameP"),
    Transp=c("HardLin", "TrSoftLin"),
    EnSoft=c("EnSoftLin", "Seismic"),
    EnHard=c("Mine", "Well", "Industrial"),
    Urb=c("Rural", "Urban"),
    For=c("ForHarv"))


## max boot runs or NULL
BMAX <- 1#NULL

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

library(mefa4)
library(intrval)
library(raster)
library(rgdal)
library(sp)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
Xn <- get_model_matrix(en$DAT, en$mods)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

cfn <- list(
    hab=c("vegcBSpr", "vegcCrop", "vegcGraminoidFen",
        "vegcGrassHerb", "vegcIndustrial", "vegcLarch", "vegcMarsh",
        "vegcMine", "vegcMixedwood", "vegcPine", "vegcRoughP", "vegcRural",
        "vegcShrub", "vegcSpruce", "vegcSwamp", "vegcTameP", "vegcUrban",
        "wtAge", "wtAge2", "wtAge05", "fCC2",
        "isCon:wtAge", "isCon:wtAge2", "isUpCon:wtAge", "isBSLarch:wtAge",
        "isUpCon:wtAge2", "isBSLarch:wtAge2", "isMix:wtAge", "isPine:wtAge",
        "isWSpruce:wtAge", "isMix:wtAge2", "isPine:wtAge2", "isWSpruce:wtAge2",
        "isCon:wtAge05", "isUpCon:wtAge05", "isBSLarch:wtAge05", "isMix:wtAge05",
        "isPine:wtAge05", "isWSpruce:wtAge05"),
    modif=c("mWell", "mSoft", "mEnSft", "mTrSft", "mSeism"),
    nuisance=c("ROAD", "CMETHODSM", "CMETHODRF"),
    spclim=c("pWater_KM", "pWater2_KM", "xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
        "xMCMT", "xY", "xX", "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY"),
    ssh=c("SSH_KM", "SSH05_KM", "THF_KM",
        "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", "Cult_KM",
        "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", "Noncult2_KM"),
    yr=c("YR"))
pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
    "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05,
    "CMETHODSM"=1, "CMETHODRF"=1)
setdiff(colnames(Xage), cfn$hab) # should be intercept only

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y


xclim <- data.frame(
    transform_clim(kgrid),
    pAspen=kgrid$pAspen,
    pWater_KM=kgrid$pWater,
    pWater2_KM=kgrid$pWater^2)
## this has pAspen for the south, otherwise all the same
Xclim <- model.matrix(as.formula(paste0("~-1+", paste(cfn$spclim, collapse="+"))), xclim)
colnames(Xclim) <- fix_names(colnames(Xclim))

## ch2soil ch2veg trSoil trVeg
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
tv$ao <- as.factor(paste0(as.character(tv[, "UseInAnalysisFine"]), ifelse(tv[, "MatureOld"], "O", "")))

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2veg$rf2 <- tv$UseInAnalysisFineAge[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$UseInAnalysisFineAge[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]
ch2veg$rf3 <- tv$ao[match(ch2veg$rf, rownames(tv))]
ch2veg$cr3 <- tv$ao[match(ch2veg$cr, rownames(tv))]

ch2veg$sector <- as.character(ch2veg$sector)
ch2veg$sector2 <- as.character(ch2veg$sector)
ch2veg$sector2[ch2veg$cr3 %in% sectors_list$EnSoft] <- "EnSoft"
ch2veg$sector2[ch2veg$cr3 %in% sectors_list$EnHard] <- "EnHard"
## note: predictions treated Rural/Industrial as Energy and not RuralUrban, peat mine also Energy
ch2veg$sector[ch2veg$sector=="Misc" & ch2veg$sector2=="EnHard"] <- "Energy"
ch2veg$sector[ch2veg$sector=="RuralUrban" & ch2veg$sector2=="EnHard"] <- "Energy"
table(ch2veg$sector2, ch2veg$sector)


EXCL <- c("HWater", "SoilUnknown", "SoilWater", "Water")
MODIF <- c("SoftLin", "Well", "EnSoftLin", "TrSoftLin", "Seismic")

keepn <- rownames(ch2veg)[!(ch2veg$cr2 %in% EXCL) & !(ch2veg$rf2 %in% EXCL)]
trVeg <- trVeg[,keepn]
ch2veg <- ch2veg[keepn,]
ch2veg$modif <- ch2veg$cr2 %in% MODIF
rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn

CN <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")

## define OSR
ogrListLayers("d:/spatial/Oilsands-Boundaries.gdb")
pl <- readOGR("d:/spatial/Oilsands-Boundaries.gdb", "OilsandRegionDissolve10TM")
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(pl))
o <- over(xy, pl)
#plot(xy, pch=".", col=ifelse(is.na(o$FIELDCODE), 1, 4))

ss <- !is.na(o$FIELDCODE)
pveg_mine <- 0.44 # proportion of vegetated mines in OSR from NDVI

trVegSS <- trVeg[ss,]
AVegSS <- colSums(trVegSS)

SECS <- c("All", "Agr", "Transp", "Energy", "EnS", "EnH", "Urban", "For")
for (SEC in SECS) {

    ## BF_THIS is the list of HF that is to be backfilled
    ## sector that is NOT being backfilled
    BF_THIS <- unname(unlist(sectors_list))
    if (SEC == "All") {
        ## keep all the HF
        BF_THIS <- character(0) # all current
    }
    if (SEC == "None") {
        BF_THIS <- unname(unlist(sectors_list)) # backfill all
    }
    if (SEC == "Agr") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% sectors_list$Agr)]
    }
    if (SEC == "Transp") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% sectors_list$Transp)]
    }
    if (SEC == "Energy") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% c(sectors_list$EnSoft, sectors_list$EnHard))]
    }
    if (SEC == "EnS") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% sectors_list$EnSoft)]
    }
    if (SEC == "EnH") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% sectors_list$EnHard)]
    }
    if (SEC == "Urban") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% sectors_list$Urb)]
    }
    if (SEC == "For") {
        BF_THIS <- BF_THIS[!(BF_THIS %in% sectors_list$For)]
    }


    ## define SSH based on actuall partial backfill here based on trVeg

    SSH <- en$SSH
    compare_sets(colnames(en$SSH), levels(ch2veg$cr3))
    setdiff(colnames(en$SSH), levels(ch2veg$cr3))
    SSH_EXCL <- setdiff(levels(ch2veg$cr3), colnames(en$SSH)) # do not count for SSH
    ZERO <- c("Bare", "SnowIce", "HWater", "Water") # non-habitat
    ALL_HF <- c("Crop", "RoughP", "TameP",
        "HardLin", "TrSoftLin",
        "EnSoftLin", "Seismic",
        "Mine", "Well",
        "Industrial", "Rural", "Urban",
        "ForHarv") # backfill these
    LIN_HF <- c(    "HardLin", "TrSoftLin", "EnSoftLin", "Seismic")
    Alien_HF <- c("Crop", "TameP", #"RoughP",
        "HardLin", "Mine", "Well", "Industrial", "Rural", "Urban")

    ## fully backfilled
    bf0 <- as.character(ch2veg$rf3)
    SSH0 <- row_std(groupSums(trVeg[ss,], 2, bf0))

    ## this does not separate forestry --> for SSH_KM variable
    pbf <- as.character(ch2veg$cr3)
    bfi <- pbf %in% BF_THIS
    pbf[bfi] <- as.character(ch2veg$rf3[bfi])
    SSH <- row_std(groupSums(trVeg[ss,], 2, pbf))
    ## this does separate forestry --> for HF_KM variables
    pbfFor <- as.character(ch2veg$cr3)
    pbfFor[ch2veg$sector == "Forestry"] <- "ForHarv"
    pbfFor[bfi] <- as.character(ch2veg$rf3[bfi])
    SSHfor <- row_std(groupSums(trVeg[ss,], 2, pbfFor))

    ## placeholder matrix
    dd <- data.frame(SSH_KM = rep(0, sum(ss)))
    dd$SSH05_KM <- 0
    dd$THF_KM <- rowSums(SSHfor[,colnames(SSHfor) %in% ALL_HF,drop=FALSE])
    dd$Lin_KM <- rowSums(SSHfor[,colnames(SSHfor) %in% LIN_HF,drop=FALSE])
    ## note: no abandoned or rough pasture here
    dd$Cult_KM <- rowSums(SSHfor[,colnames(SSHfor) %in% c("Crop", "TameP"),drop=FALSE])
    dd$Alien_KM <- rowSums(SSHfor[,colnames(SSHfor) %in% Alien_HF,drop=FALSE])

    dd$Nonlin_KM <- dd$THF_KM - dd$Lin_KM
    dd$Noncult_KM <- dd$THF_KM - dd$Cult_KM
    dd$Succ_KM <- dd$THF_KM - dd$Alien_KM
    dd$THF2_KM <- dd$THF_KM^2
    dd$Succ2_KM <- dd$Succ_KM^2
    dd$Alien2_KM <- dd$Alien_KM^2
    dd$Noncult2_KM <- dd$Noncult_KM^2
    dd$Nonlin2_KM <- dd$Nonlin_KM^2

    Xssh <- model.matrix(as.formula(paste0("~-1+", paste(cfn$ssh, collapse="+"))), dd)
    colnames(Xssh) <- fix_names(colnames(Xssh))
    Xssh0 <- Xssh
    Xssh0[] <- 0

    for (spp in SPP) {
        cat(SEC, spp)

        resn <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))

        ## north estimates
        names(en$mods)
        ESTN <- suppressWarnings(get_coef(resn, Xn, stage=STAGE, na.out=FALSE))

        b <- nrow(ESTN)
        if (!is.null(BMAX))
            b <- min(BMAX, b)

        CR <- matrix(0, sum(ss), b)
        rownames(CR) <- rownames(kgrid)[ss]
        RF <- CR
        HABCR <- matrix(0, nrow(ch2veg), b)
        rownames(HABCR) <- rownames(ch2veg)
        HABRF <- HABCR

        for (i in seq_len(b)) {
            ir <- i %% round(b/10)
            if (is.na(ir) || ir == 0) {
                cat(".")
                flush.console()
            }
            estn <- ESTN[i,]

            ## surrounding SSH and HF
            ssh <- resn[[as.integer(i)]]$ssh
            Xssh0[,"SSH_KM"] <- rowSums(SSH0[,colnames(SSH0) %in% ssh$labels,drop=FALSE]) # reference
            Xssh0[,"SSH05_KM"] <- sqrt(Xssh0[,"SSH_KM"])
            Xssh[,"SSH_KM"] <- rowSums(SSH[,colnames(SSH) %in% ssh$labels,drop=FALSE])
            Xssh[,"SSH05_KM"] <- sqrt(Xssh[,"SSH_KM"])
            munSsh0 <- drop(Xssh0[,cfn$ssh] %*% estn[cfn$ssh])
            munSsh <- drop(Xssh[,cfn$ssh] %*% estn[cfn$ssh])
            munCl <- drop(Xclim[ss,cfn$spclim] %*% estn[cfn$spclim])

            ## let climate include SSH and HF too
            munClim0 <- munCl + munSsh0 # reference (fully backfilled)
            munClim <- munCl + munSsh # current/partial backfilled

            ## sof linear modifiers & habitat
            if (estn["mSoft"] != 0 & estn["mEnSft"] == 0) {
                estn["mEnSft"] <- estn["mSoft"]
                estn["mTrSft"] <- estn["mSoft"]
                estn["mSeism"] <- estn["mSoft"]
            }
            munHab <- drop(Xage %*% estn[colnames(Xage)])
            munMod <- structure(numeric(length(cfn$modif)), names=cfn$modif)
            for (k in cfn$modif)
                munMod[k] <- log(linexp(1, estn[k], pm[k]))
            munHab <- c(munHab,
                HardLin=-1000,
                Bare=-1000,
                SnowIce=-1000,
                Well=unname(munMod["mWell"]),
                EnSoftLin=unname(munMod["mEnSft"]),
                TrSoftLin=unname(munMod["mTrSft"]),
                Seismic=unname(munMod["mSeism"]))
            #munHab["Mine"] <- -1000
            ## expand coefficients for north
            prnCr <- munHab[match(ch2veg$cr2, names(munHab))]
            prnRf <- munHab[match(ch2veg$rf2, names(munHab))]
            prnCr[ch2veg$modif] <- prnRf[ch2veg$modif] + prnCr[ch2veg$modif]

            ## put pieces together for north
            ## multiplying with row normalized area gives the weighted averaging
            ADnCr <- 100 * t(exp(prnCr) * t(trVegSS)) * exp(munClim) # males / km cell
            ADnRf <- 100 * t(exp(prnRf) * t(trVegSS)) * exp(munClim0)
            #ADnCrHab <- groupSums(ADnCr, 2, ch2veg$cr2)

            ## Mines: 70 is vegetated
            ADnCr[,ch2veg$cr == "MineSite"] <- pveg_mine * ADnCr[,ch2veg$cr == "MineSite"]

            ## quantiles not applied -- look at that post hoc before summing up
            CR[,i] <- rowSums(ADnCr) # no pair adjustment applied, just ha to km
            RF[,i] <- rowSums(ADnRf) # no pair adjustment applied, just ha to km
            if (i == 1L) {
                SCR <- groupSums(ADnCr, 2, ch2veg$sector2)
                SRF <- groupSums(ADnRf, 2, ch2veg$sector2)
            }
            HABCR[colnames(ADnCr),i] <- colSums(ADnCr)
            HABRF[colnames(ADnRf),i] <- colSums(ADnRf)
        }
        CR <- apply(CR, 1, median)
        RF <- apply(RF, 1, median)
        DIRO <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/sector/", spp)
        if (!dir.exists(DIRO))
            dir.create(DIRO)
        save(
            SCR, SRF,
            file=paste0(DIRO, "/", spp, "-", toupper(STAGE), "-", toupper(SEC), "-pixel.RData"))
    #    save(
    #        CR, RF,
    #        HABRF, HABCR,
    #        file=paste0(DIRO, "/", spp, "-", toupper(STAGE), "-", toupper(SEC), ".RData"))
        cat("OK\n")
    }

}

