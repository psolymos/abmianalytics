library(mefa4)
#library(intrval)
library(raster)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
#ROOT <- "~/GoogleWork/tmp"

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

#xwalk <- read.csv("~/repos/abmianalytics/lookup/veg-v61-crosswalk.csv")
#rownames(xwalk) <- xwalk[,1]

cfs <- list(
    hab=c("soilcClay", "soilcCrop", "soilcRapidDrain",
        "soilcRoughP", "soilcSaline", "soilcTameP", "soilcUrbInd"),
    modif=c("mWell", "mSoft"),
    nuisance=c("ROAD", "CMETHODSM", "CMETHODRF"),
    spclim=c("pAspen", "pWater_KM",
        "pWater2_KM", "xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
        "xMCMT", "xY", "xX", "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY"),
    ssh=c("SSH_KM", "SSH05_KM", "THF_KM",
        "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", "Cult_KM",
        "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", "Noncult2_KM"),
    yr=c("YR"))

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

load("d:/abmi/AB_data_v2019/misc/overlap/OverlapReg.RData")
rownames(OverlapReg) <- OverlapReg$LinkID
OverlapReg$pAspen <- kgrid[rownames(OverlapReg), "pAspen"]
OverlapReg$wN <- OverlapReg$pAspen / (OverlapReg$pAspen + (1-OverlapReg$pForest))
kgrid$wN <- ifelse(kgrid$NRNAME == "Grassland", 0, 1)
kgrid[rownames(OverlapReg), "wN"] <- OverlapReg$wN

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth | Lookup$ModelSouth,])
rownames(tax) <- tax$Code

xclim <- data.frame(
    transform_clim(kgrid),
    pAspen=kgrid$pAspen,
    pWater_KM=kgrid$pWater,
    pWater2_KM=kgrid$pWater^2)
## this has pAspen for the south, otherwise all the same
Xclim <- model.matrix(as.formula(paste0("~-1+", paste(cfs$spclim, collapse="+"))), xclim)
colnames(Xclim) <- fix_names(colnames(Xclim))

## ch2soil ch2veg trSoil trVeg
#load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v61hf2016v3WildFireUpTo2016.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]

compare_sets(ch2soil$cr, rownames(ts))
setdiff(ch2soil$cr, rownames(ts))
setdiff(rownames(ts), ch2soil$cr)

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2soil$rf2 <- ts$UseInAnalysisCoarse[match(ch2soil$rf, rownames(ts))]
ch2soil$cr2 <- ts$UseInAnalysisCoarse[match(ch2soil$cr, rownames(ts))]
ch2soil$sector <- ts$Sector61[match(ch2soil$cr, rownames(ts))]
ch2veg$rf2 <- tv$UseInAnalysisFineAge[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$UseInAnalysisFineAge[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

str(ch2soil)
str(ch2veg)

## checking if sectors are assigned correctly
if (FALSE) {
isN <- kgrid$NRNAME != "Grassland" & kgrid$NRNAME != "Rocky Mountain" &
    kgrid$NRNAME != "Parkland" & kgrid$NSRNAME != "Dry Mixedwood"
isS <- kgrid$NRNAME == "Grassland" | kgrid$NRNAME == "Parkland" |
    kgrid$NSRNAME == "Dry Mixedwood"

load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3noDistVeg.Rdata")
s1 <- groupSums(dd_kgrid[[1]], 2, tv$Sector61[match(colnames(dd_kgrid[[1]]), rownames(tv))])
s2 <- groupSums(dd_kgrid[[3]], 2, ts$Sector61[match(colnames(dd_kgrid[[3]]), rownames(ts))])

s3 <- groupSums(trVeg, 2, ch2veg$sector[match(colnames(trVeg), rownames(ch2veg))])
s4 <- groupSums(trSoil, 2, ch2soil$sector[match(colnames(trSoil), rownames(ch2soil))])

load("d:/abmi/AB_data_v2018/data/analysis/checks/Veg transtion and sector lookup.RData")

ch2veg$VegTSL <- VegTSL[rownames(ch2veg),"Sector"]
table(Peter=ch2veg$sector, Ermias=ch2veg$VegTSL)

aa <- data.frame(VegKg=100*colSums(s1[isN,])[levels(tv$Sector61)]/sum(s1[isN,]),
    VegTr=100*colSums(s3[isN,])[levels(tv$Sector61)]/sum(s3[isN,]),
    SoilKg=100*colSums(s2[isS,])[levels(tv$Sector61)]/sum(s2[isS,]),
    SoilTr=100*colSums(s4[isS,])[levels(tv$Sector61)]/sum(s4[isS,]))
colSums(aa)
round(aa, 3)

}



EXCL <- c("HWater", "SoilUnknown", "SoilWater", "Water")
MODIF <- c("SoftLin", "Well", "EnSoftLin", "TrSoftLin", "Seismic")

keeps <- rownames(ch2soil)[!(ch2soil$cr2 %in% EXCL) & !(ch2soil$rf2 %in% EXCL) ]
keepn <- rownames(ch2veg)[!(ch2veg$cr2 %in% EXCL) & !(ch2veg$rf2 %in% EXCL)]
trSoil <- trSoil[,keeps]
trVeg <- trVeg[,keepn]
ch2soil <- ch2soil[keeps,]
ch2veg <- ch2veg[keepn,]
ch2soil$modif <- ch2soil$cr2 %in% MODIF
ch2veg$modif <- ch2veg$cr2 %in% MODIF
rss <- rowSums(trSoil)
rss[rss==0] <- 1
trSoil <- trSoil / rss
rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn

summary(ch2soil)
summary(ch2veg)

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

CN <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")
stemp <- matrix(0, nrow(kgrid), nlevels(ch2veg$sector))
dimnames(stemp) <- list(rownames(kgrid), CN)
stemp <- as(stemp, "dgCMatrix")

## spp runs

#spp <- "ALFL"
#SPP <- rownames(tax)

bbb <- read.csv("s:/sppweb2018/c4i/tables/StandardizedOutput-birds-final-lookup-withChecks.csv")
bbb <- bbb[is.na(bbb$Exclude),]
bbb <- bbb[bbb$ModelNorth | bbb$ModelSouth,]
rownames(bbb) <- bbb$Code.1

SPP <- rownames(bbb)
DONE <- gsub(".RData", "", list.files("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-09-20/"))
SPP <- SPP[!(SPP %in% DONE)]

b <- 100

#SPP <- SPP[1:30]
#SPP <- SPP[31:60]
#SPP <- SPP[61:90]
#SPP <- SPP[91:126]

for (spp in SPP) {
    gc()
    cat(spp)
    flush.console()
    CURRB <- REFB <- matrix(0, nrow(kgrid), b)
    rownames(CURRB) <- rownames(REFB) <- rownames(kgrid)
    RUN_OK <- logical(b) # FALSE

    TYPE <- "C" # combo
    if (bbb[spp, "ModelSouth"] && !bbb[spp, "ModelNorth"])
        TYPE <- "S"
    if (!bbb[spp, "ModelSouth"] && bbb[spp, "ModelNorth"])
        TYPE <- "N"

    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))

    #i <- 1 # boot run
    for (i in seq_len(b)) {
        gc()
        if (i %% 5 == 0) {
            cat(".")
            flush.console()
        }

        if (TYPE != "N") {
            #ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
            ## south estimates
            ests <- suppressWarnings(get_coef(ress, Xs, stage="Space", na.out=FALSE))[i,]
            musClim <- drop(Xclim[,cfs$spclim] %*% ests[cfs$spclim])
            musHab <- c(ests[1], ests[1]+ests[cfs$hab])
            names(musHab) <- c("Productive", gsub("soilc", "", cfs$hab))
            musMod <- structure(numeric(length(cfs$modif)), names=cfs$modif)
            for (k in cfs$modif)
                musMod[k] <- log(linexp(1, ests[k], pm[k]))
            musHab <- c(musHab,
                HardLin=-1000,
                HFor=unname(musHab["Productive"]),
                SoftLin=unname(musMod["mSoft"]))
            if (all(is.finite(exp(musHab)))) {
                ## expand coefficients for south
                prsCr <- musHab[match(ch2soil$cr2, names(musHab))]
                prsRf <- musHab[match(ch2soil$rf2, names(musHab))]
                prsCr[ch2soil$modif] <- prsRf[ch2soil$modif] + prsCr[ch2soil$modif]
                ## put pieces together for south
                ADsCr <- t(exp(prsCr) * t(trSoil)) * exp(musClim)
                ADsRf <- t(exp(prsRf) * t(trSoil)) * exp(musClim)
                ## add up by sector for south
                ADsCrSect <- groupSums(ADsCr, 2, ch2soil$sector)
                ADsRfSect <- groupSums(ADsRf, 2, ch2soil$sector)
                RUN_OK[i] <- TRUE
            }
        } else {
            ADsCrSect <- stemp
            ADsRfSect <- stemp
        }

        if (TYPE != "S") {
            #resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
            ## north estimates
            estn <- suppressWarnings(get_coef(resn, Xn, stage="Space", na.out=FALSE))[i,]
            munClim <- drop(Xclim[,cfn$spclim] %*% estn[cfn$spclim])
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
            munHab["Mine"] <- -1000
            if (all(is.finite(exp(munHab)))) {
                ## expand coefficients for north
                prnCr <- munHab[match(ch2veg$cr2, names(munHab))]
                prnRf <- munHab[match(ch2veg$rf2, names(munHab))]
                prnCr[ch2veg$modif] <- prnRf[ch2veg$modif] + prnCr[ch2veg$modif]
                ## put pieces together for north
                ## multiplying with row normalized area gives the weighted averaging
                ADnCr <- t(exp(prnCr) * t(trVeg)) * exp(munClim)
                ADnRf <- t(exp(prnRf) * t(trVeg)) * exp(munClim)
                ## add up by sector for north
                ADnCrSect <- groupSums(ADnCr, 2, ch2veg$sector)
                ADnRfSect <- groupSums(ADnRf, 2, ch2veg$sector)
                RUN_OK[i] <- TRUE
            }
        } else {
            ADnCrSect <- stemp
            ADnRfSect <- stemp
        }

        if (RUN_OK[i]) {
            ## weighted average
            #wS <- 1-kgrid$pAspen
            wS <- 1-kgrid$wN
            if (TYPE == "S")
                wS[] <- 1
            if (TYPE == "N")
                wS[] <- 0
            wS[kgrid$useS] <- 1
            wS[kgrid$useN] <- 0
            #wS <- wS * kgrid$pSoil
            Curr <- wS * ADsCrSect[,CN] + (1-wS) * ADnCrSect[,CN]
            Ref <- wS * ADsRfSect[,CN] + (1-wS) * ADnRfSect[,CN]

            if (i == 1) {
                CURR1 <- Curr
                REF1 <- Ref
                CURR <- Curr
                REF <- Ref
            } else {
                ## this updates the numbers for eventually give the mean
                px <- sum(RUN_OK-1) / sum(RUN_OK)
                CURR <- CURR * px + Curr * (1-px)
                REF <- REF * px + Ref * (1-px)
            }
            CURRB[,i] <- rowSums(Curr)
            REFB[,i] <- rowSums(Ref)
        } else {
            CURRB[,i] <- NA
            REFB[,i] <- NA
        }
    }

#    save(CURR, REF, CURR1, REF1, CURRB, REFB,
#        file=paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-09-20/", spp, ".RData"))
    save(CURR, REF, CURR1, REF1, CURRB, REFB,
        file=paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/tmp/", spp, ".RData"))

    cat("DONE\n")
}


## north models with bootstrap to get pop sizes in BCR6
DONE <- gsub(".RData", "", list.files("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/"))
SPP <- setdiff(rownames(tax)[tax$ModelNorth & rownames(tax) %in% colnames(en$OFF)], DONE)
ss <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"

for (spp in SPP) {
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    ## north estimates
    ESTN <- suppressWarnings(get_coef(resn, Xn, stage="Space", na.out=FALSE))
    b <- nrow(ESTN)
    cat(spp, "with", b, "runs ")
    flush.console()
    CR <- matrix(0, sum(ss), b)
    rownames(CR) <- rownames(kgrid)[ss]
    HAB <- matrix(0, nlevels(ch2veg$cr2), b)
    rownames(HAB) <- levels(ch2veg$cr2)

    for (i in seq_len(b)) {
        if (i %% round(b/10) == 0) {
            cat(".")
            flush.console()
        }
        estn <- ESTN[i,]
        munClim <- drop(Xclim[ss,cfn$spclim] %*% estn[cfn$spclim])
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
        munHab["Mine"] <- -1000
        ## expand coefficients for north
        prnCr <- munHab[match(ch2veg$cr2, names(munHab))]
        prnRf <- munHab[match(ch2veg$rf2, names(munHab))]
        prnCr[ch2veg$modif] <- prnRf[ch2veg$modif] + prnCr[ch2veg$modif]
        ## put pieces together for north
        ## multiplying with row normalized area gives the weighted averaging
        ADnCr <- 100 * t(exp(prnCr) * t(trVeg[ss,])) * exp(munClim) # males / km cell
#        ADnRf <- t(exp(prnRf) * t(trVeg)) * exp(munClim)
        ADnCrHab <- groupSums(ADnCr, 2, ch2veg$cr2)

        ## quantiles not applied -- look at that post hoc before summing up
        CR[,i] <- rowSums(ADnCrHab) # no pair adjustment applied, just ha to km
        HAB[colnames(ADnCrHab),i] <- colSums(ADnCrHab)

    }
    save(CR, HAB,
        file=paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/", spp, ".RData"))
    cat(" OK\n")

}


## mapping

library(mefa4)
library(raster)

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth | Lookup$ModelSouth,])
rownames(tax) <- tax$Code

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}


Rmaskn <- make_raster(as.integer(1-kgrid$useS), kgrid, rt)
values(Rmaskn)[values(Rmaskn) == 0] <- NA
Rmasks <- make_raster(as.integer(1-kgrid$useN), kgrid, rt)
values(Rmasks)[values(Rmasks) == 0] <- NA
#Rmaskm <- make_raster(as.integer(kgrid$NRNAME == "Rocky Mountain"), kgrid, rt)
#values(Rmaskm)[values(Rmaskm) == 0] <- NA
Rw <- make_raster(as.integer(kgrid$pWater > 0.99), kgrid, rt)
values(Rw)[values(Rw) == 0] <- NA

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
col3 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221"))(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude

## checking results
if (FALSE) {
cn <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")
rn <- rownames(kgrid)[kgrid$NRNAME != "Grassland"]
Aveg <- groupSums(trVeg, 2, ch2veg$sector)[,cn]
KA <- 100*colMeans(Aveg[rn,])
#Asoil <- 100*colMeans(groupSums(trSoil, 2, ch2soil$sector)[rn,cn])

#spp <- "BTNW" # species
#load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-04-01/", spp, ".RData"))


CS <- colSums(Curr[rn,cn])
RS <- colSums(Ref[rn,cn])
NC <- sum(CS)
NR <- sum(RS)
Sector_Total <- (100 * (CS - RS) / NR)[-1]
Sector_UnderHF <- (100 * (CS - RS) / RS)[-1]
Sector_Area <- (100 * KA / sum(KA))[names(Sector_Total)]
Sector_Unit <- 100 * Sector_Total / Sector_Area

cat(spp, "\n")
round(cbind(CS=CS, RS=RS, Df=CS-RS), 2)
round(cbind(Total=Sector_Total,
    Under=Sector_UnderHF,
    Unit=Sector_Unit), 2)


Dcr <- rowSums(Curr)
q <- quantile(Dcr, 0.99)
Dcr[Dcr > q] <- q
summary(Dcr)
Drf <- rowSums(Ref)
q <- quantile(Drf, 0.99)
Drf[Drf > q] <- q
summary(Drf)
MAX <- max(Dcr, Drf)

Rcr <- make_raster(Dcr, kgrid, rt)
Rrf <- make_raster(Drf, kgrid, rt)
plot(sqrt(Rcr), col=col1)
plot(sqrt(Rrf), col=col1)

x <- make_raster(as.numeric(Aveg[,"Forestry"]), kgrid, rt)
plot(x, col=col1)

}



## detections
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)

ddd <- ee$dd
## subset based on analysis data -- same filtering applied
ddd <- ddd[unique(c(rownames(en$DAT), rownames(es$DAT))),]
ddd <- nonDuplicated(ddd, ddd$SS, TRUE)
yyy <- groupSums(ee$yy, 1, ee$dd$SS)[rownames(ddd),]
yyy[yyy > 0] <- 1
ss <- !is.na(ddd$X) & !is.na(ddd$NRNAME)
ddd <- ddd[ss,]
yyy <- yyy[ss,]

xy <- SpatialPoints(as.matrix(ddd[,c("X","Y")]))
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))
rt10 <- aggregate(rt, fact=10)
sam0 <- rasterize(xy, rt10, field=1, fun='sum')
values(sam0)[!is.na(values(sam0))] <- 1

rnr <- make_raster(as.integer(kgrid$NRNAME), kgrid, rt)
cnr <- c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae')
cnr <- cnr[c(5,6,1,2,4,3)]

## 10k level detections
if (FALSE) {
xyall <- xyFromCell(sam0, c(1:ncell(sam0)), spatial=TRUE)
tmp <- spTransform(xyall, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xyall <- as(xyall, "SpatialPointsDataFrame")
xyall@data <- data.frame(coordinates(tmp), surveyed=ifelse(is.na(values(sam0)), 0, 1))
#plot(xyall, col=xyall@data$surveyed+1, pch=".")
xyall2 <- xyall
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    ## only non ABMI detections here
    xy1 <- try(SpatialPoints(as.matrix(ddd[yyy[,spp] > 0 &
        !(ddd$PCODE %in% c("BU_ABMI", "BU_OG-ABMI", "ABMIRF", "ABMISM")) &
        startsWith(as.character(ddd$PCODE), "BBS"), c("X","Y")])))
    if (!inherits(xy1, "try-error")) {
        proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        xy1 <- spTransform(xy1, proj4string(rt))
        sam1 <- rasterize(xy1, rt10, field=1, fun='last')
        xyall@data[[as.character(tax[spp, "SpeciesID"])]] <- ifelse(is.na(values(sam1)), 0, 1)
    } else {
        xyall@data[[as.character(tax[spp, "SpeciesID"])]] <- 0
    }

    xy1 <- try(SpatialPoints(as.matrix(ddd[yyy[,spp] > 0 &
            !(ddd$PCODE %in% c("BU_ABMI", "BU_OG-ABMI", "ABMIRF", "ABMISM")) &
            !startsWith(as.character(ddd$PCODE), "BBS"), c("X","Y")])))
    if (!inherits(xy1, "try-error")) {
        proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        xy1 <- spTransform(xy1, proj4string(rt))
        sam1 <- rasterize(xy1, rt10, field=1, fun='last')
        xyall2@data[[as.character(tax[spp, "SpeciesID"])]] <- ifelse(is.na(values(sam1)), 0, 1)
    } else {
        xyall2@data[[as.character(tax[spp, "SpeciesID"])]] <- 0
    }
}
#summary(xyall@data)
xyall <- xyall@data # BBS
xyall <- xyall[xyall$surveyed == 1,]
xyall <- xyall[rowSums(xyall[,-(1:3)]) > 0,]
xyall2 <- xyall2@data # BAM
xyall2 <- xyall2[xyall2$surveyed == 1,]
xyall2 <- xyall2[rowSums(xyall2[,-(1:3)]) > 0,]
dim(xyall)
dim(xyall2)

write.csv(xyall[,-3], row.names=FALSE, file="d:/abmi/AB_data_v2018/www/detections10k_BBS.csv")
write.csv(xyall2[,-3], row.names=FALSE, file="d:/abmi/AB_data_v2018/www/detections10k_BAM.csv")


#dbWriteTable(con, "detections10k_nonABMI", xyall, overwrite=TRUE, row.names=FALSE) # birds
dbWriteTable(con, "detections10k_BBS", xyall, overwrite=TRUE, row.names=FALSE) # birds
dbWriteTable(con, "detections10k_BAM", xyall2, overwrite=TRUE, row.names=FALSE) # birds
}

library(DBI)
source("~/.ssh/postgres")
data.frame(..postgres_science)
con <- dbConnect(
    odbc::odbc(),
    driver   = "PostgreSQL Unicode(x64)",
    server   = ..postgres_science$host,
    database = ..postgres_science$database,
    uid      = ..postgres_science$username,
    pwd      = ..postgres_science$password,
    port     = ..postgres_science$port)
(dbl <- dbListTables(con))

#tmp <- data.frame(x=1:3, y=5:7)
#dbWriteTable(con, "test", tmp)
#dbSendQuery(con, "drop table test")

#d0 <- dbReadTable(con, "map_YellowheadedBlackbird")
#dbSendQuery(con, "DROP TABLE map_YellowheadedBlackbird;")


#dbDisconnect(con)

## saving BOOT objects for c4i

NAfound <- NULL
for (spp in rownames(bbb)) {

    cat(spp, "\n");flush.console()
    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-09-20/", spp, ".RData"))

    TYPE <- "C" # combo
    if (bbb[spp, "ModelSouth"] && !bbb[spp, "ModelNorth"])
        TYPE <- "S"
    if (!bbb[spp, "ModelSouth"] && bbb[spp, "ModelNorth"])
        TYPE <- "N"

    for (i in 1:100) {
        if (any(is.na(CURRB[,i]))) {
            CURRB[,i] <- CURRB[,i-1]
            cat(" NAcr")
        }
        if (any(is.na(REFB[,i]))) {
            REFB[,i] <- REFB[,i-1]
            cat(" NArf")
        }

        q <- quantile(CURRB[,i], 0.99, na.rm=TRUE)
        CURRB[CURRB[,i] > q,i] <- q
        q <- quantile(REFB[,i], 0.99, na.rm=TRUE)
        REFB[REFB[,i] > q,i] <- q
    }
#    Curr.Boot <- as.matrix(groupSums(CURRB, 1, kgrid[rownames(CURRB),"Row10_Col10"]))
#    Ref.Boot <- as.matrix(groupSums(REFB, 1, kgrid[rownames(CURRB),"Row10_Col10"]))
    Curr.Boot <- as.matrix(groupMeans(CURRB, 1, kgrid[rownames(CURRB),"Row10_Col10"]))
    Ref.Boot <- as.matrix(groupMeans(REFB, 1, kgrid[rownames(CURRB),"Row10_Col10"]))
    CN <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")
    SA.Curr <- as.matrix(CURR1[,CN])
    SA.Ref <- as.matrix(REF1[,CN])
    sum(SA.Curr)
    summary(colSums(Curr.Boot)*100)

    if (any(is.na(Curr.Boot)) || any(is.na(Ref.Boot)) || any(is.na(SA.Curr)) || any(is.na(SA.Ref))) {
        cat("\t--> NA found for", spp, "\n\n")
        NAfound <- c(NAfound, spp)
    } else {
        save(Curr.Boot, Ref.Boot, file=paste0("d:/abmi/reports/2018/results/birds/boot/",
            as.character(bbb[spp, "SpeciesID"]), ".RData"))
        save(SA.Curr, SA.Ref,     file=paste0("d:/abmi/reports/2018/results/birds/sector/",
            as.character(bbb[spp, "SpeciesID"]), ".RData"))
    }
}




PLOT <- TRUE
SAVE <- FALSE
for (spp in rownames(tax)[48:173]) {

    cat(spp, "\n");flush.console()
#    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-04-01/", spp, ".RData"))
    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-05-14/", spp, ".RData"))

    TYPE <- "C" # combo
    if (tax[spp, "ModelSouth"] && !tax[spp, "ModelNorth"])
        TYPE <- "S"
    if (!tax[spp, "ModelSouth"] && tax[spp, "ModelNorth"])
        TYPE <- "N"

    for (i in 1:ncol(CURRB)) {
        q <- quantile(CURRB[,i], 0.99, na.rm=TRUE)
        CURRB[CURRB[,i] > q,i] <- q
        q <- quantile(REFB[,i], 0.99, na.rm=TRUE)
        REFB[REFB[,i] > q,i] <- q
    }
    Dcr <- rowMeans(CURRB, na.rm=TRUE)
    Drf <- rowMeans(REFB, na.rm=TRUE)
    SD <- apply(CURRB, 1, sd, na.rm=TRUE)

#    Dcr <- rowSums(Curr)
#    q <- quantile(Dcr, 0.99)
#    Dcr[Dcr > q] <- q
#    summary(Dcr)
#    Drf <- rowSums(Ref)
#    q <- quantile(Drf, 0.99)
#    Drf[Drf > q] <- q
#    summary(Drf)
    MAX <- max(Dcr, Drf)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
    crsd <- 100 * SD / mean(Dcr)
    si <- 100 * pmin(Dcr, Drf)/pmax(Dcr, Drf)
    si[is.na(si)] <- 100
    si[si==0] <- 1

    if (SAVE) {
        tmp <- data.frame(
            ID=rownames(kgrid),
            Current=Dcr,
            Reference=Drf)
            #Color_Current=col1[cr],
            #Color_Reference=col1[rf],
            #Color_Difference=col3[df]

        if (TYPE == "S")
            tmp <- tmp[kgrid$useS,]
        if (TYPE == "N")
            tmp <- tmp[kgrid$useN,]
        dbWriteTable(con, "test_num", tmp,
            overwrite=TRUE, row.names=FALSE)
    }

    if (PLOT) {
        Rcr <- make_raster(cr, kgrid, rt)
        Rrf <- make_raster(rf, kgrid, rt)
        Rdf <- make_raster(df-100, kgrid, rt)
        Rsd <- make_raster(crsd, kgrid, rt)
        Rsi <- make_raster(si, kgrid, rt)
        if (TYPE == "S")
            Msk <- Rmasks
        if (TYPE == "N")
            Msk <- Rmaskn
        if (TYPE != "C") {
            Rcr <- mask(Rcr, Msk)
            Rrf <- mask(Rrf, Msk)
            Rdf <- mask(Rdf, Msk)
            Rsd <- mask(Rsd, Msk)
            Rsi <- mask(Rsi, Msk)
        }
        ## add here mask for Rockies if needed

        xy1 <- SpatialPoints(as.matrix(ddd[yyy[,spp] > 0,c("X","Y")]))
        proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        xy1 <- spTransform(xy1, proj4string(rt))
        sam1 <- rasterize(xy1, rt10, field=1, fun='last')
        #xyall@data[[spp]] <- ifelse(is.na(values(sam1)), 0, 1)

        png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/figs/maps2/", spp, ".png"),
            height=1500*2, width=1000*3, res=300)
        op <- par(mfrow=c(2,3), mar=c(2,1,2,3))
        plot(rt, col=CE, axes=FALSE, box=FALSE, main="Reference", legend=FALSE)
        plot(Rrf, add=TRUE, col=col1[1:max(rf)])
        plot(Rw, add=TRUE, col=CW, legend=FALSE)
        plot(rt, col=CE, axes=FALSE, box=FALSE, main="Current", legend=FALSE)
        plot(Rcr, add=TRUE, col=col1[1:max(cr)])
        plot(Rw, add=TRUE, col=CW, legend=FALSE)
        plot(rt, col=CE, axes=FALSE, box=FALSE, main="Difference", legend=FALSE)
        plot(Rdf, add=TRUE, col=col3[min(df):max(df)])
        plot(Rw, add=TRUE, col=CW, legend=FALSE)
        plot(rt, col=CE, axes=FALSE, box=FALSE, main="Intactness", legend=FALSE)
        plot(Rsi, add=TRUE, col=col2[min(si):max(si)])
        plot(Rw, add=TRUE, col=CW, legend=FALSE)
        plot(rnr,col=cnr, axes=FALSE, box=FALSE, main="Detections", legend=FALSE)
        plot(sam0,add=TRUE, col="#ffffffaa", legend=FALSE)
        plot(sam1,add=TRUE, col="red4", legend=FALSE)
        plot(rt, col=CE, axes=FALSE, box=FALSE, main="Current SD", legend=FALSE)
        plot(Rsd, add=TRUE, col=rev(col2)[1:max(cr)])
        plot(Rw, add=TRUE, col=CW, legend=FALSE)
        par(op)
        dev.off()
    }
}

if (SAVE)
    dbDisconnect(con)

## compare

library(cure4insect)
set_options(path = "d:/abmi/reports")
load_common_data()

info <- droplevels(get_species_table("birds"))

for (spp in rownames(tax)) {
    if (tax[spp, "SpeciesID"] %in% info$SpeciesID) {
        cat(spp, "\n");flush.console()
        species <- as.character(tax[spp, "SpeciesID"])
        y <- load_species_data(species)
        TYPE <- "C" # combo
        if (info[species, "model_south"] && !info[species, "model_north"])
            TYPE <- "S"
        if (!info[species, "model_south"] && info[species, "model_north"])
            TYPE <- "N"

        Dcr <- rowSums(y$SA.Curr[match(rownames(kgrid), rownames(y$SA.Curr)),])
        Dcr[is.na(Dcr)] <- 0
        q <- quantile(Dcr, 0.99)
        Dcr[Dcr > q] <- q
        summary(Dcr)
        Drf <- rowSums(y$SA.Ref[match(rownames(kgrid), rownames(y$SA.Ref)),])
        Drf[is.na(Drf)] <- 0
        q <- quantile(Drf, 0.99)
        Drf[Drf > q] <- q
        summary(Drf)
        MAX <- max(Dcr, Drf)

        df <- (Dcr-Drf) / MAX
        df <- sign(df) * abs(df)^0.5
        df <- pmin(200, ceiling(99 * df)+100)
        df[df==0] <- 1
        cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
        rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
        #si <- 100 * pmin(Dcr, Drf)/pmax(Dcr, Drf)
        #si[is.na(si)] <- 100
        #si[si==0] <- 1

        Rcr <- make_raster(cr, kgrid, rt)
        Rrf <- make_raster(rf, kgrid, rt)
        Rdf <- make_raster(df-100, kgrid, rt)
        #Rsi <- make_raster(si, kgrid, rt)
        if (TYPE == "S")
            Msk <- Rmasks
        if (TYPE == "N")
            Msk <- Rmaskn
        if (TYPE != "C") {
            Rcr <- mask(Rcr, Msk)
            Rrf <- mask(Rrf, Msk)
            Rdf <- mask(Rdf, Msk)
            #Rsi <- mask(Rsi, Msk)
        }
        ## add here mask for Rockies if needed

        png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/figs/maps/00-compare-", spp, ".png"),
            height=1500*2, width=1000*3, res=300)
        op <- par(mfrow=c(2,3), mar=c(2,1,2,3))
        plot(Rcr, col=col1[1:max(cr)], axes=FALSE, box=FALSE, main="Current old", legend=FALSE)
        plot(Rrf, col=col1[1:max(rf)], axes=FALSE, box=FALSE, main="Reference old", legend=FALSE)
        plot(Rdf, col=col3[min(df):max(df)], axes=FALSE, box=FALSE, main="Difference old", legend=FALSE)




        load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-04-01/", spp, ".RData"))
        TYPE <- "C" # combo
        if (tax[spp, "ModelSouth"] && !tax[spp, "ModelNorth"])
            TYPE <- "S"
        if (!tax[spp, "ModelSouth"] && tax[spp, "ModelNorth"])
            TYPE <- "N"

        Dcr <- rowSums(Curr)
        q <- quantile(Dcr, 0.99)
        Dcr[Dcr > q] <- q
        summary(Dcr)
        Drf <- rowSums(Ref)
        q <- quantile(Drf, 0.99)
        Drf[Drf > q] <- q
        summary(Drf)
        MAX <- max(Dcr, Drf)

        df <- (Dcr-Drf) / MAX
        df <- sign(df) * abs(df)^0.5
        df <- pmin(200, ceiling(99 * df)+100)
        df[df==0] <- 1
        cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
        rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
        #si <- 100 * pmin(Dcr, Drf)/pmax(Dcr, Drf)
        #si[is.na(si)] <- 100
        #si[si==0] <- 1

        Rcr <- make_raster(cr, kgrid, rt)
        Rrf <- make_raster(rf, kgrid, rt)
        Rdf <- make_raster(df-100, kgrid, rt)
        #Rsi <- make_raster(si, kgrid, rt)
        if (TYPE == "S")
            Msk <- Rmasks
        if (TYPE == "N")
            Msk <- Rmaskn
        if (TYPE != "C") {
            Rcr <- mask(Rcr, Msk)
            Rrf <- mask(Rrf, Msk)
            Rdf <- mask(Rdf, Msk)
            #Rsi <- mask(Rsi, Msk)
        }
        ## add here mask for Rockies if needed

        plot(Rcr, col=col1[1:max(cr)], axes=FALSE, box=FALSE, main="Current new", legend=FALSE)
        plot(Rrf, col=col1[1:max(rf)], axes=FALSE, box=FALSE, main="Reference new", legend=FALSE)
        plot(Rdf, col=col3[min(df):max(df)], axes=FALSE, box=FALSE, main="Difference new", legend=FALSE)

        par(op)
        dev.off()


    }
}

## sector effects

library(mefa4)
library(raster)

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth | Lookup$ModelSouth,])
rownames(tax) <- tax$Code


## species tables for data portal

library(cure4insect)
set_options(path="d:/abmi/reports")
load_common_data()
SP <- get_species_table()
SPP <- rownames(SP)
KT <- get_id_table()

for (spp in SPP) {

    cat(spp, which(SPP==spp), "/", length(SPP), "\n")
    flush.console()

    tx <- as.character(SP[spp, "taxon"])

    TYPE <- "C" # combo
    if (SP[spp, "model_south"] && !SP[spp, "model_north"])
        TYPE <- "S"
    if (!SP[spp, "model_south"] && SP[spp, "model_north"])
        TYPE <- "N"

    y <- load_species_data(spp, boot=FALSE)
    Dcr <- rowSums(y$SA.Curr[match(rownames(KT), rownames(y$SA.Curr)),])
    q <- quantile(Dcr, 0.99, na.rm=TRUE)
    Dcr[!is.na(Dcr) & Dcr > q] <- q
    Drf <- rowSums(y$SA.Ref[match(rownames(KT), rownames(y$SA.Ref)),])
    q <- quantile(Drf, 0.99, na.rm=TRUE)
    Drf[!is.na(Drf) & Drf > q] <- q
    MAX <- max(Dcr, Drf, na.rm=TRUE)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
    d <- data.frame(Row_Col=rownames(KT), Curr=cr, Ref=rf, Diff=df)
    i <- !is.na(cr) & !is.na(rf)
    d <- d[i,]
    write.csv(d, row.names = FALSE,
        file=paste0("s:/BDQT/species-tables-1km/", tx, "/", spp, ".csv"))
}
