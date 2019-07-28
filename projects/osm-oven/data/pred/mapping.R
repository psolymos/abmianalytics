## common stuff -----------------------------------------------------------

library(mefa4)
library(intrval)
library(raster)
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

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth,])
rownames(tax) <- tax$Code

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

## for predictions -----------------------------------------------------------------

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


## north models with bootstrap -------------------------------------

PROJ <- "north"
spp <- "OVEN"

## define OSR
library(rgdal)
library(sp)

ogrListLayers("d:/spatial/Oilsands-Boundaries.gdb")
pl <- readOGR("d:/spatial/Oilsands-Boundaries.gdb", "OilsandRegionDissolve10TM")
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(pl))
o <- over(xy, pl)
#plot(xy, pch=".", col=ifelse(is.na(o$FIELDCODE), 1, 4))

ss <- !is.na(o$FIELDCODE)

trVegSS <- trVeg[ss,]
AVegSS <- colSums(trVegSS)

resn <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))

## north estimates
names(en$mods)
#STAGE <- "Space"
STAGE <- "HF"
ESTN <- suppressWarnings(get_coef(resn, Xn, stage=STAGE, na.out=FALSE))

b <- nrow(ESTN)
CR <- matrix(0, sum(ss), b)
rownames(CR) <- rownames(kgrid)[ss]
RF <- CR
HABCR <- matrix(0, nrow(ch2veg), b)
rownames(HABCR) <- rownames(ch2veg)
HABRF <- HABCR

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

## partial backfilled
## keep all the HF
#BF_THIS <- character(0) # current
#SEC <- "All" # what HF was not backfilled

## keep energy only
BF_THIS <- c(
    "EnSoftLin", "Seismic",
    "Mine", "Well",
    "Industrial") # backfill these
SEC <- "Energy"

#BF_THIS <- c("Crop", "RoughP", "TameP",
#    "HardLin", "TrSoftLin",
#    "EnSoftLin", "Seismic",
#    "Mine", "Well",
#    "Industrial", "Rural", "Urban",
#    "ForHarv") # backfill these
#SEC <- "None"

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

## modifying effects: use with Space stage
if (FALSE) {
lamMod <- matrix(NA, b, length(cfn$modif))
colnames(lamMod) <- cfn$modif
for (i in seq_len(b)) {
        estn <- ESTN[i,]
        if (estn["mSoft"] != 0 & estn["mEnSft"] == 0) {
            estn["mEnSft"] <- estn["mSoft"]
            estn["mTrSft"] <- estn["mSoft"]
            estn["mSeism"] <- estn["mSoft"]
        }
        munMod <- structure(numeric(length(cfn$modif)), names=cfn$modif)
        for (k in cfn$modif)
            lamMod[i,k] <- linexp(1, estn[k], pm[k])
}
save(lamMod, ESTN, file="~/repos/abmianalytics/projects/osm-oven/data/pred/lamMod.RData")
}


for (i in seq_len(b)) {
        if (i %% round(b/10) == 0) {
            cat(".")
            flush.console()
        }
        estn <- ESTN[i,]

        ## surrounding SSH and HF
        ssh <- resn[[as.integer(i)]]$ssh
        Xssh0[,"SSH_KM"] <- rowSums(SSH0[,colnames(SSH0) %in% ssh$labels]) # reference
        Xssh0[,"SSH05_KM"] <- sqrt(Xssh0[,"SSH_KM"])
        Xssh[,"SSH_KM"] <- rowSums(SSH[,colnames(SSH0) %in% ssh$labels])
        Xssh[,"SSH05_KM"] <- sqrt(Xssh[,"SSH_KM"])
        if (i == 1) # save 1st run
            XSSH <- Xssh[,c("SSH_KM", "THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM",
                "Noncult_KM", "Cult_KM")]
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

        ## quantiles not applied -- look at that post hoc before summing up
        CR[,i] <- rowSums(ADnCr) # no pair adjustment applied, just ha to km
        RF[,i] <- rowSums(ADnRf) # no pair adjustment applied, just ha to km
        HABCR[colnames(ADnCr),i] <- colSums(ADnCr)
        HABRF[colnames(ADnRf),i] <- colSums(ADnRf)
}
cat("\n")
save(CR, RF, HABRF, HABCR, ch2veg, tv, AVegSS, XSSH,
    file=paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/oven/",
    spp, "-", tolower(STAGE), "-", tolower(SEC), ".RData"))


## explore joint ssh landscape

ssh_labs <- list()
for (i in seq_len(b)) {
    ssh_labs[[i]] <- resn[[as.integer(i)]]$ssh$labels
}
lab <- names(table(unlist(ssh_labs))[table(unlist(ssh_labs)) > 0.5*b])

#d:\abmi\AB_data_v2018\data\analysis\birds\pred\oven\OVEN-hf-all.RData
#d:\abmi\AB_data_v2018\data\analysis\birds\pred\oven\OVEN-hf-energy.RData
#d:\abmi\AB_data_v2018\data\analysis\birds\pred\oven\OVEN-space-all.RData
#d:\abmi\AB_data_v2018\data\analysis\birds\pred\oven\OVEN-space-energy.RData

EN <- c("EnSoftLin", "Seismic", "Mine", "Well", "Industrial")
jd <- data.frame(XSSH)
colnames(jd) <- gsub("_KM", "", colnames(jd))
jd$D <- rowMeans(CR) / 100 # per ha density
jd$logD <- log(jd$D+0.5)
summary(jd)
sub <- sample(nrow(jd), 5000)
kk <- row_std(groupSums(trVeg[ss,], 2, ch2veg$cr3))
jd <- data.frame(jd, as.matrix(kk[,EN]))
jd$SSH <- rowSums(kk[,colnames(kk) %in% lab,drop=FALSE])
jd$Energy <- as.numeric(rowSums(kk[,EN]))
jd$Soft <- as.numeric(rowSums(kk[,c("EnSoftLin", "Seismic")]))
jd$Ui <- as.numeric(rowSums(kk[,c("Mine", "Well", "Industrial")]))

plot(D ~ SSH, jd[sub,], pch=19, col="#00000022")
plot(D ~ THF, jd[sub,], pch=19, col="#00000022")

# bivariate using GAM
m <- mgcv::gam(D ~ s(SSH, THF), jd, family=gaussian)
plot(m, scheme=2, main=spp)

m2 <- mgcv::gam(D ~ s(Lin), jd, family=gaussian)
plot(m2)
m2 <- mgcv::gam(D ~ s(Alien), jd, family=gaussian)
plot(m2)
m2 <- mgcv::gam(D ~ s(Succ), jd, family=gaussian)
plot(m2)
m2 <- mgcv::gam(D ~ s(Cult), jd, family=gaussian)
plot(m2)

m2 <- mgcv::gam(D ~ s(Energy, k=5), jd, family=gaussian)
plot(m2, rug=TRUE)

m2 <- mgcv::gam(D ~ s(Seismic), jd, family=gaussian)
plot(m2, rug=TRUE)

m2 <- mgcv::gam(D ~ s(EnSoftLin), jd, family=gaussian)
plot(m2, rug=TRUE)

m2 <- mgcv::gam(D ~ s(Soft), jd, family=gaussian)
plot(m2, rug=TRUE)

m2 <- mgcv::gam(D ~ s(Ui), jd, family=gaussian)
plot(m2, rug=TRUE)

m <- mgcv::gam(D ~ s(Ui, Soft), jd, family=gaussian)
plot(m, scheme=2, main=spp)

## take


## sector effects

file <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/oven/",
    spp, "-", tolower(STAGE), "-", tolower(SEC), ".RData")

load(file)


levs <- c(
    "MineSite",
    "Pipeline",
    "TransmissionLine",
    "RailHardSurface",
    "RailVegetatedVerge",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge",
    "SeismicLineNarrow",
    "SeismicLineWide",
    "IndustrialSiteRural",
    "UrbanIndustrial",
    "WellSite")

sect <- as.character(ch2veg$sector)
sect[ch2veg$cr %in% levs] <- as.character(ch2veg$cr)[ch2veg$cr %in% levs]
sect[sect %in% c("SeismicLineNarrow", "SeismicLineWide")] <- "SeismicLine"
sect[sect %in% c("Pipeline", "TransmissionLine")] <- "PipeTransLine"
sect[sect %in% c("RailHardSurface",
    "RailVegetatedVerge",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge")] <- "RoadRailVerge"
sect[sect %in% c("IndustrialSiteRural", "UrbanIndustrial")] <- "Industrial"
sect[sect == "Energy"] <- "OtherEnergy"
data.frame(x=table(sect))

HCR <- groupSums(HABCR, 1, sect)
HRF <- groupSums(HABRF, 1, sect)
A <- groupSums(matrix(AVegSS, ncol=1), 1, sect)

save(HCR, HRF, A, file=paste0(
    "~/repos/abmianalytics/projects/osm-oven/data/pred/osm-oven-",
    tolower(STAGE), "-", tolower(SEC), ".RData"))


## mine site: take NDVI and use weighted average

## NDVI was from 2017-2018.
## Each pixel was the median value from June 1 - Aug 15th.
## So it's safe to say the NDVI is from June/July either 2017 or 2018.
## from Sentinel-2 10 m data, top of atmosphere corrected.
#' Sentinel - NDVI (Normalized Difference Vegetation Index)
#'
#' This most known and used vegetation index is a simple, but effective VI for
#' quantifying green vegetation. It normalizes green leaf scattering in the
#' Near Infra-red wavelength and chlorophyll absorption in the red wavelength.
#'
#' Values description: The value range of an NDVI is -1 to 1.
#' Negative values of NDVI (values approaching -1) correspond to water.
#' Values close to zero (-0.1 to 0.1) generally correspond to barren areas of rock, sand, or snow.
#' Low, positive values represent shrub and grassland (approximately 0.2 to 0.4),
#' while high values indicate temperate and tropical rainforests (values approaching 1).




#' This predicts with or without surrounding effects,
#' output is `mu` matrix (PKEY x B) that is on log scale
#' and has no offsets added to it
predict_with_SSH <- function(res, X, SSH=NULL, stage=NULL) {
    est <- suppressWarnings(get_coef(res, X, stage=stage, na.out=FALSE))
#    OK <- !sapply(res, inherits, "try-error")
#    ii <- sapply(res[OK], "[[", "iteration")
#    notNA <- which(OK)
#    est <- est[notNA,,drop=FALSE]
    c1 <- colSums(abs(est)) > 0
    if (any(c1[c("SSH_KM", "SSH05_KM")])) {
        if (is.null(SSH))
            stop("provide SSH")
        essh <- est[,c("SSH_KM", "SSH05_KM")]
        c1[c("SSH_KM", "SSH05_KM")] <- FALSE # drop SSH
        mu <- X[,c1,drop=FALSE] %*% t(est[,c1,drop=FALSE])
        mussh <- mu
        mussh[] <- 0 # put SSH effects here
        for (i in rownames(est)) {
            ssh <- res[[as.integer(i)]]$ssh
            v <- rowSums(SSH[,ssh$labels])
            mussh[,i] <- essh[1,"SSH_KM"]*v + essh[1,"SSH05_KM"]*sqrt(v)
        }
        mu <- mu + mussh # add them up
    } else {
        mu <- X[,c1,drop=FALSE] %*% t(est[,c1,drop=FALSE])
    }
    mu
}
