## common stuff -----------------------------------------------------------

PROJ <- "north"
spp <- "CAWA"
#STAGE <- "Space"
STAGE <- "HF"

## partial backfilled
## keep all the HF
BF_THIS <- character(0) # current
SEC <- "All" # what HF was not backfilled

## keep energy only
#BF_THIS <- c(
#    "EnSoftLin", "Seismic",
#    "Mine", "Well",
#    "Industrial") # backfill these
#SEC <- "Energy"

#BF_THIS <- c("Crop", "RoughP", "TameP",
#    "HardLin", "TrSoftLin",
#    "EnSoftLin", "Seismic",
#    "Mine", "Well",
#    "Industrial", "Rural", "Urban",
#    "ForHarv") # backfill these
#SEC <- "None"

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
pveg_mine <- 0.44 # proportion of vegetated mines in OSR from NDVI

trVegSS <- trVeg[ss,]
AVegSS <- colSums(trVegSS)

resn <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))

## north estimates
names(en$mods)
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

        ## Mines: 70 is vegetated
        ADnCr[,ch2veg$cr == "MineSite"] <- pveg_mine * ADnCr[,ch2veg$cr == "MineSite"]

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


## explore joint ssh landscape effects

espall <- new.env()
ehfall <- new.env()
#ehfen <- new.env()
load("d:/abmi/AB_data_v2018/data/analysis/birds/pred/oven/OVEN-space-all.RData", envir=espall)
load("d:/abmi/AB_data_v2018/data/analysis/birds/pred/oven/OVEN-hf-all.RData", envir=ehfall)
#load("d:/abmi/AB_data_v2018/data/analysis/birds/pred/oven/OVEN-hf-energy.RData", envir=ehfen)


ssh_labs <- list()
for (i in seq_len(b)) {
    ssh_labs[[i]] <- resn[[as.integer(i)]]$ssh$labels
}
lab <- names(table(unlist(ssh_labs))[table(unlist(ssh_labs)) > 0.5*b])

EN <- c("EnSoftLin", "Seismic", "Mine", "Well", "Industrial")
jd <- data.frame(XSSH)
colnames(jd) <- gsub("_KM", "", colnames(jd))

jd$D1 <- rowMeans(espall$CR) / 100 # per ha density
jd$D2 <- rowMeans(ehfall$CR) / 100 # per ha density
jd$dD <- jd$D2 - jd$D1
#jd$logD <- log(jd$D+0.5)
summary(jd)
sub <- sample(nrow(jd), 5000)
kk <- row_std(groupSums(trVeg[ss,], 2, ch2veg$cr3))
jd <- data.frame(jd, as.matrix(kk[,EN]))
jd$SSH <- rowSums(kk[,colnames(kk) %in% lab,drop=FALSE])
jd$Energy <- as.numeric(rowSums(kk[,EN]))
jd$Soft <- as.numeric(rowSums(kk[,c("EnSoftLin", "Seismic")]))
jd$Ui <- as.numeric(rowSums(kk[,c("Mine", "Well", "Industrial")]))


m1 <- mgcv::gam(D1 ~ s(SSH, THF), jd, family=gaussian)
m2 <- mgcv::gam(D2 ~ s(SSH, THF), jd, family=gaussian)
#m3 <- mgcv::gam(dD ~ s(SSH, THF), jd, family=gaussian)

op <- par(mfrow=c(1,2))
plot(m1, scheme=2, main="Space")
plot(m2, scheme=2, main="HF")
#plot(m3, scheme=2, main="Diff")
par(op)

SSH <- seq(0, 1, 0.01)
THF <- seq(0, 1, 0.01)
pr <- expand.grid(SSH=SSH, THF=THF)
pr$D1 <- predict(m1, pr)
pr$D2 <- predict(m2, pr)

img1 <- list(x=SSH, y=THF, z=matrix(pr$D1, length(SSH), length(THF)))
img2 <- list(x=SSH, y=THF, z=matrix(pr$D2, length(SSH), length(THF)))
MAX <- c(max(img1$z), max(img2$z))

levs <- seq(0.2, 1.2, 0.2)
col <- hcl.colors(100, "YlOrRd", rev = TRUE)
op <- par(mfrow=c(1,2))
image(img1, xlab="SSH", ylab="THF", main="Space", col=col[seq_len(ceiling(100*MAX[1]/max(MAX)))])
contour(img1, add=TRUE, levels=levs)
box()
image(img2, xlab="SSH", ylab="THF", main="HF", col=col[seq_len(ceiling(100*MAX[2]/max(MAX)))])
contour(img2, add=TRUE, levels=levs)
box()
par(op)


val <- jd$SSH
val <- jd$THF
val <- jd$Seismic
val <- jd$EnSoftLin
val <- jd$Ui
val <- jd$Soft

m31 <- mgcv::gam(D1 ~ s(val, k=4), jd, family=gaussian)
m32 <- mgcv::gam(D2 ~ s(val, k=4), jd, family=gaussian)

pr2 <- data.frame(val=seq(0, quantile(val, 0.999), length.out=200))

pp <- qnorm(0.975)
tmp <- predict(m31, pr2, se.fit=TRUE)
pr2$D1 <- tmp$fit
pr2$D1cl <- tmp$se.fit*pp
tmp <- predict(m32, pr2, se.fit=TRUE)
pr2$D2 <- tmp$fit
pr2$D2cl <- tmp$se.fit*pp

plot(D1 ~ val, pr2, type="n", col=4, lwd=2, ylim=c(0, max(pr2$D1, pr2$D2)))
polygon(c(pr2$val, rev(pr2$val)), c(pr2$D1+pr2$D1cl, rev(pr2$D1-pr2$D1cl)), border=NA, col="#0000ff44")
polygon(c(pr2$val, rev(pr2$val)), c(pr2$D2+pr2$D2cl, rev(pr2$D2-pr2$D2cl)), border=NA, col="#ff000044")
lines(D1 ~ val, pr2, col=4, lwd=2)
lines(D2 ~ val, pr2, col=2, lwd=2)


## take


## sector effects

spp <- "OVEN"
STAGE <- "Space"
SEC <- "All"
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
#' This most known and used vegetation index is a simple, but effective vegetation index for
#' quantifying green vegetation. It normalizes green leaf scattering in the
#' Near Infra-red wavelength and chlorophyll absorption in the red wavelength.
#'
#' Values description: The value range of an NDVI is -1 to 1.
#' Negative values of NDVI (values approaching -1) correspond to water.
#' Values close to zero (-0.1 to 0.1) generally correspond to barren areas of rock, sand, or snow.
#' Low, positive values represent shrub and grassland (approximately 0.2 to 0.4),
#' while high values indicate temperate and tropical rainforests (values approaching 1).

x <- read.csv("d:/abmi/AB_data_v2018/data/raw/HFI_Mines_NDVI.csv")
summary(x)
hist(x$NDVImean)
table(x$Mines_FEATURE_TY, is.na(x$NDVImean), useNA="a")
x <- x[!is.na(x$NDVImean),]

m <- lm(NDVImean ~ Mines_FEATURE_TY-1, x)
summary(m)

boxplot(NDVImean ~ Mines_FEATURE_TY-1, x)
abline(h=c(-0.1,0, 0.1), col=2)

a <- aggregate(x$NDVImean, list(FTY=x$Mines_FEATURE_TY), quantile, c(0.5, 0, 1))
a <- a[order(a$x[,"50%"]),]
a

x$FEATURE_TY <- as.character(x$Mines_FEATURE_TY)
x$FEATURE_TY <- factor(x$FEATURE_TY, as.character(a$FTY))
#levels(x$FEATURE_TY) <- paste(LETTERS[1:nlevels(x$FEATURE_TY)], x$FEATURE_TY)

library(ggplot2)

p <- ggplot(x, aes(y=NDVImean, x=FEATURE_TY, fill=FEATURE_TY, color=FEATURE_TY)) +
  geom_violin() + coord_flip() + theme(legend.position = "none") +
  geom_hline(yintercept = 0) + geom_hline(yintercept = -0.1, lty=2) + geom_hline(yintercept = 0.1, lty=2)
p

## calculating average mine NDVI proportions (all)
x$vegetated <- x$NDVImean > 0.1
sum(x$area_ha[x$vegetated]) / sum(x$area_ha) # 0.7

## calculating average mine NDVI proportions (OSR)
library(rgdal)
library(sp)
pl <- readOGR("d:/spatial/Oilsands-Boundaries.gdb", "OilsandRegionDissolve10TM")
xy <- x[,c("x", "y")]
coordinates(xy) <- ~ x + y
proj4string(xy) <- proj4string(pl)
xy <- spTransform(xy, proj4string(pl))
o <- over(xy, pl)
plot(xy, pch=".", col=ifelse(is.na(o$FIELDCODE), 1, 4))
plot(pl, add=T)
x$inOSR <- !is.na(o$FIELDCODE)
sum(x$area_ha[x$vegetated & x$inOSR]) / sum(x$area_ha[x$inOSR]) # 0.44



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
