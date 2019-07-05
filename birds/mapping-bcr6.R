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

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2veg$rf2 <- tv$UseInAnalysisFineAge[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$UseInAnalysisFineAge[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

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


## north models with bootstrap to get pop sizes in BCR6 -------------------------------------

#PROJ <- "north" # all years
PROJ <- "north2" # years matching PIF/BBS years 2006-2015

#DONE <- gsub(".RData", "", list.files("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/"))
#SPP <- setdiff(rownames(tax)[tax$ModelNorth & rownames(tax) %in% colnames(en$OFF)], DONE)

SPP <- gsub(".RData", "", list.files(file.path(ROOT, "out", PROJ)))

ss <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"

for (spp in SPP) {
    resn <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
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
    #fout <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/", spp, ".RData")
    fout <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6-subset/", spp, ".RData")
    save(CR, HAB, file=fout)
    cat(" OK\n")

}


## mapping -------------------------------------------------------------------------



Rmaskn <- make_raster(as.integer(1-kgrid$useS), kgrid, rt)
values(Rmaskn)[values(Rmaskn) == 0] <- NA
Rmasks <- make_raster(as.integer(1-kgrid$useN), kgrid, rt)
values(Rmasks)[values(Rmasks) == 0] <- NA
#Rmaskm <- make_raster(as.integer(kgrid$NRNAME == "Rocky Mountain"), kgrid, rt)
#values(Rmaskm)[values(Rmaskm) == 0] <- NA

## out of BCR6 mask
Rout <- make_raster(as.integer(kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"), kgrid, rt)
values(Rout)[values(Rout) == 0] <- NA
## Water
Rw <- make_raster(as.integer(kgrid$pWater > 0.99), kgrid, rt)
values(Rw)[values(Rw) == 0] <- NA

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
col3 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221"))(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude

## detections
ddd <- nonDuplicated(en$DAT, en$DAT$SS, TRUE)
yyy <- groupSums(en$YY, 1, en$DAT$SS)[rownames(ddd),]
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

rnr <- make_raster(as.integer(kgrid$BCRCODE), kgrid, rt)
cnr <- c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae')
cnr <- cnr[c(5,6,1,2,4,3)]
cnr <- cnr[1:5]

#SPP <- rownames(tax)[tax$ModelNorth & rownames(tax) %in% colnames(en$OFF)]
SPP <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATTW", "BAOR", "BARS",
    "BAWW", "BBMA", "BBWA", "BBWO", "BCCH", "BHCO", "BHVI", "BLJA",
    "BLPW", "BOCH", "BRBL", "BRCR", "BTNW", "CAWA", "CCSP", "CEDW",
    "CHSP", "CMWA", "CONW", "CORA", "COYE", "CSWA", "DEJU", "DOWO",
    "DUFL", "EAKI", "EAPH", "EUST", "EVGR", "FOSP", "GCKI", "GRAJ",
    "GRCA", "HAWO", "HETH", "HOLA", "HOSP", "HOWR", "LCSP", "LEFL",
    "LISP", "MAWA", "MAWR", "MOBL", "MODO", "MOWA", "NOFL", "NOWA",
    "OCWA", "OSFL", "OVEN", "PAWA", "PHVI", "PIGR", "PISI", "PIWO",
    "PUFI", "RBGR", "RBNU", "RCKI", "RECR", "REVI", "ROPI", "RUBL",
    "RUGR", "RWBL", "SAVS", "SOSP", "SPGR", "SWSP", "SWTH", "TEWA",
    "TOWA", "TRES", "VATH", "VEER", "VESP", "WAVI", "WBNU", "WCSP",
    "WETA", "WEWP", "WIWA", "WIWR", "WTSP", "WWCR", "YBFL", "YBSA",
    "YEWA", "YRWA")

SPP <- gsub(".RData", "", list.files(file.path(ROOT, "out", "north2")))

for (spp in SPP) {

    gc()
    cat(spp, "\n");flush.console()

    png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/mapsSubset/", spp, ".png"),
        height=1600*2, width=1000*2, res=300)
    op <- par(mfrow=c(2,2), mar=c(6,0.5,5,0.5))

    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/", spp, ".RData"))

    for (i in 1:ncol(CR)) {
        q <- quantile(CR[,i], 0.99)
        CR[CR[,i] > q,i] <- q
    }
    tab <- t(apply(CR, 1, quantile, c(0.5, 0.025, 0.975)))
    tab <- tab / 100

    cr <- tab[match(rownames(kgrid), rownames(CR)), 1]
    cr[is.na(cr)] <- -1
    Rcr <- make_raster(cr, kgrid, rt)
    values(Rcr)[values(Rcr) < 0] <- NA

    SD <- tab[,3]-tab[,2] # IQR
    SD <- SD[match(rownames(kgrid), rownames(CR))]
    SD[is.na(SD)] <- -1
    Rsd <- make_raster(SD, kgrid, rt)
    values(Rsd)[values(Rsd) < 0] <- NA

    #xy1 <- SpatialPoints(as.matrix(ddd[yyy[,spp] > 0,c("X","Y")]))
    #proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    #xy1 <- spTransform(xy1, proj4string(rt))
    #sam1 <- rasterize(xy1, rt10, field=1, fun='last')


#    plot(rnr,col=cnr, axes=FALSE, box=FALSE, legend=FALSE,
#        main=paste0(as.character(tax[spp, "CommonName"]), "\n(", as.character(tax[spp, "ScientificName"]), ")"))
#    plot(sam0,add=TRUE, col="#ffffffaa", legend=FALSE)
#    plot(sam1,add=TRUE, col="red4", legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, legend=FALSE,
        main=paste0(as.character(tax[spp, "CommonName"]), "\n(", as.character(tax[spp, "ScientificName"]), ")"))
    plot(Rcr, add=TRUE, col=col1, legend=FALSE)
    plot(Rw, add=TRUE, col=CE, legend=FALSE)
    plot(Rcr, col=col1, legend.only=TRUE, horizontal = TRUE,
        legend.args = list(text="Median Density (inds./ha)", side = 1, line = 2, cex=0.5))

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="", legend=FALSE)
    plot(Rsd, add=TRUE, col=rev(col2), legend=FALSE)
    plot(Rw, add=TRUE, col=CE, legend=FALSE)
    plot(Rsd, col=rev(col2), legend.only=TRUE, horizontal = TRUE,
        legend.args = list(text="95% CI Range", side = 1, line = 2, cex=0.5))


    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6-subset/", spp, ".RData"))

    for (i in 1:ncol(CR)) {
        q <- quantile(CR[,i], 0.99)
        CR[CR[,i] > q,i] <- q
    }
    tab <- t(apply(CR, 1, quantile, c(0.5, 0.025, 0.975)))
    tab <- tab / 100

    cr <- tab[match(rownames(kgrid), rownames(CR)), 1]
    cr[is.na(cr)] <- -1
    Rcr <- make_raster(cr, kgrid, rt)
    values(Rcr)[values(Rcr) < 0] <- NA

    SD <- tab[,3]-tab[,2] # IQR
    SD <- SD[match(rownames(kgrid), rownames(CR))]
    SD[is.na(SD)] <- -1
    Rsd <- make_raster(SD, kgrid, rt)
    values(Rsd)[values(Rsd) < 0] <- NA

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Subset", legend=FALSE)
    plot(Rcr, add=TRUE, col=col1, legend=FALSE)
    plot(Rw, add=TRUE, col=CE, legend=FALSE)
    plot(Rcr, col=col1, legend.only=TRUE, horizontal = TRUE,
        legend.args = list(text="Median Density (inds./ha)", side = 1, line = 2, cex=0.5))

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="", legend=FALSE)
    plot(Rsd, add=TRUE, col=rev(col2), legend=FALSE)
    plot(Rw, add=TRUE, col=CE, legend=FALSE)
    plot(Rsd, col=rev(col2), legend.only=TRUE, horizontal = TRUE,
        legend.args = list(text="95% CI Range", side = 1, line = 2, cex=0.5))

    par(op)
    dev.off()

}

# habitat assoc figures

library(cure4insect)
set_options(path = "s:/reports")
load_common_data()

for (spp in SPP) {
    png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/mapsMedian/", spp, "-habitats.png"),
        height=700, width=1400, res=150)
    plot_abundance(as.character(TAB[spp, "SpeciesID"]), "veg_coef")
    dev.off()
}

