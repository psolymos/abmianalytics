#' ---
#' title: "Bird data packaging"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Sep 23, 2020"
#' output: pdf_document
#' ---
#'
#' # Preamble
#'
#' In this script we filter, transform, mutate, and get some offsets calculated.
#'
library(mefa4)
library(intrval)
library(maptools)
source("~/repos/abmianalytics/birds/00-functions.R")
knitr::opts_chunk$set(eval=FALSE)
#load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData")
load("d:/abmi/AB_data_v2020/data/analysis/species/birds/ab-birds-all-2020-09-23.RData")
if (FALSE) {
    ## quick hack to get reference data set for prediction
    vc1[] <- 0
    vc1[,colnames(vr1)] <- vr1
    vc2[] <- 0
    vc2[,colnames(vr2)] <- vr2
    sc1[] <- 0
    sc1[,colnames(sr1)] <- sr1
    sc2[] <- 0
    sc2[,colnames(sr2)] <- sr2
}
#'
#' # Remove missing data
#'
#' Check NAs first
#'
## checking NAs in veg+soil+hf summaries
dd$vshf <- ifelse(rowSums(is.na(vc1)) > 0, NA, 0)
(aa <- data.frame(n_of_NAs=colSums(is.na(dd))))
#'
#' - `ROAD`: NAs are expected because we need hard linear amount for non-BBS surveys
#' - `DATE` and `DATI`: we use use constant singing rate where these are NA
#' - `X` and `Y`: should be dropped
dd <- droplevels(dd[!is.na(dd$X) & !is.na(dd$vshf),])
#' - climate: out-of-AB-bound issue, to be dropped
dd <- dd[!is.na(dd$pAspen) & !is.na(dd$NRNAME),]
#' - `YEAR`: must be in this range
dd <- dd[dd$YEAR %[]% c(1993, 2019),]
#' - date and time: constrain the ranges to match QPAD expectations
dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
dd$start <- as.POSIXlt(dd$DATI)$hour + as.POSIXlt(dd$DATI)$min / 60
keep <- is.na(dd$DATI) # these will be constant phi
keep[dd$JULIAN %[]% c(125, 200)] <- TRUE
keep[dd$start %[]% c(3, 12)] <- TRUE
dd <- droplevels(dd[keep,])
#' Normalized ordinal day
dd$JDAY <- dd$JULIAN / 365
#' Normalized time since local sunrise
Coor <- as.matrix(dd[,c("X", "Y")])
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dd$srise <- NA
dd$srise[subset] <- sr
dd$TSSR <- (dd$start - dd$srise) / 24
#' Quadratic terms
dd$JDAY2 <- dd$JDAY^2
dd$TSSR2 <- dd$TSSR^2
#' Maximum counting distance in 100 m units (area in ha)
dd$MAXDIS <- dd$MAXDIS / 100
table(dd$MAXDIS)
#'
#' # Subset species data
#'
#' Keep species with at least 20 detections
#'
yy <- yy[rownames(dd), colSums(yy > 0) >= 20]
#'
#' # Derive predictor variables
#'
#' Year relative to start year 1993
dd$YR <- dd$YEAR - min(dd$YEAR)
#' Point level intersection (veg+soil+HF)
dd <- data.frame(dd, vs0[rownames(dd),])
#' Read lookup tables
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv0) <- tv0[,1]
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v2020.csv")
rownames(tv) <- tv[,1]
tv0 <- tv0[rownames(tv),]
tv <- data.frame(tv, tv0)
tv$UseInAnalysisNoAge <- gsub("CC", "", gsub("[[:digit:]]", "", as.character(tv$UseInAnalysis)))
tv$UseInAnalysisNoAge[endsWith(tv$UseInAnalysisNoAge, "R")] <-
    substr(tv$UseInAnalysisNoAge[endsWith(tv$UseInAnalysisNoAge, "R")], 1,
        nchar(tv$UseInAnalysisNoAge[endsWith(tv$UseInAnalysisNoAge, "R")])-1)
ts0 <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts0) <- ts0[,1]
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(ts) <- ts[,1]
ts0 <- ts0[rownames(ts),]
ts <- data.frame(ts, ts0)
#' 7 ha (150 m radius) level proportions used downsream
dd$pWater <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_water"]])
dd$pRoad <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_road"]])
dd$pRoadVeg <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_road_veg"]])
dd$pClosed <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_closed"]])
dd$pHarest <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_harvest"]])
dd$pHF <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_HF"]])
#' 4 and 2 level land cover classification for QPAD offsets
lcc4 <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"LCC4"]))
tmp <- find_max(lcc4)
dd$LCC4 <- tmp$index
dd$LCC2 <- dd$LCC4
levels(dd$LCC2) <- c("OpenWet", "Forest", "Forest", "OpenWet")
#' 1 km$^2$ (564 m radius) level proportions
dd$pWater_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_water"]])
dd$pWater2_KM <- dd$pWater_KM^2
dd$pWet_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv[colnames(vc2), "is_wet"]])
dd$pWetWater_KM <-dd$pWater_KM + dd$pWet_KM
#' Placeholder for Surrounding Suitable Habitat (SSH) at the 1 km$^2$ level
dd$SSH_KM <- 0
dd$SSH05_KM <- sqrt(dd$SSH_KM)
#' Surrounding footprint at the 1 km$^2$ level
dd$THF_KM <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc2), "is_HF"]])
dd$Lin_KM <- rowSums(row_std(vc1[rownames(dd),])[,
    tv[colnames(vc2), "is_HF"] & tv[colnames(vc2), "is_linear"]])
## note: no abandoned or rough pasture here
dd$Cult_KM <- rowSums(row_std(vc1[rownames(dd),])[,c("CultivationCrop",
    "CultivationTamePasture", "HighDensityLivestockOperation")])
dd$Alien_KM <- rowSums(row_std(vc1[rownames(dd),])[,
    tv[colnames(vc2), "is_HF"] & tv[colnames(vc2), "is_alien"]])
dd$Nonlin_KM <- dd$THF_KM - dd$Lin_KM
dd$Noncult_KM <- dd$THF_KM - dd$Cult_KM
dd$Succ_KM <- dd$THF_KM - dd$Alien_KM
dd$THF2_KM <- dd$THF_KM^2
dd$Succ2_KM <- dd$Succ_KM^2
dd$Alien2_KM <- dd$Alien_KM^2
dd$Noncult2_KM <- dd$Noncult_KM^2
dd$Nonlin2_KM <- dd$Nonlin_KM^2
#' Climate variable and lat/long transformations
dd <- data.frame(dd, transform_clim(dd))
#' Road is defined based on hard linear threshold where it is missing (non-BBS)
dd$ROAD[is.na(dd$ROAD) & dd$pRoad > 0.04] <- 1
dd$ROAD[is.na(dd$ROAD)] <- 0
table(dd$ROAD, dd$pRoad > 0.04, useNA="a")
#'
#' ## Land cover reclass
#'
vc1r <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"UseInAnalysisNoAge"]))
vr1r <- row_std(groupSums(vr1[rownames(dd),], 2, tv[colnames(vr1),"UseInAnalysisNoAge"]))
vc2r <- row_std(groupSums(vc2[rownames(dd),], 2, tv[colnames(vc2),"UseInAnalysisNoAge"]))
vr2r <- row_std(groupSums(vr2[rownames(dd),], 2, tv[colnames(vr2),"UseInAnalysisNoAge"]))

sc1r <- row_std(groupSums(sc1[rownames(dd),], 2, ts[colnames(sc1),"UseInAnalysis"]))
sr1r <- row_std(groupSums(sr1[rownames(dd),], 2, ts[colnames(sr1),"UseInAnalysis"]))
sc2r <- row_std(groupSums(sc2[rownames(dd),], 2, ts[colnames(sc2),"UseInAnalysis"]))
sr2r <- row_std(groupSums(sr2[rownames(dd),], 2, ts[colnames(sr2),"UseInAnalysis"]))
#' Point intersection based reclassified variables (w/o forest age)
dd$vegpt <- as.factor(tv$UseInAnalysisNoAge[match(dd$VEGHFAGEclass, rownames(tv))])
dd$soilpt <- as.factor(ts$UseInAnalysis[match(dd$SOILHFclass, rownames(ts))])
#' Dominant land cover type: veg+HF.
#' Calculate proportion without classes we do not use, because:
#' - it is not a stratum we are interested in (water, snow/ice is 0 density),
#' - its is too small of a feature to make up a full 7-ha buffer.
tmp <- find_max(vc1r[,colnames(vc1r) %nin% c("Water","HWater", "SnowIce", "Bare",
    "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
#' Keep levels consistent
dd$vegc <- factor(as.character(tmp$index), levels(dd$vegpt))
dd$vegv <- tmp$value
#'
#' Finding harvest areas
cc <- paste0(ifelse(tv[colnames(vc1), "is_harvest"], "CC", ""),
    as.character(tv[colnames(vc1), "UseInAnalysisNoAge"]))
tmp <- row_std(groupSums(vc1[rownames(dd),], 2, cc))
tmp <- find_max(tmp[,colnames(tmp) %ni% c("Water","HWater", "SnowIce", "Bare",
    "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")])
dd$vegccc <- factor(as.character(tmp$index), c(levels(dd$vegc),
    paste0("CC", c(c("Spruce","Decid","Mixedwood","Pine")))))
dd$vegvcc <- tmp$value
table(cc=dd$vegccc, not=dd$vegc)
#levels(dd$vegc) <- levels(dd$vegccc)
#' Indicator variable for harvest area
dd$isCC <- startsWith(as.character(dd$vegccc), "CC")
tmp <- dd$vegccc
levels(tmp) <- gsub("CC", "", levels(dd$vegccc))
dd$vegc[dd$isCC] <- tmp[dd$isCC]
dd$vegw[dd$isCC] <- dd$vegvcc[dd$isCC]
#' Weight is a function of proportion:
#' - 0 below 0.25 (too small to be considered dominant)
#' - 1 above 0.75 (it is large enough to consider it dominant)
#' - 0-1 in between
dd$vegw <- pmax(0, pmin(1, 2*dd$vegv-0.5))
#' Check how often we get the dominant class at the center
a <- table(pt=dd$vegpt,bf=droplevels(dd$vegc))
a <- a[colnames(a),]
sum(diag(a))/sum(a) # pretty good
#' Finally: drop unused levels and relevel to have Decid as reference
dd$vegc <- droplevels(dd$vegc)
dd$vegc <- relevel(dd$vegc, "Decid")
dd$vegccc <- droplevels(dd$vegccc)
dd$vegccc <- relevel(dd$vegccc, "Decid")
data.frame(table(dd$vegc))

#' Indicator variables for stand types (don't use age for TreedSwamp)
dd$isMix <- ifelse(dd$vegc == "Mixedwood", 1L, 0L)
dd$isWSpruce <- ifelse(dd$vegc == "Spruce", 1L, 0L)
dd$isPine <- ifelse(dd$vegc == "Pine", 1L, 0L)
dd$isBog <- ifelse(dd$vegc == "TreedBog", 1L, 0L)
dd$isFen <- ifelse(dd$vegc == "TreedFen", 1L, 0L)
dd$isBogFen <- ifelse(dd$vegc %in% c("TreedBog", "TreedFen"), 1L, 0L)
dd$isUpCon <- ifelse(dd$vegc %in% c("Spruce", "Pine"), 1L, 0L)
dd$isCon <- ifelse(dd$vegc %in% c("TreedBog", "TreedFen",
    "Spruce", "Pine"), 1L, 0L)
#'
#' ## Weighted age calculation
#'
tv <- tv[colnames(vc1),]
ac <- as.character(tv[, "AGE"])
ac[is.na(ac)] <- ""
ac[tv$UseInAnalysisNoAge == "TreedSwamp"] <- ""
vc1age <- row_std(groupSums(vc1[rownames(dd),], 2, ac))
## exclude unknown (0) and non-forest (blank)
AgePtCr <- t(vc1age[,c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")])
AgeMin <- structure(c(0,10,20,40,60,80,100,120,140,160)/200,
    names=c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
dd$wtAge <- colSums(AgePtCr * AgeMin) / colSums(AgePtCr)
dd$wtAge[is.na(dd$wtAge)] <- 0
dd$isFor <- dd$vegc %in% c("Spruce","Decid","Mixedwood","Pine","TreedBog", "TreedFen")
dd$wtAge[!dd$isFor] <- 0
dd$wtAge2 <- dd$wtAge^2
dd$wtAge05 <- sqrt(dd$wtAge)
#'
#' ## Forestry convergence
#'
## fCC1: linear
MAXFOR <- 50/200
dd$fCC1 <- 0
dd$fCC1[dd$isCC==1] <- pmax(0, 1 - (dd$isCC * dd$wtAge/MAXFOR)[dd$isCC==1])
#plot(fCC1 ~ wtAge, dd[dd$isCC==1,])
## fCC2: Dave Huggard's recovery trajectories
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3,
    82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85,
    88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100
#data.frame(Age=age*200, wtAge=age, fCC2_conif=conif, fCC2_decid=decid)
dd$fCC2 <- 0
tmp1 <- approxfun(age, decid)(dd$wtAge)
tmp1[is.na(tmp1)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- dd$isFor & !dd$isCon & dd$isCC==1
dd$fCC2[ii] <- tmp1[ii]
tmp2 <- approxfun(age, conif)(dd$wtAge)
tmp2[is.na(tmp2)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- dd$isCon & dd$isCC==1
dd$fCC2[ii] <- tmp2[ii]
#plot(dd$wtAge, dd$fCC1, col=2, pch=".")
#points(dd$wtAge, dd$fCC2, col=4, pch=".")
sum(is.na(dd$fCC2))
by(dd$wtAge*200, list(veg=interaction(dd$isCC,dd$vegc,drop=TRUE)), fstat, level=1)
#'
#' ## Vegetation and footprint based modifier variables
#'
#' These variables will modify the effects given some adjacent habitat around them.
#' Adjacent habitat is the dominant land cover, e.g. deciduous, or crop around a wellpad.
#'
#' Modifiers used only in the north
dd$mEnSft <- vc1r[,"EnSoftLin"]
dd$mTrSft <- vc1r[,"TrSoftLin"]
dd$mSeism <- vc1r[,"EnSeismic"]
#' Modifiers used in the north and the south
dd$mWell <- vc1r[,"Well"]
dd$mHard <- vc1r[,"HardLin"] # optional, use ROAD instead
dd$mSoft <- dd$mSeism + dd$mEnSft + dd$mTrSft
#'
#' Data subset to be used in the north:
dd$useNorth <- TRUE
#' - do not consider the Grassland natural region
dd$useNorth[dd$NRNAME == "Grassland"] <- FALSE
#' - do not consider sites that are mostly open water
#'   (makes dominant land cover designation questionable)
dd$useNorth[dd$useNorth & (dd$pWater > 0.5 | dd$vegpt == "Water")] <- FALSE
#' - drop observations that were assigned 0 weight (will not contribute to likelihood)
dd$useNorth[dd$vegw == 0] <- FALSE
#'
#' Surrounding habitat and classification
dd$vegca <- dd$vegc
levels(dd$vegca) <- c(levels(dd$vegca), "SpruceO","DecidO","MixedwoodO","PineO","TreedBogO", "TreedFenO")
ii <- dd$vegc %in% c("Spruce","Pine","TreedBog", "TreedFen") & dd$wtAge*200 >= 80
dd$vegca[ii] <- paste0(as.character(dd$vegc[ii]), "O")
ii <- dd$vegc %in% c("Decid","Mixedwood") & dd$wtAge*200 >= 50
dd$vegca[ii] <- paste0(as.character(dd$vegc[ii]), "O")
table(dd$vegca, dd$isCC)

ao <- paste0(as.character(tv[colnames(vc1), "UseInAnalysisNoAge"]),
    ifelse(tv[colnames(vc1), "MatureOld"], "O", ""))
SSH_veg <- row_std(groupSums(vc2[rownames(dd),], 2, ao))
SSH_veg <- SSH_veg[,colnames(SSH_veg) %ni% c("Water","HWater", "SnowIce", "Bare",
    "EnSeismic", "HardLin", "TrSoftLin", "EnSoftLin", "Well")]
compare_sets(levels(dd$vegca), colnames(SSH_veg))
#'
#' ## Soil reclass for the south
#'
#' Dominant land cover type: soil+HF (similarly to the veg counterpart).
tmp <- find_max(sc1r[,colnames(sc1r) %ni% c("SoilWater", "SoilUnknown", "HWater",
    "EnSeismic", "EnSoftLin", "TrSoftLin", "HardLin", "Well")])
dd$soilc <- factor(as.character(tmp$index), levels(dd$soilpt))
dd$soilv <- tmp$value
#' Use backfilled soil type where HFor is dominant
ii <- dd$soilc == "HFor"
tmp2 <- find_max(sr1r[,colnames(sr1r) %ni% c("SoilWater", "SoilUnknown", "HWater")])
dd$soilc[ii] <- tmp2$index[ii]
dd$soilv[ii] <- tmp2$value[ii]
dd$soilw <- pmax(0, pmin(1, 2*dd$soilv-0.5))
a <- table(pt=dd$soilpt,bf=droplevels(dd$soilc))
a <- a[colnames(a),]
sum(diag(a))/sum(a) # pretty good
dd$soilc <- droplevels(dd$soilc)
dd$soilc <- relevel(dd$soilc, "Loamy")
#'
#' Data subset to be used in the south:
dd$useSouth <- FALSE
#' - use the grassland and Parkland natural regions
dd$useSouth[dd$NRNAME %in% c("Grassland", "Parkland")] <- TRUE
# ' - plus the Dry Mixedwood subregion within the Boreal
dd$useSouth[dd$NSRNAME %in% c("Dry Mixedwood")] <- TRUE
#' - but only below the magical latitude limit of 56.7 degrees
dd$useSouth[dd$useSouth & dd$Y > 56.7] <- FALSE
#' - do not consider sites that are mostly open water
#'   (makes dominant land cover designation questionable)
dd$useSouth[dd$useSouth & (dd$pWater > 0.5 | dd$vegpt == "Water")] <- FALSE
#' - do not consider sites where we have no soil info
dd$useSouth[dd$useSouth & sr1r[,"SoilUnknown"] > 0] <- FALSE
#' - drop observations that were assigned 0 weight (will not contribute to likelihood)
dd$useSouth[dd$soilw == 0] <- FALSE
#'
#' Surrounding land cover in the south
SSH_soil <- sc2r[rownames(dd),]
#' ## Simplistic exploration of land cover types using opticut
#'
if (FALSE) {
    library(opticut)
    library(parallel)
    ncl <- 8

    ## south
    ys <- groupSums(yy[dd$useSouth,], 1, dd$SS[dd$useSouth])
    ys[ys > 0] <- 1
    ys <- ys[,colSums(ys) >= 20]
    ds <- nonDuplicated(dd, SS, TRUE)[rownames(ys),]

    table(ds$soilc, ds$WELL)
    hist(ds$soilv)

    cl <- makeCluster(ncl)
    o <- opticut(as.matrix(ys) ~ ROAD, data=ds, strata=ds$soilc,
        dist="binomial:cloglog", weights=ds$soilv, cl=cl)
    stopCluster(cl)
    plot(o, cex.axis=0.5)

    ## north
    yn <- groupSums(yy[dd$useNorth,], 1, dd$SS[dd$useNorth])
    yn[yn > 0] <- 1
    yn <- yn[,colSums(yn) >= 100]
    dn <- nonDuplicated(dd, SS, TRUE)[rownames(yn),]

    data.frame(table(dn$vegc))
    hist(dn$vegv)

    cl <- makeCluster(ncl)
    o <- opticut(as.matrix(yn) ~ ROAD, data=dn, strata=dn$vegc,
        dist="binomial:cloglog", weights=dn$vegw, cl=cl)
    stopCluster(cl)
    plot(o, cex.axis=0.5)
}
#'
#' # QPAD offsets
#'
#' Design matrices for singing rates (`Xp`) and for EDR (`Xq`)
Xp <- cbind("(Intercept)"=1, as.matrix(dd[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=dd$pClosed,
    LCC2OpenWet=ifelse(dd$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(dd$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(dd$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(dd$LCC4=="Wet", 1, 0))
summary(Xp)
summary(Xq)
#' Use version 3 of BAM QPAD estimates
library(QPAD)
load_BAM_QPAD(version=3)
sppp <- intersect(colnames(yy), getBAMspecieslist())
#' Save values in a matrix
off <- matrix(NA, nrow(dd), length(sppp))
rownames(off) <- rownames(dd)
colnames(off) <- sppp
#' Loop for species
for (spp in sppp) {
    ## print out where we are
    cat(spp, "\n");flush.console()
    p <- rep(NA, nrow(dd))
    A <- q <- p
    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model
    mi <- bestmodelBAMspecies(spp, type="BIC",
        model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    ## design matrices matching the coefs
    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0
    ## calculate p, q, and A based on constant phi and tau for the respective NAs
    p[!OKp] <- sra_fun(dd$MAXDUR[!OKp], cf0[1])
    unlim <- ifelse(dd$MAXDIS[!OKq] == Inf, TRUE, FALSE)
    A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dd$MAXDIS[!OKq]^2)
    q[!OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[!OKq], cf0[2]))
    ## calculate time/lcc varying phi and tau for non-NA cases
    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    p[OKp] <- sra_fun(dd$MAXDUR[OKp], phi1)
    unlim <- ifelse(dd$MAXDIS[OKq] == Inf, TRUE, FALSE)
    A[OKq] <- ifelse(unlim, pi * tau1^2, pi * dd$MAXDIS[OKq]^2)
    q[OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[OKq], tau1))
    ## log(0) is not a good thing, apply constant instead
    ii <- which(p == 0)
    p[ii] <- sra_fun(dd$MAXDUR[ii], cf0[1])
    ## store, next
    off[,spp] <- log(p) + log(A) + log(q)
}
## sanity checks
(Ra <- apply(off, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))
off_mean <- log(rowMeans(exp(off)))
#'
#' # Bootstrap blocking units
#'
#' Spatial blocking units are based on unique locations (SS), we aime for a well balanced
#' data set
cn=c("PCODE", "SS", "SSYR", "PKEY", "YEAR", "DATE", "DATI", "MAXDUR",
    "MAXDIS", "CMETHOD", "ROAD", "X", "Y", "NRNAME", "NSRNAME", "LUF_NAME", "useNorth", "useSouth")
ddd=dd[,cn]
#save(ddd, file="~/GoogleWork/abmi/bird-data.RData")
#load("~/GoogleWork/abmi/bird-data.RData")
#with(ddd, plot(X, Y, col=NRNAME, pch="."))
tmp <- nonDuplicated(ddd, SS, TRUE)
cx <- cut(tmp$X, c(-121, -116, -112,-109))
cy <- cut(tmp$Y, c(48, 51, 54, 57, 61))
ct <- cut(tmp$YEAR, c(1992, 2001, 2009, 2013, 2019))
table(cy, cx)
table(ct)
ftable(ct, cy, cx)

dd$BLOCK_X <- cut(dd$X, c(-121, -116, -112,-109))
dd$BLOCK_Y <- cut(dd$Y, c(48, 51, 54, 57, 61))
dd$BLOCK_T <- cut(dd$YEAR, c(1992, 2001, 2009, 2013, 2018))
dd$BLOCK_XY <- interaction(droplevels(dd$BLOCK_X), droplevels(dd$BLOCK_Y), sep="::", drop=TRUE)
dd$BLOCK_XYT <- interaction(dd$BLOCK_XY, dd$BLOCK_T, sep="::", drop=TRUE)
ftable(dd$BLOCK_T, dd$BLOCK_Y, dd$BLOCK_X)
#'
#' Random quantiles: these are also based on SS
set.seed(1)
tmp$RND <- sample.int(100, nrow(tmp), replace=TRUE)
dd$RND <- tmp$RND[match(dd$SS, tmp$SS)]
#'
#' # Model subsets
#'
library(mefa4)
library(opticut)
source("~/repos/abmianalytics/birds/00-functions.R")
NMIN <- 20
B <- 256 # 240
#'
#' ## South
#'
source("~/repos/abmianalytics/birds/models-soil.R")
setdiff(get_terms(mods_soil, "list"), colnames(dd))
rm(DAT, YY, OFF, OFFmean, SSH, BB, mods)

cn2 <- c(cn, get_terms(mods_soil, "list"), "soilw")
DAT <- droplevels(dd[dd$useSouth & dd$RND > 10, cn2])
YY <- yy[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= NMIN]
#YY <- YY[,colSums(groupSums(YY, 1, DAT$SS) > 0) >= NMIN]
OFF <- off[rownames(DAT), intersect(colnames(off), colnames(YY))]
OFFmean <- off_mean[rownames(DAT)]
mods <- mods_soil
#mods$SSH <- NULL
SSH <- SSH_soil[rownames(DAT),]

BB <- pbapply::pbsapply(1:B, bfun, DAT$SS, DAT$BLOCK_XYT)
nrow(DAT)
max(BB)
(lu <- length(unique(as.numeric(BB))))
stopifnot(all(BB <= nrow(DAT)))
stopifnot(lu <= nrow(DAT))
nrow(BB)

z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="soilw", ssh_class="soilc", ssh_fit="Space")
z$timer
cat("Estimate for", ncol(YY), "species and", B, "runs is", ceiling(unname(ncol(YY)*B*z$timer[3])/(60*60)), "hrs\n")

save(DAT, YY, OFF, OFFmean, SSH, BB, mods,
    file="d:/abmi/AB_data_v2020/data/analysis/species/birds/data/ab-birds-south-2020-09-23.RData")
if (FALSE) {
## update S models
library(mefa4)
load("d:/abmi/AB_data_v2020/data/analysis/species/birds/data/ab-birds-south-2020-09-23.RData")

mods$Hab[[3]] <- . ~ . + soilc2
mods$Hab[[4]] <- . ~ . + soilc2 + pAspen
mods$Hab[[5]] <- . ~ . + soilc1
mods$Hab[[6]] <- . ~ . + soilc1 + pAspen

str(DAT)
DAT$soilc2 <- as.character(DAT$soilc)
DAT$soilc1 <- as.character(DAT$soilc)
DAT$soilc1[DAT$soilc1 %in% c("Loamy", "ClaySub", "SandyLoam",
    "Other", "ThinBreak", "RapidDrain", "Blowout")] <- "SoilNative"
DAT$soilc1[DAT$soilc1 %in% c("Industrial", "Mine", "Urban", "Rural")] <- "UrbIndRur"
DAT$soilc2[DAT$soilc2 %in% c("Loamy", "ClaySub", "SandyLoam")] <- "ClaySubLoamSand"
DAT$soilc2[DAT$soilc2 %in% c("Other", "ThinBreak", "RapidDrain", "Blowout")] <- "OtherBlowThinRapid"
DAT$soilc2[DAT$soilc2 %in% c("Industrial", "Mine", "Urban")] <- "UrbInd"
DAT$soilc1 <- as.factor(DAT$soilc1)
DAT$soilc1 <- relevel(DAT$soilc1, "SoilNative")
DAT$soilc2 <- as.factor(DAT$soilc2)
DAT$soilc2 <- relevel(DAT$soilc2, "ClaySubLoamSand")
addmargins(table(DAT$soilc, DAT$soilc2))
addmargins(table(DAT$soilc, DAT$soilc1))

save(DAT, YY, OFF, OFFmean, SSH, BB, mods,
    file="d:/abmi/AB_data_v2020/data/analysis/species/birds/data/ab-birds-south-2020-12-04.RData")

}
#'
#' ## North
#'
source("~/repos/abmianalytics/birds/models-veg.R")
setdiff(get_terms(mods_veg, "list"), colnames(dd))
rm(DAT, YY, OFF, OFFmean, SSH, BB, mods)

cn2 <- c(cn, get_terms(mods_veg, "list"), "vegw", "vegca")
DAT <- dd[dd$useNorth & dd$RND > 10, cn2]
YY <- yy[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= NMIN]
#YY <- YY[,colSums(groupSums(YY, 1, DAT$SS) > 0) >= NMIN]
OFF <- off[rownames(DAT), intersect(colnames(off), colnames(YY))]
OFFmean <- off_mean[rownames(DAT)]
mods <- mods_veg
#mods$SSH <- NULL
SSH <- SSH_veg[rownames(DAT),]

BB <- pbapply::pbsapply(1:B, bfun, DAT$SS, DAT$BLOCK_XYT)
nrow(DAT)
max(BB)
(lu <- length(unique(as.numeric(BB))))
stopifnot(all(BB <= nrow(DAT)))
stopifnot(lu <= nrow(DAT))
nrow(BB)

z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="vegw", ssh_class="vegca", ssh_fit="Space")
z$timer
cat("Estimate for", ncol(YY), "species and", B, "runs is", ceiling(unname(ncol(YY)*B*z$timer[3])/(60*60)), "hrs\n")

save(DAT, YY, OFF, OFFmean, SSH, BB, mods,
    file="d:/abmi/AB_data_v2020/data/analysis/species/birds/data/ab-birds-north-2020-09-23.RData")
#'
#' ## Validation subsets
#'
source("~/repos/abmianalytics/birds/models-veg.R")
setdiff(get_terms(mods_veg, "list"), colnames(dd))
rm(DAT, YY, OFF, OFFmean, SSH, BB, mods)

## find center point not that is in grassland
zz <- droplevels(dd[dd$CMETHOD=="RF" & substr(as.character(dd$SS),1,2) != "OG",])
tmp <- strsplit(as.character(zz$SS), "_")
zz$ABMIsite <- sapply(tmp, "[[", 1)
zz$ABMIbirdpt <- sapply(tmp, "[[", 2)
n <- table(zz$ABMIsite)
n <- n[n==9]
zz$All9 <- zz$ABMIsite %in% names(n)
nn <- sum_by(zz$NRNAME == "Grassland" | zz$Y < 50, zz$ABMIsite)
zz$NotGr <- zz$ABMIsite %in% rownames(nn)[nn[,"x"] == 0]

table(ngr=zz$NotGr, a=zz$All9)/9

nam <- unique(zz[zz$NotGr & zz$All9, "ABMIsite"])
set.seed(1)
nam500 <- sample(nam, 500)

dd$ABMIsite <- zz$ABMIsite[match(rownames(dd), rownames(zz))]
dd$ABMIsite[is.na(dd$ABMIsite)] <- ""
dd$validation <- dd$ABMIsite %in% nam500
dd$validation[dd$PCODE == "BU_BG"] <- TRUE
table(validation=dd$validation,north=dd$useNorth)

with(dd[!dd$validation,], plot(X, Y, pch=19, cex=0.2, col=ifelse(Y >= 50 & useNorth, "black", "grey")))
with(dd[dd$validation,], points(X, Y, pch=19, cex=0.4, col=ifelse(PCODE == "BU_BG", 4, 2)))

cn2 <- c(cn, get_terms(mods_veg, "list"), "vegw", "vegca", "ABMIsite")
DAT <- dd[dd$Y >= 50 & dd$useNorth & !dd$validation, cn2]

YY <- yy[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= NMIN & colnames(YY) %in% colnames(off)]
OFF <- off[rownames(DAT), colnames(YY)]
mods <- mods_veg
#mods$SSH <- NULL
SSH <- SSH_veg[rownames(DAT),]

BB <- pbapply::pbsapply(1:B, bfun, DAT$SS, DAT$BLOCK_XYT)
nrow(DAT)
max(BB)
(lu <- length(unique(as.numeric(BB))))
stopifnot(all(BB <= nrow(DAT)))
stopifnot(lu <= nrow(DAT))
nrow(BB)

z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="vegw", ssh_class="vegca", ssh_fit="Space")
z$timer
cat("Estimate for", ncol(YY), "species and", B, "runs is", ceiling(unname(ncol(YY)*B*z$timer[3])/(60*60)), "hrs\n")

DATv <- dd[dd$validation, cn2]
YYv <- yy[rownames(DATv),colnames(YY)]
OFFv <- off[rownames(DATv), colnames(YY)]
SSHv <- SSH_veg[rownames(DATv),]


save(DAT, YY, OFF, SSH, BB, mods, DATv, YYv, OFFv, SSHv,
    file="d:/abmi/AB_data_v2018/data/analysis/birds/data/ab-birds-validation-2019-01-30.RData")

#' The End



## set root directory
ROOT <- "~/Dropbox/courses/st-johns-2017"
library(raster)
library(rgdal)
library(rgeos)
library(sp)

## load natural regions shape file
## dissolve polygons
setwd(file.path(ROOT, "data", "NatRegAB"))
AB <- readOGR(".",
    "Natural_Regions_Subregions_of_Alberta") # rgdal
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province


COL <- c('#e6f5c9','#f4cae4','#b3e2cd','#fff2ae','#fdcdac','#cbd5e8')

## sampling grid as spatial points data frame
x <- dd[,c("X","Y")]
coordinates(x) <- c("X", "Y") # string
proj4string(x) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
x <- spTransform(x, proj4string(AB))

png("aa.png", height=2000, width=1200, bg=NA)
op <- par(mar=c(0,0,1,0))
plot(ABnr, col=COL, border=COL)
plot(x,add=TRUE)
par(op)
dev.off()

## export all the BG data (no temporal resampling/subsetting)

source("~/repos/abmianalytics/birds/models-veg.R")
setdiff(get_terms(mods_veg, "list"), colnames(dd))
rm(DAT, YY, OFF, OFFmean, SSH, BB, mods)

dd$isBG <- dd$PCODE == "BU_BG"

cn2 <- c(cn, get_terms(mods_veg, "list"), "vegw", "vegca")
DAT <- dd[dd$isBG, cn2]

DATv <- DAT

bgu <- read.csv("d:/abmi/AB_data_v2018/data/analysis/birds/validation-BGgroups.csv")
bgu$xv_yv <- as.factor(paste0(bgu$xv, "_", bgu$yv))

Grain <- expand.grid(xv=1:10, yv=1:10)
rownames(Grain) <- paste0(Grain$xv, "_", Grain$yv)
Grain$x0 <- factor(rep("x", 100), c("x", ""))
Grain$x1 <- Grain$x0
Grain$x1[] <- ifelse(Grain$xv %in% c(1,3,5,7,9) & Grain$yv %in% c(1,3,5,7,9), "x", "")
Grain$x2 <- Grain$x0
Grain$x2[] <- ifelse(Grain$xv %in% c(1,4,7,10) & Grain$yv %in% c(1,4,7,10), "x", "")
Grain$x3 <- Grain$x0
Grain$x3[] <- ifelse(Grain$xv %in% c(1,5,9) & Grain$yv %in% c(1,5,9), "x", "")
Grain$x4 <- Grain$x0
Grain$x4[] <- ifelse(Grain$xv %in% c(1,10) & Grain$yv %in% c(1,10), "x", "")
#with(Grain, plot(xv, yv, pch=c(19, 21)[as.integer(x4)]))
Grain <- Grain[match(bgu$xv_yv, rownames(Grain)),]

op <- par(mfrow=c(2,3))
with(Grain, plot(xv, yv, pch=c(19, 21)[as.integer(x0)], main="0"))
with(Grain, plot(xv, yv, pch=c(19, 21)[as.integer(x1)], main="1"))
with(Grain, plot(xv, yv, pch=c(19, 21)[as.integer(x2)], main="2"))
with(Grain, plot(xv, yv, pch=c(19, 21)[as.integer(x3)], main="3"))
with(Grain, plot(xv, yv, pch=c(19, 21)[as.integer(x4)], main="4"))
par(op)

bgu$x0 <- Grain$x0
bgu$x1 <- Grain$x1
bgu$x2 <- Grain$x2
bgu$x3 <- Grain$x3
bgu$x4 <- Grain$x4

## extent
bgu$SS2 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, g2, sep="::", drop=TRUE)))
tmp <- table(bgu$SS2)
for (i in names(tmp))
    if (tmp[i] < 4)
        bgu$SS2[bgu$SS2 == i] <- ""
bgu$SS2 <- as.factor(bgu$SS2)
bgu$SS3 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, g3, sep="::", drop=TRUE)))
tmp <- table(bgu$SS3)
for (i in names(tmp))
    if (tmp[i] < 9)
        bgu$SS3[bgu$SS3 == i] <- ""
bgu$SS3 <- as.factor(bgu$SS3)
bgu$SS4 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, g4, sep="::", drop=TRUE)))
tmp <- table(bgu$SS4)
for (i in names(tmp))
    if (tmp[i] < 16)
        bgu$SS4[bgu$SS4 == i] <- ""
bgu$SS4 <- as.factor(bgu$SS4)
bgu$SS5 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, g5, sep="::", drop=TRUE)))
tmp <- table(bgu$SS5)
for (i in names(tmp))
    if (tmp[i] < 25)
        bgu$SS5[bgu$SS5 == i] <- ""
bgu$SS5 <- as.factor(bgu$SS5)
bgu$SS10 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, g10, sep="::", drop=TRUE)))
bgu$SS10 <- as.factor(bgu$SS10)
levels(bgu$SS10) <- c(levels(bgu$SS10), "")

## grain size

bgu$xx0 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, x0, sep="::", drop=TRUE)))
bgu$xx0[bgu$x0 == ""] <- ""
bgu$xx0 <- as.factor(bgu$xx0)
levels(bgu$xx0) <- c(levels(bgu$xx0), "")

bgu$xx1 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, x1, sep="::", drop=TRUE)))
bgu$xx1[bgu$x1 == ""] <- ""
bgu$xx1 <- as.factor(bgu$xx1)

bgu$xx2 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, x2, sep="::", drop=TRUE)))
bgu$xx2[bgu$x2 == ""] <- ""
bgu$xx2 <- as.factor(bgu$xx2)

bgu$xx3 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, x3, sep="::", drop=TRUE)))
bgu$xx3[bgu$x3 == ""] <- ""
bgu$xx3 <- as.factor(bgu$xx3)

bgu$xx4 <- as.character(with(bgu, interaction(ProjectID, Cluster, SITE, x4, sep="::", drop=TRUE)))
bgu$xx4[bgu$x4 == ""] <- ""
bgu$xx4 <- as.factor(bgu$xx4)


## extent: sampling intensity (/ unit area) stays, extent grows grows
DATv$EX2x2 <- bgu$SS2[match(DATv$SS, bgu$SS)]
DATv$EX2x2[is.na(DATv$EX2x2)] <- ""
DATv$EX3x3 <- bgu$SS3[match(DATv$SS, bgu$SS)]
DATv$EX3x3[is.na(DATv$EX3x3)] <- ""
DATv$EX4x4 <- bgu$SS4[match(DATv$SS, bgu$SS)]
DATv$EX4x4[is.na(DATv$EX4x4)] <- ""
DATv$EX5x5 <- bgu$SS5[match(DATv$SS, bgu$SS)]
DATv$EX5x5[is.na(DATv$EX5x5)] <- ""
DATv$EX10x10 <- bgu$SS10[match(DATv$SS, bgu$SS)]
DATv$EX10x10[is.na(DATv$EX10x10)] <- ""

## 'grain' size: extent stays save but sampling intensity (/ unit area) changes
DATv$GR2x2 <- bgu$xx4[match(DATv$SS, bgu$SS)]
DATv$GR2x2[is.na(DATv$GR2x2)] <- ""
DATv$GR3x3 <- bgu$xx3[match(DATv$SS, bgu$SS)]
DATv$GR3x3[is.na(DATv$GR3x3)] <- ""
DATv$GR4x4 <- bgu$xx2[match(DATv$SS, bgu$SS)]
DATv$GR4x4[is.na(DATv$GR4x4)] <- ""
DATv$GR5x5 <- bgu$xx1[match(DATv$SS, bgu$SS)]
DATv$GR5x5[is.na(DATv$GR5x5)] <- ""
DATv$GR10x10 <- bgu$xx0[match(DATv$SS, bgu$SS)]
DATv$GR10x10[is.na(DATv$GR10x10)] <- ""

bgu$bgid <- as.factor(paste0("BGID", bgu$SITE))
DATv$BGID <- bgu$bgid[match(DATv$SS, bgu$SS)]
bgu$bgxy <- as.factor(paste0("BGID", bgu$SITE, "_", bgu$xv_yv))
DATv$BGXY <- bgu$bgxy[match(DATv$SS, bgu$SS)]
DATv$XXYY <- bgu$xv_yv[match(DATv$SS, bgu$SS)]

## randomly pick one visit for the spatial groups
DATv <- DATv[sample(nrow(DATv)),]
#DATv$RND1 <- ifelse(duplicated(DATv$SS), 0, 1)
DATv <- DATv[order(rownames(DATv)),]

DATv <- DATv[DATv$YEAR >= 2015 & !is.na(DATv$BGID),]
#data.frame(nl=apply(DATv,2,function(z) nlevels(as.factor(z))), nl2=apply(droplevels(DATv),2,function(z) nlevels(as.factor(z))))

## dealing with visits
set.seed(1)
DATv$VISIT <- NA
DATv$VISITRND <- NA
for (i in levels(DATv$BGXY)) {
    ii <- DATv$BGXY == i
    tmp <- DATv[ii,,drop=FALSE]
    tmp$VISIT <- seq_len(nrow(tmp))
    tmp$VISITRND <- tmp$VISIT
    if (nrow(tmp) > 1)
        tmp$VISITRND <- sample(tmp$VISIT)
    DATv$VISIT[ii] <- tmp$VISIT
    DATv$VISITRND[ii] <- tmp$VISITRND
}
## temporal pools: seqential visits
DATv$TS1 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TS1[DATv$VISIT > 1] <- ""
DATv$TS2 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TS2[DATv$VISIT > 2] <- ""
DATv$TS3 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TS3[DATv$VISIT > 3] <- ""
DATv$TS4 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TS4[DATv$VISIT > 4] <- ""
## temporal pools: random visits
DATv$TR1 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TR1[DATv$VISITRND > 1] <- ""
DATv$TR2 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TR2[DATv$VISITRND > 2] <- ""
DATv$TR3 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TR3[DATv$VISITRND > 3] <- ""
DATv$TR4 <- paste0(DATv$BGXY, "_", DATv$VISIT)
DATv$TR4[DATv$VISITRND > 4] <- ""

DATv$EX2x2[DATv$VISITRND > 1] <- ""
DATv$EX3x3[DATv$VISITRND > 1] <- ""
DATv$EX4x4[DATv$VISITRND > 1] <- ""
DATv$EX5x5[DATv$VISITRND > 1] <- ""
DATv$EX10x10[DATv$VISITRND > 1] <- ""

DATv$GR2x2[DATv$VISITRND > 1] <- ""
DATv$GR3x3[DATv$VISITRND > 1] <- ""
DATv$GR4x4[DATv$VISITRND > 1] <- ""
DATv$GR5x5[DATv$VISITRND > 1] <- ""
DATv$GR10x10[DATv$VISITRND > 1] <- ""


YY <- yy[rownames(DATv),]
YY <- YY[,colSums(YY>0) >= NMIN & colnames(YY) %in% colnames(off)]
OFF <- off[rownames(DATv), colnames(YY)]
mods <- mods_veg
#mods$SSH <- NULL
SSH <- SSH_veg[rownames(DATv),]

YYv <- YY
OFFv <- OFF
SSHv <- SSH

save(#DAT, YY, OFF, SSH, BB,
    mods, DATv, YYv, OFFv, SSHv,
    file="d:/abmi/AB_data_v2018/data/analysis/birds/data/ab-birds-biggrids-2019-05-13.RData")


## FORSITE mixedwood redefinitions

#tmp <- read.csv("d:/bam/BAM_data_v2019/forsite/stands2PSpoints.csv")
tmp <- read.csv("d:/bam/BAM_data_v2019/forsite/stands2ALLpoints.csv")

compare_sets(dd$PKEY, tmp$PKEY)

dd$mw <- tmp$StandType[match(dd$PKEY, tmp$PKEY)]
table(dd$mw, dd$vegc)
table(tmp$StandType)
data.frame(table(dd$mw, useNA="a"))
table(dd$mw, dd$vegc, useNA="a")

dd$mw[dd$mw == "Unidentified or no data"] <- NA
dd$mw <- droplevels(dd$mw)
dd <- dd[!is.na(dd$mw),]

table(dd$mw, grepl("Mixedwood", as.character(dd$mw), fixed=TRUE))
dd$isMix <- ifelse(grepl("Mixedwood", as.character(dd$mw), fixed=TRUE), 1, 0)
dd$isWSpruce <- ifelse(dd$mw %in% c("pureWhiteSpruce"), 1, 0)
dd$isPine <- ifelse(dd$mw %in% c("purePine"), 1, 0)
dd$isBSpruce <- ifelse(dd$mw %in% c("pureBlackSpruce"), 1, 0)
dd$isLarch <- 0
dd$isBSLarch <- pmax(dd$isLarch, dd$isBSpruce)
dd$isUpCon <- ifelse(dd$mw %in% c("pureConifer_allcon", "pureConifer_Pl",
    "pureConifer_Sx", "purePine", "pureWhiteSpruce"), 1, 0)
dd$isCon <- ifelse(dd$mw %in% c("pureConifer_allcon", "pureConifer_Pl",
    "pureConifer_Sx", "purePine", "pureWhiteSpruce", "pureBlackSpruce"), 1, 0)
dd$vegc <- dd$mw
dd$vegw <- 1


source("~/repos/abmianalytics/birds/models-veg.R")
#mods_veg$Hab <- list(.~. + mw)


setdiff(get_terms(mods_veg, "list"), colnames(dd))
rm(DAT, YY, OFF, OFFmean, SSH, BB, mods)

cn2 <- c(cn, get_terms(mods_veg, "list"), "vegw", "vegca")
DAT <- dd[, cn2]
YY <- yy[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= NMIN]
#YY <- YY[,colSums(groupSums(YY, 1, DAT$SS) > 0) >= NMIN]
OFF <- off[rownames(DAT), intersect(colnames(off), colnames(YY))]
YY <- YY[,colnames(OFF)]
mods <- mods_veg
mods$SSH <- NULL
mods$HF <- NULL
mods$ARU <- NULL
#SSH <- SSH_veg[rownames(DAT),]

BB <- pbapply::pbsapply(1:B, bfun, DAT$SS, DAT$BLOCK_XYT)
nrow(DAT)
max(BB)
(lu <- length(unique(as.numeric(BB))))
stopifnot(all(BB <= nrow(DAT)))
stopifnot(lu <= nrow(DAT))
nrow(BB)

z <- .run_path1(1, "AMRO", mods, CAICalpha=1)
z$timer
cat("Estimate for", ncol(YY), "species and", B, "runs is", ceiling(unname(ncol(YY)*B*z$timer[3])/(60*60)), "hrs\n")

save(DAT, YY, OFF, BB, mods,
    file="d:/abmi/AB_data_v2018/data/analysis/birds/data/ab-birds-mixedwood-2019-10-31.RData")
