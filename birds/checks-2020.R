## comparing 2018 (1) and 2020 (2) version of northern estimates

library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT1 <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
ROOT2 <- "d:/abmi/AB_data_v2020/data/analysis/species/birds" # change this bit


en1 <- new.env()
load(file.path(ROOT1, "data", "ab-birds-north-2019-01-30.RData"), envir=en1)

en2 <- new.env()
load(file.path(ROOT2, "data", "ab-birds-north-2020-09-23.RData"), envir=en2)

Xn1 <- get_model_matrix(en1$DAT, en1$mods)
Xn2 <- get_model_matrix(en2$DAT, en2$mods)

Xage1 <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage1) <- colnames(Xn1)[match(colnames(Xage1), make.names(colnames(Xn1)))]
Xage2 <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v2020.csv"))
colnames(Xage2) <- colnames(Xn2)[match(colnames(Xage2), make.names(colnames(Xn2)))]

SPP <- intersect(colnames(en1$YY), colnames(en2$YY))

get_coef_north <- function(resn, STAGE="ARU", subset=NULL, new=TRUE, ...) {
    if (new) {
        Xn <- Xn2
        Xage <- Xage2
        en <- en2
    } else {
        Xn <- Xn1
        Xage <- Xage1
        en <- en1
    }
    OK <- !sapply(resn, inherits, "try-error")
    spp <- resn[[which(OK)[1]]]$species
    estn <- suppressWarnings(get_coef(resn, Xn, stage=STAGE, na.out=FALSE))
    if (is.null(subset))
        subset <- seq_len(nrow(estn))
    estn <- estn[subset,,drop=FALSE]

    mu <- Xage %*% t(estn[,colnames(Xage),drop=FALSE])
    lam1 <- exp(mu)
    lam1 <- lam1[!grepl("9", rownames(lam1)),,drop=FALSE]
    lam1 <- t(apply(lam1, 1, quantile, c(0.5, 0.05, 0.95)))
    lamCC <- lam1[grepl("CC", rownames(lam1)),,drop=FALSE]

    MOD <- c("ROAD", "mWell", "mSoft",
        "mEnSft", "mTrSft", "mSeism", "CMETHODSM", "CMETHODRF")
    Z <- exp(estn[,MOD,drop=FALSE])
    isSoft <- estn[,"mSoft"] != 0 & estn[,"mEnSft"] == 0
    #isSoft2 <- get_mid(resn)[,"Contrast"] == 3
    estn[isSoft,"mEnSft"] <- estn[isSoft,"mSoft"]
    estn[isSoft,"mTrSft"] <- estn[isSoft,"mSoft"]
    estn[isSoft,"mSeism"] <- estn[isSoft,"mSoft"]
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
        "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05,
        "CMETHODSM"=1, "CMETHODRF"=1)
    for (i in MOD)
        Z[,i] <- linexp(1, estn[,i], pm[i])

    HFc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")

    # not in HFc and not forestry!
    Xn2 <- Xn[en$DAT$mWell > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nWell <- nrow(Xn2)
    lamWell <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    estWell <- quantile(lamWell * Z[,"mWell"], c(0.5, 0.05, 0.95))

    Xn2 <- Xn[en$DAT$mEnSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nEnSoft <- nrow(Xn2)
    lamEnSft <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    estEnSft <- quantile(lamEnSft * Z[,"mEnSft"], c(0.5, 0.05, 0.95))

    ## TrSft incorporates ROAD effect as well?
    Xn2 <- Xn[en$DAT$mTrSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nTrSoft <- nrow(Xn2)
    lamTrSft <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    #estTrSft <- quantile(lamTrSft * Z[,"mTrSft"] * Z[,"ROAD"], c(0.5, 0.05, 0.95))
    estTrSft <- quantile(lamTrSft * Z[,"mTrSft"], c(0.5, 0.05, 0.95))

    ESMAX <- apply(rbind(lamCC[endsWith(rownames(lamCC), "_0-10"),], EnSft=estEnSft, TrSft=estTrSft), 2, max)
    Xn2 <- Xn[en$DAT$mSeism > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nSeism <- nrow(Xn2)
    lamSeism <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    estSeism <- quantile(lamSeism * Z[,"mSeism"], c(0.5, 0.05, 0.95))
    if (estSeism[1] > ESMAX[1])
        estSeism <- ESMAX

    lam1 <- rbind(lam1,
        #UrbInd=estUI,
        SoftLin=(nSeism * estSeism + nEnSoft * estEnSft) / (nSeism + nEnSoft),
        HardLin=estTrSft)
    ## add here what's needed for final output as defined (seismic, ensoft, trsoft, hardlin, well)
    colnames(lam1) <- c("Estimate", "Lower", "Upper")
    lam1
}

COEFS1 <- list()
COEFS2 <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    resn1 <- load_species(file.path(ROOT1, "out", "north", paste0(spp, ".RData")))
    resn2 <- load_species(file.path(ROOT2, "out", "north", paste0(spp, ".RData")))
    tmp1 <- try(get_coef_north(resn1, STAGE="ARU", new=FALSE))
    tmp2 <- try(get_coef_north(resn2, STAGE="ARU", new=TRUE))
    if (!inherits(tmp1, "try-error"))
        COEFS1[[spp]] <- tmp1
    if (!inherits(tmp2, "try-error"))
        COEFS2[[spp]] <- tmp2
}

compare_sets(names(COEFS1), names(COEFS2))

m <- structure(list(l1 = c("SpruceR", "Spruce1", "Spruce2", "Spruce3",
    "Spruce4", "Spruce5", "Spruce6", "Spruce7", "Spruce8", "PineR",
    "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6", "Pine7",
    "Pine8", "DecidR", "Decid1", "Decid2", "Decid3", "Decid4", "Decid5",
    "Decid6", "Decid7", "Decid8", "MixedwoodR", "Mixedwood1", "Mixedwood2",
    "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7",
    "Mixedwood8", "BSprR", "BSpr1", "BSpr2", "BSpr3", "BSpr4", "BSpr5",
    "BSpr6", "BSpr7", "BSpr8", "LarchR", "Larch1", "Larch2", "Larch3",
    "Larch4", "Larch5", "Larch6", "Larch7", "Larch8", "GrassHerb",
    "Shrub", "GraminoidFen", "Marsh", "RoughP", "TameP", "Crop",
    "Industrial", "Mine", "Rural", "Urban", "CCSpruceR", "CCSpruce1",
    "CCSpruce2", "CCSpruce3", "CCSpruce4", "CCPineR", "CCPine1",
    "CCPine2", "CCPine3", "CCPine4", "CCDecidR", "CCDecid1", "CCDecid2",
    "CCDecid3", "CCDecid4", "CCMixedwoodR", "CCMixedwood1", "CCMixedwood2",
    "CCMixedwood3", "CCMixedwood4", "Swamp", "Swamp", "Swamp", "Swamp",
    "SoftLin", "HardLin"), l2 = c("SpruceR", "Spruce1", "Spruce2",
    "Spruce3", "Spruce4", "Spruce5", "Spruce6", "Spruce7", "Spruce8",
    "PineR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
    "Pine7", "Pine8", "DecidR", "Decid1", "Decid2", "Decid3", "Decid4",
    "Decid5", "Decid6", "Decid7", "Decid8", "MixedwoodR", "Mixedwood1",
    "Mixedwood2", "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6",
    "Mixedwood7", "Mixedwood8", "TreedBogR", "TreedBog1", "TreedBog2",
    "TreedBog3", "TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7",
    "TreedBog8", "TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3",
    "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8",
    "GrassHerb", "Shrub", "GraminoidFen", "Marsh", "RoughP", "TameP",
    "Crop", "Industrial", "Mine", "Rural", "Urban", "CCSpruceR",
    "CCSpruce1", "CCSpruce2", "CCSpruce3", "CCSpruce4", "CCPineR",
    "CCPine1", "CCPine2", "CCPine3", "CCPine4", "CCDecidR", "CCDecid1",
    "CCDecid2", "CCDecid3", "CCDecid4", "CCMixedwoodR", "CCMixedwood1",
    "CCMixedwood2", "CCMixedwood3", "CCMixedwood4", "ShrubbySwamp",
    "ShrubbyBog", "ShrubbyFen", "TreedSwamp", "SoftLin", "HardLin"
    )), class = "data.frame", row.names = c(NA, -91L))

spp <- "ALFL"

V <- list()
for (spp in SPP) {
    l1 <- data.frame(COEFS1[[spp]][m$l1,])
    l1$r <- l1$Upper - l1$Lower
    l2 <- data.frame(COEFS2[[spp]][m$l2,])
    l2$r <- l2$Upper - l2$Lower

    V[[spp]] <- c(COR = c(Estimate=cor(l1$Estimate, l2$Estimate, method="spearman"),
        Range=cor(l1$r, l2$r, method="spearman")),
        SCALE = c(Estimate=median(l2$Estimate/l1$Estimate), Range=median(l2$r/l1$r)))

}

V <- data.frame(do.call(rbind, V))
summary(V)
VV <- V[rowSums(is.na(V))==0 & is.finite(V$SCALE.Estimate),]

hist(VV$COR.Estimate)
hist(VV$SCALE.Estimate)

head(VV[order(VV$COR.Estimate),], 20)
head(VV[order(VV$SCALE.Estimate),], 20)

with(V, plot(COR.Estimate ~ SCALE.Estimate, ylim=c(-1, 1), xlim=c(0,2)))
abline(h=0.8, v=c(0.8, 1.2), lty=2)

## ESTIMATES -----------------
addmargins(with(V, table(COR=cut(COR.Estimate, c(-1, 0, 0.5, 0.8, 1)),
    SCALE=cut(SCALE.Estimate, c(0, 0.5, 0.8, 1.2, 2, Inf)), useNA="a")))
## RANGE (upper - lower) -----------------
addmargins(with(V, table(COR=cut(COR.Range, c(-1, 0, 0.5, 0.8, 1)),
    SCALE=cut(SCALE.Range, c(0, 0.5, 0.8, 1.2, 2, Inf)), useNA="a")))
