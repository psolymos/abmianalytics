#' Preliminiary coefs for all-in-one reporting
#' Non-bir taxa comes here
#'

ROOT <- "s:/AB_data_v2020/Results/Results from Ermias/Preliminary/"

TAXA <- c("lichens", "mites", "mosses", "vplants")
taxon <- "mites"

COEFS <- list()
B <- 100

for (taxon in TAXA) {

    if (taxon == "lichens") {
        fn <- "Lichen_North Coefficient tables 2020.Rdata"
        fs <- "Lichen_South Coefficient tables 2020.Rdata"
    }
    if (taxon == "mites") {
        fn <- "Mite_North Coefficient tables 2020.Rdata"
        fs <- "Mite_South Coefficient tables 2020.Rdata"
    }
    if (taxon == "mosses") {
        fn <- "Moss_North Coefficient tables 2020.Rdata"
        fs <- "Moss_South Coefficient tables 2020.Rdata"
    }
    if (taxon == "vplants") {
        fn <- "VPlant_North Coefficient tables 2020.Rdata"
        fs <- "VPlant_South Coefficient tables 2020.Rdata"
    }

    #' load North and south coefs into their respective environments
    en <- new.env()
    load(paste0(ROOT, fn), envir=en)
    es <- new.env()
    load(paste0(ROOT, fs), envir=es)

    #' # North
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    cfn <- en$Coef
    colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
    cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=cfn[,"Grass"])
    cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass"))]
    xCFn <- cbind(qlogis(cfn), en$Sc.coef)
    #' array: bootstrap x species x coefs
    CFn <- array(0, c(nrow(xCFn), ncol(xCFn), B),
        dimnames=list(rownames(xCFn), colnames(xCFn), NULL))
    for (i in 1:B) {
        cfn <- en$Coef
        colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
        cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=cfn[,"Grass"])
        cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass"))]
        ses <- en$Coef.se
        colnames(ses) <- gsub("BlackSpruce", "TreedBog", colnames(ses))
        ses <- cbind(ses,
            TreedFenR=ses[,"TreedFen"],
            TreedFen1=ses[,"TreedFen"],
            TreedFen2=ses[,"TreedFen"],
            TreedFen3=ses[,"TreedFen"],
            TreedFen4=ses[,"TreedFen"],
            TreedFen5=ses[,"TreedFen"],
            TreedFen6=ses[,"TreedFen"],
            TreedFen7=ses[,"TreedFen"],
            TreedFen8=ses[,"TreedFen"],
            Rural=ses[,"UrbInd"],
            Urban=ses[,"UrbInd"],
            Industrial=ses[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=ses[,"Grass"])
        ses <- ses[,!(colnames(ses) %in% c("TreedFen", "UrbInd", "Grass"))]
        if (i > 1) {
            se <- ses
            se[se > 5] <- 5
            eps <- rnorm(length(se), 0, se)
        } else {
            eps <- 0
        }
        CFn[,,i] <- cbind(qlogis(cfn) + eps, en$Sc.coef)
    }


    #' # South
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    cfs <- es$Coef
    cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
    cfs <- cfs[,colnames(cfs) != "UrbInd"]
    xCFs <- cbind(qlogis(cfs), pAspen=es$Coef.pAspen, es$Sc.coef)
    #' array: bootstrap x species x coefs
    CFs <- array(0, c(nrow(xCFs), ncol(xCFs), B),
        dimnames=list(rownames(xCFs), colnames(xCFs), NULL))
    for (i in 1:B) {
        cfs <- es$Coef
        cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        cfs <- cfs[,colnames(cfs) != "UrbInd"]
        ses <- es$Coef.se
        ses <- cbind(ses,
            Rural=ses[,"UrbInd"],
            Urban=ses[,"UrbInd"],
            Industrial=ses[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        ses <- ses[,colnames(ses) != "UrbInd"]
        if (i > 1) {
            se <- ses
            se[se > 5] <- 5
            eps <- rnorm(length(se), 0, se)
        } else {
            eps <- 0
        }
        CFs[,,i] <- cbind(qlogis(cfs) + eps, pAspen=es$Coef.pAspen, es$Sc.coef)
    }
    #range(CFn)
    #range(CFs)
    CFn[CFn < -10^4] <- -10^4
    CFn[CFn > 10^4] <- 10^4
    CFs[CFs < -10^4] <- -10^4
    CFs[CFs > 10^4] <- 10^4

    COEFS[[taxon]] <- list(north=CFn, south=CFs)

}

save(COEFS, file="s:/AB_data_v2020/Results/COEFS-EA.RData")

#' Birds
#'
#' Organize coefs: compare land cover classes N/S
#' Organize spclim standardizations
#' Organize marginal & joint versions

library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds"

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)
bt <- ee$tax
rm(ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2020-09-23.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2020-09-23.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v2020.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

load("s:/AB_data_v2020/Results/COEFS-EA.RData")

get_coef_north <- function(resn, STAGE="ARU", subset=1, ...) {
    estn <- suppressWarnings(get_coef(resn, Xn, stage=STAGE, na.out=FALSE))
    estn <- estn[subset,]

    mu <- drop(Xage %*% estn[colnames(Xage)])
    lam1 <- exp(mu)
    lam1 <- lam1[!grepl("9", names(lam1))]
    lamCC <- lam1[grepl("CC", names(lam1))]

    MOD <- c("ROAD", "mWell", "mSoft",
        "mEnSft", "mTrSft", "mSeism")
    Z <- exp(estn[MOD])
    isSoft <- estn["mSoft"] != 0 & estn["mEnSft"] == 0
    #isSoft2 <- get_mid(resn)[,"Contrast"] == 3
    if (isSoft) {
        estn["mEnSft"] <- estn["mSoft"]
        estn["mTrSft"] <- estn["mSoft"]
        estn["mSeism"] <- estn["mSoft"]
    }
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
        "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05)
    for (i in MOD)
        Z[i] <- linexp(1, estn[i], pm[i])

    HFc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")

    FUN <- mean

    # not in HFc and not forestry!
    Xn2 <- Xn[en$DAT$mWell > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nWell <- nrow(Xn2)
    lamWell <- exp(Xn2 %*% estn[colnames(Xn2)])
    estWell <- FUN(lamWell * Z["mWell"])

    Xn2 <- Xn[en$DAT$mEnSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nEnSoft <- nrow(Xn2)
    lamEnSft <- exp(Xn2 %*% estn[colnames(Xn2)])
    estEnSft <- FUN(lamEnSft * Z["mEnSft"])

    ## TrSft incorporates ROAD effect as well?
    Xn2 <- Xn[en$DAT$mTrSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nTrSoft <- nrow(Xn2)
    lamTrSft <- exp(Xn2 %*% estn[colnames(Xn2)])
    #estTrSft <- quantile(lamTrSft * Z[,"mTrSft"] * Z[,"ROAD"], c(0.5, 0.05, 0.95))
    estTrSft <- FUN(lamTrSft * Z["mTrSft"])

    ESMAX <- max(lamCC[endsWith(names(lamCC), "R")], EnSft=estEnSft, TrSft=estTrSft)
    Xn2 <- Xn[en$DAT$mSeism > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nSeism <- nrow(Xn2)
    lamSeism <- exp(Xn2 %*% estn[colnames(Xn2)])
    estSeism <- FUN(lamSeism * Z["mSeism"])
    if (estSeism > ESMAX)
        estSeism <- ESMAX

    ## UrbInd
    #Xn2 <- Xn[en$DAT$vegc %ni% c("Industrial", "Mine", "Rural", "Urban"), colnames(Xage)]
    #nUI <- nrow(Xn2)
    #lamUI <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    #estUI <- quantile((nUI * lamUI + nWell * lamWell) / (nUI + nWell),
    #    c(0.5, 0.05, 0.95))

    lam1 <- c(lam1,
        Wellsites=estWell,
        EnSeismic=estSeism,
        EnSoftLin=estEnSft,
        TrSoftLin=estTrSft,
        HardLin=0,
        Water=0,
        MineV=unname(lam1["Mine"]))
    lam1["Mine"] <- 0
    names(lam1) <- gsub("Spruce", "WhiteSpruce", names(lam1))
    names(lam1) <- gsub("Decid", "Deciduous", names(lam1))

        #UrbInd=estUI,
        #SoftLin=(nSeism * estSeism + nEnSoft * estEnSft) / (nSeism + nEnSoft),
        #HardLin=estTrSft)
    lam1
}

get_coef_south <- function(ress, STAGE="ARU", subset=1, ...) {
    ests <- suppressWarnings(get_coef(ress, Xs, stage=STAGE, na.out=FALSE))
    ests <- ests[subset,]
    LCC <- c("soilcBlowout", "soilcClaySub", "soilcCrop",
        "soilcIndustrial", "soilcMine", "soilcOther", "soilcRapidDrain",
        "soilcRoughP", "soilcRural", "soilcSandyLoam", "soilcTameP",
        "soilcThinBreak", "soilcUrban")
    muLCC <- c(ests[1], ests[1]+ests[LCC])
    names(muLCC) <- levels(es$DAT$soilc)
    lamLCC0 <- exp(muLCC) # nontreed
    o <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other",
        "RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")
    lamLCC0 <- lamLCC0[o]

    MOD <- c("ROAD", "mWell", "mSoft")
    Z <- exp(ests[MOD])
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2)
    for (i in MOD)
        Z[i] <- linexp(1, ests[i], pm[i])
    lamMOD <- Z
    names(lamMOD) <- c("Road", "Well", "Soft")

    HFc <- c("RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")

    FUN <- mean

    Xs2 <- Xs[es$DAT$mSoft > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nSoft <- nrow(Xs2)
    lamSoft <- exp(Xs2 %*% ests[colnames(Xs2)])
    estSoft <- FUN(lamSoft * Z["mSoft"])

    Xs2 <- Xs[es$DAT$mWell > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nWell <- nrow(Xs2)
    lamWell <- exp(Xs2 %*% ests[colnames(Xs2)])
    estWell <- FUN(lamWell * Z["mWell"])

    Xs2 <- Xs[es$DAT$ROAD > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nROAD <- nrow(Xs2)
    lamROAD <- exp(Xs2 %*% ests[colnames(Xs2)])
    estROAD <- FUN(lamROAD * Z["ROAD"])

    lam0 <- c(lamLCC0,
        EnSeismic=estSoft,
        EnSoftLin=estSoft,
        TrSoftLin=estSoft,
        HardLin=estROAD,
        Wellsites=estWell,
        MineV=unname(lamLCC0["Mine"]),
        Water=0)
    lam0["Mine"] <- 0
    attr(lam0, "pAspen") <- ests["pAspen"]
    lam0
}

## all species
SPPn <- substr(list.files(file.path(ROOT, "out", "north")), 1, 4)
names(SPPn) <- bt[SPPn, "sppid"]
SPPs <- substr(list.files(file.path(ROOT, "out", "south")), 1, 4)
names(SPPs) <- bt[SPPs, "sppid"]
## same B as for other taxa
B <- dim(COEFS[[1]]$north)[3]


cfn <- list()
for (spp in SPPn) {
    cat(spp, "\n")
    flush.console()

    res <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    est1 <- suppressWarnings(get_coef(res, Xn, stage="ARU", na.out=FALSE))
    est2 <- suppressWarnings(get_coef(res, Xn, stage="Space", na.out=FALSE))

    BB <- min(B, nrow(est1))
    cf1 <- sapply(1:BB, function(i) get_coef_north(res, STAGE="ARU", subset=i))
    cf2 <- sapply(1:BB, function(i) get_coef_north(res, STAGE="Space", subset=i))

    cfn[[spp]] <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}
save(cfn, file="s:/AB_data_v2020/Results/BIRDS-North.RData")

cfs <- list()
for (spp in SPPs) {
    cat(spp, "\n")
    flush.console()

    res <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    est1 <- suppressWarnings(get_coef(res, Xs, stage="ARU", na.out=FALSE))
    est2 <- suppressWarnings(get_coef(res, Xs, stage="Space", na.out=FALSE))

    BB <- min(B, nrow(est1))
    cf1 <- sapply(1:BB, function(i) get_coef_south(res, STAGE="ARU", subset=i))
    cf2 <- sapply(1:BB, function(i) get_coef_south(res, STAGE="Space", subset=i))

    cfs[[spp]] <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}

save(cfs, file="s:/AB_data_v2020/Results/BIRDS-South.RData")

## -- load coefs and make figures

load("s:/AB_data_v2020/Results/COEFS-EA.RData")
load("s:/AB_data_v2020/Results/BIRDS-North.RData")
load("s:/AB_data_v2020/Results/BIRDS-South.RData")

## organizing bird coefs

cns <- rownames(COEFS$mites$south[1,,])
cns0 <- cns[1:(which(cns=="Intercept")-2)]
cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]


cnb <- c("pWater_KM", "pWater2_KM", "xPET",
    "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT", "xMCMT", "xY", "xX",
    "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY")

# south
spp <- "ALFL"
B <- min(100, nrow(cfs[[spp]]$estSpace))
cfsm <- rbind(log(cfs[[spp]]$coefARU), pAspen=cfs[[spp]]$estARU[1:B,"pAspen"])
cfsj <- rbind(log(cfs[[spp]]$coefSpace), t(cfs[[spp]]$estSpace[1:B,c("pAspen", cnb)]))
cfsm[cfsm > 10^4] <- 10^4
cfsm[cfsm < -10^4] <- -10^4
cfsj[cfsj > 10^4] <- 10^4
cfsj[cfsj < -10^4] <- -10^4
if (ncol(cfsm) < 100) {
    b <- sample.int(ncol(cfsm), 100-ncol(cfsm), replace=TRUE)
    cfsm <- cbind(cfsm, cfsm[,b])
    cfsj <- cbind(cfsj, cfsj[,b])
}
compare_sets(rownames(cfsm),dimnames(COEFS$lichens$south)[[2]])

CFsm <- array(0, c(length(cfs), nrow(cfsm), 100))
dimnames(CFsm) <- list(names(cfs), rownames(cfsm), NULL)
CFsj <- array(0, c(length(cfs), nrow(cfsj), 100))
dimnames(CFsj) <- list(names(cfs), rownames(cfsj), NULL)

for (spp in names(cfs)) {
    B <- min(100, nrow(cfs[[spp]]$estSpace))
    cfsm <- rbind(log(cfs[[spp]]$coefARU), pAspen=cfs[[spp]]$estARU[1:B,"pAspen"])
    cfsj <- rbind(log(cfs[[spp]]$coefSpace), t(cfs[[spp]]$estSpace[1:B,c("pAspen", cnb)]))
    cfsm[cfsm > 10^4] <- 10^4
    cfsm[cfsm < -10^4] <- -10^4
    cfsj[cfsj > 10^4] <- 10^4
    cfsj[cfsj < -10^4] <- -10^4
    if (ncol(cfsm) < 100) {
        b <- sample.int(ncol(cfsm), 100-ncol(cfsm), replace=TRUE)
        cfsm <- cbind(cfsm, cfsm[,b])
        cfsj <- cbind(cfsj, cfsj[,b])
    }
    CFsm[spp,,] <- cfsm
    CFsj[spp,,] <- cfsj
}

# north
spp <- "ALFL"
B <- min(100, nrow(cfn[[spp]]$estSpace))
cfnm <- rbind(log(cfn[[spp]]$coefARU))
cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[1:B,cnb]))
cfnm[cfnm > 10^4] <- 10^4
cfnm[cfnm < -10^4] <- -10^4
cfnj[cfnj > 10^4] <- 10^4
cfnj[cfnj < -10^4] <- -10^4
if (ncol(cfnm) < 100) {
    b <- sample.int(ncol(cfnm), 100-ncol(cfnm), replace=TRUE)
    cfnm <- cbind(cfnm, cfnm[,b])
    cfnj <- cbind(cfnj, cfnj[,b])
}
compare_sets(rownames(cfnm),dimnames(COEFS$lichens$north)[[2]])

CFnm <- array(0, c(length(cfn), nrow(cfnm), 100))
dimnames(CFnm) <- list(names(cfn), rownames(cfnm), NULL)
CFnj <- array(0, c(length(cfn), nrow(cfnj), 100))
dimnames(CFnj) <- list(names(cfn), rownames(cfnj), NULL)

for (spp in names(cfn)) {
    B <- min(100, nrow(cfn[[spp]]$estSpace))
    cfnm <- rbind(log(cfn[[spp]]$coefARU))
    cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[1:B,cnb]))
    cfnm[cfnm > 10^4] <- 10^4
    cfnm[cfnm < -10^4] <- -10^4
    cfnj[cfnj > 10^4] <- 10^4
    cfnj[cfnj < -10^4] <- -10^4
    if (ncol(cfnm) < 100) {
        b <- sample.int(ncol(cfnm), 100-ncol(cfnm), replace=TRUE)
        cfnm <- cbind(cfnm, cfnm[,b])
        cfnj <- cbind(cfnj, cfnj[,b])
    }
    CFnm[spp,,] <- cfnm
    CFnj[spp,,] <- cfnj
}

COEFS$birds <- list(
    north=list(marginal=CFnm, joint=CFnj),
    south=list(marginal=CFsm, joint=CFsj))

save(COEFS,
    file="s:/AB_data_v2020/Results/COEFS-ALL.RData")


## mapping and sector effects

load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData")


ksc <- dd_2018$soil_current
ksc <- ksc / ifelse(rowSums(ksc) > 0, rowSums(ksc), 1)

cns <- rownames(COEFS$mites$south[1,,])
cns0 <- cns[1:(which(cns=="Intercept")-2)]
cns1 <- cns[which(cns=="Intercept"):length(cns)]
cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]
cnn1 <- cnn[which(cnn=="Intercept"):length(cnn)]
all(cns1==cnn1)

LTs <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(LTs) <- LTs[,1]
LTs$Map <- as.character(LTs$UseInAnalysis)
LTs$Map[LTs$Map == "Well"] <- "Wellsites"
LTs$Map[LTs$Map %in% c("SoilWater", "HWater")] <- "Water"
LTs$Map[LTs$Map %in% c("SoilUnknown", "HFor")] <- "UNK"
LTs <- LTs[,c("Map", "Sector")]


make_clim <- function(x, birds=FALSE) {
    if (birds) {
        z <- with(x, cbind(
            pWater_KM=pWater,
            pWater2_KM=pWater^2,
            xPET=(PET - 0) / 800,
            xMAT=(MAT - 0) / 6,
            xAHM=(AHM - 0) / 50,
            xFFP=(FFP - 0) / 130,
            xMAP=(MAP - 0) / 2300,
            xMWMT=(MWMT - 0) / 20,
            xMCMT=(MCMT - 0) / 25,
            xY=(POINT_Y - 54.1) / 2.43,
            xX=(POINT_X - (-113.3)) / 2.12))
        z <- cbind(z,
            xY2=z[,"xY"]^2,
            xX2=z[,"xX"]^2,
            `xFFP:xMAP`=z[,"xFFP"]*z[,"xMAP"],
            `xMAP:xPET`=z[,"xMAP"]*z[,"xPET"],
            `xAHM:xMAT`=z[,"xAHM"]*z[,"xMAT"],
            `xX:xY`=z[,"xX"]*z[,"xY"])
    } else {
        z <- with(x, cbind(
            Intercept=1,
            Lat=x$POINT_Y,
            Long=x$POINT_X,
            AHM=x$AHM,
            PET=x$PET,
            FFP=x$FFP,
            MAP=x$MAP,
            MAT=x$MAT,
            MCMT=x$MCMT,
            MWMT=x$MWMT,
            Lat2=x$POINT_Y^2,
            Long2=x$POINT_X^2,
            LatLong=x$POINT_X*x$POINT_Y,
            MAPPET=x$MAP*x$PET,
            MATAHM=x$MAT*x$AHM,
            MAPFFP=x$MAP*x$FFP,
            MAT2=x$MAT^2,
            MWMT2=x$MWMT^2))
    }
    z
}

## south - non birds

lt <- list(
    south=structure(list(Map = c("ClaySub", "Other", "Other", "Other",
        "Other", "Other", "RapidDrain", "Loamy", "Other", "Other", "Other",
        "Other", "Other", "Other", "ClaySub", "SandyLoam", "RapidDrain",
        "RapidDrain", "RapidDrain", "RapidDrain", "RapidDrain", "ThinBreak",
        "Blowout", "Other", "Other", "Blowout", "UNK", "Water", "Urban",
        "Urban", "Rural", "Industrial", "Industrial", "Rural", "Mine",
        "Mine", "Wellsites", "EnSoftLin", "EnSoftLin", "EnSeismic", "EnSeismic",
        "HardLin", "HardLin", "TrSoftLin", "TrSoftLin", "TrSoftLin",
        "Crop", "RoughP", "RoughP", "TameP", "Industrial", "Water", "Water",
        "Water", "Water", "UNK"), Sector = c("Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "RuralUrban", "RuralUrban",
        "RuralUrban", "RuralUrban", "Energy", "Misc", "Energy", "Misc",
        "Energy", "Energy", "Energy", "Energy", "Energy", "Transportation",
        "Transportation", "Transportation", "Transportation", "Transportation",
        "Agriculture", "Agriculture", "Agriculture", "Agriculture", "Agriculture",
        "Misc", "Misc", "Misc", "Misc", "Forestry")), row.names = c("Cy",
        "Len", "LenS", "LenSP", "LenT", "LenW", "Li", "Lo", "Ltc", "LtcC",
        "LtcD", "LtcH", "LtcS", "LtcR", "Sb", "Sy", "BdL", "CS", "Gr",
        "Sa", "SwG", "TB", "BlO", "LenA", "Ov", "SL", "UNK", "Water",
        "UrbanIndustrial", "UrbanResidence", "RuralResidentialIndustrial",
        "IndustrialSiteRural", "WindGenerationFacility", "OtherDisturbedVegetation",
        "MineSite", "PeatMine", "WellSite", "Pipeline", "TransmissionLine",
        "SeismicLineNarrow", "SeismicLineWide", "RoadHardSurface", "RailHardSurface",
        "RoadTrailVegetated", "RoadVegetatedVerge", "RailVegetatedVerge",
        "CultivationCrop", "CultivationAbandoned", "CultivationRoughPasture",
        "CultivationTamePasture", "HighDensityLivestockOperation", "BorrowpitsDugoutsSumps",
        "MunicipalWaterSewage", "Reservoirs", "Canals", "CutBlocks"), class = "data.frame"),
    north=list())


## -- this is not for sector effects!
##
## cfs: coefficients (soil type x boot)
## cl: climate variables (cell x var)
## pa: pAspen values (vector of length #cell)
## ps: soil class areas or proportions
## i: bootstrap iteration (1 by default)
## returns a #cell x sector matrix
pred_south_nonbird <- function(cfs, cl, pa, ps, i=1) {
    CL <- make_clim(x)
    LTs <- lt$south
    sclass <- rownames(LTs)[!(LTs$Map %in% c("UNK", "Water"))]
    sclass <- intersect(colnames(ps), sclass)
    LTsx <- LTs[sclass,]
    w <- groupSums(ps[,sclass], 2, LTs[sclass,"Map"])
    rs <- rowSums(w)
    w <- w / ifelse(rs > 0, rs, 1)

    mucl <- drop(CL %*% cfs[colnames(CL),i])
    mupa <- pa * cfs["pAspen",i]
    mu <- as.matrix(w)
    mu[] <- mucl + mupa
    mu <- t(t(mu) + cfs[colnames(mu),i])
    rabunds <- rowSums(w * plogis(mu))



}

ss <- which(kgrid$NRNAME=="Grassland")[1:100]
i <- 1 # bootstrap run
## spp coef table with boot
cfs <- COEFS$vplants$south[1,,1:10]
## climate info
cl <- kgrid[ss,]
## soil class proportions
ps <- ksc[ss,]
## pAspen
pa <- kgrid[ss,"pAspen"]




## plots


## south
.get_stats <- function(x) {
    q <- apply(x, 1, quantile, c(0.5, 0.05, 0.95))
    data.frame(Label=factor(rownames(x), rownames(x)),
        First=x[,1], Mean=rowMeans(x), SD=apply(x, 1, sd),
        Median=q[1,], Lower=q[2,], Upper=q[3,])
}
cns <- c("Loamy", "SandyLoam", "RapidDrain", "ClaySub", "ThinBreak", "Blowout", "Other",
    "Crop", "TameP", "RoughP",
    "EnSeismic", "EnSoftLin", "TrSoftLin", "HardLin",
    "Wellsites", "Rural", "Urban", "Industrial",
    "Mine", "MineV", "Water")
cnn <- c(
    "WhiteSpruceR", "WhiteSpruce1", "WhiteSpruce2", "WhiteSpruce3",
    "WhiteSpruce4", "WhiteSpruce5", "WhiteSpruce6", "WhiteSpruce7",
    "WhiteSpruce8", "PineR", "Pine1", "Pine2", "Pine3", "Pine4",
    "Pine5", "Pine6", "Pine7", "Pine8", "DeciduousR", "Deciduous1",
    "Deciduous2", "Deciduous3", "Deciduous4", "Deciduous5", "Deciduous6",
    "Deciduous7", "Deciduous8", "MixedwoodR", "Mixedwood1", "Mixedwood2",
    "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7",
    "Mixedwood8", "TreedBogR", "TreedBog1", "TreedBog2", "TreedBog3",
    "TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7", "TreedBog8",
    "TreedSwamp", "ShrubbySwamp", "ShrubbyBog", "ShrubbyFen", "GraminoidFen",
    "Marsh", "Shrub", "CCWhiteSpruceR", "CCWhiteSpruce1", "CCWhiteSpruce2",
    "CCWhiteSpruce3", "CCWhiteSpruce4", "CCPineR", "CCPine1", "CCPine2",
    "CCPine3", "CCPine4", "CCDeciduousR", "CCDeciduous1", "CCDeciduous2",
    "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR", "CCMixedwood1",
    "CCMixedwood2", "CCMixedwood3", "CCMixedwood4", "Crop", "TameP",
    "RoughP", "Wellsites", "EnSeismic", "EnSoftLin", "TrSoftLin",
    "HardLin", "TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3",
    "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8",
    "Rural", "Urban", "Industrial", "Mine", "MineV", "Water", "GrassHerb",
    "Intercept", "Lat", "Long", "AHM", "PET", "FFP", "MAP", "MAT",
    "MCMT", "MWMT", "Lat2", "Long2", "LatLong", "MAPPET", "MATAHM",
    "MAPFFP", "MAT2", "MWMT2")

taxon <- "vplants"
spp <- "Achillea.borealis"
cf <- COEFS[[taxon]]$south[spp,,]


FUN <- plogis

ra <- rbind(data.frame(pAspen="Non-treed", get_stats(FUN(cf[cns,]))),
    data.frame(pAspen="Treed", get_stats(FUN(t(t(cf[cns,]) + cf["pAspen",])))))

library(ggplot2)

p <- ggplot(ra, aes(x=Label, y=First, fill=Label, group=pAspen)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  width=0.2, color="#105A73") +
    ylab("Relative abundance") +
    xlab("Land cover class") +
    facet_wrap(~pAspen) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_minimal() +
    labs(title=paste0(spp, " (", taxon, ")")) +
    theme(legend.position="none")



labs <- list(
    veg=list(
        native=c(`White spruce 0-9` = "WhiteSpruceR", `White spruce 10-19` = "WhiteSpruce1",
            `White spruce 20-39` = "WhiteSpruce2", `White spruce 40-59` = "WhiteSpruce3",
            `White spruce 60-79` = "WhiteSpruce4", `White spruce 80-99` = "WhiteSpruce5",
            `White spruce 100-119` = "WhiteSpruce6", `White spruce 120-139` = "WhiteSpruce7",
            `White spruce 140+` = "WhiteSpruce8", `Pine 0-9` = "PineR", `Pine 10-19` = "Pine1",
            `Pine 20-39` = "Pine2", `Pine 40-59` = "Pine3", `Pine 60-79` = "Pine4",
            `Pine 80-99` = "Pine5", `Pine 100-119` = "Pine6", `Pine 120-139` = "Pine7",
            `Pine 140+` = "Pine8", `Deciduous 0-9` = "DeciduousR", `Deciduous 10-19` = "Deciduous1",
            `Deciduous 20-39` = "Deciduous2", `Deciduous 40-59` = "Deciduous3",
            `Deciduous 60-79` = "Deciduous4", `Deciduous 80-99` = "Deciduous5",
            `Deciduous 100-119` = "Deciduous6", `Deciduous 120-139` = "Deciduous7",
            `Deciduous 140+` = "Deciduous8", `Mixedwood 0-9` = "MixedwoodR",
            `Mixedwood 10-19` = "Mixedwood1", `Mixedwood 20-39` = "Mixedwood2",
            `Mixedwood 40-59` = "Mixedwood3", `Mixedwood 60-79` = "Mixedwood4",
            `Mixedwood 80-99` = "Mixedwood5", `Mixedwood 100-119` = "Mixedwood6",
            `Mixedwood 120-139` = "Mixedwood7", `Mixedwood 140+` = "Mixedwood8",
            `Treed bog 0-9` = "TreedBogR", `Treed bog 10-19` = "TreedBog1",
            `Treed bog 20-39` = "TreedBog2", `Treed bog 40-59` = "TreedBog3",
            `Treed bog 60-79` = "TreedBog4", `Treed bog 80-99` = "TreedBog5",
            `Treed bog 100-119` = "TreedBog6", `Treed bog 120-139` = "TreedBog7",
            `Treed bog 140+` = "TreedBog8",
            `Treed fen 0-9`="TreedFenR", `Treed fen 10-19`="TreedFen1", `Treed fen 20-39`="TreedFen2",
            `Treed fen 40-59`="TreedFen3",`Treed fen 60-79`= "TreedFen4", `Treed fen 80-99`="TreedFen5",
            `Treed fen 100-119`="TreedFen6", `Treed fen 120-139`="TreedFen7", `Treed fen 140+`="TreedFen8",
            `Treed swamp` = "TreedSwamp", `Shrubby swamp` = "ShrubbySwamp",
            `Shrubby bog` = "ShrubbyBog", `Shrubby fen` = "ShrubbyFen", `Graminoid fen` = "GraminoidFen",
            `Marsh` = "Marsh", `Shrub` = "Shrub", `Grass/herb` = "GrassHerb"),
        hffor=c(`Harvest, White spruce 0-9` = "CCWhiteSpruceR",
            `Harvest, white spruce 10-19` = "CCWhiteSpruce1",
            `Harvest, white spruce 20-39` = "CCWhiteSpruce2",
            `Harvest, white spruce 40-59` = "CCWhiteSpruce3",
            `Harvest, white spruce 60-79` = "CCWhiteSpruce4",
            `Harvest, pine 0-9` = "CCPineR", `Harvest, pine 10-19` = "CCPine1",
            `Harvest, pine 20-39` = "CCPine2", `Harvest, pine 40-59` = "CCPine3",
            `Harvest, pine 60-79` = "CCPine4",
            `Harvest, deciduous 0-9` = "CCDeciduousR", `Harvest, deciduous 10-19` = "CCDeciduous1",
            `Harvest, deciduous 20-39` = "CCDeciduous2", `Harvest, deciduous 40-59` = "CCDeciduous3",
            `Harvest, deciduous 60-79` = "CCDeciduous4",
            `Harvest, mixedwood 0-9` = "CCMixedwoodR",
            `Harvest, mixedwood 10-19` = "CCMixedwood1", `Harvest, mixedwood 20-39` = "CCMixedwood2",
            `Harvest, mixedwood 40-59` = "CCMixedwood3", `Harvest, mixedwood 60-79` = "CCMixedwood4"),
        hfpoly=c(
            "Cropland"="Crop",
            "Tame pasture"="TameP",
            "Rough pasture"="RoughP",
            "Well sites"="Wellsites",
            "Rural residential"="Rural",
            "Urban/Industrial"="Urban",
            "Industrial (rural)"="Industrial"),
        hfsoft=c(
            "Seismic lines"="EnSeismic",
            "Vegetated linear (energy)"="EnSoftLin",
            "Vegetated linear (transport)"="TrSoftLin"),
        hfhard=c(
            "Non-vegetated linear"="HardLin")
    ),
    soil=list(
        native=c(
            "Loamy"="Loamy",
            "Sandy/loamy"="SandyLoam",
            "Rapid drain"="RapidDrain",
            "Clay"="ClaySub",
            "Thin break"="ThinBreak",
            "Blowout"="Blowout",
            "Other soil types"="Other"),
        hfpoly=c(
            "Cropland"="Crop",
            "Tame pasture"="TameP",
            "Rough pasture"="RoughP",
            "Well sites"="Wellsites",
            "Rural residential"="Rural",
            "Urban/Industrial"="Urban",
            "Industrial (rural)"="Industrial"),
        hfsoft=c(
            "Seismic lines"="EnSeismic",
            "Vegetated linear (energy)"="EnSoftLin",
            "Vegetated linear (transport)"="TrSoftLin"),
        hfhard=c(
            "Non-vegetated linear"="HardLin"),
        extras=c(
            "Mine (non-vegetated)"="Mine",
            "Mine (vegetated)"="MineV",
            "Open water"="Water"))
    )

.plot_abundance_soil_2020 <- function(cf, plot=TRUE, paspen=0, ylim, main, ylab,
    bird=FALSE, ...) {
    if (paspen %)(% c(0, 1))
        stop("paspen must be in [0, 1]")
    lab <- c(labs$soil$native, labs$soil$hfpoly)

    FUN <- if (bird)
        exp else plogis
    p <- t(t(cf[lab,]) + paspen*cf["pAspen",])
    p <- .get_stats(FUN(p))

    out <- p[,c("Median", "Lower", "Upper")]
    colnames(out) <- c("Estimate", "LCL", "UCL")
    rownames(out) <- names(lab)
    out[out < 10^-5] <- 10^-6
    if (plot) {
        op <- par(mai=c(1.5, 1, 0.2, 0.3))
        on.exit(par(op))

        lci <- out$LCL
        uci <- out$UCL
        y1 <- out$Estimate
        x <- 1:14

        if (missing(main))
            main <- ""
        if (missing(ylab))
            ylab <- "Relative abundance"
        if (missing(ylim)) {
            ymax <- max(min(max(uci[x]), 2*max(y1)), y1*1.02)
            ylim <- c(0, ymax)
        } else {
            ymax <- max(ylim)
        }

        space <- c(1, x[-1] -x[-length(x)]) - 0.9
        col <- c(
            hcl.colors(7,"Dark 3"),
            "#333300", "#666600", "#999900",
            "#909090","#555555","#404040","#222222")

        x1 <- barplot(y1[x], space=space, border="white", col=col,
            ylim=ylim, xlim=c(-0.5,max(x)+1.5), xaxs="i", yaxt="n", ylab=ylab,
            col.lab="grey50", cex.lab=1.2,axisnames=FALSE)[,1]
        ax <- axis(side=2, cex.axis=0.9, col.axis="grey50", col.ticks="grey50", las=2)
        abline(h=ax, col="grey80")
        x1 <- barplot(y1[x], space=space, border="white", col=col,
            xaxs="i", yaxt="n", ylab=ylab,
            col.lab="grey50", cex.lab=1.2,axisnames=FALSE, add=TRUE)[,1]
        box(bty="l", col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[i], y1[i]), col="grey90")
            lines(rep(x1[i],2), c(uci[i], y1[i]), col=col[i])
        }
        mtext(side=1, at=x1, line=0.7,
            names(lab),
            col=col, las=2)
        mtext(side=3, at=0, adj=0, main, col="grey30")
    }
    invisible(out)
}


## input is cf: LC x Boot matrix for a single species
.plot_abundance_lin_2020 <- function(cf, plot=TRUE, veg=TRUE, ylim, main, xlab, ylab,
    bird=FALSE, ...) {
    FUN <- if (bird)
        exp else plogis
    if (veg) {
        x0 <- FUN(cf[labs$veg$native,,drop=FALSE])
        xs <- FUN(cf[labs$veg$hfsoft,,drop=FALSE])
        xh <- FUN(cf[labs$veg$hfhard,,drop=FALSE])
    } else {
        x0 <- FUN(cf[labs$soil$native,,drop=FALSE])
        xs <- FUN(cf[labs$soil$hfsoft,,drop=FALSE])
        xh <- FUN(cf[labs$soil$hfhard,,drop=FALSE])
    }
    # apply land cover weights (proportions here if needed)
    tab <- rbind(AverageCoef=colMeans(x0),
        SoftLin10=0.9*colMeans(x0)+0.1*colMeans(xs),
        HardLin10=0.9*colMeans(x0)+0.1*colMeans(xh))
    out <- .get_stats(tab)

    if (plot) {
        p.mean <- out["AverageCoef", "Median"]
        p.softlin10 <- out["SoftLin10", "Median"]
        p.hardlin10 <- out["HardLin10", "Median"]

        if (missing(ylim)) {
            ymax1 <- max(p.softlin10, p.hardlin10, 2*p.mean)*1.03
            ylim <- c(0, ymax1)
        } else {
            ymax <- max(ylim)
        }
#        if (missing(main))
#            main <- species
        if (missing(main))
            main <- ""
        if (missing(xlab))
            xlab <- "Human footprint"
        if (missing(ylab))
            ylab <- "Relative abundance"

        plot(c(1,1.95,2.05), c(p.mean, p.softlin10, p.hardlin10),
            pch=c(1, 16, 15),
            col=c("grey30", "blue3", "red4"),
            xlab=xlab, ylab=ylab, xlim=c(0.8, 2.8), ylim=ylim,
            tck=0.01, yaxs="i", xaxt="n", yaxt="n", bty="l",
            cex=2, lwd=2, cex.lab=1.4, cex.axis=1.3, col.lab="grey40")
        axis(side=2, at=pretty(ylim, n=5), cex.axis=1.3, tck=0.01,
            cex.axis=1.3, col.axis="grey40", col.ticks="grey40")
        axis(side=1, at=c(1,2), labels=c("None","10% linear"),
            tck=0.01, cex.axis=1.3, col.axis="grey40", col.ticks="grey40")
        box(bty="l", col="grey40")
        lines(c(1,1.95), c(p.mean, p.softlin10), col="blue3")
        lines(c(1,2.05), c(p.mean, p.hardlin10), col="red4")
        points(c(1, 1.95, 2.05), c(p.mean, p.softlin10, p.hardlin10),
            pch=c(1,16,15), col=c("grey30", "blue3", "red4"), cex=2, lwd=2)
        lines(c(1,1), c(out["AverageCoef", "Lower"], out["AverageCoef", "Upper"]), col="grey30")
        lines(c(1.95,1.95), c(out["SoftLin10", "Lower"], out["SoftLin10", "Upper"]), col="blue3")
        lines(c(2.05,2.05), c(out["HardLin10", "Lower"], out["HardLin10", "Upper"]), col="red4")
        ly <- c(p.softlin10, p.hardlin10)
        if (abs(ly[2]-ly[1]) < ymax1/20)
            ly <- c(mean(ly)+ymax1/40*sign(ly[1]-ly[2]), mean(ly)+ymax1/40*sign(ly[2]-ly[1]))
        text(c(2.15, 2.15), ly, c("Soft linear","Hard linear"),
            col=c("blue3", "red4"), cex=1.3, adj=0)
        mtext(side=3, at=0.8, adj=0, main, col="grey30", cex=1.3)
    }
    invisible(out)
}



.plot_abundance_veg_2020 <- function(cf, plot=TRUE, ylim, main, ylab, bw=FALSE,
    bird=FALSE, ...){

    lab <- c(labs$veg$native, labs$veg$hfpoly, labs$veg$hffor)
    # collapse treed fen as average

    FUN <- if (bird)
        exp else plogis
    p <- cf[lab,]
    ii <- grep("TreedFen", lab)
    p[ii[1],] <- colMeans(p[ii,])
    p <- p[-ii[-1],]
    lab[ii[1]] <- "TreedFen"
    names(lab)[ii[1]] <- "Treed fen"
    lab <- lab[-ii[-1]]
    rownames(p) <- lab
    p <- .get_stats(FUN(p))

    out <- p[,c("Median", "Lower", "Upper")]
    colnames(out) <- c("Estimate", "LCL", "UCL")
    rownames(out) <- names(lab)
    out[out < 10^-5] <- 10^-6

    if (plot) {

        lci <- out$LCL
        uci <- out$UCL
        y1 <- out$Estimate
        names(y1) <- rownames(out)

        if (missing(main))
            main <- ""
        if (missing(ylab))
            ylab <- "Relative abundance"
        if (missing(ylim)) {
            ymax <- min(max(uci), 2 * max(y1))
            ylim <- c(0, ymax)
        } else {
            ymax <- max(ylim)
        }

        x <- c(rep(1:9, 5) + rep(seq(0, 40, 10), each=9),
            51, 53, 55, 57, 59, 61, 63, 65, 67,
            70, 72, 74, 76, 78, 80, 82)
        space <- c(1,x[-1]-x[-length(x)])-0.99

        op <- par(mai=c(1.5,1,0.2,0.3))
        on.exit(par(op))

        if (!bw) {
            col.r <- c(rep(0,9), seq(0.3,0.6,length.out=9),
                seq(0.5,1,length.out=9),
                seq(0.8,0.9,length.out=9), rep(0,9), 0, #rep(0,9),
                0.8,0.2,0,0,0,0,0,0,
                0.14,0.29,0.6,0.4,0,0,0)
            col.g <- c(seq(0.5,1,length.out=9), seq(0.4,0.8,length.out=9),
                seq(0.1,0.2,length.out=9), seq(0.4,0.8,length.out=9),
                seq(0.4,0.7,length.out=9), 0.5, #seq(0.15,0.5,length.out=9),
                0.8,0.8,0,0,0,0,0,0,
                0.14,0.29,0.6,0.40,0,0)
            col.b <- c(rep(0,9),rep(0,9),rep(0,9),seq(0.2,0.4,length.out=9),
                seq(0.2,0.6,length.out=9), 0.7, #seq(0.4,0.7,length.out=9),
                0,0,1,1,0,0,0,0,
                0,0,0,0.4,0,0,0)
        } else {
            col.r <- c(rep(seq(0.7,0.2,length.out=9), 5), rep(0.3,15))
            col.b <- col.g <- col.r
        }
        COL <- rgb(col.r, col.g, col.b)
        if (!bw) {
            COL[47:61] <- c(hcl.colors(8,"Dark 3"),
                "#333300", "#666600", "#999900",
                "#909090","#555555","#404040","#222222")

        }
        idx <- 1:length(x)
        x1 <- barplot(y1[idx],
            space=space,
            border="white",
            col=COL,
            ylim=ylim,
            xlim=c(-0.5,max(x)+1.5),
            xaxs="i", yaxt="n",
            ylab=ylab,
            col.lab="grey50",
            cex.lab=1.2, axisnames=FALSE)[,1]
        ax <- axis(side=2, cex.axis=0.9, col.axis="grey50",
            col.ticks="grey50", las=2)
        abline(h=ax, col="grey80")
        x1 <- barplot(y1[idx],
            space=space,
            border="white", col=COL,
            ylim=ylim, xaxs="i", yaxt="n",
            col.lab="grey50",
            cex.lab=1.2, axisnames=FALSE, add=TRUE)[,1]
        box(bty="l", col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[idx][i], y1[idx][i]),col="grey90")
            lines(rep(x1[i],2), c(uci[idx][i], y1[idx][i]),col=COL[i])
        }
        mtext(side=1, at=x1[c(5,14,23,32,41)], line=1.4,
            c("Upland Spruce","Pine","Deciduous","Mixedwood","Black Spruce"),
            col=rgb(col.r[c(5,14,23,32,41)], col.g[c(5,14,23,32,41)],
            col.b[c(5,14,23,32,41)]),las=1)
        at1<-rep(seq(1,9,2),5) + rep(c(0,9,18,27,36), each=5)
        mtext(side=1,at=x1[at1]-0.3,rep(c("0","20","60","100","140"),5),
            line=0.2, adj=0.5, cex=0.8, col=COL[at1])
        mtext(side=1, at=-0.25, adj=1, line=0.2, "Age:", col="grey40", cex=0.8)
        mtext(side=3, at=0, adj=0, main, col="grey30")

        ii <- match(c("TreedFen", "TreedSwamp", "ShrubbySwamp", "ShrubbyBog", "ShrubbyFen",
            "GraminoidFen", "Marsh", "Shrub", "GrassHerb"), lab)
        mtext(side=1, at=x1[ii],
            names(lab)[ii],
            col=COL[ii],
            las=2, adj=1.1)
        ii <- match(c("Crop", "TameP", "RoughP", "Wellsites", "Rural", "Urban", "Industrial"), lab)
        mtext(side=1, at=x1[ii],
            names(lab)[ii],
            col=COL[ii], las=2, adj=1.1)

        ## Add cutblock trajectories - upland conifer
        i1 <- grep("CCWhiteSpruce", lab)
        x2<-x1[1:5]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="White pruce 80-99")
        lines(c(x2[1:5], x1[x3]), y1[c(i1, x3)], col="grey30", lty=2)
        points(x2[1:5], y1[i1], pch=18, cex=1, col="grey30")
        points(x2[1:5], y1[i1], pch=5, cex=0.7, col="grey10")
        ## Pine
        i1 <- grep("CCPine", lab)
        x2<-x1[10:15]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Pine 80-99")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Deciduous
        i1 <- grep("CCDeciduous", lab)
        x2<-x1[19:24]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Deciduous 80-99")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Mixed
        i1 <- grep("CCMixedwood", lab)
        x2<-x1[28:33]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Mixedwood 80-99")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
    }
    invisible(out)
}


## etc --

CFn <- array(0, c(length(SPPn), dim(COEFS[[1]]$north)[2], B))
dimnames(CFn) <- list(
    names(SPPn),
    dimnames(COEFS[[1]]$north)[[2]],
    NULL)
CFs <- array(0, c(length(SPPs), dim(COEFS[[1]]$south)[2], B))
dimnames(CFs) <- list(
    names(SPPs),
    dimnames(COEFS[[1]]$south)[[2]],
    NULL)

spp <- "ALFL"
resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
estn1 <- suppressWarnings(get_coef(resn, Xn, stage="ARU", na.out=FALSE))
ests1 <- suppressWarnings(get_coef(ress, Xs, stage="ARU", na.out=FALSE))
estn2 <- suppressWarnings(get_coef(resn, Xn, stage="Space", na.out=FALSE))
ests2 <- suppressWarnings(get_coef(ress, Xs, stage="Space", na.out=FALSE))

cf1n <- get_coef_north(resn, STAGE="ARU", subset=1)
cf2n <- get_coef_north(resn, STAGE="Space", subset=1)
cf1s <- get_coef_south(ress, STAGE="ARU", subset=1)
cf2s <- get_coef_south(ress, STAGE="Space", subset=1)


cnn0 <- dimnames(COEFS$vplants$north)[[2]]
cns0 <- dimnames(COEFS$vplants$south)[[2]]
cnn <- cnn0[1:(which(cnn0 == "Intercept")-1L)]
cns <- cns0[1:(which(cns0 == "pAspen")-1L)]

compare_sets(cnn, names(cf1n))
intersect(cnn, names(cf1n))
setdiff(cnn, names(cf1n))
setdiff(names(cf1n), cnn)

compare_sets(cns, names(cf1s))
intersect(cns, names(cf1s))
setdiff(cns, names(cf1s))
setdiff(names(cf1s), cns)

