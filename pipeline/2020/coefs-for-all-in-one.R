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

    cfn <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
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

    cfs <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}

save(cfs, file="s:/AB_data_v2020/Results/BIRDS-South.RData")


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


## plots

load("s:/AB_data_v2020/Results/COEFS-EA.RData")

## south
get_stats <- function(x) {
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

