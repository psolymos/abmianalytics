## organizing results for mites, lichens, mosses, vplants

ROOT <- "s:/AB_data_v2020/Results/Results from Ermias/Preliminary/"

TAXA <- c("lichens", "mites", "mosses", "vplants")
#taxon <- "mites"

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

## bootstrap based coef tables (i.e. not normal SEs)

ROOT <- "s:/AB_data_v2020/Results/Results from Ermias/Boot coef/"

TAXA <- c("lichens", "mites", "mosses", "vplants")
#taxon <- "mites"

COEFS <- list()
B <- 250

for (taxon in TAXA) {

    if (taxon == "lichens") {
        fn <- "Lichen_North Bootstrap coefficents 2020.Rdata"
        fs <- "Lichen_South Bootstrap coefficents 2020.Rdata"
    }
    if (taxon == "mites") {
        fn <- "Mite_North Bootstrap coefficents 2020.Rdata"
        fs <- "Mite_South Bootstrap coefficents 2020.Rdata"
    }
    if (taxon == "mosses") {
        fn <- "Moss_North Bootstrap coefficents 2020.Rdata"
        fs <- "Moss_South Bootstrap coefficents 2020.Rdata"
    }
    if (taxon == "vplants") {
        fn <- "VPlant_North Bootstrap coefficents 2020.Rdata"
        fs <- "VPlant_South Bootstrap coefficents 2020.Rdata"
    }

    #' load North and south coefs into their respective environments
    en <- new.env()
    load(paste0(ROOT, fn), envir=en)
    es <- new.env()
    load(paste0(ROOT, fs), envir=es)

    #' # North
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    i <- 1
    cfn <- en$Coef.bs[,,i]
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
    xCFn <- cbind(qlogis(cfn), en$Sc.coef.bs[,,i])
    #' array: bootstrap x species x coefs
    CFn <- array(0, c(nrow(xCFn), ncol(xCFn), B),
        dimnames=list(rownames(xCFn), colnames(xCFn), NULL))
    for (i in 1:B) {
        cfn <- en$Coef.bs[,,i]
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
        CFn[,,i] <- cbind(qlogis(cfn), en$Sc.coef.bs[,,i])
    }


    #' # South
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    i <- 1
    cfs <- es$Coef.bs[,,i]
    cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
    cfs <- cfs[,colnames(cfs) != "UrbInd"]
    xCFs <- cbind(qlogis(cfs), pAspen=es$Coef.pAspen.bs[[i]], es$Sc.coef.bs[,,i])
    #' array: bootstrap x species x coefs
    CFs <- array(0, c(nrow(xCFs), ncol(xCFs), B),
        dimnames=list(rownames(xCFs), colnames(xCFs), NULL))
    for (i in 1:B) {
        cfs <- es$Coef.bs[,,i]
        cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        cfs <- cfs[,colnames(cfs) != "UrbInd"]
        CFs[,,i] <- cbind(qlogis(cfs), pAspen=es$Coef.pAspen.bs[[i]], es$Sc.coef.bs[,,i])
    }
    CFn[CFn < -10^4] <- -10^4
    CFn[CFn > 10^4] <- 10^4
    CFs[CFs < -10^4] <- -10^4
    CFs[CFs > 10^4] <- 10^4

    COEFS[[taxon]] <- list(north=CFn, south=CFs)

}

save(COEFS, file="s:/AB_data_v2020/Results/COEFS-EAboot.RData")


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

load("s:/AB_data_v2020/Results/COEFS-EAboot.RData")

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

## combine birds and other taxa

library(mefa4)

## species names etc
ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds" # change this bit
ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)
TAX <- ee$tax
rm(ee)

#load("s:/AB_data_v2020/Results/COEFS-EA.RData")
load("s:/AB_data_v2020/Results/COEFS-EAboot.RData")
load("s:/AB_data_v2020/Results/BIRDS-North.RData")
load("s:/AB_data_v2020/Results/BIRDS-South.RData")

if (FALSE) { # testing CC effects
spp <- "OSFL"
h <- c("CCWhiteSpruceR", "CCWhiteSpruce1",
"CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "CCPineR",
"CCPine1", "CCPine2", "CCPine3", "CCPine4", "CCDeciduousR", "CCDeciduous1",
"CCDeciduous2", "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR",
"CCMixedwood1", "CCMixedwood2", "CCMixedwood3", "CCMixedwood4")
h1 <- gsub("CC", "", h)
m  <- matrix(rowMeans(cfn[[spp]]$coefARU)[h], 5, 4)
m1 <- matrix(rowMeans(cfn[[spp]]$coefARU)[h1], 5, 4)
plot(0, type="n", xlim=c(0,4), ylim=c(0,max(m,m1)))
for (i in 1:4) {
    lines(0:4, m[,i], col=i, lty=2)
    lines(0:4, m1[,i], col=i, lty=1)
}
}

SPPn <- as.character(TAX$sppid[match(names(cfn), TAX$code)])
SPPs <- as.character(TAX$sppid[match(names(cfs), TAX$code)])
rownames(TAX) <- TAX$sppid
TAX <- TAX[sort(union(SPPn,SPPs)),]
names(cfs) <- SPPs
names(cfn) <- SPPn

## organizing bird coefs

cns <- rownames(COEFS$mites$south[1,,])
cns0 <- cns[1:(which(cns=="Intercept")-2)]
cns1 <- cns[(which(cns=="Intercept")-1):length(cns)]
cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]
cnn1 <- cnn[which(cnn=="Intercept"):length(cnn)]

cnb <- c("pWater_KM", "pWater2_KM", "xPET",
    "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT", "xMCMT", "xY", "xX",
    "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY")

BMAX <- 250

# south
spp <- 1
B <- min(BMAX, nrow(cfs[[spp]]$estSpace))
cfsm <- rbind(log(cfs[[spp]]$coefARU), pAspen=cfs[[spp]]$estARU[1:B,"pAspen"])
cfsj <- rbind(log(cfs[[spp]]$coefSpace), t(cfs[[spp]]$estSpace[1:B,c("pAspen", cnb)]))
cfsm[cfsm > 10^4] <- 10^4
cfsm[cfsm < -10^4] <- -10^4
cfsj[cfsj > 10^4] <- 10^4
cfsj[cfsj < -10^4] <- -10^4
if (ncol(cfsm) < BMAX) {
    b <- sample.int(ncol(cfsm), BMAX-ncol(cfsm), replace=TRUE)
    cfsm <- cbind(cfsm, cfsm[,b])
    cfsj <- cbind(cfsj, cfsj[,b])
}
compare_sets(rownames(cfsm),dimnames(COEFS$lichens$south)[[2]])

CFsm <- array(0, c(length(cfs), nrow(cfsm), BMAX))
dimnames(CFsm) <- list(names(cfs), rownames(cfsm), NULL)
CFsj <- array(0, c(length(cfs), nrow(cfsj), BMAX))
dimnames(CFsj) <- list(names(cfs), rownames(cfsj), NULL)

for (spp in names(cfs)) {
    B <- min(BMAX, nrow(cfs[[spp]]$estSpace))
    cfsm <- rbind(log(cfs[[spp]]$coefARU), pAspen=cfs[[spp]]$estARU[1:B,"pAspen"])
    cfsj <- rbind(log(cfs[[spp]]$coefSpace), t(cfs[[spp]]$estSpace[1:B,c("pAspen", cnb)]))
    cfsm[cfsm > 10^4] <- 10^4
    cfsm[cfsm < -10^4] <- -10^4
    cfsj[cfsj > 10^4] <- 10^4
    cfsj[cfsj < -10^4] <- -10^4
    if (ncol(cfsm) < BMAX) {
        b <- sample.int(ncol(cfsm), BMAX-ncol(cfsm), replace=TRUE)
        cfsm <- cbind(cfsm, cfsm[,b])
        cfsj <- cbind(cfsj, cfsj[,b])
    }
    CFsm[spp,,] <- cfsm
    CFsj[spp,,] <- cfsj
}

# north
spp <- 1
B <- min(BMAX, nrow(cfn[[spp]]$estSpace))
cfnm <- rbind(log(cfn[[spp]]$coefARU))
cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[1:B,cnb]))
cfnm[cfnm > 10^4] <- 10^4
cfnm[cfnm < -10^4] <- -10^4
cfnj[cfnj > 10^4] <- 10^4
cfnj[cfnj < -10^4] <- -10^4
if (ncol(cfnm) < BMAX) {
    b <- sample.int(ncol(cfnm), BMAX-ncol(cfnm), replace=TRUE)
    cfnm <- cbind(cfnm, cfnm[,b])
    cfnj <- cbind(cfnj, cfnj[,b])
}
compare_sets(rownames(cfnm),dimnames(COEFS$lichens$north)[[2]])

CFnm <- array(0, c(length(cfn), nrow(cfnm), BMAX))
dimnames(CFnm) <- list(names(cfn), rownames(cfnm), NULL)
CFnj <- array(0, c(length(cfn), nrow(cfnj), BMAX))
dimnames(CFnj) <- list(names(cfn), rownames(cfnj), NULL)

for (spp in names(cfn)) {
    B <- min(BMAX, nrow(cfn[[spp]]$estSpace))
    cfnm <- rbind(log(cfn[[spp]]$coefARU))
    cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[1:B,cnb]))
    cfnm[cfnm > 10^4] <- 10^4
    cfnm[cfnm < -10^4] <- -10^4
    cfnj[cfnj > 10^4] <- 10^4
    cfnj[cfnj < -10^4] <- -10^4
    if (ncol(cfnm) < BMAX) {
        b <- sample.int(ncol(cfnm), BMAX-ncol(cfnm), replace=TRUE)
        cfnm <- cbind(cfnm, cfnm[,b])
        cfnj <- cbind(cfnj, cfnj[,b])
    }
    CFnm[spp,,] <- cfnm
    CFnj[spp,,] <- cfnj
}

## bare & ice -- this is not really working well, no time to debog and not needed
if (FALSE) {
for (j in names(COEFS)) {
    zn <- COEFS[[j]]$north
    xn <- array(NA, dim(zn)+c(0,2,0))
    dimnames(xn) <- list(dimnames(zn)[[1]],
        c(cnn0, "Bare", "SnowIce", cnn1),
        NULL)
    xn[,cnn,] <- zn
    xn[is.na(xn)] <- -10^4
    COEFS[[j]]$north <- xn
}

zn <- CFnm
xn <- array(NA, dim(zn)+c(0,2,0))
dimnames(xn) <- list(dimnames(zn)[[1]],
    c(cnn0, "Bare", "SnowIce"),
    NULL)
xn[,cnn0,] <- zn
xn[is.na(xn)] <- -10^4
CFnm <- xn

zn <- CFnj
bbb <- dimnames(zn)[[2]]
xn <- array(NA, dim(zn)+c(0,2,0))
dimnames(xn) <- list(dimnames(zn)[[1]],
    c(cnn0, "Bare", "SnowIce", bbb[!(bbb %in% cnn0)]),
    NULL)
xn[,bbb,] <- zn
xn[is.na(xn)] <- -10^4
CFnj <- xn
}

COEFS$birds <- list(
    north=list(marginal=CFnm, joint=CFnj),
    south=list(marginal=CFsm, joint=CFsj),
    species=TAX)

save(COEFS,
    file="s:/AB_data_v2020/Results/COEFS-ALL.RData")

if (FALSE) { # testing CC effects
spp <- "AlderFlycatcher"
h <- c("CCWhiteSpruceR", "CCWhiteSpruce1",
"CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "CCPineR",
"CCPine1", "CCPine2", "CCPine3", "CCPine4", "CCDeciduousR", "CCDeciduous1",
"CCDeciduous2", "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR",
"CCMixedwood1", "CCMixedwood2", "CCMixedwood3", "CCMixedwood4")
h1 <- gsub("CC", "", h)

m  <- matrix(rowMeans(exp(COEFS$birds$north$marginal[spp,,]))[h], 5, 4)
m1 <- matrix(rowMeans(exp(COEFS$birds$north$marginal[spp,,]))[h1], 5, 4)
plot(0, type="n", xlim=c(0,4), ylim=c(0,max(m,m1)))
for (i in 1:4) {
    lines(0:4, m[,i], col=i, lty=2)
    lines(0:4, m1[,i], col=i, lty=1)
}

m  <- matrix(rowMeans(exp(CFnm[spp,,]))[h], 5, 4)
m1 <- matrix(rowMeans(exp(CFnm[spp,,]))[h1], 5, 4)
plot(0, type="n", xlim=c(0,4), ylim=c(0,max(m,m1)))
for (i in 1:4) {
    lines(0:4, m[,i], col=i, lty=2)
    lines(0:4, m1[,i], col=i, lty=1)
}

m  <- matrix(rowMeans(cfn[[spp]]$coefARU)[h], 5, 4)
m1 <- matrix(rowMeans(cfn[[spp]]$coefARU)[h1], 5, 4)
plot(0, type="n", xlim=c(0,4), ylim=c(0,max(m,m1)))
for (i in 1:4) {
    lines(0:4, m[,i], col=i, lty=2)
    lines(0:4, m1[,i], col=i, lty=1)
}

}


## TODO
##
## Coefs
## - Bare, SnowIce: 0 and placeholder
## - add Mammals
## - heads up for 150m mammal and amphibian models
##
## Preds
## - develop alternative code (current only) for bootstrap
## -
