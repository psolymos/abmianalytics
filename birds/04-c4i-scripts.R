library(mefa4)
library(intrval)
#library(raster)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
#ROOT <- "~/GoogleWork/tmp"

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)
TAX <- ee$tax

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

xwalk <- read.csv("~/repos/abmianalytics/lookup/veg-v61-crosswalk.csv")
rownames(xwalk) <- xwalk[,1]

tax <- read.csv("d:/abmi/AB_data_v2018/data/raw/species/taxonomy.csv")

bbn <- unique(sort(as.numeric(en$BB)))
bbs <- unique(sort(as.numeric(es$BB)))

#spp <- "ALFL"
#resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
#ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))

## lookup table
if (FALSE) {
Lookup <- data.frame(
    SpeciesID=TAX$sppid,
    ScientificName=TAX$scinam,
    TSNID=NA,
    CommonName=TAX$species,
    ModelNorth=FALSE,
    ModelSouth=FALSE,
    UseavailNorth=FALSE,
    UseavailSouth=FALSE,
    Nonnative=FALSE,
    LinkHabitat="poisson_log",
    LinkSpclim="poisson_log",
    AUCNorth=NA,
    AUCSouth=NA,
    Code=TAX$code,
    Order=TAX$order,
    Family=TAX$family,
    Comments=""
)
rownames(Lookup) <- Lookup$SpeciesID
Lookup$TSNID <- tax$TSNID[match(Lookup$ScientificName, tax$SCIENTIFIC_NAME)]
Lookup$Comments <- as.character(Lookup$Comments)

tmp1 <- c("Accipitridae", "Falconidae", "Pandionidae")
tmp2 <- "Raptor species with large home ranges are inadequately sampled by point counts and ARUs."
Lookup$Comments[Lookup$Family %in% tmp1] <- tmp2

tmp1 <- c("Anatidae", "Gaviidae", "Gruidae", "Laridae", "Pelecanidae", "Phalacrocoracidae",
    "Podicipedidae", "Rallidae", "Recurvirostridae")
tmp2 <- "Aquatic dependent species species are inadequately sampled by terrestrial point counts and ARUs."
Lookup$Comments[Lookup$Family %in% tmp1] <- tmp2

tmp1 <- c("Caprimulgidae", "Strigidae", "Tytonidae")
tmp2 <- "Owls and other nocturnal species are inadequately sampled by point counts and ARUs with traditional timing of early morning."
Lookup$Comments[Lookup$Family %in% tmp1] <- tmp2

#Phasianidae
#Game birds are not modeled.
}

bbb <- read.csv("d:/abmi/sppweb2018/c4i/tables/StandardizedOutput-birds-final-lookup-withChecks.csv")
bbb <- bbb[is.na(bbb$Exclude),]
rownames(bbb) <- bbb$SpeciesID
bbb$Code <- bbb$Code.1
Lookup <- bbb[,c("SpeciesID", "ScientificName", "TSNID", "CommonName", "ModelNorth",
    "ModelSouth", "UseavailNorth", "UseavailSouth", "Nonnative",
    "LinkHabitat", "LinkSpclim", "AUCNorth", "AUCSouth", "Code",
    "Order", "Family", "Comments")]
Lookup <- Lookup[Lookup$UseavailNorth+Lookup$UseavailSouth+Lookup$ModelNorth+Lookup$ModelSouth > 0,]
Lookup <- droplevels(Lookup)

## Use-Availability south

yys <- ee$yy[rownames(es$DAT),]
yys <- yys[,colSums(yys>0) > 3]
vcs <- ee$sc2[rownames(yys),] # 1km^2 buffer used here
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]
vcs <- groupSums(vcs, 2, ts[colnames(vcs), "UseInAnalysisCoarse"])
vcs <- vcs[rownames(yys),
    c("Productive",
    "Clay",
    "Saline",
    "RapidDrain",
    "Crop",
    "TameP",
    "RoughP",
    "UrbInd",
    "HardLin",
    "SoftLin")]
vcs <- vcs / rowSums(vcs)
UseavailSouth <- t(apply(as.matrix(yys), 2, function(z) wrsi(z, vcs)$rWRSI))
colnames(UseavailSouth) <- colnames(vcs)
rownames(UseavailSouth) <- TAX[rownames(UseavailSouth), "sppid"]

#Lookup$UseavailSouth[Lookup$SpeciesID %in% rownames(UseavailSouth)] <- TRUE

## Use-Availability north

yyn <- ee$yy[rownames(en$DAT),]
yyn <- yyn[,colSums(yyn>0) > 3]
vcn <- ee$vc2[rownames(yyn),] # 1km^2 buffer used here
vt <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(vt) <- vt$ID
vcn <- groupSums(vcn, 2, vt[colnames(vcn), "UseAvail"])
vcn <- vcn[,c(
    "Deciduous",
    "Mixedwood",
    "WhiteSpruce",
    "Pine",
    "BlackSpruce",
    "TreedFen",
    "Open",
    "Wetland",
    "Water",
    "Bare",
    "HFor",
    "Crop",
    "TameP",
    "RoughP",
    "UrbInd",
    "HardLin",
    "SoftLin")]
UseavailNorth <- t(apply(as.matrix(yyn), 2, function(z) wrsi(z, vcn)$rWRSI))
colnames(UseavailNorth) <- colnames(vcn)
rownames(UseavailNorth) <- TAX[rownames(UseavailNorth), "sppid"]

#Lookup$UseavailNorth[Lookup$SpeciesID %in% rownames(UseavailNorth)] <- TRUE

## coefs North

get_coef_north <- function(resn, STAGE="ARU", subset=NULL, ...) {
    OK <- !sapply(resn, inherits, "try-error")
    spp <- resn[[which(OK)[1]]]$species
    estn <- suppressWarnings(get_coef(resn, Xn, stage=STAGE, na.out=FALSE))
    if (is.null(subset))
        subset <- seq_len(nrow(estn))
    estn <- estn[subset,,drop=FALSE]

    mu <- Xage %*% t(estn[,colnames(Xage),drop=FALSE])
    lam1 <- exp(mu)
    lam1 <- lam1[!grepl("9", rownames(lam1)),,drop=FALSE]
    #lamUI <- lam1[xwalk[rownames(lam1),2] == "UrbInd",]
    #lam1 <- groupMeans(lam1, 1, xwalk[rownames(lam1),2])
    #lam1 <- lam1[rownames(lam1) != "UrbInd",,drop=FALSE]
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

    ## UrbInd
    #Xn2 <- Xn[en$DAT$vegc %ni% c("Industrial", "Mine", "Rural", "Urban"), colnames(Xage)]
    #nUI <- nrow(Xn2)
    #lamUI <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    #estUI <- quantile((nUI * lamUI + nWell * lamWell) / (nUI + nWell),
    #    c(0.5, 0.05, 0.95))

    lam1 <- rbind(lam1,
        #UrbInd=estUI,
        SoftLin=(nSeism * estSeism + nEnSoft * estEnSft) / (nSeism + nEnSoft),
        HardLin=estTrSft)
    ## add here what's needed for final output as defined (seismic, ensoft, trsoft, hardlin, well)
    colnames(lam1) <- c("Estimate", "Lower", "Upper")
    lam1
}

COEFS <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelNorth])) {
    cat(spp, "\n")
    flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    tmp <- try(get_coef_north(resn, STAGE="ARU"))
    if (!inherits(tmp, "try-error"))
        COEFS[[spp]] <- tmp
}

COEFSmarginal <- list()
COEFSjoint <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelNorth])) {
    cat(spp, "\n")
    flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    tmp <- try(get_coef_north(resn, STAGE="ARU", subset=1))
    if (!inherits(tmp, "try-error"))
        COEFSmarginal[[spp]] <- tmp
    tmp <- try(get_coef_north(resn, STAGE="Space", subset=1))
    if (!inherits(tmp, "try-error"))
        COEFSjoint[[spp]] <- tmp
}

plot_coef_north <- function(lam1, ...) {

    if (length(lam1) == 1 && is.character(lam1)) {
        spp <- lam1
        spp2 <- rownames(Lookup)[Lookup$Code == spp]
        lam1 <- cbind(
            Estimate=CoefNorth[spp2,,1],
            Lower=LowerNorth[spp2,],
            Upper=UpperNorth[spp2,])
        lam1 <- lam1[rownames(lam1) != "Water",]
        lam1 <- lam1[rownames(lam1) != "Bare",]
    }

    lam <- lam1[!grepl("CC", rownames(lam1)),]
    #lam <- lam[!(rownames(lam) %in% c("SoftLin", "HardLin")),]
    rownames(lam) <- gsub("_", " ", rownames(lam))
    lamCC <- lam1[grepl("CC", rownames(lam1)),]

    op <- par(las=2, mar=c(10,4,3,2), cex.axis=0.9)
    on.exit(par(op))
    col <- c(rep(RColorBrewer::brewer.pal(8, "Accent")[c(1,2,3,5,6)], each=9),
        RColorBrewer::brewer.pal(8, "Dark2"),
        RColorBrewer::brewer.pal(8, "Set3")[c(1,3,4)])
    k <- 2
    b <- barplot(unname(lam[,1]), col=col,
        ylab="Relative abundance",
        ylim=c(0, min(k*max(lam1[,1]), max(lam1))), ...)
    mtext(rownames(lam), col=col, side=1, at=b, cex=0.8, line=1)
    o <- b[2]-b[1]
    for (i in 1:4) {
        xi <- b[(((i-1)*9+1):(i*9))][1:6]
        yi <- rbind(lamCC[((i-1)*5+1):(i*5),], lam[(((i-1)*9+1):(i*9))[6],,drop=FALSE])
        lines(xi-0.2*o, yi[,1], lwd=1, lty=1, col="darkgrey")
        points(xi[-6]-0.2*o, yi[-6,1], col=1, lwd=2, cex=0.6)
        segments(x0=xi[-6]-0.2*o, y0=yi[-6,2], y1=yi[-6,3], lwd=1, col=1)
    }
    segments(x0=b, y0=lam[,2], y1=lam[,3], lwd=1, col=1)

    invisible(lam1)
}

## also look at current and previous results
library(cure4insect)
set_options(path = "s:/reports", version="2017")
load_common_data()
v17spp <- get_all_species("birds", "north")

pdf(file.path(ROOT, "explore-coefNorth2.pdf"), onefile=TRUE, width=10, height=12)
for (spp in names(COEFS)) {
    cat(spp, "\n");flush.console()
    y <- en$YY[bbn,spp]
    npk1 <- sum(y>0)
    npkT <- length(y)
    yss <- sum_by(y, en$DAT$SS[bbn])[,"x"]
    nss1 <- sum(yss>0)
    nssT <- length(yss)
    NAM <- as.character(TAX[spp, "species"])
    LAB <- paste0(spp, " ", NAM, ", North\n", npk1,"/",npkT, " pts, ", nss1, "/", nssT, " loc")
    ID <- as.character(TAX[spp, "sppid"])
    op <- par(mfrow=c(2,1))
    plot_coef_north(COEFS[[spp]], main=LAB)
    if (ID %in% v17spp)
        plot_abundance(ID, "veg_coef") else plot.new()
    par(op)
}
dev.off()

cn1 <- c(
    "WhiteSpruce_0-10",
    "WhiteSpruce_10-20",
    "WhiteSpruce_20-40",
    "WhiteSpruce_40-60",
    "WhiteSpruce_60-80",
    "WhiteSpruce_80-100",
    "WhiteSpruce_100-120",
    "WhiteSpruce_120-140",
    "WhiteSpruce_140+",
    "Pine_0-10",
    "Pine_10-20",
    "Pine_20-40",
    "Pine_40-60",
    "Pine_60-80",
    "Pine_80-100",
    "Pine_100-120",
    "Pine_120-140",
    "Pine_140+",
    "Deciduous_0-10",
    "Deciduous_10-20",
    "Deciduous_20-40",
    "Deciduous_40-60",
    "Deciduous_60-80",
    "Deciduous_80-100",
    "Deciduous_100-120",
    "Deciduous_120-140",
    "Deciduous_140+",
    "Mixedwood_0-10",
    "Mixedwood_10-20",
    "Mixedwood_20-40",
    "Mixedwood_40-60",
    "Mixedwood_60-80",
    "Mixedwood_80-100",
    "Mixedwood_100-120",
    "Mixedwood_120-140",
    "Mixedwood_140+",
    "BlackSpruce_0-10",
    "BlackSpruce_10-20",
    "BlackSpruce_20-40",
    "BlackSpruce_40-60",
    "BlackSpruce_60-80",
    "BlackSpruce_80-100",
    "BlackSpruce_100-120",
    "BlackSpruce_120-140",
    "BlackSpruce_140+",
    "TreedFen",
    "Shrub",
    "Grass",
    "Bare",
    "Water",
    "TreeShrubSwamp",
    "NonTreeFenMarsh")
cn2 <- c("CCWhiteSpruce_0-10",
    "CCWhiteSpruce_10-20",
    "CCWhiteSpruce_20-40",
    "CCWhiteSpruce_40-60",
    "CCWhiteSpruce_60-80",
    "CCPine_0-10",
    "CCPine_10-20",
    "CCPine_20-40",
    "CCPine_40-60",
    "CCPine_60-80",
    "CCDeciduous_0-10",
    "CCDeciduous_10-20",
    "CCDeciduous_20-40",
    "CCDeciduous_40-60",
    "CCDeciduous_60-80",
    "CCMixedwood_0-10",
    "CCMixedwood_10-20",
    "CCMixedwood_20-40",
    "CCMixedwood_40-60",
    "CCMixedwood_60-80",
    "Crop",
    "TameP",
    "RoughP",
    "UrbInd",
    "SoftLin",
    "HardLin")
cn <- c(cn1, cn2)

CoefNorth <- t(sapply(COEFS, function(z) {
    c(z[,"Estimate"], Bare=NA, Water=NA)
}))[,cn]
LowerNorth <- t(sapply(COEFS, function(z) {
    c(z[,"Lower"], Bare=NA, Water=NA)
}))[,cn]
UpperNorth <- t(sapply(COEFS, function(z) {
    c(z[,"Upper"], Bare=NA, Water=NA)
}))[,cn]
rownames(CoefNorth) <- rownames(LowerNorth) <- rownames(UpperNorth) <-
    Lookup$SpeciesID[match(rownames(CoefNorth), Lookup$Code)]
#Lookup$ModelNorth[Lookup$SpeciesID %in% rownames(CoefNorth)] <- TRUE

CoefNorth <- array(CoefNorth, c(dim(CoefNorth), 1), dimnames=c(dimnames(CoefNorth), NULL))

CoefNorthJoint <- t(sapply(COEFSjoint, function(z) {
    c(z[,"Estimate"], Bare=NA, Water=NA)
}))[,cn]
CoefNorthMarginal <- t(sapply(COEFSmarginal, function(z) {
    c(z[,"Estimate"], Bare=NA, Water=NA)
}))[,cn]


## Linear North

LinearNorth <- data.frame(AverageCoef=rowMeans(CoefNorth[,cn1,1], na.rm=TRUE))
LinearNorth$SoftLin10 <- LinearNorth$AverageCoef * 0.9 + CoefNorth[,"SoftLin",1] * 0.1
LinearNorth$HardLin10 <- LinearNorth$AverageCoef * 0.9 + CoefNorth[,"HardLin",1] * 0.1
LinearNorth <- as.matrix(LinearNorth)


## Coef South

cn1 <- c(
    "Productive",
    "Clay",
    "Saline",
    "RapidDrain",
    "Water")
cn2 <- c("Crop",
    "TameP",
    "RoughP",
    "UrbInd",
    "HardLin",
    "SoftLin")
cn <- c(cn1, cn2)

get_coef_south <- function(ress, STAGE="ARU", subset=NULL, ...) {
    ests <- suppressWarnings(get_coef(ress, Xs, stage=STAGE, na.out=FALSE))
    if (is.null(subset))
        subset <- seq_len(nrow(ests))
    ests <- ests[subset,,drop=FALSE]

#    LCC <- c("soilcClay", "soilcCrop", "soilcRapidDrain",
#        "soilcRoughP", "soilcSaline", "soilcTameP", "soilcUrbInd")
    LCC <- c("soilcBlowout", "soilcClaySub", "soilcCrop",
        "soilcIndustrial", "soilcMine", "soilcOther", "soilcRapidDrain",
        "soilcRoughP", "soilcRural", "soilcSandyLoam", "soilcTameP",
        "soilcThinBreak", "soilcUrban")
    muLCC <- t(cbind(ests[,1,drop=FALSE], ests[,1]+ests[,LCC,drop=FALSE]))
    rownames(muLCC) <- levels(es$DAT$soilc)
    qfun <- function(x)
        c(quantile(x, c(0.5, 0.05, 0.95)))
    lamLCC0 <- t(apply(exp(muLCC), 1, qfun)) # nontred
    # RoughP: Abandoned and RoughPasture (closer to natural grass)
    # TameP: TamePasture + HD livestock operations (closer to crop)
#    o <- c("Productive", "Clay", "Saline", "RapidDrain", "RoughP", "TameP", "Crop", "UrbInd")
    o <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other",
        "RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")
    lamLCC0 <- lamLCC0[o,]

    lamLCC1 <- t(apply(exp(muLCC + ests[,"pAspen"]), 1, qfun)) # treed
    lamLCC1 <- lamLCC1[o,]

    MOD <- c("ROAD", "mWell", "mSoft")
    Z <- exp(ests[,MOD,drop=FALSE])
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2)
    for (i in MOD)
        Z[,i] <- linexp(1, ests[,i], pm[i])
    lamMOD <- t(apply(Z, 2, qfun))
    rownames(lamMOD) <- c("Road", "Well", "Soft")

    HFc <- c("RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")

    Xs2 <- Xs[es$DAT$mSoft > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nSoft <- nrow(Xs2)
    lamSoft <- apply(exp(Xs2 %*% t(ests[,colnames(Xs2),drop=FALSE])), 2, median)
    estSoft <- quantile(lamSoft * Z[,"mSoft"], c(0.5, 0.05, 0.95))

    Xs2 <- Xs[es$DAT$mWell > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nWell <- nrow(Xs2)
    lamWell <- apply(exp(Xs2 %*% t(ests[,colnames(Xs2),drop=FALSE])), 2, median)
    estWell <- quantile(lamWell * Z[,"mSoft"], c(0.5, 0.05, 0.95))

    Xs2 <- Xs[es$DAT$ROAD > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nROAD <- nrow(Xs2)
    lamROAD <- apply(exp(Xs2 %*% t(ests[,colnames(Xs2),drop=FALSE])), 2, median)
    estROAD <- quantile(lamROAD * Z[,"mSoft"], c(0.5, 0.05, 0.95))

    lam0 <- rbind(lamLCC0,
        EnSeismic=estSoft,
        EnSoftLin=estSoft,
        TrSoftLin=estSoft,
        HardLin=estROAD,
        Well=estWell)
    colnames(lam0) <- c("Estimate", "Lower", "Upper")
    attr(lam0, "pAspen") <- mean(exp(ests[,"pAspen"]))
    attr(lam0, "Treed") <- lamLCC1
    lam0
}

COEFS2 <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelSouth])) {
    cat(spp, "\n")
    flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    tmp <- try(get_coef_south(ress, STAGE="ARU"))
    if (!inherits(tmp, "try-error"))
        COEFS2[[spp]] <- tmp
}

COEFS2joint <- list()
COEFS2marginal <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelSouth])) {
    cat(spp, "\n")
    flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    tmp <- try(get_coef_south(ress, STAGE="ARU", subset=1))
    if (!inherits(tmp, "try-error"))
        COEFS2marginal[[spp]] <- tmp
    tmp <- try(get_coef_south(ress, STAGE="Space", subset=1))
    if (!inherits(tmp, "try-error"))
        COEFS2joint[[spp]] <- tmp
}


plot_coef_south <- function(lam0, ...) {

    if (length(lam0) == 1 && is.character(lam0)) {
        spp <- lam0
        spp2 <- rownames(Lookup)[Lookup$Code == spp]
        lam0 <- cbind(
            Estimate=CoefSouth[spp2,,1],
            Lower=LowerSouth[spp2,],
            Upper=UpperSouth[spp2,])
        lam0 <- lam0[rownames(lam0) != "Water",]
        attr(lam0, "pAspen") <- mean(exp(CoefSouthBootlist[[spp]][,"pAspen"]))
    }

#    op <- par(mfrow=c(1,2), las=2, cex.axis=0.9)
#    on.exit(par(op))
    lam1 <- lam0 * attr(lam0, "pAspen")
    k <- 2
    b <- barplot(lam0[,1], ylim=c(0, min(k*max(lam0[,1], lam1[,1]), max(lam0, lam1))),
        ylab="Relative abundance",
        col=RColorBrewer::brewer.pal(8, "Dark2"), ...)
    segments(x0=b, y0=lam0[,2], y1=lam0[,3], lwd=2, col=1)

    b <- barplot(lam1[,1], ylim=c(0, min(k*max(lam0[,1], lam1[,1]), max(lam0, lam1))),
        ylab="Relative abundance",
        col=RColorBrewer::brewer.pal(8, "Dark2"), main="Treed")
    segments(x0=b, y0=lam1[,2], y1=lam1[,3], lwd=2, col=1)

    #lam1x <- rbind(attr(lam0, "Treed"), lam0[9:10,] * attr(lam0, "pAspen"))
    #b <- barplot(lam1x[,1], ylim=c(0, min(k*max(lam0[,1], lam1[,1]), max(lam0, lam1))),
    #    ylab="Relative abundance",
    #    col=RColorBrewer::brewer.pal(8, "Dark2"), main="Treed, joint")
    #segments(x0=b, y0=lam1x[,2], y1=lam1x[,3], lwd=2, col=1)

    invisible(lam0)
}

s17spp <- get_all_species("birds", "south")

pdf(file.path(ROOT, "explore-coefSouth2.pdf"), onefile=TRUE, width=20, height=10)
for (spp in names(COEFS2)) {
    cat(spp, "\n");flush.console()
    y <- es$YY[bbs,spp]
    npk1 <- sum(y>0)
    npkT <- length(y)
    yss <- sum_by(y, es$DAT$SS[bbs])[,"x"]
    nss1 <- sum(yss>0)
    nssT <- length(yss)
    NAM <- as.character(TAX[spp, "species"])
    LAB <- paste0(spp, " ", NAM, ", South\n", npk1,"/",npkT, " pts, ", nss1, "/", nssT, " loc")
    ID <- as.character(TAX[spp, "sppid"])
    #op <- par(mfrow=c(2,1))
    op <- par(mfrow=c(2,2), cex.axis=0.7)
    plot_coef_south(COEFS2[[spp]], main=LAB)
    if (ID %in% s17spp) {
        p0 <- plot_abundance(ID, "soil_coef", paspen=0, plot=FALSE)
        p1 <- plot_abundance(ID, "soil_coef", paspen=1, plot=FALSE)
        plot_abundance(ID, "soil_coef", paspen=0, ylim=c(0, max(p0, p1)))
        plot_abundance(ID, "soil_coef", paspen=1, ylim=c(0, max(p0, p1)))
    }
    par(op)

}
dev.off()

CoefSouth <- t(sapply(COEFS2, function(z) {
    c(z[,"Estimate"], Water=NA)
}))[,cn]
LowerSouth <- t(sapply(COEFS2, function(z) {
    c(z[,"Lower"], Water=NA)
}))[,cn]
UpperSouth <- t(sapply(COEFS2, function(z) {
    c(z[,"Upper"], Water=NA)
}))[,cn]
rownames(CoefSouth) <- rownames(LowerSouth) <- rownames(UpperSouth) <-
    Lookup$SpeciesID[match(rownames(CoefSouth), Lookup$Code)]
#Lookup$ModelSouth[Lookup$SpeciesID %in% rownames(CoefSouth)] <- TRUE

CoefSouth <- array(CoefSouth, c(dim(CoefSouth), 1), dimnames=c(dimnames(CoefSouth), NULL))

CoefSouthJoint <- t(sapply(COEFS2joint, function(z) {
    c(z[,"Estimate"], Water=NA)
}))[,cn]
CoefSouthMarginal <- t(sapply(COEFS2marginal, function(z) {
    c(z[,"Estimate"], Water=NA)
}))[,cn]


## Linear South

LinearSouth <- data.frame(AverageCoef=rowMeans(CoefSouth[,cn1,1], na.rm=TRUE))
LinearSouth$SoftLin10 <- LinearSouth$AverageCoef * 0.9 + CoefSouth[,"SoftLin",1] * 0.1
LinearSouth$HardLin10 <- LinearSouth$AverageCoef * 0.9 + CoefSouth[,"HardLin",1] * 0.1
LinearSouth <- as.matrix(LinearSouth)


## AUC based on stage=Space
Lookup$AUCNorth <- NA
for (spp in names(COEFS)) {
    cat(spp, "N\n");flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    yn <- as.numeric(en$YY[,spp])
    off <- if (spp %in% colnames(en$OFF))
        en$OFF[,spp] else en$OFFmean
    lamn <- exp(predict_with_SSH(resn, Xn, en$SSH, stage="Space") + off)
    rocn <- simple_roc(ifelse(yn > 0, 1, 0), rowMeans(lamn))
    aucn <- simple_auc(rocn)
    Lookup$AUCNorth[Lookup$Code == spp] <- aucn
}

Lookup$AUCSouth <- NA
for (spp in names(COEFS2)) {
    cat(spp, "S\n");flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    ys <- as.numeric(es$YY[,spp])
    off <- if (spp %in% colnames(es$OFF))
        es$OFF[,spp] else es$OFFmean
    lams <- exp(predict_with_SSH(ress, Xs, es$SSH, stage="Space") + off)
    rocs <- simple_roc(ifelse(ys > 0, 1, 0), rowMeans(lams))
    aucs <- simple_auc(rocs)
    Lookup$AUCSouth[Lookup$Code == spp] <- aucs
}

## additions

Lookup$SizeNorth <- colSums(en$YY[bbn,] > 0)[match(Lookup$Code, colnames(en$YY))]
Lookup$SizeSouth <- colSums(es$YY[bbs,] > 0)[match(Lookup$Code, colnames(es$YY))]

tmpn <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelNorth])) {
    cat(spp, "\n");flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    estn <- suppressWarnings(get_coef(resn, Xn, stage="ARU", na.out=FALSE))
    tmpn[[spp]] <- estn
    #CoefNorthBoot[spp,,] <- estn[1:101,dimnames(CoefNorthBoot)[[2]]]
}
tmps <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelSouth])) {
    cat(spp, "\n");flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    ests <- suppressWarnings(get_coef(ress, Xs, stage="ARU", na.out=FALSE))
    tmps[[spp]] <- ests
    #CoefSouthBoot[spp,,1:min(nrow(ests),101)] <- ests[1:min(nrow(ests),101),dimnames(CoefSouthBoot)[[2]]]
}
CoefNorthBootlistARU <- tmpn
CoefSouthBootlistARU <- tmps

tmpn <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelNorth])) {
    cat(spp, "\n");flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    estn <- suppressWarnings(get_coef(resn, Xn, stage="Space", na.out=FALSE))
    tmpn[[spp]] <- estn
    #CoefNorthBoot[spp,,] <- estn[1:101,dimnames(CoefNorthBoot)[[2]]]
}
tmps <- list()
for (spp in as.character(Lookup$Code[Lookup$ModelSouth])) {
    cat(spp, "\n");flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    ests <- suppressWarnings(get_coef(ress, Xs, stage="Space", na.out=FALSE))
    tmps[[spp]] <- ests
    #CoefSouthBoot[spp,,1:min(nrow(ests),101)] <- ests[1:min(nrow(ests),101),dimnames(CoefSouthBoot)[[2]]]
}
CoefNorthBootlistSpace <- tmpn
CoefSouthBootlistSpace <- tmps

toSave <- c("Lookup",
    "CoefNorth", "CoefSouth",
    "CoefNorthJoint", "CoefSouthJoint",
    "CoefNorthMarginal", "CoefSouthMarginal",
    "CoefNorthBootlistSpace", "CoefSouthBootlistSpace",
    "CoefNorthBootlistARU", "CoefSouthBootlistARU",
    #"SpclimNorth", "SpclimSouth",
    "UseavailNorth", "UseavailSouth",
    "LinearNorth", "LinearSouth",
    "UpperNorth", "UpperSouth",
    "LowerNorth", "LowerSouth")
for (i in toSave) {
    cat("\n\n-------------\n", i, "\n\n")
    print(str(get(i)))
}

table(rowSums(Lookup[,c("ModelNorth", "ModelSouth", "UseavailNorth", "UseavailSouth")]))
colSums(Lookup[,c("ModelNorth", "ModelSouth", "UseavailNorth", "UseavailSouth")])
table(N=Lookup$ModelNorth, S=Lookup$ModelSouth)
summary(Lookup$AUCNorth)
summary(Lookup$AUCSouth)

save(list=toSave,
    file="d:/abmi/sppweb2018/c4i/tables/StandardizedOutput-birds.RData")

## --- checking results


library(mefa4)
ROOT <- "~/GoogleWork/tmp"

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
bbn <- unique(sort(as.numeric(en$BB)))
bbs <- unique(sort(as.numeric(es$BB)))

load(file.path(ROOT, "tables", "StandardizedOutput-birds.RData"))

## N 109, S 101
b <- read.csv("~/repos/abmispecies/_data/birds.csv")

LIN <- c("SoftLin", "HardLin")
#dCIn <- rowMeans((UpperNorth - LowerNorth) / (0.00001 + CoefNorth[,,1]), na.rm=TRUE)
dCIn <- apply(UpperNorth[,!(colnames(UpperNorth) %in% LIN)], 1, max, na.rm=TRUE) /
    apply(CoefNorth[,!(colnames(UpperNorth) %in% LIN), 1], 1, max, na.rm=TRUE)
dCIn <- structure(pmin(25, round(dCIn, 1)), names=names(dCIn))
#dCIs <- rowMeans((UpperSouth - LowerSouth) / (0.00001 + CoefSouth[,,1]), na.rm=TRUE)
dCIs <- apply(UpperSouth[,!(colnames(UpperSouth) %in% LIN)], 1, max, na.rm=TRUE) /
    apply(CoefSouth[,!(colnames(UpperSouth) %in% LIN),1], 1, max, na.rm=TRUE)
dCIs <- structure(pmin(25, round(dCIs, 1)), names=names(dCIs))

nn <- Lookup[names(dCIn), "SizeNorth"]
ns <- Lookup[names(dCIs), "SizeSouth"]
SPPn <- as.character(Lookup[names(dCIn), "Code"])
SPPs <- as.character(Lookup[names(dCIs), "Code"])
dfn <- sapply(SPPn, function(i) sum(colSums(abs(CoefNorthBootlist[[i]])) > 0))
dfs <- sapply(SPPs, function(i) sum(colSums(abs(CoefSouthBootlist[[i]])) > 0))
ncn <- Lookup[names(dCIn), "Comments"] == ""
ncs <- Lookup[names(dCIs), "Comments"] == ""
names(ncn) <- names(nn) <- names(dfn) <- names(SPPn) <- names(dCIn)
names(ncs) <- names(ns) <- names(dfs) <- names(SPPs) <- names(dCIs)
summary(dfn)
summary(dfs)

plot(dCIn, log(nn), #xlim=c(0,20),
    col=ifelse(Lookup[names(dCIn), "Comments"] == "", 1, 2))
abline(h=log(60*5), v=3)

plot(dCIs, log(ns), #xlim=c(0,20),
    col=ifelse(Lookup[names(dCIs), "Comments"] == "", 1, 2))
abline(h=log(30*5), v=5)

table(CI=dCIn < 3, NoComm=ncn)
table(CI=nn >= 300, NoComm=ncn)

#OKn <- nn >= 300 & ncn & dCIn < 20
#OKs <- ns >= 150 & ncs & dCIs < 20
OKn <- nn >= dfn*5 & dCIn < 20
OKs <- ns >= dfs*5 & dCIs < 20

pdf(file.path(ROOT, "tables", "explore-coefNorthOK.pdf"),
    onefile=TRUE, width=10, height=8)
for (spp in SPPn[OKn]) {
    #cat(spp, "\n");flush.console()
    NAM <- names(SPPn[SPPn == spp])
    LAB <- paste0(spp, " ", NAM, "\n",
        ifelse(nn[NAM] < dfn[NAM]*5, "nnn ", ""),
        ifelse(dCIn[NAM] > 3, "iii ", ""),
        ifelse(Lookup[NAM, "Comments"] == "", "",  "ccc "),
        "ndet=", nn[NAM], " df=", dfn[NAM])
    plot_coef_north(spp, main=LAB)
}
dev.off()
pdf(file.path(ROOT, "tables", "explore-coefNorthNotOK.pdf"),
    onefile=TRUE, width=10, height=8)
for (spp in SPPn[!OKn]) {
    #cat(spp, "\n");flush.console()
    NAM <- names(SPPn[SPPn == spp])
    LAB <- paste0(spp, " ", NAM, "\n",
        ifelse(nn[NAM] < dfn[NAM]*5, "nnn ", ""),
        ifelse(dCIn[NAM] > 3, "iii ", ""),
        ifelse(Lookup[NAM, "Comments"] == "", "",  "ccc "),
        "ndet=", nn[NAM], " df=", dfn[NAM])
    plot_coef_north(spp, main=LAB)
}
dev.off()

pdf(file.path(ROOT, "tables", "explore-coefSouthOK.pdf"),
    onefile=TRUE, width=10, height=5)
for (spp in SPPs[OKs]) {
    #cat(spp, "\n");flush.console()
    NAM <- names(SPPs[SPPs == spp])
    LAB <- paste0(spp, " ", NAM, "\n",
        ifelse(ns[NAM] < dfs[NAM]*5, "nnn ", ""),
        ifelse(dCIs[NAM] > 3, "iii ", ""),
        ifelse(Lookup[NAM, "Comments"] == "", "",  "ccc "),
        "ndet=", ns[NAM], " df=", dfs[NAM])
    plot_coef_south(spp, main=LAB)
}
dev.off()
pdf(file.path(ROOT, "tables", "explore-coefSouthNotOK.pdf"),
    onefile=TRUE, width=10, height=5)
for (spp in SPPs[!OKs]) {
    #cat(spp, "\n");flush.console()
    NAM <- names(SPPs[SPPs == spp])
    LAB <- paste0(spp, " ", NAM, "\n",
        ifelse(ns[NAM] < dfs[NAM]*5, "nnn ", ""),
        ifelse(dCIs[NAM] > 3, "iii ", ""),
        ifelse(Lookup[NAM, "Comments"] == "", "",  "ccc "),
        "ndet=", ns[NAM], " df=", dfs[NAM])
    plot_coef_south(spp, main=LAB)
}
dev.off()

OKn <- nn >= dfn*5 & dCIn < 20
OKs <- ns >= dfs*5 & dCIs < 20

ii <- c("FRGU", "MALL")
OKn[names(SPPn[SPPn %in% ii])] <- FALSE

ii <- c("BADO", "BLBW", "BBWO", "BOOW", "CSWA", "DUFL", "EAPH",
    "EUST", "GADW", "MAWR", "MERL", "MOBL", "MOCH", "PIGR",
    "RECR", "SPGR", "SPPI", "TOWA", "VEER", "WIFL", "YERA",
    "YHBL")
OKn[names(SPPn[SPPn %in% ii])] <- TRUE

ii <- c("FEHA", "LESC")
OKs[names(SPPs[SPPs %in% ii])] <- FALSE

ii <- c("BOBO", "LASP", "PUMA", "RNPH")
OKs[names(SPPs[SPPs %in% ii])] <- TRUE

Lookup$ModelNorth <- rownames(Lookup) %in% names(OKn[OKn])
Lookup$ModelSouth <- rownames(Lookup) %in% names(OKs[OKs])

save(Lookup, file=file.path(ROOT, "tables", "lookup-birds.RData"))


