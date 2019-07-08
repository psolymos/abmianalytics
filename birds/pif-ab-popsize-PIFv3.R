## NEW STUFF FOR v3 ===============================================================================

## find intersect between PIX3 and PIF3


library(mefa4)
library(intrval)
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds"
source("~/repos/abmianalytics/birds/00-functions.R")

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
Xn <- get_model_matrix(en$DAT, en$mods)

## find common subset

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth,])
rownames(tax) <- tax$Code

fl <- list.files("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/")
Spp <- gsub(".RData", "", fl)
tax <- droplevels(tax[Spp,])

pif3 <- read.csv("~/GoogleWork/bam/PIF-AB/popBCR-6AB_v3_27-Feb-2019.csv") #v3
pif3$English.Name <- as.character(pif3$English.Name)
pif3$English.Name[pif3$English.Name == "LeConte's Sparrow"] <- "Le Conte's Sparrow"
pif3$English.Name[pif3$English.Name == "Canada Jay"] <- "Gray Jay"

mefa4::compare_sets(tax$CommonName, pif3$English.Name)
setdiff(tax$CommonName, pif3$English.Name)
tmp <- droplevels(pif3[!(pif3$English.Name %in% tax$CommonName),])
mefa4::compare_sets(tax$ScientificName, tmp$Scientific.Name)

## big table to store all the results

TAB <- droplevels(tax[tax$CommonName %in% pif3$English.Name,])
#TAB$AUC <- TAB$AUCNorth
TAB$AUCNorth <- NULL
TAB$nDet <- TAB$SizeNorth
TAB$SizeNorth <- NULL
TAB$ModelNorth <- TAB$ModelSouth <- TAB$UseavailNorth <- TAB$UseavailSouth <- NULL
TAB$Nonnative <- TAB$LinkHabitat <- TAB$LinkSpclim <- TAB$AUCSouth <- NULL
TAB$SizeSouth <- TAB$Comments <- TAB$TSNID <- NULL
#TAB$SpeciesID <- NULL
summary(TAB)
SPP <- rownames(TAB)

## PIF v3 numbers

pif3 <- droplevels(pif3[match(TAB$CommonName, pif3$English.Name),])
TAB$Npif <- pif3$Population.Estimate..unrounded. / 10^6
TAB$Npif95lower <- pif3$Lower.95..bound..unrounded. / 10^6
TAB$Npif95upper <- pif3$Upper.95..bound..unrounded. / 10^6
TAB$MDD <- pif3$Detection.Distance.Category..m.
TAB$PairAdj <- pif3$Pair.Adjust.Category
TAB$TimeAdj <- pif3$Time.Adjust.Mean

## AUC
if (FALSE) {

    AUC <- numeric(length(SPP))
    names(AUC) <- SPP
    for (spp in SPP) {
        cat(spp, "N\n");flush.console()
        resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
        yn <- as.numeric(en$YY[,spp])
        off <- en$OFF[,spp]
        lamn <- exp(predict_with_SSH(resn, Xn, en$SSH, stage="Space") + off)
        rocn <- simple_roc(ifelse(yn > 0, 1, 0), apply(lamn, 1, median))
        aucn <- simple_auc(rocn)
        AUC[spp] <- aucn
    }
    save(AUC, file="~/GoogleWork/bam/PIF-AB/results-v3/AUC.RData")
}
load("~/GoogleWork/bam/PIF-AB/results-v3/AUC.RData") # AUC
TAB$AUC <- AUC[match(rownames(TAB), names(AUC))]

## PIX pop sizes
if (FALSE) {
    load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
    k <- kgrid[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",]
    A <- k$Area_km2 - k$Area_km2*k$pWater

    NPIX <- matrix(NA, length(SPP), 240)
    rownames(NPIX) <- SPP

    for (spp in SPP) {
        cat(spp, "\n")
        flush.console()
        e <- new.env()
        load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/", spp, ".RData"), envir=e)
        CR <- e$CR[rownames(k),]
        for (i in 1:ncol(CR)) {
            q <- quantile(CR[,i], 0.99)
            CR[CR[,i] > q,i] <- q
        }
        CR <- CR * A
        NPIX[spp,1:ncol(CR)] <- colSums(CR)
    }
    save(NPIX, file="~/GoogleWork/bam/PIF-AB/results-v3/NPIX.RData")
}
## PIX pop sizes: 2006-2015 subset
if (FALSE) {
    load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
    k <- kgrid[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",]
    A <- k$Area_km2 - k$Area_km2*k$pWater

    SPP <- gsub(".RData", "", list.files(file.path(ROOT, "out", "north2")))

    NPIX2 <- matrix(NA, length(SPP), 240)
    rownames(NPIX2) <- SPP

    for (spp in SPP) {
        cat(spp, "\n")
        flush.console()
        e <- new.env()
        load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6-subset/", spp, ".RData"), envir=e)
        CR <- e$CR[rownames(k),]
        for (i in 1:ncol(CR)) {
            q <- quantile(CR[,i], 0.99)
            CR[CR[,i] > q,i] <- q
        }
        CR <- CR * A
        NPIX2[spp,1:ncol(CR)] <- colSums(CR)
    }
    save(NPIX2, file="~/GoogleWork/bam/PIF-AB/results-v3/NPIX-subset-2006-2015.RData")

    ## compare the two
    library(intrval)
    load("~/GoogleWork/bam/PIF-AB/results-v3/NPIX.RData") # NPIX
    load("~/GoogleWork/bam/PIF-AB/results-v3/NPIX-subset-2006-2015.RData") # NPIX
    NPIX <- NPIX[rownames(NPIX2),]

    tab <- t(apply(NPIX, 1, quantile, c(0.5, 0.025, 0.975), na.rm=TRUE))
    tab2 <- t(apply(NPIX2, 1, quantile, c(0.5, 0.025, 0.975), na.rm=TRUE))
    table(tab[,2:3] %[o]% tab2[,2:3])
    which(!(tab[,2:3] %[o]% tab2[,2:3]))
    data.frame(tab, tab2)["WTSP",]

    op <- par(mfrow=c(1,3))
    plot(tab[,1], tab2[,1], xlab="All (1993-2017)", ylab="Subset (2006-2015)",
        xlim=c(0, max(tab[,1])), ylim=c(0, max(tab[,1])), main=expression(N[PIX]))
    abline(0,1, lty=2)
    segments(x0=tab[,1], y0=tab2[,2], y1=tab2[,3], col="grey")
    segments(y0=tab2[,1], x0=tab[,2], x1=tab[,3], col="grey")
    points(tab[,1], tab2[,1], pch=19)
    mm <- lm(tab2[rownames(tab) != "PISI",1] ~ tab[rownames(tab) != "PISI",1])
    abline(mm, col=2)
    summary(mm)

    hist((tab2[,1]/tab[,1])[rownames(tab) != "PISI"], xlab="subset/all", main="Ratio")

    plot(tab[,3]-tab[,2], tab2[,3]-tab2[,2], xlim=c(0, 10^7), ylim=c(0, 10^7),
        xlab="All (1993-2017)", ylab="Subset (2006-2015)", main="CI range")
    abline(0,1)
    par(op)


}

load("~/GoogleWork/bam/PIF-AB/results-v3/NPIX.RData") # NPIX
tab <- t(apply(NPIX, 1, quantile, c(0.5, 0.025, 0.975), na.rm=TRUE))
TAB$Npix <- tab[rownames(TAB),1] * TAB$PairAdj / 10^6
TAB$Npix95lower <- tab[rownames(TAB),2] * TAB$PairAdj / 10^6
TAB$Npix95upper <- tab[rownames(TAB),3] * TAB$PairAdj / 10^6

## offset components
if (FALSE) {

    library(maptools)
    library(intrval)
    source("~/repos/abmianalytics/birds/00-functions.R")

    dd <- en$DAT
    dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
    dd$start <- as.POSIXlt(dd$DATI)$hour + as.POSIXlt(dd$DATI)$min / 60
    keep <- is.na(dd$DATI) # these will be constant phi
    keep[dd$JULIAN %[]% c(125, 200)] <- TRUE
    keep[dd$start %[]% c(3, 12)] <- TRUE
    dd <- droplevels(dd[keep,])
    dd$JDAY <- dd$JULIAN / 365
    Coor <- as.matrix(dd[,c("X", "Y")])
    JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
    subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
    sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
    dd$srise <- NA
    dd$srise[subset] <- sr
    dd$TSSR <- (dd$start - dd$srise) / 24
    dd$JDAY2 <- dd$JDAY^2
    dd$TSSR2 <- dd$TSSR^2

    ee <- new.env()
    load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData", envir=ee)
    vc1 <- ee$vc1
    tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
    rownames(tv) <- tv[,1]
    lcc4 <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"LCC4"]))
    tmp <- find_max(lcc4)
    dd$LCC4 <- tmp$index
    dd$LCC2 <- dd$LCC4
    levels(dd$LCC2) <- c("OpenWet", "Forest", "Forest", "OpenWet")
    dd$pClosed <- rowSums(row_std(vc1[rownames(dd),])[,tv[colnames(vc1), "is_closed"]])

    Xp <- cbind("(Intercept)"=1, as.matrix(dd[,c("TSSR","JDAY","TSSR2","JDAY2")]))
    Xq <- cbind("(Intercept)"=1, TREE=dd$pClosed,
        LCC2OpenWet=ifelse(dd$LCC2=="OpenWet", 1, 0),
        LCC4Conif=ifelse(dd$LCC4=="Conif", 1, 0),
        LCC4Open=ifelse(dd$LCC4=="Open", 1, 0),
        LCC4Wet=ifelse(dd$LCC4=="Wet", 1, 0))
    summary(Xp)
    summary(Xq)

    library(QPAD)
    load_BAM_QPAD(version=3)
    sppp <- SPP

    off <- matrix(NA, nrow(dd), length(sppp))
    rownames(off) <- rownames(dd)
    colnames(off) <- sppp

    dd0 <- dd
    TAU <- matrix(NA, length(SPP), 240)
    rownames(TAU) <- SPP
    PHI <- TAU

    for (spp in sppp) {
        ## print out where we are
        cat(spp, "\n");flush.console()
        for (i in 1:240) {
            ss <- en$BB[,i]
            dd <- dd0[ss,]
            p <- rep(NA, nrow(dd))
            A <- q <- p
            cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
            mi <- bestmodelBAMspecies(spp, type="BIC",
                model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
            cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
            Xp2 <- Xp[ss,names(cfi$sra),drop=FALSE]
            OKp <- rowSums(is.na(Xp2)) == 0
            Xq2 <- Xq[ss,names(cfi$edr),drop=FALSE]
            OKq <- rowSums(is.na(Xq2)) == 0
            p[!OKp] <- sra_fun(dd$MAXDUR[!OKp], cf0[1])
            unlim <- ifelse(dd$MAXDIS[!OKq] == Inf, TRUE, FALSE)
            A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dd$MAXDIS[!OKq]^2)
            q[!OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[!OKq], cf0[2]))
            phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
            tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
            PHI[spp,i] <- mean(phi1)
            TAU[spp,i] <- mean(tau1) * 100
        }
    }
    save(PHI, TAU, file="~/GoogleWork/bam/PIF-AB/results-v3/PHIandTAU.RData")
}
load("~/GoogleWork/bam/PIF-AB/results-v3/PHIandTAU.RData") # PHI, TAU
tmp1 <- t(apply(TAU, 1, quantile, c(0.5, 0.025, 0.975)))
tmp2 <- t(apply(PHI, 1, quantile, c(0.5, 0.025, 0.975)))
TAB$EDR <- tmp1[,1]
#TAB$tau95lower <- tmp1[,2]
#TAB$tau95upper <- tmp1[,3]
TAB$phi <- tmp2[,1]
TAB$p3 <- 1-exp(-3 * TAB$phi)
#TAB$phi95lower <- tmp2[,2]
#TAB$phi95upper <- tmp2[,3]
rm(tmp1, tmp2)

## roadside count effects
if (FALSE) {

    ROAD0 <- matrix(NA, length(SPP), 240)
    rownames(ROAD0) <- SPP
    ROAD1 <- ROAD0
    #bb <- unique(as.numeric(en$BB))
    #ddr <- en$DAT[bb,]
    #ddr <- ddr[ddr$ROAD==1,]
    ddr <- en$DAT
    ddr <- ddr[ddr$ROAD==1,]

    Xn1 <- get_model_matrix(ddr, en$mods) # has ROAD and mTrSft as >0
    Xn0 <- Xn1
    Xn0[,"ROAD"] <- 0
    Xn0[,"mTrSft"] <- 0
    summary(Xn0[,c("ROAD","mTrSft")])
    summary(Xn1[,c("ROAD","mTrSft")])

    for (spp in SPP) {
        cat(spp, "N\n");flush.console()
        resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))

        lam0 <- exp(predict_with_SSH(resn, Xn0, en$SSH, stage="Space"))
        lam1 <- exp(predict_with_SSH(resn, Xn1, en$SSH, stage="Space"))
        ROAD0[spp,1:ncol(lam0)] <- colMeans(lam0)
        ROAD1[spp,1:ncol(lam1)] <- colMeans(lam1)
    }
    save(ROAD0, ROAD1, file="~/GoogleWork/bam/PIF-AB/results-v3/ROAD.RData")
}
load("~/GoogleWork/bam/PIF-AB/results-v3/ROAD.RData") # ROAD0, ROAD1
tmp0 <- t(apply(ROAD0, 1, quantile, c(0.5, 0.025, 0.975), na.rm=TRUE))
tmp1 <- t(apply(ROAD1, 1, quantile, c(0.5, 0.025, 0.975), na.rm=TRUE))
TAB$Y0 <- tmp0[,1]
TAB$Y0_95lower <- tmp0[,2]
TAB$Y0_95upper <- tmp0[,3]
TAB$Y1 <- tmp1[,1]
TAB$Y1_95lower <- tmp1[,2]
TAB$Y1_95upper <- tmp1[,3]
#TAB$Y1[TAB$Y1 < 0.0001] <- 0.0001 # impacts 1 species, CSWA
TAB$Y0_Y1 <- sign(TAB$Y0-TAB$Y1)
TAB$Y0_Y1[tmp0[,2:3] %[o]% tmp1[,2:3]] <- 0
table(TAB$Y0_Y1)
#TAB[TAB$Y0_Y1 > 0,]
#TAB[TAB$Y0_Y1 < 0,]


## habitat representation stuff
if (FALSE) {

    load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
    k <- kgrid[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",]
    A <- k$Area_km2 - k$Area_km2*k$pWater
    ## ch2soil ch2veg trSoil trVeg
    load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
    stopifnot(all(rownames(kgrid) == rownames(trVeg)))
    stopifnot(all(rownames(kgrid) == rownames(trSoil)))

    tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
    rownames(tv) <- tv[,1]
    tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])

    ch2veg$g <- tv$UseInPIX[match(ch2veg$cr, rownames(tv))]
    ch2veg$f <- ifelse(startsWith(as.character(ch2veg$cr), "CC"), "CC", "not")
    ch2veg$s <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

    rn <- c(
        "Decid", "DecidO",
        "Mixed", "MixedO",
        "Pine", "PineO",
        "WSpr", "WSprO",
        "BSpr", "BSprO",
        "TrFen",
        "Shrub",
        "Grass",
        "GrFen",
        "Marsh",
        "Swamp",
        "SoftLin",
        "Roads",
        "Agr",
        "UrbInd")
    AAx <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "f"])
    AAx <- colSums(AAx[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",])/10^6
    AA <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "g"])
    AA <- colSums(AA[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",])
    AA <- AA[rn]/10^6 # in km^2


    if (FALSE) {
        ## study area
        AAx <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "s"])
        AAx <- colSums(AAx[kgrid$NRNAME != "Grassland",])/10^6
        AA <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "g"])
        AA <- colSums(AA[kgrid$NRNAME != "Grassland",])
        AA <- AA[rn]/10^6 # in km^2
        sum(kgrid$NRNAME != "Grassland")
        data.frame(p=round(100*AAx/sum(AAx),1))
        round(100*AA/sum(AA),1)
    }


    HABITAT <- list()

    for (spp in SPP) {
        cat(spp, "\n")
        flush.console()
        e <- new.env()
        load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6/", spp, ".RData"), envir=e)
        HAB <- e$HAB
        #compare_sets(tv$UseInAnalysisFineAge, rownames(HAB))
        g <- tv$UseInPIX[match(rownames(HAB), tv$UseInAnalysisFineAge)]
        NN <- groupSums(HAB, 1, g)[rn,]
        NN <- NN * TAB[spp, "PairAdj"] # inds in kn2 cell
        ## diff is due to weighting by 1-pWater for pop size, which is the right thing to do
        d <- TAB[spp,"Npix"] / (median(colSums(NN))/10^6)
        NN <- NN * d # adjust to same median
        DD <- (NN / 100) / AA # inds/ha
        t(apply(NN, 1, quantile, c(0.5, 0.025, 0.975)))
        t(apply(DD, 1, quantile, c(0.5, 0.025, 0.975)))

        HABITAT[[spp]] <- list(N=NN, A=AA, D=DD)
    }

    vc <- as.character(en$DAT$vegc)
    table(en$DAT$vegc)
    summary(en$DAT$wtAge*200)
    vc[vc == "Decid" & en$DAT$wtAge*200 > 60] <- "DecidO"
    vc[vc == "Mixedwood" & en$DAT$wtAge*200 > 60] <- "MixedwoodO"
    vc[vc == "Spruce" & en$DAT$wtAge*200 > 80] <- "SpruceO"
    vc[vc == "Pine" & en$DAT$wtAge*200 > 80] <- "PineO"
    vc[vc == "BSpr" & en$DAT$wtAge*200 > 80] <- "BSprO"

    nn <- c(
        "Decid"="Decid",
        "DecidO"="DecidO",
        "MixedwoodO"="MixedO",
        "Mixedwood"="Mixed",
        "Spruce"="WSpr",
        "Pine"="Pine",
        "SpruceO"="WSprO",
        "Swamp"="Swamp",
        "Shrub"="Shrub",
        "BSpr"="BSpr",
        "RoughP"="Agr",
        "Crop"="Agr",
        "PineO"="PineO",
        "Industrial"="UrbInd",
        "BSprO"="BSprO",
        "Rural"="UrbInd",
        "GrassHerb"="Grass",
        "GraminoidFen"="GrFen",
        "Marsh"="Marsh",
        "Larch"="TrFen",
        "Mine"="UrbInd",
        "TameP"="Agr",
        "Urban"="UrbInd")
    for (i in names(nn))
        vc[vc == i] <- nn[i]
    data.frame(table(vc)[rn])
    ndet <- table(vc)
    ndet <- structure(as.numeric(ndet[rn]), names=rn)

    ## determine BBS in BCR6
    r <- cure4insect::.read_raster_template()
    od <- setwd("d:/spatial/bcr/")
    BCR <- readOGR(".", "BCR_Terrestrial_master") # rgdal
    BCR <- spTransform(BCR, proj4string(r))
    setwd(od)
    xypt <- en$DAT
    coordinates(xypt) <- ~ X + Y
    proj4string(xypt) <-
        CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    xypt <- spTransform(xypt, proj4string(r))
    xypt2BCR <- over(xypt, BCR)

    ss <- startsWith(rownames(en$DAT), "BBS") & !duplicated(en$DAT$SS) & xypt2BCR$BCR==6
    vc2 <- vc[ss]
    wroad <- table(vc2)/sum(table(vc2))
    wroad <- structure(as.numeric(wroad[rn]), names=rn)
    avail <- AA/sum(AA)

    save(HABITAT, avail, wroad, ndet, file="~/GoogleWork/bam/PIF-AB/results-v3/HABITAT.RData")

}
## habitat representation stuff: 2006-2015 subset
if (FALSE) {

    load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
    k <- kgrid[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",]
    A <- k$Area_km2 - k$Area_km2*k$pWater
    ## ch2soil ch2veg trSoil trVeg
    load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
    stopifnot(all(rownames(kgrid) == rownames(trVeg)))
    stopifnot(all(rownames(kgrid) == rownames(trSoil)))

    tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
    rownames(tv) <- tv[,1]
    tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])

    ch2veg$g <- tv$UseInPIX[match(ch2veg$cr, rownames(tv))]
    ch2veg$f <- ifelse(startsWith(as.character(ch2veg$cr), "CC"), "CC", "not")
    ch2veg$s <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

    rn <- c(
        "Decid", "DecidO",
        "Mixed", "MixedO",
        "Pine", "PineO",
        "WSpr", "WSprO",
        "BSpr", "BSprO",
        "TrFen",
        "Shrub",
        "Grass",
        "GrFen",
        "Marsh",
        "Swamp",
        "SoftLin",
        "Roads",
        "Agr",
        "UrbInd")
    AAx <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "f"])
    AAx <- colSums(AAx[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",])/10^6
    AA <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "g"])
    AA <- colSums(AA[kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS",])
    AA <- AA[rn]/10^6 # in km^2


    if (FALSE) {
        ## study area
        AAx <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "s"])
        AAx <- colSums(AAx[kgrid$NRNAME != "Grassland",])/10^6
        AA <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "g"])
        AA <- colSums(AA[kgrid$NRNAME != "Grassland",])
        AA <- AA[rn]/10^6 # in km^2
        sum(kgrid$NRNAME != "Grassland")
        data.frame(p=round(100*AAx/sum(AAx),1))
        round(100*AA/sum(AA),1)
    }

    SPP <- gsub(".RData", "", list.files(file.path(ROOT, "out", "north2")))

    HABITAT2 <- list()

    for (spp in SPP) {
        gc()
        cat(spp, "\n")
        flush.console()
        e <- new.env()
        load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/bcr6-subset/", spp, ".RData"), envir=e)
        HAB <- e$HAB
        #compare_sets(tv$UseInAnalysisFineAge, rownames(HAB))
        g <- tv$UseInPIX[match(rownames(HAB), tv$UseInAnalysisFineAge)]
        NN <- groupSums(HAB, 1, g)[rn,]
        NN <- NN * TAB[spp, "PairAdj"] # inds in kn2 cell
        ## diff is due to weighting by 1-pWater for pop size, which is the right thing to do
        #d <- TAB[spp,"Npix"] / (median(colSums(NN))/10^6)
        #NN <- NN * d # adjust to same median
        DD <- (NN / 100) / AA # inds/ha
        t(apply(NN, 1, quantile, c(0.5, 0.025, 0.975)))
        t(apply(DD, 1, quantile, c(0.5, 0.025, 0.975)))

        HABITAT2[[spp]] <- list(N=NN, A=AA, D=DD)
    }

    save(HABITAT2, file="~/GoogleWork/bam/PIF-AB/results-v3/HABITAT-subset-2006-2015.RData")

    load("~/GoogleWork/bam/PIF-AB/results-v3/HABITAT.RData")
    load("~/GoogleWork/bam/PIF-AB/results-v3/HABITAT-subset-2006-2015.RData")

    fun <- function(x) {
        x <- t(apply(x$D, 1, quantile, c(0.5, 0.025, 0.975)))
        colnames(x) <- c("est", "lwr", "upr")
        x
    }
    Dsub <- lapply(HABITAT2, fun)
    Dall <- lapply(HABITAT, fun)[names(Dsub)]
    library(intrval)

    pdf("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/mapsSubset/habitat.pdf", onefile=TRUE)
    #spp <- "ALFL"
    for (spp in names(Dsub)) {
        z <- data.frame(Hab=rownames(Dall[[1]]), All=Dall[[spp]], Sub=Dsub[[spp]])
        zz <- list(z$All.lwr, z$All.upr) %)o(% list(z$Sub.lwr, z$Sub.upr)
        lim <- c(0, min(3*max(z$All.est, z$Sub.est), quantile(z[,-1], 0.99, na.rm=TRUE)))
        plot(z$All.est, z$Sub.est, main=spp, xlab="D (1993-2017)", ylab="D (2006-2015)",
            xlim=lim, ylim=lim, pch=19, col="grey")
        segments(x0=z$All.est, y0=z$Sub.lwr, y1=z$Sub.upr, col="grey")
        segments(y0=z$Sub.est, x0=z$All.lwr, x1=z$All.upr, col="grey")
        abline(0, 1, lty=2)
        lines(lowess(z$Sub.est ~ z$All.est), col=4)
        if (any(zz)) {
            points(z$All.est[zz], z$Sub.est[zz], pch=19, col=1)
            segments(x0=z$All.est[zz], y0=z$Sub.lwr[zz], y1=z$Sub.upr[zz], col="1")
            segments(y0=z$Sub.est[zz], x0=z$All.lwr[zz], x1=z$All.upr[zz], col="1")
            text(z$All.est[zz], z$Sub.est[zz], z$Hab[zz], col=2, cex=0.8)
        }
    }
    dev.off()

}
load("~/GoogleWork/bam/PIF-AB/results-v3/HABITAT.RData") # HABITAT, avail, wroad

h_fun <- function(x) {
    rn <- names(avail)[!(names(avail) %in% c("SoftLin","Roads"))]
    DD <- apply(x$D, 1, median)
    sum(DD[rn] * wroad[rn]) / sum(DD[rn] * avail[rn])
}
TAB$H <- sapply(HABITAT[rownames(TAB)], h_fun)

## table 1
tab1 <- data.frame(a=100*avail,w=100*wroad, wa=wroad/avail,n=ndet)
tab1 <- tab1[!is.na(tab1$w),]
tab1 <- round(tab1, 1)
colSums(tab1)
tab1$a <- 100*tab1$a/sum(tab1$a)
tab1$w <- 100*tab1$w/sum(tab1$w)
colSums(round(tab1,1))
round(tab1,1)


## Final subset -------------------------

#Issues:
#
#DUFL - edge species, very extreme CI range, is there an adequate sample in Sb bogs and fens?
#EUST - parkland species, extreme CI range
#HOSP - parkland species, very extreme CI range
#MAWR - extreme prediction CI, extreme phi
#MOBL - edge species, very extreme CI range
#PIGR - very rare, very extreme CI range
#SPGR - very rare, very extreme CI range, extreme roadside estimate, no PIF estimate
#RUGR - no PIF estimate
#TOWA - edge & rare species, very extreme CI range
#
#Probably OK:
#
#ROPI - southern species, rare?
#GRCA - too rare? But modelled fairly well in habitat (better than some)
#RECR - too rare? Habitat model is peculiarly specific.
#       We should get multiple types here, including both a pine and a spruce specialist.
#       The emphasis on very young pine, especially very young burnt pine, seems weird,
#       and huge CI's suggest this could be spurious.
#WBNU - too rare?
#VATH - mountain species, quite large CI range, but would not say extreme
#WCSP - rare, quite large CI range, bordering extreme
#PISI - high D in foothills when taking the mean, not with median
#
#Other notes:
#
#In some species (e.g., ALFL, AMRO, HAWO) you can see discontinuity in predicted density where there
#    is a boundary between veg data from WBNP and from Phase 1 inventories
#    (northern edge, about 40% way from left). Suggests that veg is not being classified in the
#    same way in these systems. Possibly similar issue in NW corner with
#    AVIE vs Phase 1 veg data for NOWA.
#EUST variation looks better than EAPH in habitat selection figure
#Roadside bias seems rather extreme for BLPW, BTNW, CSWA, EAPH, ROPI, SPGR
#EDR seems extreme for MAWR, VESP
#Extreme CI's for BARS, CSWA, DUFL, DEJU?, EUST, HOSP, MAWR, MOBL, OSFL?,
#    PIGR, PISI, PUFI, ROPI, SPGR, TOWA, VATH, WCSP

summary(TAB)
TAB$R <- TAB$Y0/TAB$Y1
head(TAB[TAB$R > 15,])

TAB[is.na(TAB$Npif),]

#bb <- unique(as.numeric(en$BB))
#TAB$nDet <- sapply(SPP, function(spp) sum(en$YY[bb,spp]>0))


## Deltas ---
TAB$DeltaObs <- log(TAB$Npix / TAB$Npif)
TAB$DeltaR <- log(TAB$Y0/TAB$Y1)
TAB$DeltaT <- log((1/TAB$p3)/TAB$TimeAdj)
TAB$DeltaA <- log((1/TAB$EDR^2) / (1/TAB$MDD^2))
TAB$DeltaH <- log(1/TAB$H) # we take inverse because H=1 is the PIF setup
TAB$DeltaExp <- TAB$DeltaR + TAB$DeltaT + TAB$DeltaA + TAB$DeltaH
TAB$epsilon <- TAB$DeltaObs - TAB$DeltaExp

a <- read.csv("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/NACC_list_species.csv")
levels(TAB$CommonName)[levels(TAB$CommonName)=="Le Conte's Sparrow"] <- "LeConte's Sparrow"
levels(TAB$CommonName)[levels(TAB$CommonName)=="Gray Jay"] <- "Canada Jay"
compare_sets(a$common_name, TAB$CommonName)
setdiff(TAB$CommonName, a$common_name)
a$sort <- 1:nrow(a)
a <- droplevels(a[match(TAB$CommonName, a$common_name),])
TAB$ScientificName <- a$species
TAB$Order <- a$order
TAB$Family <- a$family
TAB$sort <- a$sort
TAB <- TAB[order(TAB$sort),]

TAB$Npix_Npif <- sign(TAB$Npix-TAB$Npif)
TAB$Npix_Npif[TAB[,c("Npix95lower","Npix95upper")] %[o]% TAB[,c("Npif95lower","Npif95upper")]] <- 0
table(TAB$Npix_Npif)

write.csv(TAB, row.names=FALSE,
    file=paste0("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/pifpix-v3-all-results.csv"))

## comparing v2 to v3 ---------------------------------

TAB <- read.csv("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/pifpix-v3-all-results.csv")
rownames(TAB) <- TAB$Code

pop <- read.csv("~/GoogleWork/bam/PIF-AB/draft2/Table1-estimates.csv")
rownames(pop) <- pop[,1]
sp <- intersect(rownames(TAB), rownames(pop))

op <- par(mfrow=c(2,3))
plot(pop[sp,"Npix"], TAB[sp,"Npix"], ylab="New", xlab="Old", main="Npix");abline(0,1)
plot(pop[sp,"p3"], TAB[sp,"p3"], ylab="New", xlab="Old", main="p3");abline(0,1)
plot(pop[sp,"EDR"], TAB[sp,"EDR"], ylab="New", xlab="Old", main="EDR");abline(0,1)
plot(pop[sp,"Npif"], TAB[sp,"Npif"], ylab="New", xlab="Old", main="Npif");abline(0,1)
plot(pop[sp,"TimeAdj"], TAB[sp,"TimeAdj"], ylab="New", xlab="Old", main="TimeAdj");abline(0,1)
plot(pop[sp,"MDD"], TAB[sp,"MDD"], ylab="New", xlab="Old", main="MDD");abline(0,1)
par(op)


op <- par(mfrow=c(2,3))
plot(pop[sp,"DeltaObs"], TAB[sp,"DeltaObs"], ylab="New", xlab="Old", main="OBS");abline(0,1)
plot(pop[sp,"DeltaExp"], TAB[sp,"DeltaExp"], ylab="New", xlab="Old", main="EXP");abline(0,1)
plot(pop[sp,"DeltaT"], TAB[sp,"DeltaT"], ylab="New", xlab="Old", main="T");abline(0,1)
plot(pop[sp,"DeltaA"], TAB[sp,"DeltaA"], ylab="New", xlab="Old", main="A");abline(0,1)
plot(pop[sp,"DeltaR"], TAB[sp,"DeltaR"], ylab="New", xlab="Old", main="R");abline(0,1)
plot(pop[sp,"DeltaH"], TAB[sp,"DeltaH"], ylab="New", xlab="Old", main="H");abline(0,1)
par(op)



## ----- figures for v3 -----------------------

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(MASS)
library(KernSmooth)
library(vegan)
#library(gstat)
#library(viridis)
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds"

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")

r <- cure4insect::.read_raster_template()

od <- setwd("d:/spatial/bcr/")
BCR <- readOGR(".", "BCR_Terrestrial_master") # rgdal
BCR <- spTransform(BCR, proj4string(r))
BCR <- gSimplify(BCR, tol=500, topologyPreserve=TRUE)
setwd(od)

od <- setwd("d:/spatial/NSR")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(r))
AB <- gUnaryUnion(AB, rep(1, nrow(AB))) # province
AB <- gSimplify(AB, tol=500, topologyPreserve=TRUE)
setwd(od)

od <- setwd("d:/spatial/ab-lakes")
lak <- readOGR(".", "Lakes_named") # rgdal
lak <- spTransform(lak, proj4string(r))
#lak <- gSimplify(BCR, tol=500, topologyPreserve=TRUE)
setwd(od)
hist(lak@data$Shape_Area)

lak <- lak[lak@data$Shape_Area > 100*10^6,]

BCR2AB <- gIntersection(AB, BCR, byid=TRUE)


## load numbers
pop <- read.csv("d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/pifpix-v3-all-results.csv")
rownames(pop) <- pop$Code
SPP <- rownames(pop)
EXCL <- c("DUFL", "HOSP", "MAWR", "MOBL", "PIGR", "SPGR", "RUGR", "TOWA", "RECR", "CSWA", "EUST", "EAPH", "ROPI")
SPP <- rownames(pop)[!(rownames(pop) %in% EXCL)]
pop <- droplevels(pop[SPP,])
write.csv(pop, row.names=FALSE, file="d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/pifpix-v3-final-results.csv")

## subset years of 2006-2015
if (FALSE) {
    load("~/GoogleWork/bam/PIF-AB/results-v3/NPIX-subset-2006-2015.RData") # NPIX
    tab2 <- t(apply(NPIX2, 1, quantile, c(0.5, 0.025, 0.975), na.rm=TRUE))
    pop$Npix <- tab2[rownames(pop),1] * pop$PairAdj / 10^6
    pop$Npix95lower <- tab2[rownames(pop),2] * pop$PairAdj / 10^6
    pop$Npix95upper <- tab2[rownames(pop),3] * pop$PairAdj / 10^6
    pop <- droplevels(pop[rownames(pop) != "PISI",])
}

load("~/GoogleWork/bam/PIF-AB/results-v3/HABITAT.RData") # HABITAT, avail, wroad
Dall <- data.frame(Ahab=100*avail, Whab=100*wroad, sapply(HABITAT[SPP], function(z) apply(z$D, 1, median)))
Dall <- Dall[!(rownames(Dall) %in% c("SoftLin","Roads")),]
sum(is.na(Dall))
Dall$Ahab <- 100*Dall$Ahab/sum(Dall$Ahab)
Dall$Whab <- 100*Dall$Whab/sum(Dall$Whab)
options(scipen=999)
write.csv(Dall, file="d:/abmi/AB_data_v2018/data/analysis/birds/bcr6/pifpix-v3-densities.csv")


## maps

xnss <- nonDuplicated(en$DAT, SS, TRUE)
xy <- xnss
coordinates(xy) <- ~ X + Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, proj4string(r))
xy2BCR <- over(xy, BCR)

xypt <- en$DAT
coordinates(xypt) <- ~ X + Y
proj4string(xypt) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xypt <- spTransform(xypt, proj4string(r))
xypt2BCR <- over(xypt, BCR)


pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig1-maps.pdf", width=12, height=9)
op <- par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(BCR2AB, col=c(NA, "grey", rep(NA, 11)), border=NA, main="Roadside surveys")
#plot(lak, col="white", border=NA,add=TRUE)
plot(AB, col=NA, border=1,add=TRUE)
plot(xy[xy@data$ROAD == 1,], add=TRUE, pch=19, col=1, cex=0.25)

plot(BCR2AB, col=c(NA, "grey", rep(NA, 11)), border=NA, main="Off-road surveys")
#plot(lak, col="white", border=NA,add=TRUE)
plot(AB, col=NA, border=1,add=TRUE)
plot(xy[xy@data$ROAD == 0,], add=TRUE, pch=19, col=1, cex=0.25)
par(op)
dev.off()


## pop size

pch <- 19 # ifelse(pop$Npif %[]% list(pop$NpixLo, pop$NpixHi), 21, 19)
Min <- 3*2
Siz <- 24*2
Up <- 75
pch2 <- 19 # ifelse(ppp$Npif %[]% list(ppp$NpixLo, ppp$NpixHi), 21, 19)

pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig2-popsize.pdf", width=7, height=7)
op <- par(mfrow=c(1,1), las=1, mar=c(4,4,1,2))
plot(pop[,c("Npif", "Npix")], xlab=expression(N[PIF]), ylab=expression(N[PIX]),
    type="n",
    pch=pch, col="#00000080", xlim=c(0, max(pop[,c("Npif", "Npix")])),
    ylim=c(0, max(pop[,c("Npif", "Npix")])))
abline(0,1, lty=2)
di <- sqrt(pop[,"Npif"]^2+pop[,"Npix"]^2) > Min*1.2
pp <- pop[pop[,"Npif"] < Min & pop[,"Npix"] < Min, c("Npif", "Npix")]
ppp <- pop[pop[,"Npif"] < Min & pop[,"Npix"] < Min, ]
di2 <- sqrt(pp[,"Npif"]^2+pp[,"Npix"]^2) > 1
pp <- pp*Siz/Min
pp[,1] <- pp[,1] + (Up-Siz)
lines(c(0,Up-Siz), c(0, 0), col="grey")
lines(c(0,Up-Siz), c(Min, Siz), col="grey")
rect(0, 0, Min, Min)
rect(Up-Siz, 0, Up-Siz+Min*Siz/Min, Min*Siz/Min, col="white")
lines(c(Up-Siz, Up), c(0, Siz), col=1, lty=2)
points(pp, pch=pch2, col="#00000080")
points(pop[,c("Npif", "Npix")], pch=pch, col="#00000080")
text(pop[,"Npif"]+1.2*2, pop[,"Npix"]+0, labels=ifelse(di, rownames(pop), ""), cex=0.5)
text(pp[,"Npif"]+1.2*2, pp[,"Npix"]+0, labels=ifelse(di2, rownames(pp), ""), cex=0.5)
segments(x0=c(Up-Siz+c(0, 1, 2, 3)*Siz/3), y0=rep(Siz, 4), y1=rep(Siz, 4)+0.5*2)
text(c(Up-Siz+c(0, 1, 2, 3)*Siz/3), rep(Siz, 4)+1*2.5, 2*(0:3))
dev.off()

## pop rank

pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig3-poprank.pdf", width=7, height=7)
op <- par(mfrow=c(1,1), las=1, mar=c(4,4,1,2))
rnk <- cbind(rank(pop$Npif), rank(pop$Npix))
#rnk <- 100 * rnk / max(rnk)
plot(rnk,
    #xlab="PIF rank quantile (%)", ylab="PIX rank quantile (%)",
    xlab=expression(N[PIF]~rank), ylab=expression(N[PIX]~rank),
    #xlim=c(0,100), ylim=c(0,100),
    type="n", axes=FALSE)
abline(0,1, lty=2)
text(rnk, labels=rownames(pop), cex=0.8)
axis(1, c(1, 20, 40, 60, 80, 95), c(1, 20, 40, 60, 80, 95))
axis(2, c(1, 20, 40, 60, 80, 95), c(1, 20, 40, 60, 80, 95))
box()
#abline(h=c(24,48,72),v=c(24,48,72))
par(op)
dev.off()


## components

dots_box_plot <- function(mat, lines=FALSE, method="box", ...) {
    set.seed(1)
    rnd <- runif(nrow(mat), -0.05, 0.05)
    boxplot(mat, range=0, border="white",...)
    if (lines)
        for (i in 2:ncol(mat))
            segments(x0=i+rnd-1, x1=i+rnd, y0=mat[,i-1], y1=mat[,i], col="lightgrey")
    for (i in 1:ncol(mat))
        points(i+rnd, mat[,i], pch=19, col="#00000080")

    if (method == "box") {
        #boxplot(mat, range=0, add=TRUE, col="#ff000020", names=NA)
        boxplot(mat, range=0, add=TRUE, col="#00000020", names=NA)
    } else {
        v <- 0.2
        for (i in 1:ncol(mat)) {
            xx <- sort(mat[,i])
            st <- boxplot.stats(xx)
            s <- st$stats
            if (method == "kde")
                d <- bkde(xx) # uses Normal kernel
            if (method == "fft")
                d <- density(xx) # uses FFT
            if (method == "hist") {
                h <- hist(xx, plot=FALSE)
                xv <- rep(h$breaks, each=2)
                yv <- c(0, rep(h$density, each=2), 0)
            } else {
                xv <- d$x
                yv <- d$y
                jj <- xv >= min(xx) & xv <= max(xx)
                xv <- xv[jj]
                yv <- yv[jj]
            }
            yv <- 0.4 * yv / max(yv)
            polygon(c(-yv, rev(yv))+i, c(xv, rev(xv)), col="#00000020", border="#40404080")
            polygon(c(-v,-v,v,v)+i, s[c(2,4,4,2)], col="#40404080", border=NA)
            lines(c(-v,v)+i, s[c(3,3)], lwd=2, col="#00000020")
        }
    }

    invisible(NULL)
}


pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig4-components.pdf", width=10, height=7)
op <- par(las=1)
#mat <- pop[,c("DeltaObs", "DeltaExp", "DeltaR", "DeltaT", "DeltaA", "DeltaH")]
#colnames(mat) <- c("OBS", "EXP", "R", "T", "A", "H")
mat <- pop[,c("DeltaObs", "DeltaExp", "DeltaT", "DeltaA", "DeltaR", "DeltaH")]
colnames(mat) <- c("OBS", "EXP", "T", "A", "R", "H")
par(las=1)
dots_box_plot(mat, ylab="Log Ratio", method="kde")
abline(h=0, col=1, lwd=1,lty=2)

off <- 0.2
i <- 1
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

i <- 2
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

i <- 3
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

i <- 4
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

i <- 5
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 2), cex=0.8, rownames(z)[1:2])
text(i+off, tail(z[,i], 2), cex=0.8, rownames(z)[(nrow(z)-1):nrow(z)])

i <- 6
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

for (i in 1:6) {
    zz1 <- format(mean(mat[,i]), trim = TRUE, scientific = FALSE, digits = 2)
    zz2 <- format(sd(mat[,i]), trim = TRUE, scientific = FALSE, digits = 2)
    mtext(paste("Mean =", zz1), side=1,at=i,line=3, cex=0.8)
    mtext(paste("SD =", zz2), side=1,at=i,line=4, cex=0.8)
}
par(op)
dev.off()

## components based regression: table 2

mod <- lm(DeltaObs ~ DeltaR + DeltaT + DeltaA + DeltaH, pop)
an <- anova(mod)
an$Percent <- 100 * an[["Sum Sq"]] / sum(an[["Sum Sq"]])
an <- an[c("Df", "Sum Sq", "Percent", "Mean Sq", "F value", "Pr(>F)")]
an
summary(mod)
summary(mod)$sigma^2
zval <- c(coef(summary(mod))[,1] - c(0, 1, 1, 1, 1))/coef(summary(mod))[,2]
round(cbind(coef(summary(mod))[,1:2], ph=2 * pnorm(-abs(zval))), 3)
round(an$Percent,1)

mod2 <- step(lm(DeltaObs ~ (DeltaR + DeltaT + DeltaA + DeltaH)^2, pop), trace=0)
summary(mod2)
an2 <- anova(mod2)
an2$Percent <- 100 * an2[["Sum Sq"]] / sum(an2[["Sum Sq"]])
an2 <- an2[c("Df", "Sum Sq", "Percent", "Mean Sq", "F value", "Pr(>F)")]
an2

## looking at shared variation
vpfun <- function (x, cutoff = 0, digits = 1, Xnames, showNote=FALSE, ...)
{
    x <- x$part
    vals <- x$indfract[, 3]
    is.na(vals) <- vals < cutoff
    if (cutoff >= 0)
        vals <- round(vals, digits + 1)
    labs <- format(vals, digits = digits, nsmall = digits + 1)
    labs <- gsub("NA", "", labs)
    showvarparts(x$nsets, labs, Xnames=Xnames, ...)
    if (any(is.na(vals)) && showNote)
        mtext(paste("Values <", cutoff, " not shown", sep = ""),
            1)
    invisible()
}

prt <- vegan::varpart(Y=pop$DeltaObs, ~DeltaR, ~DeltaT, ~DeltaA, ~DeltaH, data=pop)
pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig6-varpart.pdf", width=6, height=6)
vpfun(prt, cutoff = 0, digits = 2, Xnames=c("R", "T", "A", "H"), bg=1)
dev.off()


## road avoidance and ordination
DD <- as.matrix(t(Dall[,-(1:2)]))
NN <- t(t(DD) * Dall$Ahab)
NN <- t(NN / rowSums(NN))
Cex <- pop$DeltaObs
names(Cex) <- rownames(pop)
br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
#c00 <- c('#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba')
c00 <- c('#d7191c','darkgrey','#2b83ba')
#c00 <- c("red", "darkgrey", "blue")
Col0 <- colorRampPalette(c00)(7)
Col <- Col0[cut(Cex, br)]
names(Col) <- names(Cex)
o <- cca(NN)
round(100*eigenvals(o)/sum(eigenvals(o)), 1)
round(cumsum(100*eigenvals(o)/sum(eigenvals(o))), 1)

pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig7-ordination3.pdf", width=9, height=9)
op <- par(las=1)
plot(0, type="n", xlim=c(-0.8,1.2), ylim=c(-1,1), xlab="Axis 1", ylab="Axis 2")
s2 <- scores(o)$sites
s2 <- s2 / max(abs(s2))
for (i in 1:nrow(s2))
    arrows(x0=0,y0=0,x1=s2[i,1],y1=s2[i,2], angle=20, length = 0.1, col="darkgrey")
abline(h=0,v=0,lty=2)
s1 <- scores(o)$species
s1 <- s1 / max(abs(s1))
text(s1[names(Col),]*0.8, labels=colnames(NN),cex=0.75, col=Col)
text(s2*1.05, labels=rownames(NN),cex=1,col=1)
#text(s1[names(Col),]*0.8, labels=colnames(NN),cex=0.75, col=4)
for (ii in 1:200)
    lines(c(0.9, 0.96), rep(seq(-0.55, -0.95, len=200)[ii], 2), col=colorRampPalette(c00)(200)[ii])
text(c(1,1,1)+0.1, c(-0.6, -0.75, -0.9),
    c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])))
par(op)
dev.off()

pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig7-report.pdf", width=8, height=8)
op <- par(las=1)
Colg <- colorRampPalette(c("red","yellow"))(7)[cut(Cex, br)]
names(Colg) <- names(Cex)
plot(0, type="n", xlim=c(-0.8,1.2), ylim=c(-1,1), xlab="", ylab="")
s1 <- scores(o)$species
s1 <- s1 / max(abs(s1))
points(s1[names(Colg),]*0.8, cex=2, col=Colg, pch=19)
points(s1[names(Colg),]*0.8, cex=2, col="orange")
s2 <- scores(o)$sites
s2 <- s2 / max(abs(s2))
for (i in 1:nrow(s2))
    arrows(x0=0,y0=0,x1=s2[i,1],y1=s2[i,2], angle=20, length = 0.1, col="darkgrey")
text(s2*1.05, labels=rownames(NN),cex=1,col=1)
abline(h=0,v=0,lty=2)
#text(s1[names(Col),]*0.8, labels=colnames(NN),cex=0.75, col=4)
for (ii in 1:200)
    lines(c(0.85, 0.9), rep(seq(-0.55, -0.95, len=200)[ii], 2), col=colorRampPalette(c("red","yellow"))(200)[ii])
text(c(1,1,1)+0.1, c(-0.6, -0.75, -0.9),
    c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])))
par(op)
dev.off()


pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig7-ordination-all.pdf", width=9, height=9, onefile=TRUE)
for (ii in c("DeltaObs", "DeltaExp", "DeltaR", "DeltaT", "DeltaA", "DeltaH")) {

    Cex <- pop[[ii]]
    names(Cex) <- rownames(pop)
    br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
    #c00 <- c('#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba')
    c00 <- c('#d7191c','darkgrey','#2b83ba')
    #c00 <- c("red", "darkgrey", "blue")
    Col0 <- colorRampPalette(c00)(7)
    Col <- Col0[cut(Cex, br)]
    names(Col) <- names(Cex)

    op <- par(las=1)
    plot(0, type="n", xlim=c(-0.8,1.2), ylim=c(-1,1), xlab="Axis 1", ylab="Axis 2",
        main=ii)
    s2 <- scores(o)$sites
    s2 <- s2 / max(abs(s2))
    for (i in 1:nrow(s2))
        arrows(x0=0,y0=0,x1=s2[i,1],y1=s2[i,2], angle=20, length = 0.1, col="darkgrey")
    text(s2*1.05, labels=rownames(NN),cex=1,col=1)
    abline(h=0,v=0,lty=2)
    s1 <- scores(o)$species
    s1 <- s1 / max(abs(s1))
    text(s1[names(Col),]*0.8, labels=colnames(NN),cex=0.75, col=Col)
    for (ii in 1:200)
        lines(c(0.85, 0.9), rep(seq(-0.55, -0.95, len=200)[ii], 2), col=colorRampPalette(c00)(200)[ii])
    text(c(1,1,1)+0.1, c(-0.6, -0.75, -0.9),
        c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])))
    par(op)
}
dev.off()


## roadside count/habitat bias figure


rd_sign <- pop$Y0_Y1

Cex <- pop$DeltaObs
names(Cex) <- rownames(pop)
br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
#c00 <- c('#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba')
c00 <- c('#d7191c','darkgrey','#2b83ba')
#c00 <- c("red", "darkgrey", "blue")
Col0 <- colorRampPalette(c00)(7)
Col <- Col0[cut(Cex, br)]
names(Col) <- names(Cex)

pdf("~/GoogleWork/bam/PIF-AB/draft6/Fig5-count-habitat.pdf", width=7, height=7)
op <- par(las=1)
cx <- 1 # pop$EDR/100
topl <- pop[,c("DeltaR", "DeltaH")]
plot(topl,
    xlab="R", ylab="H",
    pch=c(19, 21, 19)[rd_sign+2],
    cex=cx,
    col=Col)
abline(h=0,v=0,lty=2)
#points(topl, cex=cx, pch=21)
text(topl[,1], topl[,2]-0.06,
    labels=ifelse(rd_sign != 0, rownames(topl), ""), cex=0.6)
#text(topl[,1]-0.15, topl[,2],
#    labels=ifelse(rd_sign != 0 & (rownames(topl) %in% c("TRES","AMRO")), rownames(topl), ""), cex=0.6)
for (ii in 1:200)
    lines(c(1.6, 1.7), rep(seq(-0.55, -0.95, len=200)[ii], 2)-0.4, col=colorRampPalette(c00)(200)[ii])
text(rep(2,3)+0.1, c(-0.6, -0.75, -0.9)-0.4,
    c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])),cex=0.9)
dev.off()


## story
library(intrval)


summary(pop$Npix/pop$Npif)
pop[which.max(pop$Npix/pop$Npif),]

exp(mean(log(pop$Npix/pop$Npif)))

table(pop$Npix_Npif)

all(pop$EDR < pop$MDD)

sum(pop$Npix > pop$Npif & pop[,c("Npix95lower","Npix95upper")] %)o(% pop[,c("Npif95lower","Npif95upper")])
sum(pop$Npix < pop$Npif & pop[,c("Npix95lower","Npix95upper")] %)o(% pop[,c("Npif95lower","Npif95upper")])
cat(as.character(pop$CommonName[pop$Npix < pop$Npif & pop[,c("Npix95lower","Npix95upper")] %)o(% pop[,c("Npif95lower","Npif95upper")]]), sep=", ")
sum(pop[,c("Npix95lower","Npix95upper")] %[o]% pop[,c("Npif95lower","Npif95upper")])

get_p <- function(cf,se) {
    zval <- cf/se
    2 * pnorm(-abs(zval))
}
round(get_p(c(1.3, 1.6, 0.37, 1.1, 0.087, -0.028), c(1.9,1.2,0.38,0.33,0.72,0.58)),2)


round(100*table(pop$Y0_Y1)/nrow(pop),1)
table(sign(pop$Y0 - pop$Y1))
dots_box_plot(log(pop[,c("Y0","Y1")]), lines=TRUE)

summary(pop$Y0/pop$Y1)

cor.test(pop$DeltaH,pop$DeltaR)

o
eigenvals(o)/sum(eigenvals(o))


round(pop[c("CAWA","WEWP","OSFL","BTNW","RUBL","BRCR"),c("Npix","Npix95lower","Npix95upper","Npif","Npif95lower","Npif95upper")],4)




