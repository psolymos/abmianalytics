library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"

e <- new.env()
#load(file.path(ROOT, "data", "data-full-withrevisit.Rdata"), envir=e)
load(file.path(ROOT, "out", "birds", "data", "data-wrsi.Rdata"), envir=e)
TAX <- droplevels(e$TAX)
TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")
rm(e)

e <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-josmshf.Rdata"), envir=e)
xn <- e$DAT
mods <- e$mods
OFF <- e$OFF
yy <- e$YY
BB <- e$BB
rm(e)

fln <- list.files(file.path(ROOT, "out", "birds", "results", "josmshf"))
fln <- sub("birds_abmi-josmshf_", "", fln)
fln <- sub(".Rdata", "", fln)

SPP <- fln
tax <- droplevels(TAX[SPP,])

tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv0$Sector2 <- factor(ifelse(is.na(tv0$Sector), "NATIVE", as.character(tv0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
tv0$HAB <- paste0(tv0$Type, ifelse(tv0$AGE %in% c("5", "6", "7", "8", "9"), "O", ""))
tv0$HAB[tv0$Combined %in% c("RoadTrailVegetated", "RoadVegetatedVerge",
    "RailVegetatedVerge")] <- "Verges"
#tv0[tv0$HAB=="SoftLin",c("Combined", "HAB")]

if (FALSE) {
source("~/Dropbox/courses/st-johns-2017/R/diagnostics-functions.R")
xn$counts <- as.numeric(yy[,"ALFL"])
xn$counts01 <- ifelse(xn$counts>0, 1, 0)
plot(counts ~ wtAge + pAspen + ClosedCanopy + pWater + X + Y +
    xPET + xAHM + Succ_KM + Alien_KM, xn, B=0)
siplot(counts01 ~ wtAge + pAspen + ClosedCanopy + pWater + X + Y +
    xPET + xAHM + Succ_KM + Alien_KM, xn, B=0)
}

## --- calculating total pop size based on predictions ---

STAGE <- list(veg = 7) # hab=5, hab+clim=6, hab+clim+shf=7

OUTDIR1 <- paste0("e:/peter/josm/2017/stage", STAGE$veg, "/pred1")
OUTDIRB <- paste0("e:/peter/josm/2017/stage", STAGE$veg, "/predB")


load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
#source("~/repos/bragging/R/glm_skeleton.R")
#source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/abmianalytics/R/maps_functions.R")
regs <- levels(kgrid$LUFxNSR)
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$useBCR6 <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"

AREA_ha <- (1-kgrid$pWater) * kgrid$Area_km2 * 100
#AREA_ha <- kgrid$Area_km2 * 100
AREA_ha <- AREA_ha[kgrid$useBCR6]
names(AREA_ha) <- rownames(kgrid)[kgrid$useBCR6]

PREDS <- matrix(0, sum(kgrid$useBCR6), length(SPP))
rownames(PREDS) <- rownames(kgrid)[kgrid$useBCR6]
colnames(PREDS) <- SPP
PREDS0 <- PREDS

for (spp in SPP) {
    cat(spp, "--------------------------------------\n");flush.console()
    fl <- list.files(file.path(OUTDIR1, spp))
    ssRegs <- gsub("\\.Rdata", "", fl)
    pxNcr <- NULL
    #pxNrf <- NULL
    for (i in ssRegs) {
        cat(spp, i, "\n");flush.console()
        load(file.path(OUTDIR1, spp, paste0(i, ".Rdata")))
        rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
        pxNcr <- rbind(pxNcr, pxNcr1)
        #pxNrf <- rbind(pxNrf, pxNrf1)
    }
    PREDS[,spp] <- pxNcr[rownames(PREDS),]
    #PREDS0[,spp] <- pxNrf[rownames(PREDS0),]
}
N <- colSums(PREDS*AREA_ha) / 10^6
save(AREA_ha, N, PREDS, file=file.path(OUTDIR1, "predictions.Rdata"))

## getting habitat stuf

fl <- list.files(file.path(OUTDIR1, SPP[1]))
ssRegs <- gsub("\\.Rdata", "", fl)
Aveg <- NULL
for (i in ssRegs) {
    load(file.path(ROOT, "out", "transitions", paste0(i, ".Rdata")))
    tmp <- trVeg[rownames(trVeg) %in% names(AREA_ha),]
    tmp <- (1-kgrid[rownames(tmp), "pWater"]) * tmp
    Aveg <- rbind(Aveg, colSums(tmp))
}
rownames(Aveg) <- ssRegs
Aveg <- Aveg / 10^4

ch2veg <- t(sapply(strsplit(colnames(trVeg), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2veg <- data.frame(ch2veg)
colnames(ch2veg) <- c("rf","cr")
rownames(ch2veg) <- colnames(Aveg)
ch2veg$isHF <- ch2veg$cr %in% c("BorrowpitsDugoutsSumps",
    "Canals", "CCConif0", "CCConif1", "CCConif2", "CCConif3", "CCConif4",
    "CCConifR", "CCDecid0", "CCDecid1", "CCDecid2", "CCDecid3", "CCDecid4",
    "CCDecidR", "CCMixwood0", "CCMixwood1", "CCMixwood2", "CCMixwood3",
    "CCMixwood4", "CCMixwoodR", "CCPine0", "CCPine1", "CCPine2",
    "CCPine3", "CCPine4", "CCPineR",
    "CultivationCropPastureBareground", "HighDensityLivestockOperation",
    "IndustrialSiteRural",
    "MineSite",
    "MunicipalWaterSewage", "OtherDisturbedVegetation",
    "PeatMine", "Pipeline", "RailHardSurface",
    "RailVegetatedVerge", "Reservoirs", "RoadHardSurface", "RoadTrailVegetated",
    "RoadVegetatedVerge", "RuralResidentialIndustrial", "SeismicLine",
    "TransmissionLine", "Urban", "WellSite",
    "WindGenerationFacility")
ch2veg$HAB <- tv0$HAB[match(ch2veg$cr, tv0$Combined)]
ch2veg$HAB[is.na(ch2veg$HAB)] <- "Swamp"
ch2veg$HAB0 <- tv0$HAB[match(ch2veg$rf, tv0$Combined)]
ch2veg$HAB0[is.na(ch2veg$HAB0)] <- "Swamp"
AvegH <- groupSums(Aveg, 2, ch2veg$HAB)
sum(colSums(AvegH))/100
sum(AREA_ha)/100
ch2veg$Area_ha <- colSums(Aveg)

Nhab <- list()
for (spp in SPP) {
    cat(spp, "--------------------------------------\n");flush.console()
#    fl <- list.files(file.path(OUTDIR1, spp))
#    ssRegs <- gsub("\\.Rdata", "", fl)
    hbNcr <- 0
    for (i in ssRegs) {
        cat(spp, i, "\n");flush.console()
        aa <- AvegH[i,]
        ee <- new.env()
        load(file.path(OUTDIR1, spp, paste0(i, ".Rdata")), envir=ee)
        e <- new.env()
        load(file.path(OUTDIRB, spp, paste0(i, ".Rdata")), envir=e)

if (FALSE) {
sum(aa)
sum(AREA_ha[names(e$Cells)])

str(ee$pxNcr1)
str(e$pxNcrB)
sum(ee$pxNcr1[,1], na.rm=TRUE)
sum(e$pxNcrB[,1], na.rm=TRUE)

str(ee$hbNcr1)
str(e$hbNcrB)
sum(ee$hbNcr1[,1], na.rm=TRUE)
sum(e$hbNcrB[,1], na.rm=TRUE)

sum(ee$pxNcr1[,1]*AREA_ha[names(ee$Cells)], na.rm=TRUE)
sum(e$pxNcrB[,1]*AREA_ha[names(e$Cells)], na.rm=TRUE)
sum(ee$hbNcr1[,1]*Aveg[i,], na.rm=TRUE)
sum(e$hbNcrB[,1]*Aveg[i,], na.rm=TRUE)
}

        tmp <- e$hbNcrB
        tmp[is.na(tmp)] <- 0
        tmp <- groupSums(tmp*Aveg[i,], 1, ch2veg$HAB)
        hbNcr <- hbNcr + tmp
    }
    Nhab[[spp]] <- hbNcr
}
save(AvegH, Nhab, file=file.path(OUTDIRB, "predictions_HAB.Rdata"))

## getting CIs for Stage 7

ssRegs <- gsub("\\.Rdata", "", list.files(file.path(OUTDIR1, SPP[1])))
PREDSCI <- array(0, c(length(ssRegs), length(SPP), 100))
dimnames(PREDSCI) <- list(ssRegs, SPP, NULL)
#PREDSCI0 <- PREDSCI
## do not need to subtract water: it was not accounted for in predictions
#AA <- (1-kgrid$pWater) * kgrid$Area_km2 * 100
AA <- kgrid$Area_km2 * 100
names(AA) <- rownames(kgrid)

for (i in ssRegs) {
    cat(i, "--------------------------------------\n");flush.console()
    for (spp in SPP) {
        cat(i, spp, "\n");flush.console()
        e <- new.env()
        load(file.path(OUTDIRB, spp, paste0(i, ".Rdata")), envir=e)
        pxNcr <- e$pxNcrB
        #pxNrf <- e$pxNrfB
        rownames(pxNcr) <- names(e$Cells)
        #rownames(pxNrf) <- names(e$Cells)
        PREDSCI[i,spp,] <- colSums(pxNcr*AA[names(e$Cells)])
        #PREDSCI0[i,spp,] <- colSums(pxNrf*AA[names(e$Cells)])
    }
}
save(PREDSCI, file=file.path(OUTDIRB, "predictionsCI.Rdata"))

## looking at results

e <- new.env()
load("e:/peter/josm/2017/stage5/pred1/predictions.Rdata", envir=e)
N5 <- e$N
e <- new.env()
load("e:/peter/josm/2017/stage6/pred1/predictions.Rdata", envir=e)
N6 <- e$N
e <- new.env()
load("e:/peter/josm/2017/stage7/pred1/predictions.Rdata", envir=e)
N7 <- e$N
rm(e)

NN <- data.frame(N5=N5, N6=N6, N7=N7)
NN$spp <- rownames(NN)
write.csv(NN, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/results/PopSize567.csv")

summary(NN)
pairs(NN[NN$N5 < 20 & NN$N6 < 20 & NN$N7 < 20,])

## --- evaluating AUC based on external data ---

pr_fun_for_gof <- function(est, X, off=0) {
    if (is.null(dim(est))) {
        mu0 <- drop(X %*% est)
        exp(mu0 + off)
    } else {
        mu0 <- apply(est, 1, function(z) X %*% z)
        rowMeans(exp(mu0 + off))
    }
}
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}

## out-of-sample set
bunique <- unique(BB)
INTERNAL <- 1:nrow(xn) %in% bunique
ss <- which(INTERNAL)
ss1 <- which(!INTERNAL)
bid <- xn$bootid
levels(bid) <- sapply(strsplit(levels(bid), "\\."), "[[", 1)
up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
#    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()
Terms <- getTerms(mods, "list")
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))

#spp <- "OVEN"
res_AUC <- list()
for (spp in fln) {
    cat(spp, "\n");flush.console()
    y1sp <- yy[,spp]
    y10sp <- ifelse(y1sp > 0, 1, 0)
    off1sp <- OFF[,spp]
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est0 <- getEst(resn, stage=0, na.out=FALSE, Xn)
    est5 <- getEst(resn, stage=5, na.out=FALSE, Xn)
    est6 <- getEst(resn, stage=6, na.out=FALSE, Xn)
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn)
    j <- 1:nrow(est7)
    pr0o <- pr_fun_for_gof(est0[j,], Xn[ss1,], off=off1sp[ss1])
    pr5o <- pr_fun_for_gof(est5[j,], Xn[ss1,], off=off1sp[ss1])
    pr6o <- pr_fun_for_gof(est6[j,], Xn[ss1,], off=off1sp[ss1])
    pr7o <- pr_fun_for_gof(est7[j,], Xn[ss1,], off=off1sp[ss1])
    pr0i <- pr_fun_for_gof(est0[j,], Xn[ss,], off=off1sp[ss])
    pr5i <- pr_fun_for_gof(est5[j,], Xn[ss,], off=off1sp[ss])
    pr6i <- pr_fun_for_gof(est6[j,], Xn[ss,], off=off1sp[ss])
    pr7i <- pr_fun_for_gof(est7[j,], Xn[ss,], off=off1sp[ss])
    auc0o <- simple_auc(simple_roc(y10sp[ss1], pr0o))
    auc5o <- simple_auc(simple_roc(y10sp[ss1], pr5o))
    auc6o <- simple_auc(simple_roc(y10sp[ss1], pr6o))
    auc7o <- simple_auc(simple_roc(y10sp[ss1], pr7o))
    auc0i <- simple_auc(simple_roc(y10sp[ss], pr0i))
    auc5i <- simple_auc(simple_roc(y10sp[ss], pr5i))
    auc6i <- simple_auc(simple_roc(y10sp[ss], pr6i))
    auc7i <- simple_auc(simple_roc(y10sp[ss], pr7i))
    res_AUC[[spp]] <- c(
        auc0o=auc0o, auc5o=auc5o, auc6o=auc6o, auc7o=auc7o,
        auc0i=auc0i, auc5i=auc5i, auc6i=auc6i, auc7i=auc7i)
}

AUC <- data.frame(do.call(rbind, res_AUC))
AUC$k5 <- (AUC$auc5i-AUC$auc0i) / (1-AUC$auc0i)
AUC$k6 <- (AUC$auc6i-AUC$auc0i) / (1-AUC$auc0i)
AUC$k7 <- (AUC$auc7i-AUC$auc0i) / (1-AUC$auc0i)
AUC$spp <- rownames(AUC)

write.csv(AUC, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/results/AUC.csv")
write.csv(NN, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/results/PopSize567.csv")


## road effects

pr_fun_for_road <- function(est, X0, X1) {
    lam1 <- exp(apply(est, 1, function(z) X1 %*% z))
    Lam1 <- unname(colMeans(lam1))
    lam0 <- exp(apply(est, 1, function(z) X0 %*% z))
    Lam0 <- unname(colMeans(lam0))
    out <- rbind(Lam1, Lam0, Ratio=Lam1 / Lam0)
    cbind(mean=rowMeans(out), t(apply(out, 1, quantile, c(0.5, 0.025, 0.975))))
}
roadside_bias <- list()
Xn_onroad <- Xn_offroad <- Xn[Xn[,"ROAD01"] == 1,]
Xn_offroad[,c("ROAD01","habCl:ROAD01")] <- 0
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn)
    roadside_bias[[spp]] <- pr_fun_for_road(est7, X0=Xn_offroad, X1=Xn_onroad)
}
#rsb <- data.frame(do.call(rbind, roadside_bias))
#rsb$spp <- SPP
#write.csv(rsb, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/results/roadside_bias.csv")
#save(roadside_bias, file=file.path(ROOT, "josmshf", "roadside_bias.Rdata"))

## roadside avoidance index
rai_data <- data.frame(HAB=xn$hab1ec, ROAD=xn$ROAD01)
rai_pred <- matrix(0, nrow(Xn), nrow(tax))
rownames(rai_pred) <- rownames(rai_data) <- rownames(xn)
colnames(rai_pred) <- SPP
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn)
    rai_pred[,spp] <- pr_fun_for_gof(est7, Xn, off=0)
}
save(roadside_bias, rai_pred, rai_data,
    file="e:/peter/josm/2017/roadside_avoidance.Rdata")

## --- unifying the bits and pieces ---

## PIF table
pif <- read.csv("~/GoogleWork/bam/PIF-AB/popBCR-6AB_v2_22-May-2013.csv")
mefa4::compare_sets(tax$English_Name, pif$Common_Name)
setdiff(tax$English_Name, pif$Common_Name)
pif <- pif[match(tax$English_Name, pif$Common_Name),]
rownames(pif) <- rownames(tax)

AUC <- read.csv("~/GoogleWork/bam/PIF-AB/results/AUC.csv")
rownames(AUC) <- AUC$spp
AUC$spp <- NULL
AUC <- AUC[rownames(tax),]

if (FALSE) {
NN <- read.csv("~/GoogleWork/bam/PIF-AB/results/PopSize567.csv")
rownames(NN) <- NN$spp
NN$spp <- NULL
NN <- NN[rownames(tax),]

load("e:/peter/josm/2017/stage7/predB/predictionsCI.Rdata")
rm(PREDSCI0)
PREDSCI <- PREDSCI[,rownames(tax),] / 10^6
N7B <- apply(PREDSCI, 2, colSums)
N7CI <- t(apply(N7B, 2, quantile, c(0.5, 0.025, 0.975)))
N7mean <- colMeans(N7B)
}

load("e:/peter/AB_data_v2016/out/birds/data/mean-qpad-estimates.Rdata")
qpad_vals <- qpad_vals[rownames(tax),]

## roadside_bias, rai_pred, rai_data
load("e:/peter/josm/2017/roadside_avoidance.Rdata")
rsb <- t(sapply(roadside_bias, function(z) z[,1]))

## AvegH, Nhab
load("e:/peter/josm/2017/stage7/predB/predictions_HAB.Rdata")

## bootstrap averaged pop size estimate
b_fun <- function(h) {
    N7tb <- c(Mean=mean(colSums(h) / 10^6),
        quantile(colSums(h) / 10^6, c(0.5, 0.025, 0.975)))
    N7hb <- t(apply(h, 1, function(z)
        c(Mean=mean(z / 10^6), quantile(z / 10^6, c(0.5, 0.025, 0.975)))))
    N7b <- rbind(N7hb, TOTAL=N7tb)
    N7b
}
NestAll <- lapply(Nhab, b_fun)
NestTot <- t(sapply(NestAll, function(z) z["TOTAL",]))

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(gstat)
library(raster)
#library(viridis)
load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))

r <- raster(file.path("~/Dropbox/courses/st-johns-2017",
    "data", "ABrasters", "dem.asc"))
slope <- terrain(r, opt="slope")
aspect <- terrain(r, opt="aspect")
hill <- hillShade(slope, aspect, 40, 270)

od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/bcr/")
BCR <- readOGR(".", "BCR_Terrestrial_master") # rgdal
BCR <- spTransform(BCR, proj4string(r))
BCR <- gSimplify(BCR, tol=500, topologyPreserve=TRUE)
setwd(od)

od <- setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(r))
AB <- gUnaryUnion(AB, rep(1, nrow(AB))) # province
AB <- gSimplify(AB, tol=500, topologyPreserve=TRUE)
setwd(od)

BCR2AB <- gIntersection(AB, BCR, byid=TRUE)

For <- c("Decid", "Mixwood", "Pine", "Conif", "BSpr", "Larch")
xn$HAB <- paste0(xn$hab1, ifelse(xn$hab1 %in% For & xn$wtAge*200 >= 80, "O", ""))
compare_sets(colnames(AvegH), xn$HAB)
setdiff(colnames(AvegH), xn$HAB)
xnss <- nonDuplicated(xn, SS, TRUE)
xy <- xnss
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, proj4string(r))
xy2BCR <- over(xy, BCR)

tmp <- xnss[!is.na(xy2BCR) & xy2BCR==24,]
tmp <- tmp[!(tmp$PCODE != "BBSAB" & tmp$ROAD01 > 0),]
tab <- table(tmp$HAB, tmp$ROAD01)

#Ahab <- tab[,"0"] / sum(tab[,"0"])
Ahab <- colSums(AvegH[,rownames(tab)]) / sum(AvegH[,rownames(tab)])
Whab <- tab[,"1"] / sum(tab[,"1"])
AWhab <- data.frame(Ahab, Whab)
NAM <- names(Ahab)

h_fun <- function(h) {
    NN <- h[NAM, "50%"] * 10^6 # back to individuals
    DD <- NN / colSums(AvegH)[NAM] # density: males / ha
    sum(DD * Whab) / sum(DD * Ahab)
}
H <- sapply(NestAll, h_fun)

d_fun <- function(h) {
    NN <- h[NAM, "50%"] * 10^6 # back to individuals
    DD <- NN / colSums(AvegH)[NAM] # density: males / ha
    DD
}

## taxonomy etc ---
pop <- tax[,c("Species_ID", "English_Name", "Scientific_Name", "Spp")]
## pop size estimates ---
#pop$Npix <- NestTot[,"Mean"] # M males
pop$Npix <- NestTot[,"50%"]
pop$NpixLo <- NestTot[,"2.5%"]
pop$NpixHi <- NestTot[,"97.5%"]
pop$Npif <- (pif$Population_Estimate_unrounded / pif$Pair_Adjust) / 10^6 # M males
## roadside related metrics
pop$Y1 <- rsb[,"Lam1"]
pop$Y0 <- rsb[,"Lam0"]
## Tadj and EDR/MDD ---
pop$TimeAdj <- pif$Time_Adjust
#pop$p3 <- 1-exp(-3 * qpad_vals$phi0)
pop$p3 <- 1-exp(-3 * qpad_vals$phi)
pop$MDD <- pif$Detection_Distance_m
#pop$EDR <- qpad_vals$tau0 * 100
pop$EDR <- qpad_vals$tau * 100
## QAQC ---
pop$AUCin <- AUC$auc7i
pop$AUCout <- AUC$auc7o
pop$k7 <- AUC$k7
pop$DataQ <- pif$Data_Quality_Rating
pop$BbsVar <- pif$BBS_Variance_Rating
pop$SpSamp <- pif$Species_Sample_Rating
pop$H <- H
## Deltas ---
pop$DeltaObs <- log(pop$Npix / pop$Npif)
pop$DeltaR <- log(pop$Y0/pop$Y1)
pop$DeltaT <- log((1/pop$p3)/pop$TimeAdj)
pop$DeltaA <- log((1/pop$EDR^2) / (1/pop$MDD^2))
pop$DeltaH <- log(1/H) # we take inverse because H=1 is the PIF setup
pop$DeltaExp <- pop$DeltaR + pop$DeltaT + pop$DeltaA + pop$DeltaH
pop$epsilon <- pop$DeltaObs - pop$DeltaExp
## subset ---
pop <- droplevels(pop[rowSums(is.na(pop))==0,])
pop <- pop[sort(rownames(pop)),]
Dall <- data.frame(Ahab=100*Ahab, Whab=100*Whab, sapply(NestAll, d_fun)[,rownames(pop)])

write.csv(pop, row.names=FALSE, file="~/GoogleWork/bam/PIF-AB/draft1/Table1-estimates.csv")
write.csv(100*Dall[,1:2], file="~/GoogleWork/bam/PIF-AB/draft1/Table2-habitats.csv")
write.csv(Dall, file="~/GoogleWork/bam/PIF-AB/draft1/Table3-densities.csv")

## --- making the figures ---

## maps
pdf("~/GoogleWork/bam/PIF-AB/draft1/Fig1-maps.pdf", width=12, height=9)
op <- par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(BCR2AB, col=c(NA, "grey", rep(NA, 11)), border=NA, main="Roadside surveys")
plot(AB, col=NA, border=1,add=TRUE)
plot(xy[xy@data$PCODE=="BBSAB",], add=TRUE, pch=19, col=1, cex=0.25)

plot(BCR2AB, col=c(NA, "grey", rep(NA, 11)), border=NA, main="Off-road surveys")
plot(AB, col=NA, border=1,add=TRUE)
plot(xy[xy@data$PCODE!="BBSAB",], add=TRUE, pch=19, col=1, cex=0.25)
par(op)
dev.off()

xylonlat <- spTransform(xy, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
plot(density(coordinates(xylonlat)[,2]), type="n")
lines(density(coordinates(xylonlat)[xylonlat@data$PCODE!="BBSAB",2]), col=2)
lines(density(coordinates(xylonlat)[xylonlat@data$PCODE=="BBSAB",2]), col=4)

dots_box_plot <- function(mat, lines=FALSE, ...) {
    set.seed(1)
    rnd <- runif(nrow(mat), -0.1, 0.1)
    boxplot(mat, range=0, ...)
    if (lines)
        for (i in 2:ncol(mat))
            segments(x0=i+rnd-1, x1=i+rnd, y0=mat[,i-1], y1=mat[,i], col="lightgrey")
    for (i in 1:ncol(mat))
        points(i+rnd, mat[,i], pch=19, col="#00000080")
    #boxplot(mat, range=0, add=TRUE, col="#ff000020", names=NA)
    boxplot(mat, range=0, add=TRUE, col="#00000020", names=NA)
    invisible(NULL)
}

pdf("~/GoogleWork/bam/PIF-AB/draft1/Fig2-popsize-old.pdf", width=8, height=8)
op <- par(mfrow=c(2,2), las=1, mar=c(4,4,1,2))
dots_box_plot(pop[,c("Npif", "Npix")], lines=TRUE,
    ylab="Population Size (M singing inds.)", names=c("PIF", "PIX"))
dots_box_plot(log(pop[,c("Npif", "Npix")]), lines=TRUE,
    ylab="log Population Size (M singing inds.)", names=c("PIF", "PIX"))
plot(pop[,c("Npif", "Npix")], xlab=expression(N[PIF]), ylab=expression(N[PIX]),
    pch=19, col="#00000080", xlim=c(0, max(pop[,c("Npif", "Npix")])),
    ylim=c(0, max(pop[,c("Npif", "Npix")])))
abline(0,1)
plot(log(pop[,c("Npif", "Npix")]), xlab=expression(log(N[PIF])), ylab=expression(log(N[PIX])),
    pch=19, col="#00000080", xlim=range(log(pop[,c("Npif", "Npix")])),
    ylim=range(log(pop[,c("Npif", "Npix")])))
abline(0,1)
par(op)
dev.off()


pdf("~/GoogleWork/bam/PIF-AB/draft1/Fig2-popsize.pdf", width=7, height=7)
op <- par(mfrow=c(1,1), las=1, mar=c(4,4,1,2))
plot(pop[,c("Npif", "Npix")], xlab=expression(N[PIF]), ylab=expression(N[PIX]),
    type="n",
    pch=19, col="#00000080", xlim=c(0, max(pop[,c("Npif", "Npix")])),
    ylim=c(0, max(pop[,c("Npif", "Npix")])))
abline(0,1, lty=2)
Min <- 3
Siz <- 18
di <- sqrt(pop[,"Npif"]^2+pop[,"Npix"]^2) > Min
pp <- pop[pop[,"Npif"] < Min & pop[,"Npix"] < Min, c("Npif", "Npix")]
di2 <- sqrt(pp[,"Npif"]^2+pp[,"Npix"]^2) > 1
pp <- pp*Siz/Min
pp[,1] <- pp[,1] + (30-Siz)
lines(c(0,30-Siz), c(0, 0), col="grey")
lines(c(0,30-Siz), c(Min, Siz), col="grey")
rect(0, 0, Min, Min)
rect(30-Siz, 0, 30-Siz+Min*Siz/Min, Min*Siz/Min, col="white")
lines(c(30-Siz, 30), c(0, Siz), col=1, lty=2)
points(pp, pch=19, col="#00000080")
points(pop[,c("Npif", "Npix")], pch=19, col="#00000080")
text(pop[,"Npif"]+1.2, pop[,"Npix"]+0, labels=ifelse(di, rownames(pop), ""), cex=0.4)
text(pp[,"Npif"]+1.2, pp[,"Npix"]+0, labels=ifelse(di2, rownames(pp), ""), cex=0.4)
segments(x0=c(30-Siz+c(0, 1, 2, 3)*Siz/3), y0=rep(Siz, 4), y1=rep(Siz, 4)+0.5)
text(c(30-Siz+c(0, 1, 2, 3)*Siz/3), rep(Siz, 4)+1, 0:3)
par(op)
dev.off()


pdf("~/GoogleWork/bam/PIF-AB/draft1/Fig3-components.pdf", width=8, height=6)
op <- par(las=1)
mat <- pop[,c("DeltaObs", "DeltaExp", "DeltaR", "DeltaT", "DeltaA", "DeltaH")]
colnames(mat) <- c("OBS", "EXP", "R", "T", "A", "H")
par(las=1)
dots_box_plot(mat, ylab="Log Ratio")
abline(h=0, col=1, lwd=1,lty=2)
par(op)
dev.off()

mod <- lm(DeltaObs ~ DeltaR + DeltaT + DeltaA + DeltaH, pop)
summary(mod)
an <- anova(mod)
an$Percent <- 100 * an[["Sum Sq"]] / sum(an[["Sum Sq"]])
an <- an[c("Df", "Sum Sq", "Percent", "Mean Sq", "F value", "Pr(>F)")]
an
summary(mod)$sigma^2

mod2 <- step(lm(DeltaObs ~ (DeltaR + DeltaT + DeltaA + DeltaH)^2, pop), trace=0)
summary(mod2)
an2 <- anova(mod2)
an2$Percent <- 100 * an2[["Sum Sq"]] / sum(an2[["Sum Sq"]])
an2 <- an2[c("Df", "Sum Sq", "Percent", "Mean Sq", "F value", "Pr(>F)")]
an2

## road avoidance and ordination
library(vegan)
DD <- as.matrix(t(Dall[,-(1:2)]))
NN <- t(t(DD) * Dall$Ahab)
NN <- t(NN / rowSums(NN))
Cex <- pop$DeltaObs
names(Cex) <- rownames(pop)
br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
Col0 <- colorRampPalette(c("red", "darkgrey", "blue"))(7)
Col <- Col0[cut(Cex, br)]
names(Col) <- names(Cex)
o <- cca(NN)

pdf("~/GoogleWork/bam/PIF-AB/draft1/Fig4-ordination.pdf", width=9, height=9)
op <- par(las=1)
plot(0, type="n", xlim=c(-0.8,1.2), ylim=c(-1,1), xlab="Axis 1", ylab="Axis 2")
s2 <- scores(o)$sites
s2 <- s2 / max(abs(s2))
for (i in 1:nrow(s2))
    arrows(x0=0,y0=0,x1=s2[i,1],y1=s2[i,2], angle=20, length = 0.1, col="darkgrey")
text(s2*1.05, labels=rownames(NN),cex=1,col=1)
abline(h=0,v=0,lty=2)
s1 <- scores(o)$species
s1 <- s1 / max(abs(s1))
text(s1[names(Col),]*0.8, labels=colnames(NN),cex=0.75, col=Col)
legend("bottomright", bty="n", fill=Col0[c(1,4,7)], border=NA,
    legend=c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])))
par(op)
dev.off()



# ---


plot(Ahab, Whab, type="n", ylim=c(0,max(AWhab)), xlim=c(0,max(AWhab)))
abline(0,1)
text(Ahab, Whab, names(Ahab), cex=0.7)



summary(round(rowSums(mat[,-(1:2)]) - mat[,1], 12))
summary(round(mat[,2]+mat[,"epsilon"] - mat[,1], 12))

mod <- lm(DeltaObs ~ DeltaR + DeltaT + DeltaA + DeltaH, pop)
an <- anova(mod)
an$Percent <- 100 * an[["Sum Sq"]] / sum(an[["Sum Sq"]])
an <- an[c("Df", "Sum Sq", "Percent", "Mean Sq", "F value", "Pr(>F)")]
summary(mod)
an

cf <- coef(mod)
mat2 <- t(t(model.matrix(mod)) * cf)
mat2 <- cbind(Obs=pop$DeltaObs, mat2[,-1], eps=pop$DeltaObs-rowSums(mat2))
par(mfrow=c(2,1))
dots_box_plot(mat2)
abline(h=0, col=2, lwd=2)
dots_box_plot(mat)
abline(h=0, col=2, lwd=2)


par(las=1)
mat2 <- pop[,c("Npif", "Npix")]
colnames(mat2) <- c("PIF", "'Pixel'")
dots_box_plot(mat2, col="grey", "Population Size (M singing inds.)")

library(intrval)
table(OVER <- pop$Npif %[]% pop[,c("NpixLo", "NpixHi")])
rownames(pop)[OVER]

#v <- tanh(pop$RAI * 0.5)
#v <- plogis(pop$RAI)*2-1
#v <- pop$RAI
Cex <- pop$DeltaObs
names(Cex) <- rownames(pop)
br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
Col <- colorRampPalette(c("red", "black", "blue"))(7)[cut(Cex, br)]
names(Col) <- rownames(pop)

with(pop, plot(log(H), DeltaObs-DeltaExp, type="n", #xlim=c(-0.4, 0.4), ylim=c(-10,10),
    ylab=expression(Delta[OBS]-Delta[EXP]), xlab=expression(-Delta[H])))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(lm(I(DeltaObs-DeltaExp) ~ I(log(H)), pop), col=1)
with(pop, text(log(H), DeltaObs-DeltaExp, rownames(pop),
    cex=0.8, col=Col))
legend("topright", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))

with(pop, plot(DeltaObs-log(H), DeltaExp, type="n", #xlim=c(-0.4, 0.4), ylim=c(-10,10),
    ylab=expression(Delta[OBS]-Delta[H]), xlab=expression(Delta[R]+Delta[T]+Delta[A]+epsilon)))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(0, 1)
with(pop, text(DeltaObs-log(H), DeltaExp, rownames(pop), cex=0.8, col=Col))
legend("bottomright", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))

with(pop, plot(DeltaObs, log(H)))
with(pop, plot(epsilon, log(H)))

with(pop, plot(PropRd, log(DeltaRes), type="n",
    ylab=expression(log(Delta[Res])), xlab="Road Avoidance Index"))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(lm(log(DeltaRes) ~ PropRd, pop), col=1)
with(pop, text(PropRd, log(DeltaRes), rownames(pop),
    cex=0.8, col=Col))
legend("topleft", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))


o <- cca(rai)
plot(0,type="n", xlim=c(-1,1), ylim=c(-1,1))
s1 <- scores(o)$species
s1 <- s1 / max(abs(s1))
text(s1[names(Col),], labels=colnames(rai),cex=0.75, col=Col)
s2 <- scores(o)$sites
s2 <- s2 / max(abs(s2))
text(s2, labels=rownames(rai),cex=1,col=3)
points(s1[1,,drop=F],pch=3,col=4,cex=5)
abline(h=0,v=0,lty=2)
legend("topleft", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))

plot(Ahab, Whab, type="n", ylim=c(0,max(AWhab)), xlim=c(0,max(AWhab)))
abline(0,1)
text(Ahab, Whab, names(Ahab), cex=0.7)





op <- par(mar=c(1,1,1,1))
plot(hill, col=grey(0:100/100), legend=FALSE, bty="n",
    box=FALSE, axes=FALSE)
plot(r, legend=FALSE, col=topo.colors(50, alpha=0.35)[26:50], add=TRUE)
plot(BCR2AB, col=c("#00000060", NA, rep("#00000060", 11)), add=TRUE)
plot(xy[xy@data$PCODE!="BBSAB" & !(!is.na(xy2BCR) & xy2BCR==24),], add=TRUE, pch=19, col="white", cex=0.5)
plot(xy[xy@data$PCODE=="BBSAB" & !(!is.na(xy2BCR) & xy2BCR==24),], add=TRUE, pch=19, col="lightblue", cex=0.5)
plot(xy[xy@data$PCODE!="BBSAB" & !is.na(xy2BCR) & xy2BCR==24,], add=TRUE, pch=19, cex=0.5)
plot(xy[xy@data$PCODE=="BBSAB" & !is.na(xy2BCR) & xy2BCR==24,], add=TRUE, pch=19, col=4, cex=0.5)
legend("bottomleft", title="Surveys", pch=c(19, 19, 19, 21), bty="n",
    col=c("blue", "lightblue", "black", "black"),
    legend=c("BBS in BCR 6", "BBS outside", "Off road in BCR 6", "Off road outside"))
par(op)

tab0 <- table(xy@data$HAB)
ii1 <- xy@data$PCODE=="BBSAB" & !is.na(xy2BCR) & xy2BCR==24
tab1 <- table(xy@data$HAB, ii1)
summary(xnss$POINT_Y[ii1 & xnss$POINT_Y < 58]) # 56.51 latitude
ii2 <- xnss$POINT_Y < 56.51
tab2 <- table(xy@data$HAB, ii2)
df <- data.frame(tab0/sum(tab0), w1=tab1[,"TRUE"]/sum(tab1[,"TRUE"]),
    w2=tab2[,"TRUE"]/sum(tab2[,"TRUE"]))
all(NAM==rownames(df))

tmp <- t(as.matrix(df[,-1]))
rownames(tmp)[1] <- "a"
tmp <- tmp[,order(tmp[1,])]

par(las=1, mar=c(5,6,2,2))
barplot(100*tmp[1:2,], beside=TRUE, horiz=TRUE, xlab="% representation",
    legend.text=TRUE)


h_fun <- function(h, a, w) {
    NN <- h[NAM, "50%"] * 10^6 # back to individuals
    DD <- NN / colSums(AvegH)[NAM] # density: males / ha
    sum(DD * w) / sum(DD * a)
}
H1 <- sapply(NestAll, h_fun, a=df$Freq, w=df$w1)
H2 <- sapply(NestAll, h_fun, a=df$Freq, w=df$w2)
H3 <- sapply(NestAll, h_fun, a=df$Freq, w=df$w1*df$w2/sum(df$w1*df$w2))
names(H1) <- names(H2) <- names(H3) <- names(NestAll)

par(mfrow=c(2,2))
plot(H1,H2)
plot(log(H1[rownames(pop)]), pop$epsilon)
plot(log(H2[rownames(pop)]), pop$epsilon)
plot(log(H3[rownames(pop)]), pop$epsilon)


cor(cbind(pop$epsilon, H1[rownames(pop)], H2[rownames(pop)], H3[rownames(pop)]))

par(mfrow=c(1,3))
logH <- log(H1[rownames(pop)])
with(pop, plot(logH, epsilon, type="n", #xlim=c(-0.4, 0.4), ylim=c(-10,10),
    ylab=expression(Delta[RES]), xlab=expression(Delta[H])))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(lm(epsilon ~ logH, pop), col=1)
with(pop, text(logH, epsilon, rownames(pop),
    cex=0.8, col=Col))
legend("topleft", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))

logH <- log(H2[rownames(pop)])
with(pop, plot(logH, epsilon, type="n", #xlim=c(-0.4, 0.4), ylim=c(-10,10),
    ylab=expression(Delta[RES]), xlab=expression(Delta[H])))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(lm(epsilon ~ logH, pop), col=1)
with(pop, text(logH, epsilon, rownames(pop),
    cex=0.8, col=Col))
legend("topleft", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))

logH <- log(H3[rownames(pop)])
with(pop, plot(logH, epsilon, type="n", #xlim=c(-0.4, 0.4), ylim=c(-10,10),
    ylab=expression(Delta[RES]), xlab=expression(Delta[H])))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(lm(epsilon ~ logH, pop), col=1)
with(pop, text(logH, epsilon, rownames(pop),
    cex=0.8, col=Col))
legend("topleft", bty="n", fill=c(4,2), border=NA, legend=c("QPAD > PIF", "QPAD < PIF"))

par(mfrow=c(1,3))
plot(pop$Npix, pop$Npif, ylim=c(0,7));abline(0,1)
plot(pop$Npix, pop$Npif/H1[rownames(pop)], ylim=c(0,7));abline(0,1)
plot(pop$Npix, pop$Npif/H2[rownames(pop)], ylim=c(0,7));abline(0,1)

plot(log(pop$Npix/pop$Npif), log(H1[rownames(pop)]));abline(h=0,v=0)
plot(pop$epsilon, log(H1[rownames(pop)]));abline(h=0,v=0)
plot(pop$epsilon+log(H1[rownames(pop)]), log(H1[rownames(pop)]));abline(h=0,v=0)

pop$DeltaH <- log(H1[rownames(pop)])
pop$Res1 <- pop$DeltaObs - pop$DeltaT
pop$Res2 <- pop$DeltaObs - pop$DeltaT - pop$DeltaA
pop$Res3 <- pop$DeltaObs - pop$DeltaT - pop$DeltaA - pop$DeltaR
pop$Res4 <- pop$DeltaObs - pop$DeltaT - pop$DeltaA - pop$DeltaR - pop$DeltaH

par(mfrow=c(2,1))
mat <- pop[,c("DeltaObs", "DeltaR", "DeltaT", "DeltaA", "DeltaH", "epsilon")]
dots_box_plot(mat, col="grey", ylab=expression(Delta))
abline(h=0, col=2, lwd=2)

mat <- pop[,c("DeltaObs", "Res1", "Res2", "Res3", "Res4")]
dots_box_plot(mat, col="grey", ylab=expression(Delta))
abline(h=0, col=2, lwd=2)

yyy <- as.matrix(t(groupMeans(yy, 1, xn$ROAD01)))
r2 <- yyy[,"0"]/yyy[,"1"]
r2 <- r2[rownames(pop)]
plot(r2, exp(pop$DeltaR))

