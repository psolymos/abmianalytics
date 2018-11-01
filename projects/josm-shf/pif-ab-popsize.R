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

## making pretty maps


library(sp)
library(rgeos)
library(raster)
library(cure4insect)
opar <- set_options(path = "w:/reports")
load_common_data()

stopifnot(all(rownames(KT)==rownames(kgrid)))
KT <- cure4insect:::.c4if$KT
rt <- cure4insect:::.read_raster_template()
#INSIDE <- KT$reg_nr!="Grassland" & coordinates(cure4insect:::.c4if$XY)[,2] > 50
INSIDE <- kgrid$useBCR6
r0 <- cure4insect:::.make_raster(ifelse(INSIDE, 1, 0), KT, rt)
v <- values(r0)
values(r0)[!is.na(v) & v==0] <- NA
cf <- function(n) rev(viridis::viridis(2*n)[(1+n):(n*2)])

pdf(file.path(ROOT, "out", "birds", "josmshf", "BCR6-maps.pdf"),
    onefile=TRUE, width=6, height=9)
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
    PREDS <- pxNcr[rownames(KT[INSIDE,]),]
    PREDS <- PREDS[match(rownames(KT), names(PREDS))]
    PREDS[is.na(PREDS)] <- 0
    r1 <- .make_raster(PREDS, KT, rt)
    r1 <- mask(r1, r0)
    plot(rt, col="darkgrey", axes=FALSE, box=FALSE, legend=FALSE, main=spp)
    #plot(r1, col=cf(100), axes=FALSE, box=FALSE, main=spp)
    plot(r1, col=cf(100), add=TRUE)
}
dev.off()


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
xn$BCR6
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

## not conclusive: but habCl effect pulls the contrast to more neutral for forest species
roadside_coef <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn)
    roadside_coef[[spp]] <- est7[,c("ROAD01","habCl:ROAD01")]
}
rcf <- t(exp(sapply(roadside_coef, function(z) c(RdOp=mean(z[,1]), RdCl=mean(rowSums(z))))))
boxplot(rcf,ylim=c(0,10))
abline(h=1,col=2)
hist(rcf[,2]/rcf[,1])

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

source("~/repos/bamanalytics/R/makingsense_functions.R")

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

#load("e:/peter/AB_data_v2016/out/birds/data/mean-qpad-estimates.Rdata")
#qpad_vals <- qpad_vals[rownames(tax),]

## roadside_bias, rai_pred, rai_data
#load("e:/peter/josm/2017/roadside_avoidance.Rdata")
#rsb <- t(sapply(roadside_bias, function(z) z[,1]))

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

xypt <- xn
coordinates(xypt) <- ~ POINT_X + POINT_Y
proj4string(xypt) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xypt <- spTransform(xypt, proj4string(r))
xypt2BCR <- over(xypt, BCR)

tabpt <- table(xn$HAB)


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

## forest classes
NAM_for <- c("BSpr", "BSprO", "Conif", "ConifO", "Decid", "DecidO",
    "Larch", "LarchO", "Mixwood", "MixwoodO", "Pine", "PineO")
pref_fun <- function(h) {
    NN <- h[NAM, "50%"] * 10^6 # back to individuals
    NN_for <- h[NAM[NAM %in% NAM_for], "50%"] * 10^6 # back to individuals
    sum(NN_for) / sum(NN)
}
pref <- sapply(NestAll, pref_fun)
quantile(pref,seq(0,1,0.25))

## roadside count effect within BCR 6
source("~/repos/abmianalytics/R/results_functions.R")
pr_fun_for_road <- function(est, X0, X1) {
    lam1 <- exp(apply(est, 1, function(z) X1 %*% z))
    Lam1 <- unname(colMeans(lam1))
    lam0 <- exp(apply(est, 1, function(z) X0 %*% z))
    Lam0 <- unname(colMeans(lam0))
    out <- rbind(Lam1, Lam0, Ratio=Lam1 / Lam0)
    cbind(mean=rowMeans(out), t(apply(out, 1, quantile, c(0.5, 0.025, 0.975))))
}
roadside_bias <- list()
tmp <- xn[!is.na(xypt2BCR) & xypt2BCR==24,]
tmp <- tmp[tmp$PCODE == "BBSAB",]
Xn_offroad <- model.matrix(getTerms(mods, "formula"), tmp)
colnames(Xn_offroad) <- fixNames(colnames(Xn_offroad))
Xn_onroad <- Xn_offroad
Xn_offroad[,c("ROAD01","habCl:ROAD01")] <- 0
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn_offroad)
    roadside_bias[[spp]] <- pr_fun_for_road(est7, X0=Xn_offroad, X1=Xn_onroad)
}
rsb <- t(sapply(roadside_bias, function(z) z[,1]))

## offsets

load("e:/peter/AB_data_v2016/out/birds/data/data-offset-covars.Rdata")
tmp <- xn[!is.na(xypt2BCR) & xypt2BCR==24,]
offdat <- offdat[rownames(offdat) %in% rownames(tmp),]
Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
    LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))
library(QPAD)
load_BAM_QPAD(3)
getBAMversion()
meanphi <- meantau <- meanphi0 <- meantau0 <- structure(rep(NA, length(SPP)), names=SPP)
for (spp in SPP) {
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    mi <- bestmodelBAMspecies(spp, type="BIC", model.sra=0:8)
    cat(spp, unlist(mi), "\n");flush.console()
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0
    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    meanphi[spp] <- mean(phi1)
    meantau[spp] <- mean(tau1)
    meanphi0[spp] <- cf0[1]
    meantau0[spp] <- cf0[2]
}
qpad_vals <- data.frame(Species=SPP, phi0=meanphi0, tau0=meantau0,
    phi=meanphi, tau=meantau)

## taxonomy etc ---
pop <- tax[,c("Species_ID", "English_Name", "Scientific_Name", "Spp")]

## LH stuff
library(lhreg)
data(lhreg_data)
compare_sets(lhreg_data$spp, rownames(pop))
setdiff(rownames(pop), lhreg_data$spp)
pop <- data.frame(pop, lhreg_data[match(rownames(pop), lhreg_data$spp),
    c("Mig", "Mig2", "Hab2", "Hab3", "Hab4")])

## pop size estimates ---
#pop$Npix <- NestTot[,"Mean"] # M males
#pop$Npix <- NestTot[,"50%"]
#pop$NpixLo <- NestTot[,"2.5%"]
#pop$NpixHi <- NestTot[,"97.5%"]
#pop$Npif <- (pif$Population_Estimate_unrounded / pif$Pair_Adjust) / 10^6 # M males
pop$Npix <- NestTot[,"50%"] * pif$Pair_Adjust # M inds
pop$NpixLo <- NestTot[,"2.5%"] * pif$Pair_Adjust
pop$NpixHi <- NestTot[,"97.5%"] * pif$Pair_Adjust
pop$Npif <- pif$Population_Estimate_unrounded / 10^6 # M inds
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
#pop <- droplevels(pop[!is.na(pop$Npif),])
pop <- pop[sort(rownames(pop)),]
Dall <- data.frame(Ahab=100*Ahab, Whab=100*Whab, sapply(NestAll, d_fun)[,rownames(pop)])

write.csv(pop, row.names=FALSE, file="~/GoogleWork/bam/PIF-AB/draft2/Table1-estimates.csv")
write.csv(Dall, file="~/GoogleWork/bam/PIF-AB/draft2/Table3-densities.csv")
write.csv(cbind(Dall[,1:2], n=tabpt), file="~/GoogleWork/bam/PIF-AB/draft2/Table2-habitats.csv")

## --- making the figures ---

## maps
pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig1-maps.pdf", width=12, height=9)
op <- par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(BCR2AB, col=c(NA, "grey", rep(NA, 11)), border=NA, main="Roadside surveys")
plot(AB, col=NA, border=1,add=TRUE)
plot(xy[xy@data$ROAD01 == 1,], add=TRUE, pch=19, col=1, cex=0.25)

plot(BCR2AB, col=c(NA, "grey", rep(NA, 11)), border=NA, main="Off-road surveys")
plot(AB, col=NA, border=1,add=TRUE)
plot(xy[xy@data$ROAD01 == 0,], add=TRUE, pch=19, col=1, cex=0.25)
par(op)
dev.off()

xylonlat <- spTransform(xy, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
plot(density(coordinates(xylonlat)[,2]), type="n")
lines(density(coordinates(xylonlat)[xylonlat@data$PCODE!="BBSAB",2]), col=2)
lines(density(coordinates(xylonlat)[xylonlat@data$PCODE=="BBSAB",2]), col=4)

library(MASS)
library(KernSmooth)
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

pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig2-popsize-old.pdf", width=8, height=8)
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


pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig2-popsize.pdf", width=7, height=7)
op <- par(mfrow=c(1,1), las=1, mar=c(4,4,1,2))
pch <- 19 # ifelse(pop$Npif %[]% list(pop$NpixLo, pop$NpixHi), 21, 19)
plot(pop[,c("Npif", "Npix")], xlab=expression(N[PIF]), ylab=expression(N[PIX]),
    type="n",
    pch=pch, col="#00000080", xlim=c(0, max(pop[,c("Npif", "Npix")])),
    ylim=c(0, max(pop[,c("Npif", "Npix")])))
abline(0,1, lty=2)
Min <- 3*2
Siz <- 18*2
di <- sqrt(pop[,"Npif"]^2+pop[,"Npix"]^2) > Min
pp <- pop[pop[,"Npif"] < Min & pop[,"Npix"] < Min, c("Npif", "Npix")]
ppp <- pop[pop[,"Npif"] < Min & pop[,"Npix"] < Min, ]
pch2 <- 19 # ifelse(ppp$Npif %[]% list(ppp$NpixLo, ppp$NpixHi), 21, 19)
di2 <- sqrt(pp[,"Npif"]^2+pp[,"Npix"]^2) > 1
pp <- pp*Siz/Min
pp[,1] <- pp[,1] + (60-Siz)
lines(c(0,60-Siz), c(0, 0), col="grey")
lines(c(0,60-Siz), c(Min, Siz), col="grey")
rect(0, 0, Min, Min)
rect(60-Siz, 0, 60-Siz+Min*Siz/Min, Min*Siz/Min, col="white")
lines(c(60-Siz, 60), c(0, Siz), col=1, lty=2)
points(pp, pch=pch2, col="#00000080")
points(pop[,c("Npif", "Npix")], pch=pch, col="#00000080")
text(pop[,"Npif"]+1.2*2, pop[,"Npix"]+0, labels=ifelse(di, rownames(pop), ""), cex=0.5)
text(pp[,"Npif"]+1.2*2, pp[,"Npix"]+0, labels=ifelse(di2, rownames(pp), ""), cex=0.5)
segments(x0=c(60-Siz+c(0, 1, 2, 3)*Siz/3), y0=rep(Siz, 4), y1=rep(Siz, 4)+0.5)
text(c(60-Siz+c(0, 1, 2, 3)*Siz/3), rep(Siz, 4)+1, 0:3)
dev.off()

pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig3-poprank.pdf", width=7, height=7)
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


pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig4-components.pdf", width=10, height=7)
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
text(i+off, tail(z[,i], 2), cex=0.8, rownames(z)[(nrow(z)-1):nrow(z)])

i <- 2
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

i <- 5
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 2), cex=0.8, rownames(z)[1:2])
text(i+off, tail(z[,i], 3), cex=0.8, rownames(z)[(nrow(z)-2):nrow(z)])

i <- 4
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

i <- 6
z <- mat[order(mat[,i], decreasing=TRUE),]
text(i+off, head(z[,i], 1), cex=0.8, rownames(z)[1])
#text(i+off, tail(z[,i], 1), cex=0.8, rownames(z)[nrow(z)])

for (i in 1:6) {
    zz1 <- format(mean(mat[,i]), trim = TRUE, scientific = FALSE, digits = 2)
    zz2 <- format(sd(mat[,i]), trim = TRUE, scientific = FALSE, digits = 2)
    mtext(paste("Mean =", zz1), side=1,at=i,line=3, cex=0.8)
    mtext(paste("SD =", zz2), side=1,at=i,line=4, cex=0.8)
}
par(op)
dev.off()

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
pdf("~/GoogleWork/bam/PIF-AB/draft3/FigX-varpart.pdf", width=6, height=6)
vpfun(prt, cutoff = 0, digits = 2, Xnames=c("R", "T", "A", "H"))
dev.off()


## road avoidance and ordination
library(vegan)
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

pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig6-ordination3.pdf", width=9, height=9)
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
#text(s1[names(Col),]*0.8, labels=colnames(NN),cex=0.75, col=4)
for (ii in 1:200)
    lines(c(0.85, 0.9), rep(seq(-0.55, -0.95, len=200)[ii], 2), col=colorRampPalette(c00)(200)[ii])
text(c(1,1,1)+0.1, c(-0.6, -0.75, -0.9),
    c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])))
par(op)
dev.off()

pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig6-ordination-all.pdf", width=9, height=9, onefile=TRUE)
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

library(intrval)
rd_over <- sapply(roadside_bias[rownames(pop)], function(z) z[1,3:4] %[o]% z[2,3:4])
rd_roadhigher <- sapply(roadside_bias[rownames(pop)], function(z) z[1,3:4] %[<o]% z[2,3:4])
rd_sign <- ifelse(rd_roadhigher, 1, -1)
rd_sign[rd_over] <- 0
table(rd_over, rd_roadhigher)
table(rd_sign)
plot(pop$DeltaR, pop$DeltaH, col=rd_sign+2)
abline(h=0, v=0)

Cex <- pop$DeltaObs
names(Cex) <- rownames(pop)
br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
#c00 <- c('#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba')
c00 <- c('#d7191c','darkgrey','#2b83ba')
#c00 <- c("red", "darkgrey", "blue")
Col0 <- colorRampPalette(c00)(7)
Col <- Col0[cut(Cex, br)]
names(Col) <- names(Cex)

pdf("~/GoogleWork/bam/PIF-AB/draft3/Fig5-count-habitat.pdf", width=7, height=7)
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
    labels=ifelse(pop$DeltaR %][% c(-1.5, 0.9), rownames(topl), ""), cex=0.6)
    for (ii in 1:200)
        lines(c(1.7, 1.9), rep(seq(-0.55, -0.95, len=200)[ii], 2), col=colorRampPalette(c00)(200)[ii])
    text(c(2.5,2.5,2.5)+0.1, c(-0.6, -0.75, -0.9),
        c(expression(N[PIX] < N[PIF]),expression(N[PIX] == N[PIF]),expression(N[PIX] > N[PIF])))
dev.off()



## ranking

pdf("~/GoogleWork/bam/PIF-AB/draft2/Fig2-popsize.pdf", width=7, height=7)
op <- par(mfrow=c(1,1), las=1, mar=c(4,4,1,2))
plot(rank(pop$Npif), rank(pop$Npix), xlab="PIF rank", ylab="PIX rank",
    type="n",
    pch=19, col="#00000080",
    ylim=c(0, nrow(pop)), xlim=c(0, nrow(pop)))
abline(0,1, lty=2)
text(rank(pop$Npif), rank(pop$Npix), labels=rownames(pop), cex=0.4)
par(op)
dev.off()


## results numbers
summary(exp(pop$DeltaObs))
table(pop$Npif < pop$NpixLo)
table(pop$Npif > pop$NpixHi)
table(pop$Npif %[]% list(pop$NpixLo, pop$NpixHi))

# AMCR AMGO AMRO BARS BBMA BRBL CCSP EAPH EUST HOSP HOWR SAVS SOSP VESP
rownames(pop)[pop$Npif > pop$NpixHi]
# BANS BAOR CLSW EAKI GRCA LCSP MODO NESP RWBL TRESVEER YHBL

rownames(pop)[pop$Npif %[]% list(pop$NpixLo, pop$NpixHi)]





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

