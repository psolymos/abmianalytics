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
AREA_ha <- AREA_ha[kgrid$useBCR6]

PREDS <- matrix(0, sum(kgrid$useBCR6), length(SPP))
rownames(PREDS) <- rownames(kgrid)[kgrid$useBCR6]
colnames(PREDS) <- SPP
PREDS0 <- PREDS

for (spp in SPP) {
    cat(spp, "--------------------------------------\n");flush.console()
    fl <- list.files(file.path(OUTDIR1, spp))
    ssRegs <- gsub("\\.Rdata", "", fl)
    pxNcr <- pxNrf <- NULL
    for (i in ssRegs) {
        cat(spp, i, "\n");flush.console()
        load(file.path(OUTDIR1, spp, paste0(i, ".Rdata")))
        rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
        pxNcr <- rbind(pxNcr, pxNcr1)
        pxNrf <- rbind(pxNrf, pxNrf1)
    }
    PREDS[,spp] <- pxNcr[rownames(PREDS),]
    PREDS0[,spp] <- pxNrf[rownames(PREDS0),]
}
N <- colSums(PREDS*AREA_ha) / 10^6
save(AREA_ha, N, PREDS, PREDS0, file=file.path(OUTDIR1, "predictions.Rdata"))

## getting CIs for Stage 7

ssRegs <- gsub("\\.Rdata", "", list.files(file.path(OUTDIR1, SPP[1])))
PREDSCI <- array(0, c(length(ssRegs), length(SPP), 100))
dimnames(PREDSCI) <- list(ssRegs, SPP, NULL)
PREDSCI0 <- PREDSCI
AA <- (1-kgrid$pWater) * kgrid$Area_km2 * 100
names(AA) <- rownames(kgrid)

for (i in ssRegs) {
    cat(i, "--------------------------------------\n");flush.console()
    for (spp in SPP) {
        cat(i, spp, "\n");flush.console()
        e <- new.env()
        load(file.path(OUTDIRB, spp, paste0(i, ".Rdata")), envir=e)
        pxNcr <- e$pxNcrB
        pxNrf <- e$pxNrfB
        rownames(pxNcr) <- rownames(pxNrf) <- names(e$Cells)
        PREDSCI[i,spp,] <- colSums(pxNcr*AA[names(e$Cells)])
        PREDSCI0[i,spp,] <- colSums(pxNrf*AA[names(e$Cells)])
    }
}
save(PREDSCI, PREDSCI0, file=file.path(OUTDIRB, "predictionsCI.Rdata"))


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
pif <- read.csv("~/Dropbox/bam/PIF-AB/popBCR-6AB_v2_22-May-2013.csv")
mefa4::compare_sets(tax$English_Name, pif$Common_Name)
setdiff(tax$English_Name, pif$Common_Name)
pif <- pif[match(tax$English_Name, pif$Common_Name),]
rownames(pif) <- rownames(tax)

AUC <- read.csv("~/Dropbox/bam/PIF-AB/results/AUC.csv")
rownames(AUC) <- AUC$spp
AUC$spp <- NULL
AUC <- AUC[rownames(tax),]

NN <- read.csv("~/Dropbox/bam/PIF-AB/results/PopSize567.csv")
rownames(NN) <- NN$spp
NN$spp <- NULL
NN <- NN[rownames(tax),]

load("e:/peter/josm/2017/stage7/predB/predictionsCI.Rdata")
rm(PREDSCI0)
PREDSCI <- PREDSCI[,rownames(tax),] / 10^6
N7B <- apply(PREDSCI, 2, colSums)
N7CI <- t(apply(N7B, 2, quantile, c(0.5, 0.025, 0.975)))

load("e:/peter/AB_data_v2016/out/birds/data/mean-qpad-estimates.Rdata")
qpad_vals <- qpad_vals[rownames(tax),]

## roadside_bias, rai_pred, rai_data
load("e:/peter/josm/2017/roadside_avoidance.Rdata")

#rsb <- read.csv("~/Dropbox/bam/PIF-AB/results/roadside_bias.csv")
#rownames(rsb) <- rsb$spp
#rsb$spp <- NULL
#rsb <- rsb[rownames(tax),]
rsb <- t(sapply(roadside_bias, function(z) z[,1]))
rsbCI <- t(sapply(roadside_bias, function(z) z["Ratio",]))

tmp <- cbind(ROAD=rai_data$ROAD, rai_pred)
rai <- groupSums(tmp[BB[,1],], 1, rai_data$HAB[BB[,1]], TRUE)
rai <- t(t(rai) / colSums(rai))
RAI0 <- 1 - colSums(rai[,1] * rai)
summary(RAI0)
RAI <- log(RAI0)-log(RAI0["ROAD"])
RAIc <- RAI0-RAI0["ROAD"]

## taxonomy etc
pop <- tax[,c("Species_ID", "English_Name", "Scientific_Name", "Spp")]
## roadside related metrics
pop$RAI <- RAI[match(rownames(pop), names(RAI))]
pop$RAIc <- RAIc[match(rownames(pop), names(RAI))]
pop$Don <- rsb[,"Lam1"]
pop$Doff <- rsb[,"Lam0"]
pop$DeltaRoad <- rsb[,"Ratio"]
pop$DeltaRoadLo <- rsbCI[,"2.5%"]
pop$DeltaRoadHi <- rsbCI[,"97.5%"]
## pop size estimates
pop$Nqpad7mean <- NN$N7 # M males
pop$Nqpad7median <- N7CI[,"50%"]
pop$Nqpad7lo <- N7CI[,"2.5%"]
pop$Nqpad7hi <- N7CI[,"97.5%"]
pop$Npif <- (pif$Population_Estimate_unrounded / pif$Pair_Adjust) / 10^6 # M males
## Tadj and EDR/MDD
pop$TimeAdj <- pif$Time_Adjust
pop$MDD <- pif$Detection_Distance_m
pop$p3 <- 1-exp(-3 * qpad_vals$phi0)
pop$EDR <- qpad_vals$tau * 100
#pop$p3 <- 1-exp(-3 * qpad_vals$phi0)
#pop$EDR <- qpad_vals$tau * 100
## Deltas
pop$DeltaObs <- pop$Nqpad7mean / pop$Npif
pop$DeltaTime <- (1/pop$p3)/pop$TimeAdj
pop$DeltaDist <- pop$MDD^2 / pop$EDR^2
pop$DeltaExp <- pop$DeltaRoad * pop$DeltaTime * pop$DeltaDist
pop$DeltaRes <- pop$DeltaObs / pop$DeltaExp
## QAQC
pop$AUCin <- AUC$auc7i
pop$AUCout <- AUC$auc7o
pop$k7 <- AUC$k7
pop$DataQ <- pif$Data_Quality_Rating
pop$BbsVar <- pif$BBS_Variance_Rating
pop$SpSamp <- pif$Species_Sample_Rating

pop <- droplevels(pop[rowSums(is.na(pop))==0,])
pop <- pop[sort(rownames(pop)),]
summary(pop)
write.csv(pop, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/qpad-pif-results.csv")

## --- making the figures ---

boxplot(log(pop[,c("DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")]))
abline(h=0, col=2)

boxplot(log(pop[,c("DeltaObs", "DeltaExp")]))
abline(h=0, col=2)


dots_box_plot <- function(mat, ylab="", ...) {
    rnd <- runif(nrow(mat), -0.1, 0.1)
    boxplot(mat, range=0, ylab=ylab)
    for (i in 2:ncol(mat))
        segments(x0=i+rnd-1, x1=i+rnd, y0=mat[,i-1], y1=mat[,i], col="lightgrey")
    for (i in 1:ncol(mat))
        points(i+rnd, mat[,i], pch=19, ...)
    boxplot(mat, range=0, add=TRUE, col="#ff000020")
    invisible(NULL)
}

Col <- c("beige", "yellow", "orange")[pop$DataQ]
mat <- log(pop[,c("DeltaObs", "DeltaExp", "DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")])
colnames(mat) <- c("Obs", "Exp", "Road", "Time", "Dist", "Res")
dots_box_plot(mat, col="grey", ylab=expression(log(Delta)))
abline(h=0, col=2, lwd=2)

mat2 <- pop[,c("Npif", "Nqpad7mean")]
colnames(mat2) <- c("PIF", "QPAD")
dots_box_plot(mat2, col="grey", "Population Size (M singing inds.)")

library(intrval)
table(OVER <- pop$Npif %[]% pop[,c("Nqpad7lo", "Nqpad7hi")])
rownames(pop)[OVER]

#v <- tanh(pop$RAI * 0.5)
#v <- plogis(pop$RAI)*2-1
#v <- pop$RAI
Cex <- log(pop$DeltaObs)
names(Cex) <- rownames(pop)
br <- c(-Inf, -2, -1, -0.1, 0.1, 1, 2, Inf)
Col <- colorRampPalette(c("red", "black", "blue"))(7)[cut(Cex, br)]
names(Col) <- rownames(pop)

with(pop, plot(RAI, log(DeltaRes), type="n", xlim=c(-0.4, 0.4), ylim=c(-10,10),
    ylab=expression(log(Delta[Res])), xlab="Road Avoidance Index"))
abline(h=0, v=0, col=1, lwd=1, lty=2)
abline(lm(log(DeltaRes) ~ RAI, pop), col=1)
with(pop, text(RAI, log(DeltaRes), rownames(pop),
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
