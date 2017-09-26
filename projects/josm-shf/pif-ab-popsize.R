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

## --- calculating total pop size based on predictions ---

STAGE <- list(veg = 5) # hab=5, hab+clim=6, hab+clim+shf=7

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

PREDS <- matrix(0, sum(kgrid$useBCR6), length(SPP))
rownames(PREDS) <- rownames(kgrid)[kgrid$useBCR6]
colnames(PREDS) <- SPP
PREDS0 <- PREDS

AREA_ha <- (1-kgrid$pWater) * kgrid$Area_km2 * 100
AREA_ha <- AREA_ha[kgrid$useBCR6]

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


## road effects

roadside_bias <- list()
Xn_onroad <- Xn_offroad <- Xn[Xn[,"ROAD01"] == 1,]
Xn_offroad[,c("ROAD01","habCl:ROAD01")] <- 0
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    off1sp <- OFF[,spp]
    off1sp <- off1sp[Xn[,"ROAD01"] == 1]
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn)
    pr_on <- pr_fun_for_gof(est7, Xn_onroad, off=off1sp)
    pr_off <- pr_fun_for_gof(est7, Xn_offroad, off=off1sp)
    roadside_bias[[spp]] <- c(on=mean(pr_on), off=mean(pr_off), onoff=mean(pr_on)/mean(pr_off))
}
rsb <- data.frame(do.call(rbind, roadside_bias))
rsb$spp <- SPP
write.csv(rsb, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/results/roadside_bias.csv")
#save(roadside_bias, file=file.path(ROOT, "josmshf", "roadside_bias.Rdata"))

## roadside avoidance index
rai_data <- data.frame(HAB=xn$hab1ec, ROAD=xn$ROAD01)
rai_pred <- matrix(0, nrow(Xn), nrow(tax))
rownames(rai_pred) <- rownames(rai_data) <- rownames(xn)
colnames(rai_pred) <- SPP
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    off1sp <- OFF[,spp]
    resn <- loadSPP(file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xn)
    rai_pred[,spp] <- pr_fun_for_gof(est7, Xn, off=off1sp)
}
save(roadside_bias, rai_pred, rai_data,
    file="e:/peter/josm/2017/roadside_avoidance.Rdata")


## ---


## PIF table
pif <- read.csv("~/Dropbox/bam/PIF-AB/popBCR-6AB_v2_22-May-2013.csv")
mefa4::compare_sets(tax$English_Name, pif$Common_Name)
setdiff(tax$English_Name, pif$Common_Name)
pif <- pif[match(tax$English_Name, pif$Common_Name),]


## roadside_bias
load(file.path(ROOT, "out", "birds", "josmshf", "roadside_bias.Rdata"))

load(file.path(ROOT, "out", "birds", "data", "mean-qpad-estimates.Rdata"))
qpad_vals <- qpad_vals[rownames(tax),]

## roadside avoidance
library(mefa4)
load(file.path(ROOT, "out", "birds", "josmshf", "roadside_avoidance.Rdata"))
tmp <- cbind(ROAD=rai_data$ROAD, rai_pred)
rai <- groupSums(tmp[BBn[,1],], 1, rai_data$HAB[BBn[,1]], TRUE)
rai <- t(t(rai) / colSums(rai))
RAI <- 1 - colSums(rai[,1] * rai)
summary(RAI)
RAIc <- RAI-RAI["ROAD"]

#yy <- cbind(ALL=1, ROAD=xnn[BBn[,1],"ROAD01"],
#    ifelse(as.matrix(yyn[BBn[,1],]) > 0, 1, 0))
#rai <- groupSums(yy, 1, xnn$hab1[BBn[,1]], TRUE)
#n <- rai[,"ALL"]
#rai <- rai[,-1]
#rai <- t(t(rai) / colSums(rai))
#sai <- groupSums(yy, 1, xnn$hab1[BBn[,1]], TRUE)
#RAI <- 1 - colSums(rai[,1] * rai)

pop <- tax[,c("Species_ID", "English_Name", "Scientific_Name", "Spp")]
pop$RAI <- RAI[match(rownames(pop), names(RAI))]
pop$RAIc <- RAIc[match(rownames(pop), names(RAIc))]
pop$RAIroad <- RAI["ROAD"]
pop$Don <- roadside_bias[rownames(pop), "on"]
pop$Doff <- roadside_bias[rownames(pop), "off"]
pop$DeltaRoad <- roadside_bias[rownames(pop), "onoff"]
pop$Nqpad <- colSums(PREDS*AREA_ha) / 10^6 # M males
pop$Nqpad[pop$Nqpad > 1000] <- NA
pop$Npif <- (pif$Pop_Est / pif$Pair_Adjust) / 10^6 # M males
pop$DeltaObs <- pop$Nqpad / pop$Npif
pop$TimeAdj <- pif$Time_Adjust
pop$MDD <- pif$Detection_Distance_m
pop$p3 <- 1-exp(-3 * qpad_vals$phi0)
pop$EDR <- qpad_vals$phi0 * 100
pop$DeltaTime <- (1/pop$p3)/pop$TimeAdj
pop$DeltaDist <- pop$MDD^2 / pop$EDR^2
pop$DeltaExp <- pop$DeltaRoad * pop$DeltaTime * pop$DeltaDist
pop$DeltaRes <- pop$DeltaObs / pop$DeltaExp
pop <- pop[rowSums(is.na(pop))==0,]

#write.csv(pop, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/qpad-pif-results.csv")

boxplot(log(pop[,c("DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")]))
abline(h=0, col=2)

boxplot(log(pop[,c("DeltaObs", "DeltaExp")]))
abline(h=0, col=2)

mat <- log(pop[,c("DeltaObs", "DeltaExp", "DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")])
rnd <- runif(nrow(pop), -0.1, 0.1)
boxplot(mat, range=0)
for (i in 2:ncol(mat))
    segments(x0=i+rnd-1, x1=i+rnd, y0=mat[,i-1], y1=mat[,i], col="lightgrey")
for (i in 1:ncol(mat))
    points(i+rnd, mat[,i], col="darkgrey", pch=19)
abline(h=0, col=2, lwd=2)
boxplot(mat, range=0, add=TRUE)

with(pop, plot(RAI, log(DeltaRes), type="n"))
abline(h=0, v=RAI["ROAD"], col=2, lwd=2)
with(pop, text(RAI, log(DeltaRes), rownames(pop), cex=0.75))

boxplot(pop[,c("Npif", "Nqpad")], ylim=c(0,10))

## --- making the figures ---

