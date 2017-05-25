library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"

level <- 0.9

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

## names etc
e <- new.env()
#load(file.path(ROOT, "data", "data-full-withrevisit.Rdata"), envir=e)
load(file.path(ROOT, "data", "data-wrsi.Rdata"), envir=e)
TAX <- droplevels(e$TAX)
TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")

en <- new.env()
load(file.path(ROOT, "data", "data-josmshf.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY
tax <- droplevels(TAX[colnames(yyn),])

#rm(e, en)

## model for species
fln <- list.files(file.path(ROOT, "results", "josmshf"))
fln <- sub("birds_abmi-josmshf_", "", fln)
fln <- sub(".Rdata", "", fln)

SPP <- colnames(yys)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))

names(modsn)
stage_hab_n <- 5

## spp specific output

spp <- "BTNW"

## all N & S data for plots
res_coef <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    res_coef[[spp]] <- list(veg=NULL, max=NA)
    resn <- loadSPP(file.path(ROOT, "results", "north",
        paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    tmp1 <- pred_veghf(estn_hab, Xnn, burn_included=FALSE)
    res_coef[[spp]]$veg <- tmp1
    res_coef[[spp]]$max[1] <- fig_veghf_ymax(tmp1)
}
#save(res_coef, file=file.path(ROOT, "tables", "res_coef.Rdata"))

## get all kinds of linear coefs

get_lin <- function(tp) {
    pr <- unname(attr(tp$veg, "linear"))
    c(AverageCoef=pr[1],
        SoftLin10=pr[1] * exp(0.1*log(pr[2])),
        HardLin10=pr[1] * pr[5],
        SoftLin	=pr[2],
        SoftLin.LCI	=pr[3],
        SoftLin.UCI	=pr[4],
        HardLin	=pr[5],
        HardLin.LCI	=pr[6],
        HardLin.UCI=pr[7])
}

lin1_n <- list()
for (spp in rownames(tax)) {
    lin1_n[[spp]] <- get_lin(res_coef[[spp]])
}

## plots for N & S

fname <- file.path(ROOT, "josmshf", "habitat-assoc.pdf")
pdf(file=fname,width=16,height=7,onefile=TRUE)
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    tp <- res_coef[[spp]]
    NAM <- as.character(tax[spp, "English_Name"])
    if (max(tp$max, na.rm=TRUE) > 3*min(tp$max, na.rm=TRUE)) {
        MAXn <- tp$max[1]
    } else {
        MAXn <- max(tp$max, na.rm=TRUE)
    }
    NDAT <- sum(yyn[,spp] > 0)
    fig_veghf(tp$veg, paste0(NAM, " (n = ", NDAT, " detections)"), ymax=MAXn)
}
dev.off()

## surroundinghf-north

fname <- file.path(ROOT, "josmshf", "surrounding-hf.pdf")
pdf(file=fname, width=7.5, height=5.7, onefile=TRUE)
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(tax[spp, "English_Name"])
    resn <- loadSPP(file.path(ROOT, "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    estn_sp <- getEst(resn, stage=stage_hab_n + 2, na.out=FALSE, Xnn)
    op <- par(mai=c(0.9,1,0.2,0.3))
    fig_hf_noremn(estn_sp, Xnn, LAB=paste0(NAM, ", North"))
    par(op)
}
dev.off()

## GoF AUC ROC R2 and stuff

library(pROC)
library(ResourceSelection)
#en <- new.env()
#load(file.path(ROOT, "data", "data-north.Rdata"), envir=en)
OFFn <- en$OFF
OFFmn <- en$OFFmean
BBn <- en$BB
DATn <- en$DAT
#rm(en)
## out-of-sample set
bunique <- unique(BBn)
INTERNAL <- 1:nrow(DATn) %in% bunique
ss <- which(INTERNAL)
ss1 <- which(!INTERNAL)
bid <- DATn$bootid
levels(bid) <- sapply(strsplit(levels(bid), "\\."), "[[", 1)

pr_fun_for_gof <- function(est, X, off=0) {
    if (is.null(dim(est))) {
        mu0 <- drop(X %*% est)
        exp(mu0 + off)
    } else {
        mu0 <- apply(est, 1, function(z) X %*% z)
        rowMeans(exp(mu0 + off))
    }
}

## road effects
roadside_bias <- list()
Xnn_onroad <- Xnn_offroad <- Xnn[Xnn[,"ROAD01"] == 1,]
Xnn_offroad[,c("ROAD01","habCl:ROAD01")] <- 0
for (spp in rownames(tax)) {
cat(spp, "\n");flush.console()
off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
off1sp <- off1sp[Xnn[,"ROAD01"] == 1]
resn <- loadSPP(file.path(ROOT, "results", "josmshf",
    paste0("birds_abmi-josmshf_", spp, ".Rdata")))
est7 <- getEst(resn, stage=7, na.out=FALSE, Xnn)
pr_on <- pr_fun_for_gof(est7, Xnn_onroad, off=off1sp)
pr_off <- pr_fun_for_gof(est7, Xnn_offroad, off=off1sp)
roadside_bias[[spp]] <- c(on=mean(pr_on), off=mean(pr_off), onoff=mean(pr_on)/mean(pr_off))
}
roadside_bias <- do.call(rbind, roadside_bias)
save(roadside_bias, file=file.path(ROOT, "josmshf", "roadside_bias.Rdata"))

## roadside avoidance index
rai_data <- data.frame(HAB=xnn$hab1ec, ROAD=xnn$ROAD01)
rai_pred <- matrix(0, nrow(Xnn), nrow(tax))
rownames(rai_pred) <- rownames(rai_data) <- rownames(xnn)
colnames(rai_pred) <- rownames(tax)
for (spp in rownames(tax)) {
cat(spp, "\n");flush.console()
off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
resn <- loadSPP(file.path(ROOT, "results", "josmshf",
    paste0("birds_abmi-josmshf_", spp, ".Rdata")))
est7 <- getEst(resn, stage=7, na.out=FALSE, Xnn)
rai_pred[,spp] <- pr_fun_for_gof(est7, Xnn, off=off1sp)
}
save(rai_pred, rai_data, file=file.path(ROOT, "josmshf", "roadside_avoidance.Rdata"))


spp <- "OVEN"
all_acc <- list()
for (spp in rownames(tax)[colSums(yyn[ss1,]>0) > 0]) {
cat(spp, "\n");flush.console()

## Deviance based pseudo R^2 (only internal)
y1sp <- yyn[,spp]
y10sp <- ifelse(y1sp > 0, 1, 0)
off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
resn <- loadSPP(file.path(ROOT, "results", "josmshf",
    paste0("birds_abmi-josmshf_", spp, ".Rdata")))
est5 <- getEst(resn, stage=5, na.out=FALSE, Xnn)
est6 <- getEst(resn, stage=6, na.out=FALSE, Xnn)
est7 <- getEst(resn, stage=7, na.out=FALSE, Xnn)
## null model need to represent ARU as well
m0 <- glm(y1sp[ss] ~ ARU3, xnn[ss,], family=poisson, offset=off1sp[ss])
pr0 <- pr_fun_for_gof(coef(m0), model.matrix(~ARU3, xnn), off=off1sp)
## 1st run
pr11 <- pr_fun_for_gof(est5[1,], Xnn, off=off1sp)
pr12 <- pr_fun_for_gof(est6[1,], Xnn, off=off1sp)
pr13 <- pr_fun_for_gof(est7[1,], Xnn, off=off1sp)
## smooth
prf1 <- pr_fun_for_gof(est5, Xnn, off=off1sp)
prf2 <- pr_fun_for_gof(est6, Xnn, off=off1sp)
prf3 <- pr_fun_for_gof(est7, Xnn, off=off1sp)
## Null model: intercept and offset
ll0 <- sum(dpois(y1sp[ss], pr0[ss], log=TRUE))
## Saturated: one parameter per observation
lls <- sum(dpois(y1sp[ss], y1sp[ss], log=TRUE))
## Full: our smoothed prediction
llf1 <- sum(dpois(y1sp[ss], prf1[ss], log=TRUE))
llf2 <- sum(dpois(y1sp[ss], prf2[ss], log=TRUE))
llf3 <- sum(dpois(y1sp[ss], prf3[ss], log=TRUE))
## Full: our 1st prediction
ll11 <- sum(dpois(y1sp[ss], pr11[ss], log=TRUE))
ll12 <- sum(dpois(y1sp[ss], pr12[ss], log=TRUE))
ll13 <- sum(dpois(y1sp[ss], pr13[ss], log=TRUE))
R2f1 <- 1-(lls-llf1)/(lls-ll0)
R2f2 <- 1-(lls-llf2)/(lls-ll0)
R2f3 <- 1-(lls-llf3)/(lls-ll0)
R211 <- 1-(lls-ll11)/(lls-ll0)
R212 <- 1-(lls-ll12)/(lls-ll0)
R213 <- 1-(lls-ll13)/(lls-ll0)

## ROC (only external, and 0/1)
roc0 <- roc(y10sp[ss1], pr0[ss1])
rocf1 <- roc(y10sp[ss1], prf1[ss1])
rocf2 <- roc(y10sp[ss1], prf2[ss1])
rocf3 <- roc(y10sp[ss1], prf3[ss1])
roc11 <- roc(y10sp[ss1], pr11[ss1])
roc12 <- roc(y10sp[ss1], pr12[ss1])
roc13 <- roc(y10sp[ss1], pr13[ss1])
auc0 <- as.numeric(roc0$auc)
aucf1 <- as.numeric(rocf1$auc)
aucf2 <- as.numeric(rocf2$auc)
aucf3 <- as.numeric(rocf3$auc)
auc11 <- as.numeric(roc11$auc)
auc12 <- as.numeric(roc12$auc)
auc13 <- as.numeric(roc13$auc)
kf1 <- (aucf1-auc0) / (1-auc0)
kf2 <- (aucf2-auc0) / (1-auc0)
kf3 <- (aucf3-auc0) / (1-auc0)
k11 <- (auc11-auc0) / (1-auc0)
k12 <- (auc12-auc0) / (1-auc0)
k13 <- (auc13-auc0) / (1-auc0)

spp1res <- c(R2f1=R2f1, R2f2=R2f2, R2f3=R2f3, R211=R211, R212=R212, R213=R213,
    AUC0=auc0, AUCf1=aucf1, AUCf2=aucf2, AUCf3=aucf3, AUC11=auc11, AUC12=auc12, AUC13=auc13,
    kf1=kf1, kf2=kf2, kf3=kf3, k11=k11, k12=k12, k13=k13,
    y_in=mean(y1sp[ss]),
    lam0_in=mean(pr0[ss]),
    lam1_in=mean(prf1[ss]),
    lam2_in=mean(prf2[ss]), lam3_in=mean(prf3[ss]),
    y_out=mean(y1sp[ss1]),
    lam0_out=mean(pr0[ss1]),
    lam1_out=mean(prf1[ss1]),
    lam2_out=mean(prf2[ss1]), lam3_out=mean(prf3[ss1]))

## regional ROC analysis
regres <- list()
#REG <- "Foothills"
for (REG in levels(bid)) {
    regss <- bid == REG & !INTERNAL
    yreg <- y10sp[regss]
    if (sum(yreg) == 0) {
        reg1res <- c(AUC0=NA, AUCf1=NA, AUCf2=NA, AUCf3=NA, AUC11=NA, AUC12=NA, AUC13=NA,
            kf1=NA, kf2=NA, kf3=NA, k11=NA, k12=NA, k13=NA)
    } else {
        roc0 <- roc(yreg, pr0[regss])
        rocf1 <- roc(yreg, prf1[regss])
        rocf2 <- roc(yreg, prf2[regss])
        rocf3 <- roc(yreg, prf3[regss])
        roc11 <- roc(yreg, pr11[regss])
        roc12 <- roc(yreg, pr12[regss])
        roc13 <- roc(yreg, pr13[regss])
        auc0 <- as.numeric(roc0$auc)
        aucf1 <- as.numeric(rocf1$auc)
        aucf2 <- as.numeric(rocf2$auc)
        aucf3 <- as.numeric(rocf3$auc)
        auc11 <- as.numeric(roc11$auc)
        auc12 <- as.numeric(roc12$auc)
        auc13 <- as.numeric(roc13$auc)
        kf1 <- (aucf1-auc0) / (1-auc0)
        kf2 <- (aucf2-auc0) / (1-auc0)
        kf3 <- (aucf3-auc0) / (1-auc0)
        k11 <- (auc11-auc0) / (1-auc0)
        k12 <- (auc12-auc0) / (1-auc0)
        k13 <- (auc13-auc0) / (1-auc0)
        reg1res <- c(AUC0=auc0, AUCf1=aucf1, AUCf2=aucf2, AUCf3=aucf3,
            AUC11=auc11, AUC12=auc12, AUC13=auc13,
            kf1=kf1, kf2=kf2, kf3=kf3, k11=k11, k12=k12, k13=k13)
    }
    regres[[REG]] <- reg1res
}
regres <- do.call(rbind, regres)
regres <- regres[c("NW", "NE", "SW", "SE", "Foothills", "Parkland", "Rocky Mountain"),]

all_acc[[spp]] <- list(overall=spp1res, regions=regres)
}

## OCCC metrics
library(epiR)
occc_res <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    pr <- exp(pred_veghf(estn_hab, Xnn, burn_included=FALSE, raw=TRUE))
    occc_res[[spp]] <- epi.occc(pr)
}


## means
flam <- function(ss1) {
    c(lam_obs=mean(y1sp[ss1]),
    lam_est=mean(prf3[ss1]),
    p_obs=mean(y10sp[ss1]),
    p_est=mean(1-exp(-prf3[ss1])))
}
lam_all <- list()
for (spp in fln) {
    cat(spp, "\n");flush.console()
    y1sp <- yyn[,spp]
    y10sp <- ifelse(y1sp > 0, 1, 0)
    off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
    resn <- loadSPP(file.path(ROOT, "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    #est5 <- getEst(resn, stage=5, na.out=FALSE, Xnn)
    #est6 <- getEst(resn, stage=6, na.out=FALSE, Xnn)
    est7 <- getEst(resn, stage=7, na.out=FALSE, Xnn)
    #prf1 <- pr_fun_for_gof(est5, Xnn, off=off1sp)
    #prf2 <- pr_fun_for_gof(est6, Xnn, off=off1sp)
    prf3 <- pr_fun_for_gof(est7, Xnn, off=off1sp)

    ## regional ROC analysis
    regres <- list()
    #REG <- "Foothills"
    for (REG in levels(bid)) {
        regres[[REG]] <- flam(bid == REG & !INTERNAL)
    }
    regres <- do.call(rbind, regres)
    regres <- regres[c("NW", "NE", "SW", "SE", "Foothills", "Parkland", "Rocky Mountain"),]

    lam_all[[spp]] <- rbind(All=flam(ss1), regres)
}

flam2 <- function(ss1) {
    c(lam_obs=mean(y1sp[ss1]),
    lam_est=mean(prf2[ss1]),
    p_obs=mean(y10sp[ss1]),
    p_est=mean(1-exp(-prf2[ss1])))
}
lam_all6 <- list()
for (spp in fln) {
    cat(spp, "\n");flush.console()
    y1sp <- yyn[,spp]
    y10sp <- ifelse(y1sp > 0, 1, 0)
    off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
    resn <- loadSPP(file.path(ROOT, "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata")))
    #est5 <- getEst(resn, stage=5, na.out=FALSE, Xnn)
    est6 <- getEst(resn, stage=6, na.out=FALSE, Xnn)
    #est7 <- getEst(resn, stage=7, na.out=FALSE, Xnn)
    #prf1 <- pr_fun_for_gof(est5, Xnn, off=off1sp)
    prf2 <- pr_fun_for_gof(est6, Xnn, off=off1sp)
    #prf3 <- pr_fun_for_gof(est7, Xnn, off=off1sp)

    ## regional ROC analysis
    regres <- list()
    #REG <- "Foothills"
    for (REG in levels(bid)) {
        regres[[REG]] <- flam2(bid == REG & !INTERNAL)
    }
    regres <- do.call(rbind, regres)
    regres <- regres[c("NW", "NE", "SW", "SE", "Foothills", "Parkland", "Rocky Mountain"),]

    lam_all6[[spp]] <- rbind(All=flam2(ss1), regres)
}

save(all_acc, occc_res, lam_all, lam_all6, file=file.path(ROOT, "josmshf", "res_acc.Rdata"))

load(file=file.path(ROOT, "tables", "res_acc.Rdata"))

pdet <- sapply(fln, function(z) sum(yyn[ss,z]>0)/length(ss))
overall <- t(sapply(all_acc, function(z) z$overall))
kreg <- t(sapply(all_acc, function(z) z$regions[,"kf2"]))
occc <- t(sapply(occc_res, function(z) unlist(z[1:3])))
lam <- t(sapply(lam_all, function(z) z["All", 1:2]))

plot(overall[,c("R2f2", "R212")], xlim=c(0,1), ylim=c(0,1))
abline(0,1)

overall[overall[,"R2f2"]<0,]

plot(overall[,"R2f1"], overall[,"R2f3"], xlim=c(0,1), ylim=c(0,1), col=2,
    xlab="Pseudo R^2, Habitat", ylab="Pseudo R^2, Habitat+Space (+SHF)")
segments(x0=overall[,"R2f1"], y0=overall[,"R2f2"], y1=overall[,"R2f3"], col=4)
abline(0,1, lty=2, col="grey")

plot(overall[,"R2f2"]-overall[,"R2f1"],
    overall[,"R2f3"]-overall[,"R2f1"], xlim=c(0,0.3), ylim=c(0,0.3), col=2,
    xlab="Spece-Habitat", ylab="SHF-Space")
abline(0,1, lty=2, col="grey")

plot(lam, xlim=c(0,0.5), ylim=c(0,0.5),
    xlab="Mean Observed Count", ylab="Mean Predicted # Det.")
abline(0,1, lty=2, col="grey")

plot(pdet[rownames(overall)]*length(ss), overall[,"R2f3"])

plot(overall[,c("kf1", "kf3")], xlim=c(-1,1), ylim=c(-1,1), col=2,
    xlab="CAUC, Habitat", ylab="CAUC, Habitat+Space (+SHF)")
segments(x0=overall[,"kf1"], y0=overall[,"kf2"], y1=overall[,"kf3"], col=4)
abline(0,1, lty=2, col="grey")
abline(h=0,v=0, lty=3, col="grey")

plot(occc[rownames(overall),"occc"], overall[,"kf3"])
abline(h=0,v=0.5, lty=3, col="grey")


layout(matrix(c(1,1,1,2,1,1,1,3,1,1,1,4), 3, 4, byrow=TRUE))
plot(occc[,2], occc[,3], xlim=c(0,1), ylim=c(0,1),
    #cex=1+2*sqrt(pdet), col=as.integer(typ)+1,
    xlab="Overall Precision", ylab="Overall Accuracy")
abline(0,1, lty=2, col="grey")
abline(h=0.5, v=0.5, lty=2, col="grey")
hist(occc[,1], xlim=c(0,1), main="Overall Concordance")
hist(occc[,2], xlim=c(0,1), main="Overall Precision")
hist(occc[,3], xlim=c(0,1), main="Overall Accuracy")

occcx <- occc[rowSums(is.na(occc))==0,]
occcx[occcx[,2] < 0.2 & occcx[,3] < 0.2,,drop=FALSE]
occcx[occcx[,2] < 0.4 & occcx[,3] > 0.6,,drop=FALSE]
occcx[occcx[,2] > 0.8 & occcx[,3] < 0.2,,drop=FALSE]

foccc <- function(spp, nn=10) {
    resn <- loadSPP(file.path(ROOT, "results", "north",
        paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    pr <- exp(pred_veghf(estn_hab, Xnn, burn_included=FALSE, raw=TRUE))
    ii <- sample(ncol(pr), min(nn, ncol(pr)))

    prr <- pr[order(rowMeans(pr)),ii]
    prr <- t(t(prr)/apply(prr,2,max))
    matplot(prr, type="l", col=sample(rainbow(length(ii))), lty=1,
        axes=FALSE, xlab="Land cover types", ylab="Relative abundance",
        main=as.character(tax[spp,"English_Name"]), ylim=c(0,1.2))
    text(1, 1.1, paste("Overall Concordance =", round(occc_res[[spp]]$occc, 3),
        "\nOverall Precision =", round(occc_res[[spp]]$oprec, 3),
        "\nOverall Accuracy =", round(occc_res[[spp]]$oaccu, 3)), pos=4)
    box()
}


par(mfrow=c(2,3))
ymax <- max(pr0,prf1, prf2)
ResourceSelection:::.mep(y1sp[ss], pr0[ss], link="log", type="unique", level=0.9,
    main="Internal Null", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss], prf1[ss], link="log", type="unique", level=0.9,
    main="Internall Local", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss], prf2[ss], link="log", type="unique", level=0.9,
    main="Internal Space", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss1], pr0[ss1], link="log", type="unique", level=0.9,
    main="External Null", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss1], prf1[ss1], link="log", type="unique", level=0.9,
    main="Externall Local", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss1], prf2[ss1], link="log", type="unique", level=0.9,
    main="External Space", ylim=c(0,ymax))


## mean plots

REGNAMS <- substr(rownames(lam_all[[1]])[-1], 1, pmin(nchar(rownames(lam_all[[1]])[-1]), 5))
pdf(file.path(ROOT, "josmshf", "gof-measures.pdf"), width=12, height=6, onefile=TRUE)
for (spp in names(all_acc)) {
cat(spp, "\n");flush.console()
op <- par(las=1, mfrow=c(1,2))
MAX <- max(lam_all[[spp]][,1:2])*1.2
plot(lam_all[[spp]][-1,1:2], pch=21, ylim=c(0, MAX), xlim=c(0, MAX),
    col=4, cex=2, main=as.character(tax[spp,"English_Name"]),
    xlab="Mean Observed Count", ylab="Mean Expected Number of Detections")
abline(0,1)
abline(v=lam_all[[spp]][1,1], h=lam_all[[spp]][1,2], col=4)
text(lam_all[[spp]][-1,1], lam_all[[spp]][-1,2],
    round(all_acc[[spp]]$regions[,"kf3"], 2), pos=4, cex=0.6, col=4)
text(lam_all[[spp]][-1,1], lam_all[[spp]][-1,2], REGNAMS, pos=3, cex=0.4, col=4)

points(lam_all6[[spp]][-1,1:2], pch=21, col=2, cex=2)
abline(v=lam_all6[[spp]][1,1], h=lam_all6[[spp]][1,2], col=2)
text(lam_all6[[spp]][-1,1], lam_all6[[spp]][-1,2],
    round(all_acc[[spp]]$regions[,"kf2"], 2), pos=2, cex=0.6, col=2)
text(lam_all6[[spp]][-1,1], lam_all6[[spp]][-1,2], REGNAMS, pos=1, cex=0.4, col=2)

legend("bottomright", bty="n", col=c(4,2), lty=1, legend=c("Hab+Clim+SHF", "Hab+Clim"))
STAT <- all_acc[[spp]]$overall[c("R2f2", "R2f3", "kf2", "kf3")]
text(0, MAX*0.9,
    paste0("R^2(H+C) = ", max(0, round(STAT[1],3)),
    "\nR^2(H+C+S) = ", max(0, round(STAT[2],3)),
    "\nCAUC(H+C) = ", round(STAT[3],3),
    "\nCAUC(H+C+S) = ", round(STAT[4],3)), pos=4)
foccc(spp, 25)
par(op)
}
dev.off()


