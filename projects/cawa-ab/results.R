library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"
OUTDIR <- "c:/Users/Peter/Dropbox/josm/cawa-jeff/revision"

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
tax <- droplevels(e$TAX)
tax$Fn <- droplevels(tax$English_Name)
levels(tax$Fn) <- nameAlnum(levels(tax$Fn), capitalize="mixed", collapse="")
yy <- e$YYw
yy01 <- yy
yy01[yy01>0] <- 1

load(file.path(ROOT, "data", "data-cawa.Rdata"))
y_cawa <- YY[,"CAWA"]
off_cawa <- OFF[,"CAWA"]
rm(e, OFF)

## terms and design matrices
Terms <- getTerms(mods, "list")
Xn <- model.matrix(getTerms(mods, "formula"), DAT)
colnames(Xn) <- fixNames(colnames(Xn))

names(mods)
# [1] "Hab"      "Age"      "CC"       "Contrast" "ARU"      "Space"
# [7] "Wet"      "Dec"      "HF"       "Year"
stage_hab <- 5
spp <- "CAWA"
NAM <- as.character(tax[spp, "English_Name"])

load(file.path(ROOT, "results", "cawa", "birds_abmi-cawa_CAWA.Rdata"))
#load(file.path(ROOT, "results", "north", "birds_abmi-north_CAWA.Rdata"))

(mid_summary <- getFancyMidTab(res, mods, truncate=1000))
write.csv(mid_summary, row.names=FALSE, file=file.path(OUTDIR, "cawa-tab-mid.csv"))
printCoefmat(coef_summary <- getSummary(res))
write.csv(coef_summary, row.names=TRUE, file=file.path(OUTDIR, "cawa-tab-coef.csv"))
pdf(file.path(OUTDIR, "cawa-fig-mid.pdf"))
plotMid(res, mods)
dev.off()

est_hab <- getEst(res, stage=stage_hab, na.out=FALSE, Xn)
res_coef <- pred_veghf(est_hab, Xn, burn_included=FALSE)
MAX <- fig_veghf_ymax(res_coef)
NDAT <- sum(y_cawa > 0)
#fname <- file.path(ROOT, "figs", "veghf-north",
#    paste0(as.character(tax[spp, "Spp"]), ".png"))
#png(file=fname,width=1500,height=700)
pdf(file.path(OUTDIR, "cawa-fig-hab.pdf"), width=14)
fig_veghf(res_coef, "", ymax=MAX, "Density (males / ha)")
dev.off()
## linear
#fname <- file.path(ROOT, "figs", "linear-north",
#    paste0(as.character(tax[spp, "Spp"]), ".png"))
#png(file=fname,width=350,height=400)
pdf(file.path(OUTDIR, "cawa-fig-lin.pdf"))
fig_linear(attr(res_coef, "linear"), paste0(NAM, "\nNorth (n = ", NDAT, " det.)"))
dev.off()

## climate & surrounding hf tables, climate surface maps

est_hf <- getEst(res, stage=length(mods)-1, na.out=FALSE, Xn)

pdf(file.path(OUTDIR, "cawa-fig-shf.pdf"), width=10, height=10)
par(mfrow=c(2,2))
fig_hf_noremn(est_hf, Xn, LAB=paste0(NAM, ", North, 0% Dec/Mix"),
    fillin=0, remn=NULL)
fig_hf_noremn(est_hf, Xn, LAB=paste0(NAM, ", North, 38% Dec/Mix"),
    fillin=0.38, remn="DecMixKM")
fig_hf_noremn(est_hf, Xn, LAB=paste0(NAM, ", North, 100% Dec/Mix"),
    fillin=1, remn="DecMixKM")
fig_any("DecMixKM", est_hf, Xn, "DecMixKM")
dev.off()

pdf(file.path(OUTDIR, "cawa-fig-dec.pdf"))
fig_any("DecMixKM", est_hf, Xn, LAB="", xlab="% deciduous and mixed forest in 564m buffer")
dev.off()

est_yr <- getEst(res, stage=length(mods), na.out=FALSE, Xn)

## annual % trend
fstat(100*(exp(0.1*est_yr[,"YR"])-1))
## roadside bias
fstat(exp(est_yr[,"ROAD01"]))
## EDR ratio
fstat(sqrt(exp(est_yr[,"ARU3RF"])))
fstat(sqrt(exp(est_yr[,"ARU3SM"])))

## residual trend

mu_hf <- apply(est_hf, 1, function(z) Xn %*% z)
lam_hf <- exp(mu_hf + off_cawa)
lam_hf_mean <- rowMeans(lam_hf)
#pr_hf_stat <- fstatv(pr_hf, level=level)

c_fun <- function(i) {
    yy <- y_cawa[BB[,i]]
    offf <- lam_hf[BB[,i], i]
    yr <- DAT$YR[BB[,i]]
    m <- glm(yy ~ yr, offset=offf, family="poisson")
    coef(m)
}
tr_coef <- pbsapply(1:240, c_fun)
c_fun2 <- function(i) {
    yy <- y_cawa[BB[,i]]
    offf <- lam_hf_mean[BB[,i]]
    yr <- DAT$YR[BB[,i]]
    m <- glm(yy ~ yr, offset=offf, family="poisson")
    coef(m)
}
tr_coef2 <- pbsapply(1:240, c_fun2)

DAT$isBBS <- DAT$PCODE=="BBSAB"
c_fun3 <- function(i, bbs=FALSE) {
    isBBS <- if (bbs)
        DAT$isBBS[BB[,i]] else !DAT$isBBS[BB[,i]]
    yy <- y_cawa[BB[,i]][isBBS]
    offf <- lam_hf_mean[BB[,i]][isBBS]
    yr <- DAT$YR[BB[,i]][isBBS]
    m <- glm(yy ~ yr, offset=offf, family="poisson")
    coef(m)
}
tr_coef3bbs <- pbsapply(1:240, c_fun3, bbs=TRUE)
tr_coef3not <- pbsapply(1:240, c_fun3, bbs=FALSE)

hist(100*(exp(0.1*est_yr[,"YR"])-1))
hist(100*(exp(0.1*tr_coef[2,])-1))
hist(100*(exp(0.1*tr_coef2[2,])-1))
hist(100*(exp(0.1*tr_coef3bbs[2,])-1))
hist(100*(exp(0.1*tr_coef3not[2,])-1))

fstat(100*(exp(0.1*est_yr[,"YR"])-1))
fstat(100*(exp(0.1*tr_coef[2,])-1))
fstat(100*(exp(0.1*tr_coef2[2,])-1))
fstat(100*(exp(0.1*tr_coef3bbs[2,])-1))
fstat(100*(exp(0.1*tr_coef3not[2,])-1))

100*sum(lam_hf_mean[DAT$isBBS])/sum(lam_hf_mean)

## annual indices?

## GoF AUC/ROC

library(pROC)
Y1 <- ifelse(y_cawa>0, 1, 0)

## out-of-sample set
ss1 <- which(!(1:nrow(DAT) %in% BB[,1]))
#ss <- lapply(1:240, function(i) which(!(1:nrow(DAT) %in% BB[,i])))
#ssd <- data.frame(id=unlist(ss), iter=unlist(lapply(1:240, function(i) rep(i, length(ss[[i]])))))
#ssx <- Xtab(~iter + id, ssd)
#ssi <- which(colSums(ssx) == 240)
#ssc <- as.integer(colnames(ssx)[ssi])
#compare_sets(ssc, ss1)

pr_fun_for_gof <- function(stage, off=0) {
    est0 <- getEst(res, stage=stage, na.out=FALSE, Xn)
    mu0 <- apply(est0, 1, function(z) Xn %*% z)
    lam0 <- exp(mu0 + off)
    rowMeans(lam0)
}
prAll <- sapply(0:10, pr_fun_for_gof, off=off_cawa)
#prNoOff <- pr_fun_for_gof(0, off=0)

## external only
rocAll1 <- pblapply(1:ncol(prAll), function(i) roc(Y1[ss1], prAll[ss1,i]))
names(rocAll1) <- c("NULL", names(mods))
auc <- sapply(rocAll1, function(z) as.numeric(z$auc))
#rocNoOff <- roc(Y1[ss1], prNoOff[ss1])

#Y2 <- ifelse(y_cawa>1, 1, 0)
#rocAll2 <- lapply(1:10, function(i) roc(Y2, prAll[,i]))

pdf(file.path(OUTDIR, "cawa-fig-roc.pdf"))
Show <- rocAll1[c("Wet","HF","Year")]
#Show <- rocAll1
#Col <- rev(brewer.pal(length(Show), "RdBu"))
Col <- c("blue","lightblue","red")
op <- par(las=1)
plot(Show[[1]], col=Col[1], lty=1)
for (i in 2:length(Show))
    lines(Show[[i]], col=Col[i], lty=1)
aucx <- sapply(Show, function(z) as.numeric(z$auc))
#txt <- paste0(names(aucx), " (AUC = ", round(aucx, 3), ")")
txt <- paste0(c("Local","Spatial","Year"), " (AUC = ", round(aucx, 3), ")")
legend("bottomright", bty="n", col=rev(Col),
    lty=1, lwd=2, legend=rev(txt))
dev.off()



## CTI/Wet figure

est_wet <- est_hf[,c("xCTI","WetPT","WetPT:xCTI")]

#z$xCTI <- log((x$CTI + 1) / 10)
summary(10 * exp(DAT$xCTI) - 1)
vals1 <- seq(0, 1, len=100) # WetPT
vals2 <- seq(7, 22, len=100) # CTI
dat_wet <- expand.grid(WetPT=vals1, CTI=vals2)
dat_wet$xCTI <- log((dat_wet$CTI + 1) / 10)
X_wet <- model.matrix(~ xCTI*WetPT-1, dat_wet)

pr_wet <- exp(apply(est_wet, 1, function(z) X_wet %*% z))
pr_wet_stat <- fstatv(pr_wet, level=level)

pr_mat <- matrix(pr_wet_stat$Mean, length(vals1), length(vals2))


pdf(file.path(OUTDIR, "cawa-fig-wet.pdf"))
op <- par(las=1)
image(vals1*100, vals2, pr_mat,
      col = rev(grey(seq(0.2, 1, len=25))), axes=FALSE,
      ylab="Compound topographic index", xlab="% wet habitats in 150m buffer",
      main="")
cti <- 10 * exp(DAT$xCTI) - 1
ii <- sample(which(y_cawa==0), 10000)
with(DAT[ii,], points(100*WetPT, cti[ii], pch=21, cex=0.6, col="lightblue"))
with(DAT[y_cawa>0,], points(100*WetPT, cti[y_cawa>0], pch=19, cex=0.6, col="darkblue"))
contour(vals1*100, vals2, pr_mat, add=TRUE, labcex = 0.8, levels=c(0.1, 0.25, 0.5, 1))
axis(1)
axis(2)
box()
par(op)
dev.off()

## SRA/EDR results

library(MASS)
library(QPAD)
ROOTx <- "c:/bam/May2015"
load(file.path(ROOTx, "out", "BAMCOEFS_QPAD_v3.rda"))
.BAMCOEFS$version

load(file.path(ROOT, "data", "data-offset-covars.Rdata"))
summary(offdat)

Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
    LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))

spp <- "CAWA"
p <- rep(NA, nrow(offdat))
A <- q <- p

## constant for NA cases
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
## best model
mi <- bestmodelBAMspecies(spp, type="BIC", model.sra=0:8)
cat(spp, unlist(mi), "\n");flush.console()
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
#vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0

p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])
unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * offdat$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))

phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1, pi * offdat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))

ii <- which(p == 0)
p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])

offdat$p <- p
offdat$A <- A
offdat$q <- q
offdat <- offdat[rownames(DAT),]
summary(offdat)
summary(offdat$p)
summary(100*sqrt(offdat$A[offdat$MAXDIS==Inf]/pi))
summary(offdat$A[offdat$MAXDIS<Inf])

R <- 200
#spp <- "OVEN"
level <- 0.9
version <- 3
prob <- c(0, 1) + c(1, -1) * ((1-level)/2)

jd <- seq(0.35, 0.55, 0.01) # TSSR
ts <- seq(-0.25, 0.5, 0.01) # JDAY
ls <- seq(0, 0.25, len=length(jd)) # DSLS

xp1 <- expand.grid(JDAY=jd, # ---------- CHECK !!!
    TSSR=ts)
xp1$JDAY2 <- xp1$JDAY^2
xp1$TSSR2 <- xp1$TSSR^2
xp1$Jday <- xp1$JDAY * 365
xp1$Tssr <- xp1$TSSR * 24

xp2 <- expand.grid(DSLS=ls, # ---------- CHECK !!!
    TSSR=ts)
xp2$DSLS2 <- xp2$DSLS^2
xp2$TSSR2 <- xp2$TSSR^2
xp2$Dsls <- xp2$DSLS * 365
xp2$Tssr <- xp2$TSSR * 24

Xp1 <- model.matrix(~., xp1)
#colnames(Xp1)[1] <- "INTERCEPT"
Xp2 <- model.matrix(~., xp2)
#colnames(Xp2)[1] <- "INTERCEPT"

if (getBAMversion() < 3) {
    lc <- seq(1, 5, 1)
    tr <- seq(0, 1, 0.01)
    xq <- expand.grid(LCC=as.factor(lc),
        TREE=tr)
} else {
#    lc2 <- factor(c("Forest", "OpenWet"), c("Forest", "OpenWet"))
    lc <- factor(c("DecidMixed", "Conif", "Open", "Wet"),
        c("DecidMixed", "Conif", "Open", "Wet"))
    tr <- seq(0, 1, 0.01)
    xq <- expand.grid(LCC4=lc, TREE=tr)
    xq$LCC2 <- as.character(xq$LCC4)
    xq$LCC2[xq$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
    xq$LCC2[xq$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
    xq$LCC2 <- factor(xq$LCC2, c("Forest", "OpenWet"))
}
Xq0 <- model.matrix(~., xq)

SPP <- getBAMspecieslist()
cfall <- exp(t(sapply(SPP, function(spp)
    unlist(coefBAMspecies(spp, 0, 0)))))
t <- seq(0, 10, 0.1)
r <- seq(0, 4, 0.05)
pp <- sapply(SPP, function(spp) sra_fun(t, cfall[spp,1]))
qq <- sapply(SPP, function(spp) edr_fun(r, cfall[spp,2]))

spp <- "CAWA"

## model weights
wp <- selectmodelBAMspecies(spp)$sra$wBIC
wq <- selectmodelBAMspecies(spp)$edr$wBIC
names(wp) <- rownames(selectmodelBAMspecies(spp)$sra)
names(wq) <- rownames(selectmodelBAMspecies(spp)$edr)
nsra <- selectmodelBAMspecies(spp)$sra$nobs[1]
nedr <- selectmodelBAMspecies(spp)$edr$nobs[1]

## constant model
#lphi0 <- coefBAMspecies(spp, 0, 0)$sra
#lphise0 <- sqrt(vcovBAMspecies(spp, 0, 0)$sra[1])
#lphipi0 <- quantile(rnorm(R, lphi0, lphise0), prob)
#p <- cbind(est=sra_fun(t, exp(lphi0)))#,
#           pi1=sra_fun(t, exp(lphipi0[1])),
#           pi2=sra_fun(t, exp(lphipi0[2])))
#ltau0 <- coefBAMspecies(spp, 0, 0)$edr
#ltause0 <- sqrt(vcovBAMspecies(spp, 0, 0)$edr[1])
#ltaupi0 <- quantile(rnorm(R, ltau0, ltause0), prob)
#q <- cbind(est=edr_fun(r, exp(ltau0)))#,
#           pi1=edr_fun(r, exp(ltaupi0[1])),
#           pi2=edr_fun(r, exp(ltaupi0[2])))

## covariate effects
mi <- bestmodelBAMspecies(spp, type="BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

Xp <- if (mi$sra %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
Xp <- Xp[,names(cfi$sra),drop=FALSE]
lphi1 <- drop(Xp %*% cfi$sra)
#pcf1 <- mvrnorm(R, cfi$sra, vci$sra)
#lphipi1 <- t(apply(apply(pcf1, 1, function(z) drop(Xp %*% z)),
#    1, quantile, prob))
#px <- cbind(est=exp(lphi1), exp(lphipi1))
pmat <- matrix(exp(lphi1), length(jd), length(ts))
pmax <- sra_fun(10, max(exp(lphi1)))
pmat <- sra_fun(3, pmat)
pmax <- 1

Xq <- Xq0[,names(cfi$edr),drop=FALSE]
ltau1 <- drop(Xq %*% cfi$edr)
#pcf1 <- mvrnorm(R, cfi$edr, vci$edr)
#ltaupi1 <- t(apply(apply(pcf1, 1, function(z) drop(Xq %*% z)),
#    1, quantile, prob))
#qx <- cbind(est=exp(ltau1), exp(ltaupi1))
qmat <- matrix(exp(ltau1), length(lc), length(tr))
qmax <- edr_fun(0.5, max(exp(ltau1)))
qmat <- edr_fun(1, qmat)
qmax <- 1

#library(lattice)
#levelplot(px[,1] ~ Jday * Tssr, xp)

op <- par(las=1, mfrow=c(3,2))

barplot(wp, space=0, col=grey(1-wp), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", nsra, ") v", getBAMversion()),
    ylab="Model weight", xlab="Model ID")
barplot(wq, space=0, col=grey(1-wq), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", nedr, ") v", getBAMversion()),
    ylab="Model weight", xlab="Model ID")

plot(t, pp[,spp], type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
matlines(t, pp, col="grey", lwd=1, lty=1)
lines(t, pp[,spp], col=1, lwd=2)

plot(r*100, qq[,spp], type="n", ylim=c(0,1),
     xlab="Point count radius (m)",
     ylab="Probability of detection")
#polygon(100*c(r, rev(r)), c(q[,2], rev(q[,3])),
#        col="grey", border="grey")
matlines(r*100, qq, col="grey", lwd=1, lty=1)
lines(r*100, qq[,spp], col=1, lwd=2)
abline(v=cfall[spp,2]*100, lty=2)
rug(cfall[,2]*100, side=1, col="grey")
box()

xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main=paste("Best model:", mi$sra))
box()
image(1:nlevels(xq$LCC4), tr*100, qmat,
      col = rev(grey(seq(0, qmax, len=12))), axes=FALSE,
      xlab="Land cover types", ylab="Percent tree cover",
      main=paste("Best model:", mi$edr))
if (version < 3)
    axis(1, 1:5, c("DConif","DDecid","SConif","SDecid","Open"))
if (version > 2)
    axis(1, 1:nlevels(xq$LCC4), levels(xq$LCC4))
axis(2)
box()

par(op)


range(pmat)
range(qmat)

pdf(file.path(OUTDIR, "cawa-fig-qpad.pdf"), width=12)
op <- par(las=1, mfrow=c(1,2))

xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    #col = rev(grey(seq(1-max(pmat), 1-min(pmat), len=50))),
    col = rev(grey(seq(1-0.8, 1, len=25))),
    xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main="A")#"Probability of availability")
#contour(xval, ts*24, pmat, add=TRUE, labcex = 0.8)
box()

#ii <- sample(which(y_cawa==0), 1000)
#with(offdat[ii,], points(365*JDAY, 24*TSSR, pch=21, cex=0.8, col="lightblue"))
#with(offdat[y_cawa>0,], points(365*JDAY, 24*TSSR, pch=19, cex=0.8, col="darkblue"))

image(1:nlevels(xq$LCC4), tr*100, qmat,
      #col = rev(grey(seq(1-max(qmat), 1-min(qmat), len=50))), axes=FALSE,
      col = rev(grey(seq(1-0.8, 1, len=25))), axes=FALSE,
      xlab="Land cover types", ylab="Percent tree cover",
      main="B")#"Probability of detection")
if (version < 3)
    axis(1, 1:5, c("DConif","DDecid","SConif","SDecid","Open"))
if (version > 2)
    axis(1, 1:nlevels(xq$LCC4), levels(xq$LCC4))
axis(2)
box()
par(op)
dev.off()
