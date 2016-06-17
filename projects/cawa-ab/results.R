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
rm(e, OFF, BB)

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

(mid_summary <- getFancyMidTab(res, mods))
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
fig_veghf(res_coef, paste0(NAM, " (n = ", NDAT, " detections)"), ymax=MAX)
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

est_yr <- getEst(res, stage=length(mods), na.out=FALSE, Xn)

## annual % trend
fstat(100*(exp(0.1*est_yr[,"YR"])-1))
## roadside bias
fstat(exp(est_yr[,"ROAD01"]))
## EDR ratio
fstat(sqrt(exp(est_yr[,"ARU3RF"])))
fstat(sqrt(exp(est_yr[,"ARU3SM"])))

## CTI/Wet figure

est_wet <- est_hf[,c("xCTI","WetPT","WetPT:xCTI")]

#z$xCTI <- log((x$CTI + 1) / 10)
summary(10 * exp(DAT$xCTI) - 1)
vals1 <- seq(5, 25, len=100)
vals2 <- seq(0, 1, len=100)
dat_wet <- expand.grid(CTI=vals1, WetPT=vals2)
dat_wet$xCTI <- log((dat_wet$CTI + 1) / 10)
X_wet <- model.matrix(~ xCTI*WetPT-1, dat_wet)

pr_wet <- exp(apply(est_wet, 1, function(z) X_wet %*% z))
pr_wet_stat <- fstatv(pr_wet, level=level)

pr_mat <- matrix(pr_wet_stat$Mean, length(vals1), length(vals2))


pdf(file.path(OUTDIR, "cawa-fig-wet.pdf"))
image(vals1, vals2, pr_mat,
      col = rev(grey(seq(0.2, 1, len=25))), axes=FALSE,
      xlab="Compound topographic index", ylab="Point level wet %",
      main="Density")
cti <- 10 * exp(DAT$xCTI) - 1
ii <- sample(which(y_cawa==0), 10000)
with(DAT[ii,], points(cti[ii], WetPT, pch=21, cex=0.6, col="lightblue"))
with(DAT[y_cawa>0,], points(cti[y_cawa>0], WetPT, pch=19, cex=0.6, col="darkblue"))
contour(vals1, vals2, pr_mat, add=TRUE, labcex = 0.8, levels=c(0.1, 0.25, 0.5, 1))
box()
dev.off()

## SRA/EDR results

ROOT <- "c:/bam/May2015"

R <- 200
#spp <- "OVEN"
level <- 0.9
version <- 3
prob <- c(0, 1) + c(1, -1) * ((1-level)/2)

library(MASS)
library(QPAD)
load(file.path(ROOT, "out", "BAMCOEFS_QPAD_v3.rda"))
.BAMCOEFS$version

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


pdf(file.path(OUTDIR, "cawa-fig-qpad.pdf"), width=14)
op <- par(las=1, mfrow=c(1,2))
xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    #col = rev(grey(seq(1-max(pmat), 1-min(pmat), len=50))),
    col = rev(grey(seq(1-0.8, 1, len=25))),
    xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main="Probability of availability")
#contour(xval, ts*24, pmat, add=TRUE, labcex = 0.8)
box()

image(1:nlevels(xq$LCC4), tr*100, qmat,
      #col = rev(grey(seq(1-max(qmat), 1-min(qmat), len=50))), axes=FALSE,
      col = rev(grey(seq(1-0.8, 1, len=25))), axes=FALSE,
      xlab="Land cover types", ylab="Percent tree cover",
      main="Probability of detection")
if (version < 3)
    axis(1, 1:5, c("DConif","DDecid","SConif","SDecid","Open"))
if (version > 2)
    axis(1, 1:nlevels(xq$LCC4), levels(xq$LCC4))
axis(2)
box()
par(op)
dev.off()
