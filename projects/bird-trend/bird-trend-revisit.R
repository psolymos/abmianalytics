library(mefa4)
library(intrval)
reg <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
spp <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(spp) <- spp$sppid
rf <- new.env()
load(file.path("e:/peter/AB_data_v2017",
    "data", "analysis", "species", "OUT_birdsrf_2017-04-26.Rdata"),
    envir=rf)
sm <- new.env()
load(file.path("e:/peter/AB_data_v2017",
    "data", "analysis", "species", "OUT_birdssm_2017-04-05.Rdata"),
    envir=sm)

ii <- intersect(rownames(rf$x), rownames(rf$xt))
x_rf <- data.frame(SITE=rf$x[ii,"SITE"],
    YEAR=rf$x[ii,"YEAR"],
    part="rf")
xt3_rf <- as.matrix(rf$xt3[["1"]][ii,])
#xt10_rf <- as.matrix(rf$xt[ii,])

x_sm <- data.frame(SITE=sm$x$SITE,
    YEAR=sm$x$YEAR,
    part="sm")
xt3_sm <- sm$xt

jj <- intersect(colnames(xt3_rf), colnames(xt3_sm))
compare_sets(colnames(xt3_rf), colnames(xt3_sm))


xt <- rbind(xt3_rf[,jj], xt3_sm[,jj])
x <- rbind(x_rf, x_sm)
x$SITE_YEAR <- paste0(x$SITE, "_", x$YEAR)
xx <- droplevels(nonDuplicated(x, SITE_YEAR, TRUE))
z <- as.matrix(Xtab(~SITE + YEAR, xx))
zz <- z[rowSums(z) > 1,]
table(rowSums(zz))
for (i in 1:nrow(zz)) {
    nn <- which(zz[i,] > 0)
    if (length(nn) > 2)
        zz[i,nn[1:(length(nn)-2)]] <- 0
}
table(rowSums(zz))

X0 <- XR <- matrix(NA, nrow(zz), ncol(xt))
dimnames(X0) <- dimnames(XR) <- list(rownames(zz), colnames(xt))
T <- numeric(nrow(zz))
names(T) <- rownames(zz)

for (i in rownames(zz)) {
    yrs <- as.integer(colnames(zz)[zz[i,] > 0])
    y0 <- xt[x$SITE == i & x$YEAR == yrs[1],]
    yR <- xt[x$SITE == i & x$YEAR == yrs[2],]
#    y0 <- ifelse(y0 > 0, 1, 0)
#    yR <- ifelse(yR > 0, 1, 0)
    X0[i,] <- colMeans(y0)
    XR[i,] <- colMeans(yR)
    T[i] <- diff(yrs)
}

lambda_fun <- function(ss=NULL) {
    if (is.null(ss))
        ss <- seq_len(nrow(zz))
    WM <- colSums(T[ss]*X0[ss,]) / colSums(X0[ss,])
    TC <- colSums(XR[ss,]) / colSums(X0[ss,])
    lam <- TC^(1/WM)
    100*(lam-1)
}

regz <- reg[rownames(zz),]
table(regz$NATURAL_REGIONS)

id <- which(regz$NATURAL_REGIONS %in% c("Boreal", "Canadian Shield",
    "Foothills", "Parkland"))

B <- 999
BB <- cbind(id, replicate(B, sample(id, length(id), replace=TRUE)))

res <- apply(BB, 2, lambda_fun)
res <- res[rowSums(is.na(res))==0,]
est <- cbind(Est=res[,1], t(apply(res, 1, quantile, prob=c(0.05, 0.95))))
est <- est[order(est[,1]),]
spp2 <- spp[rownames(est),]
table(0 %[]% est[,-1])

summary(est)

est <- est[!is.na(spp2$singing) & spp2$singing,] # RESQ is NA
spp3 <- spp2[rownames(est),]
lh <- read.csv("~/Dropbox/bam/lhreg2/LH-all-for-qpadv3.csv")
lh <- lh[match(spp3$scinam, lh$scientific_name),]
rownames(lh) <- rownames(spp3)
Spp <- rownames(spp3)[!is.na(lh$spp)]
est <- est[Spp,]
spp3 <- droplevels(spp3[Spp,])
lh <- droplevels(lh[Spp,])

xi <- seq_len(nrow(est))
ty <- sign(est[,1])
ty[0 %[]% est[,-1]] <- 0
COL <- c("red", "grey", "darkgreen")

plot(est[,1], xi, type="n", xlim=range(est), axes=FALSE, ann=FALSE)
abline(v=0, col="grey")
abline(v=c(-10,-5,5,10), lty=2, col="grey")
points(est[,1], xi, col=COL[ty+2], pch=19)
segments(y0=xi, y1=xi, x0=est[,2], x1=est[,3], col=COL[ty+2])
axis(1)
title(xlab="% annual change")
text(est[,1]-10, xi-0.1, rownames(est), pos=2, cex=0.6, col="darkgrey")


lh$Mig <- lh$Migr
levels(lh$Mig) <- c("LD","WR","WR","SD")
lh$Mig2 <- lh$Migr
levels(lh$Mig2) <- c("M","W","W","M")

lh$Hab4 <- lh$Hab3 <- lh$Hab2 <- lh$habitat
#levels(lh$Hab) <- c("Forest", "Grassland", "Lake/Pond", "Marsh", "Mountains",
#    "Open Woodland", "Scrub", "Shore-line", "Town")
levels(lh$Hab4) <- c("For", "Open", "Wet", "Wet", "Open",
    "Wood", "Open", "Wet", "Open")
levels(lh$Hab3) <- c("For", "Open", "Open", "Open", "Open",
    "Wood", "Open", "Open", "Open")
levels(lh$Hab2) <- c("Closed", "Open", "Open", "Open", "Open",
    "Closed", "Open", "Open", "Open")

par(mfrow=c(2,4))
boxplot(est[,1] ~ lh$Mig);abline(h=0)
boxplot(est[,1] ~ lh$food);abline(h=0)
boxplot(est[,1] ~ lh$Hab4);abline(h=0)
boxplot(est[,1] ~ lh$behavior);abline(h=0)
plot(est[,1] ~ lh$MaxFreqkHz);abline(h=0)
abline(lm(est[,1] ~ lh$MaxFreqkHz), col=4)
plot(est[,1] ~ lh$logmass);abline(h=0)
abline(lm(est[,1] ~ lh$logmass), col=4)
plot(est[,1] ~ lh$logphi);abline(h=0)
abline(lm(est[,1] ~ lh$logphi), col=4)
plot(est[,1] ~ lh$logtau);abline(h=0)
abline(lm(est[,1] ~ lh$logtau), col=4)

## todo
## check with other estimates

