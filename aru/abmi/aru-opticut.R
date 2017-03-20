#devtools::install_github("psolymos/opticut")
library(opticut)
library(mefa4)

## analyses for ms/scilett

data(birdrec)
x <- birdrec$x
y <- birdrec$y
table(colSums(y>0))
nmin <- 1
y <- y[,colSums(y>0) >= nmin]

## species accumulation

yo <- y[order(x$YDAY),]
xo <- x[order(x$YDAY),]

yall <- groupSums(yo, 1, xo$YDAY)
xall <- nonDuplicated(xo, xo$YDAY, TRUE)

ymor <- groupSums(yo[xo$TOD=="Morning",], 1, xo$YDAY[xo$TOD=="Morning"])
xmor <- nonDuplicated(xo[xo$TOD=="Morning",], xo$YDAY[xo$TOD=="Morning"], TRUE)

ymid <- groupSums(yo[xo$TOD=="Midnight",], 1, xo$YDAY[xo$TOD=="Midnight"])
xmid <- nonDuplicated(xo[xo$TOD=="Midnight",], xo$YDAY[xo$TOD=="Midnight"], TRUE)

yall <- apply(yall, 2, function(z) pmin(1, cumsum(z)))
ymor <- apply(ymor, 2, function(z) pmin(1, cumsum(z)))
ymid <- apply(ymid, 2, function(z) pmin(1, cumsum(z)))

d1all <- apply(yall, 2, function(z) xall$YDAY[which(z>0)[1]])
d1mor <- apply(ymor, 2, function(z) xmor$YDAY[which(z>0)[1]])
d1mid <- apply(ymid, 2, function(z) xmid$YDAY[which(z>0)[1]])

Sall <- rowSums(yall)
Smor <- rowSums(ymor)
Smid <- rowSums(ymid)
dall <- diff(Sall)/diff(xall$YDAY)
dmor <- diff(Smor)/diff(xmor$YDAY)
dmid <- diff(Smid)/diff(xmid$YDAY)

par(mfrow=c(2,1))
plot(xall$YDAY, Sall, type="l")
lines(xmor$YDAY, Smor, col=2)
lines(xmid$YDAY, Smid, col=4)

plot(lowess(dall ~ xall$YDAY[-1]), type="l")
lines(lowess(dmor ~ xmor$YDAY[-1]), col=2)
lines(lowess(dmid ~ xmid$YDAY[-1]), col=4)

plot(table(colSums(y>0)))

lall <- lowess(log(rowSums(y)+1) ~ x$YDAY)
lmor <- lowess(log(rowSums(y[x$TOD=="Morning",])+1) ~ x$YDAY[x$TOD=="Morning"])
lmid <- lowess(log(rowSums(y[x$TOD=="Midnight",])+1) ~ x$YDAY[x$TOD=="Midnight"])

plot(jitter(rowSums(y), factor=2) ~ x$YDAY, pch="+", cex=0.6, col="grey")
lines(lall$x, exp(lall$y)-1, col=1)
rug(d1all)
lines(lmor$x, exp(lmor$y)-1, col=2)
lines(lmid$x, exp(lmid$y)-1, col=4)


## spp detected at midnight only
colnames(y)[colSums(y01mor) == 0 & colSums(y01mid) > 0]
sum(colSums(y01mor) == 0 & colSums(y01mid) > 0)
sum(colSums(y01mor) > 0 & colSums(y01mid) == 0)
sum(colSums(y01mor) > 0 & colSums(y01mid) > 0)





## --

load("~/Dropbox/collaborations/opticut/R/abmi-aru-1min.Rdata")

SsnLab <- c("Early", "Mid", "Late")
x$Ssn <- factor(NA, SsnLab)
x$Ssn[x$ToYc %in% 1:3] <- SsnLab[1]
x$Ssn[x$ToYc %in% 4:7] <- SsnLab[2]
x$Ssn[x$ToYc %in% 8] <- SsnLab[3]

table(x$Ssn, x$ToYc)
table(x$Ssn, x$ToDc)

table(x$Ssn)
table(x$ToYc)
table(x$ToDc)

xtv <- ifelse(xt > 0, 1, 0)

cmin <- 10
#cs <- cumsum(table(colSums(xtv)))
#plot(as.integer(names(cs)), cs, type="l")
xtv <- xtv[,colSums(xtv) >= cmin]
xtv <- xtv[,colnames(xtv) != "DomesticCow"]

## binomial opticut

o1 <- opticut(xtv ~ 1, strata=x$Ssn, dist="binomial")
o1a <- opticut(xtv ~ 1, strata=x$Ssn, dist="binomial", comb="all")
o2 <- opticut(xtv ~ 1, strata=x$ToDc, dist="binomial")

## uncertainty

library(parallel)
cl <- makeCluster(3)
u1 <- uncertainty(o1, type="multi", B=99, cl=cl)
u2 <- uncertainty(o2, type="multi", B=99, cl=cl)
stopCluster(cl)

## visualization

ocoptions(theme=c("#7b3294", "#ffffbf", "#008837"))

pdf("aru-season-2.pdf", height=10, width=10)
plot(o1,sort=1,mar=c(4,10,3,3),cut=2,ylab="",xlab="Season",
    show_I=TRUE, show_S=FALSE, cex.axis=0.6, lower=0.25)
dev.off()

pdf("aru-time-2.pdf", height=10, width=10)
plot(o2,sort=1,mar=c(4,10,3,3),cut=2,ylab="",xlab="Time",show_I=FALSE, show_S=FALSE, cex.axis=0.6)
dev.off()

## week level clustering

x$COMBO <- factor(paste(x$Ssn, x$ToDc, sep=""),
    c("EarlyMorning", "MidMorning", "LateMorning", "EarlyMidnight", "MidMidnight", "LateMidnight"))
o3 <- opticut(xtv ~ 1, data=x, strata=x$COMBO, dist="binomial")

x$ISSUE <- x$WIND + x$RAIN + x$NOISE
o4 <- opticut(xtv ~ ISSUE,
    data=x, strata=x$ToYc, dist="binomial")
bp <- summary(o4)$bestpart

pdf("aru-seasontime-2.pdf", height=10, width=10)
plot(o3, mar=c(4,10,3,3), cut=-Inf, sort=1, ylab="",xlab="Weeks", cex.axis=0.6)
dev.off()

pdf("aru-weeks-2.pdf", height=10, width=10)
plot(o4,sort=1,mar=c(4,10,3,3),cut=2,ylab="",xlab="Weeks",
    show_I=TRUE, show_S=TRUE, cex.axis=0.6, lower=0.25)
dev.off()


pdf("aru-hclust-2.pdf")
plot(hclust(vegan::vegdist(t(bp),"jaccard"), "ward.D2"))
dev.off()

pdf("aru-heat-2.pdf", height=10, width=10)
heatmap(bp, scale="none", col=occolors()(3)[2:3],
    distfun=function(x) vegan::vegdist(x, "jaccard"),
    margins=c(3,8))
dev.off()

library(parallel)
cl <- makeCluster(3)
u3 <- uncertainty(o3, type="multi", B=99, cl=cl)
u4 <- uncertainty(o4, type="multi", B=99, cl=cl)
stopCluster(cl)

save(o1, o1a, o2, u1, u2, o3, o4, u4,
    file="~/Dropbox/collaborations/opticut/R/abmi-aru-ocresults-2.Rdata")
