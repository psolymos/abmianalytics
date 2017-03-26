#devtools::install_github("psolymos/opticut")
library(opticut)
library(mefa4)
library(vegan)

## analyses for ms/scilett

data(birdrec)
x <- birdrec$samp
x$toytod <- interaction(x$TOY, x$TOD)
x$xtoy <- x$TOY
levels(x$xtoy) <- gsub("[[:digit:]]", "", levels(x$xtoy))
x$toytod <- interaction(x$TOY, x$TOD)
x$xtoytod <- interaction(x$xtoy, x$TOD)
y <- birdrec$xtab
table(colSums(y>0))
nmin <- 1
y <- y[,colSums(y>0) >= nmin]
y01 <- ifelse(y > 0, 1, 0)
z <- read.csv("~/Dropbox/bam/lhreg2/aru_z.csv")
rownames(z) <- z$Species
levels(z$MigratoryBehaviour)[levels(z$MigratoryBehaviour) ==
    "Nomadic"] <- "Short distance migrant"

if (FALSE) {
names(birdrec) <- c("xtab", "samp")
birdrec$taxa <- z
birdrec$taxa$Species_ID <- NULL
save(birdrec, file="~/repos/opticut/data/birdrec.rda")
}

if (FALSE) {
t1 <- read.csv("~/Dropbox/bam/lhreg2/LH-all-for-qpadv3.csv")
t2 <- read.csv("~/Dropbox/bam/lhreg2/taxonomy.csv")
t3 <- read.csv("~/Dropbox/bam/lhreg2/tblAvianLifeHistory.csv", encoding="latin1")

z <- data.frame(Species=colnames(y))
rownames(z) <- z$Species

levels(t2$Fn)[levels(t2$Fn) == "BlackandwhiteWarbler"] <- "BlackAndWhiteWarbler"
compare_sets(rownames(z), t2$Fn)
setdiff(rownames(z), t2$Fn)
z$CommonName <- t2$English_Name[match(rownames(z), t2$Fn)]
z$ScientificName <- t2$Scientific_Name[match(rownames(z), t2$Fn)]
z$Family <- t2$Family_Sci[match(rownames(z), t2$Fn)]
z$Order <- t2$Order[match(rownames(z), t2$Fn)]

z$Species_ID <- t2$Species_ID[match(rownames(z), t2$Fn)]
compare_sets(z$Species_ID, t1$spp)
write.csv(z, row.names=FALSE, file="~/Dropbox/bam/lhreg2/aru_z2.csv")

zz <- birdrec$taxa
z0 <- z

#z <- read.csv("~/Dropbox/bam/lhreg2/aru_z.csv")
rownames(z) <- z$Species
z$Species_ID <- t2$Species_ID[match(rownames(z), t2$Fn)]
compare_sets(z$Species_ID, t3$SpeciesID)
z$MigratoryBehaviour <- t3$Migration1[match(z$Species_ID, t3$SpeciesID)]
#write.csv(z, row.names=FALSE, file="~/Dropbox/bam/lhreg2/aru_z.csv")
z$Class <- "Aves"
zzz <- z[!(rownames(z) %in% rownames(zz)),colnames(zz)]
write.csv(zzz, row.names=FALSE, file="~/Dropbox/bam/lhreg2/aru_znew.csv")

zzz <- read.csv("~/Dropbox/bam/lhreg2/aru_znew.csv")
rownames(zzz) <- zzz$Species

aa <- rbind(zz, zzz[,colnames(zz)])
aa <- aa[colnames(y),]
rownames(aa) <- colnames(y)
for (i in 1:ncol(aa))
    aa[,i] <- as.character(aa[,i])
aa[rownames(zzz),] <- zzz

birdrec$taxa <- aa
save(birdrec, file="~/repos/opticut/data/birdrec.rda")
}

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

## 1st day of detection
d1all <- apply(yall, 2, function(z) xall$YDAY[which(z>0)[1]])
d1mor <- apply(ymor, 2, function(z) xmor$YDAY[which(z>0)[1]])
d1mid <- apply(ymid, 2, function(z) xmid$YDAY[which(z>0)[1]])

Sall <- rowSums(yall)
Smor <- rowSums(ymor)
Smid <- rowSums(ymid)
dall <- diff(Sall)/diff(xall$YDAY)
dmor <- diff(Smor)/diff(xmor$YDAY)
dmid <- diff(Smid)/diff(xmid$YDAY)

lall <- lowess(log(rowSums(y)+1) ~ x$YDAY)
lmor <- lowess(log(rowSums(y[x$TOD=="Morning",])+1) ~ x$YDAY[x$TOD=="Morning"])
lmid <- lowess(log(rowSums(y[x$TOD=="Midnight",])+1) ~ x$YDAY[x$TOD=="Midnight"])


col <- occolors()(6)[c(2,5)]
par(mfrow=c(2,1), mar=c(4,4,1,1), las=1)
plot(xall$YDAY, Sall, type="n", xlab="", ylab="Number of Species", lwd=2,
    ylim=c(0,150), axes=FALSE)
lines(xmor$YDAY, Smor, col=col[1], lwd=2)
lines(xmid$YDAY, Smid, col=col[2], lwd=2)
abline(v=c(140, 180), col="grey")
box(bty="l", col="grey")
axis(1)
axis(2)
legend(85, 140, lty=1, lwd=2, col=col,
    legend=c("Morning", "Midnight"), bty="n", cex=0.8)
text(4:8*30-32, rep(150, 5), c("April", "May", "June", "July", "August"),
    cex=0.6, col="darkgrey")

plot(jitter(rowSums(y), factor=2.5) ~ jitter(x$YDAY, factor=2.5),
    pch=19, cex=1, col=rgb(0.75, 0.75, 0.75, 0.1), type="n",
    ylim=c(0,7), axes=FALSE,
    xlab="Day of Year", ylab="Detection rate / minute")
lines(lmor$x, exp(lmor$y)-1, col=col[1], lwd=2)
lines(lmid$x, exp(lmid$y)-1, col=col[2], lwd=2)
abline(v=c(140, 180), col="grey")
box(bty="l", col="grey")
axis(1)
axis(2)
text(4:8*30-32, rep(7, 5), c("April", "May", "June", "July", "August"),
    cex=0.6, col="darkgrey")


yd <- groupSums(y, 1, x$TOD)
z$sday <- factor("All", c("All", "Morning", "Midnight"))
z$sday[yd["Midnight",]==0] <- "Morning"
z$sday[yd["Morning",]==0] <- "Midnight"
table(z$sday)

dd <- apply(y, 2, function(z) rep(x$YDAY, z))

p <- c(0.05, 0.05, 0.5, 0.95, 0.95)
ddd <- t(sapply(dd, quantile, prob=p, na.rm=TRUE))
ddm <- sapply(dd, mean)
ddm <- ddm[order(ddd[,2])]
ddd <- ddd[order(ddd[,2]),]
col <- occolors()(6)[2:4]
split <- col[as.integer(z[rownames(ddd), "MigratoryBehaviour"])]
split2 <- c("grey", "white", "black")[as.integer(z[rownames(ddd), "sday"])]

par(mfrow=c(1,1), mar=c(4,4,1,1), las=1)
plot(0, type="n", axes=FALSE,
    xlim=c(75,210), ylim=c(nrow(ddd), -1), xlab="Day of Year", ylab="")
axis(1, at=seq(100, 200, 20))
for (i in 1:nrow(ddd)) {
    #lines(ddd[i,c(1,5)], c(i,i), col=split[i], lwd=1)
    polygon(ddd[i,c(2,4,4,2)], c(i-0.5,i-0.5,i+0.5,i+0.5), col=split[i], lwd=2, border=NA)
}
abline(v=c(140, 180), col="grey", lty=1)
for (i in 1:nrow(ddd)) {
    points(ddm[i], i, pch=19, cex=0.9, col=split2[i])
    points(ddm[i], i, pch=21, cex=0.9, col=1)
    text(ddd[i,2], i, z[rownames(ddd)[i], "CommonName"], cex=0.4, pos=2,
    col=1)
}
text(4:8*30-32, rep(-1, 5), c("April", "May", "June", "July", ""),
    cex=0.6, col="darkgrey")
legend(80, 135, fill=col, border=col, cex=0.6,
    bty="n", legend=levels(z$MigratoryBehaviour), title="Migratory Behaviour")
legend(110, 135, border=1, cex=0.6, fill=c("grey", "white", "black"),
    bty="n", legend=levels(z$sday), title="Diurnal Detections")

f <- function(z) c(Mean=mean(z), quantile(z, c(0.05, 0.95)))
f(ddd[z[rownames(ddd), "MigratoryBehaviour"]=="Resident",2])
f(ddd[z[rownames(ddd), "MigratoryBehaviour"]=="Short distance migrant",2])
f(ddd[z[rownames(ddd), "MigratoryBehaviour"]=="Neotropical migrant",2])



## spp detected at midnight only
colnames(y)[colSums(ymor) == 0 & colSums(ymid) > 0]
sum(colSums(ymor) == 0 & colSums(ymid) > 0)
sum(colSums(ymor) > 0 & colSums(ymid) == 0)
sum(colSums(ymor) > 0 & colSums(ymid) > 0)

#yy <- groupSums(y, 1, x$toytod)
yy <- groupSums(y, 1, x$xtoytod)
yy <- ifelse(yy > 0, 1, 0)
D <- vegdist(yy, "jaccard")
#D <- vegdist(yy, "bray")
#D <- dist(yy, "manhattan")
h <- hclust(D, "ward.D2")
plot(h)

oc <- opticut(y01 ~ 1, strata=x$xtoytod, dist="binomial")
bp <- summary(oc)$bestpart
D <- vegdist(t(bp), "jaccard")
h <- hclust(D, "ward.D2")
plot(h)

int <- range(x$YDAY)
plot(0, xlim=int, ylim=c(1,ncol(y)))
oo <- order(d1all)
for (i in 1:ncol(y01)) {
    tmp <- y01[,oo][,i]
    lines(range(x$YDAY[tmp > 0]), rep(i, 2), lwd=2, col="grey")
}
abline(v=c(140, 180))


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

CUT <- c(0, 105, 120, 140, 150, 160, 170, 180, 365)
d1 <- cbind(all=table(cut(d1all, CUT)),
    mor=table(cut(d1mor, CUT)),
    mid=table(cut(d1mid, CUT)))
barplot(t(d1), beside=TRUE)

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
