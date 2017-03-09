#devtools::install_github("psolymos/opticut")
library(opticut)
load("~/Dropbox/collaborations/opticut/R/abmi-aru-1min.Rdata")

y <- xt[,colnames(xt) != "DomesticCow"]
x$POINT <- x$SITE_LABEL
x0 <- x
x <- x[,c("PKEY", "POINT", "SITE",
    "YEAR", "MONTH", "MDAY", "HOUR", "MINUTE", "YDAY",
    "RAIN", "WIND", "INDUSTRY", "NOISE", "MICROPHONE")]


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
