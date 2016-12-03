#devtools::install_github("psolymos/opticut")
library(opticut)
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

pdf("aru-season.pdf", height=10, width=10)
plot(o1,sort=1,mar=c(4,10,3,3),cut=2,ylab="",xlab="Season",
    show_I=TRUE, show_S=FALSE, cex.axis=0.6, lower=0.25)
dev.off()

pdf("aru-time.pdf", height=10, width=10)
plot(o2,sort=1,mar=c(4,10,3,3),cut=2,ylab="",xlab="Time",show_I=FALSE, show_S=FALSE, cex.axis=0.6)
dev.off()


## week level clustering

o3 <- opticut(xtv ~ 1,
    data=x, strata=paste(x$Ssn, x$ToDc, sep=""), dist="binomial")
o4 <- opticut(xtv ~ WIND + RAIN + NOISE + INDUSTRY,
    data=x, strata=x$ToYc, dist="binomial")

plot(o3, mar=c(4,10,3,3),cut=2,ylab="",xlab="Weeks", cex.axis=0.6)

plot(o4,sort=1,mar=c(4,10,3,3),cut=2,ylab="",xlab="Weeks",
    show_I=TRUE, show_S=TRUE, cex.axis=0.6, lower=0.25)
bp <- summary(o4)$bestpart


plot(hclust(vegan::vegdist(t(bp),"jaccard"), "ward.D2"))

heatmap(bp, scale="none", col=occolors()(3)[2:3],
    distfun=function(x) dist(x, "manhattan"))

library(parallel)
cl <- makeCluster(3)
u3 <- uncertainty(o3, type="multi", B=99, cl=cl)
u4 <- uncertainty(o4, type="multi", B=99, cl=cl)
stopCluster(cl)

save(o1, o1a, o2, u1, u2, o3, o4,# u4,
    file="~/Dropbox/collaborations/opticut/R/abmi-aru-ocresults.Rdata")