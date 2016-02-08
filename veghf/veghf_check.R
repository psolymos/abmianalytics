##% Processing backfilled veg + HF (cutblock ages incorporated)
##% P Solymos
##% April 28, 2015

SAVE <- TRUE

## root directory
ROOT <- "e:/peter"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2016"
## current year
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900

library(mefa4)
source("~/repos/abmianalytics/R/veghf_functions.R")
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type.csv")
hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY

### 1K grid --------------------------------------------------------

## dd1km_pred
load(file.path("c:/p", "AB_data_v2015", "out/kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))
dd1km_old <- dd1km_pred
load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))
## kgrid
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

p <- read.csv("c:/p/Parkland_NoSoilInfo.csv")
old <- as.matrix(dd1km_old$soil_current[p$LinkID,])
new <- as.matrix(dd1km_pred$soil_current[p$LinkID,])

old <- as.matrix(dd1km_old$soil_reference[p$LinkID,])
new <- as.matrix(dd1km_pred$soil_reference[p$LinkID,])

perc1 <- cbind(SOIL=1 - rowSums(old[,c("UNK","Water")])/rowSums(old),
    UNK=old[,c("UNK")]/rowSums(old),
    Water=old[,c("Water")]/rowSums(old))
perc2 <- cbind(SOIL=1 - rowSums(new[,c("UNK","Water")])/rowSums(new),
    UNK=new[,c("UNK")]/rowSums(new),
    Water=new[,c("Water")]/rowSums(new))
    
perc <- data.frame(old=perc1, new=perc2)
write.csv(perc, file="c:/p/Parkland_NoSoilInfo_OldVsNew.csv")

summary(1 - rowSums(new[,c("UNK","Water")])/rowSums(new))


kgrid$soil <- kgrid$pSoil >= 0.95

all(rownames(kgrid) == rownames(dd1km_pred[[1]]))
all(rownames(kgrid) == rownames(dd1km_pred[[2]]))
all(rownames(kgrid) == rownames(dd1km_pred[[3]]))
all(rownames(kgrid) == rownames(dd1km_pred[[4]]))

w1 <- as.numeric(dd1km_pred$veg_current[,c("Water")]) / (kgrid$Area_km2 * 10^6)
w2 <- as.numeric(dd1km_pred$veg_reference[,c("Water")]) / (kgrid$Area_km2 * 10^6)
w3 <- as.numeric(dd1km_pred$soil_current[,c("Water")]) / (kgrid$Area_km2 * 10^6)
w4 <- as.numeric(dd1km_pred$soil_reference[,c("Water")]) / (kgrid$Area_km2 * 10^6)

w1 <- as.numeric(dd1km_pred$veg_current[,c("Water")]) / rowSums(dd1km_pred$veg_current)
w2 <- as.numeric(dd1km_pred$veg_reference[,c("Water")]) / rowSums(dd1km_pred$veg_current)
w3 <- as.numeric(dd1km_pred$soil_current[,c("Water")]) / rowSums(dd1km_pred$veg_current)
w4 <- as.numeric(dd1km_pred$soil_reference[,c("Water")]) / rowSums(dd1km_pred$veg_current)

range(w1)
range(w2)
range(w3[kgrid$soil])
range(w4[kgrid$soil])

summary(w1 - w2)

summary(w1[kgrid$soil] - w3[kgrid$soil])
summary(w1[kgrid$soil] - w4[kgrid$soil])

summary(w2[kgrid$soil] - w3[kgrid$soil])
summary(w2[kgrid$soil] - w4[kgrid$soil])

summary(w3[kgrid$soil] - w4[kgrid$soil])


summary(abs(w1 - w2))

summary(abs(w1[kgrid$soil] - w3[kgrid$soil]))
summary(abs(w1[kgrid$soil] - w4[kgrid$soil]))

summary(abs(w2[kgrid$soil] - w3[kgrid$soil]))
summary(abs(w2[kgrid$soil] - w4[kgrid$soil]))

summary(abs(w3[kgrid$soil] - w4[kgrid$soil]))


kgrid2 <- kgrid
kgrid2$wdiff <- w1-w3
kgrid2 <- kgrid2[order(abs(kgrid2$wdiff), decreasing=TRUE),]
kgrid2 <- kgrid2[kgrid2$soil,]
kgrid2 <- kgrid2[1:100,]




by(w1-w3, list(NR=kgrid$NRNAME), summary)
by(abs(w1-w3), list(NR=kgrid$NRNAME), summary)
by(abs(w1-w3), list(NSR=kgrid$NSRNAME), summary)

kgrid2 <- kgrid
kgrid2$wdiff <- w1-w3
kgrid2 <- kgrid2[order(abs(kgrid2$wdiff), decreasing=TRUE),]

plot(w1[kgrid$soil], w3[kgrid$soil])

## ------------------------------------------

e1 <- new.env()
load(file.path("c:/p", "AB_data_v2015", "out/abmi_onoff", 
    "veg-hf-clim-reg_abmi-onoff_fix-fire.Rdata"), envir=e1)
e2 <- new.env()
load(file.path(ROOT, VER, "out/abmi_onoff", 
    "veg-hf-clim-reg_abmi-onoff_fix-fire.Rdata"), envir=e2)

s1 <- e1$dd1ha[[3]]
s2 <- e2$dd1ha[[3]]
#s1 <- e1$dd1ha[[4]]
#s2 <- e2$dd1ha[[4]]
#s1 <- e1$dd150m[[3]]
#s2 <- e2$dd150m[[3]]
all(rownames(s1) == rownames(s2))

#slt <- read.csv("~/repos/abmianalytics/lookup/lookup-soil.csv")
#rownames(slt) <- slt$SOILclass
slt <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

s1x <- groupSums(s1, 2, slt[colnames(s1), "UseInAnalysis"])
s2x <- groupSums(s2, 2, slt[colnames(s2), "UseInAnalysis"])

ss <- s2x[s2x[,"SoilUnknown"] == 0, -1]
ss <- ss / rowSums(ss)
summary(ss)
summary(rowSums(ss))
ss <- ss[,!(colnames(ss) %in% c("HWater","SoilWater"))] # Soil Water excl

op <- par(mfrow=c(4,3))
for (i in 1:ncol(ss)) {
    hist(100*ss[,i], main=colnames(ss)[i], col="gold")
}
par(op)

sss <- ss[,!(colnames(ss) %in% c("SoftLin","HardLin"))] # Soil Water excl
aa <- data.frame(table(SoilType=colnames(sss)[apply(sss, 1, which.max)]))
aa$MeanPercent <- round(100*colMeans(sss[,as.character(aa$SoilType)]), 2)
aa


do 150 m, include HF
