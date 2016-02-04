##% Processing backfilled veg + HF (cutblock ages incorporated)
##% P Solymos
##% April 28, 2015

SAVE <- TRUE

## root directory
ROOT <- "c:/p"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2015"
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
load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
## kgrid
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
## transitions:
load(file.path(ROOT, VER, "out", "transitions", paste0(i, ".Rdata")))



kgrid$psoil <- 1 - dd1km_pred$soil_reference[,"UNK"] / 10^6
summary(kgrid$psoil)
kgrid$psoil[kgrid$psoil < 0] <- 0
kgrid$soil <- kgrid$psoil >= 0.95

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

