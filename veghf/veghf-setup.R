##% Processing backfilled veg + HF (cutblock ages incorporated)
##% P Solymos
##% Feb 5, 2016

#### Vegetation and HF processing (backfilled V5 and verified HF for most points)
#### (2012 w2w HF for some BAM/BBS points and prediction grid)

SAVE <- TRUE

## root directory
ROOT <- "e:/peter"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2018"
## current year
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900
#HF_YEAR <- 2014 # HF inventory update year

library(mefa4)
library(DBI)
source("~/repos/abmianalytics/R/veghf_functions.R")
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## this has pasture/cultivation etc differences
if (HF_VERSION %in% c("2016_fine")) {
    hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type-v2014.csv")
    hftypes <- droplevels(hftypes[!is.na(hftypes$HF_GROUP_COMB) &
        !duplicated(hftypes$FEATURE_TY),])
    hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class-v2014.csv")
    hflt <- hfgroups[match(hftypes$HF_GROUP_COMB, hfgroups$HF_GROUP_COMB),]
    HF_YEAR <- 2016 # HF inventory update year
}
## this has pasture/cultivation etc differences
if (HF_VERSION %in% c("2014_fine", "2014v2_fine")) {
    hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type-v2014.csv")
    hftypes <- droplevels(hftypes[hftypes$Source=="W2W_HF2014",])
    hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class-v2014.csv")
    hflt <- hfgroups[match(hftypes$HF_GROUP_COMB, hfgroups$HF_GROUP_COMB),]
    HF_YEAR <- 2014 # HF inventory update year
}
## this has NO pasture/cultivation etc differences
if (HF_VERSION %in% c("2014_coarse", "2014v2_coarse")) {
    hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type-v2014.csv")
    hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class-v2014.csv")
    hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
    HF_YEAR <- 2014 # HF inventory update year
}
## this has NO pasture/cultivation etc differences
if (HF_VERSION == "2012") {
    hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type.csv")
    hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
    hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
    HF_YEAR <- 2012 # HF inventory update year
}
## this has NO pasture/cultivation etc differences
if (HF_VERSION == "2010_coarse") {
    hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type-v2014.csv")
    hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class-v2014.csv")
    hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
    HF_YEAR <- 2010 # HF inventory update year
}
rownames(hflt) <- hftypes$FEATURE_TY

Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0",
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0",
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")

recl <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-v6-comb.csv")
recl <- recl[,1:2]

VEG_LEVS <- c(
    "Alkali",
    "AlpineLarch",
    "Bare", "Decid", "Fir", "GraminoidFen",
    "GrassHerb", "Marsh", "Mixedwood", "Pine", "Shrub", "ShrubbyBog",
    "ShrubbyFen", "ShrubbySwamp", "SnowIce", "Spruce", "TreedBog-BSpr",
    "TreedFen-BSpr", "TreedFen-Decid", "TreedFen-Larch", "TreedFen-Mixedwood",
    "TreedSwamp-Conif", "TreedSwamp-Decid", "TreedSwamp-Fir", "TreedSwamp-Forest",
    "TreedSwamp-Mixedwood", "TreedSwamp-Spruce", "TreedWetland-Mixedwood",
    "Water")

