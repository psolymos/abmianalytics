##% Processing backfilled veg + HF (cutblock ages incorporated)
##% P Solymos
##% Feb 5, 2016

#### Vegetation and HF processing (backfilled V5 and verified HF for most points)
#### (2012 w2w HF for some BAM/BBS points and prediction grid)

SAVE <- TRUE

## root directory
ROOT <- "e:/peter"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2016"
## current year
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900
HF_YEAR <- 2012 # HF inventory update year

library(mefa4)
source("~/repos/abmianalytics/R/veghf_functions.R")
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type.csv")
hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY

Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0", 
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0", 
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")
