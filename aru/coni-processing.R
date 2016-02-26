library(mefa4)
library(pbapply)

ROOT <- "e:/peter/AB_data_v2016"

int <- read.csv(file.path(ROOT, "data", "aru-coni", 
    "2015_CONI_DetectionsByRecording_1minInterval.csv"))
tms <- read.csv(file.path(ROOT, "data", "aru-coni",
    "2015_CONIPeent3.4_30_70_truepositives_details.csv"))
sit <- read.csv(file.path(ROOT, "data", "aru-coni",
    "CONImodel_ARU_sitesv2.csv"))
loc <- read.csv(file.path(ROOT, "data", "aru-coni",
    "CONImodel_ARU_locations-proper-ID.csv"))
compare_sets(sit$ID, loc$ID)
setdiff(sit$ID, loc$ID)
setdiff(loc$ID, sit$ID)

## dd150m, dd1km, points
load(file.path(ROOT, "out", "aru", 
        "veg-hf-clim-reg_aru-coni_fix-fire_fix-age0.Rdata"))

hist(int$Size)
abline(v=30000000)


compare_sets(points$POINT_ID, sit$ID)
compare_sets(points$POINT_ID, int$ID)

compare_sets(sit$ID, int$ID)
## these must be the problem sites???
setdiff(sit$ID, int$ID)

tail(sort(setdiff(points$POINT_ID, sit$ID)))
tail(sort(setdiff(sit$ID, points$POINT_ID)))

## need to process ABMI 2015 list and combine that with the rest of the CONI pts
## 150m and 1km
