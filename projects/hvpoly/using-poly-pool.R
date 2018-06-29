## using the polygon tool
library(DBI)
library(cure4insect)
library(rgdal)
library(mefa4)
load_common_data()

f <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")

db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)
#x <- dbReadTable(db, "south")
x <- dbReadTable(db, "north")
dbDisconnect(db)

for (i in 1:ncol(x))
    if (is.character(x[,i]))
        x[,i] <- as.factor(x[,i])

veg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
veg_mapping <- cbind(V6=as.character(veg$VEGHFAGE_FINE), V5=as.character(veg$PolyReclass))
x$VEGAGEclass2 <- reclass(x$VEGAGEclass, veg_mapping)
setdiff(x$VEGAGEclass2, get_levels()$veg)
x$VEGHFAGEclass2 <- reclass(x$VEGHFAGEclass, veg_mapping)
setdiff(x$VEGHFAGEclass2, get_levels()$veg)

soil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2014.csv")
soil_mapping <- cbind(In=as.character(soil$SOILHF_FINE), Out=as.character(soil$UseInAnalysis))
x$SOILclass2 <- reclass(x$SOILclass, soil_mapping)
setdiff(x$SOILclass2, get_levels()$soil)
x$SOILHFclass2 <- reclass(x$SOILHFclass, soil_mapping)
setdiff(x$SOILHFclass, soil_mapping[,1])

xy <- x[,c("xcoord", "ycoord")]
coordinates(xy) <- ~ xcoord + ycoord
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

species <- "Achillea.millefolium"
object <- load_spclim_data(species)

## reference
pred_ref <- predict(object, xy=xy, veg=x$VEGAGEclass2, soil=x$SOILclass2)
summary(pred_ref)
#table(x$VEGAGEclass2, is.na(pred_ref$veg))
#pred_ref[is.na(pred_ref)] <- 0 # water and non-veg

## current
pred_curr <- predict(object, xy=xy, veg=x$VEGHFAGEclass2, soil=x$SOILHFclass2)
summary(pred_curr)
#table(x$VEGHFAGEclass2, is.na(pred_curr$veg))
#pred_curr[is.na(pred_curr)] <- 0 # water and non-veg


