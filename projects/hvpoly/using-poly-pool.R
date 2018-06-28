options(stringsAsFactors = FALSE)
HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")

library(dplyr)
lut <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-v6_forMetaData.csv")
lut <- lut %>%
  rename(PreBackfill_Source=preBackfill_Source)
DomInEachNReg <- read.csv("~/repos/abmianalytics/lookup/DomInEachNReg.csv")
upland <- c('AlpineLarch', 'Decid', 'Fir', 'Mixedwood', 'Pine', 'Spruce')
Harvest_Area <-  c("CUTBLOCK","HARVEST-AREA")

## south
xs <- read.csv(file.path(ROOT, "AB_data_v2018", "data", "raw", "hvpoly", "Backfilled100kmtestarea",
    "south-100km-lat-long", "south-converted-attribute-table.csv"))

## restore truncated column headers
xs$Origin_Year <- xs$Origin_Yea
xs$Origin_Yea <- NULL
xs$PreBackfill_Source <- xs$PreBackfil
xs$PreBackfil <- NULL
xs$Moisture_Reg <- xs$Moisture_R
xs$Moisture_R <- NULL
xs$Pct_of_Larch <- xs$Pct_of_Lar
xs$Pct_of_Lar <- NULL
xs$Combined_ChgByCWCS <- Combine_ChgByCWCS(xs)

dfs <- make_vegHF_wide_v6(xs,
    col.label="OBJECTID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes
dfs$VEGAGEclass <- make_older(dfs$VEGAGEclass, "5")
dfs$VEGHFAGEclass <- make_older(dfs$VEGHFAGEclass, "5")

## north
xn <- read.csv(file.path(ROOT, "AB_data_v2018", "data", "raw", "hvpoly", "Backfilled100kmtestarea",
    "north-100km-lat-long", "north-converted-attribute-table.csv"))

## restore truncated column headers
xn$Origin_Year <- xn$Origin_Yea
xn$Origin_Yea <- NULL
xn$PreBackfill_Source <- xn$PreBackfil
xn$PreBackfil <- NULL
xn$Moisture_Reg <- xn$Moisture_R
xn$Moisture_R <- NULL
xn$Pct_of_Larch <- xn$Pct_of_Lar
xn$Pct_of_Lar <- NULL
xn$Combined_ChgByCWCS <- Combine_ChgByCWCS(xn)

dfn <- make_vegHF_wide_v6(xn,
    col.label="OBJECTID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes
dfn$VEGAGEclass <- make_older(dfn$VEGAGEclass, "5")
dfn$VEGHFAGEclass <- make_older(dfn$VEGHFAGEclass, "5")

cn <- c("OBJECTID", "NSRNAME", "NRNAME", "LUF_NAME",
    "Shape_Area", "xcoord", "ycoord",
    "PreBackfill_Source", "Moisture_Reg",
    "Pct_of_Larch", "Combined_ChgByCWCS",
    "HF_Year", "SampleYear", "Origin_Year",
    "HFclass", "VEGclass", "AgeRf",
    "CC_ORIGIN_YEAR", "AgeCr", "VEGAGEclass", "VEGHFclass", "VEGHFAGEclass",
    "SOILclass", "SOILHFclass")

## write to an SQLite db
library(DBI)
f <- file.path(ROOT, "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")
con <- dbConnect(RSQLite::SQLite(), f)
dbWriteTable(con, "south", dfs[,cn], overwrite = TRUE)
dbWriteTable(con, "north", dfn[,cn], overwrite = TRUE)
dbListTables(con)
dbDisconnect(con)


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


