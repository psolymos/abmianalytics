library(cure4insect)

ROOT <- "e:/peter/AB_data_v2018/data/raw/hvpoly/Backfilled100kmtestarea"

xs <- read.csv(file.path(ROOT, "south-100km-lat-long", "south-converted-attribute-table.csv"))

HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")

xs$Combined_ChgByCWCS <- "Mixedwood"
xs$Origin_Year <- xs$Origin_Yea
xs$Origin_Yea <- NULL

dfs <- make_vegHF_wide_v6(xs,
    col.label="OBJECTID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes

xn <- read.csv(file.path(ROOT, "north-100km-lat-long", "north-converted-attribute-table.csv"))
