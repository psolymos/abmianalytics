library(RSQLite)

od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)

f <- "s:/GC_eric/FromEric/20210525_SC_veghf3by7_2019.sqlite"

db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)

area <- dbReadTable(db, "SC_SevenGeneration_Area3by7_rawdata")
d0 <- dbReadTable(db, "SC_SevenGeneration_Vegetation_VHF3by7_rawdata")
dbDisconnect(db)



d0$Shape_Area <- d0$Area_m2
d0$SITE_YR <- paste0(d0$ABMI_ID, "_", d0$year_3by7)

d1 <- d0[d0$Region=="Region",]
d2 <- d0[d0$Region=="7G Lease Area",]
d3 <- d0[d0$Region=="Watershed",]

d1w <- make_vegHF_wide_v6(d1,
    col.label="SITE_YR",
    col.year="year_3by7",
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d1, d1$SITE_YR, TRUE)[rownames(d1w[[1]]),]
d1w <- fill_in_0ages_v6(d1w, dx$NSRNAME, ages_list)


d2w <- make_vegHF_wide_v6(d2,
    col.label="SITE_YR",
    col.year="year_3by7",
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d2, d2$SITE_YR, TRUE)[rownames(d2w[[1]]),]
d2w <- fill_in_0ages_v6(d2w, dx$NSRNAME, ages_list)

d3w <- make_vegHF_wide_v6(d3,
    col.label="SITE_YR",
    col.year="year_3by7",
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d3, d3$SITE_YR, TRUE)[rownames(d3w[[1]]),]
d3w <- fill_in_0ages_v6(d3w, dx$NSRNAME, ages_list)


veg <- list(
  region=as.matrix(d1w$veg_current),
  lease_area_7g=as.matrix(d2w$veg_current),
  watershed=as.matrix(d3w$veg_current))

save(area, veg, file="s:/AB_data_v2021/data/SevenGeneration_Area3by7_2021-05-26.RData")

## -norbord

db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)

area <- dbReadTable(db, "SC_Norbord_Area3by7_rawdata")
d0 <- dbReadTable(db, "SC_Norbord_Vegetation_VHF3by7_rawdata")
dbDisconnect(db)



d0$Shape_Area <- d0$Area_m2
d0$SITE_YR <- paste0(d0$ABMI_ID, "_", d0$year_3by7)

table(d0$Region)

d1 <- d0[d0$Region=="Northern Area of Analysis",]
d2 <- d0[d0$Region=="Southern Area of Analysis",]
d3 <- d0[d0$Region=="Reference Area",]


d1w <- make_vegHF_wide_v6(d1,
    col.label="SITE_YR",
    col.year="year_3by7",
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d1, d1$SITE_YR, TRUE)[rownames(d1w[[1]]),]
d1w <- fill_in_0ages_v6(d1w, dx$NSRNAME, ages_list)


d2w <- make_vegHF_wide_v6(d2,
    col.label="SITE_YR",
    col.year="year_3by7",
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d2, d2$SITE_YR, TRUE)[rownames(d2w[[1]]),]
d2w <- fill_in_0ages_v6(d2w, dx$NSRNAME, ages_list)

d3w <- make_vegHF_wide_v6(d3,
    col.label="SITE_YR",
    col.year="year_3by7",
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d3, d3$SITE_YR, TRUE)[rownames(d3w[[1]]),]
d3w <- fill_in_0ages_v6(d3w, dx$NSRNAME, ages_list)


veg <- list(
  northern_area=as.matrix(d1w$veg_current),
  sourthern_area=as.matrix(d2w$veg_current),
  reference_area=as.matrix(d3w$veg_current))

save(area, veg, file="s:/AB_data_v2021/data/Norbord_Area3by7_2021-05-26.RData")
