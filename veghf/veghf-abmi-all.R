#e:/peter/AB_data_v2018/data/raw/veghf/abmi

HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

## ABMI sites (on+off) cetre -----------------------------------------------

## point intersections
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_SiteCentre_point_rev03.csv")
d <- read.csv(f)
dd_point <- make_vegHF_wide_v6(d,
    col.label="Site_YEAR",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes

## 1 ha in 4 x 0.25ha quadrants
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_SiteCentre_1ha_4quad_rev05.csv")
d <- read.csv(f)
d$QID <- with(d, interaction(Site_YEAR, Section, sep="_", drop=TRUE))

dd <- make_vegHF_wide_v6(d,
    col.label="QID",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "1/4 ha quadrant at site centre [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, QID, TRUE)[rownames(dd[[1]]),]
dd_qha <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

#dw_qha <- Xtab(Shape_Area ~ QID + Soil_Type_1, d[d$Combined_ChgByCWCS == "Water",])
#dw_qha <- Xtab(Shape_Area ~ QID + Soil_Type_1, d)
tmp <- make_vegHF_wide_v6(d,
    col.label="QID",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes
tmp$Soil_Type_2 <- as.character(tmp$SOILHFclass)
tmp$Soil_Type_2[tmp$SOILHFclass == "Water"] <-
    as.character(tmp$Soil_Type_1[tmp$SOILHFclass == "Water"])
tmp$Soil_Type_2 <- factor(tmp$Soil_Type_2, unique(c(levels(tmp$Soil_Type_1),
    levels(tmp$SOILHFclass))))
dw_qha <- Xtab(Shape_Area ~ QID + Soil_Type_2, tmp)


## 1 ha
dd <- make_vegHF_wide_v6(d,
    col.label="Site_YEAR",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "1 ha square around site centre [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, Site_YEAR, TRUE)[rownames(dd[[1]]),]
dd_1ha <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

#dw_1ha <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_1, d[d$Combined_ChgByCWCS == "Water",])
#dw_1ha <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_1, d)
dw_1ha <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_2, tmp)

## 150m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_SiteCentre_buf150m_rev04.csv")
d <- read.csv(f)
dd <- make_vegHF_wide_v6(d,
    col.label="Site_YEAR",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "150m circle buffer around site centre [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, Site_YEAR, TRUE)[rownames(dd[[1]]),]
dd_150m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

#dw_150m <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_1, d[d$Combined_ChgByCWCS == "Water",])
#dw_150m <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_1, d)

## 564m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_SiteCentre_buf564m_rev04.csv")
d <- read.csv(f)
dd <- make_vegHF_wide_v6(d,
    col.label="Site_YEAR",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "564m circle buffer around site centre [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, Site_YEAR, TRUE)[rownames(dd[[1]]),]
dd_564m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

#dw_564m <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_1, d[d$Combined_ChgByCWCS == "Water",])
#dw_564m <- Xtab(Shape_Area ~ Site_YEAR + Soil_Type_1, d)

sapply(dd_1ha[1:4], function(z) summary(rowSums(z)))
sapply(dd_qha[1:4], function(z) summary(rowSums(z)))
sapply(dd_150m[1:4], function(z) summary(rowSums(z)))
sapply(dd_564m[1:4], function(z) summary(rowSums(z)))

rs <- 100 * rowSums(dd_qha[[1]]) / 2500
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

rs <- 100 * rowSums(dd_1ha[[1]]) / 10000
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

rs <- 100 * rowSums(dd_150m[[1]]) / (150^2*pi)
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

rs <- 100 * rowSums(dd_564m[[1]]) / (564^2*pi)
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

f1 <- file.path(ROOT, VER, "data", "raw", "clim",
    "site-center-2003-2016.csv")
f2 <- file.path(ROOT, VER, "data", "raw", "clim",
    "1_SiteCentre2017_Summary_Climate_data.csv")
xx1 <- read.csv(f1)
xx2 <- read.csv(f2)
rownames(xx2) <- xx2$Site_YEAR
rownames(xx1) <- xx1$Site_YEAR
xx2$PET <- xx2$Eref
xx2$pAspen <- xx2$Populus_tremuloides_brtpred_nofp

tmp <- strsplit(as.character(xx2$ABMI_ASSIGNED__SITE_ID), "-")
xx2$NearestOnGridSite <- sapply(tmp, function(z) {
    zz <- if (length(z) > 1) z[3] else z[1]
    as.integer(gsub("\\D+", "", zz))
})
xx2$ABMI_Assigned_Site_ID <- xx2$ABMI_ASSIGNED__SITE_ID
cn <- intersect(colnames(xx1), colnames(xx2))
xx <- rbind(xx1[,cn], xx2[,cn])
xx <- data.frame(xx, meta[match(xx$NearestOnGridSite, meta$SITE_ID),])
rownames(xx) <- xx$Site_YEAR

compare_sets(dd_point$Site_YEAR, rownames(xx))
compare_sets(rownames(dd_1ha[[1]]), rownames(xx))
compare_sets(rownames(dd_150m[[1]]), rownames(xx))
compare_sets(rownames(dd_564m[[1]]), rownames(xx))

save(dd_point, dd_qha, dd_1ha, dd_150m, dd_564m, xx,
    dw_qha, dw_1ha, #dw_150m, dw_564m,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_SiteCenter_v6verified.Rdata"))

## ABMI cam/aru -----------------------------------------------

pat <- c(paste0("-", c("NW", "NE", "SW", "SE"), "-BOTH"),
    paste0("-", c("NW", "NE", "SW", "SE"), "-CAM"),
    paste0("-", c("NW", "NE", "SW", "SE"), "-ARU"))

## point intersections
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_CamARUBird_point_rev04.csv")
d <- read.csv(f)

d$SITE <- d$Site_ID
a <- levels(d$SITE)
a[endsWith(a, "B")] <- substr(a[endsWith(a, "B")], 1, nchar(a[endsWith(a, "B")])-1)
a[endsWith(a, "-")] <- substr(a[endsWith(a, "-")], 1, nchar(a[endsWith(a, "-")])-1)
a <- Gsub("-b", "", a)
a <- Gsub(pat, "", a)
levels(d$SITE) <- a

d$Site_bird_year <- as.factor(paste0(
    as.character(d$SITE), "_",
    ifelse(is.na(d$Cam_ARU_Bird_Location), "NA", as.character(d$Cam_ARU_Bird_Location)),
    "_", as.character(d$survey_year)))
#d[is.na(d$Site_YEAR_bird),]
dd_point <- make_vegHF_wide_v6(d,
    col.label="Site_bird_year",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes


## 150m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_CamARUBird_buf150m_rev05.csv")
d <- read.csv(f)

d$SITE <- d$Site_ID
a <- levels(d$SITE)
a[endsWith(a, "B")] <- substr(a[endsWith(a, "B")], 1, nchar(a[endsWith(a, "B")])-1)
a[endsWith(a, "-")] <- substr(a[endsWith(a, "-")], 1, nchar(a[endsWith(a, "-")])-1)
a <- Gsub("-b", "", a)
a <- Gsub(pat, "", a)
levels(d$SITE) <- a

d$Site_bird_year <- as.factor(paste0(
    as.character(d$SITE), "_",
    ifelse(is.na(d$Cam_ARU_Bird_Location), "NA", as.character(d$Cam_ARU_Bird_Location)),
    "_", as.character(d$survey_year)))
dd <- make_vegHF_wide_v6(d,
    col.label="Site_bird_year",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "150m circle buffer around Camera/ARU [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, Site_bird_year, TRUE)[rownames(dd[[1]]),]
dd_150m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

## 564m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_CamARUBird_buf564m_rev05.csv")
d <- read.csv(f)

d$SITE <- d$Site_ID
a <- levels(d$SITE)
a[endsWith(a, "B")] <- substr(a[endsWith(a, "B")], 1, nchar(a[endsWith(a, "B")])-1)
a[endsWith(a, "-")] <- substr(a[endsWith(a, "-")], 1, nchar(a[endsWith(a, "-")])-1)
a <- Gsub("-b", "", a)
a <- Gsub(pat, "", a)
levels(d$SITE) <- a

d$Site_bird_year <- as.factor(paste0(
    as.character(d$SITE), "_",
    ifelse(is.na(d$Cam_ARU_Bird_Location), "NA", as.character(d$Cam_ARU_Bird_Location)),
    "_", as.character(d$survey_year)))
dd <- make_vegHF_wide_v6(d,
    col.label="Site_bird_year",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "564m circle buffer around Camera/ARU [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, Site_bird_year, TRUE)[rownames(dd[[1]]),]
dd_564m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

sapply(dd_150m[1:4], function(z) summary(rowSums(z)))
sapply(dd_564m[1:4], function(z) summary(rowSums(z)))

rs <- 100 * rowSums(dd_150m[[1]]) / (150^2*pi)
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

rs <- 100 * rowSums(dd_564m[[1]]) / (564^2*pi)
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

save(dd_point, dd_150m, dd_564m,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_CameraARU_v6verified.Rdata"))

## BAM/BBS/BU -----------------------------------------------

## point intersections
f <- file.path(ROOT, VER, "data", "raw", "veghf", "bambbs",
    "20180627_BAMBBS_1993_2017.sqlite")

db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)

d00 <- dbReadTable(db, "BU_SiteCentre")
#d10 <- dbReadTable(db, "BU_Devon_SiteCentre")
d01 <- dbReadTable(db, "BU_Buffer150m")
#d11 <- dbReadTable(db, "BU_Devon_Buffer150m")
d02 <- dbReadTable(db, "BU_Buffer564m")
#d12 <- dbReadTable(db, "BU_Devon_Buffer564m")

dbDisconnect(db)

d00 <- make_char2fact(d00)

dd_point <- make_vegHF_wide_v6(d00,
    col.label="UID",
    col.year="Survey_Year",
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes

## 150m
keep <- c("Origin_Year",
    "Pct_of_Larch",
    "NSRNAME",
    "Soil_Type_1",
    "UID",
    "Survey_Year",
    "FEATURE_TY",
    "Combined_ChgByCWCS",
    "YEAR",
    "Shape_Area")
d01 <- make_char2fact(d01[,keep])
summary(d01[d01$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA"), "YEAR"])
d01[d01$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & is.na(d01$YEAR), "YEAR"] <- 1930
summary(d01[d01$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA"), "YEAR"])

dd <- make_vegHF_wide_v6(d01,
    col.label="UID",
    col.year="Survey_Year",
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "150m circle buffer around bird point [V6 backfilled + verified HF]"
dx <- nonDuplicated(d01, UID, TRUE)[rownames(dd[[1]]),]
dd_150m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

## 564m
d02 <- make_char2fact(d02[,keep])
summary(d02[d02$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA"), "YEAR"])
d02[d02$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & is.na(d02$YEAR), "YEAR"] <- 1930
summary(d02[d02$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA"), "YEAR"])
dd <- make_vegHF_wide_v6(d02,
    col.label="UID",
    col.year="Survey_Year",
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "564m circle buffer around site centre [V6 backfilled + verified HF]"
dx <- nonDuplicated(d02, UID, TRUE)[rownames(dd[[1]]),]
dd_564m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

sapply(dd_150m[1:4], function(z) summary(rowSums(z)))
sapply(dd_564m[1:4], function(z) summary(rowSums(z)))

rs <- 100 * rowSums(dd_150m[[1]]) / (150^2*pi)
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

rs <- 100 * rowSums(dd_564m[[1]]) / (564^2*pi)
summary(rs)
table(cut(rs, c(0, 99, 101, Inf)))
data.frame(percent=rs[which(rs < 99 | rs > 101)])

save(dd_point, dd_150m, dd_564m,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_BAM-BBS-BU_v6verified.Rdata"))


## All sites -----------------------------------------------

f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "20180706_All_Sites.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)
d0 <- dbReadTable(db, "All_Sites_SiteCentre")
d1 <- dbReadTable(db, "All_Sites_1ha")
d2 <- dbReadTable(db, "All_Sites_buffer564m")
dbDisconnect(db)

d0$ABMI <- as.integer(d0$ABMI)
d1$ABMI <- as.integer(d1$ABMI)
d2$ABMI <- as.integer(d2$ABMI)

dd_point <- make_vegHF_wide_v6(d0,
    col.label="ABMI",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes

dd <- make_vegHF_wide_v6(d1,
    col.label="ABMI",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "1 ha square around site centre [V6 backfilled + 2016 HF]"
dx <- nonDuplicated(d1, ABMI, TRUE)[rownames(dd[[1]]),]
dd_1ha <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

dd <- make_vegHF_wide_v6(d2,
    col.label="ABMI",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "564m circle buffer around site centre [V6 backfilled + 2016 HF]"
dx <- nonDuplicated(d2, ABMI, TRUE)[rownames(dd[[1]]),]
dd_564m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

save(dd_point, dd_1ha, dd_564m,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_allSites_v6hfi2016.Rdata"))

x1 <- as.matrix(dd_1ha[[1]])
x1 <- round(100*x1 / rowSums(x1), 4)
write.csv(x1, file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_allSites_v6hfi2016-1ha.csv"))
x2 <- as.matrix(dd_564m[[1]])
x2 <- round(100*x2 / rowSums(x2), 4)
write.csv(x2, file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_allSites_v6hfi2016-564m.csv"))

## Wetlands ----------------------------

f <- file.path(ROOT, VER, "data", "raw", "veghf", "wetlands",
    "20180820_Wetlands_2018.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
tabs <- list()
for (i in dbListTables(db))
    tabs[[i]] <- dbReadTable(db, i)
dbDisconnect(db)
str(tabs)

nm_ci <- c("AGGC_CI_2009", "AGGC_CI_2010", "AGGC_CI_2011", "AGGC_CI_2012",
    "AGGC_CI_2013", "AGGC_CI_2014", "AGGC_CI_2015", "AGGC_CI_2016",
    "AGGC_CI_2017")
nm_clim <- c("Climate_AHM", "Climate_Eref", "Climate_FFP",
    "Climate_MAP", "Climate_MAT", "Climate_MCMT", "Climate_MWMT",
    "Climate_Populus_tremuloides_brtpred_nofp")
nm_buff <- "SummaryTableInBuffers"
nm_catch <- "SummaryTableInCatchments_rev01"

## climate stuff
t(sapply(nm_clim, function(i)
    sapply(1:3, function(j) all(tabs[[nm_clim[1]]][,j] == tabs[[i]][,j]))))
clim <- data.frame(tabs[[nm_clim[1]]][,1:3],
    sapply(tabs[nm_clim], "[[", "RASTERVALU"))
clim$SiteYear <- paste0(clim$Site_ID_Ref, "_", clim$year)
for (i in 1:ncol(clim))
    if (is.character(clim[,i]))
        clim[,i] <- as.factor(clim[,i])

## catchment level veg/hf/soil, veg61 and HFI2014
d <- tabs[[nm_catch]]
sum(is.na(d$Combined_ChgByCWCS))
#d <- d[!is.na(d$Combined_ChgByCWCS),]
d$SiteYear <- paste0(d$Site_ID, "_", d$Year)
for (i in 1:ncol(d))
    if (is.character(d[,i]))
        d[,i] <- as.factor(d[,i])

dd <- make_vegHF_wide_v6(d,
    col.label="SiteYear",
    col.year="Year",
    col.HFyear="YEAR_1",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "Catchment area around wetland [V6 backfilled + 2014v2 HFI]"
dx <- nonDuplicated(d, SiteYear, TRUE)[rownames(dd[[1]]),]
ddw_catch <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

## buffer level veg/hf/soil, 0-20, 20-100 and 100-250 m buffer
d <- tabs[[nm_buff]]
#d <- d[!is.na(d$Combined_ChgByCWCS),]
d$SiteYear <- paste0(d$Site_ID, "_", d$Year)
for (i in 1:ncol(d))
    if (is.character(d[,i]))
        d[,i] <- as.factor(d[,i])

ddw_point <- make_vegHF_wide_v6(d[d$buffer == "20m",],
    col.label="SiteYear",
    col.year="Year",
    col.HFyear="YEAR_1",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes

dd <- make_vegHF_wide_v6(d[d$buffer == "20m",],
    col.label="SiteYear",
    col.year="Year",
    col.HFyear="YEAR_1",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "0-20 m buffer area around wetland [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, SiteYear, TRUE)[rownames(dd[[1]]),]
ddw_20m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

dd <- make_vegHF_wide_v6(d[d$buffer %in% c("20m", "20_100m"),],
    col.label="SiteYear",
    col.year="Year",
    col.HFyear="YEAR_1",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "0-100 m buffer area around wetland [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, SiteYear, TRUE)[rownames(dd[[1]]),]
ddw_100m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

dd <- make_vegHF_wide_v6(d,
    col.label="SiteYear",
    col.year="Year",
    col.HFyear="YEAR_1",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "0-250 m buffer area around wetland [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, SiteYear, TRUE)[rownames(dd[[1]]),]
ddw_250m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

save(clim, ddw_catch, ddw_20m, ddw_100m, ddw_250m, ddw_point,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_wetlands_v6x.Rdata"))

## Mammal inter level transects --------------------------------

f <- file.path(ROOT, VER, "data", "raw", "veghf", "mammals",
    "20180822_MammalTransects_VEG61_2018.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
d <- dbReadTable(db, "MammalTransects_Interlevel")
dbDisconnect(db)

for (cn in c("year", "Origin_Year", "Pct_of_Larch", "year_1", "Shape_Area"))
    d[,cn] <- as.numeric(d[,cn])
for (i in 1:ncol(d))
    if (is.character(d[,i]))
        d[,i] <- as.factor(d[,i])
d$UID <- as.factor(paste0(d$ABMISite, "_", d$interLevel, "_", d$year))

dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="year",
    col.HFyear="year_1",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "0-250(?) m buffer area around inter level snow track segments [V6 backfilled + verified HF]"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd_inter <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

save(dd_inter,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_mammals_v6verifiedHF.Rdata"))

## Protected areas -------------------------------------

f <- file.path(ROOT, VER, "data", "raw", "veghf", "misc",
    "20181030_PolyTool_Park_Protected.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)
d <- dbReadTable(db, "ParkProtected_Veg61HF2016FTY_InternalUseOnly")
dbDisconnect(db)

d <- make_char2fact(d)
d$NAME_NRNAME <- interaction(d$NAME, d$NRNAME, sep="_", drop=TRUE)

dd <- make_vegHF_wide_v6(d,
    col.label="NAME_NRNAME",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "ParkProtected_Veg61HF2016FTY_InternalUseOnly"
dx <- nonDuplicated(d, NAME_NRNAME, TRUE)[rownames(dd[[1]]),]
dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

f <- file.path(ROOT, VER, "data", "raw", "veghf", "misc",
    "20181113_PolyTool_Crown_LarpOutsideProtectedAndCrown.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)
d1 <- dbReadTable(db, "CrownReservations_Veg61HF2016FTY_InternalUseOnly")
d2 <- dbReadTable(db, "LarpOutsideCrownAndProtected_Veg61HF2016FTY_InternalUseOnly")
dbDisconnect(db)

d1 <- make_char2fact(d1)
d2 <- make_char2fact(d2)
d1$NAME_NRNAME <- interaction(d1$NAME, d1$NRNAME, sep="_", drop=TRUE)
d2$NAME_NRNAME <- interaction(d2$AreaName, d2$NRNAME, sep="_", drop=TRUE)

dd1 <- make_vegHF_wide_v6(d1,
    col.label="NAME_NRNAME",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd1$scale <- "CrownReservations_Veg61HF2016FTY_InternalUseOnly"
dx1 <- nonDuplicated(d1, NAME_NRNAME, TRUE)[rownames(dd1[[1]]),]
dd1 <- fill_in_0ages_v6(dd1, dx1$NSRNAME, ages_list)
dd2 <- make_vegHF_wide_v6(d2,
    col.label="NAME_NRNAME",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd2$scale <- "LarpOutsideCrownAndProtected_Veg61HF2016FTY_InternalUseOnly"
dx2 <- nonDuplicated(d2, NAME_NRNAME, TRUE)[rownames(dd2[[1]]),]
dd2 <- fill_in_0ages_v6(dd2, dx2$NSRNAME, ages_list)

save(dd, dx, dd1, dx1, dd2, dx2,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_protected-areas.Rdata"))

## 3x7km yearly HF --------------------------------------

f <- file.path(ROOT, VER, "data", "raw", "veghf", "3x7",
    "20181105_Backfill61_3x7.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
lt <- dbListTables(db)
yr <- as.integer(substr(lt, nchar(lt)-3, nchar(lt)))

for (iyr in seq_len(length(yr))) {
    cat("\n\nYEAR", yr[iyr], "------------------------\n\n")
    flush.console()
    d <- dbReadTable(db, lt[iyr])

    for (i in 1:ncol(d))
        if (is.character(d[,i]))
            d[,i] <- as.factor(d[,i])

    summary(d[d$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA"), "YEAR"])
    sum(d[d$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & is.na(d$YEAR), "Shape_Area"])/10^6
    d[d$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & is.na(d$YEAR), "YEAR"] <- 1930
    summary(d[d$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA"), "YEAR"])
    sum(d[d$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & is.na(d$YEAR), "Shape_Area"])/10^6

    dd <- make_vegHF_wide_v6(d,
        col.label="ABMI_ID",
        col.year=yr[iyr],
        col.HFyear="YEAR",
        col.HABIT="Combined_ChgByCWCS",
        col.SOIL="Soil_Type_1",
        sparse=TRUE, HF_fine=TRUE) # use refined classes
    dd$scale <- "Backfill61 + VerifiedHF in 3x7km"
    dx <- nonDuplicated(d, ABMI_ID, TRUE)[rownames(dd[[1]]),]
    dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

    save(dd,
        file=file.path(ROOT, VER, "data", "inter", "veghf", "3x7",
        paste0("veg-hf_yearly-hf_", yr[iyr],"_bf61verifiedHF.Rdata")))

}
dbDisconnect(db)

for (iyr in seq_len(length(yr))) {
    cat("YEAR", yr[iyr], "\n")
    flush.console()
    load(file=file.path(ROOT, VER, "data", "inter", "veghf", "3x7",
        paste0("veg-hf_yearly-hf_", yr[iyr],"_bf61verifiedHF.Rdata")))
    if (iyr == 1) {
        veg_curr <- array(0, c(dim(dd[[1]]), length(yr)))
        dimnames(veg_curr) <- c(dimnames(dd[[1]]), list(yr))
        veg_ref <- array(0, c(dim(dd[[2]]), length(yr)))
        dimnames(veg_ref) <- c(dimnames(dd[[2]]), list(yr))
        soil_curr <- array(0, c(dim(dd[[3]]), length(yr)))
        dimnames(soil_curr) <- c(dimnames(dd[[3]]), list(yr))
        soil_ref <- array(0, c(dim(dd[[4]]), length(yr)))
        dimnames(soil_ref) <- c(dimnames(dd[[4]]), list(yr))
        SITES <- dimnames(veg_curr)[[1]]
    }
    veg_curr[,,iyr] <- as.matrix(dd[[1]])[SITES,]
    veg_ref[,,iyr] <- as.matrix(dd[[2]])[SITES,]
    soil_curr[,,iyr] <- as.matrix(dd[[3]])[SITES,]
    soil_ref[,,iyr] <- as.matrix(dd[[4]])[SITES,]
}
z=veg_curr[,"MineSite",]
summary(z)
matplot(yr, t(z)/(21*10^6),lty=1,type="l")

save(veg_curr, veg_ref, soil_curr, soil_ref,
    file=file.path(ROOT, VER, "data", "analysis", "yearly",
    paste0("veg-hf_yearly-hf_all-years_bf61verifiedHF.Rdata")))


## w2w 1km^2 scale ---------- short & memory hungry version ----------------------------

t0 <- proc.time()
load(file.path(ROOT, VER, "data", "analysis", "kgrid_table_km.Rdata"))
f <- file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_VEG61HFI2016v3fty.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
tn <- dbListTables(db) # "Summary_1sqkm_grid_FTY"

hd <- dbGetQuery(db, paste('SELECT * FROM', tn, 'LIMIT 5'))
d <- dbGetQuery(db,
    paste("SELECT
      Origin_Year,
      Pct_of_Larch,
      NSRNAME,
      Soil_Type_1,
      GRID_LABEL,
      FEATURE_TY,
      Combined_ChgByCWCS,
      YEAR,
      Shape_Area
    FROM
      ", tn, ";")
)
dbDisconnect(db)

d <- make_char2fact(d)

d2 <- make_vegHF_wide_v6(d,
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d2 <- d2[,c("GRID_LABEL", "NSRNAME",
    "VEGAGEclass", "VEGHFAGEclass",
    "SOILclass", "SOILHFclass",
    "Shape_Area")]

## wide format: amount
dd <- make_vegHF_wide_v6(d2,
    col.label="GRID_LABEL",
    HF_fine=TRUE, wide=TRUE, widen_only=TRUE) # use refined classes
dd$scale <- "1 km^2 areas [V6 backfilled + 2016v3 w2w HFI]"
for (i in 1:4)
    dd[[i]] <- dd[[i]][rownames(kgrid),]
all(rownames(kgrid) == rownames(dd[[1]]))
dx <- nonDuplicated(d2, GRID_LABEL, TRUE)[rownames(dd[[1]]),]
dd_kgrid <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
proc.time() - t0

save(dd,
    file=file.path(ROOT, VER, "data", "analysis", "grid",
    "veg-hf_grid_v6hf2016v3noDistVeg-unsorted.Rdata"))
save(dd_kgrid,
    file=file.path(ROOT, VER, "data", "analysis", "grid",
    "veg-hf_grid_v6hf2016v3noDistVeg.Rdata"))
save(d2,
    file=file.path(ROOT, VER, "data", "inter", "veghf", "grid",
    "veg-hf_grid_v6hf2016v3noDistVeg-long-format.Rdata"))

## wide format: transitions --------------------------------------

load(file.path(ROOT, VER, "data", "inter", "veghf", "grid",
    "veg-hf_grid_v6hf2016v3noDistVeg-long-format.Rdata"))
load(file.path(ROOT, VER, "data", "analysis", "kgrid_table_km.Rdata"))

d2$soilTr <- as.character(d2$SOILclass)
ss <- as.character(d2$SOILclass) != as.character(d2$SOILHFclass)
#table(ss)
d2$soilTr[ss] <- paste0(as.character(d2$SOILclass[ss]),
    "->", as.character(d2$SOILHFclass[ss]))
d2$soilTr <- as.factor(d2$soilTr)
trSoil <- Xtab(Shape_Area ~ GRID_LABEL + soilTr, d2)
trSoil <- trSoil[rownames(kgrid),]
range(rowSums(trSoil)/10^6)
dim(trSoil)

d2$vegTr <- as.character(d2$VEGAGEclass)
ss <- as.character(d2$VEGAGEclass) != as.character(d2$VEGHFAGEclass)
#table(ss)
d2$vegTr[ss] <- paste0(as.character(d2$VEGAGEclass[ss]),
    "->", as.character(d2$VEGHFAGEclass[ss]))
d2$vegTr <- as.factor(d2$vegTr)
trVeg <- Xtab(Shape_Area ~ GRID_LABEL + vegTr, d2)
trVeg <- trVeg[rownames(kgrid),]
range(rowSums(trVeg)/10^6)
dim(trVeg)

ch2soil <- nonDuplicated(d2, soilTr, TRUE)[,c("SOILclass", "SOILHFclass")]
colnames(ch2soil) <- c("rf","cr")
ch2veg <- nonDuplicated(d2, vegTr, TRUE)[,c("VEGAGEclass", "VEGHFAGEclass")]
colnames(ch2veg) <- c("rf","cr")

save(trVeg, trSoil, ch2veg, ch2soil,
    file=file.path(ROOT, VER, "data", "analysis", "grid",
    "veg-hf_transitions_v6hf2016v3noDistVeg.Rdata"))


## w2w 1km^2 scale ---------- long version ----------------------------

f <- file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
#    "20180821_Grid_1sqkm_VEG61HFI2016v3.sqlite")
    "20181010_Grid_1sqkm_VEG61HFI2016v3fty.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
tn <- dbListTables(db) # "Summary_1sqkm_grid_FTY"

hd <- dbGetQuery(db, paste('SELECT * FROM', tn, 'LIMIT 5'))
d <- dbGetQuery(db,
    paste("SELECT
      Origin_Year,
      Pct_of_Larch,
      NSRNAME,
      Soil_Type_1,
      GRID_LABEL,
      FEATURE_TY,
      Combined_ChgByCWCS,
      YEAR,
      Shape_Area
    FROM
      ", tn, "LIMIT 1000;")
)
d1 <- make_vegHF_wide_v6(d,
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes

d <- dbGetQuery(db, paste("SELECT GRID_LABEL FROM", tn, ";"))
d[,1] <- as.factor(d[,1])
gc()
d0 <- dbGetQuery(db, paste("SELECT NSRNAME FROM", tn, ";"))
d$NSRNAME <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Combined_ChgByCWCS FROM", tn, ";"))
d$Combined_ChgByCWCS <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Soil_Type_1 FROM", tn, ";"))
d$Soil_Type_1 <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT FEATURE_TY FROM", tn, ";"))
d$FEATURE_TY <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Origin_Year FROM", tn, ";"))
d$Origin_Year <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Pct_of_Larch FROM", tn, ";"))
d$Pct_of_Larch <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT YEAR FROM", tn, ";"))
d$YEAR <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Shape_Area FROM", tn, ";"))
d$Shape_Area <- d0[,1]
rm(d0)
gc()

dbDisconnect(db)

save(d, file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    #"20180821_Grid_1sqkm_VEG61HFI2016v3.RData"))
    "20181010_Grid_1sqkm_VEG61HFI2016v3fty.RData"))

## long format first, in chunks

load(file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    #"20180821_Grid_1sqkm_VEG61HFI2016v3.RData"))
    "20181010_Grid_1sqkm_VEG61HFI2016v3fty.RData"))

cn <- c("VEGAGEclass", "VEGHFAGEclass", "SOILclass", "SOILHFclass")
## 0-2
d1 <- make_vegHF_wide_v6(d[1:20000000,],
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d0 <- d1[,cn]
d <- d[-(1:20000000),]
rm(d1)
gc()
## 2-4
d1 <- make_vegHF_wide_v6(d[1:20000000,],
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d0 <- rbind(d0, d1[,cn])
d <- d[-(1:20000000),]
rm(d1)
gc()
## 4-6
d1 <- make_vegHF_wide_v6(d[1:20000000,],
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d0 <- rbind(d0, d1[,cn])
d <- d[-(1:20000000),]
rm(d1)
gc()
## 6-7
d1 <- make_vegHF_wide_v6(d[1:10000000,],
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d0 <- rbind(d0, d1[,cn])
d <- d[-(1:10000000),]
rm(d1)
gc()
## 7-8
d1 <- make_vegHF_wide_v6(d[1:10000000,],
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d0 <- rbind(d0, d1[,cn])
d <- d[-(1:10000000),]
rm(d1)
gc()
## >8
d1 <- make_vegHF_wide_v6(d[1:nrow(d),],
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes
d0 <- rbind(d0, d1[,cn])
rm(d1,d)
gc()

load(file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    #"20180821_Grid_1sqkm_VEG61HFI2016v3.RData"))
    "20181010_Grid_1sqkm_VEG61HFI2016v3fty.RData"))
d <- d[,c("GRID_LABEL", "NSRNAME", "Shape_Area")]

save(d, file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_common.RData"))
d0v <- d0[,1:2]
d0s <- d0[,3:4]
save(d0v, file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_veg.RData"))
save(d0s, file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_soil.RData"))

## wide format (taking preprocessed d0 and widen_only=TRUE)

load(file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_common.RData"))
load(file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_veg.RData"))
d$VEGAGEclass <- d0v$VEGAGEclass
d$VEGHFAGEclass <- d0v$VEGHFAGEclass
rm(d0v)
gc()

load(file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181010_Grid_1sqkm_soil.RData"))
d$SOILclass <- d0s$SOILclass
d$SOILHFclass <- d0s$SOILHFclass
rm(d0s)
gc()

dd <- make_vegHF_wide_v6(d,
    col.label="GRID_LABEL",
    widen_only=TRUE)
dd$scale <- "1 km^2 areas [V6 backfilled + 2016v3 w2w HFI]"
gc()
dx <- nonDuplicated(d, GRID_LABEL, TRUE)[rownames(dd[[1]]),]
dd_kgrid <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

load(file.path(ROOT, VER, "data", "analysis", "kgrid_table_km.Rdata"))
for (i in 1:4) {
    dd[[i]] <- dd[[i]][rownames(kgrid),]
    dd_kgrid[[i]] <- dd_kgrid[[i]][rownames(kgrid),]
}

save(dd,
    file=file.path(ROOT, VER, "data", "analysis", "grid",
    "veg-hf_grid_v6hf2016v3-unsorted.Rdata"))
save(dd_kgrid,
    file=file.path(ROOT, VER, "data", "analysis", "grid",
    "veg-hf_grid_v6hf2016v3.Rdata"))

## w2w poly tool --------------------------------------
## easting/northing is based on 10TM forest projection
## +proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs

f <- file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    "20181012_PolyTool.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
(tn <- dbListTables(db))

hd <- dbGetQuery(db, paste('SELECT * FROM', tn, 'LIMIT 5'))
d <- dbGetQuery(db,
    paste("SELECT
      Origin_Year,
      Pct_of_Larch,
      NSRNAME,
      Soil_Type_1,
      UID,
      FEATURE_TY,
      Combined_ChgByCWCS,
      YEAR,
      Shape_Area,
      Easting, Northing
    FROM
      ", tn, "LIMIT 1000;")
)
d1 <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE) # use refined classes

d <- dbGetQuery(db, paste("SELECT UID FROM", tn, ";"))
d[,1] <- d[,1]
gc()
d0 <- dbGetQuery(db, paste("SELECT NSRNAME FROM", tn, ";"))
d$NSRNAME <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Combined_ChgByCWCS FROM", tn, ";"))
d$Combined_ChgByCWCS <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Soil_Type_1 FROM", tn, ";"))
d$Soil_Type_1 <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT FEATURE_TY FROM", tn, ";"))
d$FEATURE_TY <- as.factor(d0[,1])
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Origin_Year FROM", tn, ";"))
d$Origin_Year <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Pct_of_Larch FROM", tn, ";"))
d$Pct_of_Larch <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT YEAR FROM", tn, ";"))
d$YEAR <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Shape_Area FROM", tn, ";"))
d$Shape_Area <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Easting FROM", tn, ";"))
d$Easting <- d0[,1]
rm(d0)
gc()
d0 <- dbGetQuery(db, paste("SELECT Northing FROM", tn, ";"))
d$Northing <- d0[,1]
rm(d0)
gc()

dbDisconnect(db)

save(d, file=file.path(ROOT, VER, "data", "raw", "veghf", "w2w",
    #"20180821_Grid_1sqkm_VEG61HFI2016v3.RData"))
    "20181010_Grid_1sqkm_VEG61HFI2016v3fty.RData"))



## add here transformations

## Cam/ARU/Bird -----------------------------------------------

## point intersections
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_CamARUBird_point_rev03.csv")
d <- read.csv(f)
d$UID <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
dd_point <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE) # use refined classes

## 150m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_CamARUBird_buf150m_rev04.csv")
d <- read.csv(f)
levels(d$Cam_ARU_Bird_Location) <- c(levels(d$Cam_ARU_Bird_Location), "UNKNOWN")
d$Cam_ARU_Bird_Location[is.na(d$Cam_ARU_Bird_Location)] <- "UNKNOWN"
d$UID <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="survey_year",
    col.HFyear="year",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
dd$scale <- "150m circle buffer around camera/ARU locations"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd_150m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
## 564m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all",
    "Summary_2003_2017_CamARUBird_buf564m_rev04.csv")
d <- read.csv(f)
levels(d$Cam_ARU_Bird_Location) <- c(levels(d$Cam_ARU_Bird_Location), "UNKNOWN")
d$Cam_ARU_Bird_Location[is.na(d$Cam_ARU_Bird_Location)] <- "UNKNOWN"
d$UID <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=TRUE) # use refined classes
dd$scale <- "564m circle buffer around camera/ARU locations"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd_564m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)



## comparing with previous years
if (FALSE) {

f1 <- file.path(ROOT, "AB_data_v2017", "data", "analysis", "site",
    "veg-hf_siteCenter_v6-fixage0.Rdata")
e1 <- new.env()
load(f1, envir=e1)

f2 <- file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_SiteCenter_v6verified.Rdata")
e2 <- new.env()
load(f2, envir=e2)

a1 <- as.matrix(e1$dd_1ha[[1]])
a1 <- a1/rowSums(a1)
a2 <- as.matrix(e2$dd_1ha[[1]])
a2 <- a2/rowSums(a2)

compare_sets(colnames(a1), colnames(a2))
setdiff(colnames(a1), colnames(a2))
setdiff(colnames(a2), colnames(a1))
cn <- colnames(a2)
names(cn) <- cn
cn[c("UrbanIndustrial","UrbanResidence")] <- "Urban"
cn[c("SeismicLineNarrow", "SeismicLineWide")] <- "SeismicLine"
cn[c("CultivationCrop", "CultivationAbandoned",
    "CultivationRoughPasture", "CultivationTamePasture")] <- "CultivationCropPastureBareground"
a2x <- groupSums(a2, 2, cn)
compare_sets(colnames(a1), colnames(a2x))
a1 <- a1[,colnames(a2x)]
compare_sets(rownames(a1), rownames(a2x))
rn <- intersect(rownames(a1), rownames(a2x))

cf <- t(sapply(colnames(a1), function(i) {
    m <- lm(new ~ old, data.frame(new=a2x[rn,i], old=a1[rn,i]))
    cf <- c(coef(m), sigma(m)^2)
    names(cf) <- c("intercept", "slope", "variance")
    cf
}))
df <- data.frame(class=colnames(a2x),
    round(cf, 3),
    mean_perc_old=round(100*colMeans(a1), 4),
    mean_perc_new=round(100*colMeans(a2x), 4))
write.csv(df, row.names=FALSE,
    file=file.path(ROOT, VER, "data", "inter", "veghf1ha-old-new-comparison.csv"))

}




xx <- read.csv(file.path(ROOT, VER, "data", "raw", "veghf",
    "site_all", "siteCenter_climate.csv"))
rownames(xx) <- xx$Site_YEAR

compare_sets(rownames(dd_1ha[[1]]), rownames(xx))
compare_sets(rownames(dd_150m[[1]]), rownames(xx))
compare_sets(rownames(dd_564m[[1]]), rownames(xx))
setdiff(rownames(xx), rownames(dd_1ha[[1]]))

xx <- droplevels(xx[rownames(dd_1ha[[1]]),])
xx$PET <- xx$Eref
xx$pAspen <- xx$Populus_tremuloides_brtpred_nofp
xx$Populus_tremuloides_brtpred_nofp <- NULL
xx$OBJECTID <- NULL
xx$NSRNAME <- dx$NSRNAME[match(rownames(xx), rownames(dx))]
xx$NRNAME <- dx$NRNAME[match(rownames(xx), rownames(dx))]
xx$LUF_NAME <- dx$LUF_NAME[match(rownames(xx), rownames(dx))]

xx$OnOffGrid <- as.factor(ifelse(xx$On_Off == "On-Grid", "IG", "OG"))
#xx$ABMI_Assigned_Site_ID[xx$On_Off != "On-Grid"]
xx$DataProvider <- sapply(strsplit(as.character(xx$ABMI_Assigned_Site_ID), "-"),
    function(z) {
        if (length(z) == 1)
            "ABMI" else z[2]
    })
xx$Label2 <- as.factor(with(xx, paste("T", OnOffGrid, DataProvider,
    ABMI_Assigned_Site_ID, survey_year, 1, sep="_")))

if (FALSE) {
e <- new.env()
load(file=file.path(ROOT, VER, "data", "analysis", "species",
    #"OUT_mites_2017-04-24.Rdata"), envir=e)
    "OUT_vplants_2017-04-13.Rdata"), envir=e)
zzz <- as.character(e$res2$Label2)
zzz <- paste0(substr(zzz, 1, nchar(zzz)-2), "_1")
compare_sets(xx$Label2, zzz)
setdiff(xx$Label2, zzz)
setdiff(zzz, xx$Label2)
}

## XY for sites based on nearest IG site
xx$NearestOnGridSite <- 0
xx$isBsite <- FALSE
for (i in 1:nrow(xx)) {
    tmp <- as.character(xx$ABMI_Assigned_Site_ID[i])
    if (xx$OnOffGrid[i] == "OG") {
        xx$NearestOnGridSite[i] <- as.integer(strsplit(tmp, "-")[[1]][3])
    } else {
        if (endsWith(tmp, "B")) {
            xx$NearestOnGridSite[i] <- as.integer(substr(tmp, 1, nchar(tmp)-1))
            xx$isBsite[i] <- TRUE
        } else {
            xx$NearestOnGridSite[i] <- as.integer(tmp)
        }
    }
}
summary(xx)
xx[is.na(xx$NearestOnGridSite),]
xx <- data.frame(xx, meta[match(xx$NearestOnGridSite, meta$SITE_ID),])
summary(xx)

#climSite$Site<-ifelse(climSite$Site=="OG-ABMI-1054-11","OG-ABMI-1054-2",climSite$Site)  # Apparent typo
#climSite$Site<-ifelse(climSite$Site=="OG-ABMI-1122-11","OG-ABMI-1122-2",climSite$Site)  # Apparent typo
#climSite$Site<-ifelse(climSite$Site=="OG-ABMI-1190-11","OG-ABMI-1190-2",climSite$Site)  # Apparent typo
#climSite$Site<-ifelse(climSite$Site=="OG-ABMI-1331-11","OG-ABMI-1331-2",climSite$Site)  # Apparent typo

save(dd_1ha, dd_150m, dd_564m, xx,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_siteCenter_v6-fixage0.Rdata"))

## ABMI Camera/ARU
## 10m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all", "CamARUBird_Buf10m.csv")
d <- read.csv(f)

levels(d$Cam_ARU_Bird_Location) <- c(levels(d$Cam_ARU_Bird_Location), "UNKNOWN")
d$Cam_ARU_Bird_Location[is.na(d$Cam_ARU_Bird_Location)] <- "UNKNOWN"
d$UID <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
table(d$Cam_ARU_Bird_Location,useNA="a")

dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="survey_year",
    col.HFyear="year",
    col.HABIT="Combined_ChgByCWCS",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "10m circle buffer around camera/ARU locations"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd_10m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
## 150m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all", "CamARUBird_Buf150m.csv")
d <- read.csv(f)
d$UID <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="survey_year",
    col.HFyear="year",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
dd$scale <- "150m circle buffer around camera/ARU locations"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd_150m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
## 564m
f <- file.path(ROOT, VER, "data", "raw", "veghf", "site_all", "CamARUBird_Buf564m.csv")
d <- read.csv(f)
d$UID <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year="survey_year",
    col.HFyear="year",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
dd$scale <- "564m circle buffer around camera/ARU locations"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd_564m <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

xx <- read.csv(file.path(ROOT, VER, "data", "raw", "veghf",
    "site_all", "CamARUBird_climate.csv"))
xx$UID <- with(xx, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, survey_year,
    sep="_", drop=TRUE))
rownames(xx) <- xx$UID

compare_sets(rownames(dd_10m[[1]]), rownames(xx))
compare_sets(rownames(dd_150m[[1]]), rownames(xx))
compare_sets(rownames(dd_564m[[1]]), rownames(xx))
setdiff(rownames(xx), rownames(dd_10m[[1]]))

xx <- droplevels(xx[rownames(dd_10m[[1]]),])
xx$PET <- xx$Eref
xx$pAspen <- xx$Populus_tremuloides_brtpred_nofp
xx$Populus_tremuloides_brtpred_nofp <- NULL
xx$OBJECTID <- NULL
xx$NSRNAME <- dx$NSRNAME[match(rownames(xx), rownames(dx))]
xx$NRNAME <- dx$NRNAME[match(rownames(xx), rownames(dx))]
xx$LUF_NAME <- dx$LUF_NAME[match(rownames(xx), rownames(dx))]

save(dd_10m, dd_150m, dd_564m, xx,
    file=file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_CamARUBird_v6-fixage0.Rdata"))

## check classes:

#load(file.path(ROOT, VER, "data", "analysis", "site",
#    "veg-hf_CamARUBird_v6-fixage0.Rdata"))
load(file.path(ROOT, VER, "data", "analysis", "site",
    "veg-hf_siteCenter_v6-fixage0.Rdata"))

zz <- data.frame(x=colMeans(dd_1ha[[1]]))
zz$x <- round(100 * zz$x / sum(zz$x), 2)


d$Site_YEAR <- with(d, interaction(ABMI_ID_WithB, survey_year, sep="_", drop=TRUE))
head(d1ha)
## 2014 site updates
f1hax <- file.path(ROOT, VER, "data/veghf/update2014", "Center1ha_2014.csv")
d1hax <- read.csv(f1hax)
d1hax$Site_YEAR <- with(d1hax, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d1hax[,colnames(d1ha)])
d1ha <- d1ha[d1ha$survey_year != 2014,]
d1ha <- rbind(d1ha, d1hax[,colnames(d1ha)])

dd1ha <- make_vegHF_wide(d1ha, col.label = "Site_YEAR",
    col.year="survey_year", col.HFyear="year_")
dd1ha$scale <- "1 ha square around site centre"

## ABMI sites (on+off) 9 bird points / site, 150 m radius buffer
f150m <- file.path(ROOT, VER, "data/veghf", "Bird150mFixFire.csv")
d150m <- read.csv(f150m)
d150m$Site_YEAR_bird <- with(d150m,
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d150m)
## 2014 site updates
f150mx <- file.path(ROOT, VER, "data/veghf/update2014", "Bird150m_2014.csv")
d150mx <- read.csv(f150mx)
d150mx$Site_YEAR_bird <- with(d150mx,
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d150mx[,colnames(d150m)])
d150m <- d150m[d150m$survey_year != 2014,]
d150m <- rbind(d150m, d150mx[,colnames(d150m)])

dd150m <- make_vegHF_wide(d150m, col.label = "Site_YEAR_bird",
    col.year="survey_year", col.HFyear="year_")
dd150m$scale <- "150 m radius circle around bird points"

## ABMI sites (on+off) 9 bird points / site, 1 km^2 buffer
f1km <- file.path(ROOT, VER, "data/veghf", "Bird564mFixFire.csv")
d1km <- read.csv(f1km)
d1km$Site_YEAR_bird <- with(d1km,
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d1km)
## 2014 site updates
f1kmx <- file.path(ROOT, VER, "data/veghf/update2014", "Bird564m_2014.csv")
d1kmx <- read.csv(f1kmx)
d1kmx$Site_YEAR_bird <- with(d1kmx,
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d1kmx[,colnames(d1km)])
d1km <- d1km[d1km$survey_year != 2014,]
d1km <- rbind(d1km, d1kmx[,colnames(d1km)])

dd1km <- make_vegHF_wide(d1km, col.label = "Site_YEAR_bird",
    col.year="survey_year", col.HFyear="year_")
dd1km$scale <- "564 m radius circle around bird points"

#### Climate and regions

## Public coordinates
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rownames(gis) <- gis$SITE_ID

## climate for all bird pts (pt=1 centre for 1ha)
clim1 <- read.csv(file.path(ROOT, VER, "data/climate", "OnOffBirds_climateLUFNrg.csv"))
colnames(clim1)[colnames(clim1) == "Eref"] <- "PET"
colnames(clim1)[colnames(clim1) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
clim1$Site_YEAR_bird <- with(clim1, interaction(Site_ID, Bird, sep="_", drop=TRUE))
rownames(clim1) <- clim1$Site_YEAR_bird
clim1 <- droplevels(clim1[!grepl("OG-ALPAC-SK", rownames(clim1)),])

tmp <- strsplit(as.character(clim1$Site_ID), "_")
clim1$Site <- as.factor(sapply(tmp, "[[", 1))
clim1$Year <- as.integer(sapply(tmp, "[[", 2))
tmp <- strsplit(as.character(clim1$Site), "-")
clim1$Nearest <- sapply(tmp, function(z) if (length(z)>1) z[3] else z)
clim1$DataProvider <- sapply(tmp, function(z) if (length(z)>1) z[2] else "ABMI")
clim1$OnOffGrid <- sapply(tmp, function(z) if (length(z)>1) z[1] else "IG")
clim1$POINT_X <- gis$PUBLIC_LONGITUDE[match(clim1$Nearest, rownames(gis))]
clim1$POINT_Y <- gis$PUBLIC_LATTITUDE[match(clim1$Nearest, rownames(gis))]

clim2 <- droplevels(clim1[clim1$Bird == 1,])
rownames(clim2) <- clim2$Site_ID

compare.sets(rownames(clim2), rownames(dd1ha$veg_current))
setdiff(rownames(clim2), rownames(dd1ha$veg_current))

compare.sets(rownames(clim1), rownames(dd150m$veg_current))
setdiff(rownames(clim1), rownames(dd150m$veg_current))

compare.sets(rownames(clim1), rownames(dd1km$veg_current))
setdiff(rownames(clim1), rownames(dd1km$veg_current))

birds <- read.csv(file.path(ROOT, VER, "out","species",
    "OUT_Birds_Species_PC-Counts_2015-05-22.csv"))
rownames(birds) <- birds$Label
mites <- read.csv(file.path(ROOT, VER, "out","species",
    "OUT_Mites_Species_Site-Binomial_2015-05-22.csv"))
rownames(mites) <- mites$Label2

clim1$Label <- with(clim1, paste0("T_", OnOffGrid, "_", DataProvider,
    "_", Site, "_", Year, "_1_PT_", Bird))
clim1$Label2 <- with(clim1, paste0("T_", OnOffGrid, "_", DataProvider,
    "_", Site, "_", Year, "_1"))

clim1 <- clim1[rownames(dd150m$veg_current),]
clim2 <- clim2[rownames(dd1ha$veg_current),]
stopifnot(all(rownames(clim1) == rownames(dd150m$veg_current)))
stopifnot(all(rownames(clim2) == rownames(dd1ha$veg_current)))

rownames(clim1) <- clim1$Label
rownames(clim2) <- clim2$Label2
rownames(dd150m[[1]]) <- rownames(dd150m[[2]]) <- rownames(clim1)
rownames(dd150m[[3]]) <- rownames(dd150m[[4]]) <- rownames(clim1)
rownames(dd1km[[1]]) <- rownames(dd1km[[2]]) <- rownames(clim1)
rownames(dd1km[[3]]) <- rownames(dd1km[[4]]) <- rownames(clim1)
rownames(dd1ha[[1]]) <- rownames(dd1ha[[2]]) <- rownames(clim2)
rownames(dd1ha[[3]]) <- rownames(dd1ha[[4]]) <- rownames(clim2)

climPoint <- clim1
climSite <- clim2

compare.sets(rownames(climPoint), rownames(birds))
setdiff(rownames(climPoint), rownames(birds))
setdiff(rownames(birds), rownames(climPoint))
compare.sets(rownames(climSite), as.character(birds$Label2))
setdiff(rownames(climSite), as.character(birds$Label2))
setdiff(as.character(birds$Label2), rownames(climSite))

compare.sets(rownames(climSite), rownames(mites))
setdiff(rownames(climSite), rownames(mites))
setdiff(rownames(mites), rownames(climSite))

## topo variables
topo <- read.csv(file.path(ROOT, VER, "data/topo", "ABMIBirdsCamARU_topo.csv"))
topo$Site_YEAR_bird <- with(topo, interaction(Site_ID, Cam_ARU_Bird_Location, sep="_", drop=TRUE))
compare.sets(climPoint$Site_YEAR_bird, topo$Site_YEAR_bird)

topo1 <- topo[match(climPoint$Site_YEAR_bird, topo$Site_YEAR_bird),]
topo2 <- topo[topo$Cam_ARU_Bird_Location == "1",]
compare.sets(climSite$Site_ID, topo$Site_ID)
topo2 <- topo2[match(climSite$Site_ID, topo2$Site_ID),]


climPoint$SLP <- topo1$slope
climPoint$ASP <- topo1$slpasp
climPoint$TRI <- topo1$tri
climPoint$CTI <- topo1$cti

climSite$SLP <- topo2$slope
climSite$ASP <- topo2$slpasp
climSite$TRI <- topo2$tri
climSite$CTI <- topo2$cti

## fix age 0 in saved files -----------------------------

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

sum(dd1ha[[1]][,Target0])
sum(dd1ha[[1]])
dd1ha <- fill_in_0ages(dd1ha, climSite$NSRNAME)
sum(dd1ha[[1]][,Target0])
sum(dd1ha[[1]])

sum(dd150m[[1]][,Target0])
dd150m <- fill_in_0ages(dd150m, climPoint$NSRNAME)
sum(dd150m[[1]][,Target0])

sum(dd1km[[1]][,Target0])
dd1km <- fill_in_0ages(dd1km, climPoint$NSRNAME)
sum(dd1km[[1]][,Target0])

if (SAVE)
    save(dd1ha, dd150m, dd1km, climSite, climPoint,
        file=file.path(ROOT, VER, "out/abmi_onoff",
        "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0.Rdata"))

## 2015 locations

#e:/peter/AB_data_v2016/data/veghf/update2015/BirdCamARU150m.csv
#e:/peter/AB_data_v2016/data/veghf/update2015/BirdCamARU564m.csv
#e:/peter/AB_data_v2016/data/veghf/update2015/BirdCamARUPoints.csv
#e:/peter/AB_data_v2016/data/veghf/update2015/Site150m.csv
#e:/peter/AB_data_v2016/data/veghf/update2015/Site1ha.csv
#e:/peter/AB_data_v2016/data/veghf/update2015/Site564m.csv
#e:/peter/AB_data_v2016/data/veghf/update2015/SitePoints.csv

## ABMI sites (on+off) cetre 1 ha
f1ha <- file.path(ROOT, VER, "data", "veghf", "update2015", "Site1ha.csv")
d1ha <- read.csv(f1ha)
d1ha$Site_YEAR <- with(d1ha, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d1ha)
dd1ha <- make_vegHF_wide(d1ha, col.label = "Site_YEAR",
    col.year="survey_year", col.HFyear="year_")
dd1ha$scale <- "1 ha square around site centre"
dd1ha_2015 <- dd1ha

## ABMI sites (on+off) 9 bird points / site, 150 m radius buffer, site center only
f150m <- file.path(ROOT, VER, "data", "veghf", "update2015", "Site150m.csv")
d150m <- read.csv(f150m)
d150m$Site_YEAR <- with(d150m, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d150m)
dd150m <- make_vegHF_wide(d150m, col.label = "Site_YEAR",
    col.year="survey_year", col.HFyear="year_")
dd150m$scale <- "150 m radius circle around site centre"
dd150mCenter_2015 <- dd150m

## ABMI sites (on+off) 9 bird points / site, 1 km^2 buffer, site center only
f1km <- file.path(ROOT, VER, "data", "veghf", "update2015", "Site564m.csv")
d1km <- read.csv(f1km)
d1km$Site_YEAR <- with(d1km, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d1km)
dd1km <- make_vegHF_wide(d1km, col.label = "Site_YEAR",
    col.year="survey_year", col.HFyear="year_")
dd1km$scale <- "564 m radius circle around site centre"
dd1kmCenter_2015 <- dd1km

## ABMI bird/camera/ARU points, 150 m radius buffer
f150m <- file.path(ROOT, VER, "data", "veghf", "update2015", "BirdCamARU150m.csv")
d150m <- read.csv(f150m)
d150m$survey_year <- 2015
d150m$Site_YEAR_PT <- with(d150m, interaction(
    ABMI_Assigned_Site_ID,
    survey_year,
    deployment,
    Cam_ARU_Bird_Location, sep="_", drop=TRUE))
head(d150m)
dd150m <- make_vegHF_wide(d150m, col.label = "Site_YEAR_PT",
    col.year="survey_year", col.HFyear="year_")
dd150m$scale <- "150 m radius circle around bird/Camera/ARU points"
dd150mPT_2015 <- dd150m

## ABMI bird/camera/ARU points, 1 km^2 buffer
f1km <- file.path(ROOT, VER, "data", "veghf", "update2015", "BirdCamARU564m.csv")
d1km <- read.csv(f1km)
d1km$survey_year <- 2015
d1km$Site_YEAR_PT <- with(d1km, interaction(
    ABMI_Assigned_Site_ID,
    survey_year,
    deployment,
    Cam_ARU_Bird_Location, sep="_", drop=TRUE))
head(d1km)
dd1km <- make_vegHF_wide(d1km, col.label = "Site_YEAR_PT",
    col.year="survey_year", col.HFyear="year_")
dd1km$scale <- "564 m radius circle around bird/Camera/ARU points"
dd1kmPT_2015 <- dd1km

rm(dd1ha, dd150m, dd1km)

load(file.path(ROOT, VER, "out/abmi_onoff",
    "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0.Rdata"))

all(rownames(dd150mPT_2015[[1]]) == rownames(dd1kmPT_2015[[1]]))
all(rownames(dd1ha_2015[[1]]) == rownames(dd150mCenter_2015[[1]]))
all(rownames(dd1ha_2015[[1]]) == rownames(dd1kmCenter_2015[[1]]))

## Public coordinates
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rownames(gis) <- gis$SITE_ID

## climate for all bird pts (pt=1 centre for 1ha)
clim1 <- read.csv(file.path(ROOT, VER, "data/climate/SitePoints_climate-2015.csv"))
colnames(clim1)[colnames(clim1) == "Eref"] <- "PET"
colnames(clim1)[colnames(clim1) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
clim1$Site_YEAR <- with(clim1, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
rownames(clim1) <- clim1$Site_YEAR

tmp <- strsplit(as.character(clim1$Site_YEAR), "_")
clim1$Site <- as.factor(sapply(tmp, "[[", 1))
clim1$Year <- as.integer(sapply(tmp, "[[", 2))
tmp <- strsplit(as.character(clim1$Site), "-")
clim1$Nearest <- sapply(tmp, function(z) if (length(z)>1) z[3] else z)
clim1$Nearest <- gsub("B", "", clim1$Nearest, fixed=TRUE) # B sites

clim1$DataProvider <- sapply(tmp, function(z) if (length(z)>1) z[2] else "ABMI")
clim1$OnOffGrid <- sapply(tmp, function(z) if (length(z)>1) z[1] else "IG")
clim1$POINT_X <- gis$PUBLIC_LONGITUDE[match(clim1$Nearest, rownames(gis))]
clim1$POINT_Y <- gis$PUBLIC_LATTITUDE[match(clim1$Nearest, rownames(gis))]

setdiff(colnames(clim1), colnames(climSite))
setdiff(colnames(climSite), colnames(clim1))

clim1$Label2 <- with(clim1, paste0("T_", OnOffGrid, "_", DataProvider,
    "_", Site, "_", Year, "_1"))
climCenter_2015 <- clim1

compare_sets(rownames(dd1ha_2015[[1]]), rownames(climCenter_2015))
climCenter_2015 <- climCenter_2015[rownames(dd1ha_2015[[1]]),]

## CAM/ARU
clim1 <- read.csv(file.path(ROOT, VER, "data/climate/BirdCamARUPoints_climate-2015.csv"))

colnames(clim1)[colnames(clim1) == "Eref"] <- "PET"
colnames(clim1)[colnames(clim1) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
clim1$Site_YEAR_CAMARU <- with(clim1, interaction(ABMI_ID_WithB,
    2015, deployment, Cam_ARU_Location, sep="_", drop=TRUE))
clim1 <- nonDuplicated(clim1, Site_YEAR_CAMARU, TRUE)
compare_sets(rownames(dd150mPT_2015[[1]]), rownames(clim1))
clim1 <- clim1[rownames(dd150mPT_2015[[1]]),]

tmp <- strsplit(as.character(clim1$Site_YEAR_CAMARU), "_")
clim1$Site <- as.factor(sapply(tmp, "[[", 1))
clim1$Year <- as.integer(sapply(tmp, "[[", 2))
tmp <- strsplit(as.character(clim1$Site), "-")
clim1$Nearest <- sapply(tmp, function(z) if (length(z)>1) z[3] else z)
clim1$Nearest <- gsub("B", "", clim1$Nearest, fixed=TRUE) # B sites

clim1$DataProvider <- sapply(tmp, function(z) if (length(z)>1) z[2] else "ABMI")
clim1$OnOffGrid <- sapply(tmp, function(z) if (length(z)>1) z[1] else "IG")
clim1$POINT_X <- gis$PUBLIC_LONGITUDE[match(clim1$Nearest, rownames(gis))]
clim1$POINT_Y <- gis$PUBLIC_LATTITUDE[match(clim1$Nearest, rownames(gis))]

clim1$Label <- with(clim1, paste0("T_", OnOffGrid, "_", DataProvider,
    "_", Site, "_", Year, "_1_", deployment, "_", Cam_ARU_Location))
clim1$Label2 <- with(clim1, paste0("T_", OnOffGrid, "_", DataProvider,
    "_", Site, "_", Year, "_1"))
climPT_2015 <- clim1

if (FALSE) {
## topo variables
topo1 <- read.csv(file.path(ROOT, VER, "data/topo/ABMIBirdsCamARU_topo.csv"))
topo1$Site_YEAR_CAMARU <- with(topo,
    interaction(Site_ID, deployment, Cam_ARU_Bird_Location, sep="_", drop=TRUE))
compare.sets(climPoint$Site_YEAR_bird, topo$Site_YEAR_bird)

topo1 <- topo[match(clim1$Site_YEAR_CAMARU, topo$Site_YEAR_CAMARU),]

clim1$SLP <- topo1$slope
clim1$ASP <- topo1$slpasp
clim1$TRI <- topo1$tri
clim1$CTI <- topo1$cti
}


sum(is.na(climCenter_2015))
sum(is.na(climPT_2015))

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

sum(dd1ha_2015[[1]][,Target0])
dd1ha_2015 <- fill_in_0ages(dd1ha_2015, climCenter_2015$NSRNAME)
sum(dd1ha_2015[[1]][,Target0])

sum(dd150mCenter_2015[[1]][,Target0])
dd150mCenter_2015 <- fill_in_0ages(dd150mCenter_2015, climCenter_2015$NSRNAME)
sum(dd150mCenter_2015[[1]][,Target0])

sum(dd1kmCenter_2015[[1]][,Target0])
dd1kmCenter_2015 <- fill_in_0ages(dd1kmCenter_2015, climCenter_2015$NSRNAME)
sum(dd1kmCenter_2015[[1]][,Target0])

sum(dd150mPT_2015[[1]][,Target0])
dd150mPT_2015 <- fill_in_0ages(dd150mPT_2015, climPT_2015$NSRNAME)
sum(dd150mPT_2015[[1]][,Target0])

sum(dd1kmPT_2015[[1]][,Target0])
dd1kmPT_2015 <- fill_in_0ages(dd1kmPT_2015, climPT_2015$NSRNAME)
sum(dd1kmPT_2015[[1]][,Target0])


if (SAVE)
    save(dd1ha, dd150m, dd1km, climSite, climPoint,
        dd1ha_2015, dd150mCenter_2015, dd1kmCenter_2015,
        dd150mPT_2015, dd1kmPT_2015, climCenter_2015, climPT_2015,
        file=file.path(ROOT, VER, "out/abmi_onoff",
        "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015.Rdata"))

## 2003-2013 revisit updates

load(file.path(ROOT, VER, "out/abmi_onoff",
    "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015.Rdata"))

## revisited ABMI sites (on+off) cetre 1 ha
f1ha <- file.path(ROOT, VER, "data/veghf/revisits", "Center1Ha_reDoYear2003_2013For2015RevisitedSites_tbl.csv")
d1ha <- read.csv(f1ha)
d1ha$Site_YEAR <- with(d1ha, interaction(ABMI_ID, survey_year, sep="_", drop=TRUE))
d1ha$Site_YEAR2 <- rownames(climSite)[match(d1ha$Site_YEAR, climSite$Site_ID)]
head(d1ha)
ddr1ha <- make_vegHF_wide(d1ha, col.label = "Site_YEAR2",
    col.year="survey_year", col.HFyear="year_")
ddr1ha$scale <- "1 ha square around site centre"
compare_sets(rownames(dd1ha[[1]]), rownames(ddr1ha[[1]]))

## revisited ABMI sites (on+off) 9 bird points / site, 150 m radius buffer
f150m <- file.path(ROOT, VER, "data/veghf/revisits",
    "birds150m_reDoYear2003_2013For2015RevisitedSites_tbl.csv")
d150m <- read.csv(f150m)
levels(d150m$Site_ID) <- gsub("-", "_", levels(d150m$Site_ID), fixed=TRUE)
d150m$Site_YEAR_bird <- with(d150m,
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
d150m$Site_YEAR_bird2 <- rownames(climPoint)[match(d150m$Site_YEAR_bird,
    climPoint$Site_YEAR_bird)]
head(d150m)
ddr150m <- make_vegHF_wide(d150m, col.label = "Site_YEAR_bird2",
    col.year="survey_year", col.HFyear="year_")
ddr150m$scale <- "150 m radius circle around bird points"
compare_sets(rownames(dd150m[[1]]), rownames(ddr150m[[1]]))

## revisited ABMI sites (on+off) 9 bird points / site, 1 km^2 buffer
f1km <- file.path(ROOT, VER, "data/veghf/revisits",
    "birds564m_reDoYear2003_2013For2015RevisitedSites_tbl.csv")
d1km <- read.csv(f1km)
levels(d1km$Site_ID) <- gsub("-", "_", levels(d1km$Site_ID), fixed=TRUE)
d1km$Site_YEAR_bird <- with(d1km,
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
d1km$Site_YEAR_bird2 <- rownames(climPoint)[match(d1km$Site_YEAR_bird,
    climPoint$Site_YEAR_bird)]
head(d1km)
ddr1km <- make_vegHF_wide(d1km, col.label = "Site_YEAR_bird2",
    col.year="survey_year", col.HFyear="year_")
ddr1km$scale <- "564 m radius circle around bird points"
compare_sets(rownames(dd1km[[1]]), rownames(ddr1km[[1]]))

## fix fire

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

sum(ddr1ha[[1]][,Target0])
sum(ddr1ha[[1]])
ddr1ha <- fill_in_0ages(ddr1ha, climSite[rownames(ddr1ha[[1]]), "NSRNAME"])
sum(ddr1ha[[1]][,Target0])
sum(ddr1ha[[1]])

sum(ddr150m[[1]][,Target0])
ddr150m <- fill_in_0ages(ddr150m, climPoint[rownames(ddr150m[[1]]), "NSRNAME"])
sum(ddr150m[[1]][,Target0])

sum(ddr1km[[1]][,Target0])
ddr1km <- fill_in_0ages(ddr1km, climPoint[rownames(ddr1km[[1]]), "NSRNAME"])
sum(ddr1km[[1]][,Target0])

for (i in 1:4) {
    cat("1ha", i, "\n")
    print(summary(colMeans(abs(dd1ha[[i]][rownames(ddr1ha[[i]]),] - ddr1ha[[i]]))))
    cat("150m", i, "\n")
    print(summary(colMeans(abs(dd150m[[i]][rownames(ddr150m[[i]]),] - ddr150m[[i]]))))
    cat("1km", i, "\n")
    print(summary(colMeans(abs(dd1km[[i]][rownames(ddr1km[[i]]),] - ddr1km[[i]]))))
}

for (i in 1:4) {
    dd1ha[[i]][rownames(ddr1ha[[i]]),] <- ddr1ha[[i]]
    dd150m[[i]][rownames(ddr150m[[i]]),] <- ddr150m[[i]]
    dd1km[[i]][rownames(ddr1km[[i]]),] <- ddr1km[[i]]
}

for (i in 1:4) {
    cat("1ha", i, "\n")
    print(summary(colMeans(abs(dd1ha[[i]][rownames(ddr1ha[[i]]),] - ddr1ha[[i]]))))
    cat("150m", i, "\n")
    print(summary(colMeans(abs(dd150m[[i]][rownames(ddr150m[[i]]),] - ddr150m[[i]]))))
    cat("1km", i, "\n")
    print(summary(colMeans(abs(dd1km[[i]][rownames(ddr1km[[i]]),] - ddr1km[[i]]))))
}

## might want to add point level intersection here?

if (SAVE)
    save(dd1ha, dd150m, dd1km, climSite, climPoint,
        dd1ha_2015, dd150mCenter_2015, dd1kmCenter_2015,
        dd150mPT_2015, dd1kmPT_2015, climCenter_2015, climPT_2015,
        file=file.path(ROOT, VER, "out/abmi_onoff",
        "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015-revisits.Rdata"))

## merging objects at site center

source("~/repos/abmianalytics/veghf/veghf-setup.R")
e <- new.env()
load(file.path(ROOT, VER, "out/abmi_onoff",
    "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015-revisits.Rdata"),
    envir=e)

climSite <- e$climSite
climCenter_2015 <- e$climCenter_2015
compare_sets(colnames(climSite), colnames(climCenter_2015))
setdiff(colnames(climSite), colnames(climCenter_2015))
setdiff(colnames(climCenter_2015), colnames(climSite))

climSite$Site_ID <- NULL
climSite$Bird <- NULL
climSite$Site_YEAR_bird <- NULL
climSite$SLP <- NULL
climSite$ASP <- NULL
climSite$TRI <- NULL
climSite$CTI <- NULL

climCenter_2015$Site_YEAR <- NULL
climCenter_2015$On_Off <- NULL
climCenter_2015$survey_year <- NULL
climCenter_2015$ABMI_Assigned_Site_ID <- NULL

climCenter_2015 <- climCenter_2015[,colnames(climSite)]

dd1ha <- e$dd1ha
dd1ha_2015 <- e$dd1ha_2015
all(rownames(climCenter_2015) == rownames(dd1ha_2015[[1]]))

dd1km <- e$dd1km
dd1kmCenter_2015 <- e$dd1kmCenter_2015
all(rownames(climCenter_2015) == rownames(dd1kmCenter_2015[[1]]))

dd1kmCenter <- dd1km
climPoint <- e$climPoint
for (i in 1:4) {
    dd1kmCenter[[i]] <- dd1km[[i]][match(rownames(climSite), climPoint$Label2),]
    rownames(dd1kmCenter[[i]]) <- rownames(climSite)
}

for (i in 1:4) {
    rownames(dd1ha_2015[[i]]) <- climCenter_2015$Label2
    rownames(dd1kmCenter_2015[[i]]) <- climCenter_2015$Label2
}
rownames(climCenter_2015) <- climCenter_2015$Label2
all(rownames(climCenter_2015) == rownames(dd1ha_2015[[1]]))
all(rownames(climCenter_2015) == rownames(dd1kmCenter_2015[[1]]))

climSite$Site <- as.character(climSite$Site)
climCenter_2015$Site <- as.character(climCenter_2015$Site)
climSite$MAP <- as.numeric(gsub(",", "", as.character(climSite$MAP)))

for (i in 1:4) {
    dd1ha[[i]] <- rBind(dd1ha[[i]], dd1ha_2015[[i]])
    dd1kmCenter[[i]] <- rBind(dd1kmCenter[[i]], dd1kmCenter_2015[[i]])
}
climSite <- rbind(climSite, climCenter_2015)
all(rownames(climSite) == rownames(dd1ha[[1]]))
all(rownames(climSite) == rownames(dd1kmCenter[[1]]))

save(dd1ha, dd1kmCenter, climSite,
    file=file.path(ROOT, VER, "out/abmi_onoff",
    "veg-hf-clim-reg_abmi-onoff_siteCentre_incl2015.Rdata"))

## merging objects at bird points and ARU locations

source("~/repos/abmianalytics/veghf/veghf-setup.R")
e <- new.env()
load(file.path(ROOT, VER, "out/abmi_onoff",
    "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015-revisits.Rdata"),
    envir=e)
climPT_2015 <- e$climPT_2015
climPT_2015$Label0 <- climPT_2015$Label

topo1 <- read.csv(file.path(ROOT, VER, "data", "topo", "ABMIBirdsCamARU_topo.csv"))
topo1$ID <- with(topo1, paste0(Site_ID, "_", deployment, "_", Cam_ARU_Bird_Location))
compare_sets(topo1$ID, rownames(climPT_2015))
topo1 <- topo1[match(rownames(climPT_2015), topo1$ID),]
climPT_2015$SLP <- topo1$slope
climPT_2015$ASP <- topo1$slpasp
climPT_2015$TRI <- topo1$tri
climPT_2015$CTI <- topo1$cti

## RiverForks points
climPT1 <- droplevels(climPT_2015[climPT_2015$deployment == "Bird",])
climPT1$Label <- gsub("_Bird_", "_PT_", climPT1$Label)
## ARU points
climPT2 <- droplevels(climPT_2015[climPT_2015$deployment %in% c("AppCenter", "ARU",
    "BOTH", "SciCenter"),])
climPT2$Label <- gsub("_AppCenter_", "_STATION_", climPT2$Label)
climPT2$Label <- gsub("_ARU_", "_STATION_", climPT2$Label)
climPT2$Label <- gsub("_BOTH_", "_STATION_", climPT2$Label)
climPT2$Label <- gsub("_SciCenter_", "_STATION_", climPT2$Label)
climPT2$Label <- gsub("T_IG_ABMI_388B_2015_1", "T_IG_ABMI_388_2015_1", climPT2$Label)
## Camera points
climPT3 <- droplevels(climPT_2015[climPT_2015$deployment %in% c("CAM", "BOTH"),])

## climate etc for RiverForks
e$climPoint$Label0 <- e$climPoint$Label
compare_sets(colnames(e$climPoint), colnames(climPT1))
setdiff(colnames(e$climPoint), colnames(climPT1))
setdiff(colnames(climPT1), colnames(e$climPoint))
(cn <- intersect(colnames(e$climPoint), colnames(climPT1)))
climRF <- rbind(e$climPoint[,cn], climPT1[,cn])
## climate etc for SM units
climSM <- climPT2[,cn]

dd150m_RF <- list()
dd150m_SM <- list()
dd1km_RF <- list()
dd1km_SM <- list()

for (i in 1:4) {
    dd150m_RF[[i]] <- rbind(e$dd150m[[i]], e$dd150mPT_2015[[i]])
    dd150m_RF[[i]] <- dd150m_RF[[i]][rownames(climRF), ]
    rownames(dd150m_RF[[i]]) <- climRF$Label

    dd150m_SM[[i]] <- e$dd150mPT_2015[[i]]
    dd150m_SM[[i]] <- dd150m_SM[[i]][rownames(climSM), ]
    rownames(dd150m_SM[[i]]) <- climSM$Label

    dd1km_RF[[i]] <- rbind(e$dd1km[[i]], e$dd1kmPT_2015[[i]])
    dd1km_RF[[i]] <- dd1km_RF[[i]][rownames(climRF), ]
    rownames(dd1km_RF[[i]]) <- climRF$Label

    dd1km_SM[[i]] <- e$dd1kmPT_2015[[i]]
    dd1km_SM[[i]] <- dd1km_SM[[i]][rownames(climSM), ]
    rownames(dd1km_SM[[i]]) <- climSM$Label
}
dd150m_RF$scale <- dd150m_SM$scale <- e$dd150m$scale
dd1km_RF$scale <- dd1km_SM$scale <- e$dd1km$scale

rownames(climRF) <- climRF$Label
rownames(climSM) <- climSM$Label

save(dd150m_RF, dd150m_SM, dd1km_RF, dd1km_SM, climRF, climSM,
    file=file.path(ROOT, VER, "out/abmi_onoff",
    "veg-hf-clim-reg_abmi-onoff_Birds-RF-SM_incl2015.Rdata"))

load("e:/peter/AB_data_v2016/data/species/OUT_birdssm_2016-05-27.Rdata")

zz <- rownames(climSM)
zz <- gsub("_Center", "_1", zz)
compare_sets(zz, x$SITE_LABEL)
setdiff(zz, x$SITE_LABEL)
setdiff(x$SITE_LABEL, zz)

## -- requests --

## 2018-09-28
## a) Verified Human Footprint (in the year of sampling) summarized tabular
## data for a 564m buffer zone (1 km2) surrounding sampled terrestrial ABMI
## systematic grid site centres in the Oil Sands Monitoring Region.
## b) Verified Human Footprint (in the year of sampling) summarized tabular
## data for a 250m buffer zone surrounding open water at sampled wetland ABMI
## systematic grid sites in the Oils Sands Monitoring Region.

## dd_564m, xx
load(file.path(ROOT, VER, "data", "analysis", "site", "veg-hf_SiteCenter_v6verified.Rdata"))

id <- readLines(file.path(ROOT, VER, "data", "raw", "xy", "ABMISites_InOSM.txt"))
ss <- as.character(xx$ABMI_Assigned_Site_ID) %in% id

xxs <- xx
xxs$inOSR <- ss
xxs$NRNAME_real <- dd_point$NRNAME[match(xxs$Site_YEAR, dd_point$Site_YEAR)]
xxs$NSRNAME_real <- dd_point$NSRNAME[match(xxs$Site_YEAR, dd_point$Site_YEAR)]
xxs$LUF_NAME_real <- dd_point$LUF_NAME[match(xxs$Site_YEAR, dd_point$Site_YEAR)]

vhs <- as.matrix(dd_564m[[1]][,setdiff(colnames(dd_564m[[1]]), colnames(dd_564m[[2]]))])
vhs <- vhs / rowSums(dd_564m[[1]])
For <- rowSums(vhs[,startsWith(colnames(vhs), "CC")])
vhs <- cbind(vhs[,!startsWith(colnames(vhs), "CC")], "ForestHarvest"=For)
vhs_terr <- vhs
vhs_terr <- data.frame(xxs[,c("inOSR", "NRNAME_real", "NSRNAME_real", "LUF_NAME_real")], vhs_terr)

## clim, ddw_250m
load(file.path(ROOT, VER, "data", "analysis", "site", "veg-hf_wetlands_v6x.Rdata"))

rownames(clim) <- clim$SiteYear
xx <- clim[intersect(rownames(clim), rownames(ddw_250m[[1]])),]
ss <- as.character(xx$Site_ID_Ref) %in% paste0("W", id)

xxs <- xx
xxs$inOSR <- ss
xxs$NRNAME_real <- ddw_point$NRNAME[match(xxs$SiteYear, ddw_point$SiteYear)]
xxs$NSRNAME_real <- ddw_point$NSRNAME[match(xxs$SiteYear, ddw_point$SiteYear)]
xxs$LUF_NAME_real <- ddw_point$LUF_NAME[match(xxs$SiteYear, ddw_point$SiteYear)]

vhs <- as.matrix(ddw_250m[[1]][,setdiff(colnames(ddw_250m[[1]]), colnames(ddw_250m[[2]]))])
vhs <- vhs / rowSums(ddw_250m[[1]])
For <- rowSums(vhs[,startsWith(colnames(vhs), "CC")])
vhs <- cbind(vhs[,!startsWith(colnames(vhs), "CC")], "ForestHarvest"=For)
vhs_wet <- vhs[rownames(xxs),]
vhs_wet <- data.frame(xxs[,c("inOSR", "NRNAME_real", "NSRNAME_real", "LUF_NAME_real")], vhs_wet)

write.csv(vhs_terr, file=file.path(ROOT, VER, "data", "misc", "requests",
    "2018-09-28", "Terrestrial-ABMI-all-sites-HF-564mBuffer.csv"))
write.csv(vhs_wet, file=file.path(ROOT, VER, "data", "misc", "requests",
    "2018-09-28", "Wetland-ABMI-all-sites-HF-250mBuffer.csv"))


## 2018-11-30
##
## Hello.  I'm doing the 3x7 veg+HF trend summary (with forestry recovery, CI's, etc)
## for the three 7G areas: lease, watershed and five watershed areas.
## I have the updated 3x7's from Peter - the same ones as for the OSR report.
## But I don't have the 2016 wall-to-wall summarized for the 7G area, which I will
## use to standardardize the 3x7 values.
## The attached is the 2016 W2W summary that Peter sent for each area in the OSR.
## Could I get the exact same thing for the three 7G areas?
## I think Peter did this from the 1km2 rasters -
## do you have the info on which of these rasters are in each of the three 7G
## summary areas?  Or, Katherine, could you work with Peter to get him those?

library(mefa4)
library(readxl)
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3noDistVeg.Rdata")
f <- "d:/abmi/AB_data_v2018/data/raw/xy/GridValues7GAreas_Grid_1km.xlsx"
x1 <- read_excel(f, sheet="LeaseArea7G_Grid_1km")
x2 <- read_excel(f, sheet="Watershed07GB_Grid_1km")
x3 <- read_excel(f, sheet="FiveWatershedin7G_Grid_1km")

out <- rbind(
    LeaseArea7G=colSums(dd_kgrid[[1]][x1$GRID_LABEL[x1[["Shape_Area (m^2)"]] > 0.5*10^6],]),
    Watershed07GB=colSums(dd_kgrid[[1]][x2$GRID_LABEL[x2[["Shape_Area (m^2)"]] > 0.5*10^6],]),
    FiveWatershedin7G=colSums(dd_kgrid[[1]][x3$GRID_LABEL[x3[["Shape_Area (m^2)"]] > 0.5*10^6],]))

write.csv(as.data.frame(out), file=file.path("d:/abmi/AB_data_v2018", "data", "misc", "requests",
    "2018-11-30", "VegHF2016-GridValues7GAreas.csv"))
