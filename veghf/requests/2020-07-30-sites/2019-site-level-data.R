## putting together sile level data up to 2019
library(mefa4)
library(raster)
library(sp)

## 1ha and qha 2003-2018:
f1 <- "d:/abmi/AB_data_v2019/data/analysis/site/veg-hf_SiteCenter_Veg61-vHF.Rdata"
e1 <- new.env()
load(f1, envir=e1)

## climate 2003-2016
f2 <- "d:/abmi/AB_data_v2017/data/raw/veghf/site_all/siteCenter_climate.csv"
cl0316 <- read.csv(f2)
rownames(cl0316) <- paste0(cl0316$ABMI_Assigned_Site_ID, "_", cl0316$survey_year)

## all up to 2017
f3 <- "d:/abmi/AB_data_v2018/data/analysis/site/veg-hf_SiteCenter_v6verified.Rdata"
e3 <- new.env()
load(f3, envir=e3)

## 2019 updates
f4 <- "s:/AB_data_v2020/data/analysis/site/veg-hf_SITE1HA-2019_Veg61-vHF.Rdata"
e4 <- new.env()
load(f4, envir=e4)

xy <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

f5 <- "d:/abmi/AB_data_v2018/data/analysis/site/veg-hf_CameraARU_v6verif_2017-2018-sites.Rdata"
e5 <- new.env()
load(f5, envir=e5)


table(e1$xx$survey_year)

## site IDs and regions
SITES <- e1$xx[,1:6]
s19 <- nonDuplicated(e4$d_long, UID_old, TRUE)
s19 <- s19[,c("UID_old", "ABMI_ID_WithB", "survey_year", "NSRNAME", "NRNAME", "LUF_NAME")]
colnames(s19) <- colnames(SITES)
SITES <- rbind(SITES, s19)
table(SITES$survey_year)

SITES$offgrid <- grepl("OG-", SITES$site_id) | grepl("OGW-", SITES$site_id)
table(SITES$offgrid)
SITES$nearest <- SITES$site_id

tmp <- SITES$site_id[SITES$offgrid]
tmp <- gsub("Confidential-", "", tmp)
tmp <- strsplit(tmp, "-")
SITES$nearest[SITES$offgrid] <- sapply(tmp, "[[", 3)

SITES$bsite <- endsWith(SITES$nearest, "B")
table(SITES$bsite)
SITES$nearest[SITES$bsite] <- gsub("B", "", SITES$nearest[SITES$bsite])
SITES$nearest <- as.integer(SITES$nearest)
summary(SITES$nearest)

SITES$X <- xy$PUBLIC_LONGITUDE[match(SITES$nearest, xy$SITE_ID)]
SITES$Y <- xy$PUBLIC_LATTITUDE[match(SITES$nearest, xy$SITE_ID)]


## available climate
str(e4$clim)
head(e4$clim)
str(cl0316)
cl0316$PET <- cl0316$Eref
cl0316$pAspen <- cl0316$Populus_tremuloides_brtpred_nofp
cl <- rbind(e4$clim, cl0316[,colnames(e4$clim)])
rownames(cl) <- gsub("Confidential_", "Confidential-", rownames(cl))

compare_sets(rownames(SITES), rownames(cl))
setdiff(rownames(SITES), rownames(cl))
setdiff(rownames(cl), rownames(SITES))

SITES <- data.frame(SITES, cl[match(rownames(SITES), rownames(cl)),])

summary(SITES)

## try to use previously visited info to be used for 2018 sites
tmp <- nonDuplicated(SITES[!is.na(SITES$AHM),], nearest)
SITES[is.na(SITES$AHM), colnames(cl)] <- SITES[match(SITES$nearest[is.na(SITES$AHM)], tmp$nearest), colnames(cl)]

SITES$clim_source <- ifelse(is.na(SITES$AHM), "public", "exact")
table(SITES$survey_year, SITES$clim_source)

## extract clim values for public xy

rr <- raster("d:/abmi/AB_data_v2016/data/kgrid/AHM1k.tif")
r <- list()
for (i in colnames(cl0316)[6:13]) {
    r[[i]] <- raster(paste0("s:/Base_shapefiles/ab-climate/", i, ".asc"))
}
r <- stack(r)
proj4string(r) <- proj4string(rr)

cr <- SITES[is.na(SITES$AHM),]
coordinates(cr) <- ~ X + Y
proj4string(cr) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
cr <- spTransform(cr, proj4string(rr))

v <- extract(r, cr, method="bilinear")
summary(v)
v <- data.frame(v)
v$PET <- v$Eref
v$pAspen <- v$Populus_tremuloides_brtpred_nofp
SITES[is.na(SITES$AHM), colnames(cl)] <- v[,colnames(cl)]
summary(SITES)

## ----------- veghf 1ha level

names(e1)
names(e3)
names(e4)

i <- 3
all(colnames(e1$dd_1ha[[i]]) == colnames(e3$dd_1ha[[i]]))
all(colnames(e1$dd_1ha[[i]]) == colnames(e4$d_wide_1ha[[i]]))
compare_sets(colnames(e1$dd_1ha[[i]]), colnames(e4$d_wide_1ha[[i]]))

table(e1$dd_1ha$sample_year)
table(e3$dd_1ha$sample_year)
table(e4$d_wide_1ha$sample_year)

## 1ha
dd_1ha <- e1$dd_1ha
for (i in 1:4)
    dd_1ha[[i]] <- rbind(
        e1$dd_1ha[[i]],
        e4$d_wide_1ha[[i]][,colnames(e1$dd_1ha[[i]])])
dd_1ha$scale <- NULL
dd_1ha$sample_year <- c(e1$dd_1ha[[5]], e4$d_wide_1ha[[5]])
compare_sets(rownames(SITES), rownames(dd_1ha[[1]]))
clim_1ha <- SITES[rownames(dd_1ha[[1]]),]
all(rownames(clim_1ha) == rownames(dd_1ha[[1]]))

## qha
dd_qha <- e1$dd_qha
for (i in 1:4)
    dd_qha[[i]] <- rbind(
        e1$dd_qha[[i]],
        e4$d_wide_qha[[i]][,colnames(e1$dd_qha[[i]])])
dd_qha$scale <- NULL
dd_qha$sample_year <- c(e1$dd_qha[[5]], e4$d_wide_qha[[5]])

clim_qha <- rbind(
    data.frame(SITES, quadrant="NE"),
    data.frame(SITES, quadrant="NW"),
    data.frame(SITES, quadrant="SE"),
    data.frame(SITES, quadrant="SW"))
rownames(clim_qha) <- paste0(clim_qha$UID, "_", clim_qha$quadrant)
compare_sets(rownames(clim_qha), rownames(dd_qha[[1]]))
clim_qha <- clim_qha[rownames(dd_qha[[1]]),]

## 1km
dd_564m <- e3$dd_564m
for (i in 1:4)
    dd_564m[[i]] <- rbind(
        e3$dd_564m[[i]],
        e4$d_wide_1km[[i]][,colnames(e3$dd_564m[[i]])])
dd_564m$scale <- NULL
dd_564m$sample_year <- c(e3$dd_564m[[5]], e4$d_wide_1km[[5]])

compare_sets(rownames(SITES), rownames(dd_564m[[1]]))


## use 2018 w2w results???






od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)


f <- "s:/AB_data_v2018/data/raw/veghf/site_all/20180706_All_Sites.sqlite"

f <- "d:/abmi/AB_data_v2019/data/raw/veghf/site_all/20190318_SummaryTables_1ha_TerrestrialSites_Veg61_vHF_LandFacets_SurveyYear_2003_2018.sqlite"
f <- "d:/abmi/AB_data_v2019/data/raw/veghf/site_all/20190129_SummaryTables_CAMARU_2017_2018_Veg61_vHFSPOT2017.sqlite"


f <- "s:/GC_eric/FromEric/Sites_summaries/Round2020/20200602_SC_Sites_visit2019_V61HFI2018v1_SummaryTables_Round_2020_batch04.sqlite"

u <- "s:/GC_eric/FromEric/Sites_summaries/Round2020/"
f <- "20200224_CAMARU_SurveyYear_2019_Buffers_facet_batch01"
f <- "20200224_SC_Sites_SummaryTables_Round_2020.sqlite"
f <- "20200429_SC_Sites_SummaryTables_Round_2020_batch01.sqlite"

f <- "20200429_SC_Sites_SummaryTables_Round_2020_batch02.sqlite"
f <- "20200506_SC_Sites_SummaryTables_Round_2020_batch03.sqlite"
f <- "20200602_SC_Sites_visit2019_V61HFI2018v1_SummaryTables_Round_2020_batch04.sqlite"
f <- "20200702_Bird_Sites.sqlite"
db <- dbConnect(RSQLite::SQLite(), paste0(u,f))
dbListTables(db)
dbDisconnect(db)


db <- dbConnect(RSQLite::SQLite(), f)
tt <- dbListTables(db)
tt
d0 <- dbReadTable(db, "Summary_1ha")

d1 <- dbReadTable(db, "All_Sites_1ha")
d2 <- dbReadTable(db, "All_Sites_SiteCentre")
d3 <- dbReadTable(db, "All_Sites_buffer564m")
d1 <- make_char2fact(d1)
d2 <- make_char2fact(d2)
d3 <- make_char2fact(d3)

## climate 2003-2016
d4 <- read.csv("d:/abmi/AB_data_v2017/data/raw/veghf/site_all/siteCenter_climate.csv",
    stringsAsFactors = TRUE)


dbDisconnect(db)


## 2003-2017
f <- "s:/AB_data_v2018/data/analysis/site/veg-hf_SiteCenter_v6verified.Rdata"
#f <- "s:/AB_data_v2018/data/analysis/site/veg-hf_allSites_v6hfi2016.Rdata"
load(f)

ls()
table(dd_point$SampleYear)




## read in data
if (endsWith(tolower(FILE), ".csv")) {
  cat("Reading CSV file:\n", FILE, "... ")
  d <- read.csv(FILE)
} else {
  cat("Connecting to SQLite database:\n", FILE)
  db <- dbConnect(RSQLite::SQLite(), FILE)
  cat("\n\nFound the following tables:\n")
  cat(paste(dbListTables(db), collapse="\n"))
  cat("\n\nLoading table:\n", TABLE)
  d <- dbReadTable(db, TABLE)
  cat("\n\nDisconnecting ... ")
  dbDisconnect(db)
  d <- make_char2fact(d)
}

## take a subset if needed
if (!is.null(SUB_COL)) {
  cat("OK\n\nTaking subset ... ")
  d <- d[d[[SUB_COL]] %in% SUB_VAL,]
}

if (AREA_COL != "Shape_Area") {
  d[["Shape_Area"]] <- d[[AREA_COL]]
  d[[AREA_COL]] <- NULL
}
cat("OK\n\n")



