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
f4 <- "d:/abmi/AB_data_v2020/data/analysis/site/veg-hf_SITE1HA-2019_Veg61-vHF.Rdata"
e4 <- new.env()
load(f4, envir=e4)

## 2017-2018 updates
f5 <- "d:/abmi/AB_data_v2020/data/analysis/site/veg-hf_SITES-2017-2018_Veg61-vHF.Rdata"
e5 <- new.env()
load(f5, envir=e5)

xy <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")


table(e1$xx$survey_year)

compare_sets(rownames(e1$dd_1ha$veg_current), rownames(e5$d_wide_1ha$veg_current))

## site IDs and regions
SITES <- e1$xx[,1:6]
s19 <- nonDuplicated(e4$d_long, UID_old, TRUE)
s19 <- s19[,c("UID_old", "ABMI_ID_WithB", "survey_year", "NSRNAME", "NRNAME", "LUF_NAME")]
colnames(s19) <- colnames(SITES)
SITES <- rbind(SITES, s19)
table(SITES$survey_year)
table(e5$clim$survey_year)
s78 <- e5$clim[,c("survey_year","NSRNAME","NRNAME","LUF_NAME")]
tmp <- strsplit(rownames(s78), "_")
s78$UID <- rownames(s78)
s78$site_id <- sapply(tmp, "[[", 1)
s78 <- s78[,colnames(SITES)]
SITES <-SITES[SITES$survey_year %notin% c(2017, 2018),]
SITES <- rbind(SITES, s78)

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

cl2 <- e5$clim[,colnames(cl)]
compare_sets(rownames(cl), rownames(cl2))
cl <- rbind(cl, cl2)

SITES <- data.frame(SITES, cl[match(rownames(SITES), rownames(cl)),])

summary(SITES)

if (FALSE) { # not necessary, has been extracted properly
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
}

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

## qha
dd_qha <- e1$dd_qha
for (i in 1:4)
    dd_qha[[i]] <- rbind(
        e1$dd_qha[[i]],
        e4$d_wide[[i]][,colnames(e1$dd_qha[[i]])])
dd_qha$scale <- NULL
dd_qha$sample_year <- c(e1$dd_qha[[5]], e4$d_wide[[5]])

table(dd_qha$sample_year)

# drop 2017-2018 years
ss <- dd_qha$sample_year %notin% c(2017, 2018)
for (i in 1:4)
    dd_qha[[i]] <- rbind(dd_qha[[i]][ss,],
        e5$d_wide_qha[[i]][,colnames(dd_qha[[i]])])
dd_qha$sample_year <- c(dd_qha$sample_year, e5$d_wide_qha$sample_year)

ss <- dd_1ha$sample_year %notin% c(2017, 2018)
for (i in 1:4)
    dd_1ha[[i]] <- rbind(dd_1ha[[i]][ss,],
        e5$d_wide_1ha[[i]][,colnames(dd_1ha[[i]])])
dd_1ha$sample_year <- c(dd_1ha$sample_year, e5$d_wide_1ha$sample_year)

## 1km
dd_564m <- e3$dd_564m
for (i in 1:4)
    dd_564m[[i]] <- rbind(
        e3$dd_564m[[i]],
        e4$d_wide_1km[[i]][,colnames(e3$dd_564m[[i]])])
dd_564m$scale <- NULL
dd_564m$sample_year <- c(e3$dd_564m[[5]], e4$d_wide_1km[[5]])

for (i in 1:4)
    rownames(dd_564m[[i]]) <- gsub("Confidential_", "Confidential-", rownames(dd_564m[[i]]))

ss <- dd_564m$sample_year %notin% c(2017, 2018)
for (i in 1:4)
    dd_564m[[i]] <- rbind(dd_564m[[i]][ss,],
        e5$d_wide_1km[[i]][,colnames(dd_564m[[i]])])
dd_564m$sample_year <- c(dd_564m$sample_year, e5$d_wide_1km$sample_year)

compare_sets(rownames(SITES), rownames(dd_564m[[1]]))
data.frame(x=setdiff(rownames(dd_564m[[1]]), rownames(SITES)))

if (FALSE) { # no need, fixed now
## use 2018 w2w results
library(cure4insect)
load("s:/AB_data_v2020/data/analysis/veghf/w2w_veghf2018_grid1sqkm.RData")
load_common_data()
TB <- get_id_table()
XY <- get_id_locations()
XY <- spTransform(XY, proj4string(rr))
XY <- coordinates(XY)

cr <- SITES[!(rownames(SITES) %in% rownames(dd_564m[[1]])),]
coordinates(cr) <- ~ X + Y
proj4string(cr) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
cr <- spTransform(cr, proj4string(rr))
cr <- coordinates(cr)

id <- character(nrow(cr))
names(id) <- rownames(cr)

for (i in 1:nrow(cr)) {
    Dis <- sqrt((XY[,1] - cr[i,1])^2 + (XY[,2] - cr[i,2])^2)
    id[[i]] <- rownames(XY)[which.min(Dis)]
}
compare_sets(rownames(dd18[[1]]), id)

for (i in 1:4) {
    tmp <- dd18[[i]][id,]
    rownames(tmp) <- names(id)
    dd_564m[[i]] <- rbind(
        dd_564m[[i]],
        tmp)
}
dd_564m$sample_year <- c(dd_564m[[5]], rep(2018, length(id)))
names(dd_564m$sample_year) <- rownames(dd_564m[[1]])

for (i in 1:4)
    dd_564m[[i]] <- dd_564m[[i]][rownames(SITES),]
dd_564m[[5]] <- dd_564m[[5]][rownames(SITES)]

compare_sets(rownames(SITES), rownames(dd_564m[[1]]))

SITES$km_source <- ifelse(rownames(SITES) %in% names(id), "w2w_HFI2018", "exact")
}

## some messed up site names
rn <- matrix(c(
  c("OG-ABMI-448-2_2013", "OG-ABMI-448-1_2013",
  "OG-ABMI-1057-1_2014", "OG-ABMI-1057-11_2014",
  "OG-ABMI-468-11_2014", "OG-ABMI-468-1_2014",
  "OG-ABMI-531-11_2014", "OG-ABMI-531-1_2014",
  "OG-ABMI-543-11_2014", "OG-ABMI-543-1_2014",
  "OG-ABMI-562-11_2014", "OG-ABMI-562-1_2014",
  "OG-ABMI-590-11_2014", "OG-ABMI-590-1_2014",
  "OG-ABMI-653-11_2014", "OG-ABMI-653-1_2014",
  "OG-ABMI-660-11_2014", "OG-ABMI-660-1_2014",
  "OG-ABMI-687-11_2014", "OG-ABMI-687-1_2014",
  "OG-ABMI-720-11_2014", "OG-ABMI-720-1_2014",
  "OG-ABMI-783-11_2014", "OG-ABMI-783-1_2014",
  "OG-ABMI-927-11_2014", "OG-ABMI-927-1_2014",
  "OG-ABMI-910-11_2014", "OG-ABMI-910-3_2014",
  "OG-ABMI-943-11_2014", "OG-ABMI-943-2_2014")
  ), ncol=2, byrow=TRUE)
colnames(rn) <- c("Species", "GIS")

names(dd_1ha$sample_year) <- rownames(dd_1ha$veg_current)
names(dd_qha$sample_year) <- rownames(dd_qha$veg_current)
names(dd_564m$sample_year) <- rownames(dd_564m$veg_current)

for (i in 1:4) {
  dd_1ha[[i]] <- dd_1ha[[i]][rownames(SITES),]
  dd_564m[[i]] <- dd_564m[[i]][rownames(SITES),]
}
dd_1ha[[5]] <- dd_1ha[[5]][rownames(SITES)]
dd_564m[[5]] <- dd_564m[[5]][rownames(SITES)]


SITES$WrongUID <- SITES$UID
for (i in 1:nrow(rn))
    SITES$UID[SITES$UID == rn[i,"GIS"]] <- rn[i,"Species"]
rownames(SITES) <- SITES$UID

for (i in 1:4) {
  rownames(dd_1ha[[i]]) <- SITES$UID
  rownames(dd_564m[[i]]) <- SITES$UID
}
names(dd_1ha[[5]]) <- SITES$UID
names(dd_564m[[5]]) <- SITES$UID


clim_1ha <- SITES[rownames(dd_1ha[[1]]),]
all(rownames(clim_1ha) == rownames(dd_1ha[[1]]))

clim_qha <- rbind(
    data.frame(SITES, quadrant="NE"),
    data.frame(SITES, quadrant="NW"),
    data.frame(SITES, quadrant="SE"),
    data.frame(SITES, quadrant="SW"))
rownames(clim_qha) <- paste0(clim_qha$WrongUID, "_", clim_qha$quadrant)
compare_sets(rownames(clim_qha), rownames(dd_qha[[1]]))
clim_qha <- clim_qha[rownames(dd_qha[[1]]),]

rownames(clim_qha) <- paste0(clim_qha$UID, "_", clim_qha$quadrant)
for (i in 1:4) {
  rownames(dd_qha[[i]]) <- rownames(clim_qha)
}
names(dd_qha[[5]]) <-rownames(clim_qha)
all(rownames(clim_qha) == rownames(dd_qha[[1]]))


save(
#    SITES,
    clim_1ha, clim_qha,
    dd_1ha, dd_qha, dd_564m,
    file="s:/AB_data_v2020/data/analysis/site/veg-hf_SITE-all-years-combined.RData")
save(
#    SITES,
    clim_1ha, clim_qha,
    dd_1ha, dd_qha, dd_564m,
    file="d:/abmi/AB_data_v2020/data/analysis/site/veg-hf_SITE-all-years-combined.RData")


