source("~/repos/abmianalytics/veghf/veghf-setup.R")

### ABMI on+off grid sites --------------------------------------------------

## ABMI sites (on+off) cetre 1 ha
f1ha <- file.path(ROOT, VER, "data/veghf", "Center1haFixFire.csv")
d1ha <- read.csv(f1ha)
d1ha$Site_YEAR <- with(d1ha, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
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


## clim/topo/nsr table
## fix age0
## merge
## save


if (SAVE)
    save(dd1ha, dd150m, dd1km, climSite, climPoint,
        dd1ha_2015, dd150mCenter_2015, dd1kmCenter_2015,
        dd150mPT_2015, dd1kmPT_2015, climCenter_2015, climPT_2015,
        file=file.path(ROOT, VER, "out/abmi_onoff", 
        "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015.Rdata"))





