source("~/repos/abmianalytics/veghf/veghf-setup.csv")

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
gis <- read.csv(file.path("y:/Oracle_access_2015", "data", "sitemetadata.csv"))
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
clim2$Label2 <- with(clim2, paste0("T_", OnOffGrid, "_", DataProvider, 
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

if (SAVE)
    save(dd1ha, dd150m, dd1km, climSite, climPoint,
        file=file.path(ROOT, VER, "out/abmi_onoff", 
        "veg-hf-clim-reg_abmi-onoff_fix-fire.Rdata"))

## fix age 0 in saved files -----------------------------

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

## dd1ha, dd150m, dd1km, climSite, climPoint
load(file.path(ROOT, VER, "out/abmi_onoff", 
    "veg-hf-clim-reg_abmi-onoff_fix-fire.Rdata"))
#dd1hav <- fill_in_0ages(dd1ha, climSite$NSRNAME)
#round(data.frame(w0=100*colSums(dd1ha$veg_current)/sum(dd1ha$veg_current), 
#    wo0=100*colSums(dd1hav$veg_current)/sum(dd1hav$veg_current)), 4)
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

save(dd1ha, dd150m, dd1km, climSite, climPoint,
    file=file.path(ROOT, VER, "out/abmi_onoff", 
    "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0.Rdata"))

