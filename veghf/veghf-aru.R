source("~/repos/abmianalytics/veghf/veghf-setup.R")

### ARU CONI

## CONI ARU sites, 150 m radius buffer
f150m <- file.path(ROOT, VER, "data/veghf/aru", "CONImodel_ARUs_Buf150m_VegHFV5.csv")
d150m <- read.csv(f150m)
head(d150m)
dd150m <- make_vegHF_wide(d150m, col.label = "POINT_ID", 
    col.year="Year_", col.HFyear="CutYear")
dd150m$scale <- "150 m radius circle around bird points"

## CONI ARU sites, 1 km^2 buffer
f1km <- file.path(ROOT, VER, "data/veghf/aru", "CONImodel_ARUs_Buf564m_VegHFV5.csv")
d1km <- read.csv(f1km)
head(d1km)
dd1km <- make_vegHF_wide(d1km, col.label = "POINT_ID", 
    col.year="Year_", col.HFyear="CutYear")
dd1km$scale <- "564 m radius circle around bird points"

stopifnot(all(rownames(dd150m[[1]]) == rownames(dd1km[[1]])))
clim <- read.csv("e:/peter/AB_data_v2016/data/veghf/aru/CONImodel_ARUs_climate.csv")
clim <- nonDuplicated(clim, POINT_ID, TRUE)
clim <- clim[rownames(dd150m[[1]]),]
topo <- read.csv("e:/peter/AB_data_v2016/data/veghf/aru/CONImodel_ARUs_topo.csv")
topo <- nonDuplicated(topo, POINT_ID, TRUE)
topo <- topo[rownames(dd150m[[1]]),]

pts <- read.csv("e:/peter/AB_data_v2016/data/veghf/aru/CONImodel_ARUs_VegHFV5.csv")
pts <- nonDuplicated(pts, POINT_ID, TRUE)
pts <- pts[rownames(dd150m[[1]]),]

points <- data.frame(pts[,c("POINT_ID", "NRNAME","NSRNAME")],
    clim[,c("ProjectID", "Cluster", "SITE", "STATION", 
    "ARU", "UTM_Zone", "EASTING", "NORTHING", "Year_", "DEPLOY_DATE", 
    "RETRIEVE_DATE", "AHM", "bFFP", "CMD", "eFFP", "EMT", 
    "Eref", "FFP", "MAP", "MAT", "MCMT", "MSP", "MWMT", "NFFD", "PAS", 
    "PPT01", "PPT02", "PPT03", "PPT04", "PPT05", "PPT06", "PPT07", 
    "PPT08", "PPT09", "PPT10", "PPT11", "PPT12", "SHM", "Tave01", 
    "Tave02", "Tave03", "Tave04", "Tave05", "Tave06", "Tave07", "Tave08", 
    "Tave09", "Tave10", "Tave11", "Tave12", "TD", "TIMEFIREINV", 
    "Tmax01", "Tmax02", "Tmax03", "Tmax04", "Tmax05", "Tmax06", "Tmax07", 
    "Tmax08", "Tmax09", "Tmax10", "Tmax11", "Tmax12", "Tmin01", "Tmin02", 
    "Tmin03", "Tmin04", "Tmin05", "Tmin06", "Tmin07", "Tmin08", "Tmin09", 
    "Tmin10", "Tmin11", "Tmin12")],
    topo[,c("aggregates", "cti", "dunes", "geol_bedcode", 
    "geol_surf", "slope", "slpasp", "spi", "tpi2km", "tri", "tri_1", 
    "vrm11x11", "vrm5x5")])

## fix age 0 in saved files -----------------------------

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

sum(dd150m[[1]][,Target0])
dd150m <- fill_in_0ages(dd150m, points$NSRNAME)
sum(dd150m[[1]][,Target0])

sum(dd1km[[1]][,Target0])
dd1km <- fill_in_0ages(dd1km, points$NSRNAME)
sum(dd1km[[1]][,Target0])

if (SAVE)
    save(dd150m, dd1km, points,
        file=file.path(ROOT, VER, "out/aru", 
        "veg-hf-clim-reg_aru-coni_fix-fire_fix-age0.Rdata"))

