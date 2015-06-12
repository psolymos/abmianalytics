##% Processing backfilled veg + HF (cutblock ages incorporated)
##% P Solymos
##% April 28, 2015

## root directory
ROOT <- "c:/p"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2015"
## current year
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900

library(mefa4)
source("~/repos/abmianalytics/veghf/veghf_functions.R")
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

#hftypes <- read.csv(file.path(ROOT, VER, "lookup/HFtype_lookup_20150428.csv"))
#hfgroups <- read.csv(file.path(ROOT, VER, "lookup/HFclassification_20150428.csv"))
hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type.csv")
hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY

### ABMI on+off grid sites

#### Vegetation and HF processing

## ABMI sites (on+off) cetre 1 ha
f1ha <- file.path(ROOT, VER, "data/veghf", "Center1ha.csv")
d1ha <- read.csv(f1ha)
d1ha$Site_YEAR <- with(d1ha, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d1ha)
dd1ha <- make_vegHF_wide(d1ha, col.label = "Site_YEAR", col.year="survey_year")
dd1ha$scale <- "1 ha square around site centre"

## ABMI sites (on+off) 9 bird points / site, 150 m radius buffer
f150m <- file.path(ROOT, VER, "data/veghf", "Bird150m.csv")
d150m <- read.csv(f150m)
d150m$Site_YEAR_bird <- with(d150m, 
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d150m)
dd150m <- make_vegHF_wide(d150m, col.label = "Site_YEAR_bird", col.year="survey_year")
dd150m$scale <- "150 m radius circle around bird points"

## ABMI sites (on+off) 9 bird points / site, 1 km^2 buffer
f1km <- file.path(ROOT, VER, "data/veghf", "Bird564m.csv")
d1km <- read.csv(f1km)
d1km$Site_YEAR_bird <- with(d1km, 
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d1km)
dd1km <- make_vegHF_wide(d1km, col.label = "Site_YEAR_bird", col.year="survey_year")
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

save(dd1ha, dd150m, dd1km, climSite, climPoint,
    file=file.path(ROOT, VER, "out/abmi_onoff", "veg-hf-clim-reg_abmi-onoff.Rdata"))


### Snow transects

## 1 km length (250 m buffer) mammal transect (inter level)
fmi <- file.path(ROOT, VER, "data/veghf", "InterLevel.csv")
dmi <- read.csv(fmi)
dmi$Site_YEAR_tr <- with(dmi, interaction(ABMISite, survey_year, interLevel, sep="_", drop=TRUE))
ddmi <- make_vegHF_wide(dmi, col.label = "Site_YEAR_tr", col.year="survey_year")
ddmi$scale <- "inter level mammal transects"

## 9-10 km length (250 m buffer) mammal transect (full transect level)
fmt <- file.path(ROOT, VER, "data/veghf", "TransectLevel.csv")
dmt <- read.csv(fmt)
dmt$Site_YEAR <- with(dmt, interaction(ABMISite, survey_year, sep="_", drop=TRUE))
## strange site issue: "394-2005_2005" --> "394-2005_2006"
levels(dmt$Site_YEAR)[levels(dmt$Site_YEAR)=="394-2005_2005"] <- "394-2005_2006"
ddmt <- make_vegHF_wide(dmt, col.label = "Site_YEAR", col.year="survey_year")
ddmt$scale <- "full transect level mammal transects"

## Transect segment labels
seg <- nonDuplicated(dmi, Site_YEAR_tr, TRUE)
tmp <- strsplit(rownames(seg), "_")
seg$Site <- as.factor(sapply(tmp, "[[", 1))
seg$Inter <- as.integer(sapply(tmp, "[[", 3))
seg$Site_Inter <- with(seg, interaction(Site, Inter, sep="_", drop=TRUE))
seg$OnOffGrid <- as.factor(ifelse(substr(sapply(tmp, "[[", 1), 1, 2) == "OG", "OG", "IG"))

## mammal stuff
clim3 <- read.csv(file.path(ROOT, VER, "data/climate", 
    "mamTrack_interLevel_latLong_climate_naturalReg_V2.csv"))
colnames(clim3)[colnames(clim3) == "Eref"] <- "PET"
colnames(clim3)[colnames(clim3) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
clim3$Site_Inter <- with(clim3, interaction(ABMISite, interLevel, sep="_", drop=TRUE))
clim3$Year <- seg$survey_year[match(clim3$Site_Inter, seg$Site_Inter)]
clim3$Site_Year_Inter <- with(clim3, 
    interaction(ABMISite, Year, interLevel, sep="_", drop=TRUE))
clim3$Site_Inter[is.na(clim3$Site_Year_Inter)]
clim3 <- droplevels(clim3[!is.na(clim3$Site_Year_Inter),])
clim3$Site_Year <- as.factor(paste(clim3$ABMISite, clim3$Year, sep="_"))
rownames(clim3) <- clim3$Site_Year_Inter

compare.sets(rownames(clim3), rownames(ddmi$veg_current))

seg2 <- nonDuplicated(clim3, Site_Year, TRUE)
seg2 <- seg2[,c("ABMISite","interLevel","Site_Inter","Year","Site_Year_Inter","Site_Year")]
clim4 <- read.csv(file.path(ROOT, VER, "data/climate", 
    "mamTrack_latLong_climate_naturalReg_V2.csv"))
colnames(clim4)[colnames(clim4) == "Eref"] <- "PET"
colnames(clim4)[colnames(clim4) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
clim4$Year <- seg2$Year[match(clim4$ABMISite, seg2$ABMISite)]
rownames(clim4) <- paste(clim4$ABMISite, clim4$Year, sep="_")

compare.sets(clim4$ABMISite, seg2$ABMISite)
compare.sets(rownames(clim4), rownames(ddmt$veg_current))
setdiff(rownames(clim4), rownames(ddmt$veg_current))
setdiff(rownames(ddmt$veg_current), rownames(clim4))

inter <- read.csv(file.path(ROOT, VER, "out","species",
    "OUT_Mammals_Species_InterSegment_2015-06-01.csv"))
tran <- read.csv(file.path(ROOT, VER, "out","species",
    "OUT_Mammals_Species_Transect-Binomial-Length-DSS_2015-06-01.csv"))
climInter <- clim3[rownames(ddmi$veg_current),]
climTr <- clim4[rownames(ddmt$veg_current),]
all(rownames(climInter) == rownames(ddmi[[1]]))
all(rownames(climTr) == rownames(ddmt[[1]]))
rm(clim3, clim4)

tran$SiteFunny <- tran$Site
levels(tran$SiteFunny)[levels(tran$SiteFunny) == "1"] <- "001"
levels(tran$SiteFunny)[levels(tran$SiteFunny) == "2"] <- "002"
levels(tran$SiteFunny)[levels(tran$SiteFunny) == "3"] <- "003"
tran$YearFunny <- paste0("-", as.character(tran$Year))
tran$YearFunny[grep("-ILM-", as.character(tran$Site))] <- ""
tran$SiteYear <- paste0(tran$SiteFunny, tran$YearFunny)
compare.sets(tran$SiteYear, climTr$ABMISite)
## Transect level is all good
setdiff(tran$SiteYear, climTr$ABMISite)
setdiff(climTr$ABMISite, tran$SiteYear)
climTr$label_tr <- as.character(tran$label_tr)[match(climTr$ABMISite,
    tran$SiteYear)]
climTr$label_tr[is.na(climTr$label_tr)] <- as.character(climTr$ABMISite)[is.na(climTr$label_tr)]
climTr$label_tr <- as.factor(climTr$label_tr)


inter$SiteFunny <- inter$Site
levels(inter$SiteFunny)[levels(inter$SiteFunny) == "1"] <- "001"
levels(inter$SiteFunny)[levels(inter$SiteFunny) == "2"] <- "002"
levels(inter$SiteFunny)[levels(inter$SiteFunny) == "3"] <- "003"
inter$YearFunny <- paste0("-", as.character(inter$Year))
inter$YearFunny[grep("-ILM-", as.character(inter$Site))] <- ""
inter$SiteYear <- paste0(inter$SiteFunny, inter$YearFunny)
compare.sets(inter$SiteYear, climInter$ABMISite)
inter$Site_Year_Inter <- with(inter, interaction(SiteYear, Year, InterSegID, 
    sep="_", drop=TRUE))
compare.sets(inter$Site_Year_Inter, climInter$Site_Year_Inter)
setdiff(inter$Site_Year_Inter, climInter$Site_Year_Inter)
setdiff(climInter$Site_Year_Inter, inter$Site_Year_Inter)

climInter$label_int <- as.character(inter$label_int)[match(climInter$Site_Year_Inter,
    inter$Site_Year_Inter)]
climInter$label_int[is.na(climInter$label_int)] <- 
    as.character(climInter$Site_Year_Inter)[is.na(climInter$label_int)]

climInter$label_tr <- as.character(inter$label_tr)[match(climInter$Site_Year_Inter,
    inter$Site_Year_Inter)]
climInter$label_tr[is.na(climInter$label_tr)] <- 
    as.character(climInter$Site_Year)[is.na(climInter$label_tr)]
climInter$label_int <- as.factor(climInter$label_int)
climInter$label_tr <- as.factor(climInter$label_tr)

rownames(climInter) <- climInter$label_int
rownames(ddmi[[1]]) <- rownames(ddmi[[2]]) <- rownames(climInter)
rownames(ddmi[[3]]) <- rownames(ddmi[[4]]) <- rownames(climInter)
rownames(climTr) <- climTr$label_tr
rownames(ddmt[[1]]) <- rownames(ddmt[[2]]) <- rownames(climTr)
rownames(ddmt[[3]]) <- rownames(ddmt[[4]]) <- rownames(climTr)
all(rownames(climInter) == rownames(ddmi[[1]]))
all(rownames(climTr) == rownames(ddmt[[1]]))

save(ddmi, ddmt, climInter, climTr,
    file=file.path(ROOT, VER, "out/abmi_onoff", "veg-hf-clim-reg_mammals-onoff.Rdata"))


### BAM+BBS bird points, 150 m radius buffer

## processing csv files
fl <- list.files(file.path(ROOT, VER, "data", "veghf", "bammbbs150m"))
tmplist <- list()
for (fn in fl) {
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "veghf", "bammbbs150m", fn)
    d <- read.csv(f)
    dd <- make_vegHF_wide(d, col.label="PKEY", col.year="YEAR_", sparse=TRUE)
    tmplist[[fn]] <- dd
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
}

## assembling return object
dd150m_bambbs <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year=NA,
    scale = "150 m radius circle around bird points")

### BAM+BBS bird points, 1 km^2 buffer

## processing csv files
fl <- list.files(file.path(ROOT, VER, "data", "veghf", "bammbbs564m"))
tmplist <- list()
for (fn in fl) {
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "veghf", "bammbbs564m", fn)
    d <- read.csv(f)
    dd <- make_vegHF_wide(d, col.label="PKEY", col.year="YEAR_", sparse=TRUE)
    tmplist[[fn]] <- dd
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
}

## assembling return object
dd1km_bambbs <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year=NA,
    scale = "564 m radius circle around bird points")

save(dd150m_bambbs, dd1km_bambbs,
    file=file.path(ROOT, VER, "out/bambbs", "veg-hf_bambbs.Rdata"))

### 1K grid

## Sample year is current year, so that forest ages are relative to present
## and not relative to HF or veg inventory year.

fl <- list.files(file.path(ROOT, VER, "data", "kgrid", "tiles"))

## test feature types
if (FALSE) {
NEW <- character(0)
for (fn in fl) {
    cat("checking", which(fl == fn), "/", length(fl));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
    d <- read.csv(f)
    diff <- setdiff(levels(d$FEATURE_TY), rownames(hflt))
    if (length(diff))
        NEW <- union(NEW, diff)
    cat("\t", length(NEW), "new types found\n")
}

## tracking the strange area mismatch
natrack <- list()
for (fn in fl) {
    cat("checking", which(fl == fn), "/", length(fl), "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
    d <- read.csv(f)
    dd <- make_vegHF_wide(d, col.label="Row_Col", col.year=NULL, wide=FALSE)
    tmp <- colSums(is.na(dd[,c("VEGAGEclass",
        "VEGHFAGEclass","SOILclass","SOILHFclass")]))
    natrack[[fn]] <- tmp
    if (tmp[1] > 0)
        break
}

}

## processing csv files in batches of 50

Start <- c(1, 51, 101, 151, 201, 251, 301, 351, 401, 451,
    501, 551, 601, 651, 701, 751, 802)
tmplist <- list()

for (s in 1:(length(Start)-1)) {

    gc()
    fn <- fl[Start[s]]
    cat("\n\n------------- batch", s, "----------------\n")
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
    d <- read.csv(f)
    dd <- make_vegHF_wide(d, col.label="Row_Col", col.year=NULL, sparse=TRUE)
    veg_current <- dd$veg_current
    veg_reference <- dd$veg_reference
    soil_current <- dd$soil_current
    soil_reference <- dd$soil_reference
    sample_year <- dd$sample_year[1]

#lapply(dd[1:4], sum)

    for (i in (Start[s]+1):(Start[s+1]-1)) {

        fn <- fl[i]
        cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
        f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
        d <- read.csv(f)
        dd <- make_vegHF_wide(d, col.label="Row_Col", col.year=NULL, sparse=TRUE)
        veg_current <- bind_fun2(veg_current, dd$veg_current)
        veg_reference <- bind_fun2(veg_reference, dd$veg_reference)
        soil_current <- bind_fun2(soil_current, dd$soil_current)
        soil_reference <- bind_fun2(soil_reference, dd$soil_reference)

    }
    tmplist[[s]] <- list(
        veg_current = veg_current,
        veg_reference = veg_reference,
        soil_current = soil_current,
        soil_reference = soil_reference,
        sample_year = sample_year,
        scale = "1 km x 1 km prediction grid cells")
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
}

## assembling return object
dd1km_pred <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year = tmplist[[1]]$sample_year,
    scale = "1 km x 1 km prediction grid cells")

kgrid <- read.csv(
    file.path(ROOT, VER, "data", "kgrid", 
    "Grid1km_template_final_clippedBy_ABBound_with_atts_to_Peter.csv"))
rownames(kgrid) <- kgrid$Row_Col

compare.sets(rownames(dd1km_pred$veg_current), rownames(kgrid))

## NSR x LUF regions used as prediction regions in sector effects etc.
kgrid$nsr_luf <- with(kgrid, paste(as.integer(NSRNAME), as.integer(LUF_NAME), sep="_"))
colnames(kgrid)[colnames(kgrid) == "col"] <- "Col"
colnames(kgrid)[colnames(kgrid) == "Eref"] <- "PET"
colnames(kgrid)[colnames(kgrid) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"

## 10 x 10 km grid
kgrid$Row10 <- 1 + kgrid$Row %/% 10
kgrid$Col10 <- 1 + kgrid$Col %/% 10
kgrid$Row10_Col10 <- interaction(kgrid$Row10, kgrid$Col10, sep="_", drop=TRUE)

## random pick from 10K grid
tmp <- as.integer(kgrid$Row10_Col10)
kgrid$Rnd10 <- integer(length(tmp))
set.seed(1234)
for (i in seq_len(max(tmp))) {
    lg <- tmp == i
    kgrid$Rnd10[lg] <- sample.int(sum(lg))
}

dd1km_pred$veg_current <- dd1km_pred$veg_current[rownames(kgrid),]
dd1km_pred$veg_reference <- dd1km_pred$veg_reference[rownames(kgrid),]
dd1km_pred$soil_current <- dd1km_pred$soil_current[rownames(kgrid),]
dd1km_pred$soil_reference <- dd1km_pred$soil_reference[rownames(kgrid),]

## check area diff
range(sapply(dd1km_pred[1:4], sum) / 10^6)

## proportion of water -- for mapping purposes
kgrid$pWater <- dd1km_pred$veg_current[,"Water"] / 10^6

## veg based area < soil based area, thus using the max
kgrid$Area_km2 <- rowSums(dd1km_pred$soil_reference) / 10^6

## UTM projection for fake maps
library(raster)
library(sp)
library(rgdal)
XYlatlon <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(XYlatlon) <- ~ POINT_X + POINT_Y
proj4string(XYlatlon) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
XY <- as.data.frame(spTransform(XYlatlon, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))
kgrid$X <- XY$POINT_X
kgrid$Y <- XY$POINT_Y

kgrid$NEAR_DIST <- NULL

save(dd1km_pred, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid.Rdata"))
save(kgrid,
    file=file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

#tab_veg <- data.frame(Label=colnames(dd1km_pred[[2]]),
#    prop_cr=colSums(dd1km_pred[[2]]) / sum(dd1km_pred[[2]]))
#write.csv(tab_veg, file=file.path(ROOT, VER, "out", "tab_veg.csv"))

ltveg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg.csv")
ltsoil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil.csv")

tveg <- data.frame(VEGHFAGE=colnames(dd1km_pred$veg_current))
tveg$HF <- hfgroups$HF_GROUP[match(tveg$VEGHFAGE, hfgroups$HF_GROUP)]
tveg$HF[substr(as.character(tveg$VEGHFAGE), 1, 2) == "CC"] <- "CutBlocks"
tveg$VEGAGE <- ltveg$VEGAGE[match(tveg$VEGHFAGE, ltveg$VEGAGE)]
tveg$VEGAGE <- as.character(tveg$VEGAGE)
tveg$VEGHFAGE <- as.character(tveg$VEGHFAGE)
tveg$VEGAGE[substr(tveg$VEGHFAGE, 1, 2) == "CC"] <- 
    substr(tveg$VEGHFAGE, 3, nchar(tveg$VEGHFAGE))[substr(tveg$VEGHFAGE, 1, 2) == "CC"]
tveg <- data.frame(tveg, 
    ltveg[match(tveg$VEGAGE, ltveg$VEGAGE),-1],
    hfgroups[match(tveg$HF, hfgroups$HF_GROUP),-1])
tveg$VEGAGE <- as.factor(tveg$VEGAGE)
rownames(tveg) <- tveg$VEGHFAGE
colnames(tveg)[colnames(tveg)=="Type.1"] <- "HFtype"

tsoil <- data.frame(SOILHF=colnames(dd1km_pred$soil_current))
tsoil$HF <- hfgroups$HF_GROUP[match(tsoil$SOILHF, hfgroups$HF_GROUP)]
tsoil$SOIL <- ltsoil$SOILclass[match(tsoil$SOILHF, ltsoil$SOILclass)]
tsoil <- data.frame(tsoil, 
    ltsoil[match(tsoil$SOILHF, ltsoil$SOILclass),-1],
    hfgroups[match(tsoil$SOILHF, hfgroups$HF_GROUP),-1])
rownames(tsoil) <- tsoil$SOILHF

write.csv(tveg, file="~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
write.csv(tsoil, file="~/repos/abmianalytics/lookup/lookup-soil-hf.csv")


### Transition for 1K grid




## QS for transition -------------------------------------------------------

source("c:/Dropbox/abmi/intactness/dataproc/data_proc_common_2014.R")

fl <- list.files(file.path(ROOT, VER, "data/veghf/qs/csv"))

cc <- c("LinkID","VEGAGEclass","VEGHFAGEclass","SOILclass","SOILHFclass","Shape_Area")

Start <- c(0:79*10+1, 802)

d <- read.csv(file.path(ROOT, VER, "data/veghf/qs/csv", fl[1]))
colnames(d)[colnames(d) == "year"] <- "HF_Year"
dd <- make_vegHF_wide2(d, col.label="LinkID", col.year=NULL, wide=FALSE)
ddd0 <- dd[character(0),cc]
xddd0 <- ddd0

for (s in 1:(length(Start)-1)) {
    cat("----------------------\nStarting block", s, "\n")
    for (i in Start[s]:(Start[s+1]-1)) {
        cat(i, "of", length(fl), "-", fl[i], "\t")
        flush.console()
        if (i == Start[s]) {
            d <- read.csv(file.path(ROOT, VER, "data/veghf/qs/csv", fl[i]))
            colnames(d)[colnames(d) == "year"] <- "HF_Year"
            dd <- make_vegHF_wide2(d, col.label="LinkID", col.year=NULL, wide=FALSE)
            dd0 <- dd[,cc]
        } else {
            d <- read.csv(file.path(ROOT, VER, "data/veghf/qs/csv", fl[i]))
            colnames(d)[colnames(d) == "year"] <- "HF_Year"
            dd <- make_vegHF_wide2(d, col.label="LinkID", col.year=NULL, wide=FALSE)
            dd0 <- rbind(dd0, dd[,cc])
        }
        cat("OK", nrow(dd0), "\n")
    }
    ddd0 <- rbind(ddd0, dd0)
    cat("\nFinished block", s, "dim:", nrow(ddd0), "\n")
    if (i %in% c(100, 200, 300, 400, 500, 600, 700, 801)) {
        save(ddd0, file=file.path(ROOT, VER, "data/veghf/qs/long", paste0("Long-part", i, ".Rdata")))
        ddd0 <- xddd0
        gc()
    }
}


load(file.path(ROOT, VER, "R/xy_clim_regions_QSlevel.Rdata"))
lu <- read.csv(file.path(ROOT, VER, "lookup/VEG_HF_interim.csv"))
su <- read.csv(file.path(ROOT, VER, "lookup/SOIL_HF_interim.csv"))

fl3 <- list.files(file.path(ROOT, VER, "data/veghf/qs/long"))

resv <- list()
ress <- list()

for (j in 1:length(fl3)) {
    cat("\n", j);flush.console()

    load(file.path(ROOT, VER, "data/veghf/qs/long", fl3[j]))
    ddd0 <- ddd0[ddd0$Shape_Area >= 0,]
    ddd0$NR <- QS$NR[match(ddd0$LinkID, rownames(QS))]
    ddd0$NSR <- QS$NSR[match(ddd0$LinkID, rownames(QS))]
    ddd0 <- ddd0[rowSums(is.na(ddd0))==0,]
    zz <- ddd0$NR != ""
    cat("\tblanks:", sum(!zz), "\n")
    ddd0 <- ddd0[zz,]
    ddd0$NR <- droplevels(ddd0$NR)
    ddd0$NSR <- droplevels(ddd0$NSR)
#    levels(ddd0$NR)[levels(ddd0$NR) %in% c("Boreal","Canadian Shield")] <- "BorSh"

    ## define region here
    #ddd0$REG <- ddd0$NSR

    ddd0$REG <- ddd0$NR
    levels(ddd0$REG)[levels(ddd0$REG) %in% c("Grassland","Parkland")] <- "South"
    levels(ddd0$REG)[levels(ddd0$REG) %in% c("Boreal","Foothills","Rocky Mountain","Canadian Shield")] <- "North"
    ddd0$REG[ddd0$NSR %in% "Dry Mixedwood"] <- "South"

    levels(ddd0$VEGAGEclass) <- lu$Levels4[match(levels(ddd0$VEGAGEclass), lu$VEGHFAGE)]
    levels(ddd0$VEGHFAGEclass) <- lu$Levels4[match(levels(ddd0$VEGHFAGEclass), lu$VEGHFAGE)]
    levels(ddd0$SOILclass) <- su$Levels3[match(levels(ddd0$SOILclass), su$SOILclass)]
    levels(ddd0$SOILHFclass) <- su$Levels3[match(levels(ddd0$SOILHFclass), su$SOILclass)]

    for (type in c("veg", "soil")) {
#    for (type in c("veg")) {
        if (type=="veg") {
            xt0 <- Xtab(Shape_Area ~ VEGHFAGEclass + VEGAGEclass, ddd0)
            xt <- Xtab(Shape_Area ~ VEGHFAGEclass + VEGAGEclass + REG, ddd0)
        }
        if (type=="soil") {
            xt0 <- Xtab(Shape_Area ~ SOILHFclass + SOILclass, ddd0)
            xt <- Xtab(Shape_Area ~ SOILHFclass + SOILclass + REG, ddd0)
        }

        xtv <- vector("list", nlevels(ddd0$REG)+1)
        names(xtv) <- c("All", levels(ddd0$REG))
        xtv$All <- Melt(xt0)
        xtv$All$gr <- factor("All", levels=names(xtv))
        for (i in names(xtv)[-1]) {
            if (i %in% names(xt)) {
                xtv[[i]] <- Melt(xt[[i]])
                xtv[[i]]$gr <- factor(i, levels=names(xtv))
            } else {
                tmp <- xt0
                tmp[] <- 0
                xtv[[i]] <- Melt(tmp)
                xtv[[i]]$gr <- factor(character(0), levels=names(xtv))
            }
        }
        if (type=="veg")
            resv[[j]] <- do.call(rbind, xtv)
        if (type=="soil")
            ress[[j]] <- do.call(rbind, xtv)
    }
}
resvv <- do.call(rbind, resv)
resss <- do.call(rbind, ress)

save(resvv, resss, file=file.path(ROOT, VER, "R/transitions_NorthSouth.Rdata"))
#save(resvv, resss, file=file.path(ROOT, VER, "R/transitions_NSR.Rdata"))

mat_fun <- function(what="All", type="veg") {
    d <- switch(type, "veg"=resvv, "soil"=resss)
    tr0 <- as.matrix(Xtab(value ~ rows + cols, d, subset=d$gr %in% what))
    i <- rownames(tr0) %in% colnames(tr0)
    tr0 <- tr0[c(rownames(tr0)[i], rownames(tr0)[!i]),]
    #tr0 <- as.matrix(Xtab(Shape_Area ~ cr + rf, dd))
    tr <- round(100*t(t(tr0)/colSums(tr0)), 2)
    #tr[tr<1] <- 0
    c <- round(100*colSums(tr0)/sum(tr0),2)
    r <- c(round(100*rowSums(tr0)/sum(tr0), 2), Total=100)
    tr2 <- cbind(rbind(tr, Total=c), Total=r)
    #write.csv(tr2, file.path(ROOT, VER, "R/transitions.csv"))
    tr2
}

a1 <- mat_fun("All", "veg")
a2 <- mat_fun("South", "soil")
a3 <- mat_fun("North", "veg")

write.csv(a1, file.path(ROOT, VER, "R/transitions_veg_allProvince.csv"))
#write.csv(a1, file.path(ROOT, VER, "R/transitions_veg_allProvince_allClasses.csv"))
write.csv(a2, file.path(ROOT, VER, "R/transitions_soil_South.csv"))
write.csv(a3, file.path(ROOT, VER, "R/transitions_veg_North.csv"))


con1 <- file(file.path(ROOT, VER, "R/transitions_veg_byNSR.csv"), "w")
con2 <- file(file.path(ROOT, VER, "R/transitions_soil_byNSR.csv"), "w")
for (i in levels(resvv$gr)) {
    cat("Natural Subregion:", i, "\n", file=con1)
    cat("Natural Subregion:", i, "\n", file=con2)
    b1 <- mat_fun(i, "veg")
    b1[is.na(b1)] <- 0
    b2 <- mat_fun(i, "soil")
    b2[is.na(b2)] <- 0
    write.csv(b1, file=con1)
    write.csv(b2, file=con2)
    cat("\n", file=con1)
    cat("\n", file=con2)
}
close(con1)
close(con2)


zz <- file("ex.data", "w")  # open an output file connection
cat("TITLE extra line", "2 3 5 7", "", "11 13 17", file = zz, sep = "\n")
cat("One more line\n", file = zz)
close(zz)
readLines("ex.data")
unlink("ex.data")



## checks

library(mefa4)
## root directory
ROOT <- "c:/p"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2014"

load(file.path(ROOT, VER, "R/veg_current_QSlevel.Rdata"))
load(file.path(ROOT, VER, "R/veg_reference_QSlevel.Rdata"))
load(file.path(ROOT, VER, "R/soil_current_QSlevel.Rdata"))
load(file.path(ROOT, VER, "R/soil_reference_QSlevel.Rdata"))
load(file.path(ROOT, VER, "R/xy_clim_regions_QSlevel.Rdata"))

lu <- read.csv(file.path(ROOT, VER, "lookup/VEG_HF_interim.csv"))
su <- read.csv(file.path(ROOT, VER, "lookup/SOIL_HF_interim.csv"))


plot_subset <- function(z, ...) {
    points(QS$POINT_X[z], QS$POINT_Y[z], pch=15, cex=0.3, col=2, ...)
    invisible(NULL)
}
plot_base <- function(main="", ...) {
    plot(QS$POINT_X, QS$POINT_Y, pch=15, cex=0.3, col="lightgrey", 
        ann=FALSE, axes=FALSE, ...)
    title(main=main)
    invisible(NULL)
}

plot_brdr <- function() {
    nr1 <- QS$NEAR_DIST_GRA_km<1 & QS$NEAR_DIST_GRA_km>0
    nr1[is.na(nr1)] <- FALSE
    plot_subset(nr1)
    nr1 <- QS$NEAR_PARK_GRA_km<1 & QS$NEAR_PARK_GRA_km>0
    nr1[is.na(nr1)] <- FALSE
    plot_subset(nr1)
    nr1 <- QS$NEAR_FOOTHILL_GRA_km<1 & QS$NEAR_FOOTHILL_GRA_km>0
    nr1[is.na(nr1)] <- FALSE
    plot_subset(nr1)
    nr1 <- QS$NEAR_DIST_BoSh_km<1 & QS$NEAR_DIST_BoSh_km>0
    nr1[is.na(nr1)] <- FALSE
    plot_subset(nr1)


plot_all <- function(z, main="", ...) {
    br <- seq(0, max(z, na.rm=TRUE), len=11)
    txt <- paste0(format(br[-length(br)], digits=3, nsmall=1), "-", format(br[-1], digits=3, nsmall=1))
    br[1] <- br[1] - 1
    cz <- cut(z, br)
    #col <- rev(terrain.colors(nlevels(cz)+5)[1:nlevels(cz)])
    col <- rev(terrain.colors(nlevels(cz)))
    plot(QS$POINT_X, QS$POINT_Y, pch=15, cex=0.3, col=col[cz], 
        ann=FALSE, axes=FALSE)
    legend("bottomleft", fill=rev(col), legend=rev(txt), bty="n")
    title(main=main)
    invisible(NULL)
}

setwd(file.path(ROOT, VER, "data/veghf/qs/plots"))

png("Climate-missing.png", width=500, height=1000)
op <- par(mar=c(0, 1, 3, 0) + 0.1)
plot_base(main="Missing climate data")
plot_subset(is.na(QS$MAP))
par(op)
dev.off()

png("NatReg-missing.png", width=500, height=1000)
op <- par(mar=c(0, 1, 3, 0) + 0.1)
plot_base(main="Missing NR")
plot_subset(is.na(QS$NR))
par(op)
dev.off()

png("Aspen-missing.png", width=500, height=1000)
op <- par(mar=c(0, 1, 3, 0) + 0.1)
plot_base(main="Missing Aspen data")
plot_subset(is.na(QS$pAspen_mean))
par(op)
dev.off()

cn <- as.character(lu$Levels3)[match(colnames(veg_current), lu$VEGHFAGE)]
dat <- groupSums(veg_current, 2, cn)
dat <- as.matrix(dat)
for (i in colnames(dat)) {
    cat(i, "\n");flush.console()
    png(paste0("Current-veg-", i, ".png"), width=500, height=1000)
    op <- par(mar=c(0, 1, 3, 0) + 0.1)
    plot_all(dat[,i], main=paste0("Current-veg-", i))
    par(op)
    dev.off()
}

cn <- as.character(lu$Levels3)[match(colnames(veg_reference), lu$VEGHFAGE)]
dat <- groupSums(veg_reference, 2, cn)
dat <- as.matrix(dat)
for (i in colnames(dat)) {
    cat(i, "\n");flush.console()
    png(paste0("Reference-veg-", i, ".png"), width=500, height=1000)
    op <- par(mar=c(0, 1, 3, 0) + 0.1)
    plot_all(dat[,i], main=paste0("Reference-veg-", i))
    par(op)
    dev.off()
}

cn <- as.character(su$Levels2)[match(colnames(soil_current), su$SOILclass)]
dat <- groupSums(soil_current, 2, cn)
dat <- as.matrix(dat)
for (i in colnames(dat)) {
    cat(i, "\n");flush.console()
    png(paste0("Current-soil-", i, ".png"), width=500, height=1000)
    op <- par(mar=c(0, 1, 3, 0) + 0.1)
    plot_all(dat[,i], main=paste0("Current-soil-", i))
    par(op)
    dev.off()
}

cn <- as.character(su$Levels2)[match(colnames(soil_reference), su$SOILclass)]
dat <- groupSums(soil_reference, 2, cn)
dat <- as.matrix(dat)
for (i in colnames(dat)) {
    cat(i, "\n");flush.console()
    png(paste0("Reference-soil-", i, ".png"), width=500, height=1000)
    op <- par(mar=c(0, 1, 3, 0) + 0.1)
    plot_all(dat[,i], main=paste0("Reference-soil-", i))
    par(op)
    dev.off()
}




## maps




## -------------- soil stuff

dd <- make_vegHF_wide2(d1km, col.label = "Site_YEAR", col.year="year", wide=FALSE)
tt <- read.csv("c:/p/AB_data_v2014/lookup/VEG_HF_interim.csv")
dd$cr <- dd$VEGHFAGEclass
levels(dd$cr) <- as.character(tt$Levels3)[match(levels(dd$cr), tt$VEGHFAGE)]
dd$rf <- dd$VEGAGEclass
levels(dd$rf) <- as.character(tt$Levels3)[match(levels(dd$rf), tt$VEGHFAGE)]

tr0 <- as.matrix(Xtab(Shape_Area ~ cr + rf, dd))
tr <- round(100*t(t(tr0)/colSums(tr0)), 2)
tr[tr<1] <- 999
c <- round(100*colSums(tr0)/sum(tr0),2)
r <- c(round(100*rowSums(tr0)/sum(tr0), 2), Total=100)
tr2 <- cbind(rbind(tr, Total=c), Total=r)
write.csv(tr2, file.path(ROOT, VER, "R/transitions.csv"))


x <- dd1km$soil_reference
SITES$soil <- rowSums(x[,!(colnames(x) %in% c("UNK","Water","Len","LenW","Ltc","LtcR"))])
SITES$nonwtr <- rowSums(x[,!(colnames(x) %in% c("Water","Len","LenW","Ltc","LtcR"))])
SITES$psoil <- with(SITES, soil / nonwtr)
SITES$psoil[is.na(SITES$psoil)] <- 0
cc <- c("Grassland","Shrubland","CultivationCropPastureBareground","HighDensityLivestockOperation")
xx <- dd1km$veg_current
SITES$grass <- rowSums(xx[,cc])
plot(SITES$grass, SITES$psoil)


by(SITES[,c("soil","grass")], SITES$NR, summary)
by(SITES[SITES$NR=="Boreal",c("soil","grass")], droplevels(SITES$NSR[SITES$NR=="Boreal"]), summary)
plot(REG[,6:5],pch=ifelse(REG$NATURAL_SUBREGIONS=="Dry Mixedwood",19,21))
## incorporate soil %
## separate the 2 blobs by diagonal (lat*long function)

rn <- rownames(dd1km$current)
all(rn == rownames(dd1ha$current))

pdf(file.path(ROOT, VER, "results/tmp.pdf"), onefile=TRUE)
for (i in colnames(dd1ha$reference)) {
    cr1km <- dd1km$current[,i]
    cr1ha <- dd1ha$current[,i]
    rf1km <- dd1km$reference[,i]
    rf1ha <- dd1ha$reference[,i]
    boxplot(cbind(cr1km,cr1ha,rf1km,rf1ha), main=i)
    #plot(cr1km, rf1km)
    #points(cr1ha, rf1ha, col=2)
}
for (i in setdiff(colnames(dd1ha$current), colnames(dd1ha$reference))) {
    cr1km <- dd1km$current[,i]
    cr1ha <- dd1ha$current[,i]
    boxplot(cbind(cr1km,cr1ha), main=i)
    #plot(cr1km, rf1km)
    #points(cr1ha, rf1ha, col=2)
}
for (i in colnames(dd1ha$soil)) {
    rf1km <- dd1km$soil[,i]
    rf1ha <- dd1ha$soil[,i]
    boxplot(cbind(rf1km,rf1ha), main=i)
    #plot(cr1km, rf1km)
    #points(cr1ha, rf1ha, col=2)
}
dev.off()



plot(dd1km$current[,""]


## --
library(mefa4)
load("veg_current_QSlevel.Rdata")
lt <- read.csv("lookup/VEG_HF_interim.csv")
rownames(lt) <- lt[,1]
cn <- as.character(lt$Levels4)[match(colnames(veg_current), rownames(lt))]
## leave it as sparse, that speeds up things a bit
vhf <- groupSums(veg_current, 2, cn)


## HF 2012

library(mefa4)
ROOT <- "c:/p"
VER <- "AB_data_v2014"
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900

hftypes <- read.csv(file.path(ROOT, VER, "lookup/HFtype_lookup_20140514.csv"))
hfgroups <- read.csv(file.path(ROOT, VER, "lookup/HFclassification_20140514.csv"))
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY

load(file.path(ROOT, VER, "R/xy_clim_regions_QSlevel.Rdata"))

setdiff(d10$FEATURE_TY, rownames(hflt))
setdiff(d12$FEATURE_TY, rownames(hflt))

d10 <- read.csv(file.path(ROOT, VER, "data/veghf", "HF2010_int_QS.csv"))
d12 <- read.csv(file.path(ROOT, VER, "data/veghf", "HF2012_int_QS.csv"))

d10$HF <- hflt$HF_GROUP[match(d10$FEATURE_TY, rownames(hflt))]
levels(d10$LinkID) <- c(levels(d10$LinkID), setdiff(rownames(QS), levels(d10$LinkID)))
d12$HF <- hflt$HF_GROUP[match(d12$FEATURE_TY, rownames(hflt))]
levels(d12$LinkID) <- c(levels(d12$LinkID), setdiff(rownames(QS), levels(d12$LinkID)))

hf10 <- Xtab(Shape_Area ~ LinkID + HF, d10)
hf10 <- hf10[rownames(QS),]

hf12 <- Xtab(Shape_Area ~ LinkID + HF, d12)
hf12 <- hf12[rownames(QS),]

AAm <- QS$Area_km2*10^6

range(AAm)
range(hf10)
range(hf12)

px10 <- hf10 / AAm
px12 <- hf12 / AAm

px10=hf10
px12=hf12


setwd("c:/p/AB_data_v2014/R")

for (i in 1:ncol(px12)) {
    cat(i, "\n");flush.console()
    dat <- cbind(px12[,i], px10[,i])
    #png(paste0("w2wHFcomp-", colnames(px12)[i], "-old.png"))
    png(paste0("Area-w2wHFcomp-", colnames(px12)[i], "-awa.png"))
    plot(dat, main=colnames(px12)[i], 
        #ylim=c(0,1), xlim=c(0,1),
        ylim=c(0,max(dat)), xlim=c(0,max(dat)),
        xlab="w2w HF 2012", ylab="w2w HF 2010")
    abline(0,1)
    dev.off()
}

save(hf10,hf12,file="c:/p/AB_data_v2014/R/w2wHF_2010vs2012.Rdata")

res <- list()
for (i in 1:ncol(px12)) {
    j <- colnames(px12)[i]
    zz <- round(cbind(HF2010=by(100*px10[,i], QS$NR, mean),
        HF2012=by(100*px12[,i], QS$NR, mean)), 2)
    res[[j]] <- data.frame(HFtype=j, NR=rownames(zz), zz)
}
res <- do.call(rbind, res)
rownames(res) <- NULL
res$diff <- res$HF2012 - res$HF2010
res$abs_diff <- abs(res$diff)

res2 <- list()
for (i in 1:ncol(px12)) {
    j <- colnames(px12)[i]
    zz <- round(cbind(HF2010=by(100*px10[,i], QS$NR, mean),
        HF2012=by(100*px12[,i], QS$NSR, mean)), 2)
    res2[[j]] <- data.frame(HFtype=j, NSR=rownames(zz), zz)
}
res2 <- do.call(rbind, res2)
rownames(res2) <- NULL
res2$diff <- res2$HF2012 - res2$HF2010
res2$abs_diff <- abs(res2$diff)

res[order(res$abs_diff),]
res2[order(res2$abs_diff),]

write.csv(res, file="w2wHF_2010-2012_comparison_NR-level.csv", row.names=FALSE)
write.csv(res2, file="w2wHF_2010-2012_comparison_NSR-level.csv", row.names=FALSE)

op <- par(mar=c(4,16,4,1)+0.1, las=1, mfcol=c(2,2))
boxplot(diff~HFtype,res, horizontal=TRUE, main="NR")
abline(v=0, col=2)
boxplot(diff~HFtype,res2, horizontal=TRUE, main="NSR")
abline(v=0, col=2)
boxplot(diff~NR,res, horizontal=TRUE, main="NR")
abline(v=0, col=2)
boxplot(diff~NSR,res2, horizontal=TRUE, main="NSR")
abline(v=0, col=2)
par(op)


aa <- rowSums(x) / 656000
names(aa) <- rownames(x)
aa <- sort(aa, decreasing=TRUE)

bb <- px12[,"CultivationCropPastureBareground"] - px10[,"CultivationCropPastureBareground"]
bb <- sort(bb, decreasing=TRUE)

cc <- px12[,"CutBlocks"] - px10[,"CutBlocks"]
cc <- sort(cc, decreasing=TRUE)

m <- matrix(c(1, 0.97,
    0.8, 0.7,
    -0.7, -0.8,
    -0.9, -1), 2, 4)
db <- apply(m, 2, function(z) names(bb)[bb < z[1] & bb > z[2]][1:10])
dc <- apply(m, 2, function(z) names(cc)[cc < z[1] & cc > z[2]][1:10])
colnames(db) <- colnames(dc) <- c("12only","12>10","10>12","10only")


## HF 2010 for HF letter

library(mefa4)
ROOT <- "c:/p"
VER <- "AB_data_v2014"
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900

hftypes <- read.csv(file.path(ROOT, VER, "lookup/HFtype_lookup_20140514.csv"))
hfgroups <- read.csv(file.path(ROOT, VER, "lookup/HFclassification_20140514.csv"))
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY
rownames(hfgroups) <- hfgroups$HF_GROUP

load(file.path(ROOT, VER, "R/xy_clim_regions_QSlevel.Rdata"))
load(file.path(ROOT, VER, "R/veg_current_QSlevel.Rdata"))

NR <- QS$NR
levels(NR)[levels(NR) %in% c("Boreal","Canadian Shield")] <- "Boreal & Canadian Shield"

x <- veg_current * rs_veg_current
x <- as.matrix(x)
x <- groupSums(x, 1, NR)
x <- x[rownames(x) != "",]
Ar <- rowSums(x)
AA <- sum(x)
cn <- colnames(x)
cn[substr(cn, 1, 2) == "CC"] <- "CutBlocks"
x <- groupSums(x, 2, cn)
#Ac <- colSums(x)

x <- x[,colnames(x) %in% rownames(hfgroups)]
cn <- as.character(hfgroups$UseInReporting)[match(colnames(x), rownames(hfgroups))]
x <- groupSums(x, 2, cn)

xx <- colSums(x)
xxx <- rowSums(x)
xxxx <- sum(x)

x <- 100*x/Ar
xx <- 100*xx/AA

xxx <- 100*xxx/Ar
xxxx <- 100*xxxx/AA

x <- rbind(x, Total=xx)
x <- cbind(x, Total=c(xxx, xxxx))
x <- t(x)
x <- round(x, 2)

write.csv(x, file="c:/Dropbox/abmi/w2wHF2010_byNR_forHFletter.csv")


## AlPac 3x7 yearly stuff

ROOT <- "c:/p"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2014"
## current year
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900

library(mefa4)
#setwd("c:/Dropbox/abmi/intactness/dataproc")
source("c:/Dropbox/abmi/intactness/dataproc/data_proc_common_2014.R")

hftypes <- read.csv(file.path(ROOT, VER, "lookup/HFtype_lookup_20140514.csv"))
hfgroups <- read.csv(file.path(ROOT, VER, "lookup/HFclassification_20140514.csv"))
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY

fl <- list.files(file.path(ROOT, VER, "data/veghf/alpac"))

#i <- 1
for (i in 1:length(fl)) {
    cat(i, "\n")
    flush.console()
    fp <- file.path(ROOT, VER, "data/veghf/alpac", fl[i])

    d <- read.csv(fp)
    print(colnames(d))
    cat("\n")
    flush.console()

    if (!("year" %in% colnames(d)))
        d$year <- d$YEAR
    if (any(is.na(d$year)))
        d$year[is.na(d$year)] <- 0
    if (any(d$year==1))
        d$year[d$year==1] <- 0
    d$HF_Year <- d$year
    yr <- as.integer(strsplit(fl[i], "_")[[1]][4])
    d$YEAR <- yr
    type <- strsplit(fl[i], "_")[[1]][5]
    d$Site_YEAR <- interaction(d$ABMI, yr, drop=TRUE, sep="_")

    dd <- make_vegHF_wide2(d, col.label = "Site_YEAR", col.year="YEAR")
    dd$year <- yr
    dd$bound_type <- type

    save(dd, file=paste0(file.path(ROOT, VER, "data/veghf/alpac", fl[i]), ".Rdata"))
}

## all provice 3x7 yearly stuff

fl <- list.files(file.path(ROOT, VER, "data/veghf/allprov3x7"))

#i <- 1
for (i in 1:length(fl)) {
    gc()
    cat(i, "\n")
    flush.console()
    fp <- file.path(ROOT, VER, "data/veghf/allprov3x7", fl[i])

    d <- read.csv(fp)
    print(colnames(d))
    cat("\n")
    flush.console()

    if (!("year" %in% colnames(d)))
        d$year <- d$YEAR
    if (any(is.na(d$year)))
        d$year[is.na(d$year)] <- 0
    if (any(d$year==1))
        d$year[d$year==1] <- 0
    d$HF_Year <- d$year
    yr <- as.integer(substr(strsplit(fl[i], "_")[[1]][4], 1, 4))
    d$YEAR <- yr
    d$Site_YEAR <- interaction(d$ABMI, yr, drop=TRUE, sep="_")

    dd <- make_vegHF_wide2(d, col.label = "Site_YEAR", col.year="YEAR")
    dd$year <- yr

    save(dd, file=paste0(file.path(ROOT, VER, "data/veghf/allprov3x7", fl[i]), ".Rdata"))
}

## mammals
fl <- c("c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2001.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2002.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2003.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2004.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2005.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2006.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2007.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2008.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2009.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2010.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2011.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2012.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_2013.csv")

full_res <- list()
#i <- 1
for (i in 1:length(fl)) {
    cat(i, "\n");flush.console()
    d <- read.csv(fl[[i]])
    yr <- as.integer(substr(strsplit(fl[i], "_")[[1]][7], 1, 4))
    d$YEAR <- yr
    d$Site_YEAR <- interaction(d$ABMISite, yr, drop=TRUE, sep="_")
    dd <- make_vegHF_wide2(d, col.label = "Site_YEAR", col.year="YEAR")
    dd$sample_year <- yr
    full_res[[as.character(yr)]] <- dd
}



fl <- c("c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2001.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2002.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2003.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2004.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2005.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2006.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2007.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2008.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2009.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2010.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2011.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2012.csv",
    "c:/p/AB_data_v2014/data/veghf/mammal/veg_hf_mammal_transect_interLevel2013.csv")

int_res <- list()
#i <- 1
for (i in 1:length(fl)) {
    cat(i, "\n");flush.console()
    d <- read.csv(fl[[i]])
    yr <- as.integer(substr(strsplit(fl[i], "_")[[1]][7], 11, 14))
    d$YEAR <- yr
    d$Site_YEAR_inter <- interaction(d$ABMISite_1, yr, d$interLevel, drop=TRUE, sep="_")
    dd <- make_vegHF_wide2(d, col.label = "Site_YEAR_inter", col.year="YEAR")
    dd$sample_year <- yr
    int_res[[as.character(yr)]] <- dd
}

mammal_transects <- full_res[[1]]
mammal_transects$sample_year <- NA
for (i in 2:length(full_res)) {
    mammal_transects$veg_current <- rBind(mammal_transects$veg_current, full_res[[i]]$veg_current)
    mammal_transects$veg_reference <- rBind(mammal_transects$veg_reference, full_res[[i]]$veg_reference)
    mammal_transects$soil_current <- rBind(mammal_transects$soil_current, full_res[[i]]$soil_current)
    mammal_transects$soil_reference <- rBind(mammal_transects$soil_reference, full_res[[i]]$soil_reference)
}
sapply(full_res,function(z) dim(z[[1]])[1])

mammal_inter <- int_res[[1]]
mammal_inter$sample_year <- NA
for (i in 2:length(int_res)) {
    mammal_inter$veg_current <- rBind(mammal_inter$veg_current, int_res[[i]]$veg_current)
    mammal_inter$veg_reference <- rBind(mammal_inter$veg_reference, int_res[[i]]$veg_reference)
    mammal_inter$soil_current <- rBind(mammal_inter$soil_current, int_res[[i]]$soil_current)
    mammal_inter$soil_reference <- rBind(mammal_inter$soil_reference, int_res[[i]]$soil_reference)
}
sapply(int_res,function(z) dim(z[[1]])[1])

save(mammal_transects, mammal_inter, 
    file=paste0(file.path(ROOT, VER, "results/mammals_veghf"), ".Rdata"))

## wetlands for MC

fl <- c("c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2007.csv",
    "c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2008.csv",
    "c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2009.csv",
    "c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2010.csv",
    "c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2011.csv",
    "c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2012.csv",
    "c:/p/AB_data_v2014/data/veghf/wetland/veg_hf_on_wetland_site_2013.csv")
ba <- read.csv("c:/p/AB_data_v2014/data/veghf/wetland/wetland_sites_bufferArea_allYears.csv")
zc <- read.csv("c:/p/AB_data_v2014/data/veghf/wetland/ZoneCategories_on_wetland_sites_allYears.csv")

w <- list()
#i <- 1
for (i in 1:length(fl)) {
    cat(i, "\n");flush.console()
    d <- read.csv(fl[[i]])
    yr <- as.integer(substr(strsplit(fl[i], "_")[[1]][8], 1, 4))
    #d$YEAR <- yr
    levels(d$FEATURE_TY) <- sub(' +$', '', levels(d$FEATURE_TY))
    d$Site_Year_Dist <- interaction(d$Site_ID_DPan, d$year_, d$distance, drop=TRUE, sep="_")
    dd <- try(make_vegHF_wide2(d, col.label = "Site_Year_Dist", col.year="year_"))
    dd$sample_year <- yr
    w[[as.character(yr)]] <- dd
}
sapply(w,function(z) dim(z[[1]])[1])

ww <- w[[1]]
ww$sample_year <- NA
for (i in 2:length(w)) {
    ww$veg_current <- rBind(ww$veg_current, w[[i]]$veg_current)
    ww$veg_reference <- rBind(ww$veg_reference, w[[i]]$veg_reference)
    ww$soil_current <- rBind(ww$soil_current, w[[i]]$soil_current)
    ww$soil_reference <- rBind(ww$soil_reference, w[[i]]$soil_reference)
}

r <- strsplit(rownames(ww[[1]]), "_")
rr <- data.frame(site=sapply(r, "[[", 1), year=sapply(r, "[[", 2),
    dist=sapply(r, "[[", 3))
rr$site_year <- interaction(rr$site, rr$year, sep="_", drop=TRUE)
labels <- levels(rr$site_year)

ww1 <- ww2 <- ww3 <- ww12 <- ww123 <- ww
for (i in 1:4) {
    ww1[[i]] <- ww[[i]][rr$dist=="20",]
    rownames(ww1[[i]]) <- rr$site_year[rr$dist=="20"]
    ww1[[i]] <- ww1[[i]][labels,]

    ww2[[i]] <- ww[[i]][rr$dist=="100",]
    rownames(ww2[[i]]) <- rr$site_year[rr$dist=="100"]
    ww2[[i]] <- ww2[[i]][labels,]

    ww3[[i]] <- ww[[i]][rr$dist=="250",]
    rownames(ww3[[i]]) <- rr$site_year[rr$dist=="250"]
    ww3[[i]] <- ww3[[i]][labels,]

    ww12[[i]] <- ww1[[i]] + ww2[[i]]
    ww123[[i]] <- ww1[[i]] + ww2[[i]] + ww3[[i]]

    ww1[[i]] <- ww1[[i]]/rowSums(ww1[[i]])
    ww12[[i]] <- ww12[[i]]/rowSums(ww12[[i]])
    ww123[[i]] <- ww123[[i]]/rowSums(ww123[[i]])
}

ww1$buffer <- "0-20 m buffer around wetlands"
ww12$buffer <- "0-100 m buffer around wetlands"
ww123$buffer <- "0-250 m buffer around wetlands"

#ba$site_year <- interaction(ba$Site_ID_DPan, ba$year_, sep="_", drop=TRUE)
#rownames(ba) <- ba$site_year

save(ww1, ww12, ww123, file=file.path(ROOT, VER, "R/veghf_abmiWetlands_allbuffers.Rdata"))



