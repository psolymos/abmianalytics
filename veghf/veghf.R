##% Processing backfilled veg + HF (cutblock ages incorporated)
##% P Solymos
##% April 28, 2015

SAVE <- TRUE

## root directory
ROOT <- "c:/p"
## version (structure is still in change, so not really useful)
VER <- "AB_data_v2015"
## current year
THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900

library(mefa4)
source("~/repos/abmianalytics/R/veghf_functions.R")
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

hftypes <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type.csv")
hfgroups <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
hflt <- hfgroups[match(hftypes$HF_GROUP, hfgroups$HF_GROUP),]
rownames(hflt) <- hftypes$FEATURE_TY

#### Vegetation and HF processing

### ABMI on+off grid sites --------------------------------------------------

## ABMI sites (on+off) cetre 1 ha
f1ha <- file.path(ROOT, VER, "data/veghf", "Center1haFixFire.csv")
d1ha <- read.csv(f1ha)
d1ha$Site_YEAR <- with(d1ha, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d1ha)
dd1ha <- make_vegHF_wide(d1ha, col.label = "Site_YEAR", 
    col.year="survey_year", col.HFyear="year_")
dd1ha$scale <- "1 ha square around site centre"

## ABMI sites (on+off) 9 bird points / site, 150 m radius buffer
f150m <- file.path(ROOT, VER, "data/veghf", "Bird150mFixFire.csv")
d150m <- read.csv(f150m)
d150m$Site_YEAR_bird <- with(d150m, 
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d150m)
dd150m <- make_vegHF_wide(d150m, col.label = "Site_YEAR_bird", 
    col.year="survey_year", col.HFyear="year_")
dd150m$scale <- "150 m radius circle around bird points"

## ABMI sites (on+off) 9 bird points / site, 1 km^2 buffer
f1km <- file.path(ROOT, VER, "data/veghf", "Bird564mFixFire.csv")
d1km <- read.csv(f1km)
d1km$Site_YEAR_bird <- with(d1km, 
    interaction(Site_ID, Bird, sep="_", drop=TRUE))
head(d1km)
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


### Snow transects -------------------------------------------------------

## 1 km length (250 m buffer) mammal transect (inter level)
fmi <- file.path(ROOT, VER, "data/veghf", "InterLevel_SRDFireFix.csv")
dmi <- read.csv(fmi)
dmi$Site_YEAR_tr <- with(dmi, interaction(ABMISite, survey_year, interLevel, sep="_", drop=TRUE))
head(dmi)
ddmi <- make_vegHF_wide(dmi, col.label = "Site_YEAR_tr", 
    col.year="survey_year", col.HFyear="year_")
ddmi$scale <- "inter level mammal transects"

## 9-10 km length (250 m buffer) mammal transect (full transect level)
fmt <- file.path(ROOT, VER, "data/veghf", "TransectLevel_SRDFireFix.csv")
dmt <- read.csv(fmt)
dmt$Site_YEAR <- with(dmt, interaction(ABMISite, survey_year, sep="_", drop=TRUE))
## strange site issue: "394-2005_2005" --> "394-2005_2006"
levels(dmt$Site_YEAR)[levels(dmt$Site_YEAR)=="394-2005_2005"] <- "394-2005_2006"
head(dmt)
ddmt <- make_vegHF_wide(dmt, col.label = "Site_YEAR", 
    col.year="survey_year", col.HFyear="year_")
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

if (SAVE)
    save(ddmi, ddmt, climInter, climTr,
        file=file.path(ROOT, VER, "out/abmi_onoff", 
        "veg-hf-clim-reg_mammals-onoff_fix-fire.Rdata"))


### BAM+BBS bird points, 150 m radius buffer --------------------------------

## processing csv files
fl <- list.files(file.path(ROOT, VER, "data", "veghf", "bammbbs150m"))
tmplist <- list()
for (fn in fl) {
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "veghf", "bammbbs150m", fn)
    d <- read.csv(f)
    hfc <- "year"
    if (!(hfc %in% colnames(d)))
        hfc <- "YEAR"
    dd <- make_vegHF_wide(d, col.label="PKEY", 
        col.year="YEAR_", col.HFyear=hfc, sparse=TRUE)
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
    hfc <- "year"
    if (!(hfc %in% colnames(d)))
        hfc <- "YEAR"
    dd <- make_vegHF_wide(d, col.label="PKEY", 
        col.year="YEAR_", col.HFyear=hfc, sparse=TRUE)
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

climPoint_bambbs <- read.csv(
    file.path(ROOT, VER, "data/climate", "AllBird_fromPeter_april2015_climates.csv"))
colnames(climPoint_bambbs)[colnames(climPoint_bambbs) == "Eref"] <- "PET"
colnames(climPoint_bambbs)[colnames(climPoint_bambbs) == 
    "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
climPoint_bambbs$OBJECTID <- NULL
colnames(climPoint_bambbs)[colnames(climPoint_bambbs) == "YEAR_"] <- "YEAR"
rownames(climPoint_bambbs) <- climPoint_bambbs$PKEY

compare.sets(rownames(dd150m_bambbs[[1]]), rownames(dd1km_bambbs[[1]]))
compare.sets(rownames(climPoint_bambbs), rownames(dd150m_bambbs[[1]]))

all(rownames(dd150m_bambbs[[1]]) == rownames(dd1km_bambbs[[1]]))
climPoint_bambbs <- climPoint_bambbs[rownames(dd150m_bambbs[[1]]),]
all(rownames(dd150m_bambbs[[1]]) == rownames(climPoint_bambbs))

if (SAVE)
    save(dd150m_bambbs, dd1km_bambbs, climPoint_bambbs,
        file=file.path(ROOT, VER, "out/bambbs", "veg-hf_bambbs_fix-fire.Rdata"))

## wetland zones ----------------------------------------------------

fw <- file.path(ROOT, VER, "data/veghf/wetlands", 
    "VerifiedHF_Veg_onWetSitesBuffferRings_allYear_July14_2015.csv")
dw250m <- read.csv(fw)
dw250m$Site_YEAR <- with(dw250m, 
    interaction(Pin_Wetland_ID, Year_survey, sep="_", drop=TRUE))
head(dw250m)
table(dw250m$BUFF_DIST)
setdiff(levels(dw250m$FEATURE_TY), levels(hftypes$FEATURE_TY))
levels(dw250m$FEATURE_TY) <- gsub(" ", "", levels(dw250m$FEATURE_TY))
levels(dw250m$FEATURE_TY) <- toupper(levels(dw250m$FEATURE_TY))
setdiff(levels(dw250m$FEATURE_TY), levels(hftypes$FEATURE_TY))

dw20m <- dw250m[dw250m$BUFF_DIST <= 20,]
dw100m <- dw250m[dw250m$BUFF_DIST <= 100,]

ddw20m <- make_vegHF_wide(dw20m, col.label = "Site_YEAR", 
    col.year="Year_survey", col.HFyear="YEAR_cut")
ddw20m$scale <- "0-20 m buffer around wetlands"

ddw100m <- make_vegHF_wide(dw100m, col.label = "Site_YEAR", 
    col.year="Year_survey", col.HFyear="YEAR_cut")
ddw100m$scale <- "0-100 m buffer around wetlands"

ddw250m <- make_vegHF_wide(dw250m, col.label = "Site_YEAR", 
    col.year="Year_survey", col.HFyear="YEAR_cut")
ddw250m$scale <- "0-250 m buffer around wetlands"

all(rownames(ddw20m[[1]]) == rownames(ddw100m[[1]]))
all(rownames(ddw20m[[1]]) == rownames(ddw250m[[1]]))
all(rownames(ddw100m[[1]]) == rownames(ddw250m[[1]]))

climWet <- read.csv(file.path(ROOT, VER, "data/climate", 
    "climates_on_wetlandPin.csv"))
colnames(climWet)[colnames(climWet) == "Eref"] <- "PET"
colnames(climWet)[colnames(climWet) == 
    "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
climWet$OBJECTID <- NULL
climWet$Site_YEAR <- with(climWet, 
    interaction(Pin_Wetland_ID, year, sep="_", drop=TRUE))

compare.sets(rownames(ddw250m[[1]]), levels(climWet$Site_YEAR))
setdiff(rownames(ddw250m[[1]]), levels(climWet$Site_YEAR))
setdiff(levels(climWet$Site_YEAR), rownames(ddw250m[[1]]))


source("~/repos/abmianalytics/species/00globalvars_wetland.R")
sort(REJECT)

#totalA <- read.csv(file.path(ROOT, VER, "data/veghf/wetlands", 
#    "BufferRings_all_year_July14_2015.csv"))

ii <- intersect(levels(climWet$Site_YEAR), rownames(ddw250m[[1]]))
ii <- ii[ii != "W-213_2013"] # outside of AB bound

for (i in 1:4) {
    ddw20m[[i]] <- ddw20m[[i]][ii,]
    ddw100m[[i]] <- ddw100m[[i]][ii,]
    ddw250m[[i]] <- ddw250m[[i]][ii,]
}
rownames(climWet) <- climWet$Site_YEAR
climWet <- droplevels(climWet[ii,])

fsw <- file.path(ROOT, VER, "data/veghf/wetlands", 
    "sketch_inter_BufRings_allYearMerged.csv")
dsw <- read.csv(fsw)
dsw$Site_YEAR <- with(dw250m, 
    interaction(Pin_Wetland_ID, Year_survey, sep="_", drop=TRUE))



if (SAVE)
    save(ddw20m, ddw100m, ddw250m, climWet,
        file=file.path(ROOT, VER, "out/wetlands", "veg-hf_wetlands_fix-fire.Rdata"))

### 1K grid --------------------------------------------------------

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
    dd <- make_vegHF_wide(d, col.label="Row_Col", 
        col.year=NULL, col.HFyear="CutYear", wide=FALSE)
    tmp <- colSums(is.na(dd[,c("VEGAGEclass",
        "VEGHFAGEclass","SOILclass","SOILHFclass")]))
    natrack[[fn]] <- tmp
    if (tmp[1] > 0)
        break
}

## check blank HABIT cases
blank_n <- numeric(length(fl)) # no. of cases
blank_a <- numeric(length(fl)) # total area
for (i in 1:length(fl)) {
    cat("\nchecking", i, "/", length(fl));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fl[i])
    d <- read.csv(f)
    tmp <- d[d$HABIT == "",,drop=FALSE]
    if (nrow(tmp)) {
        blank_n[i] <- nrow(tmp)
        blank_a[i] <- sum(tmp$Shape_Area)
        cat("\tfound:", nrow(tmp))
    }
}

(blanks <- which(blank_n > 0))
sum(blank_a) # 1.262801 < 2 m^2
blank_a[blanks]

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
    dd <- make_vegHF_wide(d, col.label="Row_Col", 
        col.year=NULL, col.HFyear="CutYear", sparse=TRUE)
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
        dd <- make_vegHF_wide(d, col.label="Row_Col", 
            col.year=NULL, col.HFyear="CutYear", sparse=TRUE)
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

## fill-in NA values with nearest

lnas <- is.na(kgrid[,"pAspen"])
wnas <- which(!lnas)
for (i in which(lnas)) {
    j <- wnas[which.min(sqrt((kgrid$X[!lnas] - kgrid$X[i])^2 +
        (kgrid$Y[!lnas] - kgrid$Y[i])^2))]
    kgrid[i,"pAspen"] <- kgrid[j,"pAspen"]
}

cvs <- c("AHM", "PET", "FFP", "MAP", "MAT", "MCMT", "MWMT")
lnas <- is.na(kgrid[,cvs[1]])
wnas <- which(!lnas)
for (i in which(lnas)) {
    j <- wnas[which.min(sqrt((kgrid$X[!lnas] - kgrid$X[i])^2 +
        (kgrid$Y[!lnas] - kgrid$Y[i])^2))]
    kgrid[i,cvs] <- kgrid[j,cvs]
}

sum(is.na(kgrid))

kgrid$LUFxNSR <- interaction(kgrid$LUF_NAME, kgrid$NSRNAME, drop=TRUE, sep="_")
levels(kgrid$LUFxNSR) <- gsub(" ", "", levels(kgrid$LUFxNSR))

if (SAVE) {
save(dd1km_pred, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
save(kgrid,
    file=file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
}

## compiling lookup tables -------------------------

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

## compiling hab-age tables by NSR ------------------------

lt <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

summary(rowSums(dd1km_pred$veg_current)/10^6)
vhf <- groupSums(dd1km_pred$veg_current, 1, kgrid$NSRNAME)
vhf2 <- groupSums(dd1km_pred$veg_current, 1, kgrid$NRNAME)
veg <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NSRNAME)
veg2 <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NRNAME)
stopifnot(all(colnames(vhf) == rownames(lt)))

lt$AGE[lt$AGE == ""] <- NA
labs_to_keep <- is.na(lt$HF) & !is.na(lt$AGE)
cols_to_keep <- colnames(vhf)[labs_to_keep]
known_ages <- cols_to_keep[!grepl("0", cols_to_keep)]

## no CC included here
Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0", 
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0", 
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")
Target <- gsub("0", "", Target0)
Ages <- c("0", "R", as.character(1:9))

NSRs <- levels(kgrid$NSRNAME)
NRs <- nonDuplicated(kgrid, NSRNAME, TRUE)
NRs <- as.character(NRs[NSRs, "NRNAME"])
names(NRs) <- NSRs

ages_cr <- array(NA, c(length(Target), length(Ages), length(NSRs)),
    list(Target, Ages, NSRs))
ages_rf <- ages_cr

for (What in c("cr","rf")) {
    cat("\n-------------", What, "---------------")
    for (nsr in NSRs) {
        for (i in Target) {
            tmp <- if (What == "cr")
                vhf[nsr, paste0(i, Ages)] else veg[nsr, paste0(i, Ages)]
            if (tmp[1] > 0 && sum(tmp[-1]) <= 0) {## these are all <1%, just ignore
                cat("\nproblem:", nsr, "\t", i, "\t", 
                    round(100 * sum(tmp) / sum(vhf[nsr, cols_to_keep]), 4), 
                    ifelse(sum(vhf[nsr, known_ages]) > 0, "OK", "zero !!!"))
            } else {
                #cat("\nOK:", nsr, "\t", i)
            }
            if (tmp[1] > 0 && sum(tmp[-1]) <= 0) {
                ## we do NOT have known ages to work with
                if (sum(vhf[nsr, known_ages]) <= 0) {
                ## NO know age forest in NSR --> use NR
                    allages <- if (What == "cr")
                        vhf2[NRs[nsr], known_ages] else veg2[NRs[nsr], known_ages]
                } else {
                ## we have known age forest in NSR
                    allages <- if (What == "cr")
                        vhf[nsr, known_ages] else veg[nsr, known_ages]
                }
                dat <- data.frame(A=allages,
                    VEG=substr(names(allages), 1, nchar(names(allages))-1),
                    AGE=substr(names(allages), nchar(names(allages)),
                        nchar(names(allages))))
                mat <- as.matrix(Xtab(A ~ VEG + AGE, dat))
                mat <- mat[Target, Ages[-1]]
                tmp <- colSums(mat)
                tmp <- c("0"=0, tmp)
            } else {
                ## we have known ages to work with
                tmp[1] <- 0
            }
            #stopifnot(all(!is.na(tmp / sum(tmp))))
            tmp <- tmp / ifelse(sum(tmp) <= 0, 1, sum(tmp))
            stopifnot(all(!is.na(tmp)))
            if (What == "cr")
                ages_cr[i,,nsr] <- tmp
            if (What == "rf")
                ages_rf[i,,nsr] <- tmp
        }
    }
    cat("\n\n")
}

sum(is.na(ages_cr))
sum(is.na(ages_rf))

## we need to know availability for combining forest classes
nsr_cr <- groupSums(dd1km_pred$veg_current, 1, kgrid$NSRNAME)
nsr_rf <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NSRNAME)
cn1 <- as.character(lt$VEG)
cn1[!is.na(lt$HF)] <- "HF"
cn2 <- cn1[is.na(lt$HF)]
nsr_cr <- as.matrix(groupSums(nsr_cr, 2, cn1))
nsr_rf <- as.matrix(groupSums(nsr_rf, 2, cn2))
nsr_cr <- nsr_cr[,dimnames(ages_cr)[[1]]]
nsr_rf <- nsr_rf[,dimnames(ages_cr)[[1]]]

AvgAges <- list(current=ages_cr, reference=ages_rf,
    area_cr=nsr_cr, area_rf=nsr_rf)

if (SAVE)
save(AvgAges, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

## fix age 0 in saved files -----------------------------

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))
Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0", 
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0", 
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")

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

## ddmi, ddmt, climInter, climTr
load(file.path(ROOT, VER, "out/abmi_onoff", 
    "veg-hf-clim-reg_mammals-onoff_fix-fire.Rdata"))

sum(ddmi[[1]][,Target0])
ddmi <- fill_in_0ages(ddmi, climInter$NSRNAME)
sum(ddmi[[1]][,Target0])

sum(ddmt[[1]][,Target0])
ddmt <- fill_in_0ages(ddmt, climTr$NSRNAME)
sum(ddmt[[1]][,Target0])

save(ddmi, ddmt, climInter, climTr,
    file=file.path(ROOT, VER, "out/abmi_onoff", 
    "veg-hf-clim-reg_mammals-onoff_fix-fire_fix-age0.Rdata"))


## dd150m_bambbs, dd1km_bambbs -- need NSR from previous climate table

load(file.path(ROOT, VER, "out/bambbs", "veg-hf_bambbs_fix-fire.Rdata"))

sum(dd150m_bambbs[[1]][,Target0])
dd150m_bambbs <- fill_in_0ages(dd150m_bambbs, climPoint_bambbs$NSRNAME)
sum(dd150m_bambbs[[1]][,Target0])

sum(dd1km_bambbs[[1]][,Target0])
dd1km_bambbs <- fill_in_0ages(dd1km_bambbs, climPoint_bambbs$NSRNAME)
sum(dd1km_bambbs[[1]][,Target0])

save(dd150m_bambbs, dd1km_bambbs, climPoint_bambbs,
    file=file.path(ROOT, VER, "out/bambbs", "veg-hf_bambbs_fix-fire_fix-age0.Rdata"))

## Wetlands

sum(ddw20m[[1]][,Target0])
ddw20m <- fill_in_0ages(ddw20m, climWet$NSRNAME)
sum(ddw20m[[1]][,Target0])

sum(ddw100m[[1]][,Target0])
ddw100m <- fill_in_0ages(ddw100m, climWet$NSRNAME)
sum(ddw100m[[1]][,Target0])

sum(ddw250m[[1]][,Target0])
ddw250m <- fill_in_0ages(ddw250m, climWet$NSRNAME)
sum(ddw250m[[1]][,Target0])

save(ddw20m, ddw100m, ddw250m, climWet,
    file=file.path(ROOT, VER, "out/wetlands", 
    "veg-hf_wetlands_fix-fire_fix-age0.Rdata"))

## 1 km grid
load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

sum(dd1km_pred[[1]][,Target0])
sum(dd1km_pred[[2]][,Target0])
sum(dd1km_pred[[1]])
sum(dd1km_pred[[2]])
dd1km_pred <- fill_in_0ages(dd1km_pred, kgrid$NSRNAME)
sum(dd1km_pred[[1]][,Target0])
sum(dd1km_pred[[2]][,Target0])
sum(dd1km_pred[[1]])
sum(dd1km_pred[[2]])

save(dd1km_pred, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))


### Transition for 1K grid ------------------------------------------------

## this is based on the fix-fire fix-age0 version
## label collapsing as desired (swamp/wetland, ages?)

source("~/repos/abmianalytics/R/veghf_functions.R")

load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
fl <- list.files(file.path(ROOT, VER, "data", "kgrid", "tiles"))

cc <- c("Row_Col","VEGAGEclass","VEGHFAGEclass","SOILclass","SOILHFclass","Shape_Area")

Start <- c(0:79*10+1, 802)


d <- read.csv(file.path(ROOT, VER, "data", "kgrid", "tiles", fl[1]))
dd <- make_vegHF_wide(d, col.label="Row_Col", 
    col.year=NULL, col.HFyear="CutYear", wide=FALSE)
ddd0 <- dd[character(0),cc]
xddd0 <- ddd0

for (s in 1:(length(Start)-1)) {
    cat("----------------------\nStarting block", s, "\n")
    for (i in Start[s]:(Start[s+1]-1)) {
        cat(i, "of", length(fl), "-", fl[i], "\t")
        flush.console()
        d <- read.csv(file.path(ROOT, VER, "data", "kgrid", "tiles", fl[i]))
        dd <- make_vegHF_wide(d, col.label="Row_Col", 
            col.year=NULL, col.HFyear="CutYear", wide=FALSE)
        if (i == Start[s]) {
            dd0 <- dd[,cc]
        } else {
            dd0 <- rbind(dd0, dd[,cc])
        }
        cat("OK", nrow(dd0), "\n")
    }
    ddd0 <- rbind(ddd0, dd0)
    cat("\nFinished block", s, "dim:", nrow(ddd0), "\n")
    if (i %in% c(100, 200, 300, 400, 500, 600, 700, 801)) {
        save(ddd0, file=file.path(ROOT, VER, 
            "data", "kgrid", "long", paste0("Long-part", i, ".Rdata")))
        ddd0 <- xddd0
        gc()
    }
}

## -- works on pre-saved chunks

load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
lu <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
su <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

lu$use_tr <- as.character(lu$VEGAGE_use)
lu$use_tr[!is.na(lu$HF)] <- as.character(lu$VEGHFAGE[!is.na(lu$HF)])
lu$use_tr[lu$use_tr == "WetBare"] <- "NonVeg"
allVegTr <- unique(c(lu$use_tr[is.na(lu$HF)], 
    paste0(rep(lu$use_tr[is.na(lu$HF)], sum(!is.na(lu$HF))), "->",
    rep(lu$use_tr[!is.na(lu$HF)], each=sum(is.na(lu$HF))))))
lu$use_tr <- as.factor(lu$use_tr)

su$use_tr <- as.character(su$Levels1)
su$use_tr[!is.na(su$HF)] <- as.character(su$SOILHF[!is.na(su$HF)])
allSoilTr <- unique(c(su$use_tr[is.na(su$HF)], 
    paste0(rep(su$use_tr[is.na(su$HF)], sum(!is.na(su$HF))), "->",
    rep(su$use_tr[!is.na(su$HF)], each=sum(is.na(su$HF))))))
su$use_tr <- as.factor(su$use_tr)

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))
Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0", 
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0", 
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")

recl <- list(
    bf=c("Conif", "Decid", "Mixwood", "Pine", "Swamp-Conif", "Swamp-Decid", 
        "Swamp-Mixwood", "Swamp-Pine", "Wetland-BSpr", "Wetland-Decid", "Wetland-Larch"),
    target=c("Conif0", "Decid0", "Mixwood0", "Pine0", "Swamp-Conif0", "Swamp-Decid0", 
        "Swamp-Mixwood0", "Swamp-Pine0", "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0"),
    reclass=c("Conif0", "Decid0", "Mixwood0", "Pine0", "Conif0", "Decid0", 
        "Mixwood0", "Pine0", "BSpr0", "Decid0", "Larch0"))
    
fl3 <- list.files(file.path(ROOT, VER, "data", "kgrid", "long"))


## do one LUFxNSR class at a time and save it
#i <- "UpperAthabasca_CentralMixedwood"
for (ii in 1:nlevels(kgrid$LUFxNSR)) {
    i <- levels(kgrid$LUFxNSR)[ii]
    cat("\n---------", i)
    #j <- 4
    units <- list()
    sunits <- list()
    for (j in 1:length(fl3)) {
        cat("\n", j);flush.console()
        load(file.path(ROOT, VER, "data", "kgrid", "long", fl3[j]))
        flush.console()
        ddd0$LUFxNSR <- kgrid$LUFxNSR[match(ddd0$Row_Col, kgrid$Row_Col)]
        nsr <- as.character(kgrid[which(kgrid$LUFxNSR==i)[1], "NSRNAME"])
        
        if (any(ddd0$LUFxNSR == i)) {
            cat(" processing ... ")

            xx <- ddd0[ddd0$LUFxNSR == i,,drop=FALSE]
            xx$Row_Col <- droplevels(xx$Row_Col)
            xx$LUFxNSR <- NULL

            xx$soil <- su$use_tr[match(xx$SOILclass, rownames(su))]
            xx$shf <- su$use_tr[match(xx$SOILHFclass, rownames(su))]

            xx$veg <- lu$use_tr[match(xx$VEGAGEclass, rownames(lu))]
            xx$vhf <- lu$use_tr[match(xx$VEGHFAGEclass, rownames(lu))]

            xx$soilTr <- ifelse(as.character(xx$soil) == as.character(xx$shf),
                as.character(xx$soil), paste0(as.character(xx$soil),
                "->", as.character(xx$shf)))
            
            xx$vegTr <- ifelse(as.character(xx$veg) == as.character(xx$vhf),
                as.character(xx$veg), paste0(as.character(xx$veg),
                "->", as.character(xx$vhf)))

            sxt <- Xtab(Shape_Area ~ Row_Col + soilTr, xx)
            sxxx <- Melt(sxt)
            colnames(sxxx) <- c("Row_Col", "soilTr", "Shape_Area")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")
            xxx0 <- xxx[xxx$vegTr %in% Target0,,drop=FALSE]
            if (nrow(xxx0)>0) {
                cat("age0")
                xxx1 <- xxx[!(xxx$vegTr %in% Target0),,drop=FALSE]
                xxx0$vegTr <- as.character(xxx0$vegTr)
                xxx0$veg <- sapply(strsplit(as.character(xxx0$vegTr), "->"), "[[", 1)
                xxx0$vhf <- sapply(strsplit(as.character(xxx0$vegTr), "->"), 
                    function(z) z[length(z)])
                xxx0$vhf[xxx0$vhf == xxx0$veg] <- ""

                ## needs to sum to 1, include availability
                ages <- AvgAges$reference[,,nsr]
                areas <- AvgAges$area_rf[nsr,]
                bf0 <- groupMeans(ages * areas, 1, recl$reclass)[,-1]
                bf0 <- bf0 / rowSums(bf0)

                tmp <- list()
                for (k in 1:10) {
                    tmpv <- xxx0
                    target <- substr(tmpv$veg, 1, nchar(tmpv$veg)-1)
                    tmpv$Shape_Area <- tmpv$Shape_Area * bf0[match(tmpv$veg, rownames(bf0)),k]
                    tmpv$veg <- paste0(target, colnames(bf0)[k])
                    tmpv$vegTr <- ifelse(tmpv$vhf == "", tmpv$veg,
                        paste0(tmpv$veg, "->", tmpv$vhf))
                    tmp[[k]] <- tmpv[,colnames(xxx1)]
                }
                xxx0v <- do.call(rbind, tmp)
                xxx <- rbind(xxx1, xxx0v)
            }
            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xxx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")
            units[[j]] <- xxx
            sunits[[j]] <- sxxx
        } else cat(" onto the next chunk")
    }
    units <- do.call(rbind, units)
    levels(units$vegTr) <- c(levels(units$vegTr),
        setdiff(allVegTr, levels(units$vegTr)))
    sunits <- do.call(rbind, sunits)
    levels(sunits$soilTr) <- c(levels(sunits$soilTr),
        setdiff(allSoilTr, levels(sunits$soilTr)))
    
    trVeg <- Xtab(Shape_Area ~ Row_Col + vegTr, units)
    trVeg <- trVeg[,allVegTr]
    trSoil <- Xtab(Shape_Area ~ Row_Col + soilTr, sunits)
    trSoil <- trSoil[rownames(trVeg),allSoilTr]
    range(rowSums(trVeg)/10^6)
    range(rowSums(trSoil)/10^6)
    
    save(trVeg, trSoil, file=file.path(ROOT, VER, "out", "transitions", 
        paste0(i, ".Rdata")))
}


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




