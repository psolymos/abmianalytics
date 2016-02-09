source("~/repos/abmianalytics/veghf/veghf-setup.R")

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

## fix age 0 in saved files -----------------------------

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

