HF_VERSION <- "2014_coarse" # load 2014 defs
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))

## wetland zones ----------------------------------------------------

## Buffers

f <- file.path(ROOT, VER, "data", "raw", "veghf",
    "wetland", "VegV6VerifiedHF_summaryOnBufRings_allYear.csv")
dw250m <- read.csv(f)

## site label issue from 2016
levels(dw250m$Pin_Wetland_ID)[levels(dw250m$Pin_Wetland_ID) == "W-936"] <- "W-956"

dw250m$Site_YEAR <- with(dw250m,
    interaction(Pin_Wetland_ID, survey_Year, sep="_", drop=TRUE))
head(dw250m)
table(dw250m$BUFF_DIST)

setdiff(levels(dw250m$FEATURE_TY), levels(hftypes$FEATURE_TY))

dw20m <- dw250m[dw250m$BUFF_DIST <= 20,]
dw100m <- dw250m[dw250m$BUFF_DIST <= 100,]

ddw20m <- make_vegHF_wide_v6(dw20m,
    col.label="Site_YEAR",
    col.year="survey_Year",
    col.HFyear="YEAR_1",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
ddw20m$scale <- "0-20 m buffer around wetlands"
dx <- nonDuplicated(dw20m, Site_YEAR, TRUE)[rownames(ddw20m[[1]]),]
ddw20m <- fill_in_0ages_v6(ddw20m, dx$NSRNAME, ages_list)

ddw100m <- make_vegHF_wide_v6(dw100m,
    col.label="Site_YEAR",
    col.year="survey_Year",
    col.HFyear="YEAR_1",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
ddw100m$scale <- "0-100 m buffer around wetlands"
dx <- nonDuplicated(dw100m, Site_YEAR, TRUE)[rownames(ddw100m[[1]]),]
ddw100m <- fill_in_0ages_v6(ddw100m, dx$NSRNAME, ages_list)

ddw250m <- make_vegHF_wide_v6(dw250m,
    col.label="Site_YEAR",
    col.year="survey_Year",
    col.HFyear="YEAR_1",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
ddw250m$scale <- "0-250 m buffer around wetlands"
dx <- nonDuplicated(dw250m, Site_YEAR, TRUE)[rownames(ddw250m[[1]]),]
ddw250m <- fill_in_0ages_v6(ddw250m, dx$NSRNAME, ages_list)

sites <- droplevels(dx[,c("Site_YEAR","Pin_Wetland_ID","survey_Year",
    "LUF_NAME", "BASIN", "NRNAME", "NSRNAME")])

## catchment

f <- file.path(ROOT, VER, "data", "raw", "veghf",
    "wetland", "VegV6VerifiedHF_summaryOnCatchment_allYear.csv")

dwCat <- read.csv(f)

## site label issue from 2016 -- already cleaned up
#levels(dwCat$Pin_Wetland_ID)[levels(dwCat$Pin_Wetland_ID) == "W-936"] <- "W-956"

dwCat$Site_YEAR <- with(dwCat,
    interaction(Pin_Wetland_ID, year, sep="_", drop=TRUE))

setdiff(levels(dwCat$FEATURE_TY), levels(hftypes$FEATURE_TY))

ddwCat <- make_vegHF_wide_v6(dwCat,
    col.label="Site_YEAR",
    col.year="year",
    col.HFyear="YEAR_1",
    sparse=TRUE, HF_fine=FALSE) # don't use refined classes
ddwCat$scale <- "Catchment around wetlands"
dx <- nonDuplicated(dwCat, Site_YEAR, TRUE)[rownames(ddwCat[[1]]),]
ddwCat <- fill_in_0ages_v6(ddwCat, dx$NSRNAME, ages_list)

all(rownames(ddw20m[[1]]) == rownames(ddw100m[[1]]))
all(rownames(ddw20m[[1]]) == rownames(ddw250m[[1]]))

compare_sets(rownames(ddw20m[[1]]), rownames(ddwCat[[1]]))

save(ddw20m, ddw100m, ddw250m, ddwCat, sites,
    file=file.path(ROOT, VER, "data", "analysis", "wetland",
    "veg-hf_wetland_v6-fixage0.Rdata"))









## xy/clim/etc

climWet <- read.csv(file.path(ROOT, VER, "data", "climate",
    "climates_on_wetlandPin.csv"))
climWet2 <- read.csv(file.path(ROOT, VER, "data", "veghf", "wetlands",
    "wetlandSite2015_climates_Luf_NatReg.csv"))
climWet$ABMISite <- NULL
colnames(climWet)[colnames(climWet) == "Pin_Wetland_ID"] <- "Wetland_ID"
compare_sets(colnames(climWet), colnames(climWet2))
climWet <- rbind(climWet, climWet2[,colnames(climWet)])
levels(climWet$Wetland_ID) <- toupper(levels(climWet$Wetland_ID))

colnames(climWet)[colnames(climWet) == "Eref"] <- "PET"
colnames(climWet)[colnames(climWet) ==
    "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
climWet$OBJECTID <- NULL
climWet$Site_YEAR <- with(climWet,
    interaction(Wetland_ID, year, sep="_", drop=TRUE))

compare_sets(rownames(ddw250m[[1]]), levels(climWet$Site_YEAR))
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

fsw <- file.path(ROOT, VER, "data", "veghf", "wetlands",
    "sketch_inter_BufRings_allYearMerged.csv")
dsw <- read.csv(fsw)

fsw2 <- file.path(ROOT, VER, "data", "veghf", "wetlands",
    "sketch_inter_BufRings_Year2015.csv")
dsw2 <- read.csv(fsw2)

dsw$LinkID <- NULL
compare_sets(colnames(dsw), colnames(dsw2))
setdiff(colnames(dsw), colnames(dsw2))
setdiff(colnames(dsw2), colnames(dsw))
dsw <- rbind(dsw, dsw2)

dsw$Site_YEAR <- with(dsw,
    interaction(Pin_Wetland_ID, Year_, sep="_", drop=TRUE))

## fix age 0 in saved files -----------------------------
load(file.path(ROOT, VER, "out", "kgrid", "veg-hf_avgages_fix-fire.Rdata"))

sum(ddw20m[[1]][,Target0])
ddw20m <- fill_in_0ages(ddw20m, climWet$NSRNAME)
sum(ddw20m[[1]][,Target0])

sum(ddw100m[[1]][,Target0])
ddw100m <- fill_in_0ages(ddw100m, climWet$NSRNAME)
sum(ddw100m[[1]][,Target0])

sum(ddw250m[[1]][,Target0])
ddw250m <- fill_in_0ages(ddw250m, climWet$NSRNAME)
sum(ddw250m[[1]][,Target0])

if (SAVE)
    save(ddw20m, ddw100m, ddw250m, climWet,
        file=file.path(ROOT, VER, "out", "wetlands",
        "veg-hf_wetlands_fix-fire_fix-age0.Rdata"))

