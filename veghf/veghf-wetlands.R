source("~/repos/abmianalytics/veghf/veghf-setup.R")

## wetland zones ----------------------------------------------------

fw <- file.path(ROOT, VER, "data", "veghf", "wetlands", 
    "VerifiedHF_Veg_onWetSitesBuffferRings_allYear_July14_2015.csv")
dw250m <- read.csv(fw)
fw2 <- file.path(ROOT, VER, "data", "veghf", "wetlands", 
    "verifiedHF_Veg_OnWetSitesBufferRings_Yr2015_v2.csv")
dw250m2 <- read.csv(fw2)
levels(dw250m2$Pin_Wetland_ID) <- toupper(levels(dw250m2$Pin_Wetland_ID))

colnames(dw250m2)[colnames(dw250m2) == "SOURCE_1"] <- "SOURCE"
dw250m$LinkID <- NULL
compare_sets(colnames(dw250m), colnames(dw250m2))
setdiff(colnames(dw250m), colnames(dw250m2))
setdiff(colnames(dw250m2), colnames(dw250m))

dw250m <- rbind(dw250m, dw250m2)

dw250m$Site_YEAR <- with(dw250m, 
    interaction(Pin_Wetland_ID, Year_survey, sep="_", drop=TRUE))
head(dw250m)
table(dw250m$BUFF_DIST)
## 2015 FEATURE_TY labels are fine
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

