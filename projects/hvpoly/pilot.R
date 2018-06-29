## todo: linear feature stuff for birds: ise backfill with age R

## using the polygon tool
library(DBI)
library(cure4insect)
library(rgdal)
library(mefa4)
load_common_data()
set_options(path = "w:/reports")

veg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
veg_mapping <- cbind(V6=as.character(veg$VEGHFAGE_FINE), V5=as.character(veg$PolyReclass))
soil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2014.csv")
soil_mapping <- cbind(In=as.character(soil$SOILHF_FINE), Out=as.character(soil$UseInAnalysis))

f <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")
Taxa <- c("lichens", "mammals", "mites", "mosses", "vplants")# "birds")

## species density

#PilotArea <- "south"
for (PilotArea in c("south", "north")) {

    db <- dbConnect(RSQLite::SQLite(), f)
    dbListTables(db)
    if (PilotArea == "south") {
        x <- dbReadTable(db, "south")
    } else {
        x <- dbReadTable(db, "north")
    }
    dbDisconnect(db)

    for (i in 1:ncol(x))
        if (is.character(x[,i]))
            x[,i] <- as.factor(x[,i])

    x$VEGAGEclass2 <- reclass(x$VEGAGEclass, veg_mapping)
    #setdiff(x$VEGAGEclass2, get_levels()$veg)
    x$VEGHFAGEclass2 <- reclass(x$VEGHFAGEclass, veg_mapping)
    #setdiff(x$VEGHFAGEclass2, get_levels()$veg)

    x$SOILclass2 <- reclass(x$SOILclass, soil_mapping)
    #setdiff(x$SOILclass2, get_levels()$soil)
    x$SOILHFclass2 <- reclass(x$SOILHFclass, soil_mapping)
    #setdiff(x$SOILHFclass, soil_mapping[,1])

    xy <- x[,c("xcoord", "ycoord")]
    coordinates(xy) <- ~ xcoord + ycoord
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

    #taxon <- "mammals"
    for (taxon in Taxa) {
        Species <- get_all_species(taxon=taxon)
        N <- length(Species)
        SPDEN <- list()
        x$VEGHFAGEclass3 <- x$VEGHFAGEclass2
        x$SOILHFclass3 <- x$SOILHFclass2
        if (taxon == "birds") {
            stop("fix linear!")
        #    x$VEGHFAGEclass3[lin] <- x$VEGAGEclass[lin]
        #    x$SOILHFclass3[lin] <- x$SOILclass2
        }

        #species <- "Achillea.millefolium"
        OUT <- 0
        n <- 1
        for (species in Species) {
            cat(PilotArea, taxon, species, "---", n, "/", N, "\n")
            flush.console()
            object <- load_spclim_data(species)
            pred_curr <- suppressWarnings(predict(object, xy=xy,
                veg=x$VEGHFAGEclass3, soil=x$SOILHFclass3))
            pred_curr[is.na(pred_curr)] <- 0 # water and non-veg
            if (taxon == "birds")
                pred_curr <- 1-exp(-pred_curr)
            #OUT <- ((n-1)/N) * OUT  + (1/N) * pred_curr
            OUT <- OUT + pred_curr
            n <- n + 1
        }
        SPDEN[[taxon]] <- OUT
        save(OUT, file=file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
            paste0("SpeciesDensity-", PilotArea, "-", taxon, ".RData")))
    }
}
