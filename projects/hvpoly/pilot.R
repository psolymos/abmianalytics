## using the polygon tool
library(DBI)
library(cure4insect)
library(rgdal)
library(mefa4)
load_common_data()
set_options(path = "w:/reports")
#source("~/repos/abmianalytics/R/veghf_functions.R")

make_younger <- function(z) {
    for (i in rev(as.character(c(10, 20, 40, 60, 80, 100, 120, 140, 160)))) {
        z[grep(i, z)] <- gsub(i, "0", z[grep(i, z)])
    }
    z
}

#Job <- "sppden"
Job <- "ofbirds"

veg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
veg_mapping <- cbind(V6=as.character(veg$VEGHFAGE_FINE), V5=as.character(veg$PolyReclass))
soil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2014.csv")
soil_mapping <- cbind(In=as.character(soil$SOILHF_FINE), Out=as.character(soil$PolyReclass))

f <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")
if (Job == "sppden")
    Taxa <- c("lichens", "mammals", "mites", "mosses", "vplants", "birds")
if (Job == "ofbirds")
    Taxa <- "birds"

## species density

PilotAreas <- if (Job == "sppden")
    c("south", "north") else "north"
#PilotArea <- "south"
for (PilotArea in PilotAreas) {

    db <- dbConnect(RSQLite::SQLite(), f)
    #dbListTables(db)
    if (PilotArea == "south") {
        x <- dbReadTable(db, "south")
    } else {
        x <- dbReadTable(db, "north")
    }
    dbDisconnect(db)

    for (i in 1:ncol(x))
        if (is.character(x[,i]))
            x[,i] <- as.factor(x[,i])

    x$VEGAGEclass2 <- reclass(x$VEGAGEclass, veg_mapping, all=TRUE)
    #setdiff(x$VEGAGEclass2, get_levels()$veg)
    x$VEGHFAGEclass2 <- reclass(x$VEGHFAGEclass, veg_mapping, all=TRUE)
    #setdiff(x$VEGHFAGEclass2, get_levels()$veg) # Zero is 0 abundance

    x$SOILclass2 <- reclass(x$SOILclass, soil_mapping, all=TRUE)
    #setdiff(x$SOILclass2, get_levels()$soil)
    x$SOILHFclass2 <- reclass(x$SOILHFclass, soil_mapping, all=TRUE)
    #setdiff(x$SOILHFclass, soil_mapping[,1]) # Zero is 0 abundance

    xy <- x[,c("xcoord", "ycoord")]
    coordinates(xy) <- ~ xcoord + ycoord
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

    #taxon <- "mammals"
    for (taxon in Taxa) {
        Species <- if (Job == "sppden")
            get_all_species(taxon=taxon) else get_all_species(taxon, mregion="north")
        N <- length(Species)
        #SPDEN <- list()
        x$VEGHFAGEclass3 <- x$VEGHFAGEclass2
        x$SOILHFclass3 <- x$SOILHFclass2
        if (taxon == "birds") {
            lin1 <- x$VEGHFAGEclass2 %in% c("SoftLin", "HardLin")
            repl1 <- reclass(x$VEGAGEclass[lin1], veg_mapping, all=TRUE)
            repl1 <- make_younger(repl1)
            x$VEGHFAGEclass3[lin1] <- repl1

            lin2 <- x$SOILHFclass2 %in% c("SoftLin", "HardLin")
            repl2 <- reclass(x$SOILclass[lin2], soil_mapping, all=TRUE)
            x$SOILHFclass3[lin2] <- repl2
        }

        #species <- "Achillea.millefolium"
        OUT <- if (Job == "sppden")
            0 else list()
        n <- 1
        for (species in Species) {
            cat(Job, PilotArea, taxon, species, "---", n, "/", N, "\n")
            flush.console()
            object <- load_spclim_data(species)
            pred_curr <- suppressWarnings(predict(object, xy=xy,
                veg=x$VEGHFAGEclass3, soil=x$SOILHFclass3))
            pred_curr[is.na(pred_curr)] <- 0 # water and non-veg
            if (Job == "sppden") {
                if (taxon == "birds")
                    pred_curr <- 1-exp(-pred_curr)
                #OUT <- ((n-1)/N) * OUT  + (1/N) * pred_curr
                OUT <- OUT + pred_curr
            } else {
                OUT[[species]] <- pred_curr[,1,drop=FALSE]
            }
            n <- n + 1
        }
        #SPDEN[[taxon]] <- OUT
        if (Job == "sppden") {
            save(OUT, file=file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
                paste0("SpeciesDensity-", PilotArea, "-", taxon, ".RData")))
        } else {
            save(OUT, file=file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
                paste0("OFbirds-", PilotArea, ".RData")))
        }
    }
}
