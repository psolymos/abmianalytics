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

## processing raw output for final layer values
fb <- "e:/peter/AB_data_v2018/data/raw/hvpoly/OFbirds-north.RData"

fn <- c("e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-north-birds.RData",
    "e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-north-vplants.RData",
    "e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-north-mosses.RData",
    "e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-north-mites.RData",
    "e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-north-mammals.RData",
    "e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-north-lichens.RData")


## using the polygon tool
library(DBI)
library(cure4insect)
library(rgdal)
library(mefa4)
f <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")
rt <- .read_raster_template()
#rt200 <- disaggregate(rt, 5)
Taxa <- c("lichens", "mammals", "mites", "mosses", "vplants", "birds")
PilotArea <- "north"
#PilotArea <- "south"

db <- dbConnect(RSQLite::SQLite(), f)
if (PilotArea == "south") {
    x <- dbReadTable(db, "south")
} else {
    x <- dbReadTable(db, "north")
}
dbDisconnect(db)
xy <- x[,c("xcoord", "ycoord")]
coordinates(xy) <- ~ xcoord + ycoord
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, proj4string(rt))

#scaled <- TRUE
for (scaled in c(TRUE, FALSE)) {
#taxon <- "birds"
ALL <- 0
pdf(paste0("~/GoogleWork/abmi/hvpoly/results/species-density-", PilotArea,
    "-", if (scaled) "scaled" else "unscaled", ".pdf"), onefile=TRUE)
for (taxon in Taxa) {
    cat(taxon, "\n");flush.console()
    fn <- paste0("e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-",
        PilotArea, "-", taxon, ".RData")
    e <- new.env()
    load(fn, envir=e)
    OUT <- e$OUT
    is0 <- OUT[,"soil"] == 0
    OUT[is0,"soil"] <- OUT[is0,"veg"]
    OUT[is0,"comb"] <- OUT[is0,"veg"]
    rm(e)
    if (scaled) {
        for (i in 1:3)
            OUT[,i] <- round(100*OUT[,i]/max(OUT[,i]))
    }
    ALL <- ALL + OUT

    rrVeg <- trim(rasterize(xy, rt, OUT[,"veg"], fun=mean), values=NA)
    rrSoil <- trim(rasterize(xy, rt, OUT[,"soil"], fun=mean), values=NA)
    rrComb <- trim(rasterize(xy, rt, OUT[,"comb"], fun=mean), values=NA)
    List <- list(rrVeg, rrSoil, rrComb)
    names(List) <- paste(c("Veg", "Soil", "Combined"), taxon)
    plot(stack(List),
        colNA="grey", col=viridis::viridis(101), axes=FALSE, box=FALSE)
}

if (scaled) {
    ALL <- ALL / length(Taxa)
    for (i in 1:3)
        ALL[,i] <- round(100*ALL[,i]/max(ALL[,i]))
}

rrVeg <- trim(rasterize(xy, rt, ALL[,"veg"], fun=mean), values=NA)
rrSoil <- trim(rasterize(xy, rt, ALL[,"soil"], fun=mean), values=NA)
rrComb <- trim(rasterize(xy, rt, ALL[,"comb"], fun=mean), values=NA)
List <- list(rrVeg, rrSoil, rrComb)
names(List) <- paste(c("Veg", "Soil", "Combined"), "All")
plot(stack(List),
    colNA="grey", col=viridis::viridis(101), axes=FALSE, box=FALSE)

dev.off()
}

