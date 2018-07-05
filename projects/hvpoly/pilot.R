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

Job <- "sppden"
#Job <- "ofbirds"

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
    rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))
    xy <- spTransform(xy, proj4string(rpa))
    x$pAspen <- extract(rpa, xy)

    ## use N/S
    x$useN <- !(x$NRNAME %in% c("Grassland", "Parkland") | x$NSRNAME == "Dry Mixedwood")
    x$useN[x$NSRNAME == "Dry Mixedwood" & x$POINT_Y > 56.7] <- TRUE
    x$useS <- x$NRNAME == "Grassland"


    #taxon <- "mammals"
    for (taxon in Taxa) {
        Species <- if (Job == "sppden")
            get_all_species(taxon=taxon) else get_all_species(taxon, mregion="north")
        SPPS <- get_all_species(taxon, mregion="south")
        SPPN <- get_all_species(taxon, mregion="north")

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

        OUT <- if (Job == "sppden")
            0 else list()
        n <- 1
        #species <- "Achillea.millefolium"
        for (species in Species) {
            cat(Job, PilotArea, taxon, species, "---", n, "/", N, "\n")
            flush.console()
            object <- load_spclim_data(species)
            pred_curr <- suppressWarnings(predict(object, xy=xy,
                veg=x$VEGHFAGEclass3, soil=x$SOILHFclass3))
            pred_curr[is.na(pred_curr)] <- 0 # water and non-veg

            if (Job == "sppden") {
                TYPE <- "C" # combo
                if (!(species %in% SPPN))
                    TYPE <- "S"
                if (!(species %in% SPPS))
                    TYPE <- "N"

                wS <- 1-x$pAspen
                if (TYPE == "S")
                    wS[] <- 1
                if (TYPE == "N")
                    wS[] <- 0
                wS[x$useS] <- 1
                wS[x$useN] <- 0

                cr <- wS * pred_curr[,"soil"] + (1-wS) * pred_curr[,"veg"]

                if (taxon == "birds")
                    cr <- 1-exp(-cr)
                #cr[is.na(cr)] <- 0 # water and non-veg
                OUT <- OUT + cr
            } else {
                pred_curr[is.na(pred_curr)] <- 0 # water and non-veg
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
scaled <- TRUE


#PilotArea <- "north"
for (PilotArea in c("south", "north")) {

    db <- dbConnect(RSQLite::SQLite(), f)
    if (PilotArea == "south") {
        x <- dbReadTable(db, "south")
    } else {
        x <- dbReadTable(db, "north")
    }
    dbDisconnect(db)
    #xy <- x[,c("xcoord", "ycoord")]
    #coordinates(xy) <- ~ xcoord + ycoord
    #proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    #xy <- spTransform(xy, proj4string(rt))

    #taxon <- "birds"
    ALL <- 0
    TAB <- data.frame(OBJECTID=x$OBJECTID)
    for (taxon in Taxa) {
        cat(PilotArea, taxon, "\n")
        flush.console()
        fn <- paste0("e:/peter/AB_data_v2018/data/raw/hvpoly/SpeciesDensity-",
            PilotArea, "-", taxon, ".RData")
        load(fn)
        if (scaled) {
            OUT <- round(100*OUT/max(OUT))
        }
        ALL <- ALL + OUT
        TAB[[taxon]] <- if (scaled) as.integer(OUT) else OUT
    }
    if (scaled) {
        ALL <- ALL / length(Taxa)
        for (i in 1:3)
            ALL <- round(100*ALL/max(ALL))
    }
    TAB[["all"]] <- if (scaled) as.integer(ALL) else ALL


    fout <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
        "polygon-tool-pilot-results.sqlite")
    con <- dbConnect(RSQLite::SQLite(), fout)
    if (PilotArea == "south") {
        dbWriteTable(con, "spden_south", TAB, overwrite = TRUE)
    } else {
        dbWriteTable(con, "spden_north", TAB, overwrite = TRUE)
    }
    #dbListTables(con)
    dbDisconnect(con)
}


f2 <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "polygon-tool-pilot-results.sqlite")
db <- dbConnect(RSQLite::SQLite(), f2)
zs <- dbReadTable(db, "spden_south")
zn <- dbReadTable(db, "spden_north")
dbDisconnect(db)
db <- dbConnect(RSQLite::SQLite(), f)
xs <- dbReadTable(db, "south")
xn <- dbReadTable(db, "north")
dbDisconnect(db)
xys <- xs[,c("xcoord", "ycoord")]
coordinates(xys) <- ~ xcoord + ycoord
proj4string(xys) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xys <- spTransform(xys, proj4string(rt))
xyn <- xn[,c("xcoord", "ycoord")]
coordinates(xyn) <- ~ xcoord + ycoord
proj4string(xyn) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xyn <- spTransform(xys, proj4string(rt))
zs <- zs[,-1]
zn <- zn[,-1]

ls <- list()
for (i in colnames(zs)) {
    rrN <- trim(rasterize(xys, rt, zs[,i] * xs$Shape_Area, fun=sum), values=NA)
    rrA <- trim(rasterize(xys, rt, xs$Shape_Area, fun=sum), values=NA)
    rr <- rrN/rrA
    values(rr)[!is.na(values(rr)) & values(rr) == max(values(rr), na.rm=TRUE)] <- 100
    ls[[i]] <- rr
}
ls <- stack(ls)

ln <- list()
for (i in colnames(zn)) {
    rrN <- trim(rasterize(xyn, rt, zn[,i] * xn$Shape_Area, fun=sum), values=NA)
    rrA <- trim(rasterize(xyn, rt, xn$Shape_Area, fun=sum), values=NA)
    rr <- rrN/rrA
    values(rr)[!is.na(values(rr)) & values(rr) == max(values(rr), na.rm=TRUE)] <- 100
    ln[[i]] <- rr
}
ln <- stack(ln)

pdf("~/GoogleWork/abmi/hvpoly/results/species-density-south.pdf")
plot(ls, colNA="grey", col=viridis::viridis(101), axes=FALSE, box=FALSE)
dev.off()

pdf("~/GoogleWork/abmi/hvpoly/results/species-density-north.pdf")
plot(ln, colNA="grey", col=viridis::viridis(101), axes=FALSE, box=FALSE)
dev.off()

