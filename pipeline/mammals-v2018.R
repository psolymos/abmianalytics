library(cure4insect)
library(mefa4)
set_options(path = "s:/reports")
load_common_data()

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y

## ch2soil ch2veg trSoil trVeg
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]

compare_sets(ch2soil$cr, rownames(ts))
setdiff(ch2soil$cr, rownames(ts))
setdiff(rownames(ts), ch2soil$cr)

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2soil$rf2 <- ts$UseInAnalysisCoarse[match(ch2soil$rf, rownames(ts))]
ch2soil$cr2 <- ts$UseInAnalysisCoarse[match(ch2soil$cr, rownames(ts))]
ch2soil$sector <- ts$Sector61[match(ch2soil$cr, rownames(ts))]

ch2veg$rf2 <- tv$CoefTabs[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$CoefTabs[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

str(ch2soil)
str(ch2veg)

EXCL <- c("HWater", "SoilUnknown", "SoilWater", "Water", "HFor")
keeps <- rownames(ch2soil)[!(ch2soil$cr2 %in% EXCL) & !(ch2soil$rf2 %in% EXCL) ]
keepn <- rownames(ch2veg)[!(ch2veg$cr2 %in% EXCL) & !(ch2veg$rf2 %in% EXCL)]
trSoil <- trSoil[,keeps]
trVeg <- trVeg[,keepn]
ch2soil <- ch2soil[keeps,]
ch2veg <- ch2veg[keepn,]

rss <- rowSums(trSoil)
rss[rss==0] <- 1
trSoil <- trSoil / rss
rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn

summary(ch2soil)
summary(ch2veg)

stopifnot(all(rownames(trVeg) == rownames(kgrid)))
stopifnot(all(rownames(trSoil) == rownames(kgrid)))

## lat/long for centroids
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X+POINT_Y
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))

## pAspen: also part of kgrid
#xy <- spTransform(xy, proj4string(rpa))
#pAspen <- extract(rpa, xy)

## checking levels
compare_sets(ch2veg$cr2, get_levels()$veg)

compare_sets(ch2soil$cr2, get_levels()$soil)
setdiff(ch2soil$cr2, get_levels()$soil)
setdiff(get_levels()$soil, ch2soil$cr2)

(SPP <- get_all_species("mammals"))

(species <- SPP[3])
object <- load_spclim_data(species)


predict_mat2 <-
function(object, xy, veg, soil, ...)
{
    xy <- .tr_xy(xy)
    ## coefs in object are on log/logit scale, need linkinv
    fi <- if (object$taxon == "birds")
        poisson("log")$linkinv else binomial("logit")$linkinv
    if (missing(veg) && missing(soil))
        stop("veg or soil must be provided")
    xy <- spTransform(xy, proj4string(.read_raster_template()))
    if (!missing(veg)) {
        if (is.null(object$cveg)) {
            warning(sprintf("veg based estimates are unavailable for %s", object$species))
            Nveg <- NULL
        } else {
            if (nrow(veg) != nrow(coordinates(xy)))
                stop("nrow(veg) must equal number of points in xy")
            .check(as.factor(colnames(veg)), names(object$cveg))
            veg <- as(veg, "dgCMatrix") # sparse
            veg01 <- as(veg, "lMatrix") # sparse logical
            iveg <- extract(object$rveg, xy)
            ## 0 out the NAs
            iNAv <- is.na(iveg)
            iveg[iNAv] <- 0
            mveg <- object$cveg[match(colnames(veg), names(object$cveg))]
            mNAv <- is.na(mveg)
            mveg[mNAv] <- 0
            ## exploit spareness
            imatv <- iveg * veg01
            mmatv <- t(mveg * t(veg01))
            matv <- mmatv + imatv
            suppressMessages({
                matv[veg01] <- fi(matv[veg01])
            })
            Nveg <- matv * veg
            ## this blows, better to catch NAs and return or fail
            #Nveg[iNA,] <- NA
            #Nveg[,mNA] <- NA
            if (any(colnames(veg) == "SoftLin") && object$taxon == "birds" && object$version == "2017") {
                warning("veg contained SoftLin: check your assumptions")
                Nveg[,colnames(veg) == "SoftLin"] <- NA
            }
        }
    } else {
        Nveg <- NULL
    }
    if (!missing(soil)) {
        if (is.null(object$csoil)) {
            warning(sprintf("soil based estimates are unavailable for %s", object$species))
            Nsoil <- NULL
        } else {
            if (nrow(soil) != nrow(coordinates(xy)))
                stop("nrow(veg) must equal number of points in xy")
            .check(as.factor(colnames(soil)), names(object$csoil))
            soil <- as(soil, "dgCMatrix") # sparse
            soil01 <- as(soil, "lMatrix") # sparse logical
            isoil <- extract(object$rsoil, xy)
            rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))
            ipa <- extract(rpa, xy)
            ## 0 out the NAs
            iNAs <- is.na(isoil)
            isoil[iNAs] <- 0
            msoil <- object$csoil[match(colnames(soil), names(object$csoil))]
            mNAs <- is.na(msoil)
            msoil[mNAs] <- 0
            ## exploit spareness
            imats <- isoil * soil01
            mmats <- t(msoil * t(soil01))
            mats <- mmats + imats
            suppressMessages({
                mats[soil01] <- fi(mats[soil01])
            })
            Nsoil <- mats * soil
            #imats <- t(array(object$caspen * ipa + isoil, dim(soil), dimnames(soil)))
            #Nsoil <- fi(t(msoil + imats)) * soil
            if (any(colnames(soil) == "SoftLin") && object$taxon == "birds" && object$version == "2017") {
                warning("soil contained SoftLin: check your assumptions")
                Nsoil[,colnames(soil) == "SoftLin"] <- NA
            }
        }
    } else {
        Nsoil <- NULL
    }
    OUT <- list(
        veg=Nveg, veg_NA=list(row=unname(which(iNAv)), col=unname(which(mNAv))),
        soil=Nsoil, soil_NA=list(row=unname(which(iNAs)), col=unname(which(mNAs))))
    class(OUT) <- c("c4ippredmat")
    OUT
}

veg <- trVeg
colnames(veg) <- ch2veg$cr2

soil <- trSoil
colnames(soil) <- ch2soil$cr2

prm <- predict_mat2(object, xy, veg, soil)
str(prm)

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

## OF birds

library(DBI)

load("e:/peter/AB_data_v2018/data/raw/hvpoly/OFbirds-north.RData")
tab <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(tab) <- tab$sppid
ofs <- read.csv("~/repos/abmianalytics/lookup/OF-specificity.csv")
tab$ofspec <- ofs$ofspec[match(tab$AOU,ofs$species)]
tabn <- tab[tab$modelN & !is.na(tab$ofspec),]
boxplot(ofspec ~ oldforest, tabn)

SPP1 <- rownames(tabn)[tabn$oldforest == 1]
SPP2 <- rownames(tabn)[tabn$ofspec > 0.5]

load(system.file("extdata/raw_all.rda", package="cure4insect"))
MAX1 <- sapply(res, "[[", "max")[SPP1]
MAX2 <- sapply(res, "[[", "max")[SPP2]

of1 <- 0
for (spp in SPP1) {
    v <- OUT[[spp]]$veg
    v <- 100 * v / MAX1[spp]
    of1 <- of1 + v
}
#of1 <- 100 * of1 / max(of1)
of1 <- of1 / length(SPP1)

of2 <- 0
for (spp in SPP2) {
    v <- OUT[[spp]]$veg
    v <- 100 * v / MAX2[spp]
    of2 <- of2 + v
}
#of2 <- 100 * of2 / max(of2)
of2 <- of2 / length(SPP2)

ss <- sample(1:length(of1), 10^5)
plot(of1[ss], of2[ss],col="#0000FF20",pch=19,
    xlab="OF bird index (list)", ylab="OF bird index (specificity)")
abline(0,1,col=2)

fout <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
        "polygon-tool-pilot-results.sqlite")
db <- dbConnect(RSQLite::SQLite(), fout)
dbListTables(db)

tmp <- dbReadTable(db, "spden_north")
ofb <- data.frame(OBJECTID=tmp$OBJECTID, ofb_list=of1, ofb_spec=of2)
dbWriteTable(db, "ofbirds_north", ofb, overwrite = TRUE)

dbDisconnect(db)

tmp <- tabn[,c("species","oldforest","ofspec")]
tmp <- tmp[order(tmp$ofspec, decreasing=TRUE),]
write.csv(tmp,row.names=FALSE,file="OF-specificity.csv")

## responsibility --------------------------------------

#lt <- read.csv("~/repos/abmispecies/_data/birds.csv")
#rt <- read.csv("~/repos/abmianalytics/lookup/birds-global-responsibility.csv")
#rt$sppid <- lt$sppid[match(rt$Species, lt$species)]
#rt$AOU <- lt$AOU[match(rt$Species, lt$species)]
#write.csv(rt,row.names=FALSE,file="~/repos/abmianalytics/lookup/birds-global-responsibility.csv")

library(DBI)
library(cure4insect)
library(rgdal)
library(mefa4)
set_options(path = "w:/reports")
load_common_data()
#source("~/repos/abmianalytics/R/veghf_functions.R")

make_younger <- function(z) {
    for (i in rev(as.character(c(10, 20, 40, 60, 80, 100, 120, 140, 160)))) {
        z[grep(i, z)] <- gsub(i, "0", z[grep(i, z)])
    }
    z
}

rt <- read.csv("~/repos/abmianalytics/lookup/birds-global-responsibility.csv")
rownames(rt) <- rt$sppid
load(system.file("extdata/raw_all.rda", package="cure4insect"))

veg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
veg_mapping <- cbind(V6=as.character(veg$VEGHFAGE_FINE), V5=as.character(veg$PolyReclass))
soil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2014.csv")
soil_mapping <- cbind(In=as.character(soil$SOILHF_FINE), Out=as.character(soil$PolyReclass))

f <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")
Taxa <- "birds"
PilotAreas <- c("south", "north")
#PilotArea <- "south"
for (PilotArea in PilotAreas) {

    db <- dbConnect(RSQLite::SQLite(), f)
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
    x$VEGHFAGEclass2 <- reclass(x$VEGHFAGEclass, veg_mapping, all=TRUE)
    x$SOILclass2 <- reclass(x$SOILclass, soil_mapping, all=TRUE)
    x$SOILHFclass2 <- reclass(x$SOILHFclass, soil_mapping, all=TRUE)

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

    for (taxon in Taxa) {
        Species <- get_all_species(taxon=taxon)
        SPPS <- get_all_species(taxon, mregion="south")
        SPPN <- get_all_species(taxon, mregion="north")

        N <- length(Species)
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

        OUT <- 0
        n <- 1
        #species <- Species[1]
        for (species in Species) {
            cat("Responsibility", PilotArea, taxon, species, "---", n, "/", N, "\n")
            flush.console()
            object <- load_spclim_data(species)
            pred_curr <- suppressWarnings(predict(object, xy=xy,
                veg=x$VEGHFAGEclass3, soil=x$SOILHFclass3))
            pred_curr[is.na(pred_curr)] <- 0 # water and non-veg
            w <- rt[species, "Weighting"] # needs to sum to 1 over Species
            MAX <- sapply(res, "[[", "max")[species]

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

            ## provincial 1km^2 max might be lower than poly level max
            OUT <- OUT + (w * 100 * pmin(cr, MAX) / MAX)

            n <- n + 1
        }
        save(OUT, file=file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
            paste0("birds-responsibility-", PilotArea, ".RData")))
    }
}

library(DBI)
load(file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    paste0("peter-results.RData")))
fout <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
        "polygon-tool-pilot-results.sqlite")
db <- dbConnect(RSQLite::SQLite(), fout)
dbListTables(db)

e <- new.env()
load(file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    paste0("birds-responsibility-north.RData")), envir=e)
tmp <- dbReadTable(db, "spden_north")
df <- data.frame(OBJECTID=tmp$OBJECTID, resp_birds=e$OUT)
datn$resp_birds <- e$OUT
dbWriteTable(db, "resp_birds_north", df, overwrite = TRUE)

e <- new.env()
load(file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    paste0("birds-responsibility-south.RData")), envir=e)
tmp <- dbReadTable(db, "spden_south")
df <- data.frame(OBJECTID=tmp$OBJECTID, resp_birds=e$OUT)
dats$resp_birds <- e$OUT
dbWriteTable(db, "resp_birds_south", df, overwrite = TRUE)

dbDisconnect(db)

save(datn, dats, file=file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    paste0("peter-results.RData")))


## plot results ------------------------------------------------

library(DBI)
library(cure4insect)
library(rgdal)
library(mefa4)
f <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "Backfilled100kmtestarea","polygon-tool-pilot.sqlite")
rt <- .read_raster_template()
#rt200 <- disaggregate(rt, 5)

f2 <- file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "polygon-tool-pilot-results.sqlite")
db <- dbConnect(RSQLite::SQLite(), f2)
zs <- dbReadTable(db, "spden_south")
zn <- dbReadTable(db, "spden_north")
zb <- dbReadTable(db, "ofbirds_north")
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
xyn <- spTransform(xyn, proj4string(rt))
zs <- zs[,-1]
zn <- zn[,-1]
zb <- zb[,-1]

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

lb <- list()
for (i in colnames(zb)) {
    rrN <- trim(rasterize(xyn, rt, zb[,i] * xn$Shape_Area, fun=sum), values=NA)
    rrA <- trim(rasterize(xyn, rt, xn$Shape_Area, fun=sum), values=NA)
    rr <- rrN/rrA
    values(rr)[!is.na(values(rr)) & values(rr) == max(values(rr), na.rm=TRUE)] <- 100
    lb[[i]] <- rr
}
lb <- stack(lb)

pdf("~/GoogleWork/abmi/hvpoly/results/ofbirds-north.pdf", width=14)
plot(lb, colNA="grey", col=viridis::viridis(101), axes=FALSE, box=FALSE)
dev.off()

datn <- data.frame(OBJECTID=xn$OBJECTID, zn, zb)
datn$ofb_list <- as.integer(round(datn$ofb_list))
datn$ofb_spec <- as.integer(round(datn$ofb_spec))
dats <- data.frame(OBJECTID=xs$OBJECTID, zs)

save(datn, dats, file=file.path("e:/peter", "AB_data_v2018", "data", "raw", "hvpoly",
    "peter-results.RData"))



