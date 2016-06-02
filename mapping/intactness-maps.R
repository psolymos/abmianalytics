## using veg based models only

library(mefa4)
library(raster)
library(sp)
library(rgdal)

ROOT <- "c:/p/AB_data_v2015"
ROOT2 <- "e:/peter/sppweb2015-round3/Species Intactness"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))

source("~/repos/abmianalytics/R/maps_functions.R")

## raster template to use
rt <- raster(file.path(ROOT, "data", "kgrid", "AHM1k.asc"))
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
projection(rt) <- crs
mat0 <- as.matrix(rt)

## raster layers that are used for mapping
## mostly water cells
r_water <- as_Raster(kgrid$Row, kgrid$Col, kgrid$pWater, rt)
r_water[r_water <= 0.99] <- NA
## Rockies to mask out
r_mask <- as_Raster(kgrid$Row, kgrid$Col, 
    ifelse(kgrid$NRNAME == "Rocky Mountain" & kgrid$X < 800000, 1, 0), rt)
r_mask[r_mask < 1] <- NA
## combine the 2 (this is necessary due to a bug in raster plotting
## function: when >3 layers are shown there is a misterious mismatch...)
r_overlay <- r_water
r_overlay[!is.na(r_water)] <- 1
r_overlay[!is.na(r_mask)] <- 2
r_overlay <- as.factor(r_overlay)
rat <- levels(r_overlay)[[1]]
rat$levels <- c("Water", "Rockies")
rat$code <- 1:2
levels(r_overlay) <- rat
## natural regions
nr <- as.matrix(Xtab(as.integer(NRNAME) ~ Row + Col, kgrid))
nr[is.na(mat0)] <- NA
nr <- as.factor(raster(x=nr, template=rt))
rat <- levels(nr)[[1]]
rat$NaturalRegion <- levels(kgrid$NRNAME)
rat$code <- seq_len(nlevels(kgrid$NRNAME))
levels(nr) <- rat
## city coordinates
city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
city <- spTransform(city, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

H <- 1000 
W <- 600
kgrid$Row <- as.factor(kgrid$Row)
kgrid$Col <- as.factor(kgrid$Col)

## non birds

xx <- c(#birds="Birds_Species Intactness",
    lichens="Lichens_Species Intactness",
    mosses="Mosses_Species Intactness",
    mammals="Mammals_Species Intactness",
    mites="Mites_Species Intactness",
    vplants="VPlants_Species Intactness")

for (zz in 1:length(xx)) {
DIR <- xx[zz] # "Birds_Species Intactness"
TAXA <- names(xx)[zz] # "birds"

fl <- list.files(file.path(ROOT2, DIR))
names(fl) <- sapply(strsplit(fl, " "), "[[", 1)
tab <- read.csv(paste0("c:/Users/Peter/repos/abmispecies/_data/", TAXA, ".csv"))
rownames(tab) <- tab$sppid
if (is.null(tab$nonnative))
    tab$nonnative <- FALSE
use <- rownames(tab)[!tab$nonnative & tab$map.pred]

res <- matrix(NA, nrow(kgrid), length(use))
rownames(res) <- rownames(kgrid)
colnames(res) <- use
for (i in 1:length(use)) {
    cat(i, "/", length(use), "-", TAXA, use[i], "\n")
    flush.console()
    e <- new.env()
    load(file.path(ROOT2, DIR, fl[i]), envir=e)
    SI <- e$SI
    res[,i] <- SI$SI[match(rownames(kgrid), SI$LinkID)]
}

SI <- res
meanSI <- rowMeans(SI, na.rm=TRUE)
meanSI[is.na(meanSI)] <- 100

## make raster
png(file.path(ROOT2, "out", paste0("intactness-", TAXA, ".png")), width=W, height=H)
op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(meanSI, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
xxx <- switch(TAXA,
    birds="Birds",
    lichens="Lichens",
    mosses="Btyophytes",
    mammals="Mammals",
    mites="Mites",
    vplants="Vascular Plants (Native)")
title(main=paste0("Intactness - ", xxx))
par(op)
dev.off()
writeRaster(r, file.path(ROOT2, "out", paste0("intactness-", TAXA, ".asc")),
    overwrite=TRUE)

}


## birds
TAXA <- "birds"
fl <- list.files("e:/peter/sppweb2015/birds-si")
names(fl) <- sapply(strsplit(fl, "."), "[[", 1)
tab <- read.csv(paste0("c:/Users/Peter/repos/abmispecies/_data/", TAXA, ".csv"))
rownames(tab) <- tab$sppid
if (is.null(tab$nonnative))
    tab$nonnative <- FALSE
use <- rownames(tab)[!tab$nonnative & tab$map.pred]

res <- matrix(NA, nrow(kgrid), length(use))
rownames(res) <- rownames(kgrid)
colnames(res) <- use
for (i in 1:length(use)) {
    cat(i, "/", length(use), "-", TAXA, use[i], "\n")
    flush.console()
    SI <- read.csv(file.path("e:/peter/sppweb2015/birds-si", fl[i]))
    res[,i] <- SI$SI[match(rownames(kgrid), SI$LinkID)]
}

SI <- res
meanSI <- rowMeans(SI, na.rm=TRUE)
#meanSI <- rowMeans(SI)
meanSI[is.na(meanSI)] <- 0

## make raster
png(file.path(ROOT2, "out", paste0("intactness-", TAXA, ".png")), width=W, height=H)
#op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(meanSI, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
title(main="Intactness - Birds")
#par(op)
dev.off()
writeRaster(r, file.path(ROOT2, "out", paste0("intactness-", TAXA, ".asc")),
    overwrite=TRUE)




## All species: average of taxa (no weighting)

ROOT2 <- "e:/peter/sppweb2015-round3/Species Intactness"

x1 <- raster(file.path(ROOT2, "out", "intactness-birds.asc"))
x1 <- x1*100
x2 <- raster(file.path(ROOT2, "out", "intactness-mammals.asc"))
x3 <- raster(file.path(ROOT2, "out", "intactness-mites.asc"))
x4 <- raster(file.path(ROOT2, "out", "intactness-mosses.asc"))
x5 <- raster(file.path(ROOT2, "out", "intactness-vplants.asc"))
x6 <- raster(file.path(ROOT2, "out", "intactness-lichens.asc"))

x <- (x1 + x2 + x3 + x4 + x5 + x6) / 6
#plot(x)

## make raster
png(file.path(ROOT2, "out", paste0("intactness-all.png")), width=W, height=H)
#op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(x, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
title(main="Intactness - All species")
#par(op)
dev.off()
writeRaster(r, file.path(ROOT2, "out", paste0("intactness-all.asc")),
    overwrite=TRUE)

png(file.path(ROOT2, "out", paste0("intactness-birds.png")), width=W, height=H)
#op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(x1, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
title(main="Intactness - Birds")
#par(op)
dev.off()
writeRaster(r, file.path(ROOT2, "out", paste0("intactness-birds.asc")),
    overwrite=TRUE)

## uniqueness -----------------------

ROOT3 <- "e:/peter/sppweb2015-round3/fromJim/Uniqueness Maps Peter"
fl <- c("Bird uniqueness probability based_Only native habitat February 2016.Rdata",
    "Lichens uniqueness_native habitat based February 2016.Rdata",
    "Mammals uniqueness_native habitat based February 2016.Rdata",
    "Mites uniqueness February 2016_native habitat based February 2016.Rdata",
    "Moss uniqueness_native habitat based February 2016.Rdata",
    "Vascular plants uniqueness_native habitat based February 2016.Rdata")

i <- 1
for (i in 1:6) {
lab <- c("birds", "lichens", "mammals", "mites", "mosses", "vplants")[i]
Lab <- c("Birds", "Lichens", "Mammals", "Mites", "Bryophytes", "Vascular plants")[i]
load(file.path(ROOT3, fl[i]))
xx <- Alberta_Uniq[,2][match(rownames(kgrid), Alberta_Uniq[,1])]

png(file.path(ROOT3, "out", paste0("uniqueness-PSversion-", lab, ".png")), width=W, height=H)
#op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(xx, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
title(main=paste0("Uniqueness - ", Lab))
#par(op)
dev.off()
writeRaster(r, file.path(ROOT3, paste0("uniqueness-", lab, ".asc")),
    overwrite=TRUE)
}

x1 <- raster(file.path(ROOT3, "uniqueness-birds.asc"))
x2 <- raster(file.path(ROOT3, "uniqueness-mammals.asc"))
x3 <- raster(file.path(ROOT3, "uniqueness-mites.asc"))
x4 <- raster(file.path(ROOT3, "uniqueness-mosses.asc"))
x5 <- raster(file.path(ROOT3, "uniqueness-vplants.asc"))
x6 <- raster(file.path(ROOT3, "uniqueness-lichens.asc"))

x <- (x1 + x2 + x3 + x4 + x5 + x6) / 6
png(file.path(ROOT3, paste0("uniqueness-PSversion-all.png")), width=W, height=H)
#op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(x, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
title(main=paste0("Uniqueness - All species"))
#par(op)
dev.off()
writeRaster(r, file.path(ROOT3, paste0("uniqueness-all.asc")),
    overwrite=TRUE)


## richness rasters

x <- read.csv(file.path(ROOT3, "Relative Species richness for all taxa March 2016.csv"))

xx <- x[match(rownames(kgrid), x$LinkID),]

for (i in 1:7) {
    lab <- c("mammals", "birds", "mites", "vplants", "mosses", "lichens", "all")[i]
    z <- xx[,i+1]
    r <- map_fun(z, q=1, main="", colScale="intactness",
        plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
        legend=TRUE,
        mask=NULL)
    writeRaster(r, file.path(ROOT3, paste0("richness-", lab, ".asc")),
        overwrite=TRUE)
}


r1 <- map_fun(xx[,"Curr"], q=1, main="", colScale="abund",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
r2 <- map_fun(xx[,"Ref"], q=1, main="", colScale="abund",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
writeRaster(r1, "e:/peter/sppweb2015-round3/fromJim/NN plant/nNNSpp-Curr.asc",
    overwrite=TRUE)
writeRaster(r2, "e:/peter/sppweb2015-round3/fromJim/NN plant/nNNSpp-Ref.asc",
    overwrite=TRUE)