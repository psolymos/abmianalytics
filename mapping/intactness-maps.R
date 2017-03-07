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

## intactness non birds

xx <- c(#birds="Birds_Species Intactness",
    lichens="Lichens_Species Intactness",
    mosses="Mosses_Species Intactness",
    mammals="Mammals_Species Intactness",
    mites="Mites_Species Intactness",
    vplants="VPlants_Species Intactness")

nval_list <- list()

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

nval <- pmin(100, ceiling(99 * meanSI/100)+1)
nval_list[[zz]] <- nval

if (FALSE) {
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

}
names(nval_list) <- names(xx)

## intactness birds
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

nval <- pmin(100, ceiling(99 * meanSI)+1)
nval_list[["birds"]] <- nval

OF_si <- rowMeans(res[,tab[use,"oldforest"]==1], na.rm=TRUE)
OF_si[is.na(OF_si)] <- 0
OF_si <- pmin(100, ceiling(99 * OF_si)+1)


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

## intactness, all species: average of taxa (no weighting)

cnval <- do.call(cbind, nval_list)
AllSI <- rowMeans(cnval)
AllSI <- pmin(100, ceiling(99 * AllSI/100)+1)
out <- data.frame(LinkID=kgrid$Row_Col, cnval, all=AllSI, oldforestbirds=OF_si)
write.csv(out, row.names=FALSE, file="w:/normalized_data/intactness-normalized.csv")


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

xxx <- xx
for (i in 2:8)
    xxx[,i] <- pmin(100, ceiling(99 * xxx[,i])+1)
write.csv(xxx, row.names=FALSE, file="w:/normalized_data/richness-normalized.csv")

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

## normalize uniqueness


xx <- list(all="e:/peter/sppweb2015-round3/uniqueness/All Taxa Uniqueness Maps and Figures/All Taxa uniqueness North and South_native habitat based April 2016.Rdata",
birds="e:/peter/sppweb2015-round3/uniqueness/Birds Uniqueness Maps and Figures/Birds uniqueness North and South_native habitat based April 2016.Rdata",
habitats="e:/peter/sppweb2015-round3/uniqueness/Habitat Uniqueness Maps and Figures/Habitat uniqueness North and South_April 2016.Rdata",
lichens="e:/peter/sppweb2015-round3/uniqueness/Lichen Uniqueness Maps and Figures/Lichens uniqueness North and South_native habitat based April 2016.Rdata",
mammals="e:/peter/sppweb2015-round3/uniqueness/Mammals Uniqueness Maps and Figures/Mammals uniqueness North and South_native habitat based April 2016.Rdata",
mites="e:/peter/sppweb2015-round3/uniqueness/Mites Uniqueness Maps and Figures/Mites uniqueness North and South_native habitat based April 2016.Rdata",
mosses="e:/peter/sppweb2015-round3/uniqueness/Moss Uniqueness Maps and Figures/Moss uniqueness North and South_native habitat based April 2016.Rdata",
vplants="e:/peter/sppweb2015-round3/uniqueness/VPlants Uniqueness Maps and Figures/Vascular plants uniqueness North and South_native habitat based April 2016.Rdata")


out <- data.frame(LinkID=kgrid$Row_Col)
fff <- function(x) pmin(100, ceiling(99 * x)+1)
for (i in 1:8) {

load(xx[[i]])
North_out <- North_out[match(rownames(kgrid), North_out$LinkID),]
South_out <- South_out[match(rownames(kgrid), South_out$LinkID),]
out[[paste0(names(xx)[i], "_uniq_north")]] <- fff(North_out[,2])
out[[paste0(names(xx)[i], "_uniq_south")]] <- fff(South_out[,2])

}

write.csv(out, row.names=FALSE, file="w:/normalized_data/uniqueness-normalized.csv")

## new color scale for intactness

## Colour gradient for intactness map
si1<-c("#D73027","#FC8D59","#FEE090","#E0F3D8","#41DB45","#12A412")
## Function to interpolate among these colours for intactness map
si2<-colorRampPalette(si1, space = "rgb")

ff <- list("Intactness - Birds"="w:/multispecies/intactness/intactness-birds.asc",
    "Intactness - All Species"="w:/multispecies/intactness/intactness-all.asc",
    "Intactness - Vasculer Plants"="w:/multispecies/intactness/intactness-vplants.asc",
    "Intactness - Soil Mites"="w:/multispecies/intactness/intactness-mites.asc",
    "Intactness - Mammals (Snow Tracking)"="w:/multispecies/intactness/intactness-mammals.asc",
    "Intactness - Bryophytes"="w:/multispecies/intactness/intactness-mosses.asc",
    "Intactness - Lichens"="w:/multispecies/intactness/intactness-lichens.asc")

nf <- list("Intactness - Birds"="birds",
    "Intactness - All Species"="all",
    "Intactness - Vasculer Plants"="vplants",
    "Intactness - Soil Mites"="mites",
    "Intactness - Mammals (Snow Tracking)"="mammals",
    "Intactness - Bryophytes"="mosses",
    "Intactness - Lichens"="lichens")

i <- 1
for (i in 1:7) {
r <- raster(ff[[i]])

png(file.path("v:/contents/2016/multispecies/intactness", paste0("intactness-", nf[[i]], ".png")), width=W, height=H)
op <- par(mar=c(0, 0, 4, 0) + 0.1)
map_fun(r, main=names(ff)[i], colScale=si2)
par(op)
dev.off()
}
