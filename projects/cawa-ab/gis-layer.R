load("e:/peter/AB_data_v2016/out/birds/results/cawa/cawa-km-predB.Rdata")
load("e:/peter/AB_data_v2016/out/kgrid/kgrid_table.Rdata")

kg <- kgrid[rownames(km),]



## Packages
library(raster)
library(sp)
library(rgdal)
library(mefa4)
#library(rasterVis)
library(utils)
source("~/repos/abmianalytics/R/maps_functions.R")


## raster template to use
rt <- raster("e:/peter/AB_data_v2016/data/kgrid/AHM1k.asc")
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
projection(rt) <- crs
mat0 <- as.matrix(rt)

x <- kgrid$Row
y <- kgrid$Col
z <- km[match(rownames(kgrid), rownames(km)), "Median"]
mat <- as.matrix(Xtab(z ~ x + y))
mat[is.na(mat0)] <- NA
dim(mat) <- dim(mat0)
r2 <- raster(x=mat, template=rt)
writeRaster(r2, "e:/peter/AB_data_v2016/out/birds/results/cawa/cawa-km-predB.tif")

km$Run1 <- km$Aland <- NULL
km <- round(km, 6)
km$POINT_X <- kg$POINT_X
km$POINT_Y <- kg$POINT_Y
write.csv(km, row.names=FALSE, file=
    "e:/peter/AB_data_v2016/out/birds/results/cawa/cawa-km-predB.csv")
