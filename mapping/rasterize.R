library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(raster)

## grid points
load("kgrid_table_km.Rdata")
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
xy <- kgrid
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, crs)

## template raster
rt <- raster("AHM1k.asc")
projection(rt) <- crs
mat0 <- as.matrix(rt)

## load some species data
spp <- "ALFL"
load(file.path("km2", paste0(spp, ".Rdata")))
stopifnot(all(rownames(km2) == rownames(kgrid)))

## abundance to raster
val <- as.numeric(km2[,"Ref"]) # use "Curr" for current abundance
mat <- as.matrix(Xtab(val ~ Row + Col, kgrid))
mat[is.na(mat0)] <- NA
r <- raster(x=mat, template=rt)
