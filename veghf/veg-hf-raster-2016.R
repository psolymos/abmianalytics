library(mefa4)
library(raster)
library(rgdal)


load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])

load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

ch2veg$g <- tv$UseInPIX[match(ch2veg$cr, rownames(tv))]
ch2veg$f <- ifelse(startsWith(as.character(ch2veg$cr), "CC"), "CC", "not")
ch2veg$s <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

AA <- groupSums(trVeg, 2, ch2veg[colnames(trVeg), "g"])

cn <- c("Bare"="Open",
    "HWater"="Water",
    "Agr"="Agr",
    "UrbInd"="UrbInd",
    "SoftLin"="SoftLin",
    "Roads"="Roads",
    "Decid"="Decid",
    "DecidO"="Decid",
    "GrFen"="OpenWet",
    "Grass"="Open",
    "Marsh"="OpenWet",
    "Mixed"="Decid",
    "MixedO"="Decid",
    "Pine"="Conif",
    "PineO"="Conif",
    "Shrub"="Open",
    "BSpr"="ConifWet",
    "TrFen"="ConifWet",
    "Swamp"="OpenWet",
    "SnowIce"="Open",
    "WSpr"="Conif",
    "WSprO"="Conif",
    "BSprO"="ConifWet",
    "Water"="Water")

AA <- groupSums(AA[,names(cn)], 2, cn)

rs <- rowSums(AA)
AA <- AA / rs

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

rr <- list()
for (i in colnames(AA)) {
    cat(i, "\n")
    flush.console()
    rr[[i]] <- make_raster(AA[,i], kgrid, rt)
}

rr <- stack(rr)

#ply <- readOGR(dsn=system.file("extdata/OSA_bound.geojson", package="cure4insect"))
#ply <- spTransform(ply, proj4string(rt))
#rr <- crop(rr, extent(ply))
#rr <- mask(rr, ply)
rr2 <- crop(rr, extent(c(350000, 880000, 6030000, 6420000)))
rr2 <- trim(rr2, values = NA)

plot(rr2)
writeRaster(rr2, overwrite=TRUE,
    file="~/Dropbox/courses/aos-2019-anchorage/data/landcover-hfi2016.grd")

#aa=stack("~/Dropbox/courses/aos-2019-anchorage/data/pif-level-vaghf2016.grd")

