## Packages
library(raster)
library(sp)
library(rgdal)
library(mefa4)
library(rasterVis)

## root directory and version
ROOT <- "c:/p"
VER <- "AB_data_v2015"

## source functions
source("~/repos/abmianalytics/R/maps_functions.R")

## lookup tables for veg and soil classes (combined with hf)
tveg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tsoil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

## cell x vag/soil matrices and xy lookup table (climate, region, etc)
load(file.path(ROOT, VER, "out", "kgrid", "veg-hf_1kmgrid.Rdata"))
load(file.path(ROOT, VER, "out", "kgrid", "kgrid_table.Rdata"))

## raster template to use
rt <- raster(file.path(ROOT, VER, "data", "kgrid", "AHM1k.asc"))
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

## ------------- common stuff ends here -----------------

## maps: HF

hf1 <- (dd1km_pred$veg_current[,!is.na(tveg$HF)]) / rowSums(dd1km_pred$veg_current)
hft <- tveg[!is.na(tveg$HF),]
hft <- hft[colnames(hf1),]
hf2 <- cbind(TotalHF=rowSums(hf1), as.matrix(groupSums(hf1, 2, hft$UseInAnalysis)))
hf1 <- as.matrix(groupSums(hf1, 2, hft$HF))

#r_hf <- as_Raster(kgrid$Row, kgrid$Col, hf2, rt, verbose=1)
#pdf(file.path(ROOT, VER, "out", "figs", "hf", "hf-use-in-analysis.pdf"),
#    width=10, height=8)
#levelplot(r_hf, margin=FALSE, par.settings=BuRdTheme(), 
#    scales=list(draw=FALSE), contour=FALSE, main="Human Footprint Types")
#dev.off()

pdf(file.path(ROOT, VER, "out", "figs", "hf", "hf-1file.pdf"),
    width=8, height=12, onefile=TRUE)
for (i in 1:ncol(hf2)) {
    cat(i, "/", ncol(hf2), "\n");flush.console()
    x <- 100 * hf2[,i]
    map_fun(x, q=0.999, main=colnames(hf2)[i], 
        colScale="hf", maskRockies=FALSE)
}
dev.off()

## maps: soils

soil1 <- (dd1km_pred$soil_reference) / rowSums(dd1km_pred$soil_reference)
st <- tsoil[colnames(soil1),]
soil1 <- as.matrix(groupSums(soil1, 2, st$Levels1))
pSoil <- 1 - soil1[,"SoilUnknown"]
#pSoil <- (rowSums(dd1km_pred$soil_reference) - dd1km_pred$soil_reference[,"UNK"]) /
#    rowSums(dd1km_pred$soil_reference)
soil1 <- soil1[,colnames(soil1) != "SoilUnknown"]

## mask soil=UNK
r_mask2 <- as_Raster(kgrid$Row, kgrid$Col, 
    ifelse(pSoil < 0.99, 1, 0), rt)
r_mask2[r_mask2 < 1] <- NA

## replace the overlay object
r_overlay <- r_water
r_overlay[!is.na(r_water)] <- 1
r_overlay[!is.na(r_mask2)] <- 2
r_overlay <- as.factor(r_overlay)
rat <- levels(r_overlay)[[1]]
rat$levels <- c("Water", "NoSoilInfo")
rat$code <- 1:2
levels(r_overlay) <- rat

source("~/repos/abmianalytics/R/maps_functions.R")

pdf(file.path(ROOT, VER, "out", "figs", "soil", "soil-1file.pdf"),
    width=8, height=12, onefile=TRUE)
for (i in 1:ncol(soil1)) {
    cat(i, "/", ncol(soil1), "\n");flush.console()
    map_fun(100*soil1[,i], q=0.999, main=colnames(soil1)[i], 
        colScale="soil", maskRockies=TRUE, plotWater=TRUE)
}
dev.off()

## maps: veg+age

veg1 <- (dd1km_pred$veg_reference) / rowSums(dd1km_pred$veg_reference)
vt <- tveg[colnames(veg1),]
veg2 <- as.matrix(groupSums(veg1, 2, vt$AGE))
veg2 <- veg2[,colnames(veg2) != ""]
colnames(veg2) <- paste0("Forest age: ", 
    as.character(vt$AGE_Description)[match(colnames(veg2), vt$AGE)])
veg3 <- as.matrix(groupSums(veg1, 2, vt$Broad))
veg3 <- veg3[,colnames(veg3) != "NonVeg"]
veg3 <- veg3[,colnames(veg3) != "Water"]
veg4 <- as.matrix(groupSums(veg1, 2, vt$Type))
veg4[,"NonVeg"] <- veg4[,"NonVeg"] + veg4[,"Bare"]
veg4 <- veg4[,colnames(veg4) != "Bare"]
veg4 <- veg4[,colnames(veg4) != "Water"]
veg4 <- veg4[,setdiff(colnames(veg4), colnames(veg3))]
colnames(veg4) <- paste0(colnames(veg4), " (up+low)")
veg5 <- as.matrix(groupSums(veg1, 2, vt$VEG))
veg5 <- veg5[,!(colnames(veg5) %in% c("NonVeg","GrassHerb","Shrub",
    "Water","Wetland-Bare"))]

vegAll <- cbind(veg3, veg5, veg4, veg2)

pdf(file.path(ROOT, VER, "out", "figs", "veg", "veg-1file.pdf"),
    width=8, height=12, onefile=TRUE)
for (i in 1:ncol(vegAll)) {
    cat(i, "/", ncol(vegAll), colnames(vegAll)[i], "\n");flush.console()
    map_fun(100*vegAll[,i], q=0.999, main=colnames(vegAll)[i], 
        colScale="terrain", maskRockies=FALSE, plotWater=TRUE)
}
dev.off()

## this is how detections can be mapped
## xy has lat/long in it -- reprojected into UTM

## mites
spp1 <- "AchipteriaColeoptrata"
x <- read.csv(file.path(ROOT, VER, "out", "species",
    "OUT_Mites_Species_Site-Binomial_2015-05-22.csv"))
x <- x[!is.na(x$PUBLIC_X),]
y <- as.matrix(x[,which(colnames(x) == spp1):ncol(x)])
xy <- x[,c("PUBLIC_X","PUBLIC_Y")]
colnames(xy) <- c("POINT_X", "POINT_Y")
#xy <- kgrid[sample.int(nrow(kgrid), 200), c("POINT_X","POINT_Y")]

#i <- 1
pdf(file.path(ROOT, VER, "out", "figs", "mites-pa.pdf"),
    width=6, height=9, onefile=TRUE)
for (i in 1:ncol(y)) {
    cat(i, "/", ncol(y), "\n");flush.console()
    xy_map(xy, y[,i], main=colnames(y)[i])
}
dev.off()



## ------- tests

mat <- as.matrix(Xtab(MAT ~ Row + Col, kgrid))
mat[is.na(mat0)] <- NA
rr <- raster(x=mat, template=rt)

plot(rr, axes=FALSE, box=FALSE, legend=TRUE, main="Title") 

levelplot(rr, margin=FALSE, par.settings=PuOrTheme(), 
    scales=list(draw=FALSE), contour=FALSE, main="Title")

rr2 <- stack(rr, rr^2)
names(rr2) <- c("a", "b")
levelplot(rr2)
levelplot(rr2, margin=FALSE, par.settings=PuOrTheme(), 
    scales=list(draw=FALSE), contour=FALSE, main="Title")
## see levelplot examples for putting points on a map
pts <- sampleRandom(rr, size=20, sp=TRUE)
## Using +.trellis and layer from latticeExtra
levelplot(rr, margin=FALSE, par.settings=BTCTheme, 
    scales=list(draw=FALSE), contour=FALSE, 
    main="Title") + layer(sp.points(pts, col = 'red'))

colTheme <- rasterTheme(region=brewer.pal('Blues', n=9))

nr <- as.matrix(Xtab(as.integer(NRNAME) ~ Row + Col, kgrid))
nr[is.na(mat0)] <- NA
nr <- as.factor(raster(x=nr, template=rt))
rat <- levels(nr)[[1]]
rat$NaturalRegion <- levels(kgrid$NRNAME)
rat$code <- seq_len(nlevels(kgrid$NRNAME))
levels(nr) <- rat

levelplot(nr, margin=FALSE, #par.settings=PuOrTheme(), 
    scales=list(draw=FALSE), contour=FALSE, main="Title",
    col.regions=c('palegreen', 'midnightblue', 'indianred1', 'green', 'blue', 
    'red')) + layer(sp.points(pts, col = 'red'))

plot(nr, axes=FALSE, box=FALSE, legend=TRUE, add=TRUE, alpha=0.2) 

x0 <- kgrid[sample.int(nrow(kgrid), 200), c("X","Y")]
coordinates(x0) <- ~ X + Y
proj4string(x0) <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")


