## Packages
library(raster)
library(sp)
library(rgdal)
library(mefa4)
library(rasterVis)

ROOT <- "c:/p"
VER <- "AB_data_v2015"

tveg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tsoil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid.Rdata"))
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

rt <- raster(file.path(ROOT, VER, "data", "kgrid", "AHM1k.asc"))
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
projection(rt) <- crs
mat0 <- as.matrix(rt)

## HF to plot
hf1 <- (dd1km_pred$veg_current[,!is.na(tveg$HF)]) / rowSums(dd1km_pred$veg_current)
hft <- tveg[!is.na(tveg$HF),]
hf2 <- cbind(TotalHF=rowSums(hf1), as.matrix(groupSums(hf1, 2, hft$UseInAnalysis)))
hf1 <- as.matrix(groupSums(hf1, 2, hft$HF))

## raster stack from hf

r_water <- as_Raster(kgrid$Row, kgrid$Col, kgrid$pWater, rt)
r_water[r_water <= 0.99] <- NA
r_mask <- as_Raster(kgrid$Row, kgrid$Col, 
    ifelse(kgrid$NRNAME == "Rocky Mountain" & kgrid$X < 800000, 1, 0), rt)
r_mask[r_mask < 1] <- NA
r_overlay <- r_water
r_overlay[!is.na(r_water)] <- 1
r_overlay[!is.na(r_mask)] <- 2
r_overlay <- as.factor(r_overlay)
rat <- levels(r_overlay)[[1]]
rat$levels <- c("Water", "Rockies")
rat$code <- 1:2
levels(r_overlay) <- rat

nr <- as.matrix(Xtab(as.integer(NRNAME) ~ Row + Col, kgrid))
nr[is.na(mat0)] <- NA
nr <- as.factor(raster(x=nr, template=rt))
rat <- levels(nr)[[1]]
rat$NaturalRegion <- levels(kgrid$NRNAME)
rat$code <- seq_len(nlevels(kgrid$NRNAME))
levels(nr) <- rat

city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
city <- spTransform(city, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))


#x <- 100 * hf2[,1]

#r_hf <- as_Raster(kgrid$Row, kgrid$Col, rowSums(hf1), rt)
r_hf <- as_Raster(kgrid$Row, kgrid$Col, hf2, rt, verbose=1)


pdf(file.path(ROOT, VER, "out", "figs", "hf", "hf-use-in-analysis.pdf"),
    width=10, height=8)
levelplot(r_hf, margin=FALSE, par.settings=BuRdTheme(), 
    scales=list(draw=FALSE), contour=FALSE, main="Human Footprint Types")
dev.off()

pdf(file.path(ROOT, VER, "out", "figs", "hf", "hf-1file.pdf"),
    width=8, height=12, onefile=TRUE)
for (i in 1:ncol(hf2)) {
    cat(i, "\n");flush.console()
    x <- 100 * hf2[,i]
    map_fun(x, q=0.999, main=colnames(hf2)[i])
}
dev.off()


## xy has lat/long in it -- reprojected into UTM
xy <- kgrid[sample.int(nrow(kgrid), 200), c("POINT_X","POINT_Y")]
xy_map(xy, pch=19, cex=0.8, col=1:4)


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


