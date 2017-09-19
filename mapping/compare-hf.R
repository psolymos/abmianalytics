## comparing HF inventories
library(mefa4)

f14v1 <- "e:/peter/AB_data_v2017/data/analysis/veg-hf_1kmgrid_v6-fixage0.Rdata"
f10 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2010-grid_v6hf2010_coarse-fixage0.Rdata"
f14v2 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_coarse-fixage0.Rdata"

e <- new.env()
load(f14v1, envir=e)
x14v1 <- e$dd1km_pred$veg_current
e <- new.env()
load(f14v2, envir=e)
x14v2 <- e$dd1km_pred$veg_current
e <- new.env()
load(f10, envir=e)
x10 <- e$dd1km_pred$veg_current
rm(e)


library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(gstat)
library(raster)
#library(viridis)
load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))

xy <- kgrid
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

rt <- raster(file.path("e:/peter/AB_data_v2016", "data", "kgrid", "AHM1k.asc"))
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
projection(rt) <- crs
xy <- spTransform(xy, crs)
mat0 <- as.matrix(rt)
matR <- as.matrix(Xtab(ifelse(kgrid$NRNAME=="Rocky Mountain", 1, 0) ~ Row + Col, kgrid))


od <- setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(rt))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province
setwd(od)

colDiv <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B",
    "#FFFFBF","#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
colSeq <- rev(viridis::magma(100))

x10 <- x10[rownames(kgrid),]
x14v1 <- x14v1[rownames(kgrid),]
x14v2 <- x14v2[rownames(kgrid),]

stopifnot(all(rownames(x10)==rownames(kgrid)))
stopifnot(all(rownames(x14v1)==rownames(kgrid)))
stopifnot(all(rownames(x14v2)==rownames(kgrid)))

groups <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
## need to check colnames and types !!!!!!!!!!!!!!!!!!!!!!!
x10 <- x10[,rownames(groups)]
x14v2 <- x14v2[,rownames(groups)]


get_z <- function(x) {
    rs <- rowSums(x)
    # match here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gs <- groupSums(x[,groups$IS_HF], 2, groups$VEGHF[groups$IS_HF])
    100 * gs / rs
}
z10 <- get_z(x10)
z14v1 <- get_z(x14v1)
z14v2 <- get_z(x14v2)
stopifnot(all(colnames(z10)==colnames(z14v1)))
stopifnot(all(colnames(z10)==colnames(z14v2)))
stopifnot(all(colnames(z14v2)==colnames(z14v1)))

rast_ix <- function(i, x) {
    val <- as.numeric(x[,i])
    mat <- as.matrix(Xtab(val ~ Row + Col, kgrid))
    mat[is.na(mat0)] <- NA
#    mat[1] <- 0
#    mat[2] <- 100
    raster(x=mat, template=rt)
}
rast_pl <- function(i, x, r=NULL, fact=NULL, NR=TRUE, ...) {
    if (is.null(r)) {
        r <- rast_ix(i, x)
        col <- colSeq
    } else {
        col <- colDiv
    }
    if (!is.null(fact))
        r <- aggregate(r, fact, fun=mean)
    plot(r, axes=FALSE, box=FALSE, col=col, ...)
    if (NR)
        plot(ABnr, add=TRUE, border=1, lwd=0.5)
    invisible(r)
}

fact <- NULL
NR <- FALSE

pdf("e:/peter/sppweb2017/footprint-results/hf-2014v2-vs-2010-maps.pdf",
    onefile=TRUE, height=8, width=15)
## compares 2014v2 to 2010
#i <- "CCDecid"
for (i in colnames(z10)) {
cat(i, "\n");flush.console()
op <- par(mar=rep(2, 4), mfrow=c(1,3))
r1 <- rast_pl(i, z10, main=paste("2010\n", i), fact=fact, NR=NR)
r2 <- rast_pl(i, z14v2, main=paste("2014v2\n", i), fact=fact, NR=NR)
r3 <- r2-r1
r3@data@values[1] <- -100
r3@data@values[2] <- +100
rast_pl(i, r=r3, main=paste("2014-2010\n", i), NR=NR)
par(op)
}
dev.off()

## compares 2014v2 to 2014v1
pdf("e:/peter/sppweb2017/footprint-results/hf-2014v2-vs-2014v1-maps.pdf",
    onefile=TRUE, height=8, width=15)
#i <- "CCDecid"
for (i in colnames(z10)) {
cat(i, "\n");flush.console()
op <- par(mar=rep(2, 4), mfrow=c(1,3))
r1 <- rast_pl(i, z14v1, main=paste("2014v1\n", i), fact=fact, NR=NR)
r2 <- rast_pl(i, z14v2, main=paste("2014v2\n", i), fact=fact, NR=NR)
r3 <- r2-r1
r3@data@values[1] <- -100
r3@data@values[2] <- +100
rast_pl(i, r=r3, main=paste("v2-v1\n", i), NR=NR)
par(op)
}
dev.off()

## TODO: simplify v1 headers
## plot regional averages etc.

