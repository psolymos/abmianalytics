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
ABnrS <- gSimplify(ABnr, tol=500, topologyPreserve=TRUE)
object.size(ABnrS)/object.size(ABnr)

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
setdiff(rownames(groups), colnames(x14v1))
setdiff(colnames(x14v1), rownames(groups))

x14v1 <- cbind(as.matrix(x14v1),
    CultivationCropPastureBareground=rowSums(x14v1[,c("CultivationCrop",
    "CultivationAbandoned", "CultivationRoughPasture", "CultivationTamePasture")]),
    SeismicLine=rowSums(x14v1[,c("SeismicLineNarrow", "SeismicLineWide")]),
    Urban=rowSums(x14v1[,c("UrbanIndustrial", "UrbanResidence")]))
x14v1 <- x14v1[,rownames(groups)]

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
rast_pl <- function(i, x, r=NULL, fact=NULL, NR=TRUE, col=NULL, ...) {
    if (is.null(r)) {
        r <- rast_ix(i, x)
        if (is.null(col))
            col <- colSeq
    } else {
        if (is.null(col))
            col <- colDiv
    }
    if (!is.null(fact))
        r <- aggregate(r, fact, fun=mean)
    plot(r, axes=FALSE, box=FALSE, col=col, ...)
    if (NR)
        plot(ABnrS, add=TRUE, border=1, lwd=0.5)
    invisible(r)
}

fact <- NULL
NR <- TRUE

## compares 2014v2 to 2010
pdf(paste0("e:/peter/sppweb2017/footprint-results/hf-2014v2-vs-2010-maps_",
    if (is.null(fact)) 1 else fact, "km-scale.pdf"),
    onefile=TRUE, height=8, width=15)
for (i in colnames(z10)) {
cat(i, "\n");flush.console()
op <- par(mar=c(1, 1, 4, 8)+0.1, mfrow=c(1,3))
r1 <- rast_pl(i, z10, main=paste("2010\n", i), fact=fact, NR=NR)
r2 <- rast_pl(i, z14v2, main=paste("2014v2\n", i), fact=fact, NR=NR)
r3 <- r2-r1
Mx <- max(abs(r3@data@values),na.rm=TRUE)
r3@data@values[1] <- -Mx
r3@data@values[2] <- +Mx
rast_pl(i, r=r3, main=paste("2014-2010\n", i), NR=NR)
par(op)
}
dev.off()

## compares 2014v2 to 2014v1
pdf(paste0("e:/peter/sppweb2017/footprint-results/hf-2014v2-vs-2014v1-maps_",
    if (is.null(fact)) 1 else fact, "km-scale.pdf"),
    onefile=TRUE, height=8, width=15)
for (i in colnames(z10)) {
cat(i, "\n");flush.console()
op <- par(mar=c(1, 1, 4, 8)+0.1, mfrow=c(1,3))
r1 <- rast_pl(i, z14v1, main=paste("2014v1\n", i), fact=fact, NR=NR)
r2 <- rast_pl(i, z14v2, main=paste("2014v2\n", i), fact=fact, NR=NR)
r3 <- r2-r1
Mx <- max(abs(r3@data@values),na.rm=TRUE)
r3@data@values[1] <- -Mx
r3@data@values[2] <- +Mx
rast_pl(i, r=r3, main=paste("v2-v1\n", i), NR=NR)
par(op)
}
dev.off()

## TODO:
## plot regional averages etc.

kgrid$NRNAMEs <- kgrid$NRNAME
levels(kgrid$NRNAMEs) <- abbreviate(levels(kgrid$NRNAMEs),3)
zz10 <- groupMeans(z10, 1, kgrid$NRNAMEs)
zz14v1 <- groupMeans(z14v1, 1, kgrid$NRNAMEs)
zz14v2 <- groupMeans(z14v2, 1, kgrid$NRNAMEs)

COL <- RColorBrewer::brewer.pal(6, "Dark2")
plot_ix <- function(i) {
    m <- cbind(zz10[,i], zz14v1[,i], zz14v2[,i])
    M <- max(m)
    op <- par(mfrow=c(1,2))
    plot(m[,1], m[,3], xlim=c(0,M), ylim=c(0, M), main=i, xlab="% 2010", ylab="% 2014v2",
        col=COL, cex=1.5, pch=3)
    abline(0,1,lty=1,col="grey")
    text(m[,1], m[,3], rownames(m), pos=2, col=COL)
    plot(m[,2], m[,3], xlim=c(0,M), ylim=c(0, M), main=i, xlab="% 2014v1", ylab="% 2014v2",
        col=COL, cex=1.5, pch=3)
    abline(0,1,lty=1,col="grey")
    text(m[,2], m[,3], rownames(m), pos=2, col=COL)
    par(op)
}

pdf("e:/peter/sppweb2017/footprint-results/hf-by-regions.pdf",
    onefile=TRUE, height=5, width=10)
for (i in colnames(zz10))
    plot_ix(i)
dev.off()

## total HF

z <- cbind(THF2010=rowSums(z10), THF2014v1=rowSums(z14v1), THF2014v2=rowSums(z14v2))
zz <- groupMeans(z, 1, kgrid$NRNAMEs)
write.csv(zz, "e:/peter/sppweb2017/footprint-results/hf-total-by-regions.csv")
pdf("e:/peter/sppweb2017/footprint-results/hf-total-by-regions.pdf", height=5, width=10)
    m <- zz
    M <- max(m)
    op <- par(mfrow=c(1,2))
    plot(m[,1], m[,3], xlim=c(0,M), ylim=c(0, M), main=i, xlab="% 2010", ylab="% 2014v2",
        col=COL, cex=1.5, pch=3)
    abline(0,1,lty=1,col="grey")
    text(m[,1], m[,3], rownames(m), pos=2, col=COL)
    plot(m[,2], m[,3], xlim=c(0,M), ylim=c(0, M), main=i, xlab="% 2014v1", ylab="% 2014v2",
        col=COL, cex=1.5, pch=3)
    abline(0,1,lty=1,col="grey")
    text(m[,2], m[,3], rownames(m), pos=2, col=COL)
    par(op)
dev.off()

val <- as.numeric(z[,1])
mat <- as.matrix(Xtab(val ~ Row + Col, kgrid))
mat[is.na(mat0)] <- NA
r1 <- raster(x=mat, template=rt)

val <- as.numeric(z[,2])
mat <- as.matrix(Xtab(val ~ Row + Col, kgrid))
mat[is.na(mat0)] <- NA
r2 <- raster(x=mat, template=rt)

val <- as.numeric(z[,3])
mat <- as.matrix(Xtab(val ~ Row + Col, kgrid))
mat[is.na(mat0)] <- NA
r3 <- raster(x=mat, template=rt)

pdf(paste0("e:/peter/sppweb2017/footprint-results/hf-total-maps_",
    if (is.null(fact)) 1 else fact, "km-scale.pdf"),
    height=16, width=15)
op <- par(mar=c(1, 1, 4, 8)+0.1, mfrow=c(2,3))

rr1 <- rast_pl("Total HF", r=r1, main=paste("2010\n", "Total HF"), fact=fact, NR=NR, col=colSeq)
rr2 <- rast_pl("Total HF", r=r3, main=paste("2014v2\n", "Total HF"), fact=fact, NR=NR, col=colSeq)
rr <- rr2-rr1
Mx <- max(abs(rr@data@values),na.rm=TRUE)
rr@data@values[1] <- -Mx
rr@data@values[2] <- +Mx
rast_pl("Total HF", r=rr, main=paste("2014-2010\n", "Total HF"), NR=NR)

rr1 <- rast_pl("Total HF", r=r2, main=paste("2014v1\n", "Total HF"), fact=fact, NR=NR, col=colSeq)
rr2 <- rast_pl("Total HF", r=r3, main=paste("2014v2\n", "Total HF"), fact=fact, NR=NR, col=colSeq)
rr <- rr2-rr1
Mx <- max(abs(rr@data@values),na.rm=TRUE)
rr@data@values[1] <- -Mx
rr@data@values[2] <- +Mx
rast_pl("Total HF", r=rr, main=paste("v2-v1\n", "Total HF"), NR=NR)

par(op)
dev.off()

