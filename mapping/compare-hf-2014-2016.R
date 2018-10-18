library(cure4insect)
set_options(path = "w:/reports")
load_common_data()

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
#library(gstat)
library(raster)
#library(viridis)
load(file.path("e:/peter/AB_data_v2017", "data", "analysis", "kgrid_table_km.Rdata"))
load("e:/peter/AB_data_v2018/data/analysis/ages-by-nsr.Rdata")
source("~/repos/abmianalytics/R/veghf_functions.R")

xy <- kgrid
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rt <- .read_raster_template()
mat0 <- as.matrix(rt)
xy <- spTransform(xy, proj4string(rt))

od <- setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(rt))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province
setwd(od)
ABnrS <- gSimplify(ABnr, tol=500, topologyPreserve=TRUE)

colDiv <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B",
    "#FFFFBF","#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
colSeq <- rev(viridis::magma(100))
COL <- RColorBrewer::brewer.pal(6, "Dark2")

rast_ix <- function(i, x)
    .make_raster(x[,i], kgrid, rt)
rast_pl <- function(i, x, r=NULL, fact=NULL, NR=TRUE, col=NULL, ...) {
    if (is.null(r)) {
        r <- .make_raster(x[,i], kgrid, rt)
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
plot_ix <- function(i) {
    m <- cbind(zz14[,i], zz16[,i])
    M <- max(m)
    plot(m[,1], m[,2], xlim=c(0,M), ylim=c(0, M), main=i, xlab="% 2014", ylab="% 2016",
        col=COL, cex=1.5, pch=3)
    abline(0,1,lty=1,col="grey")
    text(m[,1], m[,2], rownames(m), pos=2, col=COL)
}

fact <- 10
NR <- TRUE

## comparing HF inventories

TYPE <- "baseyr2014wUnknAges"
## 2014v2HFI base year 2014 with unknown ages
f14 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine_unsorted.Rdata"
## 2016v3HFI base year 2016 with unknown ages
f16 <- "e:/peter/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3-unsorted.Rdata"

TYPE <- "baseyr2014woUnknAges"
## 2014v2HFI base year 2014 without unknown ages
f14 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine-fixage0.Rdata"
## 2016v3HFI base year 2016 without unknown ages
f16 <- "e:/peter/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3.Rdata"

TYPE <- "baseyr2016wUnknAges"
## 2014v2HFI base year 2016 with unknown ages
f14 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine_unsorted-2016ref.Rdata"
## 2016v3HFI base year 2016 with unknown ages
f16 <- "e:/peter/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3-unsorted.Rdata"

TYPE <- "baseyr2016woUnknAges"
## 2014v2HFI base year 2016 without unknown ages
f14 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine-fixage0-2016ref.Rdata"
## 2016v3HFI base year 2016 without unknown ages
f16 <- "e:/peter/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3.Rdata"

## 2014v2HFI base year 2014 with unknown ages
"e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine_unsorted.Rdata"
## 2014v2HFI base year 2014 without unknown ages
"e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine-fixage0.Rdata"
## 2014v2HFI base year 2016 with unknown ages
"e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine_unsorted-2016ref.Rdata"
## 2014v2HFI base year 2016 without unknown ages
"e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_fine-fixage0-2016ref.Rdata"

## 2016v3HFI base year 2016 with unknown ages
"e:/peter/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3-unsorted.Rdata"
## 2016v3HFI base year 2016 without unknown ages
"e:/peter/AB_data_v2018/data/analysis/grid/veg-hf_grid_v6hf2016v3.Rdata"

e <- new.env()
load(f14, envir=e)
x14 <- e$dd1km_pred$veg_current
x14 <- x14[rownames(kgrid),]
stopifnot(all(rownames(kgrid) == rownames(x14)))

e <- new.env()
load(f16, envir=e)
x16 <- e$dd_kgrid$veg_current
if (is.null(x16))
    x16 <- e$dd$veg_current
x16 <- x16[rownames(kgrid),]
stopifnot(all(rownames(kgrid) == rownames(x16)))

compare_sets(colnames(x14), colnames(x16))
stopifnot(all(colnames(x14) == colnames(x16)))

x14 <- as.matrix(x14)
x16 <- as.matrix(x16)
x14 <- 100 * x14 / rowSums(x14)
x16 <- 100 * x16 / rowSums(x16)

perc <- data.frame(x2014=colMeans(x14), x2016=colMeans(x16))
round(perc, 2)
plot(perc)
abline(0,1)


cn <- colnames(x14)
cn1 <- substr(cn, 1, nchar(cn)-1)
CN <- ifelse(substr(cn, nchar(cn), nchar(cn)) %in% c("R", 0:9), cn1, cn)

xx14 <- groupSums(x14, 2, CN)
xx16 <- groupSums(x16, 2, CN)

kgrid$NRNAMEs <- kgrid$NRNAME
levels(kgrid$NRNAMEs) <- abbreviate(levels(kgrid$NRNAMEs),3)
zz14 <- groupMeans(x14, 1, kgrid$NRNAMEs)
zz16 <- groupMeans(x16, 1, kgrid$NRNAMEs)

## compares 2016 to 2014
pdf(paste0("e:/peter/sppweb2018/footprint-results/hf-2014v2-vs-2016v3-maps_",
    TYPE, "_", if (is.null(fact)) 1 else fact, "km-scale.pdf"),
    onefile=TRUE, height=8, width=15)
for (i in colnames(xx14)) {
cat(i, "\n");flush.console()
op <- par(mar=c(1, 1, 4, 8)+0.1, mfrow=c(1,3))
r1 <- rast_pl(i, xx14, main=paste("2014\n", i), fact=fact, NR=NR)
r2 <- rast_pl(i, xx16, main=paste("2016\n", i), fact=fact, NR=NR)
r3 <- r2-r1
Mx <- max(abs(r3@data@values),na.rm=TRUE)
r3@data@values[1] <- -Mx
r3@data@values[2] <- +Mx
rast_pl(i, r=r3, main=paste("2016-2014\n", i), NR=NR)
par(op)
}
dev.off()

## plot regional averages etc.

pdf(paste0("e:/peter/sppweb2018/footprint-results/hf-by-regions_", TYPE, ".pdf"),
    onefile=TRUE, height=5, width=5)
for (i in colnames(zz14))
    plot_ix(i)
dev.off()
