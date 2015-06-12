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

if (FALSE) { # test
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


## see rasterTheme help page

}

## HF to plot
hf1 <- (dd1km_pred$veg_current[,!is.na(tveg$HF)]) / rowSums(dd1km_pred$veg_current)
hft <- tveg[!is.na(tveg$HF),]
hf2 <- cbind(TotalHF=rowSums(hf1), as.matrix(groupSums(hf1, 2, hft$UseInAnalysis)))
hf1 <- as.matrix(groupSums(hf1, 2, hft$HF))

## raster stack from hf
as_Raster0 <- 
function(x, y, z, r) 
{
    mat0 <- as.matrix(r)
    mat <- as.matrix(Xtab(z ~ x + y))
    mat[is.na(mat0)] <- NA
    raster(x=mat, template=r)
}
as_Raster <- 
function(x, y, z, r, verbose=0) 
{
    if (is.null(dim(z))) {
        out <- as_Raster0(x, y, z, r)
    } else {
        if (verbose) {
            cat("rasterizing: 1\n")
            flush.console()
        }
        out <- as_Raster0(x, y, z[,1], r)
        for (i in 2:ncol(z)) {
            if (verbose) {
                cat("rasterizing:", i, "\n")
                flush.console()
            }
            tmp <- as_Raster0(x, y, z[,i], r)
            out <- stack(out, tmp)
        }
        names(out) <- if (is.null(colnames(z)))
            paste0("V", seq_len(ncol(z))) else colnames(z)
    }
    out
}
#r_hf <- as_Raster(kgrid$Row, kgrid$Col, rowSums(hf1), rt)
r_hf <- as_Raster(kgrid$Row, kgrid$Col, hf2, rt, verbose=1)

pdf(file.path(ROOT, VER, "out", "figs", "hf", "hf-use-in-analysis.pdf"),
    width=10, height=8)
levelplot(r_hf, margin=FALSE, par.settings=BuRdTheme(), 
    scales=list(draw=FALSE), contour=FALSE, main="Human Footprint Types")
dev.off()
