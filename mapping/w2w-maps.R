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

r_water <- as_Raster(kgrid$Row, kgrid$Col, kgrid$pWater, rt)
r_water[r_water <= 0.5] <- NA
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

map_fun <-
function(x, q=1, main="", colScale="terrain",
plotWater=TRUE, plotCities=TRUE, maskRockies=TRUE)
{

    ## these might need to be truncated to have pale in the middle
    
    ## Colour gradient for reference and current
    Pal1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8",
        "#91BFDB","#4575B4")))(255)
    ## Colour gradient for difference map
    Pal2 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0",
        "#A1D76A","#4D9221"))(255)

    col <- switch(colScale,
        "terrain" = rev(terrain.colors(255)),
        "heat" = rev(heat.colors(255)),
        "topo" = rev(topo.colors(255)),
        "grey" = grey(seq(1,0,len=255)),
        "abund" = Pal1,
        "diff" = Pal2)

    r <- as_Raster(kgrid$Row, kgrid$Col, x, rt)
    r[r > quantile(r, q)] <- quantile(r, q)
    ## output raster should be masked as well
    if (maskRockies)
        r[!is.na(r_mask)] <- NA

    plot(r, axes=FALSE, box=FALSE, legend=TRUE, 
        main=main, maxpixels=10^6, col=col)
    if (plotWater && !maskRockies) {
        plot(r_water, add=TRUE, 
            alpha=1, legend=FALSE, 
            col="#664DCC")
    }
    if (maskRockies && !plotWater) {
        plot(r_mask, add=TRUE, alpha=1, 
            col="#7A8B8B", legend=FALSE)
    }
    if (plotWater && maskRockies) {
        plot(r_overlay, add=TRUE, 
            alpha=1, legend=FALSE, 
            col=colorRampPalette(c("#664DCC", "#7A8B8B"))(255))
    }
    
    if (plotCities) {
        points(city, pch=18, col="grey10")
        text(city, labels=rownames(city@coords), cex=0.8, pos=3, col="grey10")
    }
    invisible(r)
}




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


## plot with some cities and water

