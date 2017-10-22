setwd("~/repos/abmianalytics/apps/hfchange")
source("globals.R")

## shape files for boundary
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/nsr")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
setwd(od)
AB <- spTransform(AB, crs)
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABnr <- gSimplify(ABnr, tol=500, topologyPreserve=TRUE)
ABnsr <- gUnaryUnion(AB, AB@data$NSRNAME) # natural regions
ABnsr <- gSimplify(ABnsr, tol=500, topologyPreserve=TRUE)
pts <- gis
coordinates(pts) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(pts) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
pts <- spTransform(pts, crs)
#ab <- spTransform(ABnr, '+proj=longlat +datum=WGS84')
od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/luf")
ABluf <- readOGR(".", "LUF Integrated Regional Plan Boundaries") # rgdal
setwd(od)
ABluf <- spTransform(ABluf, crs)
ABluf <- gUnaryUnion(ABluf, ABluf@data$LUF_NAME) # natural regions
ABluf <- gSimplify(ABluf, tol=500, topologyPreserve=TRUE)
SHP <- list(nr=ABnr, nsr=ABnsr, luf=ABluf)
#save(ABnr, ABnsr, ABluf, pts, file=file.path(DIR, "AB_NatReg_Simplified.Rdata"))
#load(file.path(DIR, "AB_NatReg_Simplified.Rdata"))


## raster manipulation

f14v2 <- "e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_coarse-fixage0.Rdata"
e <- new.env()
load(f14v2, envir=e)
x14v2 <- e$dd1km_pred$veg_current
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v6.csv")

x <- as.matrix(x14v2)[,rownames(tv)]
x <- 100 * x / rowSums(x)
x <- groupSums(x, 2, tv$ETA_UseInAnalysis_Sector)

load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))
x <- x[rownames(kgrid),c("Agriculture", "Forestry", "Energy", "RuralUrban", "Transportation", "Misc")]
colnames(x) <- c("Agriculture", "Forestry", "Energy",
        "Urban", "Transportation", "Other")

rt <- raster(file.path("e:/peter/AB_data_v2016", "data", "kgrid", "AHM1k.asc"))
projection(rt) <- crs
mat0 <- as.matrix(rt)
matNR <- as.matrix(Xtab(as.integer(kgrid$NRNAME) ~ Row + Col, kgrid))
matNSR <- as.matrix(Xtab(as.integer(kgrid$NSRNAME) ~ Row + Col, kgrid))
matLUF <- as.matrix(Xtab(as.integer(kgrid$LUF_NAME) ~ Row + Col, kgrid))

rast_ix <- function(i, x, fact=1) {
    val <- as.numeric(x[,i])
    mat <- as.matrix(Xtab(val ~ Row + Col, kgrid))
    mat[is.na(mat0)] <- NA
    r <- raster(x=mat, template=rt)
    if (fact > 1)
        r <- aggregate(r, fact, fun=mean)
    r
}

FACT <- 4
rr <- lapply(colnames(x), rast_ix, x=x, fact=FACT)
names(rr) <- colnames(x)
RHF <- stack(rr)

matNR[is.na(mat0)] <- NA
matNSR[is.na(mat0)] <- NA
matLUF[is.na(mat0)] <- NA
REG <- list(
    nr = raster(x=matNR, template=rt),
    nsr = raster(x=matNSR, template=rt),
    luf = raster(x=matLUF, template=rt))
if (FACT > 1) {
    REG$nr <- aggregate(REG$nr, FACT, fun=median)
    REG$nsr <- aggregate(REG$nsr, FACT, fun=median)
    REG$luf <- aggregate(REG$luf, FACT, fun=median)
}

## HF by year array
load(file.path("e:/peter/AB_data_v2017", "data", "analysis",
    "veg-hf_3x7_1999-2015-summaries.Rdata"))
f <- function(z) {
    z <- z[,colnames(z) != "NATIVE"]
    z <- z[,c("Agriculture", "Forestry", "Energy", "RuralUrban", "Transportation", "Misc")]
    colnames(z) <- c("Agriculture", "Forestry", "Energy",
        "Urban", "Transportation", "Other")
    z
}
hf <- list()
for (j in names(veg3x7_sector))
    hf[[j]] <- f(veg3x7_sector[[j]])
HF <- array(NA, c(dim(hf[[1]]), length(veg3x7_sector)))
dimnames(HF) <- list(rownames(hf[[1]]), colnames(hf[[1]]), names(veg3x7_sector))
for (i in 1:length(hf))
    HF[,,i] <- hf[[i]]

GIS <- gis[,c("NATURAL_REGIONS", "NATURAL_SUBREGIONS","LANDUSE_FRAMEWORK")]
colnames(GIS) <- c("nr", "nsr", "luf")

save(HF, RHF, REG, GIS, file="hfchange.rda")

d <- get_data0(1:2,1:2,T,"luf")
matplot(d$x, d$y)
get_gplot(1:2,1:2,T,"luf")
get_mape(1:2,1:2,"nsr")




get_rmap <- function(r) {

    nr <- c(4, 5, 3, 1, 2, 6)
    names(nr) <- c("Grassland", "Parkland", "Foothills", "Boreal",
        "Canadian Shield", "Rocky Mountain")
    Show <- nr %in% r

    rrr <- rr[[r[1]]]
    if (length(r) > 1)
        for (i in 2:length(r))
            rrr <- rr[[r[i]]] + rrr
    msk <- rr[[1]]
    msk[rnr %in% r] <- NA
    rrr[!(rnr %in% r)] <- NA
    op <- par(mar=c(1,1,1,1))
    on.exit(par(op))
    plot(rrr, col=rev(viridis::magma(100)), axes=FALSE, box=FALSE)
    plot(msk, col="lightgrey", add=TRUE, legend=FALSE)
    #plot(ABnr[!Show], add=TRUE, col="#808080", border=NA)
    plot(ABnr, add=TRUE, border="darkgrey")

    invisible(NULL)
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

colDiv <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B",
    "#FFFFBF","#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
colSeq <- rev(viridis::magma(100))
