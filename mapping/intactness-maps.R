## using veg based models only

library(mefa4)
library(raster)
library(sp)
library(rgdal)

ROOT <- "c:/p/AB_data_v2015"
ROOT2 <- "e:/peter/sppweb2015-round3/Species Intactness"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))

source("~/repos/abmianalytics/R/maps_functions.R")

## raster template to use
rt <- raster(file.path(ROOT, "data", "kgrid", "AHM1k.asc"))
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

H <- 1000 
W <- 600

## non birds

xx <- c(#birds="Birds_Species Intactness",
    lichens="Lichens_Species Intactness",
    mosses="Mosses_Species Intactness",
    mammals="Mammals_Species Intactness",
    mites="Mites_Species Intactness",
    vplants="VPlants_Species Intactness")

for (zz in 1:length(xx)) {
DIR <- xx[zz] # "Birds_Species Intactness"
TAXA <- names(xx)[zz] # "birds"

fl <- list.files(file.path(ROOT2, DIR))
names(fl) <- sapply(strsplit(fl, " "), "[[", 1)
tab <- read.csv(paste0("c:/Users/Peter/repos/abmispecies/_data/", TAXA, ".csv"))
rownames(tab) <- tab$sppid
if (is.null(tab$nonnative))
    tab$nonnative <- FALSE
use <- rownames(tab)[!tab$nonnative & tab$map.pred]

res <- matrix(NA, nrow(kgrid), length(use))
rownames(res) <- rownames(kgrid)
colnames(res) <- use
for (i in 1:length(use)) {
    cat(i, "/", length(use), "-", TAXA, use[i], "\n")
    flush.console()
    e <- new.env()
    load(file.path(ROOT2, DIR, fl[i]), envir=e)
    SI <- e$SI
    res[,i] <- SI$SI[match(rownames(kgrid), SI$LinkID)]
}

SI <- res
meanSI <- rowMeans(SI, na.rm=TRUE)
meanSI[is.na(meanSI)] <- 100

## make raster
png(file.path(ROOT2, "out", paste0("intactness-", TAXA, ".png")), width=W, height=H)
op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(meanSI, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
xxx <- switch(TAXA,
    birds="Birds",
    lichens="Lichens",
    mosses="Btyophytes",
    mammals="Mammals",
    mites="Mites",
    vplants="Vascular Plants (Native)")
title(main=paste0("Intactness - ", xxx))
par(op)
dev.off()
writeRaster(r, file.path(ROOT2, "out", paste0("intactness-", TAXA, ".asc")),
    overwrite=TRUE)

}


## birds
TAXA <- "birds"
fl <- list.files("e:/peter/sppweb2015/birds-si")
names(fl) <- sapply(strsplit(fl, "."), "[[", 1)
tab <- read.csv(paste0("c:/Users/Peter/repos/abmispecies/_data/", TAXA, ".csv"))
rownames(tab) <- tab$sppid
if (is.null(tab$nonnative))
    tab$nonnative <- FALSE
use <- rownames(tab)[!tab$nonnative & tab$map.pred]

res <- matrix(NA, nrow(kgrid), length(use))
rownames(res) <- rownames(kgrid)
colnames(res) <- use
for (i in 1:length(use)) {
    cat(i, "/", length(use), "-", TAXA, use[i], "\n")
    flush.console()
    SI <- read.csv(file.path("e:/peter/sppweb2015/birds-si", fl[i]))
    res[,i] <- SI$SI[match(rownames(kgrid), SI$LinkID)]
}

SI <- res
meanSI <- rowMeans(SI, na.rm=TRUE)
#meanSI <- rowMeans(SI)
meanSI[is.na(meanSI)] <- 0

## make raster
png(file.path(ROOT2, "out", paste0("intactness-", TAXA, ".png")), width=W, height=H)
#op <- par(mar=c(0, 0, 4, 0) + 0.1)
r <- map_fun(meanSI, q=1, main="", colScale="intactness",
    plotWater=TRUE, maskRockies=TRUE, plotCities=TRUE, 
    legend=TRUE,
    mask=NULL)
title(main="Intactness - Birds")
#par(op)
dev.off()
writeRaster(r, file.path(ROOT2, "out", paste0("intactness-", TAXA, ".asc")),
    overwrite=TRUE)








#source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
regs <- levels(kgrid$LUFxNSR)
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useS <- kgrid$NRNAME == "Grassland"

e <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-full-withrevisit.Rdata"), envir=e)
tax <- e$TAX
rm(e)
tax$file <- nameAlnum(as.character(tax$English_Name), "mixed", "")

load(file.path(ROOT, "out", "transitions", paste0(regs[1], ".Rdata")))
Aveg <- rbind(colSums(trVeg))
rownames(Aveg) <- regs[1]
colnames(Aveg) <- colnames(trVeg)
Asoil <- rbind(colSums(trSoil))
rownames(Asoil) <- regs[1]
colnames(Asoil) <- colnames(trSoil)

for (i in 2:length(regs)) {
    cat(regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "transitions", paste0(regs[i], ".Rdata")))
    Aveg <- rbind(Aveg, colSums(trVeg))
    rownames(Aveg) <- regs[1:i]
    Asoil <- rbind(Asoil, colSums(trSoil))
    rownames(Asoil) <- regs[1:i]
}
## m^2 to ha
Aveg <- Aveg / 10^4
Asoil <- Asoil / 10^4


city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
city <- spTransform(city, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
city <- as.data.frame(city)

cex <- 0.25
legcex <- 1.5

Col1 <- rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4"))  # Colour gradient for reference and current
Col1fun <- colorRampPalette(Col1, space = "rgb") # Function to interpolate among these colours for reference and current
C1 <- Col1fun(100)
Col2 <- c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221")  # Colour gradient for difference map
Col2fun <- colorRampPalette(Col2, space = "rgb") # Function to interpolate among these colours for difference map
C2 <- Col2fun(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude

q <- 0.99
H <- 1000 
W <- 600



