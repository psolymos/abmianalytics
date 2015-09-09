## Packages
library(raster)
library(sp)
library(rgdal)
library(mefa4)
#library(rasterVis)
library(utils)
source("~/repos/abmianalytics/R/maps_functions.R")

#ROOT <- "c:/p"
#VER <- "AB_data_v2015"
ROOT <- "c:/Users/Peter"
VER <- "www"

## cell x veg/soil matrices and xy lookup table (climate, region, etc)
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

## mammals
dir_in <- "e:/peter/sppweb2015/Mammals/Km2 summaries"
dir_out <- "e:/peter/sppweb-ftp-content/species/mammals"
lt <- read.csv("e:/peter/sppweb2015/Mammals/mammals-lookup.csv")

## mites
dir_in <- "e:/peter/sppweb2015/Mites/Km2 summaries"
dir_out <- "e:/peter/sppweb-ftp-content/species/mites"
lt <- read.csv("e:/peter/sppweb2015/Mites/mites-lookup.csv")

## vplants
dir_in <- "c:/Users/Peter/www/Plants/Km2 summaries"
dir_out <- "c:/Users/Peter/www/species/vplants"
lt <- read.csv("c:/Users/Peter/www/vplants-lookup.csv")

## birds

## mosses
dir_in <- "c:/Users/Peter/www/Moss/Km2 summaries"
dir_out <- "c:/Users/Peter/www/species/mosses"
lt <- read.csv("c:/Users/Peter/www/mosses-lookup.csv")

## lichens
dir_in <- "c:/Users/Peter/www/Lichens/Km2 summaries"
dir_out <- "c:/Users/Peter/www/species/lichens"
lt <- read.csv("c:/Users/Peter/www/lichens-lookup.csv")

#spp <- "CanadaLynx"
spplist <- as.character(lt$sppid[lt$map.pred])
for (spp in spplist) {

    gc()
    f_in <- file.path(dir_in, paste0(spp, ".csv"))
    dir.create(file.path(dir_out, spp))
    f_out <- file.path(dir_out, paste0(spp, ".zip"))
    ff <- c(paste0(spp, 
        c("-reference.asc", "-current.asc", "-reference.png", "-current.png")),
        "README.md")
    f_olist <- file.path(dir_out, spp, ff)

    tab <- read.csv(f_in)
    tab <- tab[match(rownames(kgrid), tab$LinkID),]
    tab$RefNA <- ifelse(is.na(tab$Ref), 1L, 0L)
    tab$CurrNA <- ifelse(is.na(tab$Ref), 1L, 0L)
    tab$Ref[tab$RefNA == 1L] <- 0
    tab$Curr[tab$CurrNA == 1L] <- 0

    NAM <- as.character(lt$species[lt$sppid == spp])

    na_rf <- as_Raster0(kgrid$Row, kgrid$Col, tab$RefNA, rt)
    na_cr <- as_Raster0(kgrid$Row, kgrid$Col, tab$CurrNA, rt)

    cat(spp, "\n");flush.console()
    png(f_olist[3], height=1200, width=800)
    r_rf <- map_fun(tab$Ref, 
        main=paste0(NAM, "\nreference  relative abundance"), 
        colScale="abund", q=1, 
        maskRockies=TRUE, plotWater=TRUE, mask=na_rf)
    dev.off()
    writeRaster(r_rf, f_olist[1], overwrite=TRUE)

    png(f_olist[4], height=1200, width=800)
    r_cr <- map_fun(tab$Curr, 
        main=paste0(NAM, "\ncurrent relative abundance"), 
        colScale="abund", q=1, 
        maskRockies=TRUE, plotWater=TRUE, mask=na_cr)
    dev.off()
    writeRaster(r_cr, f_olist[2], overwrite=TRUE)

    readme <- c(paste0("# ", NAM),
        "\n## Contents\n",
        paste0("* ", spp, 
            c("-reference.asc", "-current.asc", "-reference.png", "-current.png")),
        "\n## Version\n",
        "Alberta Biodiversity Monitoring Institute, species website",
        "Version 3, http://species.abmi.ca",
        "\n## Raster file information\n",
        "ASCII grid format\n",
        "CRS: '+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'")
    writeLines(readme, file.path(dir_out, spp, "README.md"))

    setwd(dir_out)
    zip(f_out, paste0("./", spp, "/", ff))
    unlink(file.path(dir_out, spp), recursive=TRUE, force=TRUE)
}


