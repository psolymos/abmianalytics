HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

library(sp)
library(rgdal)
library(rgeos)
library(cure4insect)
set_options(path = "d:/abmi/reports")
load_common_data()

ROOT <- "s:/Cure4Insect-tutorial/data/base/attribute-tables/"
#od <- setwd(ROOT)

yrs <- c(1999, 2004:2017)

SITE <- 788

ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]

# function to address origin year rounding in a reproducible manner
age_unround <- function(y) {
    v <- 1:10*10
    j <- c(-4L, 3L, 4L, 2L, -3L, 5L, 0L, -1L, 1L, -2L)
    for (i in 1:10) {
        s <- !is.na(y) & y %% v[i] == 0
        y[s] <- y[s] + j[i]
    }
    y
}

yri <- 2017

pl <- readOGR(file.path(ROOT, paste0("veg-hf-3x7-", SITE, "_", yri, ".shp")))
xx <- pl@data
xx$Origin_Year <- age_unround(xx$Origin_Yea)
xx$Soil_Type_1 <- xx$Soil_Type_
if (is.null(xx$YEAR))
    xx$YEAR <- xx$year
x <- make_vegHF_wide_v6(xx,
    col.label="ABMI_ID",
    col.year=yri,
    col.HFyear="YEAR",
    col.HABIT="Combined_C",
    col.SOIL="Soil_Type_1",
    wide=FALSE, HF_fine=TRUE) # use refined classes
is0 <- endsWith(as.character(x$VEGHFAGEclass), "0")
v <- as.character(x$VEGHFAGEclass)[is0]
x$VEGHFAGEclass[is0] <- paste0(substr(v, 1, nchar(v)-1), "6")
#table(endsWith(as.character(x$VEGHFAGEclass), "0"))
## soild class harmonization
levels(x$SOILHFclass) <- as.character(ts[levels(x$SOILHFclass), "UseInAnalysisCoarse"])
## NA will be treated as 0 (water and HFor)
x$SOILHFclass[x$SOILHFclass %in% setdiff(levels(x$SOILHFclass), get_levels()$soil)] <- NA
x$SOILHFclass <- droplevels(x$SOILHFclass)
## veg class harmonization
levels(x$VEGHFAGEclass) <- as.character(tv[levels(x$VEGHFAGEclass), "CoefTabsBDQT"])
## NA will be treated as 0 (SnowIce)
x$VEGHFAGEclass[x$VEGHFAGEclass %in% setdiff(levels(x$VEGHFAGEclass), get_levels()$veg)] <- NA
x$VEGHFAGEclass <- droplevels(x$VEGHFAGEclass)
levels(x$SOILHFclass) <- c(levels(x$SOILHFclass), "N/A")
x$SOILHFclass[is.na(x$SOILHFclass)] <- "N/A"

pl2 <- pl
pl2@data <- data.frame(label=interaction(x$VEGHFAGEclass, x$SOILHFclass, sep=" ", drop=TRUE))
pl3 <- aggregate(pl2, list(label=pl2@data$label), FUN=identity, dissolve=TRUE)
pl3@data[[2]] <- NULL
tmp <- pl3@data
tmp2 <- strsplit(as.character(pl3@data$label), " ")
tmp$veg <- as.factor(sapply(tmp2, "[[", 1))
tmp$soil <- as.factor(sapply(tmp2, "[[", 2))
tmp$label <- NULL
pl3@data <- tmp
summary(pl3@data)

library(sf)

sf3 <- st_as_sf(pl3)
plot(sf3)

pl3@data$area_sqkm <- area(pl3) / 10^6
sum(pl3@data$area_sqkm)

rgdal::writeOGR(pl3, paste0("s:/Cure4Insect-tutorial/data/base/examples/veghf3x7-", SITE, "_", yri),
    layer="current-veg-soil", driver="ESRI Shapefile", overwrite_layer=TRUE)
# check if it is correct
In <- readOGR(paste0("s:/Cure4Insect-tutorial/data/base/examples/veghf3x7-", SITE, "_", yri,
    "/current-veg-soil.shp"))

rt <- cure4insect::.read_raster_template()
#values(rt)[!is.na(values(rt))] <- -9999


sl <- list()
for (i in levels(pl3@data$soil)) {
    j <- pl3@data$soil == i
    p <- pl3[j,]
    r1 <- rasterize(p, rt, field="area_sqkm", fun=sum)
    r1 <- trim(r1)
}
r1 <- rasterize(pl3, rt, field="area", fun=sum)
r1 <- trim(r1)
plot(r1)

library(velox)
vx <- velox(rt)
vx$rasterize(pl3, field="veg", band=1, background = 0)

r1 <- rasterize(pl3, rt, field="soil", fun=sum)
r2 <- rasterize(pl3, rt, field="veg", fun=sum, getCover=TRUE)


library(fasterize)
library(sf)

sf3 <- st_as_sf(pl3)
r <- fasterize(sf3, rt, field="veg", fun=sum)

## todo
## OK - deal with unknown ages
## OK - unif rounding of ages
## OK - reclass the veg/soil labels as in BDQT
## OK - write polygon data

## - overlay kgrid
## - calculate composition
