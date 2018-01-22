## this scripts does some preprocessing
## and explains how to calculate bird densities based on
## Patchworks output using ABMI predictions and the cure4insect R package

library(cure4insect)
library(mefa4)
library(rgdal)
library(rgeos)
library(sp)

load_common_data()
veg <- get_levels()$veg
## load this file from where you put it
load("hf-an-xy-for patchworks.RData")
qsid <- rownames(coordinates(xy))

## ignor this part
if (FALSE) {
## coordinates for QSs
load("e:/peter/AB_data_v2017/data/analysis/kgrid_table_qs.Rdata")
xy <- kgrid[,c("POINT_X", "POINT_Y", "MER", "RGE", "TWP", "SEC", "QS", "TWNSHIP",
    "SECTION", "Area_km2")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- proj4string(get_id_locations())
head(xy@data)

hf <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
rownames(hf) <- tolower(hf$HF_GROUP)
save(hf, xy, file="hf-an-xy-for patchworks.RData")
}

## define species and directory with the csv input
spp <- "Ovenbird"
DIR <- "c:/Users/Peter/Downloads/Caribou scenario habitat HF combined"
fl <- list.files(DIR, pattern="bothhabitatHF") # filename must have bothhabitatHF
## load species data (spatial/climate part of predictions)
y <- load_spclim_data(spp)

for (i in 1:length(fl)) {
    fn <- fl[i]
    cat(fn, i, "/", length(fl), "\n")
    flush.console()
    d <- read.csv(file.path(DIR, fn))
    rownames(d) <- d$QSLinkID
    d$X <- NULL # rownames saved and included as 1st col
    d$QSLinkID <- NULL # tracking as rownames
    weird <- !(rownames(d) %in% qsid)
    if (any(weird)) {
        print(paste("found & excluded weid row(s):", sum(weird)))
        print(rownames(d)[weird])
        d <- d[!weird,,drop=FALSE]
    }
    colnames(d)[colnames(d) == "NA."] <- "Unknown" # some stuff was unknown
    colnames(d) <- gsub("\\.", "", colnames(d)) # remove dots from colnames
    colnames(d)[colnames(d)=="MunicipalWaterandSewage"] <- "MunicipalWaterSewage"
    ## lookup table
    tab <- data.frame(col_data=colnames(d))
    rownames(tab) <- tab$col_data
    tab$col_in <- as.character(tab$col_data)
    tab$col_in <- gsub("UP", "CC", tab$col_in)
    tab$col_coef <- as.character(veg[match(tab$col_in, veg)])
    hf2 <- hf[tolower(as.character(tab$col_in[is.na(tab$col_coef)])),]
    tab$col_coef[is.na(tab$col_coef)] <- as.character(hf2$UseInAnalysis)
    tab$col_coef[tab$col_in == "Unknown"] <- "Unknown"
    tab$col_coef <- gsub("CC80", "80", tab$col_coef) # CC converged by 60
    tab$col_coef <- gsub("CC100", "100", tab$col_coef) # CC converged by 60
    tab$col_coef <- gsub("CC120", "120", tab$col_coef) # CC converged by 60
    tab$col_coef <- gsub("CC140", "140", tab$col_coef) # CC converged by 60
    tab$col_zero <- tolower(tab$col_in) %in% tolower(c("Water",
        "BorrowPitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals",
        "RailHardSurface", "RoadHardSurface", "MineSite", "PeatMine",
        "Unknown", "Cutblock"))
    veg2 <- intersect(tab$col_coef, veg)
    ## making inputs and predict
    X <- array(1, c(nrow(d), length(veg2)), list(rownames(d), veg2))
    XY <- xy[match(rownames(d), qsid),]
    prmat <- suppressWarnings(predict_mat(y, xy=XY, veg=X)$veg)
    ## recast predictions to fit input file
    d2 <- prmat[,match(tab$col_coef, colnames(prmat))]
    colnames(d2) <- colnames(d)
    d2[,tab$col_zero] <- 0
    ## this calculation assumes that input areas are in m^2
    d3 <- d2 * d / 10^4
    ## write results
    fn_out <- paste0("Density_", spp, "_", fn)
    if (!dir.exists(file.path(DIR, spp)))
        dir.create(file.path(DIR, spp))
    write.csv(d2, file=file.path(DIR, spp, paste0("Density_", spp, "_", fn)))
    write.csv(d3, file=file.path(DIR, spp, paste0("Abundance_", spp, "_", fn)))
}
