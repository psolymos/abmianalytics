library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
#ROOT <- "~/GoogleWork/tmp"

en <- new.env()
load("d:/abmi/AB_data_v2018/data/analysis/birds/data/ab-birds-mixedwood-2019-10-31.RData", envir=en)
xn <- en$DAT
Xn <- get_model_matrix(xn, en$mods)
SPP <- colnames(en$YY)
Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]


spp <- "ALFL"
resn <- load_species(file.path(ROOT, "out", "mixedwood", paste0(spp, ".RData")))

estn <- get_coef(resn, Xn, stage="Space", na.out=FALSE)
printCoefmat(get_summary(estn))

mu <- Xn %*% t(estn[,colnames(Xn)])
lam <- t(apply(exp(mu), 1, quantile, c(0.5, 0.05, 0.95)))
lbc <- aggregate(list(D=lam[,1]), list(LandCov=xn$vegc), summary)
rownames(lbc) <- lbc[,1]
lbc[,1] <- NULL
round(lbc, 4)

RES <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    resn <- load_species(file.path(ROOT, "out", "mixedwood", paste0(spp, ".RData")))
    estn <- get_coef(resn, Xn, stage="Space", na.out=FALSE)
    mu <- Xn %*% t(estn[,colnames(Xn)])
    lam <- t(apply(exp(mu), 1, quantile, c(0.5, 0.05, 0.95)))
    lbc <- aggregate(list(D=lam[,1]), list(LandCov=xn$vegc), summary)
    rownames(lbc) <- lbc[,1]
    lbc[,1] <- NULL
    RES[[spp]] <- lbc
}
save(RES, file=file.path(ROOT, "mixedwood-summaries.RData"))



## predicting for arbitrary polygons based on mixedwood estimates


## make a raster stack of all the climate stuff
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y

xclim <- data.frame(
    transform_clim(kgrid),
    pAspen=kgrid$pAspen,
    pWater_KM=kgrid$pWater,
    pWater2_KM=kgrid$pWater^2)

library(cure4insect)
load_common_data()
r <- .read_raster_template()
l <- list()
for (i in colnames(xclim))
    l[[i]] <- .make_raster(xclim[,i], kgrid, r)
l <- stack(l)
writeRaster(l, "d:/collaborations/foresite/climate-stack.grd")


## the relevant bits start here -----------------------------------------------

library(raster)
library(mefa4)

# some functions defined here that we'll need
source("d:/collaborations/foresite/00-functions.R")

# this is the data used to fit the model
en <- new.env()
load("d:/collaborations/foresite/ab-birds-mixedwood-2019-10-31.RData", envir=en)

# this is your input data
x <- read.csv("d:/collaborations/foresite/first28rows.csv")

# this has the spatial predictors as a raster stack
l <- list()
for (i in 1:14) {
    tmp <- raster("d:/collaborations/foresite/climate-stack.grd", band=i)
    l[[names(tmp)[1]]] <- tmp
}
l <- stack(l)

# check stand types
compare_sets(x$StandType, en$DAT$vegc)
setdiff(x$StandType, xn$vegc) # class '0' won't work --> drop
x <- droplevels(x[x$StandType != "0",])
compare_sets(x$StandType, en$DAT$vegc)

# these are the terms we need
get_terms(en$mods, "list")

# stand type
x$vegc <- x$StandType

# stand type dummies for age interactions
x$isMix <- ifelse(grepl("Mixedwood", as.character(x$vegc), fixed=TRUE), 1, 0)
x$isWSpruce <- ifelse(x$vegc %in% c("pureWhiteSpruce"), 1, 0)
x$isPine <- ifelse(x$vegc %in% c("purePine"), 1, 0)
x$isBSpruce <- ifelse(x$vegc %in% c("pureBlackSpruce"), 1, 0)
x$isLarch <- 0
x$isBSLarch <- pmax(x$isLarch, x$isBSpruce)
x$isUpCon <- ifelse(x$vegc %in% c("pureConifer_allcon", "pureConifer_Pl",
    "pureConifer_Sx", "purePine", "pureWhiteSpruce"), 1, 0)
x$isCon <- ifelse(x$vegc %in% c("pureConifer_allcon", "pureConifer_Pl",
    "pureConifer_Sx", "purePine", "pureWhiteSpruce", "pureBlackSpruce"), 1, 0)

# age & origin (assumed not to be harvest)
x$wtAge <- x$AGE_2017 / 200
x$wtAge2 <- x$wtAge^2
x$wtAge05 <- sqrt(x$wtAge)
x$fCC2 <- 0 # assumes no harvest

# other modifiers that we can ignore
x$ROAD <- 0
x$mWell <- 0
x$mSoft <- 0
x$mEnSft <- 0
x$mTrSft <- 0
x$mSeism <- 0
x$YR <- mean(en$DAT$YR)

# project xy
xy <- x[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(l))

# extract space/climate values
v <- ipa <- extract(l, xy)

# get model matrix for prediction
Xn <- get_model_matrix(data.frame(x, v), en$mods)

# predict

spp <- "ALFL"
# this is the model output for the species
resn <- load_species(file.path(ROOT, "out", "mixedwood", paste0(spp, ".RData")))

estn <- get_coef(resn, Xn, stage="Space", na.out=FALSE)

mu <- Xn %*% t(estn[,colnames(Xn)])
lam <- t(apply(exp(mu), 1, quantile, c(0.5, 0.05, 0.95)))
summary(lam)

