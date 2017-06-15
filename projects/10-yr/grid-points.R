library(rgdal)
library(rgeos)
library(sp)
library(raster)
setwd("x:/toPeter/GridPoints727")

pt <- readOGR(".", "GridPoints727")
proj4string(pt)
str(pt@data)

x <- read.csv("x:/toPeter/Export_tables_NativeVeg_X_validationPoints/PointsPhotoPlotBackfillv6Ident.csv")
x1 <- read.csv("x:/toPeter/Export_tables_NativeVeg_X_validationPoints/PointsPhotoPlotIdent.csv")
x2 <- read.csv("x:/toPeter/Export_tables_NativeVeg_X_validationPoints/random_points.csv")

save(x, x1, x2, file="x:/toPeter/Export_tables_NativeVeg_X_validationPoints/data.Rdata")

## --
load("x:/toPeter/Export_tables_NativeVeg_X_validationPoints/data.Rdata")
rm(x1, x2)

library(mefa4)

table(x$duplicates)
summary(as.numeric(table(x$ABMI_SITE)))

## filter: duplicates, HF classes (where LC3 is blank)

levels(x$LC3)[levels(x$LC3) == " "] <- ""
table(lc=x$LC3 == "", lu=x$LU1_LEVEL1 == "")

x$keep <- TRUE
## LU1_LEVEL1 that are disturbed (HF), not NU or PL
Drop <- c("AG", "TR", "UB", "BU", "MI", "AQ", "CI", "RS", "TR", "SE", "IN", "FO", "RC")
x$keep[x$duplicates > 0] <- FALSE
x$keep[x$LU1_LEVEL1 %in% Drop] <- FALSE
x$keep[x$LC3 == ""] <- FALSE
x$keep[x$Combined == ""] <- FALSE

x <- droplevels(x[x$keep,])

## classify ORIGIN_YR and Origin_Year
table(x$ORIGIN_YR)
x$ORIGIN_YR[x$ORIGIN_YR == 18880] <- 1880
x$ORIGIN_YR[x$ORIGIN_YR == 0] <- 1700
table(x$Origin_Year)
x$Origin_Year[x$Origin_Year == 9999] <- 1700

x$pv_yc <- cut(2014 - x$ORIGIN_YR, c(-1, 10, 20, 40, 60, 80, 100, 120, 140, 160, Inf))
levels(x$pv_yc) <- c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")
x$bf_yc <- cut(2014 - x$Origin_Year, c(-1, 10, 20, 40, 60, 80, 100, 120, 140, 160, Inf))
levels(x$bf_yc) <- c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")
levels(x$pv_yc) <- c(levels(x$pv_yc), "")
levels(x$bf_yc) <- c(levels(x$bf_yc), "")

## finalize classes
pvl <- list(
    TUFC = c("Conif", "Forest conifer"),
    TUFD = c("Decid", "Forest deciduous"),
    TUFM = c("Mixedwood", "Forest mixed"),
    BOXC = c("ShrubbyBog", "Bog, Open, permafrost, collapse scar"),
    BOXN = c("ShrubbyBog", "Bog, Open, permafrost, no internal lawns"),
    BTXC = c("TreedBog", "Bog, Wooded, permafrost, collapse scar"),
    BFXC = c("TreedBog", "Bog, Forested, permafrost, collapse scar"),
    BTXN = c("TreedBog", "Bog, Wooded, permafrost, no internal lawns"),
    BFXN = c("TreedBog", "Bog, Forested, permafrost, no internal lawns"),
    BTNN = c("TreedBog", "Bog, Wooded, permafrost or patterning, no internal lawns"),
    BTNR = c("TreedBog", "Bog, Wooded, permafrost or patterning, internal islands of forested peat plateau"),
    BTNI = c("TreedBog", "Bog, Wooded, permafrost or patterning, internal lawns"),
    FOPN = c("GraminoidFen", "Fen, Open, patterning, no internal lawns"),
    FTPN = c("TreedFen", "Fen, Wooded, patterning, no internal lawns"),
    FONS = c("ShrubbyFen", "Fen, Open, permafrost or patterning, shrub cover"),
    FONG = c("GraminoidFen", "Fen, Open, permafrost or patterning, graminoid cover"),
    FTNN = c("TreedFen", "Fen, Wooded, permafrost or patterning, no internal lawns"),
    FTNR = c("TreedFen", "Fen, Wooded, permafrost or patterning, internal islands of forested peat plateau"),
    FTNI = c("TreedFen", "Fen, Wooded, permafrost or patterning, internal lawns"),
    MONG = c("Marsh", "Marsh, Open, permafrost or patterning, graminoid cover"),
    MOTG = c("Marsh", "Marsh, Open, temporary, graminoid cover"),
    MOSG = c("Marsh", "Marsh, Open, seasonal, graminoid cover"),
    MOQG = c("Marsh", "Marsh, Open, semi-permanent to permanent, graminoid cover"),
    MOAG = c("Marsh", "Marsh, Open, alkali, graminoid cover"),
    MOAX = c("Marsh", "Marsh, Open, alkali, non-vegetated"),
    SFNN = c("TreedSwamp", "Swamp, Forested, permafrost or patterning, no internal lawns"),
    STNN = c("TreedSwamp", "Swap, Wooded, permafrost or patterning, no internal lawns"),
    SONS = c("ShrubbySwamp", "Swamp, Open, permafrost or patterning, shrub cover"),
    SOTS = c("ShrubbySwamp", "Swamp, Open, temporary, shrub cover"),
    SOSS = c("ShrubbySwamp", "Swamp, Open, seasonal, shrub cover"),
    SOQS = c("ShrubbySwamp", "Swamp, Open, semi-permanent to permanent, shrub cover"),
    OUST = c("Shrub", "Tall shrub"),
    OUSS = c("Shrub", "Short shrub"),
    OUHG = c("GrassHerb", "Herbaceous grassland"),
    OUHF = c("GrassHerb", "Herbacesous forbs (non-wetland)"),
    OUBR = c("Bare", "Bryophyte (moss, non-wetland)"),
    OWWL = c("Water", "Lake"),
    OWWS = c("Water", "Salt water"),
    OWWR = c("Water", "River"),
    OWWA = c("HF", "Reservoir"),
    OWWW = c("Water", "Shallow open water"),
    OWWT = c("Water", "Stream"),
    SISC = c("SnowIce", "Snow cover"),
    SIGL = c("SnowIce", "Glacier"),
    ROBR = c("Bare", "Bedrock"),
    RORT = c("Bare", "Rubble, talus, blockfield"),
    ROMO = c("Bare", "Moraine"),
    ELBU = c("Burn", "Burned area"),
    ELRS = c("Bare", "River sediments"),
    ELLS = c("Bare", "Pond or lake sediments"),
    ELCC = c("HF", "Clearcut (fresh)"),
    ELRM = c("HF", "Reservoir margin"),
    ELMU = c("Bare", "Mudflat sediment"),
    ELES = c("HF", "Exposed soil or substratum"),
    ELON = c("Bare", "Other non-vegetated, undeveloped"),
    ASAS = c("HF", "artificial surface/material (including mixed surfaces, e.g. suburbia)"),
    WSNL = c("Water", "Wetland, Lentic- Seasonal"),
    WALK = c("Water", "Wetland, Lentic-Alkali"),
    WSMP = c("Water", "Wetland, Lentic-semi to permanent"),
    OUHE = c("GrassHerb", "Vegetated Open Upland Herbaceous undifferentiated"),
    WTMP = c("Water", "Wetland, Lentic- Temporary"))
compare_sets(x$LC3, names(pvl))
setdiff(x$LC3, names(pvl)) # "WSNL" "WALK" "WSMP" "OUHE" "WTMP"
pvmap <- data.frame(a=names(pvl), b=sapply(pvl, "[[", 1))
compare_sets(x$LC3, pvmap[,1])
x$pv_veg <- reclass(x$LC3, pvmap)

bfl <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-v6-comb.csv")
compare_sets(x$Combined, bfl$V6_COMB)
x$bf_veg <- reclass(x$Combined, bfl[,1:2])
levels(x$bf_veg)[levels(x$bf_veg) %in% c("Pine", "Spruce")] <- "Conif"

compare_sets(x$pv_veg, x$bf_veg)
setdiff(x$pv_veg, x$bf_veg)
setdiff(x$bf_veg, x$pv_veg)

met <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
x$NR <- met$NATURAL_REGIONS[match(x$ABMI_SITE, met$SITE_ID)]
x$NSR <- met$NATURAL_SUBREGIONS[match(x$ABMI_SITE, met$SITE_ID)]

x$pv_veg <- as.factor(as.character(x$pv_veg))
x$bf_veg <- as.factor(as.character(x$bf_veg))
Fclass <- c("Conif", "Decid", "Mixedwood", "TreedBog")
x$pv_yc[!(x$pv_veg %in% Fclass)] <- ""
x$bf_yc[!(x$bf_veg %in% Fclass)] <- ""
Levs <- c(paste0(rep(Fclass, each=10), c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")),
    levels(x$pv_veg)[!(levels(x$pv_veg) %in% Fclass)])
x$pv_vag <- factor(paste0(x$pv_veg, x$pv_yc), levels=Levs)
x$bf_vag <- factor(paste0(x$bf_veg, x$bf_yc), levels=Levs)

x$keep <- TRUE
x$keep[x$pv_veg %in% c("UNK", "Burn", "HF")] <- FALSE
x <- droplevels(x[x$keep,])
compare_sets(x$pv_veg, x$bf_veg)
all(levels(x$pv_veg)==levels(x$pv_veg))
all(levels(x$pv_vag)==levels(x$pv_vag))

xt_prov <- Xtab(~ pv_veg + bf_veg, x)
all(rownames(xt_prov)==colnames(xt_prov))
xt_site <- Xtab(~ pv_veg + bf_veg + ABMI_SITE, x)
xt_nr <- Xtab(~ pv_veg + bf_veg + NR, x)
xt_nsr <- Xtab(~ pv_veg + bf_veg + NSR, x)

xa_prov <- Xtab(~ pv_vag + bf_vag, x)
all(rownames(xt_prov)==colnames(xt_prov))
xa_site <- Xtab(~ pv_vag + bf_vag + ABMI_SITE, x)
xa_nr <- Xtab(~ pv_vag + bf_vag + NR, x)
xa_nsr <- Xtab(~ pv_vag + bf_vag + NSR, x)


save(xt_prov, xt_site, xt_nr, xt_nsr,
    xa_prov, xa_site, xa_nr, xa_nsr,
    x, file="x:/toPeter/Export_tables_NativeVeg_X_validationPoints/data-xt.Rdata")

## --
load("x:/toPeter/Export_tables_NativeVeg_X_validationPoints/data-xt.Rdata")
library(mefa4)
source("~/repos/opticut/extras/multiclass.R")
f <- function(x, digits=3)
    round(as.matrix(x/sum(x)), digits)
a <- function(x) sum(diag(x) / sum(x))
k <- function(x) kappa(x, etable(as.matrix(x)))
m <- function(x) multiclass(as.matrix(x))$average
s <- function(x) multiclass(as.matrix(x))$single

met <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rownames(met) <- met$SITE_ID
AS <- t(sapply(xt_site, k))
AA <- t(sapply(xa_site, k))
MS <- t(sapply(xt_site, m))
MA <- t(sapply(xa_site, m))
met2 <- data.frame(met[rownames(AS),], AS, age=AA)
met2 <- met2[!is.na(met2$k),]

## compare margins

xt2 <- lapply(c(All=xt_prov, xt_nr), as.matrix)
xa2 <- lapply(c(All=xa_prov, xa_nr), as.matrix)

col <- "grey"
cols <- colorRampPalette(c("blue","red"))(5)
pdf("x:/toPeter/grid-results/proportions.pdf",
    width=5.5, height=5.5, onefile=TRUE)
for (i in 1:ncol(xt2[[1]])) {
    j <- colnames(xt2[[1]])[i]
    pv <- sapply(xt2, function(z) 100*sum(z[i,])/sum(z))
    bf <- sapply(xt2, function(z) 100*sum(z[,i])/sum(z))
    lim <- c(0, max(0.1, pv, bf))
    d <- ceiling(abs(pv - bf))
    d[d < 1] <- 1
    d[d >= 5] <- 5
    plot(pv, bf, main=j,
        xlim=lim, ylim=lim,
        xlab="Photoplot (%)", ylab="Backfilled V6 (%)", col=cols[d], pch=19)
    text(pv, bf, ifelse(d > 1, names(pv), ""),
        col=cols[d], pos=3, cex=0.4, xpd=NA)
        abline(0,1,col=col,lty=1)
        abline(-1,1,col=col,lty=2)
        abline(1,1,col=col,lty=2)
        abline(-5,1,col=col,lty=3)
        abline(5,1,col=col,lty=3)
}
dev.off()

pdf("x:/toPeter/grid-results/proportions-age.pdf",
    width=5.5, height=5.5, onefile=TRUE)
for (i in 1:ncol(xa2[[1]])) {
    j <- colnames(xa2[[1]])[i]
    pv <- sapply(xa2, function(z) 100*sum(z[i,])/sum(z))
    bf <- sapply(xa2, function(z) 100*sum(z[,i])/sum(z))
    lim <- c(0, max(0.1, pv, bf))
    d <- ceiling(abs(pv - bf))
    d[d < 1] <- 1
    d[d >= 5] <- 5
    plot(pv, bf, main=j,
        xlim=lim, ylim=lim,
        xlab="Photoplot (%)", ylab="Backfilled V6 (%)", col=cols[d], pch=19)
    text(pv, bf, ifelse(d > 1, names(pv), ""),
        col=cols[d], pos=3, cex=0.4, xpd=NA)
        abline(0,1,col=col,lty=1)
        abline(-1,1,col=col,lty=2)
        abline(1,1,col=col,lty=2)
        abline(-5,1,col=col,lty=3)
        abline(5,1,col=col,lty=3)
}
dev.off()


library(rgdal)
library(rgeos)
library(sp)
library(gstat)
library(raster)

met1 <- met
coordinates(met1) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(met1) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#coordinates(met2) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
#proj4string(met2) <-
#    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- met2[,c("PUBLIC_LONGITUDE", "PUBLIC_LATTITUDE")]
coordinates(xy) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy0 <- met[,c("PUBLIC_LONGITUDE", "PUBLIC_LATTITUDE")]
coordinates(xy0) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(xy0) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

r <- raster("~/Dropbox/courses/st-johns-2017/data/ABrasters/dem.asc")
xy_tm <- spTransform(xy, proj4string(r))
setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(r))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province

setwd("x:/toPeter/grid-results")
## hill shade
slope <- terrain(r, opt='slope')
aspect <- terrain(r, opt='aspect')
hill <- hillShade(slope, aspect, 40, 270)

pdf("grid-sites.pdf",width=3.5, height=6)
op <- par(mar=c(0,0,0,0))
plot(hill, col=grey(0:100/100), legend=FALSE, bty="n", box=FALSE, axes=FALSE)
plot(r, legend=FALSE, col=topo.colors(50, alpha=0.35)[26:50], add=TRUE)
plot(xy_tm, add=TRUE, cex=0.4)
par(op)
dev.off()

## NR
pdf("grid-sites2.pdf",width=3.5, height=6)
COL <- c('#e6f5c9','#f4cae4','#b3e2cd','#fff2ae','#fdcdac','#cbd5e8')
op <- par(mar=c(0,0,0,0))
plot(ABnr, col=COL, border=COL)
plot(xy_tm, add=TRUE, cex=0.4)
par(op)
dev.off()



write.csv(as.matrix(xt_prov), file="x:/toPeter/grid-results/Prov_Landcov.csv")
write.csv(as.matrix(xa_prov), file="x:/toPeter/grid-results/Prov_LandcovAge.csv")
for (i in names(xa_nr)) {
    write.csv(as.matrix(xt_nr[[i]]),
        file=paste0("x:/toPeter/grid-results/NR_Landcov-", gsub(" ", "", i), ".csv"))
    write.csv(as.matrix(xa_nr[[i]]),
        file=paste0("x:/toPeter/grid-results/NR_LandcovAge-", gsub(" ", "", i), ".csv"))
}

zz <- data.frame(
    noage=rbind(Prov=c(k(xt_prov), m(xt_prov)), cbind(t(sapply(xt_nr, k)), t(sapply(xt_nr, m)))),
    age=rbind(Prov=c(k(xa_prov), m(xa_prov)), cbind(t(sapply(xa_nr, k)), t(sapply(xa_nr, m))))
    )
write.csv(zz, file="x:/toPeter/grid-results/Accuracy-and-kappa.csv")


## kriging
gr <- data.frame(coordinates(r)[!is.na(values(r)),])
coordinates(gr) <- ~ x + y
proj4string(gr) <- proj4string(r)
gr <- SpatialPixels(gr)

coordinates(met2) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(met2) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
met2 <- spTransform(met2, proj4string(r))
#v <- variogram(k~1, met2)
g1 <- gstat(id = "k", formula = k~1, data = met2, model=NULL)
p1 <- predict(g1, newdata=gr)

g2 <- gstat(id = "age.k", formula = age.k~1, data = met2, model=NULL)
p2 <- predict(g2, newdata=gr)

png("kriging-lc3.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(p1)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
par(op)
dev.off()

png("kriging-lc3age.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(p2)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
par(op)
dev.off()


