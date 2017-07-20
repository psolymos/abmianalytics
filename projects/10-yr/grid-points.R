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
    WSNL = c("GrassHerb", "Wetland, Lentic- Seasonal"),
    WALK = c("Water", "Wetland, Lentic-Alkali"),
    WSMP = c("Water", "Wetland, Lentic-semi to permanent"),
    OUHE = c("GrassHerb", "Vegetated Open Upland Herbaceous undifferentiated"),
    WTMP = c("GrassHerb", "Wetland, Lentic- Temporary"))
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
src <- read.csv("~/repos/abmianalytics/projects/10-yr/3by7Center_VegV6_summary.csv")
x$NR <- met$NATURAL_REGIONS[match(x$ABMI_SITE, met$SITE_ID)]
x$NSR <- met$NATURAL_SUBREGIONS[match(x$ABMI_SITE, met$SITE_ID)]
x$Source <- src$General_Source[match(x$ABMI_SITE, src$ABMI)]

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
xt_src <- Xtab(~ pv_veg + bf_veg + Source, x)

xa_prov <- Xtab(~ pv_vag + bf_vag, x)
all(rownames(xt_prov)==colnames(xt_prov))
xa_site <- Xtab(~ pv_vag + bf_vag + ABMI_SITE, x)
xa_nr <- Xtab(~ pv_vag + bf_vag + NR, x)
xa_nsr <- Xtab(~ pv_vag + bf_vag + NSR, x)
xa_src <- Xtab(~ pv_vag + bf_vag + Source, x)


save(xt_prov, xt_site, xt_nr, xt_nsr, xt_src,
    xa_prov, xa_site, xa_nr, xa_nsr, xa_src,
    x,
    file="x:/toPeter/Export_tables_NativeVeg_X_validationPoints/data-xt.Rdata")

## --
load("x:/toPeter/Export_tables_NativeVeg_X_validationPoints/data-xt.Rdata")
library(mefa4)
library(viridis)
source("~/repos/opticut/extras/multiclass.R")
library(mgcv)

f <- function(x, digits=3)
    round(as.matrix(x/sum(x)), digits)
a <- function(x) sum(diag(x) / sum(x))
k <- function(x) kappa(x, etable(as.matrix(x)))
m <- function(x) multiclass(as.matrix(x))$average
#ss <- function(x) multiclass(as.matrix(x))$single

met <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
src <- read.csv("~/repos/abmianalytics/projects/10-yr/3by7Center_VegV6_summary.csv")
met$BfSource <- droplevels(src$General_Source[match(met$SITE_ID, src$ABMI)])

rownames(met) <- met$SITE_ID
AS <- t(sapply(xt_site, k))
AA <- t(sapply(xa_site, k))
MS <- t(sapply(xt_site, m))
MA <- t(sapply(xa_site, m))
met2 <- data.frame(met[rownames(AS),], AS, age=AA)
met2 <- met2[!is.na(met2$k),]

fullxt <- do.call(rbind, lapply(names(xt_site), function(i) {
    out <- xt_site[[i]]
    rownames(out) <- paste0("Site_", i, "_Pp_", rownames(out))
    colnames(out) <- paste0("Bf_", rownames(out))
    as.matrix(out)
}))

## compare margins

names(xt_nr) <- paste0("NR=", names(xt_nr))
names(xa_nr) <- paste0("NR=", names(xa_nr))
names(xt_src) <- paste0("Src=", names(xt_src))
names(xa_src) <- paste0("Src=", names(xa_src))
xt2 <- lapply(c(Alberta=xt_prov, xt_nr, xt_src), as.matrix)
xa2 <- lapply(c(Alberta=xa_prov, xa_nr, xa_src), as.matrix)

col <- "grey"
cols <- colorRampPalette(c("blue","red"))(10)
CN <- c("Conif", "Decid", "Mixedwood",
    "TreedBog", "TreedFen", "TreedSwamp",
    "ShrubbyBog", "ShrubbyFen", "ShrubbySwamp",
    "GraminoidFen", "Marsh",
    "Shrub", "GrassHerb",
    "SnowIce", "Bare", "Water")
pdf("x:/toPeter/grid-results/proportions.pdf",
    width=5.5, height=5.5, onefile=TRUE)
for (j in CN) {
    #j <- colnames(xt2[[1]])[i]
    pv <- sapply(xt2, function(z) 100*sum(z[j,])/sum(z))
    bf <- sapply(xt2, function(z) 100*sum(z[,j])/sum(z))
    tp <- sapply(xt2, function(z) z[j,j]/sum(z[j,]))
    lim <- c(0, max(0.1, pv, bf))
    d <- pmax(1, ceiling(round(10*tp)))
    r <- pmax(1, ceiling(round(10*(1-pmin(pv, bf) / pmax(pv, bf)))))
    plot(pv, bf, main=j,
        xlim=lim, ylim=lim,
        xlab="Photoplot (%)", ylab="'Backfilled' V6 (%)",
        col=cols[d], pch=19, cex=1)
    abline(0,1,col=col,lty=1)
    abline(0,1/0.95,col=col,lty=2)
    abline(0,0.95,col=col,lty=2)
    abline(0,1/0.9,col=col,lty=3)
    abline(0,0.9,col=col,lty=3)
    text(pv, bf, ifelse(r > 0, names(pv), ""),
        col=1, pos=3, cex=0.4, xpd=NA)
    legend("topleft", bty="n", col=col, lty=1:3,
        legend=c("1:1", "1:0.95", "1:0.9"), title="Guides", cex=0.6)
    legend("bottomright", bty="n", pch=19, col=rev(cols), cex=0.6,
        legend=rev(paste0(0:9/10, "-", 1:10/10)), title="True Positive Rate")
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
        xlab="Photoplot (%)", ylab="'Backfilled' V6 (%)",
        col=cols[d], pch=19)
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
library(akima)
library(mgcv)

xy <- met2
coordinates(xy) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

r <- raster(file.path("~/Dropbox/courses/st-johns-2017", "data", "ABrasters", "dem.asc"))
xy <- spTransform(xy, proj4string(r)) # make CRS identical
## make things a bit coarser for saving some time
## comment this out for full-scale 500m resolution
r <- aggregate(r, 2) # 1km resolution
#r <- aggregate(r, 10) # 5km resolution
#r <- aggregate(r, 20) # 10km resolution
#r <- aggregate(r, 40) # 20km resolution
## crop to the range of bird data points
#r <- crop(r, extent(coordinates(x)))

## spatial grid that we use for prediction
## based on the raster coordinates
## dropping values outside of the province (NA)
g <- data.frame(coordinates(r)[!is.na(values(r)),])
## making a gridded object suitable for gstat package
gridded(g) <- ~ x + y
proj4string(g) <- proj4string(r) # projection

setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(r))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province

setwd("x:/toPeter/grid-results")

colnames(xy@data)
write.csv(xy@data, row.names=FALSE, file="site-level-accuracy.csv")
write.csv(fullxt, row.names=TRUE, file="site-level-crosstabs.csv")

#pol3 <- krige(k ~ 1, locations=xy, newdata=g, degree=3)


df <- data.frame(k=xy@data$a, kk=xy@data$age.a, coordinates(xy))
#df <- data.frame(k=xy@data$k, kk=xy@data$age.k, coordinates(xy))
colnames(df) <- c("k", "kk", "x", "y")

gam_fit <- mgcv::gam(k ~ s(x, y, k=200), df, family=gaussian)
dfpred2 <- data.frame(coordinates(g))
colnames(dfpred2) <- c("x", "y")
dfpred2$pr <- predict(gam_fit, dfpred2)
gridded(dfpred2) <- ~ x + y
#plot(dfpred2)

gam_fita <- mgcv::gam(kk ~ s(x, y, k=200), df, family=gaussian)
dfpred3 <- data.frame(coordinates(g))
colnames(dfpred3) <- c("x", "y")
dfpred3$pr <- predict(gam_fita, dfpred3)
gridded(dfpred3) <- ~ x + y
#plot(dfpred3)

## IDW

g1 <- gstat(id = "k", formula = k~1, data = xy, model=NULL)
p1 <- predict(g1, newdata=g)

g2 <- gstat(id = "age.k", formula = age.k~1, data = xy, model=NULL)
p2 <- predict(g2, newdata=g)

pch <- 1:nlevels(xy@data$BfSource)
col <- plasma(255) # viridis(255)

png("idw-lc3.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(p1, col=col)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
plot(xy, pch=pch[as.integer(xy@data$BfSource)], cex=1.5, col="white", add=TRUE)
legend("bottomleft", col=1, pch=pch, legend=levels(xy@data$BfSource), bty="n",
    title="Backfilled Source", cex=2)
par(op)
dev.off()

png("idw-lc3age.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(p2, col=col)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
plot(xy, pch=pch[as.integer(xy@data$BfSource)], cex=1.5, col="white", add=TRUE)
legend("bottomleft", col=1, pch=pch, legend=levels(xy@data$BfSource), bty="n",
    title="Backfilled Source", cex=2)
par(op)
dev.off()

#png("gam-smooth-lc3-kappa.png",width=1000, height=1500)
png("gam-smooth-lc3-accuracy.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(dfpred2, col=col)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
plot(xy, pch=pch[as.integer(xy@data$BfSource)], cex=1.5, col="white", add=TRUE)
legend("bottomleft", col=1, pch=pch, legend=levels(xy@data$BfSource), bty="n",
    title="Backfilled Source", cex=2)
par(op)
dev.off()

#png("gam-smooth-lc3age-kappa.png",width=1000, height=1500)
png("gam-smooth-lc3age-accuracy.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(dfpred3, col=col)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
plot(xy, pch=pch[as.integer(xy@data$BfSource)], cex=1.5, col="white", add=TRUE)
legend("bottomleft", col=1, pch=pch, legend=levels(xy@data$BfSource), bty="n",
    title="Backfilled Source", cex=2)
par(op)
dev.off()

setwd("x:/toPeter/grid-results")

write.csv(as.matrix(xt_prov), file="x:/toPeter/grid-results/Prov_Landcov.csv")
write.csv(as.matrix(xa_prov), file="x:/toPeter/grid-results/Prov_LandcovAge.csv")
for (i in names(xa_nr)) {
    write.csv(as.matrix(xt_nr[[i]]),
        file=paste0("x:/toPeter/grid-results/NR_Landcov-", gsub(" ", "", i), ".csv"))
    write.csv(as.matrix(xa_nr[[i]]),
        file=paste0("x:/toPeter/grid-results/NR_LandcovAge-", gsub(" ", "", i), ".csv"))
}
for (i in names(xa_src)) {
    write.csv(as.matrix(xt_src[[i]]),
        file=paste0("x:/toPeter/grid-results/NR_Landcov-", gsub(" ", "", i), ".csv"))
    write.csv(as.matrix(xa_src[[i]]),
        file=paste0("x:/toPeter/grid-results/NR_LandcovAge-", gsub(" ", "", i), ".csv"))
}

zz <- data.frame(
    noage=rbind(Prov=c(k(xt_prov), m(xt_prov)),
        cbind(t(sapply(xt_nr, k)), t(sapply(xt_nr, m))),
        cbind(t(sapply(xt_src, k)), t(sapply(xt_src, m)))
    ),
    age=rbind(Prov=c(k(xa_prov), m(xa_prov)),
        cbind(t(sapply(xa_nr, k)), t(sapply(xa_nr, m))),
        cbind(t(sapply(xa_src, k)), t(sapply(xa_src, m)))
    ))
write.csv(zz, file="x:/toPeter/grid-results/Accuracy-and-kappa.csv")

png("gam-smooth-lc3age-accuracy.png",width=1000, height=1500)
op <- par(mar=c(0,0,0,0))
plot(dfpred3, col=col)
plot(ABnr, add=TRUE, col=NA, border="grey", lwd=0.25)
plot(xy, pch=pch[as.integer(xy@data$BfSource)], cex=1.5, col="white", add=TRUE)
legend("bottomleft", col=1, pch=pch, legend=levels(xy@data$BfSource), bty="n",
    title="Backfilled Source", cex=2)
par(op)
dev.off()


COL <- c('#e6f5c9','#f4cae4','#b3e2cd','#fff2ae','#fdcdac','#cbd5e8')
library(RColorBrewer)
geocol <- brewer.pal(7, "Dark2")
xyz <- met[!(met$SITE_ID %in% met2$SITE_ID),]
coordinates(xyz) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(xyz) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xyz <- spTransform(xyz, proj4string(r)) # make CRS identical

#png("compared-sites-map.png",width=1000, height=1500)
pdf("compared-sites-map.pdf", height=8, width=5)
op <- par(mar=c(0,0,0,0))

plot(ABnr, col=paste0(COL, "60"), border=paste0(COL, "80"))
plot(ABpr, border="grey", add=TRUE)
points(xyz, col=geocol[as.integer(xyz@data$BfSource)], pch=19, cex=0.25)
points(xy, col=geocol[as.integer(xy@data$BfSource)], pch=19, cex=0.75)
legend("bottomleft", col=geocol, pch=19, legend=levels(xy@data$BfSource), bty="n",
    title="Backfilled Source", cex=1)

par(op)
dev.off()



geolabelN <- c("1 North-East","2 South-West","3 South","4 South-East","5 North-West",
    "6 Central-East","7 Central-West","8 Centre")
geolabelS <- c("1 North","2 East","3 Centre","4 West", "5 South-East","6 South-West",
    "7 North-West")

cn <- groupMeans(coordinates(sxy)[!is.na(sxy@data$Ngeo),], 1,
    as.character(sxy@data$Ngeo)[!is.na(sxy@data$Ngeo)])
cn <- cn[order(rownames(cn)),]
rownames(cn) <- geolabelN
cs <- groupMeans(coordinates(sxy)[!is.na(sxy@data$Sgeo),], 1,
    as.character(sxy@data$Sgeo)[!is.na(sxy@data$Sgeo)])
cs <- cs[order(rownames(cs)),]
rownames(cs) <- geolabelS

pdf("e:/peter/sppweb2017/geoxv-maps.pdf", onefile=TRUE, height=8, width=10)
op <- par(mar=rep(1,4), mfrow=c(1,2))
plot(ABnr, col=paste0(COL, "40"), border=paste0(COL, "80"), main="Geographic clusters, North")
points(sxy, col=paste0(geocol, "80")[sxy@data$Ngeo], pch=19)
plot(ABpr, border="grey", add=TRUE)
text(cn[,1], cn[,2], rownames(cn))
