library(mefa4)

#load("e:/peter/AB_data_v2016/data/species/OUT_mites_2016-07-28.Rdata")
#load("e:/peter/AB_data_v2016/data/species/OUT_vplants_2016-08-18.Rdata")
load("e:/peter/AB_data_v2016/out/abmi_onoff/veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015.Rdata")

## predictor data

lt <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
cn <- intersect(colnames(climCenter_2015), colnames(climSite))
sm <- rbind(climCenter_2015[,cn], climSite[,cn])
rownames(sm) <- sm$Label2
sm$s_y <- paste0(sm$Site, "_", sm$Year)

sm2 <- sm

v1 <- as.matrix(dd1ha$veg_current)
v2 <- as.matrix(dd1ha_2015$veg_current)
rownames(v2) <- sm$Label2[match(rownames(v2), sm$s_y)]
compare_sets(rownames(v1), rownames(v2))
v <- rbind(v1, v2)
compare_sets(rownames(sm), rownames(v))
v <- v[rownames(sm),]
vv <- groupSums(v, 2, lt[colnames(v), "TypeCC"])

cn <- c("Conif", "Decid", "Mixwood", "Pine", "BSpr", "Larch", "CC")
tmp <- rowSums(vv[, cn]) / rowSums(vv)
sm$pforest <- tmp

cn <- c("Conif", "Pine", "BSpr", "Larch")
cn2 <- c("CCConif0", "CCConifR", "CCConif1", "CCConif2", "CCConif3", "CCConif4",
    "CCPine0", "CCPineR", "CCPine1", "CCPine2", "CCPine3", "CCPine4")
tmp <- (rowSums(vv[, cn]) + rowSums(v[, cn2])) / rowSums(vv)
sm$pconif <- tmp

cn <- c("Swamp", "WetGrass", "WetShrub", "BSpr", "Larch")
tmp <- rowSums(vv[, cn]) / rowSums(vv)
sm$pwet <- tmp

cn <- c("Conif5", "Conif6", "Conif7", "Conif8", "Conif9",
    "Decid5", "Decid6", "Decid7", "Decid8", "Decid9",
    "Mixwood5", "Mixwood6", "Mixwood7", "Mixwood8", "Mixwood9",
    "Pine5", "Pine6", "Pine7", "Pine8", "Pine9",
    "Swamp-Conif5", "Swamp-Conif6", "Swamp-Conif7", "Swamp-Conif8", "Swamp-Conif9",
    "Swamp-Decid5", "Swamp-Decid6", "Swamp-Decid7", "Swamp-Decid8", "Swamp-Decid9",
    "Swamp-Mixwood5", "Swamp-Mixwood6", "Swamp-Mixwood7", "Swamp-Mixwood8", "Swamp-Mixwood9",
    "Swamp-Pine5", "Swamp-Pine6", "Swamp-Pine7", "Swamp-Pine8", "Swamp-Pine9",
    "Wetland-BSpr5", "Wetland-BSpr6", "Wetland-BSpr7", "Wetland-BSpr8", "Wetland-BSpr9",
    "Wetland-Decid5", "Wetland-Decid6", "Wetland-Decid7", "Wetland-Decid8", "Wetland-Decid9",
    "Wetland-Larch5", "Wetland-Larch6", "Wetland-Larch7", "Wetland-Larch8", "Wetland-Larch9")
tmp <- rowSums(v[, cn]) / rowSums(v)
sm$poldfor <- tmp

tmp <- vv[, c("Cult", "UrbInd", "SoftLin", "HardLin", "CC")] / rowSums(vv)
sm <- data.frame(sm, tmp)
rownames(sm) <- sm$s_y
sm$Alien <- sm$Cult + sm$UrbInd + sm$HardLin
sm$Succ <- sm$CC + sm$SoftLin
sm$THF <- sm$Alien + sm$Succ

#plot(sm[,c("Cult", "UrbInd", "SoftLin", "HardLin", "CC", "pforest", "pconif", "pwet", "poldfor")])

## taxa

#### vplants
load("e:/peter/AB_data_v2016/data/ermias/Compiled species data_to PS/VPlants.RData")
slt <- read.csv("~/repos/abmispecies/_data/vplants.csv")
rownames(slt) <- slt$sppid
spp <- as.character(slt$sppid[slt$map.det])
#colnames(d)[colnames(d)=="Hemerocallis..Stella.de.Oro."] <- "Hemerocallis.Stella.de.Oro."
d2 <- dd
rownames(d2) <- d2$SiteYear
compare_sets(rownames(d2), sm$s_y)
yy <- d2[intersect(rownames(d2), sm$s_y), spp]
tax <- slt[,c("sppid","scinam","rank","nonnative","species","tsnid")]
vplants <- list(y=as.matrix(yy), tax=tax)

#### mites
load("e:/peter/AB_data_v2016/data/ermias/Compiled species data_to PS/Mites.RData")
slt <- read.csv("~/repos/abmispecies/_data/mites.csv")
rownames(slt) <- slt$sppid
spp <- as.character(slt$sppid[slt$map.det])
d2 <- dd
rownames(d2) <- d2$SiteYear
compare_sets(rownames(d2), sm$s_y)
yy <- d2[intersect(rownames(d2), sm$s_y),spp]
tax <- slt[,c("sppid","scinam","rank","species","tsnid")]
mites <- list(y=as.matrix(yy), tax=tax)

#### mosses
load("e:/peter/AB_data_v2016/data/ermias/Compiled species data_to PS/Moss.RData")
slt <- read.csv("~/repos/abmispecies/_data/mosses.csv")
rownames(slt) <- slt$sppid
spp <- as.character(slt$sppid[slt$map.det])
d2 <- dd
rownames(d2) <- d2$SiteYear
compare_sets(rownames(d2), sm$s_y)
yy <- d2[intersect(rownames(d2), sm$s_y),spp]
tax <- slt[,c("sppid","scinam","species","tsnid")]
mosses <- list(y=as.matrix(yy), tax=tax)

#### lichens
load("e:/peter/AB_data_v2016/data/ermias/Compiled species data_to PS/Lichens.RData")
slt <- read.csv("~/repos/abmispecies/_data/lichens.csv")
rownames(slt) <- slt$sppid
spp <- as.character(slt$sppid[slt$map.det])
d2 <- dd
rownames(d2) <- d2$SiteYear
compare_sets(rownames(d2), sm$s_y)
yy <- d2[intersect(rownames(d2), sm$s_y),spp]
tax <- slt[,c("sppid","scinam","species","tsnid")]
lichens <- list(y=as.matrix(yy), tax=tax)

#### birds
load("e:/peter/AB_data_v2016/data/species/OUT_birdsrf_2016-05-27.Rdata")
slt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(slt) <- slt$sppid
spp <- as.character(slt$sppid[slt$map.det])
mm <- m[samp(m)$TBB_POINT_COUNT == 1,]
rownames(mm) <- paste0(samp(mm)$SiteLabel, "_", samp(mm)$YEAR)
d2 <- as.matrix(xtab(mm))
compare_sets(rownames(d2), sm$s_y)
spp <- intersect(spp, colnames(d2))
yy <- d2[intersect(rownames(d2), sm$s_y),spp]
tax <- slt[,c("sppid","scinam","family","species","singing","tsnid")]
birds <- list(y=as.matrix(yy), tax=tax)


## birds vegHF piece: dd150m, dd1km
#load("~/Dropbox/collaborations/opticut/R/abmi-data/abmi-data.Rdata")
load(file.path("e:/peter/AB_data_v2016", "out", "birds", "data", "data-full-veghf.Rdata"))

v <- as.matrix(dd150m$veg_current)
rn <- rownames(v)
tmp <- strsplit(rn, ":")
table(sapply(tmp, length)) # 1: ABMI RF, 3: ABMI SM
rn <- rn[sapply(tmp, length) == 1]
tmp <- data.frame(do.call(rbind, strsplit(rn, "_")))
rownames(tmp) <- rn
tmp$s_y <- paste0(tmp$X4, "_", tmp$X5)
tmp <- tmp[tmp$X8 == "1",] # center point
v <- v[rownames(tmp),]
rownames(v) <- tmp$s_y
compare_sets(tmp$s_y, sm2$s_y)
rownames(sm2) <- sm2$s_y
compare_sets(rownames(v), rownames(sm2))
isect <- intersect(rownames(v), rownames(sm2))
sm2 <- sm2[isect,]

v <- v[rownames(sm2),]
vv <- groupSums(v, 2, lt[colnames(v), "TypeCC"])

cn <- c("Conif", "Decid", "Mixwood", "Pine", "BSpr", "Larch", "CC")
tmp <- rowSums(vv[, cn]) / rowSums(vv)
sm2$pforest <- tmp

cn <- c("Conif", "Pine", "BSpr", "Larch")
cn2 <- c("CCConif0", "CCConifR", "CCConif1", "CCConif2", "CCConif3", "CCConif4",
    "CCPine0", "CCPineR", "CCPine1", "CCPine2", "CCPine3", "CCPine4")
tmp <- (rowSums(vv[, cn]) + rowSums(v[, cn2])) / rowSums(vv)
sm2$pconif <- tmp

cn <- c("Swamp", "WetGrass", "WetShrub", "BSpr", "Larch")
tmp <- rowSums(vv[, cn]) / rowSums(vv)
sm2$pwet <- tmp

cn <- c("Conif5", "Conif6", "Conif7", "Conif8", "Conif9",
    "Decid5", "Decid6", "Decid7", "Decid8", "Decid9",
    "Mixwood5", "Mixwood6", "Mixwood7", "Mixwood8", "Mixwood9",
    "Pine5", "Pine6", "Pine7", "Pine8", "Pine9",
    "Swamp-Conif5", "Swamp-Conif6", "Swamp-Conif7", "Swamp-Conif8", "Swamp-Conif9",
    "Swamp-Decid5", "Swamp-Decid6", "Swamp-Decid7", "Swamp-Decid8", "Swamp-Decid9",
    "Swamp-Mixwood5", "Swamp-Mixwood6", "Swamp-Mixwood7", "Swamp-Mixwood8", "Swamp-Mixwood9",
    "Swamp-Pine5", "Swamp-Pine6", "Swamp-Pine7", "Swamp-Pine8", "Swamp-Pine9",
    "Wetland-BSpr5", "Wetland-BSpr6", "Wetland-BSpr7", "Wetland-BSpr8", "Wetland-BSpr9",
    "Wetland-Decid5", "Wetland-Decid6", "Wetland-Decid7", "Wetland-Decid8", "Wetland-Decid9",
    "Wetland-Larch5", "Wetland-Larch6", "Wetland-Larch7", "Wetland-Larch8", "Wetland-Larch9")
tmp <- rowSums(v[, cn]) / rowSums(v)
sm2$poldfor <- tmp

tmp <- vv[, c("Cult", "UrbInd", "SoftLin", "HardLin", "CC")] / rowSums(vv)
sm2 <- data.frame(sm2, tmp)
sm2$Alien <- sm2$Cult + sm2$UrbInd + sm2$HardLin
sm2$Succ <- sm2$CC + sm2$SoftLin
sm2$THF <- sm2$Alien + sm2$Succ

ABMI <- list(sites=sm, sites_pc=sm2,
    species=list(birds=birds$tax, mites=mites$tax,
        vascular_plants=vplants$tax, bryophytes=mosses$tax, lichens=lichens$tax),
    detections=list(birds=birds$y, mites=mites$y,
        vascular_plants=vplants$y, bryophytes=mosses$y, lichens=lichens$y))

save(ABMI, file="~/Dropbox/collaborations/opticut/R/abmi-data/abmi-data-AB.Rdata")

## -------- raw data based ordination as per request --------------

library(mefa4)
load("~/Dropbox/collaborations/opticut/R/abmi-data/abmi-data-AB.Rdata")

lapply(lapply(ABMI$detections, rownames), function(z) table(ABMI$sites[z,"NRNAME"]))

ii <- sort(intersect(rownames(ABMI$sites),
    Reduce(intersect, lapply(ABMI$detections, rownames))))
for (j in names(ABMI$species))
    ABMI$species[[j]]$taxon <- j
jj <- Reduce(intersect, lapply(ABMI$species, colnames))
z <- do.call(rbind, lapply(ABMI$species, function(z) z[,jj]))
rownames(z) <- z$sppid
m <- Mefa(do.call(cbind, lapply(ABMI$detections, function(z) z[ii,])),
    ABMI$sites[ii,],
    z,
    join="inner")
xtab(m)[xtab(m) > 0] <- 1
m <- m[rowSums(xtab(m)) > 0, colSums(xtab(m)) > 0]

nmin <- 5
m <- m[, colSums(xtab(m)) >= nmin]
m <- m[rowSums(xtab(m)) > 0,]

library(vegan)

Y <- as.matrix(xtab(m))
X <- samp(m)
Z <- taxa(m)
Z$freq <- colSums(Y)
Z$p <- colSums(Y) / nrow(Y)
Z$taxon <- as.factor(Z$taxon)

m1 <- rda(Y ~ pforest + pconif + pwet + poldfor +
    AHM + PET + MWMT + MCMT +
    Cult + UrbInd + SoftLin + HardLin + CC, X)

m2 <- cca(Y ~ pforest + pconif + pwet + poldfor +
    AHM + PET + MWMT + MCMT +
    Cult + UrbInd + SoftLin + HardLin + CC, X)

op <- ordiplot(m1, type="n")
points(op, "species")
points(op, "biplot", col=2)

pal <- colorRampPalette(c("red","blue"))(100)
plot(op$species, pch=c(3,4,5,6)[Z$taxon], col=pal[ceiling(100*sqrt(Z$p))])

sc <- scores(m1)

m3 <- metaMDS(Y)

