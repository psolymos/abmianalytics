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
tmp <- vv / rowSums(vv)
sm <- data.frame(sm, tmp)
sm$XXX <- NULL

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
#sm <- data.frame(sm, tmp)
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

ABMI <- list(sites=sm,
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
taxa(m)$taxon <- as.factor(taxa(m)$taxon)

library(vegan)

mm <- m
#mm <- m[,taxa(m)$taxon=="vascular_plants"]

Y <- as.matrix(xtab(mm))
X <- samp(mm)
Z <- taxa(mm)
Z$freq <- colSums(Y)
Z$p <- colSums(Y) / nrow(Y)
#br <- quantile(Z$freq, seq(0,1,0.1))
br <- c(min(Z$freq),10,15,20,25,30,40,50,75,100,200,500,max(Z$freq))
q <- cut(Z$freq, br, labels=FALSE, include.lowest=TRUE)
X$MAP <- gsub(",", "", X$MAP)
X$MAP <- as.integer(X$MAP)

cn <- c("MWMT", "MCMT", "MAT", "MAP", "FFP", "PET",
    "Conif", "Decid", "Mixwood", "Pine", "GrassHerb",
    "Shrub", "Swamp", "WetGrass", "WetShrub", "BSpr", "Larch",
    "Cult", "UrbInd", "SoftLin", "HardLin", "CC", "poldfor")
mod <- cca(Y ~ ., X[,cn])
odp <- ordiplot(mod, type="n")

library(KernSmooth)
bw4bkde <- function (x, kernel = "normal", canonical = FALSE, bandwidth,
gridsize = 401L, range.x, truncate = TRUE) {
    if (!missing(bandwidth) && bandwidth <= 0)
        stop("'bandwidth' must be strictly positive")
    kernel <- match.arg(kernel, c("normal", "box", "epanech",
        "biweight", "triweight"))
    n <- length(x)
    del0 <- switch(kernel, normal = (1/(4 * pi))^(1/10), box = (9/2)^(1/5),
        epanech = 15^(1/5), biweight = 35^(1/5), triweight = (9450/143)^(1/5))
    h <- if (missing(bandwidth))
        del0 * (243/(35 * n))^(1/5) * sqrt(var(x))
    else if (canonical)
        del0 * bandwidth
    else bandwidth
    h
}
quantiles4bkde <- function(x, y) {
    bw <- c(x=bw4bkde(x), y=bw4bkde(y))
    d2 <- bkde2D(cbind(x, y), bandwidth=bw)
    o <- order(d2$fhat)
    i <- order(o)
    f <- d2$fhat
    f <- f/sum(f)
    cs <- cumsum(f[o])[i]
    dim(cs) <- dim(d2$fhat)
    d2$cs <- cs
    d2$bw <- bw
    d2
}

Biplot <- function(x, y, col=1, type="d", contour=TRUE, ...) {
    Pal <- colorRampPalette(c("#FFFFFF", col))(100)
    d <- quantiles4bkde(x, y)
    plot(x, y, type="n", xlim=range(d$x1), ylim=range(d$x2), axes=FALSE, ...)
    u <- par("usr")
    image(d$x1, d$x2, d$fhat, col=Pal[1:66], add=TRUE)
    #points(x, y, col=paste0(Pal[50], "80"), pch=".")
    if (contour)
        contour(d$x1, d$x2, d$cs, add=TRUE, col=Pal[100], levels=c(0.25, 0.5, 0.75))
    box(col=Pal[33])
    axis(1, col=Pal[100])
    axis(2, col=Pal[100])
    invisible(NULL)
}

pal <- rev(viridis::viridis(max(q)))
#pal <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba'))(max(q))
v <- 1 * max(abs(odp$species))/max(abs(odp$biplot))
c("MWMT", "MCMT", "MAT", "MAP", "FFP", "PET", "Conif", "Decid",
"Mixwood", "Pine", "GrassHerb", "Shrub", "Swamp", "WetGrass",
"WetShrub", "BSpr", "Larch", "Cult", "UrbInd", "SoftLin", "HardLin",
"CC", "poldfor")

col <- rep(2, nrow(odp$biplot))
col[rownames(odp$biplot) %in% c("MWMT", "MCMT", "MAT", "MAP", "FFP", "PET")] <- 3
col[rownames(odp$biplot) %in% c("Cult", "UrbInd", "SoftLin", "HardLin", "CC")] <- 4


op <- par(mfrow=c(2,3))

Biplot(odp$species[,1], odp$species[,2], xlab="CCA1", ylab="CCA2", main="All")
abline(h=0, v=0, lty=2)
tmp <- sapply(1:nrow(odp$biplot), function(i)
    arrows(x0=0, y0=0, x1=v*odp$biplot[i,1], y1=v*odp$biplot[i,2], col=col[i],
    lwd=1, length = 0.1, angle = 15))
par(xpd=TRUE)
tmp <- sapply(1:nrow(odp$biplot), function(i)
    text(1.1*v*odp$biplot[i,1], 1.1*v*odp$biplot[i,2], rownames(odp$biplot)[i],
    col=col[i]))
par(xpd=FALSE)
legend("topright", bty="n", fill=pal, legend=paste0(br[-length(br)], "-", br[-1]),
    title="# detections")

for (tt in levels(Z$taxon)) {#tt <- "birds"
    ss <- Z$taxon == tt
    Biplot(odp$species[,1], odp$species[,2], xlab="CCA1", ylab="CCA2", main=tt, contour=FALSE)
    abline(h=0, v=0, lty=2)
    d <- quantiles4bkde(odp$species[ss,1], odp$species[ss,2])
    points(odp$species[ss,], col=pal[q[ss]], pch=19)
    contour(d$x1, d$x2, d$cs, add=TRUE, col=1, levels=c(0.25, 0.5, 0.75))
}
par(op)

