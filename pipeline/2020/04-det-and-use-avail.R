## detection maps and use-availability figures

library(mefa4)
library(raster)
make_raster <- function(value, rc, rt) {
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

## detection maps for non birds

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid


m <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

rnr <- make_raster(as.integer(kgrid$NRNAME), kgrid, rt)
cnr <- c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae')
cnr <- cnr[c(5,6,1,2,4,3)]

f <- paste0(
    "s:/AB_data_v2020/Results/Results from Ermias/species data for detection map/",
    c("Lichen", "Mites", "Moss", "VPlants"),
    "_Site_Compiled species data.RData")
names(f) <- c("lichens", "mites", "mosses", "vplants")
# no det maps for nnplants

#taxon <- "lichens"
for (taxon in names(f)) {
    ## models for these species
    SPPn <- dimnames(COEFS[[taxon]]$north)[[1]]
    SPPs <- dimnames(COEFS[[taxon]]$south)[[1]]
    SPP <- sort(union(SPPn, SPPs))

    ## data for these species
    ex <- new.env()
    load(f[taxon], envir=ex)
    d <- ex[[names(ex)[1]]]
    d <- d[!startsWith(rownames(d), "OG-ALPAC-SK"),]
    d0 <- d[,1:5]
    yy <- d[6:ncol(d)]
    SPP2 <- colnames(yy)
    SPPnew <- setdiff(SPP2, SPP)

    site <- d0$Site
    og <- d0$OnOffGrid == "OG"
    ogs <- sapply(strsplit(site[og], "-"), "[[", 3)
    site[og] <- ogs
    site <- gsub("B", "", site)
    siten <- as.integer(site)
    #site[is.na(siten)]

    xy <- data.frame(x=m$PUBLIC_LONGITUDE, y=m$PUBLIC_LATTITUDE)[match(site, m$SITE_ID),]
    coordinates(xy) <- ~x+y
    proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    xy <- spTransform(xy, proj4string(rt))

    for (spp in SPP2) {
        cat(taxon, spp, "\n")
        flush.console()

        DIR <- file.path("s:/AB_data_v2020/Results/web", taxon, spp)
        if (!dir.exists(DIR))
            dir.create(DIR)

        xy0 <- xy[yy[,spp] == 0 & !duplicated(site),]
        xy1 <- xy[yy[,spp] > 0,]
        png(file.path(DIR, "det.png"),
            height=1500*1.5, width=1000*1.5, res=300)
        op <- par(mar=c(1,1,1,1))
        plot(rnr,col=cnr, axes=FALSE, box=FALSE, main=spp, legend=FALSE)
        plot(xy0, add=TRUE, pch=19, col="#aaaaaa88", legend=FALSE, cex=0.8)
        plot(xy1, add=TRUE, pch=19, col="red4", legend=FALSE, cex=0.8)
        par(op)
        dev.off()

    }


}

## detection maps for birds

ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds" # change this bit

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2020-09-23.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2020-09-23.RData"), envir=es)


ddd <- ee$dd
## subset based on analysis data -- same filtering applied
ddd <- ddd[unique(c(rownames(en$DAT), rownames(es$DAT))),]
ddd <- nonDuplicated(ddd, ddd$SS, TRUE)
yyy <- groupSums(ee$yy, 1, ee$dd$SS)[rownames(ddd),]
yyy[yyy > 0] <- 1
ss <- !is.na(ddd$X) & !is.na(ddd$NRNAME)
ddd <- ddd[ss,]
yyy <- yyy[ss,]
yyy <- yyy[,colSums(yyy) > 0]

xy <- SpatialPoints(as.matrix(ddd[,c("X","Y")]))
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))
rt10 <- aggregate(rt, fact=10)
sam0 <- rasterize(xy, rt10, field=1, fun='sum')
values(sam0)[!is.na(values(sam0))] <- 1

taxon <- "birds"
SPPn <- dimnames(COEFS[[taxon]]$north$joint)[[1]]
SPPs <- dimnames(COEFS[[taxon]]$south$joint)[[1]]
SPP <- sort(union(SPPn, SPPs))
SPPall <- rownames(COEFS$birds$species)
tax <- ee$tax

colnames(yyy) <- tax$sppid[match(colnames(yyy), tax$code)]
SPP2 <- colnames(yyy)
compare_sets(SPPall, SPP2)
SPP2 <- SPPall

for (spp in SPP2) {
    cat(taxon, spp, "\n");flush.console()
    xy1 <- SpatialPoints(as.matrix(ddd[yyy[,spp] > 0,c("X","Y")]))
    proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    xy1 <- spTransform(xy1, proj4string(rt))
    sam1 <- rasterize(xy1, rt10, field=1, fun='last')

    DIR <- file.path("s:/AB_data_v2020/Results/web", taxon, spp)
    if (!dir.exists(DIR))
        dir.create(DIR)

    png(file.path(DIR, "det.png"),
        height=1500*1.5, width=1000*1.5, res=300)
    op <- par(mar=c(1,1,1,1))
    plot(rnr,col=cnr, axes=FALSE, box=FALSE, main=spp, legend=FALSE)
    plot(sam0,add=TRUE, col="#aaaaaa88", legend=FALSE)
    plot(sam1,add=TRUE, col="red4", legend=FALSE)
    par(op)
    dev.off()
}

## use-avail for birds

library(mefa4)
source("~/repos/abmianalytics/birds/00-functions.R")

ls <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(ls) <- ls[,1]
lv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v2020.csv")
rownames(lv) <- lv[,1]

Cn <- c("Decid", "Mixedwood", "Pine", "Spruce", "TreedBogFen", "Swamp",
    "GrassHerb", "Shrub", "NonTreedBogFen", "Marsh",
    "Crop", "RoughP", "TameP", "UrbInd", "SoftLin",
    "HardLin", "Forestry")
Cs <- c("ClayLoamSand", "RapidDrain", "Blowout", "ThinBreak", "Other",
     "Crop", "RoughP", "TameP", "UrbInd", "SoftLin", "HardLin")

ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds" # change this bit

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2020-09-23.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2020-09-23.RData"), envir=es)

Ys <- es$YY
As <- ee$sc1[rownames(Ys),]
Yn <- en$YY
An <- ee$vc1[rownames(Yn),]

tax <- ee$tax
colnames(Ys) <- tax$sppid[match(colnames(Ys), tax$code)]
colnames(Yn) <- tax$sppid[match(colnames(Yn), tax$code)]

spt <- COEFS$birds$species
with(spt, table(ModelNorth, rownames(spt) %in% colnames(Yn)))
with(spt, table(ModelSouth, rownames(spt) %in% colnames(Ys)))
SPPsw <- rownames(spt)[rownames(spt) %in% colnames(Ys) &
        !spt$ModelSouth]
SPPnw <- rownames(spt)[rownames(spt) %in% colnames(Yn) &
        !spt$ModelNorth]
Ys <- Ys[,SPPsw]
Yn <- Yn[,SPPnw]

compare_sets(ls$ID, colnames(As))
compare_sets(lv$ID, colnames(An))

pn <- groupSums(An, 2, lv[colnames(An), "UseAvail"])
pn <- pn[,colnames(pn) != "EXCLUDE"]
pn <- pn / rowSums(pn)
pn <- pn[,Cn]
ps <- groupSums(As, 2, ls[colnames(As), "UseAvail"])
ps <- ps[,colnames(ps) != "EXCLUDE"]
ps <- ps / rowSums(ps)
ps <- ps[,Cs]

wrs <- t(pbapply::pbapply(as.matrix(Ys), 2, function(z) wrsi(z, x=ps)$rWRSI))
colnames(wrs) <- rownames(wrsi(Ys[,1], x=ps))
wrn <- t(pbapply::pbapply(as.matrix(Yn), 2, function(z) wrsi(z, x=pn)$rWRSI))
colnames(wrn) <- rownames(wrsi(Yn[,1], x=pn))


## Use avail for other taxa

Taxa <- c(lichens="Lichen", mites="Mite", mosses="Moss", vplants="VPlant")
UA <- list()
for (taxa in names(Taxa)) {
    f <- sprintf(
        paste0(
        "s:/AB_data_v2020/Results/Results from Ermias/Compiled data/",
        "Compiled %s %s analysis 2020_Quad.RData"), "South", Taxa[taxa])

    e <- new.env()
    load(f, envir=e)
    names(e)
    d <- e$d
    ys <- d[,e$SpTable.ua]
    v <- d[,which(colnames(d) == "ClaySub"):ncol(d)]
    v <- v[,!startsWith(colnames(v), "Surr")]
    tmp <- intersect(colnames(v),ls$UseInAnalysis)
    vv1 <- groupSums(as.matrix(v[,tmp]), 2, ls$UseAvail[match(tmp, ls$UseInAnalysis)])
    setdiff(ls$UseAvail, colnames(vv1))
    vv1 <- cbind(vv1, UrbInd=v[,"UrbInd"])
    vv1 <- vv1[,colnames(vv1) != "EXCLUDE"]
    psx <- vv1 / rowSums(vv1)
    psx <- psx[,colnames(ps)]


    f <- sprintf(
        paste0(
        "s:/AB_data_v2020/Results/Results from Ermias/Compiled data/",
        "Compiled %s %s analysis 2020_Quad.RData"), "North", Taxa[taxa])
    e <- new.env()
    load(f, envir=e)
    names(e)
    d <- e$d
    yn <- d[,e$SpTable.ua]
    v <- d[,which(colnames(d) == "DeciduousR"):ncol(d)]
    v <- v[,!startsWith(colnames(v), "Surr")]
    colnames(v) <- gsub("Deciduous", "Decid", colnames(v))
    colnames(v) <- gsub("WhiteSpruce", "Spruce", colnames(v))
    colnames(v) <- gsub("BlackSpruce", "TreedBog", colnames(v))
    colnames(v) <- gsub("TreedFen", "TreedFenR", colnames(v))
    colnames(v) <- gsub("TreedSwamp", "TreedSwampR", colnames(v))
    colnames(v) <- gsub("Grass", "GrassHerb", colnames(v))
    tmp <- intersect(colnames(v),lv$UseInAnalysis)
    vv1 <- groupSums(as.matrix(v[,tmp]), 2, lv$UseAvail[match(tmp, lv$UseInAnalysis)])
    setdiff(lv$UseAvail, colnames(vv1))
    vv1 <- cbind(vv1, UrbInd=v[,"UrbInd"])
    vv1 <- vv1[,colnames(vv1) != "EXCLUDE"]
    pnx <- vv1 / rowSums(vv1)
    pnx <- pnx[,colnames(pn)]

    pnx[is.na(pnx)] <- 0
    psx[is.na(psx)] <- 0

    spt <- COEFS[[taxa]]$species
    SPPnw <- rownames(spt)[spt$UseavailNorth]
    SPPsw <- rownames(spt)[spt$UseavailSouth]
    yn <- yn[,SPPnw]
    ys <- ys[,SPPsw]

    N <- t(pbapply::pbapply(as.matrix(yn), 2, function(z) wrsi(z, x=pnx)$rWRSI))
    colnames(N) <- rownames(wrsi(yn[,1], x=pnx))
    S <- t(pbapply::pbapply(as.matrix(ys), 2, function(z) wrsi(z, x=psx)$rWRSI))
    colnames(S) <- rownames(wrsi(ys[,1], x=psx))


    UA[[taxa]] <- list(north=N, south=S)

}

UA$birds <- list(north=wrn, south=wrs)

## use-avail for mammals

load("s:/AB_data_v2020/Results/UseAvail-all.RData")

load("s:/AB_data_v2020/Results/Camera mammal models revised June 2020/North/UseavailNorth.Rdata")
load("s:/AB_data_v2020/Results/Camera mammal models revised June 2020/South/UseavailSouth.Rdata")

cn <- colnames(UA$lichens$north)
mefa4::compare_sets(colnames(UseavailNorth), cn)

CN <- c(Decid="Deciduous",
    Mixedwood="Mixedwood",
    Pine="Pine",
    Spruce="WhiteSpruce",
    TreedBogFen="BlackSpruce", # BlackSpruce + TreedFen
    Swamp="Wetland",
    GrassHerb="Open",
    Shrub="Open",
    NonTreedBogFen="Wetland",
    Marsh="Wetland",
    Crop="Crop",
    RoughP="RoughP",
    TameP="TameP",
    UrbInd="RurUrbInd",
    SoftLin="SoftLin",
    HardLin="HardLin",
    Forestry="HFor")
xn <- UseavailNorth[,CN]
colnames(xn) <- names(CN)
xn[,"TreedBogFen"] <-  0.5 * (UseavailNorth[,"BlackSpruce"] + UseavailNorth[,"TreedFen"])
xn[,"UrbInd"] <-  0.8 * UseavailNorth[,"RurUrbInd"] + 0.2 * UseavailNorth[,"Well"]

cs <- colnames(UA$lichens$south)
mefa4::compare_sets(colnames(UseavailSouth), cs)
setdiff(colnames(UseavailSouth), cs)

CS <- c(ClayLoamSand="Clay",
    RapidDrain="RapidDrain",
    Blowout="Saline",
    ThinBreak="Saline",
    Other="Productive",
    Crop="Crop",
    RoughP="RoughP",
    TameP="TameP",
    UrbInd="RurUrbInd",
    SoftLin="SoftLin",
    HardLin="HardLin")
xs <- UseavailSouth[,CS]
colnames(xs) <- names(CS)
xs[,"UrbInd"] <-  0.8 * UseavailSouth[,"RurUrbInd"] + 0.2 * UseavailSouth[,"Well"]
xs[,"ThinBreak"] <-  0

load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")

spt <- COEFS2$mammals$species
SPPnw <- rownames(spt)[spt$UseavailNorth & !spt$ModelNorth]
SPPsw <- rownames(spt)[spt$UseavailSouth & !spt$ModelSouth]
xn <- xn[SPPnw,]
xs <- xs[SPPsw,]


UA$mammals <- list(north=xn, south=xs)

save(UA, file="s:/AB_data_v2020/Results/UseAvail-all.RData")


## use-avail plots N/S

library(RColorBrewer)

load("s:/AB_data_v2020/Results/UseAvail-all.RData")

col1<-brewer.pal(8, "Dark2")[c(1,1,1,1,1, 6, 7,7)]
col2<-brewer.pal(12, "Paired")[c(12,7,7)]
cols <- c(col1,col2)
names(cols) <- Cs
colsS <- cols

col1<-brewer.pal(8, "Dark2")[c(1,1,1,1,5,5,2,2,3,3, 6, 7,7)]
col2<-brewer.pal(12, "Paired")[c(12,7,7,3)]
cols <- c(col1,col2)
names(cols) <- Cn
colsN <- cols

for (taxon in names(UA)) {

    Z <- UA[[taxon]]$north
    for (spp in rownames(Z)) {
        cat(taxon, spp, "north\n")
        flush.console()

        DIR <- file.path("s:/AB_data_v2020/Results/web", taxon, spp)
        if (!dir.exists(DIR))
            dir.create(DIR)

        v <- Z[spp,names(colsN)]
#        unlink(file.path(DIR, "useavail-north.svg"))
#        svglite::svglite(file.path(DIR, "useavail-north.svg"),
#            height=4, width=7)
        png(file.path(DIR, "useavail-north.png"),
            height=480, width=600)
        op <- par(mar=c(7,4,2,2)+0.1, las=2)
        x1 <- barplot(v, horiz=FALSE, ylab="Affinity",space=NULL, col=cols, border=cols, ylim=c(-1,1), axes=FALSE,axisnames=F )
        axis(side=2)
        abline(h=0, col="red4", lwd=2)
        mtext(side=3,at=x1[1],adj=0,spp,cex=1.2,col="grey40",las=1)
        text(x=x1, y=par()$usr[3]-0.01,labels=names(colsN), srt=60, adj=1, col=cols, xpd=TRUE)
        par(op)
        dev.off()
    }

    Z <- UA[[taxon]]$south
    for (spp in rownames(Z)) {
        cat(taxon, spp, "south\n")
        flush.console()

        DIR <- file.path("s:/AB_data_v2020/Results/web", taxon, spp)
        if (!dir.exists(DIR))
            dir.create(DIR)

        v <- Z[spp,names(colsS)]
#        unlink(file.path(DIR, "useavail-south.svg"))
#        svglite::svglite(file.path(DIR, "useavail-south.svg"),
#            height=4, width=7)
        png(file.path(DIR, "useavail-south.png"),
            height=480, width=600)
        op <- par(mar=c(7,4,2,2)+0.1, las=2)
        x1 <- barplot(v, horiz=FALSE, ylab="Affinity",space=NULL, col=cols, border=cols, ylim=c(-1,1), axes=FALSE,axisnames=F )
        axis(side=2)
        abline(h=0, col="red4", lwd=2)
        mtext(side=3,at=x1[1],adj=0, spp, cex=1.2,col="grey40",las=1)
        text(x=x1, y=par()$usr[3]-0.01,labels=names(colsS), srt=60, adj=1, col=cols, xpd=TRUE)
        par(op)
        dev.off()
    }
}


## det maps for mammals

l <- list.files("s:/AB_data_v2020/Results/Camera mammal models revised June 2020/Maps/Occurrence")
l <- gsub("\\.png", "", l)

l2 <- list.dirs("s:/AB_data_v2020/Results/web/mammals", full.name=FALSE, recursive=FALSE)
mefa4::compare_sets(l, l2)

for (i in l) {
    if (!dir.exists(paste0("s:/AB_data_v2020/Results/web/mammals/", i)))
        dir.create(paste0("s:/AB_data_v2020/Results/web/mammals/", i))
    file.copy(
        paste0("s:/AB_data_v2020/Results/Camera mammal models revised June 2020/Maps/Occurrence/",
            i, ".png"),
        paste0("s:/AB_data_v2020/Results/web/mammals/", i, "/det.png")
    )
}



## mammal exact prediction maps
library(magick)

load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
taxon <- "mammals"
SPPn <- dimnames(COEFS2[[taxon]]$north)[[1]]
SPPs <- dimnames(COEFS2[[taxon]]$south)[[1]]
SPPs <- SPPs[SPPs != "Fisher"]
SPP <- sort(unique(c(SPPn, SPPs)))

spp="Coyote"

for (spp in SPP) {
f3 <- paste0(
    "s:/AB_data_v2020/Results/Camera mammal models South revised Nov 2020/Combine regions/Maps/",
    c("Reference/", "Current/", "Difference/"),
    spp, ".png")

i3 <- image_append( image_read(f3))
image_write(i3, paste0("s:/AB_data_v2020/Results/web/mammals/", spp, "/map2.png"))
}


# updated south figures
for (spp in SPPs) {
    cat(spp, "\n")
    f1 <- paste0(
        "s:/AB_data_v2020/Results/Camera mammal models South revised Nov 2020/",
        "South/Figures/Best model/",
        c("Non-treed", "Treed"), "/Soil+HF figure best model ", spp, ".png")

    i1 <- image_append( image_read(f1), stack = TRUE)
    image_write(i1,
        paste0("s:/AB_data_v2020/Results/web/mammals/", spp, "/soilhf2.png"))

    f2 <- paste0(
        "s:/AB_data_v2020/Results/Camera mammal models South revised Nov 2020/",
        "South/Maps/",
        c("Reference", "Current", "Difference"),
        "/", spp, ".jpg")

    i2 <- image_append( image_read(f2), stack = FALSE)
    image_write(i2,
        paste0("s:/AB_data_v2020/Results/web/mammals/", spp, "/map2s.png"))

}
