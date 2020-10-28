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

tax <- ee$tax

colnames(yyy) <- tax$sppid[match(colnames(yyy), tax$code)]
SPP2 <- colnames(yyy)

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




load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")

library(RColorBrewer)


tab <- OUT$Species
uan <- OUT$UseavailNorth
uas <- OUT$UseavailSouth

x <- as.matrix(uan[ ,c("Deciduous","Mixedwood","WhiteSpruce","Pine","BlackSpruce","TreedFen","Open","Wetland","HFor","Crop", "TameP", "RoughP","UrbInd","HardLin","SoftLin")])
HabLabel <- c("Deciduous","Mixedwood","Upland Spruce","Pine","Black Spruce","Treed Fen","Open Upland","Open Wetland","Forestry","Crop", "Tame Pasture", "Rough Pasture","Urban/Industry","Hard Linear","Soft Linear" )
col1<-brewer.pal(8, "Dark2")[c(1,1,1,1, 5,5, 6,7)]
col2<-brewer.pal(12, "Paired")[c(4,7,7,7,12,12,10)]
cols <- c(col1,col2)

#SPP <- rownames(tab)[tab$UseavailNorth]
SPP <- rownames(tab)[tab$UseavailNorth & tab$Group == "birds"]

for (spp in SPP) {
    gr <- tab[spp, "Group"]
    spnam <- if (is.na(tab[spp, "CommonName"])) {
        as.character(tab[spp, "ScientificName"])
    } else {
        paste0(as.character(tab[spp, "CommonName"]), " (", as.character(tab[spp, "ScientificName"]), ")")
    }
    cat(gr, spp, "\n");flush.console()
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-useavail-north.png"),
        height=480, width=600)
    op <- par(mar=c(6,4,2,2)+0.1, las=2)
    x1 <- barplot(as.vector(x [spp, ]), horiz=FALSE, ylab="Affinity",space=NULL, col=cols, border=cols, ylim=c(-1,1), axes=FALSE,axisnames=F )
    axis(side=2)
    abline(h=0, col="red4", lwd=2)
    mtext(side=3,at=x1[1],adj=0, spnam, cex=1.2,col="grey40",las=1)
    text(x=x1, y=par()$usr[3]-0.01,labels=HabLabel, srt=60, adj=1, col=cols, xpd=TRUE)
    par(op)
    dev.off()
}

x<- as.matrix(uas[ , c("Productive","Clay","Saline","RapidDrain","Crop","TameP","RoughP","UrbInd","HardLin","SoftLin")])
HabLabel <- c("Productive","Clay","Saline","Rapid Drain","Crop", "Tame Pasture", "Rough Pasture","Urban/Industry","Hard Linear","Soft Linear" )
col1<-brewer.pal(8, "Dark2")[c(7,7,7,7)]
col2<-brewer.pal(12, "Paired")[c(7,7,7,12,12,10)]
cols <- c(col1,col2)

#SPP <- rownames(tab)[tab$UseavailSouth]
SPP <- rownames(tab)[tab$UseavailSouth & tab$Group == "birds"]

for (spp in SPP) {
    gr <- tab[spp, "Group"]
    spnam <- if (is.na(tab[spp, "CommonName"])) {
        as.character(tab[spp, "ScientificName"])
    } else {
        paste0(as.character(tab[spp, "CommonName"]), " (", as.character(tab[spp, "ScientificName"]), ")")
    }
    cat(gr, spp, "\n");flush.console()
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-useavail-south.png"),
        height=480, width=600)
    op <- par(mar=c(6,4,2,2)+0.1, las=2)
    x1 <- barplot(as.vector(x[spp, ]), horiz=FALSE, ylab="Affinity",space=NULL, col=cols, border=cols, ylim=c(-1,1), axes=FALSE,axisnames=F )
    axis(side=2)
    abline(h=0, col="red4", lwd=2)
    mtext(side=3,at=x1[1],adj=0,spnam,cex=1.2,col="grey40",las=1)
    text(x=x1, y=par()$usr[3]-0.01,labels=HabLabel, srt=60, adj=1, col=cols, xpd=TRUE)
    par(op)
    dev.off()
}


