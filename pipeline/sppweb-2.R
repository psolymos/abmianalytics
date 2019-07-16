#devtools::install_github("ABbiodiversity/cure4insect@v2018")

library(cure4insect)
library(mefa4)
#set_options(path = "s:/reports")
set_options(path = "d:/abmi/reports")
load_common_data()

load("d:/abmi/sppweb2018/c4i/tables/sector-effects.RData")
load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")
ROOT <- "d:/abmi/AB_data_v2018/www"
#ROOT <- "d:/abmi/sppweb2018/www/"

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

Rmaskn <- make_raster(as.integer(1-kgrid$useS), kgrid, rt)
values(Rmaskn)[values(Rmaskn) == 0] <- NA
Rmasks <- make_raster(as.integer(1-kgrid$useN), kgrid, rt)
values(Rmasks)[values(Rmasks) == 0] <- NA
Rmaskm <- make_raster(as.integer(!kgrid$NRNAME == "Rocky Mountain"), kgrid, rt)
values(Rmaskm)[values(Rmaskm) == 0] <- NA
Rw <- make_raster(as.integer(kgrid$pWater > 0.99), kgrid, rt)
values(Rw)[values(Rw) == 0] <- NA

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
col3 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221"))(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude


for (gr in c("birds","vplants","lichens","mosses","mites")) {
#gr <- "birds"
#gr <- "vplants"
#gr <- "lichens"
#gr <- "mosses"
#gr <- "mites"

#SPP <- rownames(resn[resn$Taxon == gr,])
SPP <- rownames(OUT$Species[(OUT$Species$ModelNorth |
        OUT$Species$ModelSouth) & OUT$Species$Group == gr,])
#spp <- "AlderFlycatcher"

for (spp in SPP) {

    cat(gr, spp, "\n");flush.console()

    TYPE <- "C"
    if (OUT$Species[spp, "ModelNorth"] && !OUT$Species[spp, "ModelSouth"])
        TYPE <- "N"
    if (!OUT$Species[spp, "ModelNorth"] && OUT$Species[spp, "ModelSouth"])
        TYPE <- "S"

    NAM <- OUT$Species[spp, "CommonName"]
    if (is.na(NAM))
        NAM <- OUT$Species[spp, "ScientificName"]
    NAM <- as.character(NAM)

    if (TYPE != "S") {
        ## hab-north
        png(paste0(ROOT, "/figs/", gr, "/", spp, "-coef-north.png"),
            height=800, width=2200, res=150)
        layout(matrix(c(1,1,2), nrow=1))
        plot_abundance(spp, "veg_coef")
        par(mar=c(12,5,4,3))
        plot_abundance(spp, "veg_lin", main="")
        dev.off()

        png(paste0(ROOT, "/figs/", gr, "/", spp, "-sector-north.png"),
            height=2*500, width=2*1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(resn[spp,], "unit")
        plot_sector(resn[spp,], "regional", main="")
        plot_sector(resn[spp,], "underhf", main="")
        par(op)
        dev.off()
    }
    if (TYPE != "N") {
        ## hab-south
        png(paste0(ROOT, "/figs/", gr, "/", spp, "-coef-south.png"),
            height=800, width=2000, res=150)
        layout(matrix(c(1,2,3), nrow=1))
        p1 <- plot_abundance(spp, "soil_coef", paspen=0, plot=FALSE)
        p2 <- plot_abundance(spp, "soil_coef", paspen=1, plot=FALSE)
        plot_abundance(spp, "soil_coef", paspen=0, ylim=c(0, max(p1, p2)), main=paste0(NAM, " - non treed"))
        plot_abundance(spp, "soil_coef", paspen=1, ylim=c(0, max(p1, p2)), main="treed")
        par(mar=c(11,5,4,3))
        plot_abundance(spp, "soil_lin", main="")
        dev.off()

        png(paste0(ROOT, "/figs/", gr, "/", spp, "-sector-south.png"),
            height=2*500, width=2*1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(ress[spp,], "unit")
        plot_sector(ress[spp,], "regional", main="")
        plot_sector(ress[spp,], "underhf", main="")
        par(op)
        dev.off()
    }

if (TRUE) {
    y <- load_species_data(spp, boot=TRUE)
    ry <- rasterize_results(y)

    Curr <- y$SA.Curr[match(rownames(kgrid), rownames(y$SA.Curr)),]
    Ref <- y$SA.Ref[match(rownames(kgrid), rownames(y$SA.Ref)),]

    Dcr <- rowSums(Curr)
    q <- quantile(Dcr, 0.99)
    Dcr[Dcr > q] <- q
    Drf <- rowSums(Ref)
    q <- quantile(Drf, 0.99)
    Drf[Drf > q] <- q
    MAX <- max(Dcr, Drf)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)

    Rcr <- make_raster(cr, kgrid, rt)
    Rrf <- make_raster(rf, kgrid, rt)
    Rdf <- make_raster(df-100, kgrid, rt)
    if (TYPE == "S")
        Msk <- Rmasks
    if (TYPE == "N")
        Msk <- Rmaskn
    if (TYPE != "C") {
        Rcr <- mask(Rcr, Msk)
        Rrf <- mask(Rrf, Msk)
        Rdf <- mask(Rdf, Msk)
    }
    if (gr != "birds") {
        Rcr <- mask(Rcr, Rmaskm)
        Rrf <- mask(Rrf, Rmaskm)
        Rdf <- mask(Rdf, Rmaskm)
    }

    writeRaster(Rcr, paste0(ROOT, "/normalized-maps/", gr, "/", spp, "-cr.tif"), overwrite=TRUE)
    writeRaster(Rrf, paste0(ROOT, "/normalized-maps/", gr, "/", spp, "-rf.tif"), overwrite=TRUE)
    writeRaster(Rdf, paste0(ROOT, "/normalized-maps/", gr, "/", spp, "-df.tif"), overwrite=TRUE)
}

    ## add here mask for Rockies if needed
if (TRUE) {
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-map.png"),
        height=1500*2, width=1000*2, res=300)
    op <- par(mfrow=c(2,2), mar=c(2,1,2,3))

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Reference", legend=FALSE)
    plot(Rrf, add=TRUE, col=col1[1:max(rf)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Current", legend=FALSE)
    plot(Rcr, add=TRUE, col=col1[1:max(cr)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Difference", legend=FALSE)
    plot(Rdf, add=TRUE, col=col3[min(df):max(df)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Std. Error (Current)", legend=FALSE)
    plot(ry[["SE"]]/MAX, add=TRUE, col=rev(col2))
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    par(op)
    dev.off()

    gc()
}

}}

## use avail figures

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

## detection maps for non birds

m <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

rnr <- make_raster(as.integer(kgrid$NRNAME), kgrid, rt)
cnr <- c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae')
cnr <- cnr[c(5,6,1,2,4,3)]


ex <- new.env()
gr <- "mites"
load("s:/Result from Ermias_2018/mites/Species detection Mites 2018.RData", envir=ex)

ex <- new.env()
gr <- "lichens"
load("s:/Result from Ermias_2018/lichens/Species detection Lichens 2018.RData", envir=ex)

ex <- new.env()
gr <- "mosses"
load("s:/Result from Ermias_2018/mosses/Species detection Moss 2018.RData", envir=ex)

ex <- new.env()
gr <- "vplants"
load("s:/Result from Ermias_2018/vplants/Species detection Vascular plants 2018.RData", envir=ex)



site <- ex$dd$Site
og <- ex$dd$OnOffGrid == "OG"
ogs <- sapply(strsplit(site[og], "-"), "[[", 3)
site[og] <- ogs
site <- gsub("B", "", site)
siten <- as.integer(site)

xy <- data.frame(x=m$PUBLIC_LONGITUDE, y=m$PUBLIC_LATTITUDE)[match(site, m$SITE_ID),]
coordinates(xy) <- ~x+y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))

yy <- ex$dd[,6:ncol(ex$dd)]
compare_sets(colnames(yy), rownames(tab[tab$Group == gr,]))

SPP <- rownames(tab[tab$Group == gr,])

for (spp in SPP) {
    spnam <- if (is.na(tab[spp, "CommonName"])) {
        as.character(tab[spp, "ScientificName"])
    } else {
        paste0(as.character(tab[spp, "CommonName"]), " (", as.character(tab[spp, "ScientificName"]), ")")
    }
    cat(gr, spp, "\n");flush.console()
    xy0 <- xy[yy[,spp] == 0 & !duplicated(site),]
    xy1 <- xy[yy[,spp] > 0,]
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-det.png"),
        height=1500*1.5, width=1000*1.5, res=300)
    op <- par(mar=c(1,1,1,1))
    plot(rnr,col=cnr, axes=FALSE, box=FALSE, main=spnam, legend=FALSE)
    plot(xy0, add=TRUE, pch=19, col="#aaaaaa88", legend=FALSE, cex=0.8)
    plot(xy1, add=TRUE, pch=19, col="red4", legend=FALSE, cex=0.8)
    par(op)
    dev.off()

}

## detection maps for birds

## detections
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)

ddd <- ee$dd
## subset based on analysis data -- same filtering applied
ddd <- ddd[unique(c(rownames(en$DAT), rownames(es$DAT))),]
ddd <- nonDuplicated(ddd, ddd$SS, TRUE)
yyy <- groupSums(ee$yy, 1, ee$dd$SS)[rownames(ddd),]
yyy[yyy > 0] <- 1
ss <- !is.na(ddd$X) & !is.na(ddd$NRNAME)
ddd <- ddd[ss,]
yyy <- yyy[ss,]

xy <- SpatialPoints(as.matrix(ddd[,c("X","Y")]))
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))
rt10 <- aggregate(rt, fact=10)
sam0 <- rasterize(xy, rt10, field=1, fun='sum')
values(sam0)[!is.na(values(sam0))] <- 1

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$UseavailNorth | Lookup$UseavailSouth,])
rownames(tax) <- tax$Code

SPP <- rownames(OUT$Species[OUT$Species$Group == gr,])
tax <- tax[tax$SpeciesID %in% SPP,]

gr <- "birds"
for (spp in rownames(tax)) {
    cat(gr, spp, "\n");flush.console()
    xy1 <- SpatialPoints(as.matrix(ddd[yyy[,spp] > 0,c("X","Y")]))
    proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    xy1 <- spTransform(xy1, proj4string(rt))
    sam1 <- rasterize(xy1, rt10, field=1, fun='last')
    png(paste0("d:/abmi/AB_data_v2018/www", "/figs/", gr, "/", as.character(tax[spp, "SpeciesID"]), "-det.png"),
        height=1500*1.5, width=1000*1.5, res=300)
    op <- par(mar=c(1,1,1,1))
    plot(rnr,col=cnr, axes=FALSE, box=FALSE, main=as.character(tax[spp, "CommonName"]), legend=FALSE)
    plot(sam0,add=TRUE, col="#aaaaaa88", legend=FALSE)
    plot(sam1,add=TRUE, col="red4", legend=FALSE)
    par(op)
    dev.off()
}




## mammal camera dump

fl <- list.files("s:/Camera mammals Mar 2019/Figures North/Best model")
x <- gsub("Veg+HF figure best model ", "", fl, fixed=TRUE)
#x[endsWith(x, "Winter.png")] <- gsub("Winter.png", "-winter.png", x[endsWith(x, "Winter.png")])
#x[endsWith(x, "Summer.png")] <- gsub("Summer.png", "-summer.png", x[endsWith(x, "Summer.png")])

i1 <- endsWith(x, "Winter.png")
i2 <- endsWith(x, "Summer.png")

## combined (fine coef, summer+winter)
file.copy(paste0("s:/Camera mammals Mar 2019/Figures North/Best model/", fl[!i1 & !i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-north-combo/", x[!i1 & !i2]))
## honest models, winter
file.copy(paste0("s:/Camera mammals Mar 2019/Figures North/Best model/", fl[i1]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-north-winter/",
    gsub("Winter.png", ".png", x[i1])))
## honest models, summer
file.copy(paste0("s:/Camera mammals Mar 2019/Figures North/Best model/", fl[i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-north-summer/",
    gsub("Summer.png", ".png", x[i2])))

fl1 <- list.files("s:/Camera mammals Mar 2019/Figures South/Best model/Treed")
fl2 <- list.files("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed")
x1 <- gsub("Soil+HF figure best model ", "", fl1, fixed=TRUE)
x2 <- gsub("Soil+HF figure best model ", "", fl2, fixed=TRUE)
all(x1 == x2)

i1 <- endsWith(x1, "Winter.png")
i2 <- endsWith(x1, "Summer.png")

## combined (fine coef, summer+winter)
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Treed/", fl1[!i1 & !i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-combo-treed/", x1[!i1 & !i2]))
## honest models, winter
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Treed/", fl1[i1]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-winter-treed/",
    gsub("Winter.png", ".png", x1[i1])))
## honest models, summer
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Treed/", fl1[i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-summer-treed/",
    gsub("Summer.png", ".png", x1[i2])))

file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed/", fl2[!i1 & !i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-combo-nontreed/", x2[!i1 & !i2]))
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed/", fl2[i1]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-winter-nontreed/",
    gsub("Winter.png", ".png", x2[i1])))
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed/", fl2[i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-summer-nontreed/",
    gsub("Summer.png", ".png", x2[i2])))

## jpg to png

library(magick)
fl <- list.files("s:/Camera mammals Mar 2019/Maps/Space climate North/")
x <- gsub(".jpg", ".png", fl, fixed=TRUE)
for (i in 1:length(fl)) {
    img <- image_read(paste0("s:/Camera mammals Mar 2019/Maps/Space climate North/", fl[i]))
    image_write(img, paste0("s:/Camera mammals Mar 2019/Maps/Space climate North/", x[i]), format="png")
}
fl <- list.files("s:/Camera mammals Mar 2019/Maps/Space climate South/")
x <- gsub(".jpg", ".png", fl, fixed=TRUE)
for (i in 1:length(fl)) {
    img <- image_read(paste0("s:/Camera mammals Mar 2019/Maps/Space climate South/", fl[i]))
    image_write(img, paste0("s:/Camera mammals Mar 2019/Maps/Space climate South/", x[i]), format="png")
}


## lookup table

r <- "d:/abmi/reports/2018/images/mammals-camera/"
fl <- list.files("d:/abmi/reports/2018/images/mammals-camera/map-det")

d <- list.dirs(r, full.names=FALSE)[-1]

ls <- lapply(d, function(z) list.files(paste0(r, z)))
names(ls) <- d

spp <- sort(unique(unlist(ls)))

mefa4::compare_sets(fl, spp)
setdiff(fl, spp)
setdiff(spp, fl)
intersect(spp, fl)


SPP <- c(
    "Badger" = "Badger",
    "Beaver" = "Beaver",
    "Bighornsheep" = "Bighorn Sheep",
    "Bison" = "Bison",
    "BlackBear" = "Black Bear",
    "Bobcat" = "Bobcat",
    "CanadaLynx" = "Canada Lynx",
    "ColumbianGroundSquirrel" = "Columbian Ground Squirrel",
    "Cougar" = "Cougar",
    "Coyote" = "Coyote",
    "Deer" = "Deer",
    "Elk" = "Elk",
    "Fisher" = "Fisher",
    "Foxes" = "Foxes",
    "GoldenMantledGroundSquirrel" = "Golden Mantled Ground Squirrel",
    "GrayWolf" = "Gray Wolf",
    "Grizzlybear" = "Grizzly Bear",
    "Groundhog" = "Groundhog",
    "HoaryMarmot" = "Hoary Marmot",
    "LeastChipmunk" = "Least Chipmunk",
    "Marten" = "Marten",
    "Mink" = "Mink",
    "Moose" = "Moose",
    "Mountaingoat" = "Mountain Goat",
    "Muledeer" = "Muledeer",
    "Muskrat" = "Muskrat",
    "NorthernFlyingSquirrel" = "Northern Flying Squirrel",
    "Porcupine" = "Porcupine",
    "Pronghorn" = "Pronghorn",
    "Raccoon" = "Raccoon",
    "Redfox" = "Red Fox",
    "RedSquirrel" = "Red Squirrel",
    "RichardsonsGroundSquirrel" = "Richardson's Ground Squirrel",
    "RiverOtter" = "River Otter",
    "SnowshoeHare" = "Snowshoe Hare",
    "StripedSkunk" = "Striped Skunk",
    "VolesMiceandAllies" = "Voles, Mice and Allies",
    "WeaselsandErmine" = "Weasels and Ermine",
    "WhitetailedDeer" = "White-tailed Deer",
    "WhitetailedJackRabbit" = "Whitetailed Jack Rabbit",
    "Wolverine" = "Wolverine",
    "WolvesCoyotesandAllies" = "Wolves, Coyotes and Allies",
    "WoodlandCaribou" = "Woodland Caribou")

tab <- data.frame(SpeciesID=names(SPP), CommonName=SPP)
for (i in d) {
    tab[[i]] <- names(SPP) %in% gsub(".png", "", ls[[i]])
}


library(jsonlite)
toJSON(tab, rownames=FALSE,pretty=TRUE)




library(cure4insect)
library(mefa4)
set_options(path = "s:/reports")
load_common_data()

tab <- get_species_table()
tab$DisplayName <- paste0(tab$CommonName, " (", tab$ScientificName, ")")
toJSON(tab[c("AlderFlycatcher", "Ovenbird"),], rownames=FALSE,pretty=TRUE)


gr <- "birds"
SPP <- rownames(resn[resn$Taxon == gr,])


## API ----------------------------------
## all species: display name, group, speciesID

library(jsonlite)

load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")
tab <- OUT$Species

dot <- endsWith(rownames(tab), ".")
dotBad <- rownames(tab)[dot]
dotOK <- substr(dotBad, 1, nchar(dotBad)-1)
data.frame(Bad=dotBad,Good=dotOK, tab$Group[dot])

for (i in seq_along(dotBad)) {
    levels(tab$SpeciesID)[levels(tab$SpeciesID) == dotBad[i]] <- dotOK[i]
    rownames(tab)[rownames(tab) == dotBad[i]] <- dotOK[i]
}
any(endsWith(rownames(tab), "."))

tab$DisplayName <- ifelse(is.na(tab$CommonName), as.character(tab$ScientificName),
    paste0(tab$CommonName, " (", tab$ScientificName, ")"))

grs <- c("birds", "vplants", "mites", "lichens", "mosses", "mammals-camera")

tmp <- list()
for (g in grs[1:5]) {
    p <- file.path("d:/abmi/reports/2018", "images", g)
    f <- list.files(p)
    s <- unique(sapply(strsplit(f, "-"), "[[", 1))
    tmp[[g]] <- tab[s,c("SpeciesID","DisplayName", "Group")]
    tmp[[g]]$sppprevious <- c(rownames(tmp[[g]])[nrow(tmp[[g]])], rownames(tmp[[g]])[-nrow(tmp[[g]])])
    tmp[[g]]$sppnext <- c(rownames(tmp[[g]])[-1], rownames(tmp[[g]])[1])
}


SPP <- c(
    "Badger" = "Badger",
    "Beaver" = "Beaver",
    "Bighornsheep" = "Bighorn Sheep",
    "Bison" = "Bison",
    "BlackBear" = "Black Bear",
    "Bobcat" = "Bobcat",
    "CanadaLynx" = "Canada Lynx",
    "ColumbianGroundSquirrel" = "Columbian Ground Squirrel",
    "Cougar" = "Cougar",
    "Coyote" = "Coyote",
    "Deer" = "Deer",
    "Elk" = "Elk",
    "Fisher" = "Fisher",
    "Foxes" = "Foxes",
    "GoldenMantledGroundSquirrel" = "Golden Mantled Ground Squirrel",
    "GrayWolf" = "Gray Wolf",
    "Grizzlybear" = "Grizzly Bear",
    "Groundhog" = "Groundhog",
    "HoaryMarmot" = "Hoary Marmot",
    "LeastChipmunk" = "Least Chipmunk",
    "Marten" = "Marten",
    "Mink" = "Mink",
    "Moose" = "Moose",
    "Mountaingoat" = "Mountain Goat",
    "Muledeer" = "Muledeer",
    "Muskrat" = "Muskrat",
    "NorthernFlyingSquirrel" = "Northern Flying Squirrel",
    "Porcupine" = "Porcupine",
    "Pronghorn" = "Pronghorn",
    "Raccoon" = "Raccoon",
    "Redfox" = "Red Fox",
    "RedSquirrel" = "Red Squirrel",
    "RichardsonsGroundSquirrel" = "Richardson's Ground Squirrel",
    "RiverOtter" = "River Otter",
    "SnowshoeHare" = "Snowshoe Hare",
    "StripedSkunk" = "Striped Skunk",
    "VolesMiceandAllies" = "Voles, Mice and Allies",
    "WeaselsandErmine" = "Weasels and Ermine",
    "WhitetailedDeer" = "White-tailed Deer",
    "WhitetailedJackRabbit" = "Whitetailed Jack Rabbit",
    "Wolverine" = "Wolverine",
    "WolvesCoyotesandAllies" = "Wolves, Coyotes and Allies",
    "WoodlandCaribou" = "Woodland Caribou")

g <- grs[6]
tmp[[g]] <- data.frame(SpeciesID=names(SPP), DisplayName=SPP, Group=g)
tmp[[g]]$sppprevious <- c(rownames(tmp[[g]])[nrow(tmp[[g]])], rownames(tmp[[g]])[-nrow(tmp[[g]])])
tmp[[g]]$sppnext <- c(rownames(tmp[[g]])[-1], rownames(tmp[[g]])[1])


ALL <- do.call(rbind, tmp)

writeLines(toJSON(ALL, rownames=FALSE), file.path("d:/abmi/reports/2018", "api", "index.json"))

for (g in grs) {
    writeLines(toJSON(tmp[[g]], rownames=FALSE), file.path("d:/abmi/reports/2018", "api", g, "index.json"))
}

figs <- c("det", "useavail-north", "useavail-south",
    "coef-north", "coef-south", "map", "sector-north", "sector-south")

for (g in grs[1:5]) {
    s <- as.character(tmp[[g]]$SpeciesID)
    f <- list.files(file.path("d:/abmi/reports/2018", "images", g))
    m <- as.data.frame(matrix(FALSE, length(s), length(figs)))
    dimnames(m) <- list(s, figs)

    for (i in s) {
#        ff <- f[startsWith(f, i)]
#        ff <- gsub(".png", "", gsub(paste0(i, "-"), "", ff))
        d <- OUT$Species[i,]
        ff <- c("det",
            if (d$ModelNorth) c("sector-north", "coef-north") else NULL,
            if (d$ModelSouth) c("sector-south", "coef-south") else NULL,
            if (d$ModelNorth || d$ModelSouth) c("map") else NULL,
            if (d$UseavailNorth && !d$ModelNorth) c("useavail-north") else NULL,
            if (d$UseavailSouth && !d$ModelSouth) c("useavail-south") else NULL
        )
        m[i,] <- figs %in% ff
        v <- cbind(tab[i,], tmp[[g]][i,c("sppprevious", "sppnext")], m[i,])
        #toJSON(as.list(v), rownames=FALSE,pretty=TRUE,auto_unbox=TRUE)
        if (!dir.exists(file.path("d:/abmi/reports/2018", "api", g, i)))
            dir.create(file.path("d:/abmi/reports/2018", "api", g, i))
        writeLines(toJSON(as.list(v), rownames=FALSE, auto_unbox=TRUE,pretty=TRUE),
            file.path("d:/abmi/reports/2018", "api", g, i, "index.json"))
    }
}

figs <- list.dirs(file.path("d:/abmi/reports/2018", "images", g),full.names=FALSE)[-1]
g <- grs[6]
s <- as.character(tmp[[g]]$SpeciesID)
f <- list.files(file.path("d:/abmi/reports/2018", "images", g))
m <- as.data.frame(matrix(FALSE, length(s), length(figs)))
dimnames(m) <- list(s, figs)
for (j in figs) {
    fff <- gsub(".png", "", list.files(file.path("d:/abmi/reports/2018", "images", g, j)))
    m[fff,j] <- TRUE
}

for (i in s) {
    ff <- f[startsWith(f, i)]
    ff <- gsub(".png", "", gsub(paste0(i, "-"), "", ff))
    v <- cbind(tmp[[g]][i,], m[i,])
    #toJSON(as.list(v), rownames=FALSE,pretty=TRUE,auto_unbox=TRUE)
    if (!dir.exists(file.path("d:/abmi/reports/2018", "api", g, i)))
        dir.create(file.path("d:/abmi/reports/2018", "api", g, i))
    writeLines(toJSON(as.list(v), rownames=FALSE, auto_unbox=TRUE,pretty=TRUE),
        file.path("d:/abmi/reports/2018", "api", g, i, "index.json"))
}

## clean up dotted names
if (FALSE) {
fl <- list.files("d:/abmi/reports/2018", recursive = TRUE)
xx <- NULL
for (i in seq_along(dotBad)) {
    xx <- c(xx, grep(dotBad[i], fl))
}
xx <- sort(unique(xx))
str(xx)
fl[xx]
for (j in xx) {
    In <- Out <- fl[j]
    for (i in seq_along(dotBad)) {
        Out <- gsub(dotBad[i], dotOK[i], Out)
    }
    file.rename(file.path("d:/abmi/reports/2018", In),
        file.path("d:/abmi/reports/2018", Out))
}
}
