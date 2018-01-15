## processing common data (kgrid, areas, species lookups)

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(raster)

## HF 2012
load("e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv0$Sector2 <- factor(ifelse(is.na(tv0$Sector), "NATIVE", as.character(tv0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
VHF <- dd1km_pred[[1]]
tv0 <- tv0[colnames(VHF),]
dd1km_pred$sample_year
KA_2012 <- groupSums(VHF/10^6, 2, tv0$Sector2) # in km^2

## HF 2014
load("e:/peter/AB_data_v2017/data/analysis/veg-hf_1kmgrid_v6-fixage0.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
VHF <- dd1km_pred[[1]]
VHF <- VHF[,colnames(VHF) != "CutBlocks"]
cn <- colnames(VHF)
cn[cn %in% c("UrbanIndustrial", "UrbanResidence")] <- "Urban"
cn[cn %in% c("SeismicLineNarrow", "SeismicLineWide")] <- "SeismicLine"
cn[cn %in% c("CultivationAbandoned", "CultivationRoughPasture",
    "CultivationTamePasture", "CultivationCrop")] <- "CultivationCropPastureBareground"
VHF <- groupSums(VHF, 2, cn)
tv0 <- tv0[colnames(VHF),]
dd1km_pred$sample_year
KA_2014 <- groupSums(VHF/10^6, 2, tv0$ETA_UseInAnalysis_Sector) # in km^2

load("e:/peter/AB_data_v2017/data/analysis/kgrid_table_km.Rdata")
stopifnot(all(rownames(kgrid) == rownames(KA_2012)))
stopifnot(all(rownames(kgrid) == rownames(KA_2014)))

KT <- kgrid[,c("Row", "Col", "Row10_Col10")]
XY <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(XY) <- ~ POINT_X + POINT_Y
proj4string(XY) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

SP <- read.csv("w:/reports/2017/data/species-info.csv")
SP <- droplevels(SP[SP$map.pred,])
SP$infoOK <- TRUE
#SP$infoOKboot <- TRUE
for (i in 1:nrow(SP)) {
fn1 <- file.path("w:/reports/2017/results", as.character(SP[i,"taxon"]), "sector",
    paste0(as.character(SP[i,"SpeciesID"]), ".RData"))
#fn2 <- file.path("w:/reports/2017/results", as.character(SP[i,"taxon"]), "sector",
#    paste0(as.character(SP[i,"SpeciesID"]), ".RData"))
if (!file.exists(fn1))
    SP$infoOK[i] <- FALSE
#if (!file.exists(fn2))
#    SP$infoOKboot[i] <- FALSE
}
#table(SP$infoOK, SP$infoOKboot)
table(SP$infoOK, SP$taxon)

SP <- SP[SP$infoOK,]
SP$infoOK <- NULL
rownames(SP) <- SP$SpeciesID
save(KA_2012, KA_2014, KT, XY, SP,
    file="w:/reports/2017/data/kgrid_areas_by_sector.RData")

## no need to normalize files
#"Curr.Boot" "Ref.Boot" -- > CB/RB
#w:/reports/2017/Files from Ermias/Mites/Combine regions/Boot data/
#"SA.Curr" "SA.Ref" -- > CS/RS
#w:/reports/2017/Files from Ermias/Mites/Combine regions/Sector effects/Sector abundance summary/
#Native Misc Agriculture RuralUrban Energy Transportation Forestry


Sectors <- c("Native", "Misc", "Agriculture", "RuralUrban", "Energy", "Transportation", "Forestry")





load("e:/peter/AB_data_v2017/data/analysis/site/veg-hf_siteCenter_v6-fixage0.Rdata")
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
vhf <- as.matrix(dd_1ha$veg_current)
vhf <- vhf[,colnames(vhf) != "CutBlocks"]
a <- vhf / rowSums(vhf)
a <- a[,rownames(tv)]
a <- groupSums(a, 2, tv$ETA_UseInAnalysis_Sector)
summary(a)


## processing vplant data for sample based non-native sector effects

library(mefa4)

load("e:/peter/AB_data_v2017/data/analysis/site/veg-hf_siteCenter_v6-fixage0.Rdata")
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
vhf <- as.matrix(dd_1ha$veg_current)
vhf <- vhf[,colnames(vhf) != "CutBlocks"]
a <- vhf / rowSums(vhf)
a <- a[,rownames(tv)]
a <- groupSums(a, 2, tv$ETA_UseInAnalysis_Sector)
summary(a)

load("e:/peter/AB_data_v2017/data/misc/VPlant.RData")

str(dd)
Y <- as.matrix(dd[,6:ncol(dd)])
Y[Y>0] <- 1

X <- dd[,1:6]
compare_sets(rownames(X), rownames(a))
Sites <- intersect(rownames(X), rownames(a))

X <- X[Sites,]
Y <- Y[Sites,]
A <- a[Sites,]

sp <- read.csv("e:/peter/AB_data_v2017/data/misc/VPlant Species analysis name_April 2017.csv")
compare_sets(sp$Analysis_Name, colnames(Y))
Z <- sp[match(colnames(Y), sp$Analysis_Name),]
rownames(Z) <- Z$Analysis_Name
Z <- droplevels(Z[Z$RANK_NAME == "Species",])
Z$NN <- Z$Origin_Analysis != "Native"
Y <- Y[,rownames(Z)]

#save(X, Y, A, Z, file="e:/peter/AB_data_v2017/data/reporting/vpalnts-2017-11-16.Rdata")

hf <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
knitr::kable(hf[,c("HF_GROUP", "Sector")])

library(mefa4)
load("e:/peter/AB_data_v2017/data/reporting/Summary of VPlant abundance by sector.RData")
load("e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv0$Sector2 <- factor(ifelse(is.na(tv0$Sector), "NATIVE", as.character(tv0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
VHF <- dd1km_pred[[1]]
tv0 <- tv0[colnames(VHF),]

CS <- colSums(groupSums(VHF/10^6, 2, tv0$Sector2))
CS <- CS / sum(CS)

Curr <- VP_SectorAbundance.Curr["Achillea.millefolium",]
Ref <- VP_SectorAbundance.Ref["Achillea.millefolium",]
Area <- 100*CS

## Sector effect plot from Dave
plot_sector_1 <- function(Curr, Ref, Area, main="") {
    sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
    sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")
    c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")

    total.effect <- (100 * (Curr - Ref) / sum(Ref))[sectors]
    unit.effect <- 100 * total.effect / Area[sectors]
    ymax <- ifelse(max(abs(unit.effect))<20,20,
        ifelse(max(abs(unit.effect))<50,50,round(max(abs(unit.effect))+50,-2)))
    ymin <- ifelse(ymax>50,min(-100,round(min(unit.effect)-50,-2)),-ymax)
    ymax <- max(ymax,max(unit.effect)+0.08*(max(unit.effect)-min(unit.effect,0)))
    ymin <- min(ymin,min(unit.effect)-0.08*(max(unit.effect,0)-min(unit.effect)))
    q <- barplot(unit.effect,
        width=Area[sectors],
        space=0,col=c1,border=c1,ylim=c(ymin,ymax),
        ylab="Unit effect (%)",xlab="Area (% of region)",
        xaxt="n",cex.lab=1.3,cex.axis=1.2,tcl=0.3,
        xlim=c(0,round(sum(Area[sectors])+1,0)),
        bty="n",col.axis="grey40",col.lab="grey40",las=2)

    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray88",border="gray88")
    x.at<-pretty(c(0,sum(Area[sectors])))
    axis(side=1,tck=1,at=x.at,lab=rep("",length(x.at)),col="grey95")
    y.at<-pretty(c(ymin,ymax),n=6)
    axis(side=2,tck=1,at=y.at,lab=rep("",length(y.at)),col="grey95")
    q <- barplot(unit.effect,
        width=Area[sectors],
        space=0,col=c1,border=c1,ylim=c(ymin,ymax),
        ylab="Unit effect (%)",xlab="Area (% of region)",
        xaxt="n",cex.lab=1.3,cex.axis=1.2,tcl=0.3,
        xlim=c(0,round(sum(Area[sectors])+1,0)),
        bty="n",col.axis="grey40",col.lab="grey40",las=2,add=TRUE)
    box(bty="l",col="grey40")
    #mtext(side=1,line=2,at=x.at,x.at,col="grey40",cex=1.2)
    axis(side=1,at=x.at,tcl=0.3,lab=rep("",length(x.at)),col="grey40",
        col.axis="grey40",cex.axis=1.2,las=1)
    abline(h=0,lwd=2,col="grey40")
    mtext(side=1,at=q+c(0,0,-1,0,+1),sector.names,col=c1,cex=1.3,
        adj=0.5,line=c(0.1,0.1,1.1,0.1,1.1))
    y <- unit.effect+0.025*(ymax-ymin)*sign(unit.effect)
    if (abs(y[3]-y[4])<0.05*(ymax-ymin))
        y[3:4]<-mean(y[3:4])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[3:4])]
    if (abs(y[4]-y[5])<0.05*(ymax-ymin))
        y[4:5]<-mean(y[4:5])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[4:5])]
    text(q,y,paste(ifelse(total.effect>0,"+",""),
        sprintf("%.1f",total.effect),"%",sep=""),col="darkblue",cex=1.4)
    mtext(side=3,line=1,at=0,adj=0, main, cex=1.4,col="grey40")
    invisible(rbind(total=total.effect, unit=unit.effect, area=Area[sectors]))
}
plot_sector_2 <- function(Curr, Ref, regional=TRUE, main="") {
    sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
    sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")
    c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")
    total.effect <- if (regional)
        100 * (Curr - Ref)/sum(Ref) else 100 * (Curr - Ref)/Ref
    total.effect <- total.effect[sectors]
    off <- 0.25
    a <- 1-0.5-off
    b <- 5+0.5+off
    ymax <- ifelse(max(abs(total.effect))<20,20,
        ifelse(max(abs(total.effect))<50,50,round(max(abs(total.effect))+50,-2)))
    ymin <- ifelse(ymax>50,min(-100,round(min(total.effect)-50,-2)),-ymax)
    ymax <- max(ymax,max(total.effect)+0.08*(max(total.effect)-min(total.effect,0)))
    ymin <- min(ymin,min(total.effect)-0.08*(max(total.effect,0)-min(total.effect)))
    yax <- pretty(c(ymin,ymax))
    op <- par(las=1, xpd = TRUE)
    on.exit(par(op))
    plot(0, type="n", xaxs="i", yaxs = "i", ylim=c(ymin,ymax), xlim=c(a, b),
        axes=FALSE, ann=FALSE)
    polygon(c(a,a,b,b), c(ymin, ymax, ymax, ymin), col="grey88", border="grey88")
    segments(x0=rep(a, length(yax)), x1=rep(b,length(yax)),y0=yax, col="white")
    axis(2, yax, paste0(ifelse(yax>0, "+", ""), yax), tick=FALSE)
    rug(yax, side=2, ticksize=0.01, col="grey40", quiet=TRUE)
    lines(c(a,a), c(ymin, ymax), col="grey40", lwd=1)
    for (i in 1:5) {
        h <- total.effect[i]
        polygon(c(i-0.5, i-0.5, i+0.5, i+0.5), c(0,h,h,0), col=c1[i], border=NA)
    }
    lines(c(a,b), c(0, 0), col="grey40", lwd=2)
    title(ylab=if (regional) "Regional sector effects (%)" else "Local sector effects (%)",
        cex=1.3, col="grey40")
    mtext(side=1,at=1:5,sector.names,col=c1,cex=1.3,adj=0.5,line=0.5)

    y <- total.effect+0.025*(ymax-ymin)*sign(total.effect)
    if (abs(y[3]-y[4])<0.05*(ymax-ymin))
        y[3:4]<-mean(y[3:4])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[3:4])]
    text(1:5,y,paste(sprintf("%.1f",total.effect),"%",sep=""),col="darkblue",cex=1.2)
    mtext(side=3,line=1,at=0,adj=0, main, cex=1.4,col="grey40")
    invisible(total.effect)
}

plot_sector_1(Curr, Ref, Area, main="Common Yarrow")

op <- par(mfrow=c(1,2))
plot_sector_2(Curr, Ref, regional=TRUE, main="Common Yarrow")
plot_sector_2(Curr, Ref, regional=FALSE, main="Common Yarrow")
par(op)

png(paste0("combinedSE/", species_name, ".png")), width=600, height=600)
plot_sector_1(Curr, Ref, Area, main=species_name)
dev.off()

png(paste0("regionalSE/", species_name, ".png")), width=600, height=600)
plot_sector_2(Curr, Ref, regional=TRUE, main=species_name)
dev.off()

png(paste0("localSE/", species_name, ".png")), width=600, height=600)
plot_sector_2(Curr, Ref, regional=FALSE, main=species_name)
dev.off()

sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
SEreg <- 100 * (VP_SectorAbundance.Curr - VP_SectorAbundance.Ref) / rowSums(VP_SectorAbundance.Ref)
#A <- 100*CS
#SEunit <- t(100 * t(SEreg) / A)
SEloc <- 100 * (VP_SectorAbundance.Curr - VP_SectorAbundance.Ref) / VP_SectorAbundance.Ref
SEreg <- SEreg[,sectors]
#SEunit <- SEunit[,sectors]
SEloc <- SEloc[,sectors]

plot_sector_3 <- function(x, ylab="Sector effects (%)", type="kde",
breaks = "Sturges", ...) {
    type <- match.arg(type, c("kde", "fft", "hist"))
    if (!is.list(x))
        x <- as.data.frame(x)
    sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
    sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")
    c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")
    ymin <- -100
    ymax <- 100
    off <- 0.25
    a <- 1-0.5-off
    b <- 5+0.5+off
    v <- 0.1
    yax <- pretty(c(ymin,ymax))
    op <- par(las=1)
    on.exit(par(op))
    plot(0, type="n", xaxs="i", yaxs = "i", ylim=c(ymin,ymax), xlim=c(a, b),
        axes=FALSE, ann=FALSE)
    polygon(c(a,a,b,b), c(ymin, ymax, ymax, ymin), col="grey88", border="grey88")
    segments(x0=rep(a, length(yax)), x1=rep(b,length(yax)),y0=yax, col="white")
    axis(2, yax, paste0(ifelse(yax>0, "+", ""), yax), tick=FALSE)
    rug(yax, side=2, ticksize=0.01, col="grey40", quiet=TRUE)
    lines(c(a,a), c(ymin, ymax), col="grey40", lwd=1)
    lines(c(a,b), c(0, 0), col="grey40", lwd=2)
    out <- list()
    for (i in 1:5) {
        xx <- sort(x[[i]])
        k <- xx <= ymax
        out[[i]] <- sum(!k)
        st <- boxplot.stats(xx)
        s <- st$stats
        k[which(!k)[1]] <- TRUE
        if (type == "kde")
            d <- KernSmooth::bkde(xx[k]) # uses Normal kernel
        if (type == "fft")
            d <- density(xx[k]) # uses FFT
        if (type == "hist") {
            h <- hist(xx[k], plot=FALSE, breaks=breaks)
            xv <- rep(h$breaks, each=2)
            yv <- c(0, rep(h$density, each=2), 0)
        } else {
            xv <- d$x
            yv <- d$y
            j <- xv >= min(xx) & xv <= max(xx)
            xv <- xv[j]
            yv <- yv[j]
        }
        yv <- 0.4 * yv / max(yv)
        polygon(c(-yv, rev(yv))+i, c(xv, rev(xv)), col=c1[i], border=c1[i])
        polygon(c(-v,-v,v,v)+i, s[c(2,4,4,2)], col="#40404080", border=NA)
        lines(c(-v,v)+i, s[c(3,3)], lwd=2, col="grey30")
    }
    title(ylab=ylab, cex=1.3, col="grey40")
    mtext(side=1,at=1:5,sector.names,col=c1,cex=1.3,adj=0.5,line=0.5)
    op <- par(xpd = TRUE)
    on.exit(par(op), add=TRUE)
    out <- unlist(out)
    points(1:5, rep(105, 5), pch=19,
        cex=ifelse(out==0, 0, 0.5+2*out/max(out)), col=c1)
    invisible(x)
}
op <- par(mfrow=c(3,2))
plot_sector_3(SEreg, ylab="Regional (total) sector effects (%)", type="kde")
title(main="Kernel density")
plot_sector_3(SEloc, ylab="Local (within footprint) sector effects (%)", type="kde")

plot_sector_3(SEreg, ylab="Regional (total) sector effects (%)", type="fft")
title(main="Fast Fourier transform based density")
plot_sector_3(SEloc, ylab="Local (within footprint) sector effects (%)", type="fft")

plot_sector_3(SEreg, ylab="Regional (total) sector effects (%)", type="hist",
    breaks=20)
title(main="Binning")
plot_sector_3(SEloc, ylab="Local (within footprint) sector effects (%)", type="hist",
    breaks=20)
par(op)

## load
library(mefa4)
load("e:/peter/AB_data_v2017/data/reporting/vpalnts-2017-11-16.Rdata")

## subset here to the region in question

P <- t(sapply(colnames(Y), function(i) colMeans(A * Y[,i])))
SE <- 100 * P[,colnames(P) != "NATIVE"] / P[,"NATIVE"]
SI <- 100 * (1 - colMeans(Y))

P1 <- rowSums(P[,colnames(P) != "NATIVE"])
P0 <- P[,"NATIVE"]
SIv <- 100 * pmin(P1, P0) / pmax(P1, P0)

#plot(SI[Z$NN], SIv[Z$NN])

tmp <- rbind(data.frame(SI=SI[Z$NN], SE[Z$NN,], NN=1),
    data.frame(SI=SI[!Z$NN], SE[!Z$NN,], NN=0))
tmp <- tmp[,colnames(tmp) != "Misc"]
#write.csv(tmp, file="e:/peter/AB_data_v2017/data/reporting/nn-plant-sector-effects.csv")

NN <- VP_SectorAbundance.Curr / rowSums(VP_SectorAbundance.Curr)
DD <- t(t(NN) / CS)

SE2 <- 100 * DD[,colnames(DD) != "Native"] / DD[,"Native"]
SE2 <- SE2[rownames(SE),]

op <- par(mfrow=c(2,3))
for (i in 2:6) {
plot(SE[Z$NN,i], SE2[Z$NN,i], xlab="Data based", ylab="Prediction based",
     ylim=c(0, 4000), xlim=c(0, 500), main=colnames(SE)[i])
abline(0,1,col="grey")
abline(h=100,v=100,col="grey", lty=2)
}
par(op)

op <- par(mfrow=c(2,1))
boxplot(SE[Z$NN,], ylim=c(0,2000))
abline(h=100)
boxplot(SE2[Z$NN,], ylim=c(0,10000))
abline(h=100)
par(op)





## producing some pdf's

Area <- 100*CS
pdf("sector-effect_local-regional.pdf", onefile=TRUE, width=21, height=7)
for (spp in rownames(Z)[!Z$NN]) {
    cat(spp, "\n")
    flush.console()
    Curr <- VP_SectorAbundance.Curr[spp,]
    Ref <- VP_SectorAbundance.Ref[spp,]
    op <- par(mfrow=c(1,3))
    plot_sector_1(Curr, Ref, Area, main=paste("Combined -", spp))
    plot_sector_2(Curr, Ref, regional=TRUE, main=paste("Regional (Total) -", spp))
    plot_sector_2(Curr, Ref, regional=FALSE, main=paste("Local -", spp))
    par(op)
}
dev.off()


## making some polys:

library(rgdal)
library(rgeos)
library(sp)

#osa <- readOGR("e:/peter/AB_data_v2017/data/raw/xy/OilSandsBoundaries.gdb",layer='myFeatureClass')

## shape files for boundary
crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/nsr")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
setwd(od)
AB <- spTransform(AB, crs)
AB0 <- gUnaryUnion(AB, rep(1, nrow(AB@data)))
AB0 <- as(AB0, "SpatialPolygonsDataFrame")
AB0@data <- data.frame(Province="Alberta")
100*object.size(AB0)/object.size(AB)
writeOGR(AB0, "AB_bound.geojson", layer="Alberta", driver="GeoJSON")

crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/osa")
OSA <- readOGR(".", "JOSM_Boundary_Sept2012") # rgdal
setwd(od)
OSA <- spTransform(OSA, crs)
OSA0 <- gUnaryUnion(OSA, rep(1, nrow(OSA@data)))
OSA0 <- as(OSA0, "SpatialPolygonsDataFrame")
OSA0@data <- data.frame(Region="OSA")
writeOGR(OSA0, "OSA_bound.geojson", layer="OSA", driver="GeoJSON")


