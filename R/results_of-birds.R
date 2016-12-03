## old-forest birds guild results

library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"

level <- 0.9

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
#    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

slt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(slt) <- slt$AOU
slt$comments <- NULL

## guild and members
gname <- "old-forest-birds"
NAM <- "Old Forest Birds"
gspp <- rownames(slt)[slt$oldforest == 1 & slt$map.pred]

## coefs

load(file.path(ROOT, "tables", "res_coef.Rdata"))

tp <- combine_spp_coefs(res_coef, gspp)

library(plotrix)
par(las=2, mar=c(8,4,1,1))
ladderplot(t(tp$matN), pch=NA, col=1:ncol(tp$matN))

fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-veghf-north-all.pdf"))
pdf(file=fname,width=10,height=5,onefile=TRUE)
for (i in gspp)
fig_veghf(res_coef[[i]]$veg, i, res_coef[[i]]$max[1])
dev.off()

if (max(tp$max, na.rm=TRUE) > 3*min(tp$max, na.rm=TRUE)) {
    MAXn <- tp$max[1]
    MAXs <- tp$max[2]
} else {
    MAXn <- max(tp$max, na.rm=TRUE)
    MAXs <- max(tp$max, na.rm=TRUE)
}
if (is.na(MAXn))
    MAXn <- max(tp$max, na.rm=TRUE)
if (is.na(MAXs))
    MAXs <- max(tp$max, na.rm=TRUE)

prn <- tp$veg
NDAT <- length(tp$species$north)
## veghf
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-veghf-north.png"))
png(file=fname,width=1500,height=700)
fig_veghf(prn, paste0(NAM, " (n = ", NDAT, " species)"), ymax=MAXn)
dev.off()
## linear
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-linear-north.png"))
png(file=fname,width=350,height=400)
fig_linear(attr(prn, "linear"), paste0(NAM, "\nNorth (n = ", NDAT, " species)"))
dev.off()

prs <- tp$soil
NDAT <- length(tp$species$south)
## treed
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-soilhf-treed-south.png"))
png(file=fname,width=500,height=450)
fig_soilhf(prs$treed,
    paste0(NAM, ", South, Treed (n = ", NDAT, " species)"),
    ymax=MAXs)
dev.off()
## nontreed
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-soilhf-nontreed-south.png"))
png(file=fname,width=500,height=450)
fig_soilhf(prs$nontreed,
    paste0(NAM, ", South, Non-treed (n = ", NDAT, " species)"),
    ymax=MAXs)
dev.off()
## linear
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-linear-south.png"))
png(file=fname,width=350,height=400)
fig_linear(prs$linear,
    paste0(NAM, "\nSouth (n = ", NDAT, " species)"))
dev.off()

## sector effects

load(file.path(ROOT, "tables", "sector-effects.Rdata"))

fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-sector-north-all.pdf"))
pdf(fname, width=6, height=6, onefile=TRUE)
for (i in gspp)
    plot_seff(seff_res[[i]]$N, NAM=i, TAG="", WHERE="North")
dev.off()

gseff <- combine_spp_seff(seff_res, gspp)
gseff2 <- combine_spp_seff(seff_res, gspp, scale=TRUE)

fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-sector-north.png"))
png(fname, width=600, height=600)
plot_seff(gseff$N, NAM=NAM, TAG="", WHERE="North", CL1=gseff$Nmin, CL2=gseff$Nmax)
dev.off()

fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-sector-south.png"))
png(fname, width=600, height=600)
plot_seff(gseff$S, NAM=NAM, TAG="", WHERE="North", CL1=gseff$Smin, CL2=gseff$Smax)
dev.off()

fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-sector-north-scaled.png"))
png(fname, width=600, height=600)
plot_seff(gseff2$N, NAM=NAM, TAG="", WHERE="South", CL1=gseff$Nmin, CL2=gseff$Nmax)
dev.off()

fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-sector-south-scaled.png"))
png(fname, width=600, height=600)
plot_seff(gseff2$S, NAM=NAM, TAG="", WHERE="South", CL1=gseff$Smin, CL2=gseff$Smax)
dev.off()

## det map

source("~/repos/abmianalytics/R/wrsi_functions.R")

load(file.path(ROOT, "data", "data-wrsi.Rdata"))

TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")
DATw$NSRxLUF <- interaction(DATw$NSRNAME, DATw$LUF_NAME, sep="_", drop=TRUE)

load(file.path("e:/peter/AB_data_v2016/out", "kgrid", "kgrid_table.Rdata"))
regs <- levels(kgrid$LUFxNSR)
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

col1 <- c("#C8FBC8","#C8E6FA","#F5E6F5","#FFDCEC","#FFE6CD","#FFF1D2")[match(kgrid$NRNAME,
    c("Boreal","Foothills","Rocky Mountain","Canadian Shield","Parkland","Grassland"))]

library(raster)
library(sp)
library(rgdal)
city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS(paste0("+proj=longlat +datum=WGS84 ",
    "+ellps=WGS84 +towgs84=0,0,0"))
city <- as.data.frame(spTransform(city, CRS(paste0("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 ",
    "+x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))))
xyw <- as.matrix(kgrid[kgrid$pWater >= 0.99,c("X","Y")])
blank <- matrix(0, 0, 2)

    ## mask (e.g. for CAWA project)
    iii <- rep(TRUE, nrow(DATw))
    #iii <- DATw$POINT_Y > 50 & DATw$NRNAME != "Grassland"

    S <- rowSums(YYw[,intersect(colnames(YYw), gspp)] > 0)
    xy0 <- as.matrix(DATw[S == 0 & iii, c("X","Y")])
    xy1 <- as.matrix(DATw[S > 0 & iii, c("X","Y")])

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-det.png"))
	png(file=fname, width=600, height=1000)
    plot(kgrid$X, kgrid$Y, pch=15, cex=0.2, col=col1, axes=FALSE, ann=FALSE)
    #iiii <- !(kgrid$POINT_Y > 50 & kgrid$NRNAME != "Grassland")
    #points(kgrid$X[iiii], kgrid$Y[iiii], pch=15, cex=0.2, col="grey")
    points(xyw, pch=15, cex=0.2, col=rgb(0.3,0.45,0.9))
    #points(xy0, pch=19, cex=0.5, col="red3")
    points(xy0, pch=3, cex=0.6, col="red2")
    points(xy1, pch=19, cex=pmin(1.5, 0.2+S[S > 0]/5), col="red4")
    mtext(paste0(NAM), line=1, side=3, adj=0.5, cex=1.4, col="grey40")
    points(city, pch=18, col="grey10")
    text(city, rownames(city), cex=1, adj=-0.1, col="grey10")
    legend("bottomleft", pch=c(15,15,15,15,15,15, NA, 3,19),
        col=c("#C8FBC8","#C8E6FA","#F5E6F5","#FFDCEC","#FFE6CD","#FFF1D2", NA, "red2", "red4"),
        #col=c("#C8FBC8","#C8E6FA","#F5E6F5","#FFDCEC","#FFE6CD",NA, "grey", NA, "red2", "red4"),
        legend=c(
        "Boreal","Foothills","Rocky Mountain","Canadian Shield","Parkland","Grassland",
        ## "Boreal","Foothills","Rocky Mountain","Canadian Shield","Parkland",NA,"Outside of study area",
        NA, "Survey location","Detection"), title="Natural Regions",
        cex=1.2, pt.cex=c(4,4,4,4,4,4, NA, 1.5,2), bty="n")
	dev.off()

## load stuff -- need to bring in un-normalized results (abundance is not comparable)

cex <- 0.25
legcex <- 1.5

Col1 <- rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4"))  # Colour gradient for reference and current
Col1fun <- colorRampPalette(Col1, space = "rgb") # Function to interpolate among these colours for reference and current
C1 <- Col1fun(100)
Col2 <- c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221")  # Colour gradient for difference map
Col2fun <- colorRampPalette(Col2, space = "rgb") # Function to interpolate among these colours for difference map
C2 <- Col2fun(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude
CSI <- colorRampPalette(c("red","orange","yellow","darkgreen"), space = "rgb")(100)
CSR <- colorRampPalette(c("yellow","purple"), space = "rgb")(100)

q <- 0.99
H <- 1000
W <- 600

mcr <- matrix(0, nrow(kgrid), length(gspp))
dimnames(mcr) <- list(rownames(kgrid), gspp)
mrf <- mcr
for (spp in gspp) {
    cat(spp, "\n");flush.console()
    load(paste0("e:/peter/AB_data_v2016/out/birds/pred1cmb/", spp, ".Rdata"))
    stopifnot(all(rownames(kgrid)==rownames(km)))
    km <- data.frame(km)

    TYPE <- "C" # combo
    if (!slt[spp, "veghf.north"])
    #if (!(spp %in% fln))
        TYPE <- "S"
    if (!slt[spp, "soilhf.south"])
    #if (!(spp %in% fls))
        TYPE <- "N"

    wS <- 1-kgrid$pAspen
    if (TYPE == "S")
        wS[] <- 1
    if (TYPE == "N")
        wS[] <- 0
    wS[kgrid$useS] <- 1
    wS[kgrid$useN] <- 0

    cr <- wS * km$CurrS + (1-wS) * km$CurrN
    rf <- wS * km$RefS + (1-wS) * km$RefN

    qcr <- quantile(cr, q, na.rm=TRUE)
    cr[!is.na(cr) & cr > qcr] <- qcr
    qrf <- quantile(rf, q, na.rm=TRUE)
    rf[!is.na(rf) & rf > qrf] <- qrf

    if (TYPE=="N") {
        cr[!kgrid$useN] <- NA
        rf[!kgrid$useN] <- NA
    }
    if (TYPE=="S") {
        cr[!kgrid$useS] <- NA
        rf[!kgrid$useS] <- NA
    }

    mcr[,spp] <- cr
    mrf[,spp] <- rf
}
range(mcr, na.rm=TRUE)
range(mrf, na.rm=TRUE)

## richness map
Scr <- rowSums(1-exp(-mcr), na.rm=TRUE)
Srf <- rowSums(1-exp(-mrf), na.rm=TRUE)
Dcr <- rowSums(mcr, na.rm=TRUE)
Drf <- rowSums(mrf, na.rm=TRUE)
SItot <- 100 * pmin(Dcr, Drf, na.rm=TRUE) / pmax(Dcr, Drf, na.rm=TRUE)
SImean <- 100 * rowMeans(pmin(mcr, mrf) / pmax(mcr, mrf), na.rm=TRUE)
#plot(SItot, SImean)


#    SI <- round(100 * pmin(cr, rf) / pmax(cr, rf))
#    SI[is.na(SI)] <- 100 # 0/0 is defined as 100 intact
#    cr0 <- cr
#    rf0 <- rf
#    SI0 <- SI

SItot[is.na(SItot)] <- 100
SItot <- round(SItot)
SItot[SItot < 1] <- 1 # this is only for mapping

SImean[is.na(SImean)] <- 100
SImean <- round(SImean)
SImean[SImean < 1] <- 1 # this is only for mapping

Max <- max(Dcr, Drf, na.rm=TRUE)
ddf <- (Dcr-Drf) / Max
ddf <- sign(ddf) * abs(ddf)^0.5
ddf <- pmin(200, ceiling(99 * ddf)+100)
ddf[ddf==0] <- 1
dcr <- pmin(100, ceiling(99 * sqrt(Dcr / Max))+1)
drf <- pmin(100, ceiling(99 * sqrt(Drf / Max))+1)

Max2 <- max(Scr, Srf, na.rm=TRUE)
sdf <- (Scr-Srf) / Max2
sdf <- sign(sdf) * abs(sdf)^0.5
sdf <- pmin(200, ceiling(99 * sdf)+100)
sdf[sdf==0] <- 1
scr <- pmin(100, ceiling(99 * sqrt(Scr / Max2))+1)
srf <- pmin(100, ceiling(99 * sqrt(Srf / Max2))+1)

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-SIsummed.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=CSI[SItot], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nIntactness (summed)"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=CSI[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=CSI[i], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 0.5*(5450000 + 5700000), "50%")
    text(240000, 5700000, "100%")
    ## test NAs
#    with(kgrid[is.na(SI) & kgrid$pWater <= 0.99,], points(X, Y, col="black", pch=15, cex=cex))
    par(op)
    dev.off()

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-SIaveraged.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=CSI[SImean], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nIntactness (averaged)"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=CSI[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=CSI[i], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 0.5*(5450000 + 5700000), "50%")
    text(240000, 5700000, "100%")
    ## test NAs
#    with(kgrid[is.na(SI) & kgrid$pWater <= 0.99,], points(X, Y, col="black", pch=15, cex=cex))
    par(op)
    dev.off()

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-total-ref.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[drf], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nReference abundance"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=C1[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=C1[i], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 0.5*(5450000 + 5700000), "50%")
    text(240000, 5700000, "100%")
    par(op)
    dev.off()

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-total-curr.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[dcr], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nCurrent abundance"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=C1[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=C1[i], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 0.5*(5450000 + 5700000), "50%")
    text(240000, 5700000, "100%")
    par(op)
    dev.off()

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-total-diff.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C2[ddf], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nDifference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:200) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=C2[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/200
        segments(190000, 5450000+j, 220000, 5450000+j, col=C2[i], lwd=2, lend=2)
    }
    text(245000, 5450000, "-100%")
    text(245000, 0.5*(5450000 + 5700000), "0%")
    text(245000, 5700000, "+100%")
    par(op)
    dev.off()


    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-srich-ref.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=CSR[srf], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nSpecies Richness, Reference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=CSR[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=CSR[i], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 0.5*(5450000 + 5700000), "50%")
    text(240000, 5700000, "100%")
    par(op)
    dev.off()

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-srich-curr.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=CSR[scr], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nSpecies Richness, Current"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=CSR[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=CSR[i], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 0.5*(5450000 + 5700000), "50%")
    text(240000, 5700000, "100%")
    par(op)
    dev.off()

    fname <- file.path(ROOT, "guilds", gname,
        paste0(gname, "-map-srich-diff.png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C2[sdf], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nSpecies Richness, Difference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:200) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=C2[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/200
        segments(190000, 5450000+j, 220000, 5450000+j, col=C2[i], lwd=2, lend=2)
    }
    text(245000, 5450000, "-100%")
    text(245000, 0.5*(5450000 + 5700000), "0%")
    text(245000, 5700000, "+100%")
    par(op)
    dev.off()


