library(mefa4)

shf <- FALSE

ROOT <- "e:/peter/AB_data_v2016"

## surrounding hf
if (shf) {
    OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-shf-1/"
    OUTDIRB <- "e:/peter/sppweb2015/birds-pred-shf-B/"
} else {
    OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-1/"
    OUTDIRB <- "e:/peter/sppweb2015/birds-pred-B/"
}
#OUTDIR <- "e:/peter/AB_data_v2016/out/birds/figs"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
#source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/abmianalytics/R/maps_functions.R")
regs <- levels(droplevels(kgrid$LUFxNSR[kgrid$NRNAME != "Grassland"]))
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

e <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-full-withrevisit.Rdata"), envir=e)
dat <- e$DAT
dat <- dat[dat$useOK,]
yy <- e$YY[rownames(dat),]
tax <- droplevels(e$TAX[colnames(yy),])
rm(e)

## model for species
fln <- list.files(file.path(ROOT, "out", "birds", "results", "north"))
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- list.files(file.path(ROOT, "out", "birds", "results", "south"))
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)
## need to update these once checking is done !!!!!!!!!!!!!!!!!!

load(file.path(ROOT, "out", "transitions", paste0(regs[1], ".Rdata")))
Aveg <- rbind(colSums(trVeg))
rownames(Aveg) <- regs[1]
colnames(Aveg) <- colnames(trVeg)
Asoil <- rbind(colSums(trSoil))
rownames(Asoil) <- regs[1]
colnames(Asoil) <- colnames(trSoil)

for (i in 2:length(regs)) {
    cat(regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "transitions", paste0(regs[i], ".Rdata")))
    Aveg <- rbind(Aveg, colSums(trVeg))
    rownames(Aveg) <- regs[1:i]
    Asoil <- rbind(Asoil, colSums(trSoil))
    rownames(Asoil) <- regs[1:i]
}
Aveg <- Aveg / 10^4
Asoil <- Asoil / 10^4


library(raster)
library(sp)
library(rgdal)
city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
city <- spTransform(city, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
city <- as.data.frame(city)

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

q <- 0.99
H <- 1000
W <- 600


## sector effect

ch2veg <- t(sapply(strsplit(colnames(trVeg), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2veg <- data.frame(ch2veg)
colnames(ch2veg) <- c("rf","cr")
rownames(ch2veg) <- colnames(Aveg)
ch2veg$uplow <- as.factor(ifelse(ch2veg$rf %in% c("BSpr0", "BSpr1", "BSpr2", "BSpr3",
    "BSpr4", "BSpr5", "BSpr6",
    "BSpr7", "BSpr8", "BSpr9", "BSprR", "Larch0", "Larch1", "Larch2", "Larch3",
    "Larch4", "Larch5", "Larch6", "Larch7", "Larch8", "Larch9", "LarchR",
    "WetGrassHerb", "WetShrub"), "lowland", "upland"))

ch2soil <- t(sapply(strsplit(colnames(trSoil), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2soil <- data.frame(ch2soil)
colnames(ch2soil) <- c("rf","cr")
rownames(ch2soil) <- colnames(Asoil)

lxn <- nonDuplicated(kgrid[,c("LUF_NAME","NRNAME","NSRNAME")], kgrid$LUFxNSR, TRUE)
lxn$N <- lxn$NRNAME != "Grassland" & lxn$NRNAME != "Rocky Mountain" &
    lxn$NRNAME != "Parkland" & lxn$NSRNAME != "Dry Mixedwood"
lxn$S <- lxn$NRNAME == "Grassland" | lxn$NRNAME == "Parkland" |
    lxn$NSRNAME == "Dry Mixedwood"
table(lxn$NRNAME, lxn$N)
table(lxn$NRNAME, lxn$S)
lxn <- lxn[regs,]
all(rownames(Aveg) == regs)
all(rownames(Asoil) == regs)
AvegN <- colSums(Aveg[lxn$N,])
AvegN <- AvegN / sum(AvegN)
AsoilS <- colSums(Asoil[lxn$S,])
AsoilS <- AsoilS / sum(AsoilS)

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv <- droplevels(tv[!is.na(tv$Sector),])
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")
ts <- droplevels(ts[!is.na(ts$Sector),])



library(RColorBrewer)
br <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col <- rev(brewer.pal(10, "RdYlGn"))
Colxfun <- colorRampPalette(Col, space = "rgb")
Col2 <- Colxfun(100)

spp <- "CAWA"

load(file.path(ROOT, "out", "birds", "results", "cawa", "predB", paste0(regs[1], ".Rdata")))
rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
pxNcr0 <- pxNcrB
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "birds", "results", "cawa", "predB", paste0(regs[i], ".Rdata")))
    rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
    pxNcr0 <- rbind(pxNcr0, pxNcrB)
}
#pxNcr <- pxNcr0[match(rownames(kgrid), rownames(pxNcr0)),]
#sum(is.na(pxNcr))
#sum(is.na(pxNcr0))
#pxNcr[is.na(pxNcr)] <- 0
#for (k in 1:ncol(pxNcr)) {
#    qN <- quantile(pxNcr[is.finite(pxNcr[,k]),k], q, na.rm=TRUE)
#    pxNcr[pxNcr[,k] > qN,k] <- qN
#    qS <- quantile(pxScr[is.finite(pxScr[,k]),k], q, na.rm=TRUE)
#    pxScr[pxScr[,k] > qS,k] <- qS
#}

level <- 0.9
km <- fstatv(pxNcr0, level=level)
save(km, file=file.path(ROOT, "out", "birds", "results", "cawa", "cawa-km-predB.Rdata"))

load(file.path(ROOT, "out", "birds", "results", "cawa", "cawa-km-predB.Rdata"))
km2 <- km[match(rownames(kgrid), rownames(km)),]


    TYPE <- "N"
    NAM <- as.character(tax[spp, "English_Name"])

    cr <- km2$Mean
    SD <- km2$SD
    CoV <- SD / cr

    cr[is.na(cr)] <- 0

    qcr <- quantile(cr, q)
    cr[cr>qcr] <- qcr
    crmean <- mean(cr)

    Max <- max(qcr)
    cr0 <- cr
    #cr <- pmin(100, ceiling(99 * sqrt(cr0 / Max))+1)
    cr <- pmin(100, ceiling(99 * (cr0 / Max))+1)

    iiii <- !(kgrid$POINT_Y > 50 & kgrid$NRNAME != "Grassland")

    fname <- file.path("c:/Users/Peter/Dropbox/josm/cawa-jeff/revision",
        "cawa-map-cr.png")
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[cr], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    points(kgrid$X[iiii], kgrid$Y[iiii], pch=15, cex=0.2, col="grey")
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
#        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,"A",col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=C1[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=C1[i], lwd=2, lend=2)
    }
    pv <- as.character(round(c(0, 0.25, 0.5, 0.75, 1)*Max, 3))
    zzz <- 10000
    text(240000, 5730000, "Density (males / ha)")
    text(240000+zzz, 5450000, "0.000")
    text(240000+zzz, 5450000 + 0.25*(5700000-5450000), pv[2])
    text(240000+zzz, 0.5*(5450000 + 5700000), pv[3])
    text(240000+zzz, 5450000 + 0.75*(5700000-5450000), pv[4])
    text(240000+zzz, 5700000, pv[5])
    par(op)
    dev.off()


    covC <- CoV
    zval <- pmin(100, ceiling(99 * (covC / 2))+1)
    #zval <- as.integer(cut(covC, breaks=br))
    fname <- file.path("c:/Users/Peter/Dropbox/josm/cawa-jeff/revision",
        "cawa-map-cov.png")
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=Col2[zval], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    points(kgrid$X[is.na(zval) & !iiii], kgrid$Y[is.na(zval) & !iiii], col=Col[10], pch=15, cex=cex)
    points(kgrid$X[iiii], kgrid$Y[iiii], pch=15, cex=0.2, col="grey")
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "N")
#        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
#    if (TYPE == "S")
#        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,"B",col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    for (i in 1:100) {
        #lines(c(190000, 220000), c(5450000, 5700000), col=Col2[i], lwd=2)
        j <- i * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+j, 220000, 5450000+j, col=Col2[i], lwd=2, lend=2)
    }
    zzz <- 8000
    text(240000, 5730000, "Coefficient of variation")
    text(240000+zzz, 5450000, " 0.0")
    text(240000+zzz, 5450000 + 0.25*(5700000-5450000), " 0.5")
    text(240000+zzz, 0.5*(5450000 + 5700000), " 1.0")
    text(240000+zzz, 5450000 + 0.75*(5700000-5450000), " 1.5")
    text(240000+zzz, 5700000, ">2.0")
    par(op)
    dev.off()


#    sdMax <- quantile(SD,q,na.rm=TRUE)
    #br2 <- seq(0, 0.02, len=length(Col))
    zval <- as.integer(cut(SD/crmean, breaks=br))
    fname <- file.path("c:/Users/Peter/Dropbox/josm/cawa-jeff/revision",
        "cawa-map-sd.png")
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=Col[zval], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "SE"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    br2 <- round(br * crmean, 3)
    TEXT <- paste0(br2[-length(br2)], "-", br2[-1])
    INF <- grepl("Inf", TEXT)
    if (any(INF))
        TEXT[length(TEXT)] <- paste0(">", br2[length(br2)-1])
    TITLE <- paste0("Prediction Std. Error")
    legend("bottomleft", border=rev(Col), fill=rev(Col), bty="n", legend=rev(TEXT),
                title=TITLE, cex=legcex*0.8)
    par(op)
    dev.off()



