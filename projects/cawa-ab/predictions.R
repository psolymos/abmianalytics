library(mefa4)

shf <- FALSE

ROOT <- "d:/abmi/AB_data_v2016"

## surrounding hf
if (shf) {
    OUTDIR1 <- "d:/abmi/sppweb2015/birds-pred-shf-1/"
    OUTDIRB <- "d:/abmi/sppweb2015/birds-pred-shf-B/"
} else {
    OUTDIR1 <- "d:/abmi/sppweb2015/birds-pred-1/"
    OUTDIRB <- "d:/abmi/sppweb2015/birds-pred-B/"
}
#OUTDIR <- "d:/abmi/AB_data_v2016/out/birds/figs"

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
## area in ha
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

q <- 1
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

library(RColorBrewer)
br <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col <- rev(brewer.pal(10, "RdYlGn"))
Colxfun <- colorRampPalette(Col, space = "rgb")
Col2 <- Colxfun(100)
level <- 0.9

spp <- "CAWA"

load(file.path(ROOT, "out", "birds", "results", "cawa", "predB", paste0(regs[1], ".Rdata")))
rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
pxNcr0 <- pxNcrB
#reg_veg_B <- array(NA, c(nrow(ch2veg), 240))
hbNcrB[is.na(hbNcrB)] <- 0 # reason: area=0 div by 0
veg_B <- hbNcrB * Aveg[1,]
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "birds", "results", "cawa", "predB", paste0(regs[i], ".Rdata")))
    rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
    pxNcr0 <- rbind(pxNcr0, pxNcrB)
    hbNcrB[is.na(hbNcrB)] <- 0 # reason: area=0 div by 0
    veg_B <- veg_B + (hbNcrB * Aveg[i,])
}

km <- fstatv(pxNcr0, level=level)
summary(km$Aland)
km$Aland <- 1 - kgrid[rownames(pxNcr0), "pWater"]
km_tot <- colSums(100 * pxNcr0 * km$Aland) # ha to km^2, taking into account land area in pixel
save(km, km_tot, veg_B,
    file=file.path(ROOT, "out", "birds", "results", "cawa", "cawa-km-predB.Rdata"))

load(file.path(ROOT, "out", "birds", "results", "cawa", "cawa-km-predB.Rdata"))
km2 <- km[match(rownames(kgrid), rownames(km)),]
## pretty close after accounting for Aland
fstat(colSums(veg_B), .9)
fstat(km_tot, .9)


    TYPE <- "N"
    NAM <- as.character(tax[spp, "English_Name"])

    cr <- km2$Mean
    SD <- km2$SD
    CoV <- SD / cr

    cr[is.na(cr)] <- 0
    SD[is.na(SD)] <- 0
    SD[SD==0] <- 0.000001

    summary(cr)
    summary(SD)
    summary(CoV)

    qcr <- quantile(cr, q)
    cr[cr>qcr] <- qcr
    crmean <- mean(cr)

    Max <- max(qcr)
    cr0 <- cr
    cr <- pmin(100, ceiling(99 * sqrt(cr0 / Max))+1)
    #cr <- pmin(100, ceiling(99 * (cr0 / Max))+1)

Lc_quantile <- function (xx, probs=seq(0, 1, 0.1), type=c("L","p")) {
    xx <- xx[!is.na(xx)]
    o <- order(xx)
    x <- cumsum(xx[o]) / sum(xx)
    if (type=="L")
        q <- probs
    if (type=="p")
        q <- quantile(x, probs=probs, na.rm=TRUE)
    xxo <- xx[o]
    i <- sapply(q, function(z) min(xxo[x >= z]))
    i
}
#probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
probs <- 0:100/100
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
Col <- rev(brewer.pal(5, "RdYlBu"))

br <- Lc_quantile(cr0, probs=probs, type="L")
br <- unique(br)
C1x <- Col1fun(length(br))
zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(cr0)) else as.integer(cut(cr0, breaks=br))
zval[is.na(zval)] <- 1

    iiii <- !(kgrid$POINT_Y > 50 & kgrid$NRNAME != "Grassland")

## raster
library(raster)
source("~/repos/abmianalytics/R/maps_functions.R")
rt <- raster(file.path(ROOT, "data", "kgrid", "AHM1k.asc"))
x_cr_mean <- as_Raster(as.factor(kgrid$Row), as.factor(kgrid$Col),
    ifelse(!iiii, round(km2$Mean, 5), NA), rt)
#writeRaster(x_cr_mean, "CAWA_AB_2016-08-10.asc", overwrite=TRUE)
writeRaster(x_cr_mean, "CAWA_AB_2016-08-12.tif", overwrite=TRUE)



    fname <- file.path("c:/Users/Peter/Dropbox/josm/cawa-jeff/revision",
        "cawa-map-cr-cov.png")
    png(fname, width=W*2, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,2))
#    plot(kgrid$X, kgrid$Y, col=C1[cr], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    plot(kgrid$X, kgrid$Y, col=C1x[zval], pch=15, cex=cex, ann=FALSE, axes=FALSE)
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
    #pv <- as.character(round((c(0, 0.25, 0.5, 0.75, 1)*Max)^2, 3))
    pv <- round(c(0, br[round(seq(0.1, 1, 0.1)*length(br))]), 3)
    pv <- format(pv, nsmall = 3)
    zzz <- 10000
    text(240000, 5730000, "Density (males / ha)")
    for (i in 1:length(pv))
        text(240000+zzz, 5450000 + seq(0, 1, 0.1)[i]*(5700000-5450000), pv[i], cex=1)
    #par(op)
    #dev.off()

    tmp <- aggregate(CoV, list(kgrid$NSRNAME), mean, na.rm=TRUE)
    tmp <- tmp$x[match(kgrid$NSRNAME, tmp[,1])]
    CoV[is.na(CoV)] <- tmp[is.na(CoV)]

    covC <- CoV
    zval <- pmin(100, ceiling(99 * (covC / 2))+1)
    #zval <- as.integer(cut(covC, breaks=br))
    #fname <- file.path("c:/Users/Peter/Dropbox/josm/cawa-jeff/revision",
    #    "cawa-map-cov.png")
    #png(fname, width=W, height=H)
    #op <- par(mar=c(0, 0, 4, 0) + 0.1)
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
br <- seq(0, sqrt(max(SD)), len=9)^2
Col <- brewer.pal(9, "Reds")
    zval <- as.integer(cut(SD, breaks=br))
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
    br2 <- round(br, 3)
    br2 <- format(br2, nsmall = 3)
    TEXT <- paste0(br2[-length(br2)], "-", br2[-1])
    INF <- grepl("Inf", TEXT)
    if (any(INF))
        TEXT[length(TEXT)] <- paste0(">", br2[length(br2)-1])
    TITLE <- paste0("Prediction Std. Error")
    legend("bottomleft", border=rev(Col), fill=rev(Col), bty="n", legend=rev(TEXT),
                title=TITLE, cex=legcex*0.8)
    par(op)
    dev.off()

## pop size by habitat & HF in it

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv$age <- as.character(tv$AGE)
tv$age[!(tv$age %in% c("5","6","7","8","9"))] <- ""
tv$age[(tv$age %in% c("5","6","7","8","9"))] <- "O"
tv$hab <- as.character(interaction(tv$Type, tv$age, sep=""))
tv$hab[!is.na(tv$HF)] <- as.character(tv$Sector2[!is.na(tv$HF)])
tv$hab[rownames(tv) %in% c("NonVeg","Water",
    "BorrowpitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals",
    "RailHardSurface","RoadHardSurface")] <- "XXX"

tv$use_tr <- as.character(tv$VEGAGE_use)
tv$use_tr[!is.na(tv$HF)] <- as.character(tv$VEGHFAGE[!is.na(tv$HF)])

ch2veg$rfhab <- as.factor(tv$hab[match(ch2veg$rf, tv$use_tr)])
ch2veg$crhab <- as.factor(tv$hab[match(ch2veg$cr, tv$use_tr)])
ch2veg$A <- colSums(Aveg) # in ha
summary(ch2veg)

xthf <- as.matrix(Xtab(A ~ rfhab + crhab, ch2veg))

Nveg <- groupSums(veg_B, 1, ch2veg$crhab)
AA <- matrix(ch2veg$A, nrow(ch2veg), 240)
DD <- veg_B / AA
DD[is.na(DD)] <- 0
Aveg <- groupSums(AA, 1, ch2veg$crhab) # in ha
#Dveg <- groupMeans(DD, 1, ch2veg$crhab)
Dveg <- Nveg / Aveg

Nveg_stat <- fstatv(Nveg, 0.9)
Dveg_stat <- fstatv(Dveg, 0.9)

TAB <- data.frame(D=round(Dveg_stat[,c(1,5,6)], 4), # males / ha
    N=round(Nveg_stat[,c(1,5,6)]/1000,3), # 1K males
    A=round(Aveg[,1]/100)) # km^2
TAB$Prc <- round(100*TAB$N.Mean / sum(TAB$N.Mean), 1)
TAB <- TAB[rownames(TAB) != "XXX",]
TAB <- TAB[order(TAB$D.Mean, decreasing=TRUE),]

TAB <- rbind(TAB, Total=NA)

TAB["Total",1:3] <- round(fstat(colSums(Nveg[rownames(Dveg) != "XXX",]) / colSums(Aveg[rownames(Dveg) != "XXX",]))[c(1,3,4)], 4)
TAB["Total",4:6] <- round(fstat(colSums(Nveg[rownames(Dveg) != "XXX",]), 0.9)[c(1,3,4)]/1000, 3)
TAB["Total",7] <- round(sum(Aveg[rownames(Dveg) != "XXX",1]/100))
TAB["Total",8] <- round(100,1)


write.csv(TAB, file.path("c:/Users/Peter/Dropbox/josm/cawa-jeff/revision",
        "cawa-pop.csv"))

