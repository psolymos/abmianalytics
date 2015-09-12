library(mefa4)

shf <- TRUE

ROOT <- "c:/p/AB_data_v2015"

## surrounding hf
if (shf) {
    OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-shf-1/"
    OUTDIRB <- "e:/peter/sppweb2015/birds-pred-shf-B/"
} else {
    OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-1/"
    OUTDIRB <- "e:/peter/sppweb2015/birds-pred-B/"
}

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
#source("~/repos/bragging/R/glm_skeleton.R")
#source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
regs <- levels(kgrid$LUFxNSR)
useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
useS <- kgrid$NRNAME == "Grassland"

e <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-full-withrevisit.Rdata"), envir=e)
tax <- e$TAX
rm(e)
tax$file <- nameAlnum(as.character(tax$English_Name), "mixed", "")

## model for species
fl <- list.files(file.path(ROOT, "out", "birds", "results"))
fln <- fl[grep("-north_", fl)]
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- fl[grep("-south_", fl)]
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)
SPP <- union(fln, fls)

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

spp <- "CAWA"

for (spp in SPP) {

cat(spp, "\n");flush.console()

load(file.path(OUTDIR1, spp, paste0(regs[1], ".Rdata")))
rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
rownames(pxScr1) <- rownames(pxSrf1) <- names(Cells)
pxNcr <- pxNcr1
pxNrf <- pxNrf1
pxScr <- pxScr1
pxSrf <- pxSrf1
pSoil <- pSoil1
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(OUTDIR1, spp, paste0(regs[i], ".Rdata")))
    rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
    rownames(pxScr1) <- rownames(pxSrf1) <- names(Cells)
    pxNcr <- rbind(pxNcr, pxNcr1)
    pxNrf <- rbind(pxNrf, pxNrf1)
    pxScr <- rbind(pxScr, pxScr1)
    pxSrf <- rbind(pxSrf, pxSrf1)
    pSoil <- c(pSoil, pSoil1)
}

pxNcr <- pxNcr[rownames(kgrid),]
pxNrf <- pxNrf[rownames(kgrid),]
pxScr <- pxScr[rownames(kgrid),]
pxSrf <- pxSrf[rownames(kgrid),]
pSoil <- pSoil[rownames(kgrid)]

if (FALSE) {
wS <- 1-kgrid$pAspen
if (!NSest["north"])
    wS[] <- 1
if (!NSest["south"])
    wS[] <- 0

pxcr <- wS * pxScr + (1-wS) * pxNcr
pxcr[useN] <- pxNcr[useN]
pxcr[useS] <- pxScr[useS]
pxrf <- wS * pxSrf + (1-wS) * pxNrf
pxrf[useN] <- pxNrf[useN]
pxrf[useS] <- pxSrf[useS]

pxcr <- pxNcr
pxrf <- pxNrf

if (!NSest["north"]) {
    pxcr[useN] <- NA
    pxrf[useN] <- NA
}
if (!NSest["south"]) {
    pxcr[useS] <- NA
    pxrf[useS] <- NA
}

km <- data.frame(LinkID=kgrid$Row_Col,
    Ref=pxrf,
    Curr=pxcr)
}

#pxNcr[useS] <- NA
#pxNrf[useS] <- NA
#pxScr[useN] <- NA
#pxSrf[useN] <- NA
km <- data.frame(LinkID=kgrid$Row_Col,
    RefN=pxNrf,
    CurrN=pxNcr,
    RefS=pxSrf,
    CurrS=pxScr)
#NAM <- as.character(tax[spp, "English_Name"])
write.csv(km, row.names=FALSE,
    paste0("e:/peter/sppweb2015/birds-pred/", paste0(as.character(tax[spp, "file"]), ".csv")))
}

## plots

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

#SPP <- union(fln, fls)
SPP <- fln

#spp <- "CAWA"
for (spp in SPP) {
    km <- read.csv(paste0("e:/peter/sppweb2015/birds-pred/", 
        paste0(as.character(tax[spp, "file"]), ".csv")))

    cr <- km$CurrN
    rf <- km$RefN
    qcr <- quantile(cr, q)
    cr[cr>qcr] <- qcr
    qrf <- quantile(rf, q)
    rf[rf>qrf] <- qrf
    Max <- max(qcr, qrf)

    df <- (cr-rf) / Max
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(cr / Max))+1)
    rf <- pmin(100, ceiling(99 * sqrt(rf / Max))+1)
    range(cr)
    range(rf)
    range(df)

    NAM <- as.character(tax[spp, "English_Name"])
    TAG <- "-shf"
    
    cat("\n")
    cat(spp, "\t");flush.console()
    cat("rf\t");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-rf/",
        as.character(tax[spp, "file"]), TAG, ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[rf], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Grassland",], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nReference abundance"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    cat("cr\t");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-cr/",
        as.character(tax[spp, "file"]), TAG, ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[cr], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Grassland",], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nCurrent abundance"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    cat("df\n");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-df/",
        as.character(tax[spp, "file"]), TAG, ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C2[df], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Grassland",], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nDifference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()
}

## sector effect

ch2veg <- t(sapply(strsplit(colnames(trVeg), "->"), 
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2veg <- data.frame(ch2veg)
colnames(ch2veg) <- c("rf","cr")
rownames(ch2veg) <- colnames(Aveg)

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

for (spp in fln) {

cat(spp, "\n");flush.console()

load(file.path(OUTDIR1, spp, paste0(regs[1], ".Rdata")))
hbNcr <- hbNcr1[,1]
hbNrf <- hbNrf1[,1]
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(OUTDIR1, spp, paste0(regs[i], ".Rdata")))
    hbNcr <- rbind(hbNcr, hbNcr1[,1])
    hbNrf <- rbind(hbNrf, hbNrf1[,1])
}
dimnames(hbNcr) <- dimnames(hbNrf) <- list(regs, colnames(Aveg))
hbNcr[is.na(hbNcr)] <- 0
hbNrf[is.na(hbNrf)] <- 0
hbNcr <- hbNcr * Aveg
hbNrf <- hbNrf * Aveg
ThbNcr <- colSums(hbNcr[lxn$N,])
ThbNrf <- colSums(hbNrf[lxn$N,])
df <- (ThbNcr - ThbNrf) / sum(ThbNrf)
dA <- Xtab(AvegN ~ rf + cr, ch2veg)
dN <- Xtab(df ~ rf + cr, ch2veg)
dA <- colSums(as.matrix(groupSums(dA[,rownames(tv)], 2, tv$Sector)))
dN <- colSums(as.matrix(groupSums(dN[,rownames(tv)], 2, tv$Sector)))
U <- dN/dA
seff <- cbind(dA=dA, dN=dN, U=U)[-1,]
(sum(hbNcr)-sum(hbNrf))/sum(hbNrf)

}


