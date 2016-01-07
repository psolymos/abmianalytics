library(mefa4)

shf <- FALSE

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
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useS <- kgrid$NRNAME == "Grassland"

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
#SPP <- union(fln, fls)

slt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(slt) <- slt$AOU

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

## CoV

ks <- kgrid[kgrid$Rnd10 <= 10,]
xy10 <- groupMeans(as.matrix(kgrid[,c("X","Y")]), 1, kgrid$Row10_Col10)

library(RColorBrewer)
br <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col <- rev(brewer.pal(10, "RdYlGn"))

## csv

spp <- "ALFL"
SPP <- union(fln, fls)
SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL")

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

km <- data.frame(LinkID=kgrid$Row_Col,
    RefN=pxNrf,
    CurrN=pxNcr,
    RefS=pxSrf,
    CurrS=pxScr)
if (any(is.na(km)))
    km[is.na(km)] <- 0
#NAM <- as.character(tax[spp, "English_Name"])
write.csv(km, row.names=FALSE,
    paste0("e:/peter/sppweb2015/birds-pred/", 
    paste0(as.character(tax[spp, "file"]), ".csv")))
}

## plots

#SPP <- union(fln, fls)
#SPP <- fln

#spp <- "CAWA"
res_luf <- list()
res_nsr <- list()
SPP <- as.character(slt$AOU[slt$map.pred])
for (spp in SPP) {

    cat(spp, "\n");flush.console()
    km <- read.csv(paste0("e:/peter/sppweb2015/birds-pred/", 
        paste0(as.character(tax[spp, "file"]), ".csv")))

if (FALSE) {
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

km <- data.frame(LinkID=kgrid$Row_Col,
    RefN=pxNrf,
    CurrN=pxNcr,
    RefS=pxSrf,
    CurrS=pxScr)
}


    TYPE <- "C" # combo
    if (!slt[spp, "veghf.north"])
        TYPE <- "S"
    if (!slt[spp, "soilhf.south"])
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
#    cr <- km$CurrN
#    rf <- km$RefN
#    cr <- km$CurrS
#    rf <- km$RefS
    qcr <- quantile(cr, q)
    cr[cr>qcr] <- qcr
    qrf <- quantile(rf, q)
    rf[rf>qrf] <- qrf

if (FALSE) {
    mat <- 100 * cbind(Ncurrent=cr, Nreference=rf) # ha to km^2
    rownames(mat) <- rownames(kgrid)
    res_luf[[spp]] <- groupSums(mat, 1, kgrid$LUF_NAME)
    res_nsr[[spp]] <- groupSums(mat, 1, kgrid$NSRNAME)
}

if (TRUE) {
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
    TAG <- ""
    
    cat("\n")
    cat(spp, "\t");flush.console()
    cat("rf\t");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-rf/",
        as.character(tax[spp, "file"]), TAG, ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[rf], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
#        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nReference abundance"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    cat("cr\t");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-cr/",
        as.character(tax[spp, "file"]), TAG, ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C1[cr], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
#        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nCurrent abundance"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    cat("df\n");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-df/",
        as.character(tax[spp, "file"]), TAG, ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=C2[df], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
#        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nDifference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()
}
}


save(res_nsr, res_luf, file=file.path(ROOT, "out", "birds", "luf_Nsummaries.Rdata"))

LUF <- list()
for (spp in names(res_luf)) {
    tmp <- res_luf[[spp]] / 10^6 # M males
    tmp <- t(matrix(tmp, 2*nrow(tmp), 1))
    colnames(tmp) <- paste(rep(colnames(res_luf[[1]]), each=ncol(tmp)/2), 
        rep(rownames(res_luf[[1]]), 2))
    LUF[[spp]] <- data.frame(Species=slt[spp, "species"], tmp)        
}
NSR <- list()
for (spp in names(res_nsr)) {
    tmp <- res_nsr[[spp]] / 10^6 # M males
    tmp <- t(matrix(tmp, 2*nrow(tmp), 1))
    colnames(tmp) <- paste(rep(colnames(res_nsr[[1]]), each=ncol(tmp)/2), 
        rep(rownames(res_nsr[[1]]), 2))
    NSR[[spp]] <- data.frame(Species=slt[spp, "species"], tmp)        
}
LUF <- do.call(rbind, LUF)
NSR <- do.call(rbind, NSR)
write.csv(LUF, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_Abundance_by_LUFregions.csv")
write.csv(NSR, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_Abundance_by_NaturalSubregions.csv")


## sector effects
seff_res <- list()
#seff_luf <- list()
#seff_ns <- list()
uplow <- list()
uplow_full <- list()
uplow_luf <- list()

## stuff to exclude
## add col to lxn
## subset counter for loop


for (spp in as.character(slt$AOU[slt$map.pred])) {

cat(spp, "\n");flush.console()

load(file.path(OUTDIR1, spp, paste0(regs[1], ".Rdata")))
hbNcr <- hbNcr1[,1]
hbNrf <- hbNrf1[,1]
hbScr <- hbScr1[,1]
hbSrf <- hbSrf1[,1]
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(OUTDIR1, spp, paste0(regs[i], ".Rdata")))
    hbNcr <- rbind(hbNcr, hbNcr1[,1])
    hbNrf <- rbind(hbNrf, hbNrf1[,1])
    hbScr <- rbind(hbScr, hbScr1[,1])
    hbSrf <- rbind(hbSrf, hbSrf1[,1])
}

dimnames(hbNcr) <- dimnames(hbNrf) <- list(regs, colnames(Aveg))
hbNcr[is.na(hbNcr)] <- 0
hbNrf[is.na(hbNrf)] <- 0
hbNcr <- hbNcr * Aveg
hbNrf <- hbNrf * Aveg

dimnames(hbScr) <- dimnames(hbSrf) <- list(regs, colnames(Asoil))
hbScr[is.na(hbScr)] <- 0
hbSrf[is.na(hbSrf)] <- 0
hbScr <- hbScr * Asoil
hbSrf <- hbSrf * Asoil


## combined upland/lowland N/S
crN <- groupSums(hbNcr, 2, ch2veg$uplow)
rfN <- groupSums(hbNrf, 2, ch2veg$uplow)
crN[lxn$NRNAME=="Grassland","lowland"] <- 0 
crN[lxn$NRNAME=="Grassland","upland"] <- rowSums(hbScr[lxn$NRNAME=="Grassland",])
rfN[lxn$NRNAME=="Grassland","lowland"] <- 0 
rfN[lxn$NRNAME=="Grassland","upland"] <- rowSums(hbSrf[lxn$NRNAME=="Grassland",])
uplo <- data.frame(Current=crN, Reference=rfN)
uplow_full[[spp]] <- data.frame(sppid=spp, lxn[,1:3], uplo)

## Exclude stuff here
r0 <- lxn$NSRNAME %in% c("Alpine","Lower Foothills",
    "Montane","Subalpine","Upper Foothills")
crN[r0,] <- 0
rfN[r0,] <- 0

## upland/lowland
cr <- colSums(crN)
rf <- colSums(rfN)
cr <- c(total=sum(cr), cr)
rf <- c(total=sum(rf), rf)
si <- 100 * pmin(cr, rf) / pmax(cr, rf)
si2 <- ifelse(cr > rf, 200-si, si)
uplow[[spp]] <- c(Ref=rf, Cur=cr, SI=si, SI200=si2)

cr <- groupSums(groupSums(hbNcr, 2, ch2veg$uplow), 1, lxn$LUF_NAME)
rf <- groupSums(groupSums(hbNrf, 2, ch2veg$uplow), 1, lxn$LUF_NAME)
cr <- cbind(total=rowSums(cr), cr)
rf <- cbind(total=rowSums(rf), rf)
si <- sapply(1:3, function(i) 100 * pmin(cr[,i], rf[,i]) / pmax(cr[,i], rf[,i]))
colnames(si) <- colnames(cr)
si2 <- ifelse(cr > rf, 200-si, si)

uplow_luf[[spp]] <- data.frame(ID=spp, Ref=round(rf), Cur=round(cr), 
    SI=round(si, 2), SI200=round(si2, 2))

ThbNcr <- colSums(hbNcr[lxn$N,])
ThbNrf <- colSums(hbNrf[lxn$N,])
df <- (ThbNcr - ThbNrf) / sum(ThbNrf)
dA <- Xtab(AvegN ~ rf + cr, ch2veg)
if (FALSE) {
    tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
    tv2 <- nonDuplicated(tv,Combined,TRUE)
    dA2 <- as.matrix(groupSums(dA[,rownames(tv2)], 2, tv2$Sector3))
    tv3 <- tv2[rownames(dA2),]
    dA2 <- as.matrix(groupSums(dA2, 1, tv3$Sector3))
    dA3 <- dA2[,c(c("Agriculture","Forestry","Energy","RuralUrban","Transportation"))]
    dA3 <- round(100*t(t(dA3) / colSums(dA3)), 1)
    dA3[c("Decid", "Mixwood", "UpConif", "LoConif", "Wet", "OpenOther"),]
}

dN <- Xtab(df ~ rf + cr, ch2veg)
#dA <- colSums(as.matrix(groupSums(dA[,rownames(tv)], 2, tv$Sector2)))
#dN <- colSums(as.matrix(groupSums(dN[,rownames(tv)], 2, tv$Sector2)))
dA <- colSums(as.matrix(groupSums(dA[,rownames(tv)], 2, tv$Sector)))
dN <- colSums(as.matrix(groupSums(dN[,rownames(tv)], 2, tv$Sector)))
U <- dN/dA
seffN <- cbind(dA=dA, dN=dN, U=U)[c("Agriculture","Forestry",
    "Energy",#"EnergySoftLin","MineWell",
    "RuralUrban","Transportation"),]

ThbScr <- colSums(hbScr[lxn$S,])
ThbSrf <- colSums(hbSrf[lxn$S,])
df <- (ThbScr - ThbSrf) / sum(ThbSrf)
dA <- Xtab(AsoilS ~ rf + cr, ch2soil)
dN <- Xtab(df ~ rf + cr, ch2soil)
#dA <- colSums(as.matrix(groupSums(dA[,rownames(ts)], 2, ts$Sector2)))
#dN <- colSums(as.matrix(groupSums(dN[,rownames(ts)], 2, ts$Sector2)))
dA <- colSums(as.matrix(groupSums(dA[,rownames(ts)], 2, ts$Sector)))
dN <- colSums(as.matrix(groupSums(dN[,rownames(ts)], 2, ts$Sector)))
U <- dN/dA
seffS <- cbind(dA=dA, dN=dN, U=U)[c("Agriculture","Forestry",
    "Energy",#"EnergySoftLin","MineWell",
    "RuralUrban","Transportation"),]
seff_res[[spp]] <- list(N=seffN, S=seffS)

#(sum(hbNcr)-sum(hbNrf))/sum(hbNrf)
#(sum(km$CurrN)-sum(km$RefN))/sum(km$RefN)
#100*seff

}

#save(slt, seff_res, file=file.path(ROOT, "out", "birds", "sector-effects-e2.Rdata"))
save(slt, seff_res, uplow, uplow_luf, uplow_full, file=file.path(ROOT, "out", "birds", "sector-effects.Rdata"))
#load(file.path(ROOT, "out", "birds", "sector-effects.Rdata"))

nres <- list()
sres <- list()
for (spp in names(seff_res)) {
    nres[[spp]] <- 100*c(PopEffect=seff_res[[spp]]$N[,2], UnitEffect=seff_res[[spp]]$N[,3])
    sres[[spp]] <- 100*c(PopEffect=seff_res[[spp]]$S[,2], UnitEffect=seff_res[[spp]]$S[,3])
}
nres <- do.call(rbind, nres)
sres <- do.call(rbind, sres)
nres <- data.frame(Species=slt[rownames(nres), "species"], nres)
sres <- data.frame(Species=slt[rownames(sres), "species"], sres)

uplow <- do.call(rbind, uplow)
uplow <- data.frame(Species=slt[rownames(uplow), "species"], uplow)

summary(uplow$Cur.total)
sum(uplow$Cur.total >= 10^7)
uplow2 <- uplow[uplow$Cur.total < 10^7,]
#sum(uplow$Cur.total <= 1)
#uplow2 <- uplow2[uplow$Cur.total > 0,]

write.csv(uplow2, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_UplandLowlandIntactness.csv")

write.csv(nres, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_SectorEffects_North.csv")
write.csv(sres, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_SectorEffects_South.csv")

siLUF <- data.frame(Species=NA, LUF=rownames(uplow_luf[[1]]), do.call(rbind, uplow_luf))
siLUF$Species <- slt$species[match(siLUF$ID, slt$AOU)]
ta <- aggregate(siLUF$Cur.total, list(ID=siLUF$ID), sum)
keep <- as.character(ta$ID[ta$x < 10^7 & ta$x > 0.5])
siLUF <- siLUF[siLUF$ID %in% keep,]
siLUF$ID <- NULL
write.csv(siLUF, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_UplandLowlandIntactness-by-LUF.csv")

allluf <- do.call(rbind, uplow_full)
allluf <- data.frame(slt[match(allluf$sppid, slt$AOU),c("species","scinam")], 
    LUF_NSR=interaction(allluf$LUF_NAME, allluf$NSRNAME, drop=TRUE, sep="_"), allluf)
write.csv(allluf, row.names=FALSE,
    file="e:/peter/sppweb-html-content/species/birds/Birds_UplandLowland-by-LUF-NSR.csv")


for (spp in as.character(slt$AOU[slt$map.pred])) {
cat(spp, "\n");flush.console()

for (WHERE in c("N", "S")) {
if (slt[spp, ifelse(WHERE=="N","veghf.north", "soilhf.south")]) {
SEFF <- seff_res[[spp]][[WHERE]]

## Sector effect plot from Dave
## Sectors to plot and their order
sectors <- c("Agriculture","Forestry",
    "Energy",#"EnergySoftLin","MineWell",
    "RuralUrban","Transportation")
## Names that will fit without overlap  
sector.names <- c("Agriculture","Forestry",
    "Energy",#"EnergySL","EnergyEX",
    "RuralUrban","Transport")
## The colours for each sector above
c1 <- c("tan3","palegreen4","indianred3",#"hotpink4",
    "skyblue3","slateblue2")  
total.effect <- 100 * SEFF[sectors,"dN"]
unit.effect <- 100 * SEFF[sectors,"U"]
## Max y-axis at 20%, 50% or 100% increments 
## (made to be symmetrical with y-min, except if y-max is >100
ymax <- ifelse(max(abs(unit.effect))<20,20,
    ifelse(max(abs(unit.effect))<50,50,round(max(abs(unit.effect))+50,-2)))  
ymin <- ifelse(ymax>50,min(-100,round(min(unit.effect)-50,-2)),-ymax)
## This is to leave enough space at the top of bars for the text giving the % population change
ymax <- max(ymax,max(unit.effect)+0.08*(max(unit.effect)-min(unit.effect,0))) 
## This is to leave enough space at the bottom of negative bars for the 
## text giving the % population change
ymin <- min(ymin,min(unit.effect)-0.08*(max(unit.effect,0)-min(unit.effect))) 

NAM <- as.character(tax[spp, "English_Name"])
TAG <- ""

png(paste0("e:/peter/sppweb-html-content/species/birds/sector-",
    ifelse(WHERE=="N", "north/", "south/"),
    as.character(tax[spp, "file"]), TAG, ".png"),
    width=600, height=600)

q <- barplot(unit.effect,
    width=100 * SEFF[sectors,"dA"],
    space=0,col=c1,border=c1,ylim=c(ymin,ymax),
    ylab="Unit effect (%)",xlab="Area (% of region)",
    xaxt="n",cex.lab=1.3,cex.axis=1.2,tcl=0.3,
    xlim=c(0,round(sum(100 * SEFF[,"dA"])+1,0)),
    bty="n",col.axis="grey40",col.lab="grey40",las=2)

rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray88",border="gray88")
x.at<-pretty(c(0,sum(100 * SEFF[,"dA"])))
axis(side=1,tck=1,at=x.at,lab=rep("",length(x.at)),col="grey95")
y.at<-pretty(c(ymin,ymax),n=6)
axis(side=2,tck=1,at=y.at,lab=rep("",length(y.at)),col="grey95")
q <- barplot(unit.effect,
    width=100 * SEFF[sectors,"dA"],
    space=0,col=c1,border=c1,ylim=c(ymin,ymax),
    ylab="Unit effect (%)",xlab="Area (% of region)",
    xaxt="n",cex.lab=1.3,cex.axis=1.2,tcl=0.3,
    xlim=c(0,round(sum(100 * SEFF[,"dA"])+1,0)),
    bty="n",col.axis="grey40",col.lab="grey40",las=2,add=TRUE)
box(bty="l",col="grey40")
mtext(side=1,line=2,at=x.at,x.at,col="grey40",cex=1.2)
axis(side=1,at=x.at,tcl=0.3,lab=rep("",length(x.at)),col="grey40",
    col.axis="grey40",cex.axis=1.2,las=1)
abline(h=0,lwd=2,col="grey40")
## Set the lines so that nearby labels don't overlap
mtext(side=1,at=q+c(0,0,-1,0,+1),sector.names,col=c1,cex=1.3,
    adj=0.5,line=c(0.1,0.1,1.1,0.1,1.1))  
## Just above positive bars, just below negative ones
y <- unit.effect+0.025*(ymax-ymin)*sign(unit.effect)  
## Make sure there is no y-axis overlap in % change labels of 
## sectors that are close together on x-axis
if (abs(y[3]-y[4])<0.05*(ymax-ymin))  
    y[3:4]<-mean(y[3:4])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[3:4])]   
## Make sure there is no y-axis overlap in % change labels of sectors 
## that are close together on x-axis
if (abs(y[4]-y[5])<0.05*(ymax-ymin))  
    y[4:5]<-mean(y[4:5])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[4:5])]   
#if (abs(y[5]-y[6])<0.05*(ymax-ymin))  
#    y[5:6]<-mean(y[5:6])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[5:6])]   
text(q,y,paste(ifelse(total.effect>0,"+",""),
    sprintf("%.1f",total.effect),"%",sep=""),col="darkblue",cex=1.4)
mtext(side=3,line=1,at=0,adj=0, 
    paste0(NAM, " - ", ifelse(WHERE=="N", "North", "South")),
    cex=1.4,col="grey40")
dev.off()
}
}
}

## CoV
results10km_list <- list()

SPP <- as.character(slt$AOU[slt$map.pred])
for (spp in SPP) {

    load(file.path(OUTDIRB, spp, paste0(regs[1], ".Rdata")))
    rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
    rownames(pxScrB) <- rownames(pxSrfB) <- names(Cells)[Cells == 1]
    pxNcr0 <- pxNcrB
    #pxNrf0 <- pxNrfB
    pxScr0 <- pxScrB
    #pxSrf0 <- pxSrfB
    for (i in 2:length(regs)) {
        cat(spp, regs[i], "\n");flush.console()
        load(file.path(OUTDIRB, spp, paste0(regs[i], ".Rdata")))
        rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
        rownames(pxScrB) <- rownames(pxSrfB) <- names(Cells)[Cells == 1]
        pxNcr0 <- rbind(pxNcr0, pxNcrB)
    #    pxNrf0 <- rbind(pxNrf0, pxNrfB)
        pxScr0 <- rbind(pxScr0, pxScrB)
    #    pxSrf0 <- rbind(pxSrf0, pxSrfB)
    }

    pxNcr <- pxNcr0[rownames(ks),]
    #pxNrf <- pxNrf0[rownames(ks),]
    pxScr <- pxScr0[rownames(ks),]
    pxScr[is.na(pxScr)] <- 0
    #pxSrf <- pxSrf0[rownames(ks),]
    for (k in 1:ncol(pxNcr)) {
        qN <- quantile(pxNcr[is.finite(pxNcr[,k]),k], q, na.rm=TRUE)
        pxNcr[pxNcr[,k] > qN,k] <- qN
        qS <- quantile(pxScr[is.finite(pxScr[,k]),k], q, na.rm=TRUE)
        pxScr[pxScr[,k] > qS,k] <- qS
    }

TYPE <- "C" # combo
if (!slt[spp, "veghf.north"])
    TYPE <- "S"
if (!slt[spp, "soilhf.south"])
    TYPE <- "N"

wS <- 1-ks$pAspen
if (TYPE == "S")
    wS[] <- 1
if (TYPE == "N")
    wS[] <- 0
wS[ks$useS] <- 1
wS[ks$useN] <- 0

    cr <- wS * pxScr + (1-wS) * pxNcr

    crveg <- groupMeans(cr, 1, ks$Row10_Col10, na.rm=TRUE)

#results10km_list[[as.character(slt[spp,"sppid"])]] <- crveg
#}
#xy10km <- ks[,c("POINT_X","POINT_Y","Row10_Col10")]
#save(xy10km, results10km_list, file="w:/species/birds-provincial-10x10km-summary.Rdata")

    crvegm <- rowMeans(crveg)
    crvegsd <- apply(crveg, 1, sd)
    #crvegm <- apply(crveg, 1, median)
    #crvegsd <- apply(crveg, 1, IQR)
    covC <- crvegsd / crvegm
    #covN[is.na(covN)] <- 1

    #crsoil <- groupMeans(pxScr, 1, ks$Row10_Col10)
    #crsoilm <- rowMeans(crsoil)
    #crsoilsd <- apply(crsoil, 1, sd)
    #crsoilm <- apply(crsoil, 1, median)
    #crsoilsd <- apply(crsoil, 1, IQR)
    #covS <- crsoilsd / crsoilm
    #covS[is.na(covS)] <- 1

#px <- crveg[order(crvegm),]
#matplot(crvegm[order(crvegm)], crveg, type="l", lty=1)


    NAM <- as.character(tax[spp, "English_Name"])

    zval <- as.integer(cut(covC, breaks=br))
    zval <- zval[match(kgrid$Row10_Col10, rownames(crveg))]

    cat(spp, "saving CoV map\n\n");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-cov-cr/",
        as.character(tax[spp, "file"]), ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=Col[zval], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "CoV"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)

    TEXT <- paste0(100*br[-length(br)], "-", 100*br[-1])
    INF <- grepl("Inf", TEXT)
    if (any(INF))
        TEXT[length(TEXT)] <- paste0(">", 100*br[length(br)-1])

    TITLE <- "Coefficient of variation"
    legend("bottomleft", border=rev(Col), fill=rev(Col), bty="n", legend=rev(TEXT),
                title=TITLE, cex=legcex*0.8)
    par(op)
    dev.off()


}



