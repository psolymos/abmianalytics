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


## read in predictions and calculate guild intactness and richness

tax2 <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(tax2) <- tax2$AOU
tax3 <- read.csv("~/repos/abmianalytics/lookup/vertebrate-guilds.csv")
setdiff(tax2$AOU[tax2$map.pred], tax3$AOU.Code)
setdiff(tax2$AOU[tax2$map.pred], tax$Species_ID)
SPP <- intersect(tax2$AOU[tax2$map.pred], tax3$AOU.Code)

res_pocc <- list()
res_poccN <- list()
res_si <- list()
res_pmax <- list()
for (spp in SPP) {

    cat(spp, "\n");flush.console()
    km <- read.csv(paste0("e:/peter/sppweb2015/birds-pred/", 
        paste0(as.character(tax[spp, "file"]), ".csv")))

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
#    wS[kgrid$useS] <- 1
#    wS[kgrid$useN] <- 0

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

    res_pocc[[spp]] <- 1 - exp(-cr * (pi*1.5^2)) # 150m
    #res_poccN[[spp]] <- 1 - exp(-km$CurrN * (pi*1.5^2)) # 150m
    res_si[[spp]] <- 100 * pmin(cr, rf) / pmax(cr, rf)
    res_pmax[[spp]] <- max(cr, rf, na.rm=TRUE)
    ## NA means no prediction (rf=0 & cr=0)
}

#save(res_si, res_pocc, res_pmax, SPP, tax, tax2, tax3,
#    file=file.path(ROOT, "out", "birds", "si-and-richness-summaries.Rdata"))

SR <- do.call(cbind, res_pocc)
SI <- do.call(cbind, res_si)
#rownames(SR) <- rownames(SI) <- km$LinkID
SR <- SR[,SPP]
SI <- SI[,SPP]
tax2 <- tax2[SPP,]

meanSI_allbirds <- rowMeans(SI, na.rm=TRUE)
meanSI_oldforestbirds <- rowMeans(SI[,tax2$oldforest==1], na.rm=TRUE)
SR_allbirds <- rowSums(SR)
SR_oldforestbirds <- rowSums(SR[,tax2$oldforest==1])

out <- data.frame(
    SI_allbirds=meanSI_allbirds,
    SI_oldforestbirds=meanSI_oldforestbirds,
    SR_allbirds=SR_allbirds,
    SR_oldforestbirds=SR_oldforestbirds)
rownames(out) <- km$LinkID
save(out, file=file.path(ROOT, "out", "birds", "si-and-richness-summaries.Rdata"))

    colSI <- colorRampPalette(c("red","yellow","darkgreen"))(100)
    colSR <- colorRampPalette(c("blue","red"))(100)

    val <- SR_allbirds
    #val <- pmin(100, ceiling(99 * val / max(val,na.rm=TRUE))+1)
    val <- pmin(100, ceiling(99*val/max(val,na.rm=TRUE))+1)
    png(paste0("e:/peter/sppweb-html-content/birds-richness-all.png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=colSR[val], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,"Richness, Birds",col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    val <- SR_oldforestbirds
    #val <- pmin(100, ceiling(99 * val / max(val,na.rm=TRUE))+1)
    val <- pmin(100, ceiling(99*val/max(val,na.rm=TRUE))+1)
    png(paste0("e:/peter/sppweb-html-content/birds-richness-oldforest.png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=colSR[val], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,"Richness, Old Forest Birds",col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    val <- meanSI_oldforestbirds
    #val <- pmin(100, ceiling(99 * val / max(val,na.rm=TRUE))+1)
    val <- pmin(100, ceiling(0.99*val)+1)
    png(paste0("e:/peter/sppweb-html-content/birds-intactness-oldforest.png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=colSI[val], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,"Intactness, Old Forest Birds",col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()

    val <- meanSI_allbirds
    #val <- pmin(100, ceiling(99 * val / max(val,na.rm=TRUE))+1)
    val <- pmin(100, ceiling(0.99*val)+1)
    png(paste0("e:/peter/sppweb-html-content/birds-intactness-all.png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=colSI[val], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,"Intactness, Birds",col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    par(op)
    dev.off()
