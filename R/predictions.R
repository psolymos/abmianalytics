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

spp <- "ALFL"

for (spp in union(fln, fls)) {

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
for (spp in as.character(slt$AOU[slt$map.pred])) {
#    km <- read.csv(paste0("e:/peter/sppweb2015/birds-pred/", 
#        paste0(as.character(tax[spp, "file"]), ".csv")))

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

if (FALSE) {
wS <- 1-kgrid$pAspen
if (!slt[spp, "veghf.north"])
    wS[] <- 1
if (!slt[spp, "soilhf.south"])
    wS[] <- 0
wS[useS] <- 1
wS[useN] <- 1

pxcr <- wS * km$CurrNpxScr + (1-wS) * pxNcr
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

seff_res <- list()
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
ThbNcr <- colSums(hbNcr[lxn$N,])
ThbNrf <- colSums(hbNrf[lxn$N,])
df <- (ThbNcr - ThbNrf) / sum(ThbNrf)
dA <- Xtab(AvegN ~ rf + cr, ch2veg)
dN <- Xtab(df ~ rf + cr, ch2veg)
dA <- colSums(as.matrix(groupSums(dA[,rownames(tv)], 2, tv$Sector)))
dN <- colSums(as.matrix(groupSums(dN[,rownames(tv)], 2, tv$Sector)))
U <- dN/dA
seffN <- cbind(dA=dA, dN=dN, U=U)[-1,]

dimnames(hbScr) <- dimnames(hbSrf) <- list(regs, colnames(Asoil))
hbScr[is.na(hbScr)] <- 0
hbSrf[is.na(hbSrf)] <- 0
hbScr <- hbScr * Asoil
hbSrf <- hbSrf * Asoil
ThbScr <- colSums(hbScr[lxn$S,])
ThbSrf <- colSums(hbSrf[lxn$S,])
df <- (ThbScr - ThbSrf) / sum(ThbSrf)
dA <- Xtab(AsoilS ~ rf + cr, ch2soil)
dN <- Xtab(df ~ rf + cr, ch2soil)
dA <- colSums(as.matrix(groupSums(dA[,rownames(ts)], 2, ts$Sector)))
dN <- colSums(as.matrix(groupSums(dN[,rownames(ts)], 2, ts$Sector)))
U <- dN/dA
seffS <- cbind(dA=dA, dN=dN, U=U)[-1,]
seff_res[[spp]] <- list(N=seffN, S=seffS)

#(sum(hbNcr)-sum(hbNrf))/sum(hbNrf)
#(sum(km$CurrN)-sum(km$RefN))/sum(km$RefN)
#100*seff

}

save(slt, seff_res, file=file.path(ROOT, "out", "birds", "kgrid_table.Rdata"))

for (spp in as.character(slt$AOU[slt$map.pred])) {
cat(spp, "\n");flush.console()

for (WHERE in c("N", "S")) {
if (slt[spp, ifelse(WHERE=="N","veghf.north", "soilhf.south")]) {
SEFF <- seff_res[[spp]][[WHERE]]

## Sector effect plot from Dave
## Sectors to plot and their order
sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
## Names that will fit without overlap  
sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")  
## The colours for each sector above
c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")  
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

ks <- kgrid[kgrid$Rnd10 <= 10,]
xy10 <- groupMeans(as.matrix(kgrid[,c("X","Y")]), 1, kgrid$Row10_Col10)

load(file.path(OUTDIRB, spp, paste0(regs[1], ".Rdata")))
rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
rownames(pxScrB) <- rownames(pxSrfB) <- names(Cells)[Cells == 1]
pxNcr <- pxNcrB
pxNrf <- pxNrfB
pxScr <- pxScrB
pxSrf <- pxSrfB
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(OUTDIRB, spp, paste0(regs[i], ".Rdata")))
    rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
    rownames(pxScrB) <- rownames(pxSrfB) <- names(Cells)[Cells == 1]
    pxNcr <- rbind(pxNcr, pxNcrB)
#    pxNrf <- rbind(pxNrf, pxNrfB)
    pxScr <- rbind(pxScr, pxScrB)
#    pxSrf <- rbind(pxSrf, pxSrfB)
}

pxNcr <- pxNcr[rownames(ks),]
#pxNrf <- pxNrf[rownames(ks),]
pxScr <- pxScr[rownames(ks),]
#pxSrf <- pxSrf[rownames(ks),]

crveg <- groupMeans(pxNcr, 1, ks$Row10_Col10)
crvegm <- rowMeans(crveg)
crvegsd <- apply(crveg, 1, sd)
covN <- crvegsd / crvegm
covN[is.na(covN)] <- 1

crsoil <- groupMeans(pxScr, 1, ks$Row10_Col10)
crsoilm <- rowMeans(crsoil)
crsoilsd <- apply(crsoil, 1, sd)
covS <- crsoilsd / crsoilm
covS[is.na(covS)] <- 1



library(RColorBrewer)
br <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col <- rev(brewer.pal(10, "RdYlGn"))

NAM <- as.character(tax[spp, "English_Name"])

    zval <- as.integer(cut(covN, breaks=br))
    zval <- zval[match(kgrid$Row10_Col10, rownames(crsoil))]

    cat("df\n");flush.console()
    png(paste0("e:/peter/sppweb-html-content/species/birds/map-cov-cr/",
        as.character(tax[spp, "file"]), ".png"),
        width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=Col[zval], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,], 
        points(X, Y, col=CE, pch=15, cex=cex))
    with(kgrid[kgrid$NRNAME == "Grassland",], 
        points(X, Y, col=CE, pch=15, cex=cex))
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
                title=TITLE, cex=legcex)
    par(op)
    dev.off()



png(paste0("e:/peter/sppweb-html-content/species/birds/map-cov-cr/",
    as.character(tax[spp, "file"]), ".png"),
    width=W, height=H)



    op <- par(mar=c(0, 0, 1, 0) + 0.1)
    #with(kgrid, plot(X, Y, col="grey", pch=15, cex=0.35, ann=FALSE, axes=FALSE))
    with(kgrid, plot(X, Y, pch=15, cex=4, col=Col[zval]))
    title(main=NAM, cex.main=2)
    TEXT <- paste0(br[-length(br)], "-", br[-1])
    INF <- grepl("Inf", TEXT)
    if (any(INF))
        TEXT[length(TEXT)] <- paste0(">", br[length(br)-1])
    TITLE <- "Coefficient of variation"
    legend("bottomleft", border=rev(Col), fill=rev(Col), bty="n", legend=rev(TEXT),
                title=TITLE, cex=legcex)
    par(op)
    dev.off()

    invisible(NULL)
}



