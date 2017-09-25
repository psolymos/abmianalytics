library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"
#OUTDIR1 <- "e:/peter/AB_data_v2016/out/birds/pred1-josmshf"
#OUTDIRB <- "e:/peter/AB_data_v2016/out/birds/predB-josmshf"
STAGE <- list(veg = 6) # hab=5, hab+clim=6, hab+clim+shf=7

OUTDIR1 <- paste0("e:/peter/josm/2017/stage", STAGE$veg, "/pred1")
OUTDIRB <- paste0("e:/peter/josm/2017/stage", STAGE$veg, "/predB")


load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
#source("~/repos/bragging/R/glm_skeleton.R")
#source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/abmianalytics/R/maps_functions.R")
regs <- levels(kgrid$LUFxNSR)
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$useBCR6 <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"

e <- new.env()
#load(file.path(ROOT, "data", "data-full-withrevisit.Rdata"), envir=e)
load(file.path(ROOT, "out", "birds", "data", "data-wrsi.Rdata"), envir=e)
TAX <- droplevels(e$TAX)
TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")

en <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-josmshf.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY
BBn <- en$BB
tax <- droplevels(TAX[colnames(yyn),])

rm(e, en)

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
CSI <- colorRampPalette(c("red","yellow","green"), space = "rgb")(100)

q <- 0.99
H <- 1000
W <- 600


## csv

#spp <- "ALFL"
SPP <- rownames(tax)
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL")

PRED_DIR_IN <- "pred1-josmshf"
#PRED_DIR_IN <- "pred1"
#PRED_DIR_OUT <- "pred1cmb"

PREDS <- matrix(0, sum(kgrid$useBCR6), length(SPP))
rownames(PREDS) <- rownames(kgrid)[kgrid$useBCR6]
colnames(PREDS) <- SPP
PREDS0 <- PREDS

AREA_ha <- (1-kgrid$pWater) * kgrid$Area_km2 * 100
AREA_ha <- AREA_ha[kgrid$useBCR6]

for (spp in SPP) {

cat(spp, "--------------------------------------\n");flush.console()

load(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp, paste0(regs[1], ".Rdata")))
rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
pxNcr <- pxNcr1
pxNrf <- pxNrf1
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp, paste0(regs[i], ".Rdata")))
    rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
    pxNcr <- rbind(pxNcr, pxNcr1)
    pxNrf <- rbind(pxNrf, pxNrf1)
}

PREDS[,spp] <- pxNcr[rownames(PREDS),]
PREDS0[,spp] <- pxNrf[rownames(PREDS0),]

}

#save(PREDS, PREDS0, file=file.path(ROOT, "out", "birds", "josmshf", "predictions.Rdata"))
load(file.path(ROOT, "out", "birds", "josmshf", "predictions.Rdata"))

N <- colSums(PREDS*AREA_ha) / 10^6
#N <- N[N < max(N)]
summary(N)

## PIF table
pif <- read.csv("~/Dropbox/bam/PIF-AB/popBCR-6AB_v2_22-May-2013.csv")
mefa4::compare_sets(tax$English_Name, pif$Common_Name)
setdiff(tax$English_Name, pif$Common_Name)
pif <- pif[match(tax$English_Name, pif$Common_Name),]


## roadside_bias
load(file.path(ROOT, "out", "birds", "josmshf", "roadside_bias.Rdata"))

load(file.path(ROOT, "out", "birds", "data", "mean-qpad-estimates.Rdata"))
qpad_vals <- qpad_vals[rownames(tax),]

## roadside avoidance
library(mefa4)
load(file.path(ROOT, "out", "birds", "josmshf", "roadside_avoidance.Rdata"))
tmp <- cbind(ROAD=rai_data$ROAD, rai_pred)
rai <- groupSums(tmp[BBn[,1],], 1, rai_data$HAB[BBn[,1]], TRUE)
rai <- t(t(rai) / colSums(rai))
RAI <- 1 - colSums(rai[,1] * rai)
summary(RAI)
RAIc <- RAI-RAI["ROAD"]

#yy <- cbind(ALL=1, ROAD=xnn[BBn[,1],"ROAD01"],
#    ifelse(as.matrix(yyn[BBn[,1],]) > 0, 1, 0))
#rai <- groupSums(yy, 1, xnn$hab1[BBn[,1]], TRUE)
#n <- rai[,"ALL"]
#rai <- rai[,-1]
#rai <- t(t(rai) / colSums(rai))
#sai <- groupSums(yy, 1, xnn$hab1[BBn[,1]], TRUE)
#RAI <- 1 - colSums(rai[,1] * rai)

pop <- tax[,c("Species_ID", "English_Name", "Scientific_Name", "Spp")]
pop$RAI <- RAI[match(rownames(pop), names(RAI))]
pop$RAIc <- RAIc[match(rownames(pop), names(RAIc))]
pop$RAIroad <- RAI["ROAD"]
pop$Don <- roadside_bias[rownames(pop), "on"]
pop$Doff <- roadside_bias[rownames(pop), "off"]
pop$DeltaRoad <- roadside_bias[rownames(pop), "onoff"]
pop$Nqpad <- colSums(PREDS*AREA_ha) / 10^6 # M males
pop$Nqpad[pop$Nqpad > 1000] <- NA
pop$Npif <- (pif$Pop_Est / pif$Pair_Adjust) / 10^6 # M males
pop$DeltaObs <- pop$Nqpad / pop$Npif
pop$TimeAdj <- pif$Time_Adjust
pop$MDD <- pif$Detection_Distance_m
pop$p3 <- 1-exp(-3 * qpad_vals$phi0)
pop$EDR <- qpad_vals$phi0 * 100
pop$DeltaTime <- (1/pop$p3)/pop$TimeAdj
pop$DeltaDist <- pop$MDD^2 / pop$EDR^2
pop$DeltaExp <- pop$DeltaRoad * pop$DeltaTime * pop$DeltaDist
pop$DeltaRes <- pop$DeltaObs / pop$DeltaExp
pop <- pop[rowSums(is.na(pop))==0,]

#write.csv(pop, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/qpad-pif-results.csv")

boxplot(log(pop[,c("DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")]))
abline(h=0, col=2)

boxplot(log(pop[,c("DeltaObs", "DeltaExp")]))
abline(h=0, col=2)

mat <- log(pop[,c("DeltaObs", "DeltaExp", "DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")])
rnd <- runif(nrow(pop), -0.1, 0.1)
boxplot(mat, range=0)
for (i in 2:ncol(mat))
    segments(x0=i+rnd-1, x1=i+rnd, y0=mat[,i-1], y1=mat[,i], col="lightgrey")
for (i in 1:ncol(mat))
    points(i+rnd, mat[,i], col="darkgrey", pch=19)
abline(h=0, col=2, lwd=2)
boxplot(mat, range=0, add=TRUE)

with(pop, plot(RAI, log(DeltaRes), type="n"))
abline(h=0, v=RAI["ROAD"], col=2, lwd=2)
with(pop, text(RAI, log(DeltaRes), rownames(pop), cex=0.75))

boxplot(pop[,c("Npif", "Nqpad")], ylim=c(0,10))


## plots

res_luf <- list()
res_nsr <- list()
#SPP <- as.character(slt$AOU[slt$map.pred])
for (spp in SPP) {

    cat(spp, "\t");flush.console()
    load(file.path(ROOT, "out", "birds", "pred1cmb", paste0(spp, ".Rdata")))
    km <- data.frame(km)

    TYPE <- "C" # combo
    #if (!slt[spp, "veghf.north"])
    if (!(spp %in% fln))
        TYPE <- "S"
    #if (!slt[spp, "soilhf.south"])
    if (!(spp %in% fls))
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

    #km2 <- as.matrix(cbind(Curr=cr, Ref=rf))
    #rownames(km2) <- rownames(km)
    #save(km2, file=file.path(ROOT, "out", "birds", "pred1combined", paste0(spp, ".Rdata")))
    #cat("\n")
if (FALSE) {
    ndat <- normalize_data(rf=rf, cr=cr)
}
#    cr <- km$CurrN
#    rf <- km$RefN
#    cr <- km$CurrS
#    rf <- km$RefS
    qcr <- quantile(cr, q)
    cr[cr>qcr] <- qcr
    qrf <- quantile(rf, q)
    rf[rf>qrf] <- qrf

if (TRUE) {
    mat <- 100 * cbind(Ncurrent=cr, Nreference=rf) # ha to km^2
    rownames(mat) <- rownames(kgrid)
    res_luf[[spp]] <- groupSums(mat, 1, kgrid$LUF_NAME)
    res_nsr[[spp]] <- groupSums(mat, 1, kgrid$NSRNAME)
}

if (TRUE) {
    SI <- round(100 * pmin(cr, rf) / pmax(cr, rf))
    SI[is.na(SI)] <- 100 # 0/0 is defined as 100 intact
#    SI <- 100*as.matrix(dd1km_pred[[4]])[,"UNK"]/rowSums(dd1km_pred[[2]])
#    SI <- 100-SI
    cr0 <- cr
    rf0 <- rf
    SI0 <- SI
    SI[SI < 1] <- 1 # this is only for mapping

if (FALSE) {
library(raster)
source("~/repos/abmianalytics/R/maps_functions.R")
rt <- raster(file.path(ROOT, "data", "kgrid", "AHM1k.asc"))
r_si <- as_Raster(as.factor(kgrid$Row), as.factor(kgrid$Col), SI0, rt)
plot(r_si)
writeRaster(r_si, paste0(spp, "-intactness_2016-08-12.tif"), overwrite=TRUE)
}

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

    cat("si\t");flush.console()
    fname <- file.path(ROOT, "out", "birds", "figs", "map-si",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    png(fname, width=W, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1)
    plot(kgrid$X, kgrid$Y, col=CSI[SI], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, "\nIntactness"),col="grey30", cex=legcex)
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

if (FALSE) {
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
m0 <- as.matrix(dd1km_pred[[2]])
m0 <- 100*m0/rowSums(m0)
m0 <- m0[rf0==0 & m0[,"Water"] <= 99 & m0[,"NonVeg"] <= 99 & kgrid$NRNAME == "Grassland",]
m0 <- m0[,colSums(m0)>0]
#summary(m0)
round(colMeans(m0))

m0 <- as.matrix(dd1km_pred[[4]])
m0 <- 100*m0/rowSums(m0)
m0 <- m0[rf0==0 & m0[,"Water"] <= 99 & kgrid$NRNAME == "Grassland",]
m0 <- m0[,colSums(m0)>0]
round(colMeans(m0))

aggregate(100*as.matrix(dd1km_pred[[2]])[,"NonVeg"]/rowSums(dd1km_pred[[2]]),
    list(nr=kgrid$NRNAME, rf0=rf0==0 & kgrid$pWater < 0.9), mean)

}

    cat("rf\t");flush.console()
    fname <- file.path(ROOT, "out", "birds", "figs", "map-rf",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    png(fname, width=W, height=H)
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

    cat("cr\t");flush.console()
    fname <- file.path(ROOT, "out", "birds", "figs", "map-cr",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    png(fname, width=W, height=H)
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

    cat("df\n");flush.console()
    fname <- file.path(ROOT, "out", "birds", "figs", "map-df",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    png(fname, width=W, height=H)
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
}
}


#save(res_nsr, res_luf, file=file.path(ROOT, "out", "birds", "tables", "luf_Nsummaries.Rdata"))
load(file.path(ROOT, "out", "birds", "tables", "luf_Nsummaries.Rdata"))

LUF <- list()
for (spp in names(res_luf)) {
    tmp <- res_luf[[spp]] / 10^6 # M males
    tmp <- t(matrix(tmp, 2*nrow(tmp), 1))
    colnames(tmp) <- paste(rep(colnames(res_luf[[1]]), each=ncol(tmp)/2),
        rep(rownames(res_luf[[1]]), 2))
    LUF[[spp]] <- data.frame(Species=tax[spp, "English_Name"], tmp)
}
NSR <- list()
for (spp in names(res_nsr)) {
    tmp <- res_nsr[[spp]] / 10^6 # M males
    tmp <- t(matrix(tmp, 2*nrow(tmp), 1))
    colnames(tmp) <- paste(rep(colnames(res_nsr[[1]]), each=ncol(tmp)/2),
        rep(rownames(res_nsr[[1]]), 2))
    NSR[[spp]] <- data.frame(Species=tax[spp, "English_Name"], tmp)
}
LUF <- do.call(rbind, LUF)
NSR <- do.call(rbind, NSR)
write.csv(LUF, row.names=FALSE,
    file=file.path(ROOT, "out", "birds", "tables", "Birds_Abundance_by_LUFregions.csv"))
write.csv(NSR, row.names=FALSE,
    file=file.path(ROOT, "out", "birds", "tables", "Birds_Abundance_by_NaturalSubregions.csv"))


## sector effects
seff_res <- list()
tr_res <- list()
#seff_luf <- list()
#seff_ns <- list()
#uplow <- list()
#uplow_full <- list()
#uplow_luf <- list()

## stuff to exclude
## add col to lxn
## subset counter for loop

PRED_DIR_IN <- "pred1-shf" # "pred1-seismic-as-ES" # "pred1"

restrict_to_HF <- FALSE

TAX <- read.csv("~/repos/abmispecies/_data/birds.csv")
SPP <- as.character(TAX$AOU)[TAX$map.pred]

for (spp in SPP) {

cat(spp, "------------------------\n");flush.console()

#load(file.path(OUTDIR1, spp, paste0(regs[1], ".Rdata")))
load(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp, paste0(regs[1], ".Rdata")))
hbNcr <- hbNcr1[,1]
hbNrf <- hbNrf1[,1]
hbScr <- hbScr1[,1]
hbSrf <- hbSrf1[,1]
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    #load(file.path(OUTDIR1, spp, paste0(regs[i], ".Rdata")))
    load(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp, paste0(regs[i], ".Rdata")))
    hbNcr <- rbind(hbNcr, hbNcr1[,1])
    hbNrf <- rbind(hbNrf, hbNrf1[,1])
    hbScr <- rbind(hbScr, hbScr1[,1])
    hbSrf <- rbind(hbSrf, hbSrf1[,1])
}

if (!NSest["north"]) {
    hbNcr[] <- 0
    hbNrf[] <- 0
}
if (!NSest["south"]) {
    hbScr[] <- 0
    hbSrf[] <- 0
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
if (FALSE) {
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
}

## for HF-only sector effects, only need to adjust the total
## i.e. sum only where ch2veg$cr is HF
keep <- rownames(ch2veg) %in% rownames(ch2veg)[ch2veg$isHF]
hbNcr_HFonly <- hbNcr
hbNrf_HFonly <- hbNrf
hbNcr_HFonly[,!keep] <- 0
hbNrf_HFonly[,!keep] <- 0
ThbNcr_HFonly <- colSums(hbNcr_HFonly[lxn$N,])
ThbNrf_HFonly <- colSums(hbNrf_HFonly[lxn$N,])
ThbNcr <- colSums(hbNcr[lxn$N,])
ThbNrf <- colSums(hbNrf[lxn$N,])
Ntot_HFonly <- sum(ThbNrf_HFonly)
Ntot_All <- sum(ThbNrf)
Ntot_Use <- if (restrict_to_HF)
    Ntot_HFonly else Ntot_All
df <- (ThbNcr - ThbNrf) / Ntot_Use
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

keep <- rownames(ch2soil) %in% rownames(ch2soil)[ch2soil$isHF]
hbScr_HFonly <- hbScr
hbSrf_HFonly <- hbSrf
hbScr_HFonly[,!keep] <- 0
hbSrf_HFonly[,!keep] <- 0
ThbScr_HFonly <- colSums(hbScr_HFonly[lxn$S,])
ThbSrf_HFonly <- colSums(hbSrf_HFonly[lxn$S,])
ThbScr <- colSums(hbScr[lxn$S,])
ThbSrf <- colSums(hbSrf[lxn$S,])
Stot_HFonly <- sum(ThbSrf_HFonly)
Stot_All <- sum(ThbSrf)
Stot_Use <- if (restrict_to_HF)
    Stot_HFonly else Stot_All
df <- (ThbScr - ThbSrf) / Stot_Use
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
tr_res[[spp]] <- list(N=cbind(rf=ThbNrf, cr=ThbNcr), S=cbind(rf=ThbSrf, cr=ThbScr),
    NSest=NSest, total=c(Ntot_All=Ntot_All, Stot_All=Stot_All,
    Ntot_HFonly=Ntot_HFonly, Stot_HFonly=Stot_HFonly))
#(sum(hbNcr)-sum(hbNrf))/sum(hbNrf)
#(sum(km$CurrN)-sum(km$RefN))/sum(km$RefN)
#100*seff

}

## -new version has the HFonly pop sizes saved
## can be used to retro-fit the effects
#save(seff_res, tr_res, file=file.path(ROOT, "out", "birds", "tables", "sector-effects-new-seismic-as-bf.Rdata"))
#save(seff_res, tr_res, file=file.path(ROOT, "out", "birds", "tables", "sector-effects-new-seismic-as-ES.Rdata"))
#save(seff_res, tr_res, file=file.path(ROOT, "out", "birds", "tables", "sector-effects-new-shf.Rdata"))
load(file.path(ROOT, "out", "birds", "tables", "sector-effects-new-seismic-as-bf.Rdata"))
#load(file.path(ROOT, "out", "birds", "tables", "sector-effects-new-seismic-as-ES.Rdata"))
#load(file.path(ROOT, "out", "birds", "tables", "sector-effects-new-shf.Rdata"))

#spp <- "ALFL"
seff_loc <- list()
seff_lfull <- list()
for (spp in names(tr_res)) {
    seff_lfull[[spp]] <- groupSums(as.matrix(tr_res[[spp]]$N)[ch2veg$isHF,], 1,
        as.character(ch2veg$cr[ch2veg$isHF]))
    seff_loc[[spp]] <- groupSums(seff_lfull[[spp]], 1, tv[rownames(seff_lfull[[spp]]), "Sector"])
}
seff2 <- t(sapply(seff_loc, function(z) (z[,"cr"]-z[,"rf"])/z[,"rf"]))
seff2 <- seff2[,rownames(seff_res[[1]]$N)]
seff1 <- t(sapply(seff_res, function(z) z$N[,"dN"]))
seff2 <- cbind(seff2,
    Total=sapply(seff_loc, function(z) (sum(z[,"cr"])-sum(z[,"rf"]))/sum(z[,"rf"])))

seff2 <- seff2[!is.na(seff2[,1]),]
seff2 <- seff2[order(rownames(seff2)),]
seff2[seff2>2] <- 2
round(100*seff2,1)
#AA <- 100*seff_res[[1]]$N[,"dA"]
#a <- round(100*seff2,1)
#aa <- round(t(t(100*seff2) / AA), 1)

d1 <- apply(100*seff1, 2, density, na.rm=TRUE)
d2 <- apply(100*seff2, 2, density, na.rm=TRUE)
for (i in 1:5) {
    d1[[i]]$y <- d1[[i]]$y/max(d1[[i]]$y)
    d2[[i]]$y <- d2[[i]]$y/max(d2[[i]]$y)
}
par(mfrow=c(1,2))
plot(d1[[1]], xlim=c(-50,50), main="% change inside region", lwd=2)
for (i in 2:5) lines(d1[[i]], col=i, lwd=2)
abline(v=0)
plot(d2[[1]], xlim=c(-100,200), main="% change inside HF", lwd=2)
for (i in 2:5) lines(d1[[i]], col=i, lwd=2)
abline(v=0)
legend("topright", lty=1, col=1:5, bty="n", legend=colnames(seff2), lwd=2)

par(mfrow=c(2,3))
for (i in 1:6) {
    hist(100*seff2[,i], main=colnames(seff2)[i], col="lightblue",
        xlab="% population change inside HF", border=NA)
    abline(v=0, col=4, lty=2)
}
write.csv(round(100*seff2,1), file="sector-effects-birds-early-seral-seismic.csv")

nres <- list()
sres <- list()
for (spp in names(seff_res)) {
    nres[[spp]] <- 100*c(PopEffect=seff_res[[spp]]$N[,2], UnitEffect=seff_res[[spp]]$N[,3])
    sres[[spp]] <- 100*c(PopEffect=seff_res[[spp]]$S[,2], UnitEffect=seff_res[[spp]]$S[,3])
}
nres <- do.call(rbind, nres)
sres <- do.call(rbind, sres)
nres <- data.frame(Species=tax[rownames(nres), "English_Name"], nres)
sres <- data.frame(Species=tax[rownames(sres), "English_Name"], sres)

## keep only spp that are OK

write.csv(nres, row.names=FALSE,
    file=file.path(ROOT, "out", "birds", "tables", "Birds_SectorEffects_North.csv"))
write.csv(sres, row.names=FALSE,
    file=file.path(ROOT, "out", "birds", "tables", "Birds_SectorEffects_South.csv"))

for (spp in SPP) {
cat(spp, "\n");flush.console()

for (WHERE in c("north", "south")) {
SEFF <- seff_res[[spp]][[ifelse(WHERE=="north", "N", "S")]]

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

TOTALS <- if (WHERE=="north")
    tr_res[[spp]]$total[c("Ntot_All", "Ntot_HFonly")] else tr_res[[spp]]$total[c("Stot_All", "Stot_HFonly")]
SCALING <- if (restrict_to_HF)
    TOTALS[1]/TOTALS[2] else 1
total.effect <- 100 * SCALING * SEFF[sectors,"dN"]
#unit.effect <- 100 * SEFF[sectors,"U"]
unit.effect <- 100 * SCALING * SEFF[sectors,"dN"] / SEFF[sectors,"dA"]
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

png(file.path(ROOT, "out", "birds", "figs",
    paste0("sector-", if (restrict_to_HF) "HFonly-" else "", WHERE),
    paste0(as.character(tax[spp, "Spp"]), TAG, ".png")),
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
    paste0(NAM, " - ", ifelse(WHERE=="north", "North", "South")),
    cex=1.4,col="grey40")
dev.off()
}
}

## CoV

results10km_list <- list()
for (spp in SPP) {

    load(file.path(ROOT, "out", "birds", "predB", spp, paste0(regs[1], ".Rdata")))
    rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
    rownames(pxScrB) <- rownames(pxSrfB) <- names(Cells)[Cells == 1]
    pxNcr0 <- pxNcrB
    #pxNrf0 <- pxNrfB
    pxScr0 <- pxScrB
    #pxSrf0 <- pxSrfB
    for (i in 2:length(regs)) {
        cat(spp, regs[i], "\n");flush.console()
        load(file.path(ROOT, "out", "birds", "predB", spp, paste0(regs[i], ".Rdata")))
        rownames(pxNcrB) <- rownames(pxNrfB) <- names(Cells)[Cells == 1]
        rownames(pxScrB) <- rownames(pxSrfB) <- names(Cells)[Cells == 1]
        pxNcr0 <- rbind(pxNcr0, pxNcrB)
    #    pxNrf0 <- rbind(pxNrf0, pxNrfB)
        pxScr0 <- rbind(pxScr0, pxScrB)
    #    pxSrf0 <- rbind(pxSrf0, pxSrfB)
    }

    pxNcr <- pxNcr0[rownames(ks),]
    pxNcr[is.na(pxNcr)] <- 0
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

    TR <- FALSE # transform to prob scale
    TYPE <- "C" # combo
    #if (!slt[spp, "veghf.north"])
    if (!(spp %in% fln))
        TYPE <- "S"
    #if (!slt[spp, "soilhf.south"])
    if (!(spp %in% fls))
        TYPE <- "N"

wS <- 1-ks$pAspen
if (TYPE == "S")
    wS[] <- 1
if (TYPE == "N")
    wS[] <- 0
wS[ks$useS] <- 1
wS[ks$useN] <- 0

    cr <- wS * pxScr + (1-wS) * pxNcr
    cr <- 100*cr
#    if (TR)
#        cr <- 1-exp(-cr)

    crveg <- groupMeans(cr, 1, ks$Row10_Col10, na.rm=TRUE)

results10km_list[[as.character(tax[spp,"Spp"])]] <- crveg
}
xy10km <- ks[,c("POINT_X","POINT_Y","Row10_Col10")]
save(xy10km, results10km_list, file=file.path(ROOT, "out", "birds", "tables", "km10results.Rdata"))

slt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(slt) <- slt$AOU
slt$comments <- NULL

mcrvegsd <- matrix(0, nrow(results10km_list[[1]]), sum(slt$map.pred))
rownames(mcrvegsd) <- rownames(results10km_list[[1]])
colnames(mcrvegsd) <- as.character(tax[rownames(slt)[slt$map.pred],"Spp"])
for (spp in rownames(slt)[slt$map.pred]) {
    crveg <- results10km_list[[as.character(tax[spp,"Spp"])]]
    mcrvegsd[,as.character(tax[spp,"Spp"])] <- apply(crveg, 1, sd)
}
write.csv(mcrvegsd, file=file.path(ROOT, "out", "birds", "tables", "birds-10x10km-SD-summary.csv"))
write.csv(xy10km, file=file.path(ROOT, "out", "birds", "tables", "xy-for-10x10km-SD-summary.csv"))

TAG <- ""
for (spp in SPP) {
    crveg <- results10km_list[[as.character(tax[spp,"Spp"])]]
    crvegm <- rowMeans(crveg)
    crvegsd <- apply(crveg, 1, sd)
    crvegsd[crvegsd==0] <- 0.000001
    #crvegm <- apply(crveg, 1, median)
    #crvegsd <- apply(crveg, 1, IQR)
    covC <- crvegsd / crvegm
    covC[crvegm==0] <- mean(covC[crvegm!=0], na.rm=TRUE) # will not stick out...
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
    cat(spp)

if (FALSE) {
    cat("\tMean");flush.console()
    zval <- as.integer(cut(covC, breaks=br))
    zval <- pmin(100, ceiling(99 * (crvegm / max(crvegm, na.rm=TRUE)))+1)
    zval <- zval[match(kgrid$Row10_Col10, rownames(crveg))]
    fname <- file.path(ROOT, "out", "birds", "figs", "map-test",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    png(fname, width=W*3, height=H)
    op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,3))
    plot(kgrid$X, kgrid$Y, col=C1[zval], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
#    with(kgrid[kgrid$NRNAME == "Rocky Mountain" & kgrid$POINT_X < -112,],
#        points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "N")
        with(kgrid[kgrid$useS,], points(X, Y, col=CE, pch=15, cex=cex))
    if (TYPE == "S")
        with(kgrid[kgrid$useN,], points(X, Y, col=CE, pch=15, cex=cex))
    mtext(side=3,paste(NAM, TAG, "Mean"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)
    #par(op)
    #dev.off()
}

    cat("\tCoV");flush.console()
    zval <- as.integer(cut(covC, breaks=br))
    zval <- zval[match(kgrid$Row10_Col10, rownames(crveg))]
    fname <- file.path(ROOT, "out", "birds", "figs", "map-cov-cr",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    if (TR)
        fname <- file.path(ROOT, "out", "birds", "figs", "map-test",
            paste0(as.character(tax[spp, "Spp"]), TAG, "-CoV", ".png"))
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
    mtext(side=3,paste(NAM, TAG, "CoV"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)

    TEXT <- paste0(100*br[-length(br)], "-", 100*br[-1])
    INF <- grepl("Inf", TEXT)
    if (any(INF))
        TEXT[length(TEXT)] <- paste0(">", 100*br[length(br)-1])

    TITLE <- "Coefficient of variation"
    legend("bottomleft", border=rev(Col), fill=rev(Col), bty="n", legend=rev(TEXT),
                title=TITLE, cex=legcex*0.8)
    par(op)
    dev.off()

    cat("\tSD\n");flush.console()
    zval <- as.integer(cut(crvegsd/mean(crvegm,na.rm=TRUE), breaks=br))
    zval <- zval[match(kgrid$Row10_Col10, rownames(crveg))]
    fname <- file.path(ROOT, "out", "birds", "figs", "map-sd-cr",
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png"))
    if (TR)
        fname <- file.path(ROOT, "out", "birds", "figs", "map-test",
            paste0(as.character(tax[spp, "Spp"]), TAG, "-SD", ".png"))
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
    mtext(side=3,paste(NAM, TAG, "SE"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
#	text(378826,5774802,"Insufficient \n   data",col="white",cex=0.9)

    br2 <- round(br * mean(crvegm,na.rm=TRUE), 3)
    TEXT <- paste0(br2[-length(br2)], "-", br2[-1])
    INF <- grepl("Inf", TEXT)
    if (any(INF))
        TEXT[length(TEXT)] <- paste0(">", br2[length(br2)-1])

    TITLE <- paste0("Prediction Std. Error\n(mean = ", round(mean(crvegm,na.rm=TRUE), 3), ")")
    legend("bottomleft", border=rev(Col), fill=rev(Col), bty="n", legend=rev(TEXT),
                title=TITLE, cex=legcex*0.8)
    par(op)
    dev.off()

}



