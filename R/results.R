library(mefa4)
library(pbapply)
library(RColorBrewer)

ROOT <- "c:/p/AB_data_v2015/out/birds"

level <- 0.9

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

e <- new.env()
load(file.path(ROOT, "data", "data-full-withrevisit.Rdata"), envir=e)
dat <- e$DAT
dat <- dat[dat$useOK,]
yy <- e$YY[rownames(dat),]
tax <- droplevels(e$TAX[colnames(yy),])
pveghf <- e$pveghf[rownames(dat),]
pveghf <- data.frame(as.matrix(pveghf))
pveghf$Open <- pveghf$GrassHerb + pveghf$Shrub
pveghf <- as.matrix(pveghf[,c("Decid", "Mixwood", "Conif", "Pine", "BSpr", "Larch", 
    "Open", "Wetland", "Cult", "UrbInd", "HardLin", "SoftLin")])
colnames(pveghf) <- c("Deciduous", "Mixedwood", "White Spruce", "Pine", 
    "Black Spruce", "Larch", 
    "Open", "Wet", "Cultivated", "Urban/Industrial", "Hard Linear", "Soft Linear")
psoilhf <- as.matrix(e$psoilhf[rownames(dat),c("Productive", "Clay", 
    "Saline", "RapidDrain", "Cult", "UrbInd")])
colnames(psoilhf) <- c("Productive", "Clay", 
    "Saline", "Rapid Drain", "Cultivated", "Urban/Industrial")

en <- new.env()
load(file.path(ROOT, "data", "data-useok-north.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "data", "data-useok-south.Rdata"), envir=es)
xns <- es$DAT
modss <- es$mods
yys <- es$YY
rm(e, en, es)

## model for species
fl <- list.files(file.path(ROOT, "results"))
fln <- fl[grep("-north_", fl)]
fln <- sub("birds_bam-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- fl[grep("-south_", fl)]
fls <- sub("birds_bam-south_", "", fls)
fls <- sub(".Rdata", "", fls)

tax$ndet <- colSums(yy>0)
tax$modelN <- rownames(tax) %in% fln
tax$modelS <- rownames(tax) %in% fls
tax$ndet_n <- colSums(yyn>0)[match(colnames(yy), colnames(yyn))]
tax$ndet_s <- colSums(yys>0)[match(colnames(yy), colnames(yys))]
tax$ndet_n[is.na(tax$ndet_n)] <- 0
tax$ndet_s[is.na(tax$ndet_s)] <- 0

yy <- yy[,tax$ndet > 0]
tax <- droplevels(tax[colnames(yy),])
tax$file <- nameAlnum(as.character(tax$English_Name), "mixed", "")

pveghf <- pveghf[rownames(yyn),]
psoilhf <- psoilhf[rownames(yys),]

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

stage_hab_n <- 5
stage_hab_s <- 2

## tax placeholders for all the output
tax$map_det <- tax$ndet > 0
tax$useavail_north <- tax$ndet_n > 3
tax$useavail_south <- tax$ndet_s > 3
tax$trend_north <- tax$modelN
tax$trend_south <- tax$modelS
tax$veghf_north <- tax$modelN
tax$soilhf_nontreed_south <- tax$modelS
tax$soilhf_treed_south <- tax$modelS
tax$linear_north <- tax$modelN
tax$linear_south <- tax$modelS
tax$surroundinghf_north <- tax$modelN
tax$surroundinghf_south <- tax$modelS

yy <- yy[,colnames(yy) != "NONE"]
yys <- yys[,colnames(yys) != "NONE"]
yyn <- yyn[,colnames(yyn) != "NONE"]
tax <- droplevels(tax[colnames(yy),])

## spp specific output

spp <- "BTNW"

## useavail-north
## table: useavail-north

res_useavail_north <- list()
for (spp in rownames(tax)) {
if (tax[spp, "useavail_north"]) {
    cat(spp, "\n");flush.console()
    keep <- rowSums(pveghf) > 0
    yyy <- yyn[keep, spp]
    hhh <- pveghf[keep,]
    NAM <- as.character(tax[spp, "English_Name"])
    NDAT <- sum(yyn[,spp] > 0)
    fname <- file.path(ROOT, "figs", "useavail-north", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    png(file=fname, width=480, height=480)
    tmp <- plot_wrsi(yyy, hhh, south=FALSE)
    mtext(paste0(NAM, " (n = ", NDAT, " detections)"), adj=0, line=2,
        side=3, cex=1.2, col="grey40", las=1)
    dev.off()
    res_useavail_north[[spp]] <- tmp
}
}
wrsi_n <- t(sapply(res_useavail_north, "[[", "WRSI"))
colnames(wrsi_n) <- paste0(nameAlnum(rownames(res_useavail_north[[1]]), "mixed", ""),
    "_WRSI")
rwrsi_n <- t(sapply(res_useavail_north, "[[", "rWRSI"))
colnames(rwrsi_n) <- paste0(nameAlnum(rownames(res_useavail_north[[1]]), "mixed", ""),
    "_rWRSI")
wrsi_n <- data.frame(tax[rownames(wrsi_n), c("English_Name","Scientific_Name")],
    cbind(wrsi_n, rwrsi_n))
write.csv(wrsi_n, file=file.path(ROOT, "figs", "useavail-north.csv"))

## useavail-south
## table: useavail-south

res_useavail_south <- list()
for (spp in rownames(tax)) {
if (tax[spp, "useavail_south"]) {
    cat(spp, "\n");flush.console()
    keep <- rowSums(psoilhf) > 0
    yyy <- yys[keep, spp]
    hhh <- psoilhf[keep,]
    NAM <- as.character(tax[spp, "English_Name"])
    NDAT <- sum(yys[,spp] > 0)
    fname <- file.path(ROOT, "figs", "useavail-south", 
        paste0(as.character(tax[spp, "file"]), ".png"))
	png(file=fname, width=480, height=480)
    plot_wrsi(yyy, hhh, south=TRUE)
    mtext(paste0(NAM, " (n = ", NDAT, " detections)"), adj=0, line=2,
        side=3, cex=1.2, col="grey40", las=1)
    dev.off()
    res_useavail_south[[spp]] <- tmp
}
}
wrsi_s <- t(sapply(res_useavail_south, "[[", "WRSI"))
colnames(wrsi_s) <- paste0(nameAlnum(rownames(res_useavail_south[[1]]), "mixed", ""),
    "_WRSI")
rwrsi_s <- t(sapply(res_useavail_south, "[[", "rWRSI"))
colnames(rwrsi_s) <- paste0(nameAlnum(rownames(res_useavail_south[[1]]), "mixed", ""),
    "_WRSI")
wrsi_s <- data.frame(tax[rownames(wrsi_s), c("English_Name","Scientific_Name")],
    cbind(wrsi_s, rwrsi_s))
write.csv(wrsi_s, file=file.path(ROOT, "figs", "useavail-south.csv"))

## map-det

load(file.path("c:/p/AB_data_v2015/out", "kgrid", "kgrid_table.Rdata"))
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
for (spp in rownames(tax)) {
if (tax[spp, "map_det"]) {
    cat(spp, "\n");flush.console()
    xy0 <- as.matrix(dat[yy[,spp] == 0,c("X","Y")])
    xy1 <- as.matrix(dat[yy[,spp] > 0,c("X","Y")])
    NAM <- as.character(tax[spp, "English_Name"])
    NDAT <- sum(yy[,spp] > 0)
    fname <- file.path(ROOT, "figs", "map-det", 
        paste0(as.character(tax[spp, "file"]), ".png"))
	png(file=fname, width=600, height=1000)

    plot(kgrid$X, kgrid$Y, pch=15, cex=0.2, col=col1, axes=FALSE, ann=FALSE)
    points(xyw, pch=15, cex=0.2, col=rgb(0.3,0.45,0.9))
    points(xy0, pch="+", cex=0.5, col="red3")
    points(xy1, pch=16, cex=1.6, col="red4")
    mtext(paste0(NAM, " (n = ", NDAT, " detections)"), line=2,
        side=3, adj=0.5, cex=1.4, col="grey40")
    points(city, pch=18, col="grey10")
    text(city, rownames(city), cex=0.8, adj=-0.1, col="grey10")

	dev.off()
}
}


## veghf-north
## soilhf-treed-south
## soilhf-nontreed-south
## linear-north
## linear-south
## table: veghf-north
## table: soilhf-north

resn <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-north_", spp, ".Rdata")))
ress <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-south_", spp, ".Rdata")))
estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
ests_hab <- getEst(ress, stage=stage_hab_s, na.out=FALSE, Xns)
prn <- pred_veghf(estn_hab, Xnn)
prs <- pred_soilhf(ests_hab, Xns)

## FIXME produce plots / save tables--------------------------------------------- FIXME

## surroundinghf-north
## surroundinghf-south
## table: residual climate coefs north (with surrounding hf)
## table: residual climate coefs south (with surrounding hf)

## climate & surrounding hf
resn <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-north_", spp, ".Rdata")))
ress <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-south_", spp, ".Rdata")))
cn <- c("xPET", "xMAT", "xAHM", "xFFP", 
    "xMAP", "xMWMT", "xMCMT", "xlat", "xlong", "xlat2", "xlong2", 
    "THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", 
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", 
    "Noncult2_KM")
estn_sp <- getEst(resn, stage=stage_hab_n + 2, na.out=FALSE, Xnn)
ests_sp <- getEst(ress, stage=stage_hab_s + 2, na.out=FALSE, Xns)
sp_n <- colMeans(estn_sp[,cn])
sp_s <- colMeans(ests_sp[,cn])

## FIXME surrounding plot--------------------------------------------- FIXME
## FIXME save tables--------------------------------------------- FIXME

## trend-north
## trend-south
## table: north and south trend estimates

res_trend <- matrix(NA, nrow(tax), 10)
colnames(res_trend) <- c("Mean_North","Median_North","LCL_North","UCL_North","n_North",
    "Mean_South","Median_South","LCL_South","UCL_South","n_South")
res_trend[,5] <- tax$ndet_n
res_trend[,10] <- tax$ndet_s
rownames(res_trend) <- rownames(tax)
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
if (tax[spp, "trend_north"]) {
    resn <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-north_", spp, ".Rdata")))
    estn_yr <- getEst(resn, stage=stage_hab_n + 3, na.out=FALSE, Xnn)
    yr_n <- 100 * (exp(estn_yr[,"YR"]) - 1)
    res_trend[spp, 1:4] <- fstat(yr_n)
    NDATN <- sum(yyn[,spp] > 0)
    NN <- aggregate(yyn[,spp], list(year=xnn$YEAR), mean)
}
if (tax[spp, "trend_south"]) {
    ress <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-south_", spp, ".Rdata")))
    ests_yr <- getEst(ress, stage=stage_hab_s + 3, na.out=FALSE, Xns)
    yr_s <- 100 * (exp(ests_yr[,"YR"]) - 1)
    res_trend[spp, 6:9] <- fstat(yr_s)
    NDATS <- sum(yys[,spp] > 0)
    NS <- aggregate(yys[,spp], list(year=xns$YEAR), mean)
}
if (tax[spp, "trend_north"] | tax[spp, "trend_south"]) {
    NAM <- as.character(tax[spp, "English_Name"])
    fname <- file.path(ROOT, "figs", "trend", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    png(file=fname, width=600, height=600)
    op <- par(mfrow=c(2,2), cex=0.8)
    if (tax[spp, "trend_north"]) {
        plot(NN, ylab="Annual Mean Abundance Index", xlab="Year", 
            type="b", col=1, pch=19,
            main=paste0(NAM, ", North (n = ", NDATN, " detections)"))
        abline(lm(x ~ year, NN), col="red4", lty=1, lwd=2)
        hist(yr_n, col="gold", xlab="Decadal Trend (%)", main="")
        abline(v=fstat(yr_n)[1], col="red4", lty=1, lwd=2)
        abline(v=fstat(yr_n)[3:4], col="red4", lty=2, lwd=1)
    } else {
        plot.new()
        plot.new()
    }
    if (tax[spp, "trend_south"]) {
        plot(NS, ylab="Annual Mean Abundance Index", xlab="Year", 
            type="b", col=1, pch=19,
            main=paste0(NAM, ", South (n = ", NDATS, " detections)"))
        abline(lm(x ~ year, NS), col="red4", lty=1, lwd=2)
        hist(yr_n, col="gold", xlab="Decadal Trend (%)", main="")
        abline(v=fstat(yr_s)[1], col="red4", lty=1, lwd=2)
        abline(v=fstat(yr_s)[3:4], col="red4", lty=2, lwd=1)
    } else {
        plot.new()
        plot.new()
    }
    par(op)
    dev.off()
}
}
res_trend2 <- data.frame(tax[,c("English_Name","Scientific_Name")], res_trend)
write.csv(res_trend2, file=file.path(ROOT, "figs", "trend.csv"))

rank_fun(res_trend$Mean_North, res_trend$LCL_North, res_trend$UCL_North,
    n=res_trend$n_North, col=1, lab = rownames(res_trend))

rank_fun(res_trend$Mean_South, res_trend$LCL_South, res_trend$UCL_South,
    n=res_trend$n_South, col=1, lab = rownames(res_trend))

