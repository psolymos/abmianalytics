library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))

load(file.path(ROOT, "out", "kgrid", "kgrid_forSites.Rdata"))
kgrid0 <- kgrid[gis$closest_rowcol,]
rownames(kgrid0) <- gis$SITE_ID
yrs <- as.character(c(1999, 2001, 2002:2014))
kgrid <- kgrid0
    #rownames(kgrid) <- paste0(rownames(kgrid), ":", yrs[1])
    #kgrid$Year <- as.integer(yrs[1])
    #for (i in 2:length(yrs)) {
    #    kgridx <- kgrid0
    #    rownames(kgridx) <- paste0(rownames(kgridx), ":", yrs[i])
    #    kgridx$Year <- as.integer(yrs[i])
    #    kgrid <- rbind(kgrid, kgridx)
    #}

regs <- yrs # levels(kgrid$LUFxNSR)
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
table(gis$NATURAL_REGIONS, kgrid$NRNAME)

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


#spp <- "ALFL"
SPP <- union(fln, fls)
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL")

res_yrs <- list()
for (spp in SPP) {

cat(spp, "\n");flush.console()

load(file.path(ROOT, "out", "birds", "pred3x7", spp, paste0(regs[1], ".Rdata")))
pxNcr <- pxNcr1
pxScr <- pxScr1
for (i in 2:length(regs)) {
    load(file.path(ROOT, "out", "birds", "pred3x7", spp, paste0(regs[i], ".Rdata")))
    pxNcr <- cbind(pxNcr, pxNcr1)
    pxScr <- cbind(pxScr, pxScr1)
}
colnames(pxNcr) <- colnames(pxScr) <- yrs
rownames(pxNcr) <- rownames(pxScr) <- rownames(kgrid)
#pxNcr <- pxNcr[rownames(kgrid),]
#pxScr <- pxScr[rownames(kgrid),]
res_yrs[[spp]] <- list(vegcr=pxNcr, soilcr=pxScr)

}

save(res_yrs, file=file.path(ROOT, "out", "birds", "tables", "yearly3x7_changes.Rdata"))

plot_yrs <- function(spp) {
    x <- res_yrs[[spp]]
    #iis <- kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood"
    #iis[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- FALSE
    #iin <- kgrid$NRNAME != "Grassland"
    xx <- (kgrid$pAspen) * x$veg + (1-kgrid$pAspen) * x$soil
    xx <- groupSums(xx, 1, kgrid$NRNAME)
    xx <- 100 * xx / sum(xx[,1])

    matplot(as.integer(colnames(xx)), t(xx), type="b", pch=19, lwd=2, lty=1,
        main=tax[spp,"English_Name"])
    legend("topleft", col=1:6, lty=1, lwd=2, pch=19, legend=rownames(xx))
    xx
}

plot_yrs2 <- function(spp) {
    x <- res_yrs[[spp]]
    xx <- (kgrid$pAspen) * x$veg + (1-kgrid$pAspen) * x$soil
    xx <- colSums(xx)
    xx <- 100 * xx / sum(xx[1])

    plot(as.integer(names(xx)), xx, type="b", pch=19, lwd=2, lty=1,
        main=tax[spp,"English_Name"])
    xx
}


zz <- groupMeans(dd1km_pred$veg_current, 1, kgrid$Year)
zz <- 100 * zz / (21*10^6)
pdf("habits.pdf", onefile=TRUE)
for (i in 1:ncol(zz)){
    plot(as.integer(rownames(zz)), zz[,i], type="l", main=colnames(zz)[i])
    points(as.integer(rownames(zz)), zz[,i])
}
dev.off()
