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
SPP <- sort(union(fln, fls))
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

plot_yrs <- function(spp, all=TRUE) {
    x <- res_yrs[[spp]]
    xx <- (kgrid$pAspen) * x$veg + (1-kgrid$pAspen) * x$soil
    if (all) {
        xxx <- colSums(xx)
        xxx <- 100 * xxx / sum(xxx[1])
        LIM <- c(0, max(xxx)*1.2)
        plot(as.integer(names(xxx)), xxx, type="b", pch=19, lwd=2, lty=1,
            main=tax[spp,"English_Name"], ylim=LIM,
            axes=FALSE, ylab="Rel. abundance compared to 1999", xlab="Years")
    } else {
        xxx <- groupSums(xx, 1, kgrid$NRNAME)
        xxx <- 100 * xxx / sum(xxx[,1])
        LIM <- c(0, max(xxx)*1.2)
        matplot(as.integer(colnames(xxx)), t(xxx), type="b", pch=19, lwd=2, lty=1,
            main=tax[spp,"English_Name"], ylim=LIM,
            axes=FALSE, ylab="Rel. abundance compared to 1999", xlab="Years")
        legend("topleft", col=1:6, lty=1, lwd=2, pch=19, legend=rownames(xxx))
    }
    axis(1, at=as.integer(yrs), labels=yrs)
    axis(2)
    box()
    invisible(xxx)
}

ch_fun <- function(spp) {
    x <- res_yrs[[spp]]
    xx <- (kgrid$pAspen) * x$veg + (1-kgrid$pAspen) * x$soil
    xxx <- colSums(xx)
    xxx <- 100 * xxx / sum(xxx[1])
    n <- length(xxx)
    yr <- as.integer(names(xxx))
    unname(100 * ((xxx[n]/xxx[1])^(1/(yr[n] - yr[1])) - 1))
}

pdf("birds-yrs.pdf", onefile=TRUE)
for (spp in SPP) {
op <- par(mfrow=c(2,1))
try(plot_yrs(spp, TRUE))
try(plot_yrs(spp, FALSE))
par(op)
}
dev.off()

ch <- sapply(SPP, ch_fun)
ch <- ch[!is.na(ch) & ch < 10]
data.frame(ch=round(sort(ch),2))
## todo
- exclude spp with bad models
- calculate annual % change
