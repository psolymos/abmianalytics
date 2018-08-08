library(sp)
library(raster)
library(rgdal)

ROOT <- "e:/peter/bam/climate-preds"

## Al-Pac FMA boundary
fma <- readOGR(dsn="e:/peter/AB_data_v2018/data/raw/xy/alpac", "AlpacFMA2018")
rr <- raster("e:/peter/bam/climate-preds/OVEN/current/OVEN_currmean.asc")
fma <- spTransform(fma, proj4string(rr))
## Bird species list
SPP <- list.files(ROOT)

#spp <- "OVEN"
N <- list()
for (spp in SPP) {
    cat("\n", spp);flush.console()

fl <- structure(sprintf(c("e:/peter/bam/climate-preds/%s/current/%s_currmean.asc",
    "e:/peter/bam/climate-preds/%s/future/%s_2020mean.asc",
    "e:/peter/bam/climate-preds/%s/future/%s_2050mean.asc",
    "e:/peter/bam/climate-preds/%s/future/%s_2080mean.asc"), spp, spp),
    .Names = c("Baseline", "Yr_2020", "Yr_2050", "Yr_2080"))

    rs <- list()
    for (i in 1:4) {
        cat(".");flush.console()
        r0 <- raster(fl[i])
        rs[[names(fl)[i]]] <- trim(mask(r0, fma))
    }
    ## 4000x4000 cells: 1 cell is 40^2=1600 ha, N_cell = D males/ha * 1600 ha
    N[[spp]] <- round(sapply(rs, function(z) sum(values(z), na.rm=TRUE)) * 1600)
    rs <- stack(rs)
    writeRaster(rs, file=file.path(ROOT, spp, paste0(spp, "_clim-all_alpac-fma.tif")),
        overwrite=TRUE)

}

NN <- do.call(rbind, N)
write.csv(data.frame(Species=rownames(NN), NN), row.names=FALSE,
    file="e:/peter/bam/climate-preds/AlPacFMA_bird-climate-projections.csv")

pdf("e:/peter/bam/climate-preds/AlPacFMA_bird-climate-graphs.pdf",
    height=5, width=7, onefile=TRUE)
for (spp in SPP) {
    v <- 100*N[[spp]]/N[[spp]][1]
    plot(c(1990, 2020, 2050, 2080), v, type="n",
        ylim=c(0, max(v)+(50-floor(max(v) %% 50))),
        axes=FALSE, main=spp, xlab="Time", ylab="Abundance (% of baseline)")
    abline(h=100, col="grey", lty=2)
    lines(c(1990, 2020, 2050, 2080), v, lwd=3, pch=4, lty=1, col="lightblue")
    points(c(1990, 2020, 2050, 2080), v, pch=19, col=4, cex=1.5)
    axis(2)
    axis(1, c(1990, 2020, 2050, 2080), c(1990, 2020, 2050, 2080), tick=FALSE)
}
dev.off()


## save maps as pdf

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
pdf("e:/peter/bam/climate-preds/AlPacFMA_bird-climate-maps.pdf",
    height=10, width=18, onefile=TRUE)
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    r <- list(raster(file.path(ROOT, spp, paste0(spp, "_clim-all_alpac-fma.tif")), 1),
        raster(file.path(ROOT, spp, paste0(spp, "_clim-all_alpac-fma.tif")), 2),
        raster(file.path(ROOT, spp, paste0(spp, "_clim-all_alpac-fma.tif")), 3),
        raster(file.path(ROOT, spp, paste0(spp, "_clim-all_alpac-fma.tif")), 4))
    names(r) <- paste0(spp, c("_Baseline", "_2020", "_2050", "_2080"))
    r <- stack(r)
    MAX <- max(values(r[[1]]), values(r[[2]]), values(r[[3]]), values(r[[4]]), na.rm=TRUE)
    for (i in 1:4) {
        v <- values(r[[i]])
        v[1] <- MAX
        values(r[[i]]) <- v
    }

    op <- par(mfrow=c(2,3), mar=c(3,3,3,6))
    plot(r[[1]], col=col1, axes=FALSE, box=FALSE, main=names(r)[1], legend=TRUE)
    plot.new()
    plot.new()
    plot(r[[2]], col=col1, axes=FALSE, box=FALSE, main=names(r)[2], legend=FALSE)
    plot(r[[3]], col=col1, axes=FALSE, box=FALSE, main=names(r)[3], legend=FALSE)
    plot(r[[4]], col=col1, axes=FALSE, box=FALSE, main=names(r)[4], legend=FALSE)
    par(op)
}
dev.off()
