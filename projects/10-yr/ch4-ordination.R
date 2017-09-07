library(mefa4)
library(intrval)
library(MASS)
library(spatstat)
load("~/Dropbox/abmi/10yr/R/AllTables.Rdata")
AllIn$vplants$lt$pOcc <- AllIn$vplants$lt$nSites / 1598
source("~/repos/abmianalytics/projects/10-yr/ch4-functions.R")

colTd <- RColorBrewer::brewer.pal(6, "Dark2")
names(colTd) <- names(AllIn)
colTl <- paste0(colTd, "80")

ofun <- function(WHAT, NORTH=TRUE, SS=names(AllIn), PUT="none", LEGEND=NULL, ...) {
    mn <- ord_fun(WHAT, NORTH, alpha=0.5, scaling=2, col.text=1,
        col.pts=NA, cex.pts=0.5, cex.text=1, ...)
    s <- scores(mn, 1:3)
    TT <- if (NORTH)
        "veg" else "soilnt"
    tax <- do.call(c, lapply(names(AllIn), function(zz) rep(zz, nrow(AllIn[[zz]][[TT]]))))
    names(tax) <- do.call(c, lapply(AllIn, function(zz) as.character(zz[[TT]][,1])))
    tax <- tax[rownames(s$species)]
    xy <- lapply(structure(SS, names=SS), function(i) {
        j <- tax == i
        s$species[j,1:2]
    })
    for (i in 1:length(xy))
        points(xy[[i]], col=paste0(colTd[SS][i], "80"), pch=19, cex=1)
    chn <- lapply(structure(SS, names=SS), function(i) {
        j <- tax == i
        xx <- s$species[j,1:2]
        xx[chull(xx),]
    })
    ppn <- lapply(structure(SS, names=SS), function(i) {
        j <- tax == i
        as.ppp(s$species[j,1:2], c(-2,4,-2,2)) #c(range(s$species[,1]), range(s$species[,2])))
    })
    ell <- lapply(structure(SS, names=SS), function(i) {
        j <- tax == i
        X <- s$species[j,1:2]
        W <- weights(mn, display="species")[j]
        mat <- cov.wt(X, W)
        t <- sqrt(qchisq(0.95, 2))
        xy <- vegan:::veganCovEllipse(mat$cov, mat$center, t)
        xy
    })
    if ("chull" %in% PUT)
        for (i in 1:length(chn))
            polygon(chn[[i]], border=colTd[SS][i], col=paste0(colTd[SS][i], "10"))
    if ("ppp" %in% PUT)
        for (i in 1:length(ppn))
            contour(dfun(ppn[[i]]), levels=0.05, labels=names(ppn)[i], col=colTd[SS][i], add=TRUE)
    if ("ell" %in% PUT)
        for (i in 1:length(ell))
            lines(ell[[i]], col=colTd[SS][i])
    if (!is.null(LEGEND))
        legend(LEGEND, lty=1, lwd=2, col=colTd[SS], legend=names(chn), bty="n")
    invisible(NULL)
}

pdf("~/Dropbox/abmi/10yr/ch4/figs/ord-all-N.pdf", height=10, width=15)
par(mfrow=c(2,3), las=1)
ofun("All", NORTH=TRUE, PUT=c("chull"), SS="mammals", main="Mammals, North")
ofun("All", NORTH=TRUE, PUT=c("chull"), SS="birds", main="Birds, North")
ofun("All", NORTH=TRUE, PUT=c("chull"), SS="vplants", main="Vascular Plants, North")
ofun("All", NORTH=TRUE, PUT=c("chull"), SS="mosses", main="Bryophytes, North")
ofun("All", NORTH=TRUE, PUT=c("chull"), SS="lichens", main="Lichens, North")
ofun("All", NORTH=TRUE, PUT=c("chull"), SS="mites", main="Mites, North")
dev.off()

pdf("~/Dropbox/abmi/10yr/ch4/figs/ord-all-S.pdf", height=10, width=15)
par(mfrow=c(2,3), las=1)
ofun("All", NORTH=FALSE, PUT=c("chull"), SS="mammals", main="Mammals, South")
ofun("All", NORTH=FALSE, PUT=c("chull"), SS="birds", main="Birds, South")
ofun("All", NORTH=FALSE, PUT=c("chull"), SS="vplants", main="Vascular Plants, South")
ofun("All", NORTH=FALSE, PUT=c("chull"), SS="mosses", main="Bryophytes, South")
ofun("All", NORTH=FALSE, PUT=c("chull"), SS="lichens", main="Lichens, South")
ofun("All", NORTH=FALSE, PUT=c("chull"), SS="mites", main="Mites, South")
dev.off()

ofun1 <- function(WHAT, NORTH=TRUE, ...) {
    mn <- ord_fun(WHAT, NORTH, alpha=0.5, scaling=2, col.text=1,
        col.pts=NA, cex.pts=0.5, cex.text=1, plot=FALSE)
    s <- scores(mn, 1:3)
    xy <- s$species[,1:2]
    xy2 <- s$sites[,1:2]
    MM <- list(x=range(xy[,1], xy2[,1]), y=range(xy[,2], xy2[,2]))
    MM$x <- MM$x + 0.05*diff(MM$x)*c(-1,1)
    MM$y <- MM$y + 0.05*diff(MM$y)*c(-1,1)
    ppn <- as.ppp(xy, unlist(MM))
    k <- kde2d(x=xy[,1], y=xy[,2], n=100, lims=c(MM$x, MM$y))

    par(yaxs="i", xaxs="i", las=1)
    plot(0, ylim=MM$y, xlim=MM$x, xlab="Axis 1", ylab="Axis 2", type="n", ...)
    ## ppp
    plot(density(ppn), add=TRUE, col=colorRampPalette(c("#FFFFFF", colTd[WHAT]))(100))
    contour(dfun(ppn), col=colTd[WHAT], add=TRUE)
    ## kde
    #image(k, add=TRUE, col=colorRampPalette(c("#FFFFFF", colTd[WHAT]))(100))
    #contour(k, col=colTd[WHAT], add=TRUE)
#    cumz<-cumsum(sort(k$z,decreasing=TRUE))/sum(k$z)
#    z90<-sort(k$z,decreasing=TRUE)[which.min(abs(cumz-0.9))]
#    contour(k,levels=z90,col=colTd[WHAT],add=TRUE,labels="90%")

    points(xy, col="#00000040", pch=19, cex=1)
    text(xy2*0.8, labels=rownames(xy2), col=1, pch=19, cex=0.9)
    abline(h=0,v=0,lty=3)
    box()
}
pdf("~/Dropbox/abmi/10yr/ch4/figs/ord-one-N-ppp.pdf", height=8, width=12)
par(mfrow=c(2,3))
ofun1("mammals", NORTH=TRUE, main="Mammals, North")
ofun1("birds", NORTH=TRUE, main="Birds, North")
ofun1("vplants", NORTH=TRUE, main="Vascular Plants, North")
ofun1("mosses", NORTH=TRUE, main="Bryophytes, North")
ofun1("lichens", NORTH=TRUE, main="Lichens, North")
ofun1("mites", NORTH=TRUE, main="Mites, North")
dev.off()

pdf("~/Dropbox/abmi/10yr/ch4/figs/ord-one-S-ppp.pdf", height=8, width=12)
par(mfrow=c(2,3))
ofun1("mammals", NORTH=FALSE, main="Mammals, South")
ofun1("birds", NORTH=FALSE, main="Birds, South")
ofun1("vplants", NORTH=FALSE, main="Vascular Plants, South")
ofun1("mosses", NORTH=FALSE, main="Bryophytes, South")
ofun1("lichens", NORTH=FALSE, main="Lichens, South")
ofun1("mites", NORTH=FALSE, main="Mites, South")
dev.off()

