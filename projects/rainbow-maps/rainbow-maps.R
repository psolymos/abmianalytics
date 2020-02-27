#' ---
#' title: "Rainbow maps"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Dec 12, 2018"
#' output: word_document
#' ---
#'
knitr::opts_chunk$set(eval=FALSE)
#'
#'
#' # Appendix: R code to reproduce results
#'
#' Loading required R packages
library(mefa4) # data manipulation
library(opticut) # optimal partitioning
library(intrval) # easy interval calculations
library(parallel) # speed things up
#' Settings local path for the cure4insect package to speed things up, then load common data
library(cure4insect)
#opar <- set_options(path = "w:/reports")
opar <- set_options(path = "d:/abmi/reports")
load_common_data()
#' Load a precompiled data of ABMI species data
load("~/GoogleWork/collaborations/opticut/R/abmi-data/abmi-data-AB.Rdata")

#' ## Combine all taxa

#' Retrieve the intersect of all the rows in the different taxa tables (sites)
rn <- Reduce("intersect", lapply(ABMI$detections, rownames))
#' Make a giant site x species matrix
yy <- do.call(cbind, lapply(ABMI$detections, function(z) z[rn,]))
#' Treat entries as presence (1) ansence (0)
yy[yy > 0] <- 1
#' Drop empty columns
yy <- yy[,colSums(yy) > 0]
#' `S` is overall number of species remining in `yy`
S <- ncol(yy)
#' Get Natural Region names for each site and shroten region names
xx <- ABMI$sites[rn,]
xx$g1 <- xx$NRNAME
levels(xx$g1) <- c("B", "S", "F", "G", "P", "R")
#' Use an 8-core cluster
cl <- makeCluster(8)
#' Fit opticut model with all the species (all taxa combined):
#' - natural regions used as strata,
#' - fit logistic regression to 0/1 data.
oc1 <- opticut(yy ~ 1, strata=xx$g1, dist="binomial", cl=cl)
#' Get the optimal partition indicator matrix (species x regions)
bp1 <- summary(oc1)$bestpart

plot(oc1, cex.axis=0.1)

sh <- t(bp1) %*% bp1
#round(as.dist(100 * sh / S, 1))

plot(hclust(dist(t(bp1), "manhattan"), "ward.D2"), xlab="Natural Regions", sub="")
plot(hclust(vegan::vegdist(t(bp1), "bray"), "ward.D2"), xlab="Natural Regions", sub="")

xx$g2 <- xx$NRNAME
Class <- c("SB", "SB", "FR", "GP", "GP", "FR")
levels(xx$g2) <- Class
oc2 <- opticut(yy ~ 1, strata=xx$g2, dist="binomial", cl=cl)
bp2 <- summary(oc2)$bestpart

sh <- t(bp2) %*% bp2
round(as.dist(100 * sh / S, 1))

stopCluster(cl)

bp <- groupSums(bp1, 2, Class)
bp[bp > 0] <- 1

#' ## Each taxon separate OC

tx <- names(ABMI$detections)

i <- tx[1]
BP <- list()
OC <- list()

cl <- makeCluster(8)
for (i in tx) {
    cat(i, "\n")
    rn <- rownames(ABMI$detections[[i]])
    yy <- ABMI$detections[[i]][rn,]
#    if (i == "vascular_plants")
#        yy <- yy[,!ABMI$species[[i]]$nonnative]
    yy[yy > 0] <- 1
    yy <- yy[,colSums(yy) > 0]
    xx <- ABMI$sites[rn,]
    xx$g1 <- xx$NRNAME
    levels(xx$g1) <- c("B", "S", "F", "G", "P", "R")

    oc1 <- opticut(yy ~ 1, strata=xx$g1, dist="binomial", cl=cl)
    bp1 <- summary(oc1)$bestpart
    BP[[i]] <- bp1
    OC[[i]] <- oc1
}
stopCluster(cl)


cl <- makeCluster(8)

i <- "vascular_plants"

rn <- rownames(ABMI$detections[[i]])

yy <- ABMI$detections[[i]][rn,]
yy <- yy[,!ABMI$species[[i]]$nonnative]
yy[yy > 0] <- 1
yy <- yy[,colSums(yy) > 0]
xx <- ABMI$sites[rn,]
xx$g1 <- xx$NRNAME
levels(xx$g1) <- c("B", "S", "F", "G", "P", "R")
oc1n <- opticut(yy ~ 1, strata=xx$g1, dist="binomial", cl=cl)
bp1n <- summary(oc1n)$bestpart

yy <- ABMI$detections[[i]][rn,]
yy <- yy[,ABMI$species[[i]]$nonnative]
yy[yy > 0] <- 1
yy <- yy[,colSums(yy) > 0]
xx <- ABMI$sites[rn,]
xx$g1 <- xx$NRNAME
levels(xx$g1) <- c("B", "S", "F", "G", "P", "R")

oc1nn <- opticut(yy ~ 1, strata=xx$g1, dist="binomial", cl=cl)
bp1nn <- summary(oc1nn)$bestpart

stopCluster(cl)

bp1n <- summary(oc1n)$bestpart
bp1nn <- summary(oc1nn)$bestpart


pdf("~/GoogleWork/abmi/rainbow-maps/clusters.pdf", height=7, width=12)
op <- par(mfrow=c(2,4))
for (i in tx[c(1,2,5,4,3)])
plot(hclust(vegan::vegdist(t(BP[[i]]), "bray"), "ward.D2"), xlab="Natural Regions", sub="", main=i,
    hang=-1)
plot(hclust(vegan::vegdist(t(bp1n), "bray"), "ward.D2"), xlab="Natural Regions", sub="",
    main=paste(i, "(native)"), hang=-1)
plot(hclust(vegan::vegdist(t(bp1nn), "bray"), "ward.D2"), xlab="Natural Regions", sub="",
    main=paste(i, "(non-native)"), hang=-1)
par(op)
dev.off()

pdf(paste0("~/GoogleWork/abmi/rainbow-maps/indicators.pdf"), height=10, width=10)
for (i in tx) {
plot(OC[[i]], cex.axis=0.3, show_I=FALSE, show_S=FALSE, mar=c(4,8,2,2), xlab="", ylab="", bty="n")
title(main=i)
}
plot(oc1n, cex.axis=0.3, show_I=FALSE, show_S=FALSE, mar=c(4,8,4,2), xlab="", ylab="", bty="n")
title(main="vascular_plants (native)")
plot(oc1nn, cex.axis=0.3, show_I=FALSE, show_S=FALSE, mar=c(4,8,2,2), xlab="", ylab="", bty="n")
title(main="vascular_plants (non-native)")
dev.off()

bp <- groupSums(bp1, 2, c("B", "FSR", "FSR", "GP", "GP", "FSR"))
bp[bp > 0] <- 1

bpl <- list()
for (i in tx) {
    Class <- if (i == "mites")
        c("BS", "BS", "FR", "GP", "GP", "FR") else c("SBF", "SBF", "SBF", "GP", "GP", "R")
    bp <- groupSums(BP[[i]], 2, Class)
    bp[bp > 0] <- 1
    bpl[[i]] <- bp
}





#' Process predictions


bpl$mites <- bpl$mites[,c("BS", "GP", "FR")] # hack to reorg mites
bp <- do.call(rbind, bpl)

SppTab <- get_species_table()
SppTab$native[!SppTab$native & SppTab$taxon != "vplants"] <- TRUE
#SppTab <- SppTab[SppTab$native,]

#Spp <- intersect(SppTab$SpeciesID, unname(unlist(lapply(bpl, rownames))))
Spp <- intersect(SppTab$SpeciesID, rownames(bp))
bp <- bp[Spp,]
SppTab <- SppTab[Spp,]
rt <- .read_raster_template()
KT <- cure4insect:::.c4if$KT

rasterize_results_cr <- function (y, current=TRUE)
{
    NC <- if (current)
        rowSums(y$SA.Curr) else rowSums(y$SA.Ref)
    KT$NC <- NC[match(rownames(KT), names(NC))]
    KT$NC[is.na(KT$NC)] <- 0
    r <- .make_raster(KT$NC, rc = KT, rt = rt)
    r
}

f <- function(spp, current=TRUE) {
    y <- load_species_data(spp)
    r <- rasterize_results_cr(y, current)
    if (SppTab[spp,"taxon"] == "birds")
        r <- p_bird(r, "ha", 2)
    r
}

cl <- makeCluster(8)
system.time(f(Spp[1]))
tmp <- clusterEvalQ(cl, library(cure4insect))
tmp <- clusterEvalQ(cl, set_options(path = "d:/abmi/reports"))
tmp <- clusterEvalQ(cl, load_common_data())
tmp <- clusterExport(cl, c("rt", "KT", "rasterize_results_cr", "SppTab"))
#All1 <- pblapply(Spp, f, current=TRUE, cl=cl)
All2 <- pblapply(Spp, f, current=FALSE, cl=cl)
stopCluster(cl)
#names(All1) <- Spp
names(All2) <- Spp

All <- All2

Rall <- Reduce("+", All)
Rred   <- Reduce("+", All[bp[,"GP"]==1])
Rgreen <- Reduce("+", All[bp[,"SBF"]==1])
Rblue  <- Reduce("+", All[bp[,"R"]==1])

RRall <- stack(list(r=255*Rred/Rall, g=255*Rgreen/Rall, b=255*Rblue/Rall))
R0 <- stack(list(r=0*rt, g=0*rt, b=0*rt))

plotRGB(R0)
plotRGB(RRall, stretch="hist", add=TRUE)

h <- function(TT, what="both") {
    if (what == "native")
        ns <- SppTab$native
    if (what == "nonnative")
        ns <- !SppTab$native
    if (what == "both")
        ns <- TRUE
    ss <- SppTab$taxon == TT & ns
    Rall <- Reduce("+", All[ss])
    Rred   <- Reduce("+", All[bp[,"GP"]==1 & ss])
    Rgreen <- Reduce("+", All[bp[,"SBF"]==1 & ss])
    Rblue  <- Reduce("+", All[bp[,"R"]==1 & ss])
    if (what == "nonnative") {
        tmp <- Rgreen + Rblue
        Rblue <- Rgreen <- 0.5*tmp
    }
    RR <- stack(list(r=255*Rred/Rall, g=255*Rgreen/Rall, b=255*Rblue/Rall))
}

RRbirds <- h("birds")
RRlichens <- h("lichens")
RRmites <- h("mites")
RRmosses <- h("mosses")
RRvplants <- h("vplants")
RRvplantsNonly <- h("vplants", "native")
RRvplantsNNonly <- h("vplants", "nonnative")

hacked_plotRGB <-
function(x, r=1, g=2, b=3, scale, maxpixels=500000, stretch=NULL, ext=NULL, interpolate=FALSE, colNA='white',
alpha, bgalpha, addfun=NULL, zlim=NULL, zlimcol=NULL, axes=FALSE, xlab='', ylab='', asp=NULL, add=FALSE, ...) {

	if (missing(scale)) {
		scale <- 255
		if (! inherits(x, 'RasterStack')) {
			if ( x@data@haveminmax ) {
				scale <- max(max(x@data@max), 255)
			}
		}
	}
	scale <- as.vector(scale)[1]

	r <- sampleRegular(raster(x,r), maxpixels, ext=ext, asRaster=TRUE, useGDAL=TRUE)
	g <- sampleRegular(raster(x,g), maxpixels, ext=ext, asRaster=TRUE, useGDAL=TRUE)
	b <- sampleRegular(raster(x,b), maxpixels, ext=ext, asRaster=TRUE, useGDAL=TRUE)

	RGB <- cbind(getValues(r), getValues(g), getValues(b))

	if (!is.null(zlim)) {
		if (length(zlim) == 2) {
			zlim <- sort(zlim)
			if (is.null(zlimcol)) {
				RGB[ RGB<zlim[1] ] <- zlim[1]
				RGB[ RGB>zlim[2] ] <- zlim[2]
			} else { #if (is.na(zlimcol)) {
				RGB[RGB<zlim[1] | RGB>zlim[2]] <- NA
			}
		} else if (NROW(zlim) == 3 & NCOL(zlim) == 2) {
			for (i in 1:3) {
				zmin <- min(zlim[i,])
				zmax <- max(zlim[i,])
				if (is.null(zlimcol)) {
					RGB[RGB[,i] < zmin, i] <- zmin
					RGB[RGB[,i] > zmax, i] <- zmax
				} else { #if (is.na(zlimcol)) {
					RGB[RGB < zmin | RGB > zmax, i] <- NA
				}
			}
		} else {
			stop('zlim should be a vector of two numbers or a 3x2 matrix (one row for each color)')
		}
	}

	RGB <- stats::na.omit(RGB)

	if (!is.null(stretch)) {
		stretch = tolower(stretch)
		if (stretch == 'lin') {
			RGB[,1] <- raster:::.linStretchVec(RGB[,1])
			RGB[,2] <- raster:::.linStretchVec(RGB[,2])
			RGB[,3] <- raster:::.linStretchVec(RGB[,3])
			scale <- 255
		} else if (stretch == 'hist') {
			RGB[,1] <- raster:::.eqStretchVec(RGB[,1])
			RGB[,2] <- raster:::.eqStretchVec(RGB[,2])
			RGB[,3] <- raster:::.eqStretchVec(RGB[,3])
			scale <- 255
		} else if (stretch != '') {
			warning('invalid stretch value')
		}
	}


	naind <- as.vector( attr(RGB, "na.action") )
	if (!is.null(naind)) {
		bg <- grDevices::col2rgb(colNA)
		bg <- grDevices::rgb(bg[1], bg[2], bg[3], alpha=bgalpha, max=255)
		z <- rep( bg, times=ncell(r))
		z[-naind] <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, max=scale)
	} else {
		z <- grDevices::rgb(RGB[,1], RGB[,2], RGB[,3], alpha=alpha, max=scale)
	}

	z <- matrix(z, nrow=nrow(r), ncol=ncol(r), byrow=T)

	requireNamespace("grDevices")
	bb <- as.vector(t(bbox(r)))


	if (!add) {
		#if (!axes) graphics::par(plt=c(0,1,0,1))

		if (is.null(asp)) {
			if (couldBeLonLat(x)) {
			    ym <- mean(c(x@extent@ymax, x@extent@ymin))
				asp <- 1/cos((ym * pi)/180)
				#asp <- min(5, 1/cos((ym * pi)/180))
			} else {
				asp <- 1
			}
		}

		xlim=c(bb[1], bb[2])
		ylim=c(bb[3], bb[4])

		plot(NA, NA, xlim=xlim, ylim=ylim, type = "n", xaxs='i', yaxs='i', xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, ...)
		if (axes) {
			xticks <- graphics::axTicks(1, c(xmin(r), xmax(r), 4))
			yticks <- graphics::axTicks(2, c(ymin(r), ymax(r), 4))
			if (xres(r) %% 1 == 0) xticks = round(xticks)
			if (yres(r) %% 1 == 0) yticks = round(yticks)
			graphics::axis(1, at=xticks)
			graphics::axis(2, at=yticks, las = 1)
			#graphics::axis(3, at=xticks, labels=FALSE, lwd.ticks=0)
			#graphics::axis(4, at=yticks, labels=FALSE, lwd.ticks=0)
		}
	}
	graphics::rasterImage(z, bb[1], bb[3], bb[2], bb[4], interpolate=interpolate, ...)

	if (!is.null(addfun)) {
		if (is.function(addfun)) {
			addfun()
		}
	}
    invisible(NULL)
}

stretch="lin"
#pdf("d:/abmi/sppweb2018/rainbow-maps/results-rainbow-maps-curr.pdf")
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-ref.pdf", height=12, width=8, onefile=TRUE)
#op <- par(mfrow=c(2,4), mar=c(2,2,2,2))
op <- par(mar=c(2,2,2,2))
hacked_plotRGB(RRall, stretch=stretch, main="All")
hacked_plotRGB(RRbirds, stretch=stretch, main="Birds")
hacked_plotRGB(RRlichens, stretch=stretch, main="Lichens")
hacked_plotRGB(RRmites, stretch=stretch, main="Mites")
hacked_plotRGB(RRmosses, stretch=stretch, main="Bryophytes")
hacked_plotRGB(RRvplants, stretch=stretch, main="Vascular Plants (all)")
hacked_plotRGB(RRvplantsNonly, stretch=stretch, main="Vascular Plants (native)")
hacked_plotRGB(RRvplantsNNonly, stretch=stretch, main="Vascular Plants (non-native)")
par(op)
dev.off()


op <- par(mar=c(2,2,2,2))
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-all.pdf", height=12, width=8)
hacked_plotRGB(RRall, stretch=stretch, main="All")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-birds.pdf", height=12, width=8)
hacked_plotRGB(RRbirds, stretch=stretch, main="Birds")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-lichens.pdf", height=12, width=8)
hacked_plotRGB(RRlichens, stretch=stretch, main="Lichens")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-mites.pdf", height=12, width=8)
hacked_plotRGB(RRmites, stretch=stretch, main="Mites")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-mosses.pdf", height=12, width=8)
hacked_plotRGB(RRmosses, stretch=stretch, main="Bryophytes")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-vplants-all.pdf", height=12, width=8)
hacked_plotRGB(RRvplants, stretch=stretch, main="Vascular Plants (all)")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-vplants-native.pdf", height=12, width=8)
hacked_plotRGB(RRvplantsNonly, stretch=stretch, main="Vascular Plants (native)")
dev.off()
pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-vplants-nonnative.pdf", height=12, width=8)
hacked_plotRGB(RRvplantsNNonly, stretch=stretch, main="Vascular Plants (non-native)")
dev.off()
par(op)

pdf("~/GoogleWork/abmi/rainbow-maps/results-rainbow-maps-mosses.pdf", height=12, width=8, onefile=TRUE)
hacked_plotRGB(RRmosses, stretch=NULL, main="Bryophytes")
hacked_plotRGB(RRmosses, stretch="lin", main="Bryophytes")
hacked_plotRGB(RRmosses, stretch="hist", main="Bryophytes")
plot(RRmosses, 1)
plot(RRmosses, 2)
plot(RRmosses, 3)
dev.off()



op <- par(mfrow=c(2,3), mar=c(2,2,2,2))
hacked_plotRGB(RRall, main="All")
hacked_plotRGB(RRbirds, main="Birds")
hacked_plotRGB(RRlichens, main="Lichens")
hacked_plotRGB(RRmites, main="Mites")
#plot.new()
hacked_plotRGB(RRmosses, main="Bryophytes")
hacked_plotRGB(RRvplants, main="Vascular Plants")
par(op)

#fno <- "d:/abmi/sppweb2018/rainbow-maps/results-rainbow-maps-curr-2019-01-30.RData"
fno <- "d:/abmi/sppweb2018/rainbow-maps/results-rainbow-maps-ref-2019-01-30.RData"
save(RRall, RRbirds, RRlichens, RRmites, RRmosses, RRvplants,
    file=fno)


## post processing the results

load("s:/sppweb2018/rainbow-maps/results-rainbow-maps-ref-2019-01-30.RData")



