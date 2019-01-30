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

cl <- makeCluster(8)
for (i in tx) {
    cat(i, "\n")
    rn <- rownames(ABMI$detections[[i]])
    yy <- ABMI$detections[[i]][rn,]
    if (i == "vascular_plants")
        yy <- yy[,!ABMI$species[[i]]$nonnative]
    yy[yy > 0] <- 1
    yy <- yy[,colSums(yy) > 0]
    xx <- ABMI$sites[rn,]
    xx$g1 <- xx$NRNAME
    levels(xx$g1) <- c("B", "S", "F", "G", "P", "R")

    oc1 <- opticut(yy ~ 1, strata=xx$g1, dist="binomial", cl=cl)
    bp1 <- summary(oc1)$bestpart
    BP[[i]] <- bp1
}
stopCluster(cl)


op <- par(mfrow=c(2,3))
for (i in tx)
plot(hclust(vegan::vegdist(t(BP[[i]]), "bray"), "ward.D2"), xlab="Natural Regions", sub="", main=i)
par(op)

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

SppTab <- get_species_table()
SppTab <- SppTab[SppTab$native,]
#Spp <- intersect(SppTab$SpeciesID, unname(unlist(lapply(bpl, rownames))))
Spp <- intersect(SppTab$SpeciesID, rownames(bp))
bp <- bp[Spp,]
SppTab <- SppTab[Spp,]
rt <- .read_raster_template()
KT <- cure4insect:::.c4if$KT

rasterize_results_cr <- function (y)
{
    NC <- rowSums(y$SA.Curr)
    KT$NC <- NC[match(rownames(KT), names(NC))]
    KT$NC[is.na(KT$NC)] <- 0
    r <- .make_raster(KT$NC, rc = KT, rt = rt)
    r
}

f <- function(spp) {
    y <- load_species_data(spp)
    r <- rasterize_results_cr(y)
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
All <- pblapply(Spp, f, cl=cl)
stopCluster(cl)
names(All) <- Spp

#bpl$mites <- bpl$mites[,c("BS", "GP", "FR")] # hack to reorg mites
#bp <- do.call(rbind, bpl)[Spp,]

Rall <- Reduce("+", All)
Rred   <- Reduce("+", All[bp[,"GP"]==1])
Rgreen <- Reduce("+", All[bp[,"SB"]==1])
Rblue  <- Reduce("+", All[bp[,"FR"]==1])

RRall <- stack(list(r=255*Rred/Rall, g=255*Rgreen/Rall, b=255*Rblue/Rall))
R0 <- stack(list(r=0*rt, g=0*rt, b=0*rt))

plotRGB(R0)
plotRGB(RRall, stretch="hist")

h <- function(TT) {
    ss <- SppTab$taxon == TT
    Rall <- Reduce("+", All[ss])
    Rred   <- Reduce("+", All[bp[,"GP"]==1 & ss])
    Rgreen <- Reduce("+", All[bp[,"SB"]==1 & ss])
    Rblue  <- Reduce("+", All[bp[,"FR"]==1 & ss])
    RR <- stack(list(r=255*Rred/Rall, g=255*Rgreen/Rall, b=255*Rblue/Rall))
}

RRbirds <- h("birds")
RRlichens <- h("lichens")
RRmites <- h("mites")
RRmosses <- h("mosses")
RRvplants <- h("vplants")


hacked_plotRGB <-
function(x, r=1, g=2, b=3, scale, maxpixels=500000, stretch=NULL, ext=NULL, interpolate=FALSE, colNA='white', alpha, bgalpha, addfun=NULL, zlim=NULL, zlimcol=NULL, axes=FALSE, xlab='', ylab='', asp=NULL, add=FALSE, ...) {

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

op <- par(mfrow=c(2,3), mar=c(2,2,2,2))
hacked_plotRGB(RRall, stretch="hist", main="All")
hacked_plotRGB(RRbirds, stretch="hist", main="Birds")
hacked_plotRGB(RRlichens, stretch="hist", main="Lichens")
hacked_plotRGB(RRmites, stretch="hist", main="Mites")
#plot.new()
hacked_plotRGB(RRmosses, stretch="hist", main="Bryophytes")
hacked_plotRGB(RRvplants, stretch="hist", main="Vascular Plants")
par(op)

op <- par(mfrow=c(2,3), mar=c(2,2,2,2))
hacked_plotRGB(RRall, main="All")
hacked_plotRGB(RRbirds, main="Birds")
hacked_plotRGB(RRlichens, main="Lichens")
hacked_plotRGB(RRmites, main="Mites")
#plot.new()
hacked_plotRGB(RRmosses, main="Bryophytes")
hacked_plotRGB(RRvplants, main="Vascular Plants")
par(op)

save(RRall, RRbirds, RRlichens, RRmites, RRmosses, RRvplants,
    file="d:/abmi/sppweb2018/rainbow-maps/results-rainbow-maps-2019-01-30.RData")
