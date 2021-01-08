## load coefficients

library(intrval)
load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData") # mammals & habitat elements & nnplants

## fancy labels for the coefficients

labs <- list(
    veg=list(
        native=c(`White spruce 0-9` = "WhiteSpruceR", `White spruce 10-19` = "WhiteSpruce1",
            `White spruce 20-39` = "WhiteSpruce2", `White spruce 40-59` = "WhiteSpruce3",
            `White spruce 60-79` = "WhiteSpruce4", `White spruce 80-99` = "WhiteSpruce5",
            `White spruce 100-119` = "WhiteSpruce6", `White spruce 120-139` = "WhiteSpruce7",
            `White spruce 140+` = "WhiteSpruce8", `Pine 0-9` = "PineR", `Pine 10-19` = "Pine1",
            `Pine 20-39` = "Pine2", `Pine 40-59` = "Pine3", `Pine 60-79` = "Pine4",
            `Pine 80-99` = "Pine5", `Pine 100-119` = "Pine6", `Pine 120-139` = "Pine7",
            `Pine 140+` = "Pine8", `Deciduous 0-9` = "DeciduousR", `Deciduous 10-19` = "Deciduous1",
            `Deciduous 20-39` = "Deciduous2", `Deciduous 40-59` = "Deciduous3",
            `Deciduous 60-79` = "Deciduous4", `Deciduous 80-99` = "Deciduous5",
            `Deciduous 100-119` = "Deciduous6", `Deciduous 120-139` = "Deciduous7",
            `Deciduous 140+` = "Deciduous8", `Mixedwood 0-9` = "MixedwoodR",
            `Mixedwood 10-19` = "Mixedwood1", `Mixedwood 20-39` = "Mixedwood2",
            `Mixedwood 40-59` = "Mixedwood3", `Mixedwood 60-79` = "Mixedwood4",
            `Mixedwood 80-99` = "Mixedwood5", `Mixedwood 100-119` = "Mixedwood6",
            `Mixedwood 120-139` = "Mixedwood7", `Mixedwood 140+` = "Mixedwood8",
            `Treed bog 0-9` = "TreedBogR", `Treed bog 10-19` = "TreedBog1",
            `Treed bog 20-39` = "TreedBog2", `Treed bog 40-59` = "TreedBog3",
            `Treed bog 60-79` = "TreedBog4", `Treed bog 80-99` = "TreedBog5",
            `Treed bog 100-119` = "TreedBog6", `Treed bog 120-139` = "TreedBog7",
            `Treed bog 140+` = "TreedBog8",
            `Treed fen 0-9`="TreedFenR", `Treed fen 10-19`="TreedFen1", `Treed fen 20-39`="TreedFen2",
            `Treed fen 40-59`="TreedFen3",`Treed fen 60-79`= "TreedFen4", `Treed fen 80-99`="TreedFen5",
            `Treed fen 100-119`="TreedFen6", `Treed fen 120-139`="TreedFen7", `Treed fen 140+`="TreedFen8",
            `Treed swamp` = "TreedSwamp", `Shrubby swamp` = "ShrubbySwamp",
            `Shrubby bog` = "ShrubbyBog", `Shrubby fen` = "ShrubbyFen", `Graminoid fen` = "GraminoidFen",
            `Marsh` = "Marsh", `Shrub` = "Shrub", `Grass/herb` = "GrassHerb"),
        hffor=c(`Harvest, White spruce 0-9` = "CCWhiteSpruceR",
            `Harvest, white spruce 10-19` = "CCWhiteSpruce1",
            `Harvest, white spruce 20-39` = "CCWhiteSpruce2",
            `Harvest, white spruce 40-59` = "CCWhiteSpruce3",
            `Harvest, white spruce 60-79` = "CCWhiteSpruce4",
            `Harvest, pine 0-9` = "CCPineR", `Harvest, pine 10-19` = "CCPine1",
            `Harvest, pine 20-39` = "CCPine2", `Harvest, pine 40-59` = "CCPine3",
            `Harvest, pine 60-79` = "CCPine4",
            `Harvest, deciduous 0-9` = "CCDeciduousR", `Harvest, deciduous 10-19` = "CCDeciduous1",
            `Harvest, deciduous 20-39` = "CCDeciduous2", `Harvest, deciduous 40-59` = "CCDeciduous3",
            `Harvest, deciduous 60-79` = "CCDeciduous4",
            `Harvest, mixedwood 0-9` = "CCMixedwoodR",
            `Harvest, mixedwood 10-19` = "CCMixedwood1", `Harvest, mixedwood 20-39` = "CCMixedwood2",
            `Harvest, mixedwood 40-59` = "CCMixedwood3", `Harvest, mixedwood 60-79` = "CCMixedwood4"),
        hfpoly=c(
            "Cropland"="Crop",
            "Tame pasture"="TameP",
            "Rough pasture"="RoughP",
            "Well sites"="Wellsites",
            "Rural residential"="Rural",
            "Urban/Industrial"="Urban",
            "Industrial (rural)"="Industrial"),
        hfsoft=c(
            "Seismic lines"="EnSeismic",
            "Vegetated linear (energy)"="EnSoftLin",
            "Vegetated linear (transport)"="TrSoftLin"),
        hfhard=c(
            "Non-vegetated linear"="HardLin")
    ),
    soil=list(
        native=c(
            "Loamy"="Loamy",
            "Sandy/loamy"="SandyLoam",
            "Rapid drain"="RapidDrain",
            "Clay"="ClaySub",
            "Thin break"="ThinBreak",
            "Blowout"="Blowout",
            "Other soil types"="Other"),
        hfpoly=c(
            "Cropland"="Crop",
            "Tame pasture"="TameP",
            "Rough pasture"="RoughP",
            "Well sites"="Wellsites",
            "Rural residential"="Rural",
            "Urban/Industrial"="Urban",
            "Industrial (rural)"="Industrial"),
        hfsoft=c(
            "Seismic lines"="EnSeismic",
            "Vegetated linear (energy)"="EnSoftLin",
            "Vegetated linear (transport)"="TrSoftLin"),
        hfhard=c(
            "Non-vegetated linear"="HardLin"),
        extras=c(
            "Mine (non-vegetated)"="Mine",
            "Mine (vegetated)"="MineV",
            "Open water"="Water"))
    )

## soil/hf plots: cf is a coef table (not spp name & taxa as will be in c4i)
## coef from the S model is expected
.plot_abundance_soil_2020 <- function(cf, plot=TRUE, paspen=0, ylim, main, ylab,
    link="logit", precalc=FALSE, ...) {

    if (paspen %)(% c(0, 1))
        stop("paspen must be in [0, 1]")
    lab <- c(labs$soil$native, labs$soil$hfpoly)

    #FUN <- if (bird)
    #    exp else plogis
    FUN <- binomial(link)$linkinv
    p <- t(t(cf[lab,]) + paspen*cf["pAspen",])
    p <- .get_stats(FUN(p), precalc=precalc)

    #USE <- "Median"
    USE <- "First"
    out <- as.matrix(p[,c(USE, "Lower", "Upper")])
    colnames(out) <- c("Estimate", "LCL", "UCL")
    rownames(out) <- names(lab)
    out[is.infinite(out)] <- NA
    out[out < 10^-5] <- 10^-6
    if (plot) {
        op <- par(mai=c(1.5, 1, 0.2, 0.3))
        on.exit(par(op))

        lci <- out[,"LCL"]
        uci <- out[,"UCL"]
        y1 <- out[,"Estimate"]
        x <- 1:14

        if (missing(main))
            main <- ""
        if (missing(ylab))
            ylab <- "Relative abundance"
        if (missing(ylim)) {
            ymax <- max(min(max(uci[x], na.rm=TRUE), 2*max(y1, na.rm=TRUE)), y1*1.02)
            ylim <- c(0, ymax)
        } else {
            ymax <- max(ylim)
        }

        space <- c(1, x[-1] -x[-length(x)]) - 0.9
        col <- c(
            hcl.colors(7,"Dark 3"),
            "#333300", "#666600", "#999900",
            "#909090","#555555","#404040","#222222")

        x1 <- barplot(y1[x], space=space, border="white", col=col,
            ylim=ylim, xlim=c(-0.5,max(x)+1.5), xaxs="i", yaxt="n", ylab=ylab,
            col.lab="grey50", cex.lab=1.2,axisnames=FALSE)[,1]
        ax <- axis(side=2, cex.axis=0.9, col.axis="grey50", col.ticks="grey50", las=2)
        abline(h=ax, col="grey80")
        x1 <- barplot(y1[x], space=space, border="white", col=col,
            xaxs="i", yaxt="n", ylab=ylab,
            col.lab="grey50", cex.lab=1.2,axisnames=FALSE, add=TRUE)[,1]
        box(bty="l", col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[i], y1[i]), col="grey90")
            lines(rep(x1[i],2), c(uci[i], y1[i]), col=col[i])
        }
        mtext(side=1, at=x1, line=0.7,
            names(lab),
            col=col, las=2)
        mtext(side=3, at=0, adj=0, main, col="grey30")
    }
    invisible(as.data.frame(out))
}

## linear plots: cf is a coef table (not spp name & taxa as will be in c4i)
## need to provide N/S coef as required
.plot_abundance_lin_2020 <- function(cf, plot=TRUE, veg=TRUE, ylim, main, xlab, ylab,
    link="logit", precalc=FALSE, ...) {
#    FUN <- if (bird)
#        exp else plogis
    FUN <- binomial(link)$linkinv
    if (veg) {
        x0 <- FUN(cf[labs$veg$native,,drop=FALSE])
        xs <- FUN(cf[labs$veg$hfsoft,,drop=FALSE])
        xh <- FUN(cf[labs$veg$hfhard,,drop=FALSE])
    } else {
        x0 <- FUN(cf[labs$soil$native,,drop=FALSE])
        xs <- FUN(cf[labs$soil$hfsoft,,drop=FALSE])
        xh <- FUN(cf[labs$soil$hfhard,,drop=FALSE])
    }
    x0[is.infinite(x0)] <- NA
    xs[is.infinite(xs)] <- NA
    xh[is.infinite(xh)] <- NA
    # apply land cover weights (proportions here if needed)
    tab <- rbind(AverageCoef=colMeans(x0, na.rm=TRUE),
        SoftLin10=0.9*colMeans(x0, na.rm=TRUE)+0.1*colMeans(xs, na.rm=TRUE),
        HardLin10=0.9*colMeans(x0, na.rm=TRUE)+0.1*colMeans(xh, na.rm=TRUE))
    out <- .get_stats(tab, precalc=precalc)

    #USE <- "Median"
    USE <- "First"

    if (plot) {
        p.mean <- out["AverageCoef", USE]
        p.softlin10 <- out["SoftLin10", USE]
        p.hardlin10 <- out["HardLin10", USE]

        if (missing(ylim)) {
            ymax1 <- max(p.softlin10, p.hardlin10, 2*p.mean)*1.03
            ylim <- c(0, ymax1)
        } else {
            ymax <- max(ylim)
        }
#        if (missing(main))
#            main <- species
        if (missing(main))
            main <- ""
        if (missing(xlab))
            xlab <- "Human footprint"
        if (missing(ylab))
            ylab <- "Relative abundance"

        plot(c(1,1.95,2.05), c(p.mean, p.softlin10, p.hardlin10),
            pch=c(1, 16, 15),
            col=c("grey30", "blue3", "red4"),
            xlab=xlab, ylab=ylab, xlim=c(0.8, 2.8), ylim=ylim,
            tck=0.01, yaxs="i", xaxt="n", yaxt="n", bty="l",
            cex=2, lwd=2, cex.lab=1.4, cex.axis=1.3, col.lab="grey40")
        axis(side=2, at=pretty(ylim, n=5), cex.axis=1.3, tck=0.01,
            cex.axis=1.3, col.axis="grey40", col.ticks="grey40")
        axis(side=1, at=c(1,2), labels=c("None","10% linear"),
            tck=0.01, cex.axis=1.3, col.axis="grey40", col.ticks="grey40")
        box(bty="l", col="grey40")
        lines(c(1,1.95), c(p.mean, p.softlin10), col="blue3")
        lines(c(1,2.05), c(p.mean, p.hardlin10), col="red4")
        points(c(1, 1.95, 2.05), c(p.mean, p.softlin10, p.hardlin10),
            pch=c(1,16,15), col=c("grey30", "blue3", "red4"), cex=2, lwd=2)
        lines(c(1,1), c(out["AverageCoef", "Lower"], out["AverageCoef", "Upper"]), col="grey30")
        lines(c(1.95,1.95), c(out["SoftLin10", "Lower"], out["SoftLin10", "Upper"]), col="blue3")
        lines(c(2.05,2.05), c(out["HardLin10", "Lower"], out["HardLin10", "Upper"]), col="red4")
        ly <- c(p.softlin10, p.hardlin10)
        if (abs(ly[2]-ly[1]) < ymax1/20)
            ly <- c(mean(ly)+ymax1/40*sign(ly[1]-ly[2]), mean(ly)+ymax1/40*sign(ly[2]-ly[1]))
        text(c(2.15, 2.15), ly, c("Soft linear","Hard linear"),
            col=c("blue3", "red4"), cex=1.3, adj=0)
        mtext(side=3, at=0.8, adj=0, main, col="grey30", cex=1.3)
    }
    invisible(out)
}


## veg/hf plots: cf is a coef table (not spp name & taxa as will be in c4i)
## coef from the north model is expectes
.plot_abundance_veg_2020 <- function(cf, plot=TRUE, ylim, main, ylab, bw=FALSE,
    link="logit", precalc=FALSE, ...){

    lab <- c(labs$veg$native, labs$veg$hfpoly, labs$veg$hffor)
    # collapse treed fen as average

#    FUN <- if (bird)
#        exp else plogis
    FUN <- binomial(link)$linkinv
    p <- cf[lab,]
    ii <- grep("TreedFen", lab)
    p[ii[1],] <- colMeans(p[ii,])
    p <- p[-ii[-1],]
    lab[ii[1]] <- "TreedFen"
    names(lab)[ii[1]] <- "Treed fen"
    lab <- lab[-ii[-1]]
    rownames(p) <- lab
    p <- .get_stats(FUN(p), precalc=precalc)

    #USE <- "Median"
    USE <- "First"
    out <- as.matrix(p[,c(USE, "Lower", "Upper")])
    colnames(out) <- c("Estimate", "LCL", "UCL")
    rownames(out) <- names(lab)
    out[is.infinite(out)] <- NA

    if (plot) {

        lci <- out[,"LCL"]
        uci <- out[,"UCL"]
        y1 <- out[,"Estimate"]
        names(y1) <- rownames(out)

        if (missing(main))
            main <- ""
        if (missing(ylab))
            ylab <- "Relative abundance"
        if (missing(ylim)) {
            ymax <- min(max(uci, na.rm=TRUE), 2 * max(y1, na.rm=TRUE))
            ylim <- c(0, ymax)
        } else {
            ymax <- max(ylim)
        }

        x <- c(rep(1:9, 5) + rep(seq(0, 40, 10), each=9),
            51, 53, 55, 57, 59, 61, 63, 65, 67,
            70, 72, 74, 76, 78, 80, 82)
        space <- c(1,x[-1]-x[-length(x)])-0.99

        op <- par(mai=c(1.5,1,0.2,0.3))
        on.exit(par(op))

        if (!bw) {
            col.r <- c(rep(0,9), seq(0.3,0.6,length.out=9),
                seq(0.5,1,length.out=9),
                seq(0.8,0.9,length.out=9), rep(0,9), 0, #rep(0,9),
                0.8,0.2,0,0,0,0,0,0,
                0.14,0.29,0.6,0.4,0,0,0)
            col.g <- c(seq(0.5,1,length.out=9), seq(0.4,0.8,length.out=9),
                seq(0.1,0.2,length.out=9), seq(0.4,0.8,length.out=9),
                seq(0.4,0.7,length.out=9), 0.5, #seq(0.15,0.5,length.out=9),
                0.8,0.8,0,0,0,0,0,0,
                0.14,0.29,0.6,0.40,0,0)
            col.b <- c(rep(0,9),rep(0,9),rep(0,9),seq(0.2,0.4,length.out=9),
                seq(0.2,0.6,length.out=9), 0.7, #seq(0.4,0.7,length.out=9),
                0,0,1,1,0,0,0,0,
                0,0,0,0.4,0,0,0)
        } else {
            col.r <- c(rep(seq(0.7,0.2,length.out=9), 5), rep(0.3,15))
            col.b <- col.g <- col.r
        }
        COL <- rgb(col.r, col.g, col.b)
        if (!bw) {
            COL[47:61] <- c(hcl.colors(8,"Dark 3"),
                "#333300", "#666600", "#999900",
                "#909090","#555555","#404040","#222222")

        }
        idx <- 1:length(x)
        x1 <- barplot(y1[idx],
            space=space,
            border="white",
            col=COL,
            ylim=ylim,
            xlim=c(-0.5,max(x)+1.5),
            xaxs="i", yaxt="n",
            ylab=ylab,
            col.lab="grey50",
            cex.lab=1.2, axisnames=FALSE)[,1]
        ax <- axis(side=2, cex.axis=0.9, col.axis="grey50",
            col.ticks="grey50", las=2)
        abline(h=ax, col="grey80")
        x1 <- barplot(y1[idx],
            space=space,
            border="white", col=COL,
            ylim=ylim, xaxs="i", yaxt="n",
            col.lab="grey50",
            cex.lab=1.2, axisnames=FALSE, add=TRUE)[,1]
        box(bty="l", col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[idx][i], y1[idx][i]),col="grey90")
            lines(rep(x1[i],2), c(uci[idx][i], y1[idx][i]),col=COL[i])
        }
        mtext(side=1, at=x1[c(5,14,23,32,41)], line=1.4,
            c("Upland Spruce","Pine","Deciduous","Mixedwood","Black Spruce"),
            col=rgb(col.r[c(5,14,23,32,41)], col.g[c(5,14,23,32,41)],
            col.b[c(5,14,23,32,41)]),las=1)
        at1<-rep(seq(1,9,2),5) + rep(c(0,9,18,27,36), each=5)
        mtext(side=1,at=x1[at1]-0.3,rep(c("0","20","60","100","140"),5),
            line=0.2, adj=0.5, cex=0.8, col=COL[at1])
        mtext(side=1, at=-0.25, adj=1, line=0.2, "Age:", col="grey40", cex=0.8)
        mtext(side=3, at=0, adj=0, main, col="grey30")

        ii <- match(c("TreedFen", "TreedSwamp", "ShrubbySwamp", "ShrubbyBog", "ShrubbyFen",
            "GraminoidFen", "Marsh", "Shrub", "GrassHerb"), lab)
        mtext(side=1, at=x1[ii],
            names(lab)[ii],
            col=COL[ii],
            las=2, adj=1.1)
        ii <- match(c("Crop", "TameP", "RoughP", "Wellsites", "Rural", "Urban", "Industrial"), lab)
        mtext(side=1, at=x1[ii],
            names(lab)[ii],
            col=COL[ii], las=2, adj=1.1)

        ## Add cutblock trajectories - upland conifer
        i1 <- grep("CCWhiteSpruce", lab)
        x2<-x1[1:5]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="White pruce 80-99")
        lines(c(x2[1:5], x1[x3]), y1[c(i1, x3)], col="grey30", lty=2)
        points(x2[1:5], y1[i1], pch=18, cex=1, col="grey30")
        points(x2[1:5], y1[i1], pch=5, cex=0.7, col="grey10")
        ## Pine
        i1 <- grep("CCPine", lab)
        x2<-x1[10:15]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Pine 80-99")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Deciduous
        i1 <- grep("CCDeciduous", lab)
        x2<-x1[19:24]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Deciduous 80-99")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Mixed
        i1 <- grep("CCMixedwood", lab)
        x2<-x1[28:33]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Mixedwood 80-99")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
    }
    invisible(as.data.frame(out))
}

## get nice stats
.get_stats <- function(x, precalc=FALSE) {
    if (precalc) {
        return(
            data.frame(Label=factor(rownames(x), rownames(x)),
                First=x[,1], Mean=x[,1], SD=NA,
                Median=x[,1], Lower=x[,2], Upper=x[,3])
        )
    }
    q <- apply(x, 1, quantile, c(0.5, 0.05, 0.95))
    data.frame(Label=factor(rownames(x), rownames(x)),
        First=x[,1], Mean=rowMeans(x), SD=apply(x, 1, sd),
        Median=q[1,], Lower=q[2,], Upper=q[3,])
}

## make figures
ROOT <- "s:/AB_data_v2020/Results/web"
#ROOT <- "s:/AB_data_v2020/Results/web1"
TAXA <- c("lichens", "mites", "mosses", "vplants", "birds", "mammals", "habitats", "nnplants")
#TAXA <- c( "habitats", "mammals")
COEFS$mammals <- COEFS2$mammals
COEFS$habitats <- COEFS2$habitats
COEFS$nnplants <- COEFS2$nnplants

## store results for bio browser
RESULTS <- list(
    veghf=list(),
    linn=list(),
    lins=list(),
    soilhfnontreed=list(),
    soilhftreed=list()
)
SPECIES <- NULL

for (taxon in TAXA) {

    LINK <- switch(taxon,
        birds="log", mammals="log", nnplants="log", "logit")
    PRECALC <- taxon %in% names(COEFS2)

    if (!dir.exists(file.path(ROOT, taxon)))
        dir.create(file.path(ROOT, taxon))
    ## North
    A <- if (taxon == "birds")
        COEFS[[taxon]]$north$marginal else COEFS[[taxon]]$north
    SPPn <- dimnames(A)[[1]]
    for (spp in SPPn) {
        cat(taxon, spp, "north\n")
        flush.console()

        if (taxon=="habitats")
            LINK <- COEFS[[taxon]]$species[spp, "LinkHabitat"]

        if (!dir.exists(file.path(ROOT, taxon, spp)))
            dir.create(file.path(ROOT, taxon, spp))
        cf <- A[spp,,]
        ## veghf plot
        png(file.path(ROOT, taxon, spp, "veghf.png"), width=1000, height=500)
#        svglite(file.path(ROOT, taxon, spp, "veghf.svg"), width=10, height=5)
        vhf <- .plot_abundance_veg_2020(cf,main=spp, link=LINK, precalc=PRECALC)
        dev.off()
        RESULTS$veghf[[spp]] <- vhf
        ## linear N plot
        png(file.path(ROOT, taxon, spp, "lin-north.png"), width=500, height=500)
#        svglite(file.path(ROOT, taxon, spp, "lin-north.svg"), width=5, height=5)
        linn <- .plot_abundance_lin_2020(cf, main=spp, veg=TRUE, link=LINK, precalc=PRECALC)
        dev.off()
        RESULTS$linn[[spp]] <- linn
    }

    ## south
    A <- if (taxon == "birds")
        COEFS[[taxon]]$south$marginal else COEFS[[taxon]]$south
    SPPs <- dimnames(A)[[1]]
    for (spp in SPPs) {
        cat(taxon, spp, "south\n")
        flush.console()

        if (taxon=="habitats")
            LINK <- COEFS[[taxon]]$species[spp, "LinkHabitat"]

        if (!dir.exists(file.path(ROOT, taxon, spp)))
            dir.create(file.path(ROOT, taxon, spp))
        cf <- A[spp,,]
        ## veghf plot
        shf0 <- .plot_abundance_soil_2020(cf,main=spp, paspen=0, link=LINK, precalc=PRECALC, plot=FALSE)
        shf1 <- .plot_abundance_soil_2020(cf,main=spp, paspen=1, link=LINK, precalc=PRECALC, plot=FALSE)
        MAX <- max(shf0, shf1, na.rm=TRUE)
        MAXest <- max(shf0$Estimate, shf1$Estimate, na.rm=TRUE)
        if (MAX > 2.5*MAXest)
            MAX <- 2.5*MAXest
        png(file.path(ROOT, taxon, spp, "soilhf.png"), width=1000, height=500)
        op <- par(mfrow=c(1,2))
        .plot_abundance_soil_2020(cf,main=paste(spp, "non-treed"),
            paspen=0, link=LINK, precalc=PRECALC, plot=TRUE, ylim=c(0, MAX))
        .plot_abundance_soil_2020(cf,main=paste(spp, "treed"),
            paspen=1, link=LINK, precalc=PRECALC, plot=TRUE, ylim=c(0, MAX))
        par(op)
        dev.off()
        RESULTS$soilhfnontreed[[spp]] <- shf0
        RESULTS$soilhftreed[[spp]] <- shf1
        ## linear N plot
        png(file.path(ROOT, taxon, spp, "lin-south.png"), width=500, height=500)
        lins <- .plot_abundance_lin_2020(cf, main=spp, veg=FALSE, link=LINK, precalc=PRECALC)
        dev.off()
        RESULTS$lins[[spp]] <- lins
    }

    SPP <- data.frame(sppid=sort(unique(c(SPPn, SPPs))), taxon=taxon)
    SPP$north <- SPP$sppid %in% SPPn
    SPP$south <- SPP$sppid %in% SPPs
    SPECIES <- rbind(SPECIES, SPP)

}

save(RESULTS, SPECIES, file="s:/AB_data_v2020/Results/BB-ESTIMATES.RData")

