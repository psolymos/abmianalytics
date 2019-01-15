#' ---
#' title: "Summarizing north and south model results"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Dec 19, 2018"
#' output: pdf_document
#' ---
#'
library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
#ROOT <- "~/GoogleWork/tmp"

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)
TAX <- ee$tax
rm(ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]


spp <- "BTNW"
resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))

estn <- get_coef(resn, Xn, stage="ARU", na.out=FALSE)
printCoefmat(get_summary(estn))

explore_north <- function(resn, plot=TRUE, ...) {
    estn <- get_coef(resn, Xn, stage="ARU", na.out=FALSE)
    mu <- Xage %*% t(estn[,colnames(Xage)])
    lam1 <- t(apply(exp(mu), 1, quantile, c(0.5, 0.05, 0.95)))
    lam1 <- lam1[!grepl("9", rownames(lam1)),]
    lamCC <- lam1[grepl("CC", rownames(lam1)),]
    lam <- lam1[!grepl("CC", rownames(lam1)),]

    MOD <- c("ROAD", "mWell", "mSoft",
        "mEnSft", "mTrSft", "mSeism", "CMETHODSM", "CMETHODRF")
    lamMOD <- t(apply(exp(estn[,MOD]), 2, quantile, c(0.5, 0.05, 0.95)))
    rownames(lamMOD) <- c("Road", "Well", "Soft", "EnSoft", "TrSoft", "Seismic", "SM", "RiverFork")

    if (plot) {
        layout(matrix(c(1,1,1,1,1,1,2,2,3), 3, 3, byrow = TRUE))
        op <- par(las=2, mar=c(10,4,2,2), cex.axis=0.9)
        on.exit(par(op))
        col <- c(rep(RColorBrewer::brewer.pal(8, "Accent")[c(1,2,3,5,6,7)], each=9),
            RColorBrewer::brewer.pal(12, "Set3"))
        k <- 2
        NAM <- as.character(TAX[attr(resn, "spp"), "species"])
        b <- barplot(unname(lam[,1]), col=col,
            ylab="Relative abundance", main=paste0(NAM, ", North"),
            ylim=c(0, min(k*max(lam1[,1]), max(lam1))))
        mtext(rownames(lam), col=col, side=1, at=b, cex=0.8, line=1)
        o <- b[2]-b[1]
        for (i in 1:4) {
            xi <- b[(((i-1)*9+1):(i*9))][1:6]
            yi <- rbind(lamCC[((i-1)*5+1):(i*5),], lam[(((i-1)*9+1):(i*9))[6],,drop=FALSE])
            lines(xi-0.2*o, yi[,1], lwd=1, lty=1, col="darkgrey")
            points(xi[-6]-0.2*o, yi[-6,1], col=1, lwd=2, cex=0.6)
            segments(x0=xi[-6]-0.2*o, y0=yi[-6,2], y1=yi[-6,3], lwd=1, col=1)
        }
        segments(x0=b, y0=lam[,2], y1=lam[,3], lwd=1, col=1)
        par(las=1, mar=c(4,4,2,2))
        b <- barplot(lamMOD[1:6,1], ylim=c(0, min(k*max(lamMOD[1:6,1]), max(lamMOD[1:6,]))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(8, "Accent")[1:6], main="Modifiers, North")
        segments(x0=b, y0=lamMOD[1:6,2], y1=lamMOD[1:6,3], lwd=2,
            col=ifelse(1 %[]% lamMOD[1:6,2:3], 1, 2))
        abline(h=1, col=2)

        b <- barplot(lamMOD[7:8,1], ylim=c(0, min(k*max(lamMOD[7:8,1]), max(lamMOD[7:8,]))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(8, "Accent")[7:8], main="ARU, North")
        segments(x0=b, y0=lamMOD[7:8,2], y1=lamMOD[7:8,3], lwd=2,
            col=ifelse(1 %[]% lamMOD[7:8,2:3], 1, 2))
        abline(h=1, col=2)
    }
    invisible(list(lam=lam1, mod=lamMOD))
}

pdf(file.path(ROOT, "explore-north.pdf"), onefile=TRUE, width=10, height=10)
expln <- list()
for (i in colnames(en$YY)) {
    cat(i, "\n");flush.console()
    expln[[i]] <- explore_north(load_species(file.path(ROOT, "out", "north", paste0(i, ".RData"))))
}
dev.off()


#' Explore south

explore_south <- function(ress, plot=TRUE, ...) {
    ests <- get_coef(ress, Xs, stage="ARU", na.out=FALSE)
    LCC <- c("soilcClay", "soilcCrop", "soilcRapidDrain",
        "soilcRoughP", "soilcSaline", "soilcTameP", "soilcUrbInd")
    muLCC <- t(cbind(ests[,1], ests[,1]+ests[,LCC]))
    rownames(muLCC) <- levels(es$DAT$soilc)
    qfun <- function(x)
        c(quantile(x, c(0.5, 0.05, 0.95)), "1st"=x[1])
    lamLCC0 <- t(apply(exp(muLCC), 1, qfun)) # nontred
    lamLCC1 <- t(apply(exp(muLCC + ests[,"pAspen"]), 1, qfun)) # treed
    # RoughP: Abandoned and RoughPasture (closer to natural grass)
    # TameP: TamePasture + HD livestock operations (closer to crop)
    o <- c("Productive", "Clay", "Saline", "RapidDrain", "RoughP", "TameP", "Crop", "UrbInd")
    lamLCC0 <- lamLCC0[o,]
    lamLCC1 <- lamLCC1[o,]
    MOD <- c("ROAD", "mWell", "mSoft", "CMETHODSM", "CMETHODRF")
    lamMOD <- t(apply(exp(ests[,MOD]), 2, qfun))
    rownames(lamMOD) <- c("Road", "Well", "Soft", "SM", "RiverFork")

    if (plot) {
        k <- 2
        NAM <- as.character(TAX[attr(ress, "spp"), "species"])
        layout(matrix(c(1,2,1,2,3,4), 3, 2, byrow = TRUE))
        op <- par(mar=c(10,4,2,1), las=2, cex.axis=0.9)
        on.exit(par(op))
        b <- barplot(lamLCC0[,1], ylim=c(0, min(k*max(lamLCC0[,1], lamLCC1[,1]), max(lamLCC0, lamLCC1))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(8, "Dark2"), main=paste0(NAM, "\nNon-treed, South"), ...)
        segments(x0=b, y0=lamLCC0[,2], y1=lamLCC0[,3], lwd=2, col=1)
        #points(b, lamLCC0[,4])

        b <- barplot(lamLCC1[,1], ylim=c(0, min(k*max(lamLCC0[,1], lamLCC1[,1]), max(lamLCC0, lamLCC1))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(8, "Dark2"), main="Treed, South", ...)
        segments(x0=b, y0=lamLCC1[,2], y1=lamLCC1[,3], lwd=2, col=1)
        #points(b, lamLCC1[,4])

        par(mar=c(4,4,2,1), las=1)
        b <- barplot(lamMOD[1:3,1], ylim=c(0, min(k*max(lamMOD[1:3,1]), max(lamMOD[1:3,]))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(5, "Accent")[1:3], main="Modifiers, South", ...)
        segments(x0=b, y0=lamMOD[1:3,2], y1=lamMOD[1:3,3], lwd=2,
            col=ifelse(1 %[]% lamMOD[1:3,2:3], 1, 2))
        #points(b, lamMOD[1:3,4])
        abline(h=1, col=2)

        b <- barplot(lamMOD[4:5,1], ylim=c(0, min(k*max(lamMOD[4:5,1]), max(lamMOD[4:5,]))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(5, "Accent")[4:5], main="ARU, South", ...)
        segments(x0=b, y0=lamMOD[4:5,2], y1=lamMOD[4:5,3], lwd=2,
            col=ifelse(1 %[]% lamMOD[4:5,2:3], 1, 2))
        #points(b, lamMOD[4:5,4])
        abline(h=1, col=2)
    }
    invisible(list(lcc0=lamLCC0, lcc1=lamLCC1, mod=lamMOD))
}

pdf(file.path(ROOT, "explore-south.pdf"), onefile=TRUE, width=10)
expls <- list()
for (i in colnames(es$YY)) {
    cat(i, "\n");flush.console()
    expls[[i]] <- explore_south(load_species(file.path(ROOT, "out", "south", paste0(i, ".RData"))))
}
dev.off()

