#' ---
#' title: "Summarizing north and south model results"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Sept, 2020"
#' output: pdf_document
#' ---
#'
library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds" # change this bit
#ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
#ROOT <- "~/GoogleWork/tmp"

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)
#load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)
TAX <- ee$tax
rm(ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2020-09-23.RData"), envir=en)
#load(file.path(ROOT, "data", "ab-birds-north-2019-01-30.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2020-09-23.RData"), envir=es)
#load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

if (FALSE) {
    ## arbitrary species data
    spp <- "LISP"
    xn=en$DAT
    xn[[paste0(spp, "_count")]] <- as.numeric(en$YY[,spp])
    xn[[paste0(spp, "_offset")]] <- as.numeric(en$OFF[,spp])
    xs=es$DAT
    xs[[paste0(spp, "_count")]] <- as.numeric(es$YY[,spp])
    xs[[paste0(spp, "_offset")]] <- as.numeric(es$OFF[,spp])
    write.csv(xn, row.names = FALSE, file="~/Desktop/x/LISP-ABMI-North-2020.csv")
    write.csv(xs, row.names = FALSE, file="~/Desktop/x/LISP-ABMI-South-2020.csv")

    xxn <- data.frame(colname=colnames(xn), description="")
    xxs <- data.frame(colname=colnames(xs), description="")
    write.csv(xxn, row.names = FALSE, file="~/Desktop/x/Metadata-ABMI-North-2020.csv")
    write.csv(xxs, row.names = FALSE, file="~/Desktop/x/Metadata-ABMI-South-2020.csv")
}


#Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v2020.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

if (FALSE) {
    for (spp in colnames(en$OFF)) {
        resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))

        estn1 <- suppressWarnings(get_coef(resn, Xn, stage="ARU", na.out=FALSE))
        write.csv(estn1, file=paste0("s:/_tmp/birds-north-aru/", TAX[spp,"sppid"], "_coef.csv"))

        estn2 <- suppressWarnings(get_coef(resn, Xn, stage="Space", na.out=FALSE))
        write.csv(estn2, file=paste0("s:/_tmp/birds-north-space/", TAX[spp,"sppid"], "_coef.csv"))

        estn3 <- suppressWarnings(get_coef(resn, Xn, stage="HF", na.out=FALSE))
        #write.csv(estn3, file=paste0("_tmp/", spp, ".HFmodelcoefficients.csv"))
        write.csv(estn3, file=paste0("s:/_tmp/birds-north-hf/", TAX[spp,"sppid"], "_coef.csv"))
    }

}

spp <- "ALFL"
resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))

if (FALSE) {
    ## SSH matrix for 1st run
    m <- as.matrix(en$SSH[1:ncol(en$YY),])
    rownames(m) <- colnames(en$YY)
    m[] <- 0
    for (spp in colnames(en$YY)) {
        resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
        try(m[spp, resn[[1]]$ssh$labels] <- 1)
    }
    write.csv(m, file=file.path(ROOT, "ssh-north.csv"))
}

## --

estn <- get_coef(resn, Xn, stage="ARU", na.out=FALSE)
printCoefmat(get_summary(estn))

b <- estn[,"mWell"]
p <- 1/7
m <- exp(p*b)
z <- (m-(1-p))/p

p <- seq(0,1,0.01)
pm <- 1/7
mm <- exp(pm*b[1])
plot(p, exp(b[1]*p), type="l")
abline(v=pm, lty=2)
## need to keep m=1 at p=0: exp(0*b)=1
abline(1,(mm-1)/pm,col=2)

plot(p, exp(b[1]*p), type="l")
abline(v=pm, lty=2)
curve(linexp(x, b[1], pm), add=TRUE, col=4)

bb <- unique(sort(as.numeric(en$BB)))

explore_modifs <- function(resn, plot=TRUE, ...) {
    OK <- !sapply(resn, inherits, "try-error")
    spp <- resn[[which(OK)[1]]]$species
    y <- en$YY[bb,spp]
    npk1 <- sum(y>0)
    npkT <- length(y)
    yss <- sum_by(y, en$DAT$SS[bb])[,"x"]
    nss1 <- sum(yss>0)
    nssT <- length(yss)

    estn <- suppressWarnings(get_coef(resn, Xn, stage="ARU", na.out=FALSE))
    mu <- Xage %*% t(estn[,colnames(Xage)])
    lam1 <- t(apply(exp(mu), 1, quantile, c(0.5, 0.05, 0.95)))
    lam1 <- lam1[!grepl("9", rownames(lam1)),]
    lamCC <- lam1[grepl("CC", rownames(lam1)),]

    MOD <- c("ROAD", "mWell", "mSoft",
        "mEnSft", "mTrSft", "mSeism", "CMETHODSM", "CMETHODRF")
    Z <- exp(estn[,MOD])
    isSoft <- estn[,"mSoft"] != 0 & estn[,"mEnSft"] == 0
    #isSoft2 <- get_mid(resn)[,"Contrast"] == 3
    estn[isSoft,"mEnSft"] <- estn[isSoft,"mSoft"]
    estn[isSoft,"mTrSft"] <- estn[isSoft,"mSoft"]
    estn[isSoft,"mSeism"] <- estn[isSoft,"mSoft"]
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
        "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05,
        "CMETHODSM"=1, "CMETHODRF"=1)
    for (i in MOD)
        Z[,i] <- linexp(1, estn[,i], pm[i])

    HFc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")

    Xn2 <- Xn[en$DAT$mWell > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 > 0, colnames(Xage)]
    lamWell <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2)])), 2, median)
    estWell <- quantile(lamWell * Z[,"mWell"], c(0.5, 0.05, 0.95))

    Xn2 <- Xn[en$DAT$mEnSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 > 0, colnames(Xage)]
    lamEnSft <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2)])), 2, median)
    estEnSft <- quantile(lamEnSft * Z[,"mEnSft"], c(0.5, 0.05, 0.95))

    ## TrSft incorporates ROAD effect as well?
    Xn2 <- Xn[en$DAT$mTrSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 > 0, colnames(Xage)]
    lamTrSft <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2)])), 2, median)
    #estTrSft <- quantile(lamTrSft * Z[,"mTrSft"] * Z[,"ROAD"], c(0.5, 0.05, 0.95))
    estTrSft <- quantile(lamTrSft * Z[,"mTrSft"], c(0.5, 0.05, 0.95))

    ESMAX <- apply(rbind(lamCC[endsWith(rownames(lamCC), "R"),], EnSft=estEnSft, TrSft=estTrSft), 2, max)
    Xn2 <- Xn[en$DAT$mSeism > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 > 0, colnames(Xage)]
    lamSeism <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2)])), 2, median)
    estSeism <- quantile(lamSeism * Z[,"mSeism"], c(0.5, 0.05, 0.95))
    if (estSeism[1] > ESMAX[1])
        estSeism <- ESMAX

    lam1 <- rbind(lam1, Well=estWell, EnSft=estEnSft, TrSft=estTrSft, Seism=estSeism)
    lam <- lam1[!grepl("CC", rownames(lam1)),]

    if (plot) {
        op <- par(las=2, mar=c(10,4,3,2), cex.axis=0.9)
        on.exit(par(op))
        col <- c(rep(RColorBrewer::brewer.pal(8, "Accent")[c(1,2,3,5,6,7)], each=9),
            RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(4, "Dark2"))
        k <- 2
        NAM <- as.character(TAX[attr(resn, "spp"), "species"])
        b <- barplot(unname(lam[,1]), col=col,
            ylab="Relative abundance",
            main=paste0(NAM, ", North\n", npk1,"/",npkT, " pts, ", nss1, "/", nssT, " loc"),
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
    }
    invisible(lam)
}

pdf(file.path(ROOT, "misc", "explore-modif.pdf"), onefile=TRUE, width=10, height=8)
expln <- list()
for (i in colnames(en$YY)) {
    cat(i, "\n");flush.console()
    expln[[i]] <- explore_modifs(load_species(file.path(ROOT, "out", "north", paste0(i, ".RData"))))
}
dev.off()

explore_north <- function(resn, plot=TRUE, lin=TRUE, ...) {
    estn <- suppressWarnings(get_coef(resn, Xn, stage="ARU", na.out=FALSE))
    mu <- Xage %*% t(estn[,colnames(Xage)])
    qfun <- function(x)
        c(quantile(x, c(0.5, 0.05, 0.95)), "1st"=x[1])
    lam1 <- t(apply(exp(mu), 1, quantile, c(0.5, 0.05, 0.95)))
    lam1 <- lam1[!grepl("9", rownames(lam1)),]
    lamCC <- lam1[grepl("CC", rownames(lam1)),]
    lam <- lam1[!grepl("CC", rownames(lam1)),]

    MOD <- c("ROAD", "mWell", "mSoft",
        "mEnSft", "mTrSft", "mSeism", "CMETHODSM", "CMETHODRF")
    Z <- exp(estn[,MOD])
    isSoft <- estn[,"mSoft"] != 0 & estn[,"mEnSft"] == 0
    #isSoft2 <- get_mid(resn)[,"Contrast"] == 3
    estn[isSoft,"mEnSft"] <- estn[isSoft,"mSoft"]
    estn[isSoft,"mTrSft"] <- estn[isSoft,"mSoft"]
    estn[isSoft,"mSeism"] <- estn[isSoft,"mSoft"]
    if (lin) {
        pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
            "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05,
            "CMETHODSM"=1, "CMETHODRF"=1)
        for (i in MOD)
            Z[,i] <- linexp(1, estn[,i], pm[i])
    }
    lamMOD <- t(apply(Z, 2, qfun))
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
        ii <- c(1,2,4,5,6)
        b <- barplot(lamMOD[ii,1], ylim=c(0, min(k*max(lamMOD[ii,1]), max(lamMOD[ii,]))),
            ylab="Relative abundance",
            col=RColorBrewer::brewer.pal(8, "Accent")[ii], main="Modifiers, North")
        segments(x0=b, y0=lamMOD[ii,2], y1=lamMOD[ii,3], lwd=2,
            col=ifelse(1 %[]% lamMOD[ii,2:3], 1, 2))
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

pdf(file.path(ROOT, "misc", "explore-north.pdf"), onefile=TRUE, width=10, height=10)
expln <- list()
for (i in colnames(en$YY)) {
    cat(i, "\n");flush.console()
    expln[[i]] <- explore_north(load_species(file.path(ROOT, "out", "north", paste0(i, ".RData"))))
}
dev.off()


#' Explore south

explore_south <- function(ress, plot=TRUE, lin=TRUE, ...) {
    ests <- suppressWarnings(get_coef(ress, Xs, stage="ARU", na.out=FALSE))
    #LCC <- c("soilcClay", "soilcCrop", "soilcRapidDrain",
    #    "soilcRoughP", "soilcSaline", "soilcTameP", "soilcUrbInd")
    LCC <- c("soilcBlowout", "soilcClaySub", "soilcCrop",
        "soilcIndustrial", "soilcMine", "soilcOther", "soilcRapidDrain",
        "soilcRoughP", "soilcRural", "soilcSandyLoam", "soilcTameP",
        "soilcThinBreak", "soilcUrban")
    muLCC <- t(cbind(ests[,1], ests[,1]+ests[,LCC]))
    rownames(muLCC) <- levels(es$DAT$soilc)
    qfun <- function(x) {
        c(quantile(x, c(0.5, 0.05, 0.95)), "1st"=x[1])
    }
    lamLCC0 <- t(apply(exp(muLCC), 1, qfun)) # nontred
    lamLCC1 <- t(apply(exp(muLCC + ests[,"pAspen"]), 1, qfun)) # treed
    # RoughP: Abandoned and RoughPasture (closer to natural grass)
    # TameP: TamePasture + HD livestock operations (closer to crop)
    #o <- c("Productive", "Clay", "Saline", "RapidDrain", "RoughP", "TameP", "Crop", "UrbInd")
    o <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other",
        "RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")
    lamLCC0 <- lamLCC0[o,]
    lamLCC1 <- lamLCC1[o,]
    MOD <- c("ROAD", "mWell", "mSoft", "CMETHODSM", "CMETHODRF")
    Z <- exp(ests[,MOD])
    if (lin) {
        pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
            "CMETHODSM"=1, "CMETHODRF"=1)
        for (i in MOD)
            Z[,i] <- linexp(1, ests[,i], pm[i])
    }
    lamMOD <- t(apply(Z, 2, qfun))
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


ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(ts) <- ts[,1]
levels(ts$UseInAnalysis)

pdf(file.path(ROOT, "misc", "explore-south.pdf"), onefile=TRUE, width=10)
expls <- list()
for (i in colnames(es$YY)) {
    cat(i, "\n");flush.console()
    expls[[i]] <- explore_south(load_species(file.path(ROOT, "out", "south", paste0(i, ".RData"))))
}
dev.off()

spptab <- read.csv("~/repos/abmispecies/_data/birds.csv")
spptab <- droplevels(spptab[spptab$modelS,])
rownames(spptab) <- spptab$AOU
SPP <- rownames(spptab)
cf0 <- read.csv("d:/abmi/sppweb2017/tables/AllTaxaCombined-SoilNontreedSouth.csv")
Cult0 <- cf0[match(spptab$species, cf0$Species) ,c("Cult", "Cult.LCL", "Cult.UCL")]

expls <- list()
for (i in SPP) {
    cat(i, "\n");flush.console()
    expls[[i]] <- explore_south(load_species(file.path(ROOT, "out", "south", paste0(i, ".RData"))), plot=FALSE)
}

Crop <- t(sapply(expls, function(z) z$lcc0["Crop",]))[,1:3]
Rough <- t(sapply(expls, function(z) z$lcc0["RoughP",]))[,1:3]
Tame <- t(sapply(expls, function(z) z$lcc0["TameP",]))[,1:3]
Crop <- t(sapply(expls, function(z) z$lcc0["Crop",]))[,1:3]
Prod <- t(sapply(expls, function(z) z$lcc0["Productive",]))[,1:3]

plot(density(Cult0[,3]-Cult0[,2]))
lines(density(Crop[,3]-Crop[,2]), col=2)
lines(density(Tame[,3]-Tame[,2]), col=3)
lines(density(Rough[,3]-Rough[,2]), col=4)

o <- order(Cult0[,1], decreasing=TRUE)
plot(1:nrow(Cult0)-0.15, Cult0[o,1], type="l", ylim=c(0,max(Cult0, Crop, Tame, Rough, Prod)),
    xlim=c(1,10), axes=FALSE, ann=FALSE)
title(ylab="Relative abundance")
axis(2)
axis(1, 1:10, rownames(spptab[o,])[1:10], las=2, tick=FALSE)

segments(x0=1:nrow(Cult0)-0.15, y0=Cult0[o,2], y1=Cult0[o,3])
points(1:nrow(Cult0)-0.15, Cult0[o,1], pch=19)
points(1:nrow(Cult0)-0.05, Crop[o,1], col=2, pch=19)
segments(x0=1:nrow(Cult0)-0.05, y0=Crop[o,2], y1=Crop[o,3], col=2)
points(1:nrow(Cult0)+0.05, Tame[o,1], col=3, pch=19)
segments(x0=1:nrow(Cult0)+0.05, y0=Tame[o,2], y1=Tame[o,3], col=3)
points(1:nrow(Cult0)+0.15, Rough[o,1], col=4, pch=19)
segments(x0=1:nrow(Cult0)+0.15, y0=Rough[o,2], y1=Rough[o,3], col=4)
points(1:nrow(Cult0)+0.25, Prod[o,1], col="tan", pch=19)
segments(x0=1:nrow(Cult0)+0.25, y0=Prod[o,2], y1=Prod[o,3], col="tan")
legend("topright", pch=19, col=c(1:4, "tan"), lty=1, bty="n",
    legend=c("Cultivation (old)", "Crop", "Tame", "Rough", "Productive"))


round(cor(cbind(Crop=Crop[,1], Tame=Tame[,2], Rough=Rough[,1], Productive=Prod[,1])),2)


plot(density(Cult0[,1]))
lines(density(Crop[,1]), col=2)
lines(density(Tame[,1]), col=3)
lines(density(Rough[,1]), col=4)


## compare dominant type vs compositio in the south (no age complexity)

e <- new.env()
load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData", envir=e)
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]

sc1r <- row_std(groupSums(e$sc1[rownames(Xs),], 2, ts[colnames(e$sc1),"UseInAnalysisCoarse"]))

compare_sets(es$DAT$soilc, colnames(sc1r))

ps <- sc1r[,levels(es$DAT$soilc)]

DAT <- es$DAT
DAT$pProductive <- ps[,"Productive"]
DAT$pClay <- ps[,"Clay"]
DAT$pCrop <- ps[,"Crop"]
DAT$pRapidDrain <- ps[,"RapidDrain"]
DAT$pRoughP <- ps[,"RoughP"]
DAT$pSaline <- ps[,"Saline"]
DAT$pTameP <- ps[,"TameP"]
DAT$pUrbInd <- ps[,"UrbInd"]

mods1 <- mods2 <- es$mods
mods2$Hab <- list(
    .~. + pProductive + pClay + pCrop + pRapidDrain + pRoughP + pSaline + pTameP + pUrbInd,
    .~. + pProductive + pClay + pCrop + pRapidDrain + pRoughP + pSaline + pTameP + pUrbInd + pAspen)
k <- 3


#spp <- "ALFL"
dc <- list()
for (spp in colnames(es$OFF)) {
    cat(spp, "\n");flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    mid <- get_mid(ress)
    pieces1 <- lapply(1:k, function(i) if (mid[1,i] < 1) .~. else mods1[[i]][[mid[1,i]]])
    pieces2 <- lapply(1:k, function(i) if (mid[1,i] < 1) .~. else mods2[[i]][[mid[1,i]]])
    ff1 <- ff2 <- y ~ 1
    for (i in 1:k) {
        ff1 <- update(ff1, pieces1[[i]])
        ff2 <- update(ff2, pieces2[[i]])
    }

    dat <- DAT[es$BB[,1],]
    y <- es$YY[es$BB[,1],spp]
    off <- if (spp %in% colnames(es$OFF))
        es$OFF[es$BB[,1],spp] else es$OFFmean[es$BB[,1]]

    m1 <- glm(ff1, data=dat, offset=off, family=poisson)
    m2 <- glm(ff2, data=dat, offset=off, family=poisson)

    dc[[spp]] <- list(dom=glm_skeleton(m1), com=glm_skeleton(m2))
}


conv1 <- sapply(dc, function(z) z$dom$converge)
conv2 <- sapply(dc, function(z) z$com$converge)

aic1 <- sapply(dc, function(z) z$dom$aic)
aic2 <- sapply(dc, function(z) z$com$aic)

table(ComConverged=conv2, DomBetter=aic1-aic2 < 0)


f <- function(x, dom=TRUE) {
    m1 <- x$dom
    m2 <- x$com
    cf1 <- c(coef(m1)[1], coef(m1)[1] + coef(m1)[2:8])
    cf2 <- coef(m2)[1] + coef(m2)[2:9]
    names(cf1) <- names(cf2) <- colnames(ps)
    #cbind(Dominant=exp(cf1), Composition=exp(cf2))
    if (dom) exp(cf1) else exp(cf2)
}

cf1 <- sapply(dc, f, dom=TRUE)
cf2 <- sapply(dc, f, dom=FALSE)

plot(as.numeric(cf1[,conv2]), as.numeric(cf2[,conv2]),
    xlim=c(0,1), ylim=c(0,1))
abline(0,1)

round(data.frame(Dominant=exp(cf1), Composition=exp(cf2)),4)

plot(exp(cf1), exp(cf2))
abline(0,1)

## Bird result checking and revision

library(cure4insect)
load_common_data()
x <- read.csv("~/Dropbox/Public/birds-lookup.csv")
x$v2020 <- paste0(ifelse(x$coefnorth=="+", "N", ""), ifelse(x$coefsouth=="+", "S", ""))
x$v2020[x$v2020=="NS"] <- "N+S"
table(x$v2020)
s <- get_species_table("birds")
s$v2018 <- s$model_region
levels(s$v2018) <- c("N", "N+S", "S", "")
s <- s[match(x$id, s$SpeciesID),]
x$v2018 <- s$v2018
x$v2018[is.na(x$v2018)] <- ""

v <- read.csv("~/repos/abmispecies/_data/birds.csv")
v <- v[match(x$id, v$sppid),]
x$comments <- v$comments
x$comments[is.na(x$comments)] <- ""

addmargins(with(x, table(v2018, v2020)))
nc <- nchar(x$comments)
x$remarks <- substr(x$comments, 1, pmin(nc, 25))

xx <- x[,c("common", "order", "family", "v2018", "v2020", "remarks")]
xx <- xx[order(xx$order, xx$family, xx$common),]
write.csv(xx, row.names = FALSE, file="~/Dropbox/Public/birds-lookup-2.csv")
