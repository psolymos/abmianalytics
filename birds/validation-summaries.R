#' ---
#' title: "Summarizing validation results"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Dec 19, 2018"
#' output: pdf_document
#' ---
#'
#' # Preamble
#'
#' Load packages and functions, see repo on
#' [GitHub](https://github.com/psolymos/abmianalytics)
library(mefa4)
source("~/repos/abmianalytics/birds/00-functions.R")
#'
#' `ROOT` refers to the folder where data and results are
#' - `<ROOT>/data` hosts the data package
#' - `<ROOT>/out/<PROJ>` hosts the species specific output files
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
PROJ <- "validation"
#' Load data file that has the following objects:
#' rows of `DAT`, `YY`, `OFF`, and `SSH` match
#' (same with `DAT`, `YY`, `OFF`, and `SSH` for validation set
#' including 500 ABMI sites with 9 points and big grids)
#' - `DAT`: data frame with point count info (SS, PCODE, PKET etc) and predictors
#' - `YY`: sparse matrix with PKEY x Species data
#' - `OFF`: offsets in a PKEY x Species matrix
#' - `SSH`: surrounding suitable habitat composition (PKEY x land cover)
#' - `BB`: bootstrap PKEY index matrix
#' - `mods`: list of model terms for branching model selection
#' - `DATv`: predictors for validation data
#' - `YYv`: species matrix for validation data
#' - `OFFv`: offsets for validation data
#' - `SSHv`: surrounding suitable habitat for validation data
load(file.path(ROOT, "data", "ab-birds-biggrids-2019-05-13.RData"))
#' Make a model matrix that matches the coefficients that we estimated
#' `X` is for training data, `Xv` is for validation data set
X <- get_model_matrix(DATv, mods)
Xv <- get_model_matrix(DATv, mods)
#' We decided to use the more complex model and not to do all the stages
Stage <- "HF"
#'
#' # Aggregating up to different scales
#'
#' Functions for the magic: uses only external data
#' Groups argument defines the spatial units to aggregate to
#' must match DATv row ordering
#' unwanted rows should be blank (`""`) or `NA`
validate <- function(res, Groups, stage="HF") {
    Groups <- as.character(Groups)
    Groups[is.na(Groups)] <- ""
    ss <- Groups != ""
    muo <- predict_with_SSH(res, Xv[ss,,drop=FALSE], SSHv[ss,,drop=FALSE], stage=stage)
    offo <- OFFv[ss,res[[1]]$species]
    evo <- apply(exp(muo+offo), 1, median)
    yo <- as.numeric(YYv[ss,res[[1]]$species])
    Groups <- as.factor(Groups[ss])

    Predm <- groupSums(exp(muo+offo), 1, Groups)
    Y <- sum_by(yo, Groups)
    yobs <- Y[,"x"]
    npool <- Y[,"by"]
    ypred <- apply(Predm, 1, median)
    m1 <- get_mass(yo, evo)
    m <- get_mass(yobs, ypred)

    AUCo <- simple_auc(simple_roc(ifelse(yo>0, 1, 0), evo))
    CORo <- cor(yobs, ypred, method="spearman")
    CORUo <- .cor_uncentered(yobs, ypred)
    list(spp=res[[1]]$species,
        Groups=Groups,
        AUC=AUCo, COR=CORo, CORU=CORUo,
        y1obs=yo, lam1=evo, yobs=yobs, lam=ypred,
        cmf1=m1, cmf=m, npool=npool,
        mu=muo, off=offo)
}
#' These functions do the randomization
.randomize <- function(v) {
    g <- sample(v$Groups)
    Predm <- groupSums(exp(v$mu+v$off), 1, g)
    Y <- sum_by(v$y1obs, g)
    yobs <- Y[,"x"]
    ypred <- apply(Predm, 1, median)
    CORo <- cor(yobs, ypred, method="spearman")
    CORUo <- .cor_uncentered(yobs, ypred)
    c(COR=CORo, CORU=CORUo)
}
randomize <- function(v, B=199) {
    cbind(c(COR=v$COR, CORU=v$CORU), replicate(B, .randomize(v)))
}
summarize <- function(r) {
    q <- t(apply(r, 1, quantile, c(0.025, 0.5, 0.975)))
    p <- c(pmc(r[1,1], r[1,]), pmc(r[2,1], r[2,]))
    cbind(Est=r[,1], q, "p-value"=p)
}
#'
#' Take species with enough detections and run all grouping variables
#' with randomization where it makes sense
ALL <- list()
SPP <- colnames(YYv[,colSums(YYv>0)>100])
#spp <- "WTSP"
for (spp in SPP) {
    cat("\n", spp);flush.console()

    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))

    VV <- list()
    for (g in c("EX2x2", "EX3x3", "EX4x4", "EX5x5", "EX10x10", "GR2x2", "GR3x3", "GR4x4", "GR5x5", "GR10x10")) {
        V <- validate(res, Groups=as.factor(DATv[,g]), stage="HF")
        VV[[g]] <- randomize(V)
    }
    for (g in c("TS1", "TS2", "TS3", "TS4", "TR1", "TR2", "TR3", "TR4")) {
        v <- validate(res, Groups=as.factor(DATv[,g]), stage="HF")
        VV[[g]] <- cbind(c(COR=v$COR, CORU=v$CORU))
    }
    ALL[[spp]] <- VV
}

save(ALL, file=file.path(ROOT, "validation-BG-results-2019-05-13.RData"))


pdf(file.path(ROOT, "validation-BG-results-2019-05-13.pdf"), onefile = TRUE, width=15, height=5)

op <- par(mfrow=c(1,3))

#spp <- "WTSP"
for (spp in names(ALL)) {

VV <- ALL[[spp]]

TS <- do.call(cbind, VV[c("TS1", "TS2", "TS3", "TS4")])
TR <- do.call(cbind, VV[c("TR1", "TR2", "TR3", "TR4")])
colnames(TS) <- colnames(TR) <- paste0("Visit", 1:4)
EX <- sapply(VV[c("EX2x2", "EX3x3", "EX4x4", "EX5x5", "EX10x10")], function(z) z[,1])
GR <- sapply(VV[c("GR2x2", "GR3x3", "GR4x4", "GR5x5", "GR10x10")], function(z) z[,1])
EXr <- sapply(VV[c("EX2x2", "EX3x3", "EX4x4", "EX5x5", "EX10x10")], function(z)
    apply(z, 1, quantile, c(0.025, 0.975), na.rm=TRUE))
GRr <- sapply(VV[c("GR2x2", "GR3x3", "GR4x4", "GR5x5", "GR10x10")], function(z)
    apply(z, 1, quantile, c(0.025, 0.975), na.rm=TRUE))
rownames(EXr) <- rownames(GRr) <- c("COR_lower", "COR_upper", "CORU_lower", "CORU_upper")


plot(1:5, EX[1,], ylim=c(-1,1), type="l", col=2, xlab="Extent", ylab="Correlation",
    main=paste(spp, "Extent (intensity fixed)"), axes=FALSE)
axis(2)
box()
axis(1, 1:5, c("2x2", "3x3", "4x4", "5x5", "10x10"))
lines(1:5, EX[2,], ylim=c(0,1), type="l", col=4)
polygon(c(1:5, 5:1), c(EXr[1,], rev(EXr[2,])), border=NA, col="#FF000044")
polygon(c(1:5, 5:1), c(EXr[3,], rev(EXr[4,])), border=NA, col="#0000FF44")
legend("bottomright", col=c(2,4), lty=1, legend=c(
    "Speraman rho", "Uncentered cor."), cex=0.6, bty="n")
abline(h=0)

plot(1:5, GR[1,], ylim=c(-1,1), type="l", col=2, xlab="Intensity", ylab="Correlation",
    main="Intensity (extent fixed)", axes=FALSE)
axis(2)
box()
axis(1, 1:5, c("2x2", "3x3", "4x4", "5x5", "10x10"))
lines(1:5, GR[2,], ylim=c(0,1), type="l", col=4)
polygon(c(1:5, 5:1), c(GRr[1,], rev(GRr[2,])), border=NA, col="#FF000044")
polygon(c(1:5, 5:1), c(GRr[3,], rev(GRr[4,])), border=NA, col="#0000FF44")
legend("bottomright", col=c(2,4), lty=1, legend=c(
    "Speraman rho", "Uncentered cor."), cex=0.6, bty="n")
abline(h=0)

plot(1:4, TS[1,], ylim=c(-1,1), type="l", col=2, xlab="Number of visits", ylab="Correlation",
    main="Temporal", axes=FALSE)
axis(2)
box()
axis(1, 1:4, 1:4)
lines(1:4, TS[2,], col=4)
lines(1:4, TR[1,], col=2, lty=2)
lines(1:4, TR[2,], col=4, lty=2)
legend("bottomright", col=c(2,2,4,4), lty=c(1,2,1,2), legend=c(
    "Speraman rho, sequential", "Speraman rho, random", "Uncentered cor., sequential", "Uncentered cor., random"),
    cex=0.6, bty="n")
abline(h=0)

}

par(op)

dev.off()

#' ------------- old stuff -----------------
#'
#'
#' This code goes over all the model stages

spp <- "WTSP"
Groups <- DATv$ABMIsite
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
V <- try(c(list(validate(res, Groups=Groups, stage=0)),
    lapply(names(mods)[1:9], function(z) validate(res, Groups=Groups, stage=z))))
names(V) <- c("Null", names(mods)[1:9])
plotOne(V)


spp <- "OVEN"
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
VV <- list()
for (i in c("BGex 2x2", "BGex 3x3", "BGex 4x4", "BGex 5x5", "BGex 10x10")) {
    Groups <- switch(i,
        "BGex 2x2"=DATv$SS2,
        "BGex 3x3"=DATv$SS3,
        "BGex 4x4"=DATv$SS4,
        "BGex 5x5"=DATv$SS5,
        "BGex 10x10"=DATv$SS10)
    VV[[i]] <- validate(res, Groups=Groups, stage="HF")
}
r <- lapply(VV, randomize, B=199)
a <- lapply(r, function(z) summarize(z))

Sc <- (c(2,3,4,5,10)^2)
op <- par(mfrow=c(1,2))
j <- 1
plot(Sc, sapply(a, function(z) z[j,1]), type="l", ylim=c(0,1), col=4,
    xlab="Scale", ylab="Spearman r", main=spp)
lines(Sc, sapply(a, function(z) z[j,3]), lty=1, col="grey")
lines(Sc, sapply(a, function(z) z[j,2]), lty=2, col="grey")
lines(Sc, sapply(a, function(z) z[j,4]), lty=2, col="grey")
points(Sc, sapply(a, function(z) z[j,1]), col=4,
    pch=ifelse(sapply(a, function(z) z[j,"p-value"] < 0.05), 19, 21))
j <- 2
plot(Sc, sapply(a, function(z) z[j,1]), type="l", ylim=c(0,1), col=4,
    xlab="Scale", ylab="Uncentered r", main="")
lines(Sc, sapply(a, function(z) z[j,3]), lty=1, col="grey")
lines(Sc, sapply(a, function(z) z[j,2]), lty=2, col="grey")
lines(Sc, sapply(a, function(z) z[j,4]), lty=2, col="grey")
points(Sc, sapply(a, function(z) z[j,1]), col=4,
    pch=ifelse(sapply(a, function(z) z[j,"p-value"] < 0.05), 19, 21))
par(op)

#' Now we do it for all species

## Extent & Grain size (or is it sampling intensity?)
library(parallel)
SPP <- colnames(YYv[,colSums(YYv>0)>100])
GR <- c("ABMI 3x3",
    "BGex 2x2", "BGex 3x3", "BGex 4x4", "BGex 5x5", "BGex 10x10")#, # extent
    #"BGgr 2x2", "BGgr 3x3", "BGgr 4x4", "BGgr 5x5", "BGgr 10x10") # grain
All <- list()
cl <- makeCluster(9)
clusterEvalQ(cl, library(mefa4))
clusterEvalQ(cl, source("~/repos/abmianalytics/birds/00-functions.R"))
clusterExport(cl, c("validate", "Xv", "OFFv", "YYv", "SSHv"))
for (spp in SPP) {
    cat("\n", spp);flush.console()
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    clusterEvalQ(cl, rm(res))
    clusterExport(cl, "res")
    All[[spp]] <- list()
    for (i in GR) {
        cat(".")
        Groups <- switch(i,
            "ABMI 3x3"=DATv$ABMIsite,
            "BGex 2x2"=DATv$SS2,
            "BGex 3x3"=DATv$SS3,
            "BGex 4x4"=DATv$SS4,
            "BGex 5x5"=DATv$SS5,
            "BGex 10x10"=DATv$SS10,
            "BGgr 2x2"=DATv$xx4,
            "BGgr 3x3"=DATv$xx3,
            "BGgr 4x4"=DATv$xx2,
            "BGgr 5x5"=DATv$xx1,
            "BGgr 10x10"=DATv$xx0)
        clusterEvalQ(cl, rm(Groups))
        clusterExport(cl, "Groups")
        V <- try(c(list(validate(res, Groups=Groups, stage=0)),
            parLapply(cl, names(mods)[1:9], function(z) validate(res, Groups=Groups, stage=z))))
        if (!inherits(V, "try-error")) {
            names(V) <- c("Null", names(mods)[1:9])
            for (j in 1:length(V)) {
                V[[j]]$mu <- NULL
                V[[j]]$offset <- NULL
            }
            All[[spp]][[i]] <- V
        } else {
            All[[spp]][[i]] <- NULL
        }
    }
}
stopCluster(cl)
save(All, file=file.path(ROOT, "validation-ABMIandBG-results-2019-02-28.RData"))

## Grain size (or is it sampling intensity?)
SPP <- colnames(YYv[DATv$xx0 != "",colSums(YYv>0)>100])
All <- list()
for (spp in SPP) {
    cat("\n", spp);flush.console()
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    All[[spp]] <- list()
    for (i in 1:5) {
        cat(".")
        Groups <- switch(i,
            "1"=DATv$xx0,
            "2"=DATv$xx1,
            "3"=DATv$xx2,
            "4"=DATv$xx3,
            "5"=DATv$xx4)
        V <- try(c(list(validate(res, Groups=Groups, stage=0)),
            lapply(names(mods)[1:9], function(z) validate(res, Groups=Groups, stage=z))))
        if (!inherits(V, "try-error")) {
            names(V) <- c("Null", names(mods)[1:9])
            All[[spp]][[i]] <- V
        } else {
            All[[spp]][[i]] <- NULL
        }
    }
}
save(All, file=file.path(ROOT, "validation-BGgrain-results-2019-02-19.RData"))

COR <- matrix(NA, length(All), 5)
rownames(COR) <- names(All)
colnames(COR) <- c(
    "1"="BG 0",
    "2"="BG 1",
    "3"="BG 2",
    "4"="BG 3",
    "5"="BG 4")
for (spp in names(All)) {
    COR[spp,] <- sapply(All[[spp]], function(z) z$HF$CORU)
}
plotrix::ladderplot(COR[rowSums(is.na(COR))==0,], ylab="Uncentered correlation")
#lines(1:5, rep(0.5,5),lwd=2, col=2)
COR["CAWA","BG 3"] <- 0.475
lines(1:5, COR["CAWA",], col=2, lwd=2)
lines(1:5, COR["AMRO",], col=4, lwd=2)


Groups <- DATv$SS3
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
V <- try(c(list(validate(res, Groups=Groups, stage=0)),
    lapply(names(mods)[1:9], function(z) validate(res, Groups=Groups, stage=z))))
names(V) <- c("Null", names(mods)[1:9])
plotOne(V)




pdf(file.path(ROOT, "validation-ABMIandBG-results-2019-02-01.pdf"), onefile=TRUE, height=8, width=12)
for (spp in names(All)) {
    for (i in 1:6) {
        txt <- switch(i,
            "1"="ABMI 3x3",
            "2"="BG 2x2",
            "3"="BG 3x3",
            "4"="BG 4x4",
            "5"="BG 5x5",
            "6"="BG full")
        try(plotOne(All[[spp]][[i]], q=1, txt=txt))
    }
}
dev.off()

COR <- matrix(NA, length(All), 6)
rownames(COR) <- names(All)
colnames(COR) <- c("1"="ABMI 3x3",
    "2"="BG 2x2",
    "3"="BG 3x3",
    "4"="BG 4x4",
    "5"="BG 5x5",
    "6"="BG full")
for (spp in names(All)) {
    COR[spp,] <- sapply(All[[spp]], function(z) z$HF$CORU)
}
plotrix::ladderplot(COR[rowSums(is.na(COR))==0,-1], ylab="Uncentered correlation")

cr <- COR[rowSums(is.na(COR))==0,-1]
A <- c(2, 3, 4, 5, 10)
mn <- apply(cr, 2, mean)
cf <- coef(lm(mn[4:5] ~ A[4:5]))
E <- (1-cf[1])/cf[2]

plot(A^2, cr[1,], ylim=c(0,1), type="n", xlim=c(0,E^2), xlab="Scale (km^2)", ylab="Uncentered corr")
for (i in 1:nrow(cr))
    lines(A^2, cr[i,], col="grey")
lines(A^2, mn, lwd=2)
points(A^2, mn, pch=19)
abline(h=1, lty=2)
AA <- seq(A[5], E, length.out = 10)
rr <- cf[1]+cf[2]*AA
lines(AA^2, rr, lty=2, lwd=2)


#'
#' Below this line is just parking lot for code pieces -- DON'T USE!
#'


vv <- lapply(All, function(One) t(sapply(One, function(z) sapply(z[c("y1obs","lam1","yobs","lam")], mean))))
vvv <- t(sapply(vv, function(z) z["HF",]))

op <- par(mfrow=c(2,2))
plot(lam1 ~ y1obs, vvv, main="point", ylab="expected", xlab="observed")
abline(0,1,lty=2)
abline(lm(lam1 ~ y1obs, data.frame(vvv)), col=2)

plot(lam ~ yobs, vvv, main="landscape", ylab="expected", xlab="observed")
abline(0,1,lty=2)
abline(lm(lam ~ yobs, data.frame(vvv)), col=2)

plot(density(vvv[,"lam1"] / vvv[,"y1obs"]), main="", xlim=c(0,3), xlab="expected/observed")
abline(v=1, lty=2)
plot(density(vvv[,"lam"] / vvv[,"yobs"]), main="", xlim=c(0,3), xlab="expected/observed")
abline(v=1, lty=2)
par(op)

mui <- X %*% t(est)
pri <- apply(exp(mui), 1, median)
offi <- OFF[,spp]
evi <- apply(exp(mui+offi), 1, median)
yi <- as.numeric(YY[,spp])
table(yi)

muo <- Xv %*% t(est)
pro <- apply(exp(muo), 1, median)
offo <- OFFv[,spp]
evo <- apply(exp(muo+offo), 1, median)
yo <- as.numeric(YYv[,spp])
table(yo)

(AUCi <- simple_auc(simple_roc(ifelse(yi>0, 1, 0), evi)))
(AUCo <- simple_auc(simple_roc(ifelse(yo>0, 1, 0), evo)))

MAX <- max(quantile(evi, 0.99), quantile(evo, 0.99))
boxplot(evi ~ yi, main=spp, sub="internal", xlab="Observed count", ylab="Expected value", ylim=c(0,MAX))
boxplot(evo ~ yo, main=spp, sub="internal", xlab="Observed count", ylab="Expected value", ylim=c(0,MAX))


boxplot(pri ~ DAT$vegc, cex.axis=0.6)

Predm <- groupSums(exp(muo+offo), 1, DATv$ABMIsite)
Y <- sum_by(yo, DATv$ABMIsite)
all(rownames(Predm)==rownames(Y))

q <- 1
y9obs <- Y[rownames(Y) != "","x"]
y9pred <- apply(Predm[rownames(Predm) != "",], 1, median)
MAX <- max(quantile(y9obs, q), quantile(y9pred, q))

plot(y9obs, y9pred, xlim=c(0,MAX), ylim=c(0,MAX))
abline(0,1)
abline(lm(y9pred ~ y9obs), col=2)
abline(lm(y9pred ~ y9obs-1), col=2)

cor(y9obs, y9pred, method="spearman")




validate <- function(res, Xv, stage=NULL, offset=TRUE) {
    muo <- predict_with_SSH(res, Xv, SSHv, stage=stage)
    pro <- apply(exp(muo), 1, median)
    offo <- if (offset)
        OFFv[,spp] else 0
    evo <- apply(exp(muo+offo), 1, median)
    yo <- as.numeric(YYv[,spp])
    Predm <- groupSums(exp(muo+offo), 1, DATv$ABMIsite)
    Y <- sum_by(yo, DATv$ABMIsite)
    stopifnot(all(rownames(Predm)==rownames(Y)))
    y9obs <- Y[rownames(Y) != "","x"]
    y9pred <- apply(Predm[rownames(Predm) != "",], 1, median)
    AUCo <- simple_auc(simple_roc(ifelse(yo>0, 1, 0), evo))
    CORo <- cor(y9obs, y9pred, method="spearman")
    list(AUC=AUCo, COR9=CORo,
        y1obs=yo, lam1=evo, y9obs=y9obs, lam9=y9pred)
}
plotOne <- function(spp, q=1, All) {
    One <- All[[spp]]
    lam1 <- One$HF$lam1
    y1 <- One$HF$y1obs
    lam9 <- One$HF$lam9
    y9 <- One$HF$y9obs

    op <- par(mfrow=c(2,2))
    plot(0:9, sapply(One, "[[", "COR9"), type="l", col="#0000FF80", ylim=c(-1,1),
        xlab="Model stages", ylab="Correlation", main=spp, lwd=2)
    abline(h=0, lty=2)
    text(0:9, rep(-0.1, 10), round(sapply(One, "[[", "COR9"), 3), cex=0.6)
    plot(0:9, sapply(One, "[[", "AUC"), type="l", col="#0000FF80", ylim=c(0,1),
        xlab="Model stages", ylab="AUC", lwd=2)
    abline(h=0.5, lty=2)
    text(0:9, rep(0.45, 10), round(sapply(One, "[[", "AUC"), 3), cex=0.6)

    boxplot(lam1 ~ y1, xlab="Observed count at point", ylab="Expected value",
        ylim=c(0,quantile(lam1, q)), col="#0000FF80")

    MAX <- max(quantile(lam9, q), quantile(y9, q))
    plot(jitter(y9), lam9, xlim=c(0,MAX), ylim=c(0,MAX), xlab="Observed count at ABMI site",
        ylab="Expected value", col="#0000FF80", pch=19)
    abline(0,1, lty=2)
    abline(lm(lam9 ~ y9), col=2)
    abline(lm(lam9 ~ y9-1), col=2, lty=2)
    par(op)
    invisible(NULL)
}


spp <- "WTSP"
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
V1 <- V <- try(c(list(validate(res, Xv, 0)),
    lapply(names(mods)[1:9], function(z) validate(res, Xv, z, TRUE))))
V2 <- V <- try(c(list(validate(res, Xv, 0)),
    lapply(names(mods)[1:9], function(z) validate(res, Xv, z, FALSE))))
res2 <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
V3 <- V <- try(c(list(validate(res2, Xv, 0)),
    lapply(names(mods)[1:9], function(z) validate(res2, Xv, z, TRUE))))
names(V1) <- names(V2) <- names(V3) <- c("Null", names(mods)[1:9])
plotOne("WTSP", All=list(WTSP=V1))
plotOne("WTSP", All=list(WTSP=V2))
plotOne("WTSP", All=list(WTSP=V3))

SPP <- colnames(YYv[DATv$ABMIsite != "",colSums(YYv>0)>20])

All <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    V <- try(c(list(validate(res, Xv, 0)),
        lapply(names(mods)[1:9], function(z) validate(res, Xv, z))))
    if (!inherits(V, "try-error")) {
        names(V) <- c("Null", names(mods)[1:9])
        All[[spp]] <- V
    }
}




save(All, file="d:/abmi/AB_data_v2018/data/analysis/birds/validation-quick-results-2018-12-14.RData")

pdf("d:/abmi/AB_data_v2018/data/analysis/birds/validation-quick-results-2018-12-14.pdf",
    onefile=TRUE, height=10, width=10)
for (spp in colnames(YYv[DATv$ABMIsite != "",colSums(YYv>0)>20]))
    plotOne(spp, q=1, All=All)
dev.off()



load("d:/abmi/AB_data_v2018/data/analysis/birds/validation-quick-results-2018-12-12.RData")

coef9 <- function(spp, int=TRUE) {
    One <- All[[spp]]
    lam9 <- One$HF$lam9
    y9 <- One$HF$y9obs
    if (int)
        coef(lm(lam9 ~ y9)) else coef(lm(lam9 ~ y9-1))
}

cf9 <- t(sapply(names(All), coef9))
cf9 <- cf9[cf9[,2] < 2,]
cf9x <- sapply(rownames(cf9), coef9, int=FALSE)

hist(cf9[,2])
hist(cf9x)

library(lhreg)
data("lhreg_data")

Mass <- lhreg_data$mass[match(rownames(cf9), lhreg_data$spp)]
EDR <- exp(lhreg_data$logtau[match(rownames(cf9), lhreg_data$spp)])

plot(cf9x ~ log(Mass))
abline(lm(cf9x ~ log(Mass)))

plot(cf9[,2] ~ log(Mass))
abline(lm(cf9[,2] ~ log(Mass)))


cor(cf9x, Mass)
cor(cf9[,2], Mass)
cor(cf9x, EDR)
cor(cf9[,2], EDR)

## compare SM/RF estimates between validation and north
SPP <- colnames(YY)
stage <- "Water"
#spp <- "ALFL"
res <- matrix(0, length(SPP), 4)
dimnames(res) <- list(SPP, c("N_SM", "N_RF", "V_SM", "V_RF"))
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    res1 <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    est1 <- get_coef(res1, X, stage=stage, na.out=FALSE)
    res2 <- load_species(file.path(ROOT, "out", "validation", paste0(spp, ".RData")))
    est2 <- get_coef(res2, X, stage=stage, na.out=FALSE)
    res[spp, ] <- c(est1[1,c("CMETHODSM","CMETHODRF")], est2[1,c("CMETHODSM","CMETHODRF")])
}
res <- res[rownames(cf9),]

plot(res[,c(1,3)], xlim=c(-5,5), ylim=c(-5,5));abline(0,1)
plot(res[,c(2,4)], xlim=c(-5,5), ylim=c(-5,5));abline(0,1)

plot(cf9[,2] ~ I(res[,1]-res[,3]))
abline(lm(cf9[,2] ~ I(res[,1]-res[,3])))

plot(cf9[,2] ~ I(res[,4]-res[,2]), xlab="RF est diff: Validation-North",ylab="Validation slope (ABMI 9 pt)")
abline(lm(cf9[,2] ~ I(res[,4]-res[,2])))

plot(cf9[,2] ~ res[,2])
abline(lm(cf9[,2] ~ res[,2]))



getAUC <- function(spp) {
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    mu1 <- predict_with_SSH(res, X[BB[,1],], SSH[BB[,1],], stage="Space")
    mu2 <- predict_with_SSH(res, X[BB[,1],], SSH[BB[,1],], stage="HF")
    pr1 <- apply(exp(mu1+OFF[BB[,1],spp]), 1, median)
    pr2 <- apply(exp(mu2+OFF[BB[,1],spp]), 1, median)
    y <- as.numeric(YY[BB[,1],spp])

    c(Space=simple_auc(simple_roc(ifelse(y>0, 1, 0), pr1)),
        Landsc=simple_auc(simple_roc(ifelse(y>0, 1, 0), pr2)))
}

z <- read.csv("~/repos/abmispecies/_data/birds.csv")
SPP <- intersect(colnames(OFF), as.character(z$AOU[z$modelN]))
library(parallel)
cl <- makeCluster(4)
clusterEvalQ(cl, library(Matrix))
clusterExport(cl, c("X", "SSH", "YY", "simple_roc", "simple_auc", "predict_with_SSH", "load_species", "ROOT", "BB"))
AUC <- pbapply::pbsapply(SPP, getAUC)
AUC <- list()
for (spp in SPP) {
    cat(spp, "\n")
    AUC[[spp]] <- getAUC(spp)

}
stopCluster(cl)

getMID <- function(spp) {
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    mid <- get_mid(res)
    c(SHH=sum(mid[,"SSH"] > 0), HF=sum(mid[,"HF"] > 0))
}
MID <- pbapply::pbsapply(SPP, getMID)
summary(t(MID)/max(MID))

## process big grids

table(droplevels(DATv$PCODE))

bg <- droplevels(nonDuplicated(DATv[DATv$PCODE == "BU_BG",], SS, TRUE))
ss <- as.character(bg$SS)
tmp <- strsplit(ss, "::")
z <- data.frame(bg[,c("X", "Y")],
    ProjectID=sapply(tmp, "[[", 1),
    Cluster=sapply(tmp, "[[", 2),
    SITE=as.integer(sapply(tmp, "[[", 3)),
    STATION=as.integer(sapply(tmp, "[[", 4)))
z <- z[order(z$SITE, z$STATION),]

## problematic grids:
## 1: missing--exclude
## (7: last bit is weird, use <=100 stations)
## 14: something is not quite right up here--exclude
f <- function(i) {
    with(z[z$SITE==i & z$STATION<=100,], plot(X, Y, type="l",col=2))
    #with(z[z$SITE==i,], points(X, Y,col=4))
    with(z[z$SITE==i & z$STATION<=100,], text(X, Y, STATION,col=4))
    invisible()
}

z <- z[z$SITE %in% c(2:13) & z$STATION<=100,]

x <- data.frame(STATION=1:100, xv=rep(1:10,10), yv=rep(1:10,each=10))

c2 <- c(0,2,4,6,8,10)+0.5
x$g2 <- interaction(cut(x$xv, c2, labels=FALSE), cut(x$yv, c2, labels=FALSE), sep="_", drop=TRUE)

c3 <- c(0,3,6,9,12)+0.5
x$g3 <- interaction(cut(x$xv, c3, labels=FALSE), cut(x$yv, c3, labels=FALSE), sep="_", drop=TRUE)

c4 <- c(0,4,8,12)+0.5
x$g4 <- interaction(cut(x$xv, c4, labels=FALSE), cut(x$yv, c4, labels=FALSE), sep="_", drop=TRUE)

c5 <- c(0,5,10)+0.5
x$g5 <- interaction(cut(x$xv, c5, labels=FALSE), cut(x$yv, c5, labels=FALSE), sep="_", drop=TRUE)

x$g10 <- as.factor("1_1")

z <- data.frame(z, x[match(z$STATION, x$STATION),])

write.csv(z, row.names=TRUE, file=file.path(ROOT, "validation-BGgroups.csv"))

## looking at ABMI sites and landscape effects

library(mgcv)
SPP <- colnames(YYv[,colSums(YYv>0)>100])

#spp <- "OVEN"
pdf(file.path(ROOT, "landscape-pred.pdf"), width=12, height=12, onefile=TRUE)

for (spp in SPP) {
    cat(spp, "\n")
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    res <- res[!sapply(res, inherits, "try-error")]

    mu0 <- predict_with_SSH(res, Xv, SSHv, stage="ARU")
    mu1 <- predict_with_SSH(res, Xv, SSHv, stage="Space")
    mu2 <- predict_with_SSH(res, Xv, SSHv, stage="HF")
    lam0 <- rowMeans(exp(mu0))
    lam1 <- rowMeans(exp(mu1))
    lam2 <- rowMeans(exp(mu2))
    #y <- as.numeric(YYv[,spp])

    q <- quantile(lam0,0.99)
    lam0[lam0>q] <- q
    q <- quantile(lam1,0.99)
    lam1[lam1>q] <- q
    q <- quantile(lam2,0.99)
    lam2[lam2>q] <- q

    lc <- table(unname(unlist(lapply(res, function(z) z$ssh$labels))))/length(res)
    lc <- lc[match(colnames(SSHv), names(lc))]
    names(lc) <- colnames(SSHv)
    lc[is.na(lc)] <- 0
    LC <- names(lc)[lc >= 0.5]

    pSH <- rowSums(SSHv[,LC])
    pHF <- DATv$THF_KM

    m0 <- gam(lam0 ~ s(pSH, pHF))
    m1 <- gam(lam1 ~ s(pSH, pHF))
    m2 <- gam(lam2 ~ s(pSH, pHF))

    LEV <- pretty(c(0, fitted(m2)[fitted(m2)>0]), 6)
    LEV <- LEV[-length(LEV)]

    op <- par(mfrow=c(2,2), las=2, cex.axis=0.6)
    plot(m0, se=FALSE, scheme=2, main=paste(spp, "ARU"), rug=FALSE, levels=LEV)
    plot(m1, se=FALSE, scheme=2, main="Space", rug=FALSE, levels=LEV)
    barplot(rev(lc), horiz=TRUE, col=grey(1-rev(lc)), xlab="selection freq")
    plot(m2, se=FALSE, scheme=2, main="HF", rug=FALSE, levels=LEV)
    par(op)
}
dev.off()

