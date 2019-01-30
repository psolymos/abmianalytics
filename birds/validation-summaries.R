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
load(file.path(ROOT, "data", "ab-birds-validation-2018-12-07.RData"))
#' Make a model matrix that matches the coefficients that we estimated
#' `X` is for training data, `Xv` is for validation data set
X <- get_model_matrix(DAT, mods)
Xv <- get_model_matrix(DATv, mods)
#'
#' # Species summaries
#'
#' Let us pick a species
spp <- "WTSP"
#spp <- "OVEN"
#' Load species results
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
#' Get the coefficient matrix (bootstrap runs x possible coefficients
#' based on all the terms in `mods` list): used `X` as a template to match the
#' estimates with the column names of `X`, and uses `stage` if specified
est <- get_coef(res, X, stage=NULL, na.out=FALSE)
stopifnot(all(colnames(X)==colnames(est)))
#' Here is how you to get a sense of what results look like:
printCoefmat(get_summary(est))
#'
#' # Prediction
#'
#' Here is how to get the null model (you can use `X` and `SSH` for training data
#' and `Xv` and `SSHv` for validation data set).
#' The output is in a PKEY x bootstrap format, PKEY refers to the row names of the input.
#' The matrix represent log densities without offsets
mu0 <- predict_with_SSH(res, Xv, SSHv, stage=0)
#' Here is how to combine with offsets (use `OFF` or `OFFv` depending on previous steps).
#' This is now a PKEY x bootstrap matrix with expected abundances matching the counts
lam0 <- exp(mu0 + OFFv[,spp])
#' We can visualize the results in a boxplot using the observed counts
#' (again, use `YY` or `YYv`)
yv <- as.numeric(YYv[,spp])
#' We use `rowMeans(lambda)` as the bootstrap smoothed mean,
#' `apply(lambda, 1, median)` gives the median
#' and `apply(lambda, 1, quantile, c(0.5, 0.05, 0.95))` gives the median and 90% CI
#' (`lambda` is any matrix starting with `lam*` below)
boxplot(rowMeans(lam0) ~ yv)
#' Let us add covariates now.
#' Just to clarify, the `stage` argument can be
#' - `NULL`: use all the model stages, this is the default (but it will include Year effect!)
#' - `0`: intercept only null model
#' - any value from `names(mods)`
names(mods)
#' Pick for example `"ARU"` stage
lamARU <- exp(predict_with_SSH(res, Xv, SSHv, stage="ARU") + OFFv[,spp])
boxplot(rowMeans(lamARU) ~ yv)
#' Pick for example `"ARU"` stage for all the local effects
lamARU <- exp(predict_with_SSH(res, Xv, SSHv, stage="ARU") + OFFv[,spp])
boxplot(rowMeans(lamARU) ~ yv)
#' Pick `"HF"` stage for all local and surrounding effects
lamHF <- exp(predict_with_SSH(res, Xv, SSHv, stage="HF") + OFFv[,spp])
boxplot(rowMeans(lamHF) ~ yv)
#' Now let's compare ROC/AUC metrics (`simple_roc` is a bare bones function
#' that is blazing fast compared to pROC)
ROC0 <- simple_roc(ifelse(yv > 0, 1, 0), rowMeans(lam0))
ROCARU <- simple_roc(ifelse(yv > 0, 1, 0), rowMeans(lamARU))
ROCHF <- simple_roc(ifelse(yv > 0, 1, 0), rowMeans(lamHF))
AUC0 <- simple_auc(ROC0)
AUCARU <- simple_auc(ROCARU)
AUCHF <- simple_auc(ROCHF)
#' Plot the ROC curves
plot(ROC0[,2:1], type="l", col=2)
abline(0,1,col=1,lty=2)
lines(ROCARU[,2:1], col=3)
lines(ROCHF[,2:1], col=4)
legend("bottomright", lty=1, col=2:4, bty="n",
    legend=paste0(c("Null", "ARU", "HF"),
    " (AUC=", round(c(AUC0, AUCARU, AUCHF), 3), ")"))
#'
#' # Scaling up to ABMI grids of 9
#'
#' Functions for the magic: uses only external data
#' Groups argument defines the spatial units to aggregate
#' like ABMIsite etc, must match DATv row ordering
#' unwanted rows should be blank (`""`) or NA
validate <- function(res, Groups, stage=NULL) {
    Groups <- as.character(Groups)
    Groups[is.na(Groups)] <- ""
    muo <- predict_with_SSH(res, Xv, SSHv, stage=stage)
    pro <- apply(exp(muo), 1, median)
    offo <- OFFv[,spp]
    evo <- apply(exp(muo+offo), 1, median)
    yo <- as.numeric(YYv[,spp])

    Predm <- groupSums(exp(muo+offo), 1, Groups)
    Y <- sum_by(yo, Groups)
    yobs <- Y[rownames(Y) != "","x"]
    ypred <- apply(Predm[rownames(Predm) != "",], 1, median)
    ypred <- ypred[names(yobs)]
    m1 <- get_mass(yo, evo)
    m <- get_mass(yobs, ypred)

    AUCo <- simple_auc(simple_roc(ifelse(yo>0, 1, 0), evo))
    CORo <- cor(yobs, ypred, method="spearman")
    CORUo <- .cor_uncentered(yobs, ypred)
    list(spp=res[[1]]$species,
        AUC=AUCo, COR=CORo, CORU=CORUo,
        y1obs=yo, lam1=evo, yobs=yobs, lam=ypred,
        cmf1=m1, cmf=m)
}


plotOne <- function(One, q=1) {
    spp <- One$HF$spp
    lam1 <- One$HF$lam1
    y1 <- One$HF$y1obs
    lam <- One$HF$lam
    y <- One$HF$yobs
    p1 <- sapply(One, function(z) z$cmf1[,"expected"])
    p <- sapply(One, function(z) z$cmf[,"expected"])
    Dev1 <- sapply(One, function(z) attr(z$cmf1, "dev"))
    Dev <- sapply(One, function(z) attr(z$cmf, "dev"))

    op <- par(mfrow=c(2,3))

    plot(0:9, sapply(One, "[[", "COR"), type="l", col="#0000FF80", ylim=c(-1,1),
        xlab="Model stages", ylab="Correlation / AUC", main=spp, lwd=2)
    lines(0:9, sapply(One, "[[", "CORU"), col="#FF000080", lwd=2)
    lines(0:9, sapply(One, "[[", "AUC"), col="#00FF0080", lwd=2)
    abline(h=c(0, 0.5), lty=2)
    text(0:9, rep(-0.1, 10), round(sapply(One, "[[", "COR"), 3), cex=0.6, col="#0000FF")
    text(0:9, rep(-0.2, 10), round(sapply(One, "[[", "CORU"), 3), cex=0.6, col="#FF0000")
    text(0:9, rep(-0.3, 10), round(sapply(One, "[[", "AUC"), 3), cex=0.6, col="#00FF00")
    legend("bottomright", bty="n", lwd=2, lty=1, col=c("#0000FF80", "#FF000080", "#00FF0080"),
        legend=c("Spearman r", "Uncentered r", "AUC"))

    plot(0:9, Dev1, type="l", col="#0000FF80", ylim=c(0,max(Dev,Dev1)),
        xlab="Model stages", ylab="Goodness of fit", lwd=2)
    lines(0:9, Dev, col="#FF000080", lwd=2)
    legend("topright", bty="n", lwd=2, lty=1, col=c("#0000FF80", "#FF000080"),
        legend=c("Point", "Landscape"))

    ss <- c("Null", "ARU", "Space", "HF")
    cg <- 1:4 # RColorBrewer::brewer.pal(length(ss)+3, "Dark2")[-(1:3)]
    matplot(c(0, One$HF$cmf1[,"observed"]), rbind(0, p1[,ss]),
        type="l", lty=1, ylim=c(0,1), xlim=c(0,1), col=cg,
        xlab="Observed (point)", ylab="Expected (point)")
    legend("topleft", col=cg, lty=1, legend=ss, bty="n")
    abline(0,1,lty=2)

    boxplot(lam1 ~ y1, xlab="Observed count at point", ylab="Expected value",
        ylim=c(0,quantile(lam1, q)), col="#0000FF80")

    MAX <- max(quantile(lam, q), quantile(y, q))
    plot(jitter(y), lam, xlim=c(0,MAX), ylim=c(0,MAX), xlab="Observed count at ABMI site",
        ylab="Expected value", col="#0000FF80", pch=19)
    abline(0,1, lty=2)
    abline(lm(lam ~ y), col=2)
    abline(lm(lam ~ y-1), col=2, lty=2)


    matplot(c(0, One$HF$cmf[,"observed"]), rbind(0,p[,ss]),
        type="l", lty=1, ylim=c(0,1), xlim=c(0,1), col=cg,
        xlab="Observed (landscape)", ylab="Expected (landscape)")
    legend("topleft", col=cg, lty=1, legend=ss, bty="n")
    abline(0,1,lty=2)

    par(op)
    invisible(NULL)
}
#' This code goes over all the model stages
spp <- "WTSP"
Groups <- DATv$ABMIsite
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
V <- try(c(list(validate(res, Groups=DATv$ABMIsite, stage=0)),
    lapply(names(mods)[1:9], function(z) validate(res, Groups=Groups, stage=z))))
names(V) <- c("Null", names(mods)[1:9])
plotOne(V)
#' Now we do it for all species
SPP <- colnames(YYv[DATv$ABMIsite != "",colSums(YYv>0)>20])
All <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    V <- try(c(list(validate(res, Groups=DATv$ABMIsite, stage=0)),
        lapply(names(mods)[1:9], function(z) validate(res, Groups=Groups, stage=z))))
    if (!inherits(V, "try-error")) {
        names(V) <- c("Null", names(mods)[1:9])
        All[[spp]] <- V
    }
}
save(All, file=file.path(ROOT, "validation-quick-results-2019-01-17.RData"))

pdf(file.path(ROOT, "validation-quick-results-2019-01-17.pdf"), onefile=TRUE, height=8, width=12)
for (spp in names(All))
    plotOne(All[[spp]], q=1)
dev.off()
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