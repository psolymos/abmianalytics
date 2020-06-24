library(mefa4)
source("~/repos/abmianalytics/birds/00-functions.R")
load("d:/abmi/AB_data_v2019/data/misc/bg/bg-data-package.RData")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
X <- get_model_matrix(DAT, en$mods)
#Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
#colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

DAT$g2x2 <- paste0(DAT$Site, "_", gg$g2x2[match(DAT$Station, gg$station)])
DAT$g3x3 <- paste0(DAT$Site, "_", gg$g3x3[match(DAT$Station, gg$station)])
DAT$g4x4 <- paste0(DAT$Site, "_", gg$g4x4[match(DAT$Station, gg$station)])
DAT$g5x5 <- paste0(DAT$Site, "_", gg$g5x5[match(DAT$Station, gg$station)])

SPP <- colnames(YY[,colSums(YY>0) >= 100])
m <- structure(c(1,2,3,4,5,10), names=c(1,2,3,4,5,10))

pool <- function(d, k) {
    g <- switch(as.character(k),
        "1"="Key",
        "2"="g2x2",
        "3"="g3x3",
        "4"="g4x4",
        "5"="g5x5",
        "10"="Site")
    Y <- sum_by(d$yobs, d[[g]])
    out <- data.frame(yobs = Y[,"x"],
        yhat = sum_by(d$yhat, d[[g]])[,"x"],
        n = Y[,"by"])
    rownames(out) <- rownames(Y)
    if (g == "g2x2") {
        out <- out[out$n <= 4,]
        out$yobs <- 4 * out$yobs / out$n
        out$yhat <- 4 * out$yhat / out$n
    }
    if (g == "g3x3") {
        out <- out[out$n <= 9,]
        out$yobs <- 9 * out$yobs / out$n
        out$yhat <- 9 * out$yhat / out$n
    }
    if (g == "g4x4") {
        out <- out[out$n <= 16,]
        out$yobs <- 16 * out$yobs / out$n
        out$yhat <- 16 * out$yhat / out$n
    }
    if (g == "g5x5") {
        out <- out[out$n <= 25,]
        out$yobs <- 25 * out$yobs / out$n
        out$yhat <- 25 * out$yhat / out$n
    }
    if (g == "Site") {
        out$yobs <- 100 * out$yobs / out$n
        out$yhat <- 100 * out$yhat / out$n
    }
    COR <- cor(out$yobs, out$yhat, method="spearman")
    CORU <- .cor_uncentered(out$yobs, out$yhat)
    R2 <- summary(lm(yobs ~ yhat, out))$r.squared
    oc <- epiR::epi.occc(cbind(out$yobs, out$yhat))
    PRC <- oc$oprec
    ACC <- oc$oaccu
    list(results=out, stats=c(COR=COR, CORU=CORU, R2=R2, PRC=PRC, ACC=ACC))
}


#spp <- "OVEN"
ALL <- list()
for (spp in SPP) {

    res <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    ALL[[spp]] <- list()

    #stage <- "ARU" # local
    #stage <- "Space" # local + climate
    #stage <- "HF" # local + climate + landscape
    for (stage in c("ARU", "Space", "HF")) {
        cat(spp, stage, "\n")
        flush.console()

        mu <- predict_with_SSH(res, X, SSH, stage=stage)
        off <- OFF[,res[[1]]$species]

        d1 <- data.frame(DAT[,c("Key", "Site", "g2x2", "g3x3", "g4x4", "g5x5")],
            yhat = as.numeric(apply(exp(mu+off), 1, median)),
            yobs = as.numeric(YY[,res[[1]]$species]))
        d1 <- d1[sample(nrow(d1)),]
        d1 <- nonDuplicated(d1, Key, TRUE)


        l <- lapply(m, function(i) pool(d1, i))
        ALL[[spp]][[stage]] <- l
    }
}

#save(ALL, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALL-2020-06-23.RData")

## randomize
B <- 200
stage <- "HF"

ALLW <- list() # WS within site random
ALLB <- list() # WS between site random

for (spp in SPP) {
    cat("\n", spp)
    flush.console()

    res <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    ALLW[[spp]] <- list()
    ALLB[[spp]] <- list()

    for (i in 1:B) {
        if (i %% 10 == 0)
            cat(".")

        mu <- predict_with_SSH(res, X, SSH, stage=stage)
        off <- OFF[,res[[1]]$species]

        ## within site
        G <- data.frame(
            Key  = DAT$Key,
            Site = DAT$Site,
            g2x2 = paste0(DAT$Site, "_", sample(gg$g2x2[match(DAT$Station, gg$station)])),
            g3x3 = paste0(DAT$Site, "_", sample(gg$g3x3[match(DAT$Station, gg$station)])),
            g4x4 = paste0(DAT$Site, "_", sample(gg$g4x4[match(DAT$Station, gg$station)])),
            g5x5 = paste0(DAT$Site, "_", sample(gg$g5x5[match(DAT$Station, gg$station)])))
        d1 <- data.frame(G,
            yhat = as.numeric(apply(exp(mu+off), 1, median)),
            yobs = as.numeric(YY[,res[[1]]$species]))
        d1 <- d1[sample(nrow(d1)),]
        d1 <- nonDuplicated(d1, Key, TRUE)
        l <- lapply(m, function(i) pool(d1, i)$stats)
        ALLW[[spp]][[i]] <- l

        ## between site
        G <- data.frame(
            Key  = DAT$Key,
            Site = DAT$Site,
            g2x2 = paste0(sample(DAT$Site), "_", gg$g2x2[match(DAT$Station, gg$station)]),
            g3x3 = paste0(sample(DAT$Site), "_", gg$g3x3[match(DAT$Station, gg$station)]),
            g4x4 = paste0(sample(DAT$Site), "_", gg$g4x4[match(DAT$Station, gg$station)]),
            g5x5 = paste0(sample(DAT$Site), "_", gg$g5x5[match(DAT$Station, gg$station)]))
        d1 <- data.frame(G,
            yhat = as.numeric(apply(exp(mu+off), 1, median)),
            yobs = as.numeric(YY[,res[[1]]$species]))
        d1 <- d1[sample(nrow(d1)),]
        d1 <- nonDuplicated(d1, Key, TRUE)
        l <- lapply(m, function(i) pool(d1, i)$stats)
        ALLB[[spp]][[i]] <- l

    }
}

#save(ALLW, ALLB, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALLWB-2020-06-23.RData")


library(mefa4)
load("d:/abmi/AB_data_v2019/data/misc/bg/bg-data-package.RData")
load("d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALL.RData")
load("d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALLWB.RData")

SPP <- names(ALL)

tab <- NULL
for (spp in names(ALL)) {
    for (stage in c("ARU", "Space", "HF")) {
        for (i in names(ALL[[spp]][[stage]])) {
            tmp <- data.matrix(ALL[[spp]][[stage]][[i]]$stats)
            tab <- rbind(tab, data.frame(Species=spp, Stage=stage, Scale=i, t(tmp)))
        }
    }
}
summary(tab)

p <- colMeans(YY[,SPP]>0)
tab$pocc <- p[match(tab$Species, names(p))]

## Accuracy ~ bias
## Precision ~ variance
## COR ~ relative abundance
## R2 ~ absolute abundance

plot(tab[,4:8])

boxplot(COR ~ Stage, tab)
boxplot(R2 ~ Stage, tab)

boxplot(COR ~ Scale, tab)
boxplot(CORU ~ Scale, tab)
boxplot(R2 ~ Scale, tab)
#boxplot(CORU ~ Scale, tab[tab$Stage=="HF",])
plot(CORU ~ pocc, tab)
plot(COR ~ pocc, tab)
abline(h=0)

plot(R2 ~ COR, tab)
abline(0,1,col=2)

boxplot(R2 ~ Scale, tab)
boxplot(PRC ~ Scale, tab)
boxplot(ACC ~ Scale, tab)

## all effects move in the positive directions
m <- lm(R2 ~ Stage + Scale + pocc, tab)
m <- lm(COR ~ Stage + Scale + pocc, tab)

m <- step(m)
summary(m)

## spp ranef is really small and uncertain
#library(lme4)
#m <- lmer(COR ~ Stage + Scale + (1|Species), tab)
#summary(m)

#t(sapply(l, "[[", "stats"))

p_lot <- function(x, k, ...) {
    xx <- x[[as.character(k)]]$results
    Max <- max(xx$yobs, xx$yhat)
    plot(yhat ~ jitter(yobs), xx, ylim=c(0, Max), xlim=c(0, Max),
        pch=19, col="#65d7fc", xlab="Observed", ylab="Predicted", ...)
    abline(0,1,lty=2)
    abline(lm(yhat ~ yobs, xx))
    s <- x[[as.character(k)]]$stats
    g <- switch(as.character(k),
        "1"="1x1",
        "2"="2x2",
        "3"="3x3",
        "4"="4x4",
        "5"="5x5",
        "10"="10x10")
    txt <- paste0(c(g, paste0(names(s), " = ", round(s, 3))), collapse="\n")
    mtext(txt, at=0, adj=0, padj=1, line=-1, cex=0.75)
    invisible()
}

pdf("d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALL.pdf", onefile=TRUE, height=7, width=10)
for (spp in names(sort(p))) {
    l <- ALL[[spp]][["HF"]]

    op <- par(mfrow=c(2,3))
    p_lot(l, 1, main=spp)
    p_lot(l, 2)
    p_lot(l, 3)
    p_lot(l, 4)
    p_lot(l, 5)
    p_lot(l, 10)
    par(op)
}
dev.off()

## make big df from ALL, add spp, stage, scale and explore
## run randomized stuff, make same big big dataframe 100x and average + CI
## write up/summarize

lmfun <- function(spp) {
    l <- ALL[[spp]][["HF"]]
    lapply(m, function(k)
        coef(summary(lm(yhat ~ yobs, l[[as.character(k)]]$results))))
}





r_plot <- function(spp, k, r2=FALSE, ...) {
    s <- t(sapply(ALL[[spp]][["HF"]], "[[", "stats"))

    COR <- s[,"COR"]
    R2 <- s[,"R2"]

    CORW <- t(sapply(ALLW[[spp]], function(z) do.call(rbind, z)[,"COR"]))
    CORB <- t(sapply(ALLB[[spp]], function(z) do.call(rbind, z)[,"COR"]))

    R2W <- t(sapply(ALLW[[spp]], function(z) do.call(rbind, z)[,"R2"]))
    R2B <- t(sapply(ALLB[[spp]], function(z) do.call(rbind, z)[,"R2"]))

    if (r2) {
        XW <- R2W
        XB <- R2B
        X <- R2
    } else {
        XW <- CORW
        XB <- CORB
        X <- COR
    }

    q1 <- quantile(XW[,as.character(k)], c(0.025, 0.975))
    q2 <- quantile(XB[,as.character(k)], c(0.025, 0.975))
    d1 <- density(XW[,as.character(k)])
    d2 <- density(XB[,as.character(k)])

    d1$y <- d1$y/max(d1$y)
    d1$x[d1$x < 0] <- 0
    d1$x[d1$x > 1] <- 1
    d2$y <- d2$y/max(d2$y)
    d2$x[d2$x < 0] <- 0
    d2$x[d2$x > 1] <- 1
    plot(0, type="n", ylim=c(-0.05,1.25), xlim=c(0,1), ann=FALSE, axes=FALSE, xaxs = "i", ...)
    axis(1)
    polygon(c(d1$x[1], d1$x, d1$x[length(d1$x)]),
        0.25+c(0, d1$y, 0), col="#ff000044", border=2)
    polygon(c(d2$x[1], d2$x, d2$x[length(d2$x)]),
        c(0, d2$y, 0), col="#0000ff44", border=4)
    polygon(c(q1, rev(q1)), 0.25+c(0.05, 0.05, -0.05, -0.05), col=2, border=2)
    polygon(c(q2, rev(q2)), c(0.05, 0.05, -0.05, -0.05), col=4, border=4)
    abline(v=X[as.character(k)], lwd=3)
    invisible()
}

## red: within
## blue: between
## chop off <0 >1 regions
## express significance (interval) as a kind of rug
## retain axis for bottom only
## add text
spp <- "AMRO"
pdf("d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALLWB.pdf", onefile=TRUE, height=10, width=7)
for (spp in names(sort(p))) {
op <- par(mfrow=c(6, 2), mar=c(3,2,2,2))
for (k in m) {
    r_plot(spp, k, r2=FALSE)
    g <- switch(as.character(k),
        "1"="1x1",
        "2"="2x2",
        "3"="3x3",
        "4"="4x4",
        "5"="5x5",
        "10"="10x10")
    mtext(g, at=0, adj=0, padj=1, line=1, cex=1.25)
    if (k==1)
        mtext("r", line=-1, cex=1.25)
    r_plot(spp, k, r2=TRUE)
    if (k==1) {
        mtext(expression(R^2), line=-1, cex=1.25)
        mtext(spp, adj=1, cex=1.5)
    }
}
par(op)
}
dev.off()


r_stat <- function(spp) {
    s <- t(sapply(ALL[[spp]][["HF"]], "[[", "stats"))
    COR <- s[,"COR"]
    R2 <- s[,"R2"]
    CORW <- apply(sapply(ALLW[[spp]], function(z) do.call(rbind, z)[,"COR"]), 1, median)
    CORB <- apply(sapply(ALLB[[spp]], function(z) do.call(rbind, z)[,"COR"]), 1, median)
    R2W <- apply(sapply(ALLW[[spp]], function(z) do.call(rbind, z)[,"R2"]), 1, median)
    R2B <- apply(sapply(ALLB[[spp]], function(z) do.call(rbind, z)[,"R2"]), 1, median)
    rbind(COR=COR, CORW=CORW, CORB=CORB, R2=R2, R2W=R2W, R2B=R2B)
}

RSTAT <- lapply(SPP, r_stat)
names(RSTAT) <- SPP

boxplot(t(sapply(RSTAT, function(z) z["COR",])))
boxplot(t(sapply(RSTAT, function(z) z["CORW",])))
boxplot(t(sapply(RSTAT, function(z) z["CORB",])))

boxplot(t(sapply(RSTAT, function(z) z["R2",])))
boxplot(t(sapply(RSTAT, function(z) z["R2W",])))
boxplot(t(sapply(RSTAT, function(z) z["R2B",])))

## compare c4i predictions


pool2 <- function(spp, k, grid=TRUE, correct=FALSE) {
    g <- switch(as.character(k),
        "1"="Key",
        "2"="g2x2",
        "3"="g3x3",
        "4"="g4x4",
        "5"="g5x5",
        "10"="Site")

    Yobs <- YY[,spp]

    if (grid && !correct)
        Ypred <- DM[match(DAT$Key, rownames(DM)),spp] * 36
    if (grid && correct)
        Ypred <- DM[match(DAT$Key, rownames(DM)),spp] * exp(OFF[,spp])
    if (!grid && !correct)
        Ypred <- DMp[match(DAT$Key, rownames(DMp)),spp] * (1.5^2*pi)
    if (!grid && correct)
        Ypred <- DMp[match(DAT$Key, rownames(DMp)),spp] * exp(OFF[,spp])

    Yobs <- Yobs[!duplicated(DAT$Key)]
    Ypred <- Ypred[!duplicated(DAT$Key)]


    Yobs <- sum_by(Yobs, DAT[[g]][!duplicated(DAT$Key)])
    Ypred <- sum_by(Ypred, DAT[[g]][!duplicated(DAT$Key)])
    out <- data.frame(yobs = Yobs[,"x"],
        yhat = Ypred[,"x"],
        n = Yobs[,"by"])
    rownames(out) <- rownames(Yobs)
    if (g == "g2x2") {
        out <- out[out$n <= 4,]
        out$yobs <- 4 * out$yobs / out$n
        out$yhat <- 4 * out$yhat / out$n
    }
    if (g == "g3x3") {
        out <- out[out$n <= 9,]
        out$yobs <- 9 * out$yobs / out$n
        out$yhat <- 9 * out$yhat / out$n
    }
    if (g == "g4x4") {
        out <- out[out$n <= 16,]
        out$yobs <- 16 * out$yobs / out$n
        out$yhat <- 16 * out$yhat / out$n
    }
    if (g == "g5x5") {
        out <- out[out$n <= 25,]
        out$yobs <- 25 * out$yobs / out$n
        out$yhat <- 25 * out$yhat / out$n
    }
    if (g == "Site") {
        out$yobs <- 100 * out$yobs / out$n
        out$yhat <- 100 * out$yhat / out$n
    }
    COR <- cor(out$yobs, out$yhat, method="spearman")
    CORU <- .cor_uncentered(out$yobs, out$yhat)
    R2 <- summary(lm(yhat ~ yobs, out))$r.squared
    oc <- epiR::epi.occc(cbind(out$yobs, out$yhat))
    PRC <- oc$oprec
    ACC <- oc$oaccu
    list(results=out, stats=c(COR=COR, CORU=CORU, R2=R2, PRC=PRC, ACC=ACC))
}


#spp <- "OVEN"
sapply(lapply(m, function(k) pool2(spp, k)), "[[", "stats")

ALLC <- list()
for (spp in SPP) {
    ALLC[[spp]] <- list(
        grid=lapply(m, function(k) pool2(spp, k, grid=TRUE, correct=FALSE)),
        gridc=lapply(m, function(k) pool2(spp, k, grid=TRUE, correct=TRUE)),
        point=lapply(m, function(k) pool2(spp, k, grid=FALSE, correct=FALSE)),
        pointc=lapply(m, function(k) pool2(spp, k, grid=FALSE, correct=TRUE))
    )
}

tabc <- NULL
for (spp in names(ALLC)) {
    for (ty in c("grid", "gridc", "point", "pointc")) {
        for (i in names(ALLC[[spp]][[ty]])) {
            tmp <- data.matrix(ALLC[[spp]][[ty]][[i]]$stats)
            tabc <- rbind(tabc, data.frame(Species=spp, Type=ty, Scale=i, t(tmp)))
        }
    }
    for (i in names(ALL[[spp]][["HF"]])) {
        tmp <- data.matrix(ALL[[spp]][["HF"]][[i]]$stats)
        tabc <- rbind(tabc, data.frame(Species=spp, Type="yhat", Scale=i, t(tmp)))
    }
}
summary(tabc)

m <- lm(R2 ~ Type + Scale, tabc)
m <- step(m)
summary(m)

boxplot(COR ~ Type, tabc)
boxplot(R2 ~ Type, tabc)

abline(0,1)
abline(h=0,v=0)
plot(taba$R2, tabc$R2)
abline(0,1)
abline(h=0,v=0)

boxplot(COR ~ Scale, tabc)
boxplot(R2 ~ Scale, tabc)



p_lotc <- function(x, k, type="grid", ...) {
    xx <- x[[type]][[as.character(k)]]$results
    plot(yhat ~ jitter(yobs), xx,
        pch=19, col="#65d7fc", xlab="Observed", ylab="Predicted", ...)
    abline(lm(yhat ~ yobs, xx))
    s <- x[[type]][[as.character(k)]]$stats
    g <- switch(as.character(k),
        "1"="1x1",
        "2"="2x2",
        "3"="3x3",
        "4"="4x4",
        "5"="5x5",
        "10"="10x10")
    txt <- paste0(c(g, paste0(names(s), " = ", round(s, 3))), collapse="\n")
    mtext(txt, at=0, adj=0, padj=1, line=-1, cex=0.75)
    invisible()
}

Type <- "grid"
pdf("d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALLC.pdf", onefile=TRUE, height=7, width=10)
for (spp in names(ALLC)) {
    l <- ALLC[[spp]]
    for (ty in c("grid", "gridc", "point", "pointc")) {

        op <- par(mfrow=c(2,3))
        p_lotc(l, 1, ty, main=paste(spp, ty))
        p_lotc(l, 2, ty)
        p_lotc(l, 3, ty)
        p_lotc(l, 4, ty)
        p_lotc(l, 5, ty)
        p_lotc(l, 10, ty)
        par(op)
    }
}
dev.off()


spp <- "OVEN"
stage <- "HF"
res <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
mu <- predict_with_SSH(res, X, SSH, stage=stage)
off <- OFF[,spp]

d1 <- data.frame(
    DAT[,c("Key", "Site", "g2x2", "g3x3", "g4x4", "g5x5")],
    yhat = as.numeric(apply(exp(mu+off), 1, median)),
    yobs = as.numeric(YY[,spp]),
    mu=as.numeric(apply(exp(mu), 1, median)),
    off=as.numeric(off),
    expoff=as.numeric(exp(off)),
    Dpred=DM[match(DAT$Key, rownames(DM)),spp])
d1 <- d1[sample(nrow(d1)),]
d1 <- nonDuplicated(d1, Key, TRUE)

summary(d1)

## compare scaling: point/square, correction etc
