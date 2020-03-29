library(mefa4)
source("~/repos/abmianalytics/birds/00-functions.R")
load("d:/abmi/AB_data_v2019/data/misc/bg/bg-data-package.RData")

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
X <- get_model_matrix(DAT, en$mods)
Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

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
    R2 <- summary(lm(yhat ~ yobs, out))$r.squared
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

#save(ALL, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-predictions-ALL.RData")


t(sapply(l, "[[", "stats"))

p_lot <- function(x, k, ...) {
    xx <- x[[as.character(k)]]$results
    Max <- max(xx$yobs, xx$yhat)
    plot(yhat ~ jitter(yobs), xx, ylim=c(0, Max), xlim=c(0, Max),
        pch=19, col="grey", xlab="Observed", ylab="Predicted", ...)
    abline(0,1,lty=2)
    abline(lm(yhat ~ yobs, xx))
    invisible()
}

op <- par(mfrow=c(2,3))
p_lot(l, 1, main=spp)
mtext("Point", at=0, adj=0, line=-2)
p_lot(l, 2)
mtext("2x2", at=0, adj=0, line=-2)
p_lot(l, 3)
mtext("3x3", at=0, adj=0, line=-2)
p_lot(l, 4)
mtext("4x4", at=0, adj=0, line=-2)
p_lot(l, 5)
mtext("5x5", at=0, adj=0, line=-2)
p_lot(l, 10)
mtext("BG", at=0, adj=0, line=-2)
par(op)

## make big df from ALL, add spp, stage, scale and explore
## run randomized stuff, make same big big dataframe 100x and average + CI
## write up/summarize


