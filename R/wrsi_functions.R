## wrsi

wrsi <-
    function(Y, X)
    {
        Y <- as.integer(ifelse(Y > 0, 1L, 0L))
        X <- data.matrix(X)
        n <- length(Y)
        if (nrow(X) != n)
            stop("dimension mismatch: X and Y")
        K <- ncol(X)
        if (is.null(colnames(X)))
            colnames(X) <- paste0("V", seq_len(K))
        ## relative suitability index
        ## # of available units of type k / total # of available units (any type)
        Pavail <- colSums(X) / sum(X)
        ## # of used units of type k / total # of used units (any type)
        Xu <- X * Y
        ## sum(Xu) = sum(Y) except when rowsum != 1
        Pused <- colSums(Xu) / sum(Xu)
        ## crude weighted p-occ
        Pw <- colSums(Xu) / colSums(X)
        ## Weighted Relative Suitability Index
        WRSI <- Pused / Pavail
        Var <- (1/colSums(Xu)) - (1/sum(Xu)) + (1/colSums(X)) - (1/sum(X))
        ## fake glm
        #ymat <- cbind(colSums(Xu), colSums(X)-colSums(Xu))
        #xmat <- diag(1,K,K)
        #m <- suppressWarnings(glm(ymat ~ xmat-1, family=binomial))
        res <- data.frame(
            WRSI=WRSI,
            zWRSI=log(WRSI),
            rWRSI=(exp(2 * log(WRSI)) - 1)/(1 + exp(2 * log(WRSI))),
            Pused=Pused,
            Pavail=Pavail,
            Pw=Pw,
            u=colSums(Xu),
            a=colSums(X),
            V=Var)#,
            #coef(summary(m))[,1:2])
        rownames(res) <- colnames(X)
        class(res) <- c("wrsi", "data.frame")
        res
    }

confint.wrsi <-
    function(object, parm, level = 0.95, ...)
    {
        pnames <- rownames(object)
        if (missing(parm))
            parm <- pnames
        else if (is.numeric(parm))
            parm <- pnames[parm]
        if (is.null(parm))
            parm <- 1
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        np <- length(parm)
        pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE,
                            digits = 3), "%")
        ci <- array(NA, dim = c(np, 2L), dimnames = list(parm, pct))
        cf <- log(object$WRSI)
        fac <- qnorm(a)
        ci[] <- exp(log(object[parm, "WRSI"]) + sqrt(object[parm, "V"]) %o% fac)
        ci
    }

pwfun <- function(object, level = 0.95) {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    cf <- object[,"Estimate"]
    se <- object[,"Std..Error"]
    ci <- sapply(1:length(cf), function(i) {
        quantile(rnorm(10^4, cf[i], se[i]), a)
    })
    out <- plogis(t(ci))
    out[object$Pw==0,] <- 0
    out
}


summary.wrsi <- function(object, ...)
    data.frame(Estimate=object$WRSI, exp(confint(object, ...)),
               Pw=object$Pw, pwfun(object, ...))


load("~/Dropbox/Public/OUT_mites_2015-05-22.Rdata")
load("~/Dropbox/Public/veg-hf-clim-reg_abmi-onoff.Rdata")
tveg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")


yyy <- as.matrix(xtab(m2))
hhh <- as.matrix(dd1ha$veg_current)
iii <- intersect(rownames(yyy), rownames(hhh))
yyy <- yyy[iii,]
hhh <- hhh[iii,]
vvv <- as.character(tveg$Broad)
vvv[is.na(vvv)] <- as.character(tveg$UseInAnalysis[is.na(vvv)])
vvv[!is.na(tveg$Broad) & tveg$Broad=="Forest"] <-
    as.character(tveg$Type)[!is.na(tveg$Broad) & tveg$Broad=="Forest"]
hhh <- groupSums(hhh, 2, vvv)
hhh[,"UrbInd"] <- hhh[,"UrbInd"] + hhh[,"Mine"]
hhh <- hhh[,!(colnames(hhh) %in% c("Water","HWater","NonVeg","Mine"))]
hhh <- hhh / ifelse(rowSums(hhh) == 0, 1, rowSums(hhh))

#c("Decid", "Mixwood", "Conif", "Pine",
#  "GrassHerb", "NonVeg",
#  "Shrub", "Water", "Swamp", "Wetland", "HWater", "Cult", "UrbInd",
#  "Mine", "SoftLin", "HardLin")

spp <- 1
w <- wrsi(yyy[,spp], hhh)
rw <- w$rWRSI
names(rw) <- rownames(w)

op <- par(mar=c(4,6,2,2)+0.1, las=1)
tmp <- barplot(rw, horiz=TRUE, xlab="Affinity", main=colnames(yyy)[spp],
    xlim=c(-1,1))
box()
abline(v=0)
par(op)

op <- par(mar=c(6,4,2,2)+0.1, las=2)
tmp <- barplot(rw, horiz=FALSE, ylab="Affinity", main=colnames(yyy)[spp],
    space=NULL, width=sqrt(w$Pavail),
    ylim=c(-1,1))
box()
abline(h=0)
par(op)

w <- wrsi(YY[,spp], ZZ3)
tmp <- barplot(pw[spp,], horiz=TRUE, col=COL, xlab="Pw",
               main=paste0("#occ=",sum(YY[,spp])), xlim=c(0,1))
segments(pw1[spp,], tmp[,1], pw2[spp,], tmp[,1], col="grey", lwd=3)
axis(4, lwd=0, tick=FALSE,
     labels=paste0(round(w$Pw,3)," (", round(w$Pavail,3), ")"), at=tmp[,1])
