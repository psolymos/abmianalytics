require(RColorBrewer)

getOK <- function(res) {
    sapply(res, class) != "try-error"
}

loadSPP <- function(spp, dir=NULL, prefix=NULL) {
    if (is.null(dir))
        dir <- RESDIR
    if (is.null(prefix))
        prefix <- "birds_josm_"
    sppfn <- paste0(dir, prefix, spp, ".Rdata") 
    if (!file.exists(sppfn)) {
        cat("species results file does not exist\n")
        return(NULL)
    }
    e <- new.env()
    load(sppfn, envir=e)
    res <- e$res
#    res <- res[getOK(res)]
    res
}


fix_underscore <- function(x) {
        x <- lapply(x, strsplit, "_")
        x <- sapply(x, function(z) paste(z[[1]], collapse="\\_", sep=""))
        x
}

getTerms <- function(mods, type=c("formula", "list"), intercept=TRUE) {
    type <- match.arg(type)
    x <- unlist(lapply(unlist(mods), function(z) as.character(z)[3]))
#    x <- unname(substr(x, 5, nchar(x)))
    x <- gsub(". + ", "", x, fixed=TRUE)
    x <- unlist(strsplit(x, "+", fixed=TRUE))
    x <- unlist(strsplit(x, "*", fixed=TRUE))
    if (type == "list")
        x <- unlist(strsplit(x, ":", fixed=TRUE))
    x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
    x <- unique(x)
    if (type == "formula") {
        x <- paste("~", paste(x, collapse=" + ", sep=""))
        if (!intercept)
            x <- paste(x, "- 1")
        x <- as.formula(x)
    }
    x
}

fixNames <- function(x, sep=":") {
    unlist(lapply(x, function(z) {
        paste(sort(strsplit(z, sep)[[1]]), collapse=sep)
    }))
}

Grepl <- function(pattern, x) {
    rowSums(sapply(pattern, function(z) grepl(z, x))) > 0
}

getEst <- function(res, stage=NULL, na.out=TRUE) {
    OK <- !sapply(res, inherits, "try-error")
    if (any(!OK))
        warning(paste("try-error found:", sum(!OK)))
    ii <- sapply(res[OK], "[[", "iteration")
    est <- Xn[1:length(ii),,drop=FALSE]
    rownames(est) <- ii
    est[] <- 0
    if (is.null(stage))
        stage <- length(res[[ii[1]]]$coef)
    if (stage > 0) {
        for (i in 1:length(ii)) {
            tmp <- res[[ii[i]]]$coef[[stage]]
            names(tmp) <- fixNames(names(tmp))
            est[i,match(names(tmp), colnames(est))] <- tmp
        }
    } else {
        for (i in 1:length(ii)) {
            est[i,1] <- res[[ii[i]]]$null
        }
    }
    if (any(!OK) && na.out) {
        nas <- matrix(NA, sum(!OK), ncol(est))
        rownames(nas) <- which(!OK)
        est <- rbind(est, nas)
    }
    est
}

## add in option for offsets ???
getDataPred <- function(res, stage, X=NULL, remn=FALSE) {
    OK <- which(getOK(res))
    est <- getEst(res, stage, na.out=FALSE)
    if (is.null(X))
        X <- Xn
    tX <- t(X[,colnames(est)])
    if (remn) {
        mu <- pbsapply(1:nrow(est), function(i) {
            tX["Remn_QS",] <- rowSums(XHF[, res[[OK[i]]]$hi, drop=FALSE])
            tX["Remn2_QS",] <- tX["Remn_QS",]^2
            crossprod(tX, est[i,])
        })
    } else {
        mu <- pbapply(est, 1, function(z) crossprod(tX, z))
    }
    mu
}
#mu <- getDataPred(res)
#summary(apply(exp(mu), 1, median, na.rm=TRUE))

getCaic <- function(res, stage=NULL, run=1) {
    OK <- !sapply(res, inherits, "try-error")
    ii <- sapply(res[OK], "[[", "iteration")
    if (is.null(stage))
        stage <- length(res[[ii[1]]]$coef)
    if (stage == 0)
        return(attr(res[[run]]$caic[[1]], "StartCAIC"))
    cc <- res[[run]]$caic[[stage]]
    cc[which.min(cc)]
}

getSummary <- function(res, stage=NULL) {
    est <- getEst(res, stage=stage)
    est <- est[,colSums(abs(est), na.rm=TRUE) > 0]
    fr <- colMeans(abs(est) > 0, na.rm=TRUE)
    cf <- colMeans(est, na.rm=TRUE)
    se <- apply(est, 2, sd, na.rm=TRUE)
    z <- cf/se
    p <- 2 * pnorm(-abs(z))
    cmat <- cbind(cf, se, fr, z, p)
    colnames(cmat) <- c("Estimate", "Std. Error", "Freq.", "z value", "Pr(>|z|)")
    cmat
}
#printCoefmat(getSummary(res))

getVcov <- function(res) {
    est <- getEst(res)
    est <- est[,colSums(abs(est), na.rm=TRUE) > 0]
    cov(est)
}
getConfint <- function(res, level=0.95, type=c("tboot","quantile")) {
    type <- match.arg(type)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    s <- getSummary(res)
    parm <- rownames(s)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%", sep="")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
    if (type == "tboot") {
        fac <- qnorm(a)
        ci[] <- s[,1] + s[,2] %o% fac
    } else {
        est <- getEst(res)
        est <- est[,colSums(abs(est)) > 0]
        cii <- t(apply(est, 2, quantile, probs=a))
        rownames(cii) <- parm
        ci[] <- cii[parm,,drop=FALSE]
    }
    return(ci)
}
#getConfint(res, type="tboot")

getMidPure <- function(res, mods) {
    OK <- !sapply(res, inherits, "try-error")
    mid <- data.frame(t(sapply(res[OK], function(z) unlist(z$mid))))
    rownames(mid) <- sapply(res[OK], "[[", "iteration")
    colnames(mid) <- names(mods)
    mid
}
getMid <- function(res, mods, use_rmax=FALSE) {
    OK <- !sapply(res, inherits, "try-error")
    mid <- data.frame(t(sapply(res[OK], function(z) unlist(z$mid))))
    rownames(mid) <- sapply(res[OK], "[[", "iteration")
    colnames(mid) <- names(mods)
    mid0 <- do.call(cbind,
        lapply(1:ncol(mid), function(i) 
            apply(mid[,1:i,drop=FALSE], 1, paste, collapse=" ")))
    if (use_rmax)
        require(ade4)
    Rao <- sapply(1:ncol(mid), function(i) {
        est <- suppressWarnings(getEst(res, stage=i, na.out=FALSE))
        est0 <- ifelse(est>0, 1, 0)
    #    d <- dist(est0, "manhattan") / ncol(est)
        d <- dist(est0)
        rmax <- if (use_rmax)
            suppressWarnings(divcmax(d))$value else 1
        PI <- data.matrix(rep(1/nrow(est), nrow(est)))
        drop(t(PI) %*% as.matrix(d) %*% PI) / rmax
    })
    attr(Rao, "use_rmax") <- use_rmax
    attr(mid, "Rao") <- Rao
    Gini <- 1 - apply(mid0, 2, function(z) sum((table(z)/length(z))^2))
    attr(mid, "Gini") <- Gini
    mid
}

getMid1 <- function(allres, mods) {
    mid <- t(sapply(allres, function(z) z[[1]]$mid))
    colnames(mid) <- names(mods)
    data.frame(mid)
}
#amid <- getMid1(allres, mods)

getMid2 <- function(allres, mods) {
    f <- function(res, mods) {
        OK <- !sapply(res, inherits, "try-error")
        mid <- t(sapply(res[OK], function(z) unlist(z$mid)))
        mid
    }
    mid <- lapply(allres, f, mods=mods)
    out <- t(sapply(mid, function(z) colSums(ifelse(z>0,1,0))))
    colnames(out) <- names(mods)
    out
}
#mid2 <- getMid2(allres,mods)

getFancyMid <- function(res, mods, allspp=FALSE) {
    mid <- if (allspp)
        getMid1(res, mods) else getMidPure(res, mods)
    out <- lapply(mid, function(z) data.frame(table(Prop=z)))
    for (i in 1:length(out)) {
        tmp <- out[[i]]
        rownames(tmp) <- tmp$Prop
        Terms <- c("NULL", sapply(mods[[i]], function(z) as.character(z)[3]))
        names(Terms) <- 0:(length(Terms) - 1)
        tmp$Terms <- Terms[match(tmp$Prop, names(Terms))]
        tmp$Prop <- round(tmp$Freq/sum(tmp$Freq), 4)
        out[[i]] <- tmp
    }
    out
}
#getFancyMid(res, mods)
#getFancyMid(allres, mods, allspp=TRUE)
getFancyMidTab <- function(res, mods, ...) {
    x <- getFancyMid(res, mods, ...)
    for (i in 1:length(x)) {
        x[[i]] <- data.frame(Stage=rep(names(x)[i], nrow(x[[i]])),
            Model_ID=rownames(x[[i]]), 
            Full_ID=paste(i, rownames(x[[i]]), sep="."),
            x[[i]])
    }
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
}
getFancyModsTab <- function(mods, chmax=60) {
    out <- mods
    for (i in 1:length(out)) {
        Terms <- c("NULL", sapply(mods[[i]], function(z) as.character(z)[3]))
#        if (any(nchar(as.character(Terms)) > chmax))
#            Terms <- compress(Terms)
        tmp <- data.frame(Stage=rep(names(mods)[i], length(Terms)),
            Model_ID=0:(length(Terms) - 1), 
            Full_ID=paste(i, 0:(length(Terms) - 1), sep="."),
            Terms=Terms)
        out[[i]] <- tmp
    }
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    out
}

## Lorenz info -- this is independent of mid selection freq !!!
getLc <- function(res, type=c("stats", "hilo","lam")) {
    type <- match.arg(type)
    if (type == "stats") {
       out <- t(sapply(res, function(z) unlist(z$lc)))
    }
    if (type == "hilo") {
        tmp <- table(unlist(lapply(res, "[[", "hi")))
        out <- numeric(nlevels(samp(mm)$HabitatA3))
        names(out) <- levels(samp(mm)$HabitatA3)
        out[names(tmp)] <- tmp
    }
    if (type == "lam") {
        out <- t(sapply(res, function(z) z$habmod))
        out <- out[,colnames(out) != "ROAD01"]
        out[,-1] <- out[,1] + out[,-1]
        out <- exp(out)
        colnames(out) <- levels(samp(mm)[[ip_name]])
    }
    out
}
#res <- allres[["CAWA"]]
#summary(getLc(res, "stats"))
#data.frame(getLc(res, "hilo"))

plotLc <- function(res, allhab=TRUE) {
    if (is.character(res))
        res <- loadSPP(res)
    OK <- !sapply(res, inherits, "try-error")
    resx <- res[OK]
    op <- par(mfrow=c(1,1), las=1, mar=c(5, 8, 4, 2) + 0.1)

    all <- rep(0, nlevels(samp(mm)$hab1ec))
    names(all) <- levels(samp(mm)$hab1ec)
    x <- getLc(resx, type="hilo") / nrow(getEst(res,na.out=FALSE))
    all[names(x)] <- x
    if (allhab)
        x <- all
    o <- order(x, decreasing=FALSE)
    x <- x[o]
    x0 <- 1-x
    barplot(rbind(x,x0), space=0, horiz=TRUE, xlab="Frequency",
        main=as.character(taxa(mm)[spp,"English_Name"]))
    abline(v=c(0.25, 0.5, 0.75), col="white")
    axis(1)
    #z <- getLc(resx, type="lam")
    #boxplot(z[,o], range=0, axes=FALSE, horizontal=TRUE, col="grey")
    #box()
    #axis(1)
    #axis(2, labels=names(x), at=1:length(x))
    #title(xlab="Mean")
    par(op)
    invisible(all)
}

## bootstrap based model support

midfig <- function(mid, m=apply(mid, 2, max), ...) {
    k <- ncol(mid)
    n <- nrow(mid)
    ylim <- c(0, k+1)
    xlim <- c(-(2+max(m)/2), (2+max(m)/2))
    pt <- lapply(1:k, function(i) {
        tt <- table(mid[,i])
        tt <- tt[match(0:m[i], names(tt))]
        tt[is.na(tt)] <- 0
        names(tt) <- 0:m[i]
        tt/n
    })
    yt <- lapply(1:k, function(i) rep(i, m[i]+1))
    xt <- lapply(1:k, function(i) 0:m[i] - m[i]/2)
    wl <- lapply(1:(k-1), function(i) {
        tt <- table(interaction(mid[,i], mid[,i+1]))
        tt[tt==0] <- NA
        tt/n
    })
    yl <- lapply(1:(k-1), function(i) {
        matrix(c(i, i+1), length(wl[[i]]), 2, byrow=TRUE)
    })
    xl <- lapply(1:(k-1), function(i) {
        z <- strsplit(names(wl[[i]]), "\\.")
        z <- matrix(as.integer(unlist(z)), length(z), 2, byrow=TRUE)
        z[,1] <- z[,1] - m[i]/2
        z[,2] <- z[,2] - m[i+1]/2
        z
    })
    wl <- do.call(c, wl)
    yl <- do.call(rbind, yl)
    xl <- do.call(rbind, xl)
    rr <- order(wl)
    wl <- wl[rr]
    yl <- yl[rr,]
    xl <- xl[rr,]
    CEX <- ifelse(unlist(pt)>0,1,0)*1+2.5*unlist(pt)/max(unlist(pt))
#    CEX <- 2+2.5*unlist(pt)/max(unlist(pt))
    LWD <- 3+5*wl/max(wl,na.rm=TRUE)
    COL <- grey(ifelse(is.na(wl), 0, wl/max(wl,na.rm=TRUE)))
    plot(xlim, ylim, type="n", xlab="", ylab="", axes=FALSE, ...)
    segments(-m/2, 1:k, m/2, 1:k)
    segments(xl[,1], yl[,1],xl[,2], yl[,2], 
        col=ifelse(is.na(wl), NA, 1), lwd=LWD)
    segments(xl[,1], yl[,1],xl[,2], yl[,2], 
        col=ifelse(is.na(wl), NA, COL), lwd=LWD-3)
    points(unlist(xt), unlist(yt), 
        cex=CEX,
        pch=19, col=grey(unlist(pt)/max(unlist(pt))))
    points(unlist(xt), unlist(yt), 
        cex=CEX,
        pch=21)
    text(unlist(xt), unlist(yt), 
        unlist(lapply(m, function(z) 0:z)))
    invisible(NULL)
}

getModelVariation <- function(res, mods, use_rmax=FALSE) {
    mid <- getMid(res, mods, use_rmax)
    out <- cbind(Gini=attr(mid, "Gini"),
        Rao=attr(mid, "Rao"))
    rownames(out) <- names(mods)
    out
}

plotMid <- function(res, mods, web=TRUE, ...) {
    if (is.character(res))
        res <- allres[[res]]
    mid <- getMid(res, mods, use_rmax=FALSE)
    if (web) {
        opar <- par(mai=0.1*c(1,1,2,1))
        midfig(mid, sapply(mods, length), 
            #main=paste(taxa(mm)[res[[1]]$species,"CommonName"], " (", res[[1]]$species, ")", sep="")
            main=as.character(taxa(mm)[res[[1]]$species,"English_Name"])
            , ...)
        text(rep(-0.5-0.5*max(sapply(mods, length)), ncol(mid)), 
            seq_len(ncol(mid))+0.25, names(mods))
        par(opar)
    } else {
        rc <- c(2,2)
        if (ncol(mid) != 4)
            rc[1] <- ceiling(ncol(mid)/2)
        opar <- par(mfrow=rc, mai=0.5*c(1,5,1,1), las=1)
        for (i in 1:ncol(mid)) {
            tmp <- table(mid[,i])
            aa <- rep(0, length(mods[[i]]) + 1)
            names(aa) <- 0:length(mods[[i]])
            aa[names(tmp)] <- tmp
            names(aa) <- c(".", mods[[i]])
            col <- grey(aa/sum(aa))
            barplot(rev(aa), main=names(mods)[i], col=rev(col), horiz=TRUE)
        }
        par(opar)
    }
    invisible(NULL)
}

## gof B
getGof <- function(spp, stage) {
    off <- TRUE
    res <- loadSPP(spp)
    off0 <- if (off)
        OFF[, spp] else 0
    if (missing(stage))
        stage <- length(mods)-1
    mu_hat <- getDataPred(res, stage, Xn, remn=TRUE)
    lam_hat <- exp(mu_hat + off0)
    Y <- as.numeric(xtab(mm)[,spp])
    #E <- pbapply(lam_hat, 1, median, na.rm=TRUE)
    E <- rowMeans(lam_hat, na.rm=TRUE)
    df0 <- data.frame(y=Y, lam=E, out=iout, p=dpois(Y, E))

    rval <- list()
    for (out in c(TRUE, FALSE)) {
        df <- if (out)
            df0[df0$out,] else df0[!df0$out,]
        df$out <- NULL
            
        yy2 <- aggregate(rep(1/nrow(df), nrow(df)), list(y=df$y), sum)
        yy2$p <- NA
        pp <- sapply(yy2$y, function(z) dpois(z, df$lam))
        pp <- pp / rowSums(pp)
        pp[is.na(pp)] <- 0
        yy2$p <- colMeans(pp)
        yy2$cx <- cumsum(yy2$x)
        yy2$cp <- cumsum(yy2$p)

        dfr <- max(yy2$y) + 1 - 1
        N <- nrow(df)
        level <- 0.9
        chi2 <- sum((yy2$p*N - yy2$x*N)^2 / (yy2$p*N))
        qchi <- qchisq(level, df=dfr, lower.tail=FALSE)
        pchi <- pchisq(chi2, df=dfr, lower.tail=FALSE)
        chi <- c(Chi2=chi2, df=dfr, qchi=qchi, level=level, pchi=pchi, N=N)

        rval[[ifelse(out, "out", "in")]] <- list(
            p=yy2, 
            pred=df, 
            off=off, 
            out=out,
            chi=chi)
    }
    rval
}

getGof2 <- function(spp, cutoff=0) {
    library(pROC)

    cat(spp, "\n");flush.console()
    x0 <- getGof(spp, 0)
    x4 <- getGof(spp, 4)
    x5 <- getGof(spp, 5)
    x6 <- getGof(spp, 6)

    rval <- list()
    for (out in c(TRUE, FALSE)) {
    
        y01 <- ifelse(as.numeric(xtab(mm)[,spp])>cutoff, 1, 0)
        y01 <- if (out)
            y01[iout] else y01[!iout]

        xx0 <- x0[[ifelse(out, "out", "in")]]
        xx4 <- x4[[ifelse(out, "out", "in")]]
        xx4 <- x4[[ifelse(out, "out", "in")]]
        xx5 <- x5[[ifelse(out, "out", "in")]]
        xx6 <- x6[[ifelse(out, "out", "in")]]
        
        cat("\t", ifelse(out, "out", "in"), "\troc 0\t");flush.console()
        roc0 <- roc(y01, xx0$pred$lam)
        cat("roc 4\t");flush.console()
        roc4 <- roc(y01, xx4$pred$lam)
        cat("roc 5\t");flush.console()
        roc5 <- roc(y01, xx5$pred$lam)
        cat("roc 6\n");flush.console()
        roc6 <- roc(y01, xx6$pred$lam)

        rval[[ifelse(out, "out", "in")]] <- list(
            xx0=xx0, xx4=xx4, xx5=xx5, xx6=xx6,
            roc0=roc0, roc4=roc4, roc5=roc5, roc6=roc6)
    }
    rval$cutoff <- cutoff
    rval
}

plotGof <- function(spp) {

    tmp <- if (is.character(spp))
        getGof2(spp) else tmp

    op <- par(mfrow=c(2,3), las=1, cex=1.2)

    for (out in c(FALSE, TRUE)) {

        xx0 <- tmp[[ifelse(out, "out", "in")]]$xx0
        xx4 <- tmp[[ifelse(out, "out", "in")]]$xx4
        xx5 <- tmp[[ifelse(out, "out", "in")]]$xx5
        xx6 <- tmp[[ifelse(out, "out", "in")]]$xx6

        dd0 <- xx0$pred
        dd4 <- xx4$pred
        dd5 <- xx5$pred
        dd6 <- xx6$pred
        yy0 <- xx0$p
        yy4 <- xx4$p
        yy5 <- xx5$p
        yy6 <- xx6$p

        roc001 <- tmp[[ifelse(out, "out", "in")]]$roc0
        roc401 <- tmp[[ifelse(out, "out", "in")]]$roc4
        roc501 <- tmp[[ifelse(out, "out", "in")]]$roc5
        roc601 <- tmp[[ifelse(out, "out", "in")]]$roc6


        y <- yy0$y
        ys <- cumsum((yy0$x+0.1) / sum(yy0$x+0.1))
        ys <- ys - ys[1]
        llam0 <- aggregate(dd0$lam, list(dd0$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        llam4 <- aggregate(dd4$lam, list(dd4$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        llam5 <- aggregate(dd5$lam, list(dd5$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        llam6 <- aggregate(dd6$lam, list(dd6$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        lp0 <- aggregate(dd0$p, list(dd0$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        lp4 <- aggregate(dd4$p, list(dd4$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        lp5 <- aggregate(dd5$p, list(dd5$y), quantile, c(0.5, 0.05, 0.95))[,-1]
        lp6 <- aggregate(dd6$p, list(dd6$y), quantile, c(0.5, 0.05, 0.95))[,-1]

        Colb <- c("#00FF0015", "#DD000010", "#D0D00025", "#80008015", "#0000DD10")[c(1,2,5,4)]
        Coll <- c(rgb(0.2,0.6,0.2), "red3", "orange3", "#A000A0", "blue3")[c(1,2,5,4)]


        ## ranking
        plot(ys, llam0[,1], 
            xlim=c(0, max(ys)), ylim=c(0, max(llam0[,1],llam4[,1],llam6[,1])),
            xlab="Observed count", ylab="Expected abundance", 
            #main=spp, 
            main=paste(as.character(taxa(mm)[spp,"English_Name"]),
                ifelse(out, "(external validation)", "(internal validation)")),
            axes=FALSE, type="b", pch=3, col=Coll[1])
        polygon(c(ys, rev(ys)), c(llam0[,2], rev(llam0[,3])), col=Colb[1], border=NA)
        polygon(c(ys, rev(ys)), c(llam4[,2], rev(llam4[,3])), col=Colb[2], border=NA)
        polygon(c(ys, rev(ys)), c(llam5[,2], rev(llam5[,3])), col=Colb[3], border=NA)
        polygon(c(ys, rev(ys)), c(llam6[,2], rev(llam6[,3])), col=Colb[4], border=NA)
        lines(ys, llam4[,1], col=Coll[2], type="b", pch=4)
        lines(ys, llam5[,1], col=Coll[3], type="b", pch=4)
        lines(ys, llam6[,1], col=Coll[4], type="b", pch=21)
        axis(2)
        axis(1, at=ys, labels=yy0$y)
        box()

        ## CDF
        lim <- range(yy0[,4:5], yy4[,4:5], yy6[,4:5])
        plot(yy0$cx, yy0$cp, ylim=lim, xlim=lim,
            xlab="Empirical CDF", ylab="Fitted CDF",
            type="b", pch=3, col=Coll[1])
        abline(0,1,lty=1, col="grey")
        lines(yy4$cx, yy4$cp, col=Coll[2], type="b", pch=4)
        lines(yy5$cx, yy5$cp, col=Coll[3], type="b", pch=2)
        lines(yy6$cx, yy6$cp, col=Coll[4], type="b", pch=21)
        legend("bottomright", bty="n", pch=c(3,4,2,21), lty=1, col=Coll,
            legend=c("Constant","Habitat","Hab + Clim","Hab + Clim + QS"))

        ## ROC / AUC
        plot(roc001, main="",
            lwd=1, xlim=c(1,0), ylim=c(0,1), col=Coll[1])
        lines(roc401, col=Coll[2], lwd=1)
        lines(roc501, col=Coll[3], lwd=1)
        lines(roc601, col=Coll[4], lwd=1)
        auc <- c(as.numeric(roc001$auc), as.numeric(roc401$auc), 
            as.numeric(roc501$auc), as.numeric(roc601$auc))
        names(auc) <- c("Constant","Habitat","Hab + Clim","Hab + Clim + QS")
        txt <- paste0(names(auc), " (AUC = ", round(auc, 3), ")")
        legend("bottomright", bty="n", col=Coll,
            lty=1, lwd=1, legend=txt)

    }
    par(op)
    invisible(tmp)
}



## ---------------- predict: natural veg and age [AUPRF]

predStat <- function(X, est, level=0.95, n=0, ci=TRUE, raw=FALSE) {
    tX <- t(X)
    mu <- apply(est, 1, function(z) crossprod(tX, z))
    if (raw)
        return(mu)
    pr <- exp(mu)
    out <- matrix(NA, nrow(pr), ifelse(ci, 6, 2))
    rownames(out) <- rownames(X)
    out[,1] <- rowMeans(pr, na.rm=TRUE)
    out[,2] <- apply(pr, 1, median, na.rm=TRUE)
    if (!ci) {
        colnames(out) <- c("Mean", "Median")
    } else {
        colnames(out) <- c("Mean", "Median", "CL1n", "CL2n", "CL1q", "CL2q")
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        if (n > 0) {
            tmp <- exp(replicate(n, rnorm(nrow(mu), rowMeans(mu), apply(mu, 1, sd))))
            out[,3:4] <- t(apply(tmp, 1, quantile, probs=a, na.rm=TRUE))
        }
        out[,5:6] <- t(apply(pr, 1, quantile, probs=a, na.rm=TRUE))
    }
    out
}

getMatVegAge <- function(est, 
clim.null=FALSE, fancy.labels=TRUE, add.hf=FALSE, 
hfor=0, hard=0, soft=0) 
{
    if (PROJECT == "SRD" && add.hf)
        stop("no local HF in SRD models")
    nn <- if (add.hf)
        50 else 48
    XX <- matrix(0, nn, ncol(est))
    colnames(XX) <- colnames(est)
    if (PROJECT == "JOSM")
        cn <- c("Habitat","HabitatB","isHForC",
            "Age","Age2","isWetConif","isMix","isPine","isUplConif")
    if (PROJECT == "AUPRF")
        cn <- c("Habitat","HabitatB","isHFor",
            "HabitatC1","HabitatC2","HabitatC3",
            "Age","Age2","isWetConif","isMix","isPine","isUplConif")
    if (PROJECT == "SRD")
        cn <- c("HabitatBF",
            "Age","Age2","isWetConif","isMix","isPine","isUplConif")
    xx <- xn[1:nn,cn]
    if (fancy.labels) {
        For <- paste(rep(c("Upland Spruce", "Pine", "Mixedwood", 
            "Deciduous", "Lowland Conifer"), each=9),
            rep(0:8*20, 5))
        rownames(xx) <- if (add.hf) {
            c(For, "Grass", "Shrub", "Wetland", "Urban-Industrial", "Cultivation") 
        } else c(For, "Grass", "Shrub", "Wetland")
    } else {
        For <- paste(rep(c("Conif", "Pine", "Mixed", 
            "Decid", "WetConif"), each=9),
            rep(1:9, 5), sep="")
        rownames(xx) <- if (add.hf) {
            c(For, "Grass", "Shrub", "Wet", "UrbInd", "Cult")
        } else c(For, "Grass", "Shrub", "Wet")
    }
    tmp <- c(rep(c("Conif", "Pine", "Mixed", 
        "Decid", "WetConif"), each=9), "Grass", "Shrub", "Wet")
    if (add.hf)
        tmp <- c(tmp, "UrbInd", "Cult")
    if (PROJECT %in% c("JOSM", "AUPRF")) {
        xx$Habitat[1:nn] <- tmp
        xx$HabitatB[1:nn] <- tmp
    }
    if (PROJECT == "SRD") {
        xx$HabitatBF[1:nn] <- tmp
    }
    tmp <- c(rep(0:8*20, 5), 0, 0, 0) / 200
    if (add.hf)
        tmp <- c(tmp, 0, 0)
    xx$Age <- tmp
    xx$Age2 <- xx$Age^2
    xx$isWetConif <- ifelse(xx$Habitat %in% c("WetConif"), 1, 0)
    xx$isMix <- ifelse(xx$Habitat %in% c("Mixed"), 1, 0)
    xx$isPine <- ifelse(xx$Habitat %in% c("Pine"), 1, 0)
    xx$isUplConif <- ifelse(xx$Habitat %in% c("Conif"), 1, 0)
    tmp <- c(rep(c(rep("A", 2), rep("B", 7)), 5), rep("", 3))
    if (add.hf)
        tmp <- c(tmp, "", "")

    if (PROJECT == "AUPRF") {
        xx$HabitatC1 <- factor(paste(xx$Habitat, tmp, sep=""), levels=levels(xx$HabitatC1))
        tmp <- c(rep(c(rep("A", 4), rep("B", 5)), 5), rep("", 3))
        if (add.hf)
            tmp <- c(tmp, "", "")
        xx$HabitatC2 <- factor(paste(xx$Habitat, tmp, sep=""), levels=levels(xx$HabitatC2))
        tmp <- c(rep(c(rep("A", 2), rep("B", 2), rep("C", 5)), 5), rep("", 3))
        if (add.hf)
            tmp <- c(tmp, "", "")
        xx$HabitatC3 <- factor(paste(xx$Habitat, tmp, sep=""), levels=levels(xx$HabitatC3))
        xx$isHFor <- hfor
    }
    if (PROJECT == "JOSM") {
        MAXFOR <- 50/200
        xx$isHForC <- hfor * pmax(0, 1 - xx$Age/MAXFOR)
    }

    xx$ROAD <- hard
    xx$SoftLin_PC <- soft
#    xx$xlat <- 0
#    xx$xlong <- 0
    xx$xlat <- 0.9059 # mean(xyqs$xlat)
    xx$xlong <- 1.2012 # mean(xyqs$xlong)
    if (clim.null) {
        xx$xMAP <- 0
        xx$xPET <- 0
        xx$xCMD <- 0
        xx$xMAT <- 0
    } else {
        xx$xMAP <- 0.2109 # mean(xyqs$xMAP)
        xx$xPET <- 0.6483 # mean(xyqs$xPET)
        xx$xCMD <- 0.3128 # mean(xyqs$xCMD)
        xx$xMAT <- 0.1062 # mean(xyqs$xMAT)
    }
    XX <- model.matrix(~(.)^2, xx)
    colnames(XX) <- fixNames(colnames(XX))
    XX1 <- XX[,intersect(colnames(est), colnames(XX))]
    XX2 <- Xn[1:nrow(XX1),setdiff(colnames(est), colnames(XX))]
    XX2[] <- 0
    XX <- cbind(XX1, XX2)[,colnames(est)]
    XX
}

predVegAge <- function(res, level=0.95, stage, fancy.labels=TRUE) {
    if (missing(stage)) {
#        stage <- switch(PROJECT,    # veg,age,road,space
#            "AUPRF" = 4,
#            "SRD" = 1)
        stage <- switch(PROJECT,    # veg,age,road,space,QSHF
            "AUPRF" = 6,
            "SRD" = 3,
            "JOSM" = 3) # no QS and Space considered -- marginal hab/age
    }
    est <- getEst(res, stage=stage, na.out=TRUE) 
    est[,grepl("_QS", colnames(est))] <- 0
    out <- predStat(getMatVegAge(est, clim.null=FALSE, fancy.labels=fancy.labels), 
        est, level=level, n=0)
    ## burn
    pburn <- NULL

    OK <- !sapply(res, inherits, "try-error")
    mid <- data.frame(t(sapply(res[OK], function(z) unlist(z$mid))))
#    mid <- data.frame(t(sapply(res, function(z) unlist(z$mid))))
    if (any(mid[,1] %in% c(2,4))) {
        pburn <- predStat(matrix(1,2,2), 
            est[,c("(Intercept)","HabitatBBurn")], level=level, n=0)
        pburn <- pburn[1,]
    }
    attr(out, "burn") <- pburn
    out
}

plotVegAge <- function(res, level=0.95, ylim, srd=FALSE, stage) {
    if (is.character(res))
        res <- allres[[res]]
    x <-  predVegAge(res, level=level, stage=stage)
    rn <- gsub("[[:digit:]]", "", rownames(x))
    rn <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", rn, perl=TRUE)
    if (srd)
        x <- 1-exp(-x)
    if (srd)
        ylim <- c(0,1)
    pburn <- attr(x, "burn")
    op <- par(mar = c(7, 4, 4, 2) + 0.1)
    if (missing(ylim)) {
        maxMed <- max(x[,"Median"], pburn["Median"], na.rm=TRUE)
        maxCL <- max(x[,"CL2q"], pburn["CL2q"], na.rm=TRUE)
        ylim <- c(0, min(1.5*maxMed, maxCL))
    }
    xval <- x[,"Median"]
    col <- c(rep(brewer.pal(9,"Greens"), 5), brewer.pal(5,"Accent")[3:5])
    if (!is.null(pburn)) {
        xval <- c(xval, pburn["Median"])
        col <- c(col, "red3")
    }
    barplot(xval, beside=TRUE, space=0,
        col=col, 
        ylim=ylim,
        axisnames=FALSE)
    axis(1, at=c(1:5 * 9 - 4.5), labels=unique(rn)[1:5], tick=FALSE)
#    axis(1, at=c(46:48), labels=unique(rn)[6:8], tick=FALSE, srt=45)
    if (is.null(pburn)) {
        text(c(46:48)-0.5, par("usr")[3] - (0.04 * ylim[2]), srt = 90, adj = 1,
              labels = unique(rn)[6:8], xpd = TRUE)
    } else {
        text(c(46:49)-0.5, par("usr")[3] - (0.04 * ylim[2]), srt = 90, adj = 1,
              labels = c(unique(rn)[6:8], "Burn"), xpd = TRUE)
    }
    title(main=paste(taxa(mm)[res[[1]]$species,"CommonName"], " (", res[[1]]$species, ")", sep=""),
        ylab="Relative Abundance")
    segments(1:nrow(x) - 0.5, x[,"CL1q"], 1:nrow(x) - 0.5, x[,"CL2q"], lwd=2, col="darkgrey")
    segments(1:nrow(x) - 0.4, x[,"CL1q"], 1:nrow(x) - 0.6, x[,"CL1q"], lwd=2, col="darkgrey")
    segments(1:nrow(x) - 0.4, x[,"CL2q"], 1:nrow(x) - 0.6, x[,"CL2q"], lwd=2, col="darkgrey")
    if (!is.null(pburn)) {
        segments(nrow(x)+1 - 0.5, pburn["CL1q"], nrow(x)+1 - 0.5, pburn["CL2q"], 
            lwd=2, col="darkgrey")
        segments(nrow(x)+1 - 0.4, pburn["CL1q"], nrow(x)+1 - 0.6, pburn["CL1q"], 
            lwd=2, col="darkgrey")
        segments(nrow(x)+1 - 0.4, pburn["CL2q"], nrow(x)+1 - 0.6, pburn["CL2q"], 
            lwd=2, col="darkgrey")
    }
    par(op)
    invisible(NULL)
}

predLocalHF <- function(res, level=0.95, stage) {
    if (missing(stage)) {
        stage <- switch(PROJECT,    # veg,age,road,space,QSHF
            "AUPRF" = 6,
            "SRD" = 3,
            "JOSM" = 3) # no QS and Space considered -- marginal hab/age
    }
    est <- getEst(res, stage=stage) # veg,age,road,space
    XX <- matrix(0, 47, ncol(est))
    colnames(XX) <- colnames(est)
    xx <- xn[1:47,c("Habitat","HabitatB","isHForC",
#        "HabitatC1","HabitatC2","HabitatC3",
        "Age","Age2","isWetConif","isMix","isPine","isUplConif")]
    For <- paste(rep(c("Upland Spruce", "Pine", "Mixedwood", 
        "Deciduous", "Lowland Conifer"), each=9),
        rep(0:8*20, 5))
    rownames(xx) <- c(For, "Urban-Industrial", "Cultivation")
    xx$Habitat[1:47] <- c(rep(c("Conif", "Pine", "Mixed", 
        "Decid", "WetConif"), each=9), "UrbInd", "Cult")
    xx$HabitatB[1:47] <- xx$Habitat[1:47]
    xx$Age <- c(rep(0:8*20, 5), 0, 0) / 200
    xx$Age2 <- xx$Age^2
    xx$isWetConif <- ifelse(xx$Habitat %in% c("WetConif"), 1, 0)
    xx$isMix <- ifelse(xx$Habitat %in% c("Mixed"), 1, 0)
    xx$isPine <- ifelse(xx$Habitat %in% c("Pine"), 1, 0)
    xx$isUplConif <- ifelse(xx$Habitat %in% c("Conif"), 1, 0)

    MAXFOR <- 50/200
    xx$isHForC <- pmax(0, 1 - xx$Age/MAXFOR)
    xx$ROAD <- 0
    xx$SoftLin_PC <- 0
    xx$xlat <- 0.9059 # mean(xyqs$xlat)
    xx$xlong <- 1.2012 # mean(xyqs$xlong)
    xx$xMAP <- 0.2109 # mean(xyqs$xMAP)
    xx$xPET <- 0.6483 # mean(xyqs$xPET)
    xx$xCMD <- 0.3128 # mean(xyqs$xCMD)
    xx$xMAT <- 0.1062 # mean(xyqs$xMAT)

    XX <- model.matrix(~(.)^2, xx)
    colnames(XX) <- fixNames(colnames(XX))
    XX1 <- XX[,intersect(colnames(est), colnames(XX))]
    XX2 <- Xn[1:nrow(XX1),setdiff(colnames(est), colnames(XX))]
    XX2[] <- 0
    XX <- cbind(XX1, XX2)[,colnames(est)]
    predStat(XX, est, level=level, n=0)
}

## not for SRD
predLocalLin <- function(res, level=0.95, stage, nn=5000) {
    if (missing(stage)) {
        stage <- switch(PROJECT,    # veg,age,road
            "AUPRF" = 4,
            "JOSM" = 3)
    }
    est <- getEst(res, stage=stage) # veg,age,road,space

    Xnn <- Xn[!(xn$Habitat %in% c("Cult", "UrbInd")),]
    XX <- Xnn[sample.int(nrow(Xnn), nn),]
    XX3 <- XX2 <- XX1 <- XX
    XX[,"ROAD"] <- 0
    XX[,"SoftLin_PC"] <- 0
    XX1[,"ROAD"] <- 1 # mean(xn$ROAD)
    XX1[,"SoftLin_PC"] <- 0
    XX2[,"ROAD"] <- 0
    XX2[,"SoftLin_PC"] <- mean(xn$SoftLin_PC)
    XX3[,"ROAD"] <- 1 # mean(xn$ROAD)
    XX3[,"SoftLin_PC"] <- mean(xn$SoftLin_PC)
    if (PROJECT == "AUPRF") {
        XX[,"ROAD:SoftLin_PC"] <- 0
        XX1[,"ROAD:SoftLin_PC"] <- 0
        XX2[,"ROAD:SoftLin_PC"] <- 0
        XX3[,"ROAD:SoftLin_PC"] <- mean(xn$ROAD) * mean(xn$SoftLin_PC)
    }
    p0 <- predStat(XX, est, level=level, n=0)
    p1 <- predStat(XX1, est, level=level, n=0)
    p2 <- predStat(XX2, est, level=level, n=0)
    p3 <- predStat(XX3, est, level=level, n=0)
    rbind("None"=colMeans(p2), 
        "Hard Linear"=colMeans(p1), "Soft Linear"=colMeans(p2),
        "Both"=colMeans(p3))
}


## not for SRD
plotVegAgeHF <- function(res, level=0.95, ylim, srd=FALSE, stage) {
    if (is.character(res))
        res <- allres[[res]]
    x <-  predVegAge(res, level=level, stage=stage)
    pburn <- attr(x, "burn")
    x2 <- predLocalHF(res, level=level, stage=stage)
#    xl <- predLocalLinD(res)
    x <- rbind(x, x2[c("Urban-Industrial","Cultivation"),])
    x2 <- x2[1:45,]
    x3 <- x
    x3[1:45,] <- x2
#    xxl <- matrix(NA, 2, 6)
#    xxl[,2] <- xl[-1]
#    rownames(xxl) <- c("Hard Linear", "Soft Linear")
#    x <- rbind(x, xxl)
    rn <- gsub("[[:digit:]]", "", rownames(x))
    rn <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", rn, perl=TRUE)
    if (missing(ylim))
        ylim <- c(0, min(1.5*max(x[,"Median"], x2[,"Median"], pburn["Median"], na.rm=TRUE), 
            max(x[,"CL2q"], x3[,"CL2q"], pburn["CL2q"], na.rm=TRUE)))
    if (srd)
        x <- 1-exp(-x)
    if (srd)
        ylim <- c(0,1)
    op <- par(mar = c(9, 4, 4, 2) + 0.1)
    xval <- x[,"Median"]
    col <- c(rep(brewer.pal(9,"Greens"), 5), brewer.pal(7,"Accent")[3:7])
    if (!is.null(pburn)) {
        xval <- c(xval, pburn["Median"])
        col <- c(col, "red2")
    }
    barplot(xval, beside=TRUE, space=0,
        col=col, 
        ylim=ylim,
        axisnames=FALSE)
    axis(1, at=c(1:5 * 9 - 4.5), labels=unique(rn)[1:5], tick=FALSE)
#    axis(1, at=c(46:48), labels=unique(rn)[6:8], tick=FALSE, srt=45)
    if (is.null(pburn)) {
        text(c(46:50)-0.5, par("usr")[3] - (0.04 * ylim[2]), srt = 90, adj = 1,
            labels = unique(rn)[6:10], xpd = TRUE)
    } else {
        text(c(46:51)-0.5, par("usr")[3] - (0.04 * ylim[2]), srt = 90, adj = 1,
            labels = c(unique(rn)[6:10], "Burn"), xpd = TRUE)
    }
    title(main=paste(taxa(mm)[res[[1]]$species,"CommonName"], " (", res[[1]]$species, ")", sep=""),
        ylab="Relative Abundance")
    segments(1:nrow(x) - 0.5, x[,"CL1q"], 1:nrow(x) - 0.5, x[,"CL2q"], lwd=2, col="darkgrey")
    segments(1:nrow(x) - 0.4, x[,"CL1q"], 1:nrow(x) - 0.6, x[,"CL1q"], lwd=2, col="darkgrey")
    segments(1:nrow(x) - 0.4, x[,"CL2q"], 1:nrow(x) - 0.6, x[,"CL2q"], lwd=2, col="darkgrey")
    if (!is.null(pburn)) {
        segments(nrow(x)+1 - 0.5, pburn["CL1q"], nrow(x)+1 - 0.5, pburn["CL2q"], 
            lwd=2, col="darkgrey")
        segments(nrow(x)+1 - 0.4, pburn["CL1q"], nrow(x)+1 - 0.6, pburn["CL1q"], 
            lwd=2, col="darkgrey")
        segments(nrow(x)+1 - 0.4, pburn["CL2q"], nrow(x)+1 - 0.6, pburn["CL2q"], 
            lwd=2, col="darkgrey")
    }
    #tmp <- rep(c(TRUE, TRUE, FALSE, TRUE, FALSE), each=9)
    tmp <- c(1:5,10:14,28:32)
    lines((1:45 - 0.5)[1:5], x2[1:5,"Median"], col="tomato")
    lines((1:45 - 0.5)[10:14], x2[10:14,"Median"], col="tomato")
    lines((1:45 - 0.5)[28:32], x2[28:32,"Median"], col="tomato")
    points((1:45 - 0.5)[tmp], x2[tmp,"Median"], pch=4, col=2, cex=1)
    par(op)
    invisible(NULL)
}
#plotVegAgeHF(allres[["ALFL"]])

## not for SRD
plotLocalLin <- function(res, level=0.95, ylim, srd=FALSE, main, stage, nn=5000) {
    if (is.character(res))
        res <- allres[[res]]
    x <- predLocalLin(res, level=level, stage=stage, nn=nn)
    if (srd)
        x <- 1-exp(-x)
    rn <- rownames(x)
    if (missing(ylim))
        ylim <- c(0, min(1.5*max(x[,"Median"]), max(x[,"CL2q"])))
    if (srd)
        ylim <- c(0,1)
    barplot(x[,"Median"], beside=TRUE, space=0,
        col=brewer.pal(4,"Accent"), 
        ylim=ylim)
    if (missing(main))
        main <- paste(taxa(mm)[res[[1]]$species,"CommonName"], " (", res[[1]]$species, ")", sep="")
    ylab <- if (srd)
        "Relative Abundance" else "Relative Abundance"
    title(main=main, ylab=ylab)
    segments(1:nrow(x) - 0.5, x[,"CL1q"], 1:nrow(x) - 0.5, x[,"CL2q"], lwd=2, col="darkgrey")
    segments(1:nrow(x) - 0.4, x[,"CL1q"], 1:nrow(x) - 0.6, x[,"CL1q"], lwd=2, col="darkgrey")
    segments(1:nrow(x) - 0.4, x[,"CL2q"], 1:nrow(x) - 0.6, x[,"CL2q"], lwd=2, col="darkgrey")
    invisible(NULL)
}
#plotLocalLin(allres[["ALFL"]])


## QS level HF effect prediction

## repeated subsampling is required here with single pt
## to represent shift in intercept
predQS <- 
function(hft, est1, mu0, level=0.95, est0, single.pt=TRUE, prab=FALSE)
{
    if (hft=="IP") {
        cc <- "Remn_QS"
        cc2 <- "Remn2_QS"
    }
    if (hft=="Wet") {
        cc <- "pWet_QS"
        cc2 <- character(0)
    }
    if (hft=="WetWater") {
        cc <- "pWetWater_QS"
        cc2 <- character(0)
    }
    if (hft=="SoftLin") {
        cc <- c("Lin_QS", "Succ_QS", "Noncult_QS", "THF_QS")
        cc2 <- c("Succ2_QS", "Noncult2_QS", "THF2_QS")
    }
    if (hft=="HardLin") {
        cc <- c("Lin_QS", "Alien_QS", "Noncult_QS", "THF_QS")
        cc2 <- c("Alien2_QS", "Noncult2_QS", "THF2_QS")
    }
    if (hft=="UrbInd") {
        cc <- c("Nonlin_QS", "Alien_QS", "Noncult_QS", "THF_QS")
        cc2 <- c("Nonlin2_QS", "Alien2_QS", "Noncult2_QS", "THF2_QS")
    }
    if (hft=="Cult") {
        cc <- c("Nonlin_QS", "Alien_QS", "Cult_QS", "THF_QS")
        cc2 <- c("Nonlin2_QS", "Alien2_QS", "THF2_QS")
    }
    if (hft=="HFor") {
        cc <- c("Nonlin_QS", "Succ_QS", "Noncult_QS", "THF_QS")
        cc2 <- c("Nonlin2_QS", "Succ2_QS", "Noncult2_QS", "THF2_QS")
    }
    if (all(est1[,cc]==0))
        return(NULL)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    hfout <- seq(0, ifelse(hft %in% c("SoftLin","HardLin"), 0.24, 1), len=101)
    ## Xnqs already randomized
    HFout <- rep(hfout, ceiling(nrow(Xnqs) / length(hfout)))[1:nrow(Xnqs)]
    Xnqs2 <- Xnqs
    Xnqs2[,grepl("_QS", colnames(est1))] <- 0
    Xnqs2[,cc] <- HFout
    Xnqs2[,cc2] <- HFout^2
    if (single.pt) {
        mu0bar <- 0 # mean(mu0)
        Est <- est1
        Est[,!grepl("_QS", colnames(est1))] <- 0
        Mu <- predStat(Xnqs2, Est, raw=TRUE)
        eMu <- exp(Mu)
        if (hft %in% c("SoftLin","HardLin","Wet","WetWater")) {
            cb <- apply(eMu, 2, function(z) {
                m0 <- try(lm(log(z) ~ HFout - 1), silent=TRUE)
                c0 <- if (inherits(m0, "try-error"))
                    c(NA) else coef(m0)
                exp(mu0bar + hfout*c0[1])
            })
        } else {
            cb <- apply(eMu, 2, function(z) {
                m0 <- try(lm(log(z) ~ HFout + I(HFout^2) - 1), silent=TRUE)
                c0 <- if (inherits(m0, "try-error"))
                    c(NA, NA) else coef(m0)
                exp(mu0bar + hfout*c0[1] + hfout^2*c0[2])
            })
        }
    } else {
        mu0bar <- mean(mu0)
        Est <- est1
        Est[] <- 0
        Est[,cc] <- est1[,cc]
        Mu0 <- predStat(Xnqs, est0, raw=TRUE)
        MuHF <- predStat(Xnqs2, Est, raw=TRUE)
        eMu <- exp(Mu0 + MuHF)
        if (hft %in% c("SoftLin","HardLin","Wet","WetWater")) {
            cb <- apply(eMu, 2, function(z) {
                m0 <- try(lm(log(z) ~ HFout), silent=TRUE)
                c0 <- if (inherits(m0, "try-error"))
                    c(NA, NA) else coef(m0)
                exp(c0[1] + hfout*c0[2])
            })
        } else {
            cb <- apply(eMu, 2, function(z) {
                m0 <- try(lm(log(z) ~ HFout + I(HFout^2)), silent=TRUE)
                c0 <- if (inherits(m0, "try-error"))
                    c(NA, NA, NA) else coef(m0)
                exp(c0[1] + hfout*c0[2] + hfout^2*c0[3])
            })
        }
    }
    ci <- t(apply(cb, 1, quantile, probs=c(0.5, a), na.rm=TRUE))
    out <- cbind(x=hfout, y=ci[,1], ci[,-1])
    out[,"y"] <- rowMeans(cb)
    colnames(out)[3:4] <- c("cl1","cl2")
    if (diff(range(ci[,1])) < (exp(mu0bar) / 100)) {
        rval <- NULL
    } else {
        if (prab)
            out[,2:4] <- 1-exp(-out[,2:4])
    }
    rval <- list(summary=out, raw=cb)
    rval
}

predAllQS <- function(spp, level=0.95, stage, single.pt=TRUE, prab=FALSE) {
    if (missing(stage)) {
        stage <- switch(PROJECT,
            "AUPRF" = 6,
            "SRD" = 3,
            "JOSM" = 6) # no space involved -- marginal hab/age/contrast/QS
    }
    res <- allres[[spp]]
    est <- getEst(res, stage=stage, na.out=FALSE)
    est0 <- est1 <- est
    Cols <- colnames(est)[grepl("_QS", colnames(est))]

    est0[,Cols] <- 0
    est1[,!(colnames(est1) %in% Cols)] <- 0
    mu0 <- predStat(Xnqs, est0, raw=TRUE)
    mu0 <- apply(mu0, 1, median)
    hfts <- c("IP", "Wet","WetWater","HFor","SoftLin","HardLin","UrbInd","Cult")
    xx <- lapply(hfts, predQS, est1=est1, mu0=mu0, level=level, est0=est0,
        single.pt=single.pt, prab=prab)
    names(xx) <- hfts
    attr(xx, "args") <- list(mu0=mu0, est0=est0, est1=est1)
    xx
}

plotAllQS <- 
function(spp, level=0.95, ylim, single.pt=TRUE) 
{
    aa <- predAllQS(spp, level=level, single.pt=single.pt, prab=FALSE)
    aa0 <- attr(aa, "args")
    if (missing(ylim))
        ylim <- max(sapply(aa, function(z) {
            if (is.null(z$summary))
                return(0)
            max(z$summary[,-1], na.rm=TRUE)
        }))
    ii <- c("IP", "Wet","WetWater",
        "UrbInd","Cult","HardLin",
        "SPP","HFor","SoftLin")
    op <- par(mfrow=c(3,3), mar=c(4, 3, 2, 1) + 0.1)
    for (hft in ii) {
        if (hft=="SPP" || hft=="NONE" || is.null(aa[[hft]])) {
            plot(0,0,ann=FALSE,axes=FALSE,type="n")
            if (hft=="SPP")
                text(0,0,spp,cex=2)
            if (is.null(aa[[hft]]) && hft != "SPP" && hft != "NONE") {
                box()
                text(0,0,paste(hft, "no effect", sep=":\n"),cex=1)
            }
        } else {
            aa1 <- aa[[hft]]$raw
            aaa <- data.frame(aa[[hft]]$summary)
            plot(0, type="n", xlim=c(0,100), ylim=c(0,ylim), 
                ylab="", # ylab="Rel. abund.",
                xlab=hft)
            polygon(100*c(aaa$x, rev(aaa$x)), c(aaa[,3], rev(aaa[,4])), col="gold", border=NA)
            matlines(100*aaa$x, aa1, lty=1, col="grey")
            lines(y ~ I(100*x), aaa, type="l", lwd=3, col=2)
        }
    }
    par(op)
    invisible(aa)
}
#plotAllQS("CAWA",.9,single.pt=T)

## year effects

getYr <- function(res) {
    OK <- !sapply(res, inherits, "try-error")
    out <- t(sapply(res[OK], function(z) {
        cf <- z$coef[[length(z$coef)]]
        if (z$mid[length(z$mid)] == 0)
            out <- c(cf[1], 0, 0, 0, 0)
        if (z$mid[length(z$mid)] == 1)
            out <- c(cf[1], cf[grep("YEAR", names(cf))], 0, 0, 0)
        if (z$mid[length(z$mid)] == 2)
            out <- c(cf[1], 0, cf[grep("YR5F", names(cf))])
        out[1] <- 0
        out
    }))
    out
}

predYr <- function(res, ...) {
    Yr <- 1992:2013
    xYr <- (Yr-min(Yr))/10
    Xyr <- matrix(0, length(Yr), 5)
    Xyr[,1] <- 1
    Xyr[,2] <- xYr
    Xyr[Yr %in% 1998:2003,3] <- 1
    Xyr[Yr %in% 2004:2008,4] <- 1
    Xyr[Yr %in% 2009:2013,5] <- 1
    est <- getYr(res)
    predStat(Xyr, est, ...)
}

plotYr <- function(res, level=0.95, ...) {
    if (is.character(res))
        res <- allres[[res]]
    Yr <- 1992:2013
    pr <- predYr(res, level=level)
    plot(Yr, pr[,"Median"], ylim=c(0, max(pr[,c("Median","CL1q","CL2q")], na.rm=TRUE)), type="n",
        main=paste(taxa(mm)[res[[1]]$species,"CommonName"], " (", 
                res[[1]]$species, ")", sep=""),
        xlab="Year",ylab="Relative abundance", ...)
    polygon(c(Yr, rev(Yr)), 
        c(pr[,"CL1q"], rev(pr[,"CL2q"])), col="gold", border=NA)
    lines(Yr, pr[,"Median"], lwd=2, col="tomato")
    lines(Yr, pr[,"Mean"], lwd=1, lty=2, col="tomato")
    invisible(NULL)
}

