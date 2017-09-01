get_stuff0 <- function(TAX, TABLE, COL, north=TRUE) {
    SPP <- AllIn[[TAX]]$lt$Species
    i <- if (north)
        AllIn[[TAX]]$lt[,"veghf.north"] else AllIn[[TAX]]$lt[,"soilhf.south"]
    SPP <- as.character(SPP[i])
    z <- AllIn[[TAX]][[TABLE]]
    rownames(z) <- z$Species
    #compare_sets(rownames(z), SPP)
    z[SPP, COL]
}
rugged_mat <- function(x) {
    n <- sapply(x, length)
    out <- matrix(NA, max(n), length(x))
    colnames(out) <- names(x)
    for (i in 1:length(x))
        out[1:length(x[[i]]),i] <- x[[i]]
    out
}
get_stuff <- function(TABLE, COL, north=TRUE)
    rugged_mat(lapply(structure(names(AllIn), names=names(AllIn)),
        get_stuff0, TABLE=TABLE, COL=COL, north=north))
bounded_density <- function(x, interval=c(-Inf, Inf), ...) {
    require(intrval)
    interval <- sort(interval[1L:2L])
    if (!all(x %()% interval))
        stop("found values of x outside of open interval")
    a <- interval[1L]
    b <- interval[2L]
    bounds  <- is.finite(interval)
    if (!bounds[1L] && !bounds[2L]) { # (Inf, Inf)
        f <- finv <- function(x) x
    }
    if (bounds[1L] && !bounds[2L]) { # (a, Inf)
        f <- function(x)
            log(x-a)
        finv <- function(z)
            exp(z) + a
    }
    if (!bounds[1L] && bounds[2L]) { # (Inf, b)
        f <- function(x)
            log(b-x)
        finv <- function(z)
            b-exp(z)
    }
    if (bounds[1L] && bounds[2L]) { # (a, b)
        f <- function(x)
            qlogis((x-a) / (b-a))
        finv <- function(z)
            plogis(z) * (b-a) + a
    }
    fx <- f(x)
    d <- density(fx, ...)
    v <- d$x
    dv <- diff(v)
    h <- d$y
    n <- length(h)
    dh <- rowMeans(cbind(h[-n], h[-1L]))
    A <- dv * dh
    vinv <- finv(v)
    hinv <- A / abs(diff(vinv))
    hinv <- c(hinv[1L], hinv)
    data.frame(x=vinv, y=hinv)
}

cp1 <- function(x, p=c(0, 1), nmin=5, interval=c(-Inf, Inf), ...) {
    x2 <- x[!is.na(x)]
    q <- quantile(x2, p)
    x2 <- x2[x2 > q[1] & x2 < q[2]]
    if (length(unique(x2)) < nmin) {
        m <- mean(x2)
        s <- sd(x2)
        v <- seq(min(min(x2), m-4*s), max(max(x2), m+4*s), length.out = 512)
        v <- v %()% interval
        d <- dnorm(v, m, s)
        v <- c(v[1L], v, v[length(v)])
        d <- c(0, d, 0)
        out <- data.frame(h=v, w=d/max(d))
    } else {
        d <- bounded_density(x2, interval=interval, ...)
        out <- data.frame(h=d$x, w=d$y/max(d$y))
    }
    out
}
vp <- function(x, p=c(0,1),
col="#FF664D80", border="#FF664D",
#col="#40E0D180", border="#40E0D1",
ylim, ylab="", xlab="", main="", nmin=5, interval=c(-Inf, Inf), ...) {
    xx <- apply(x, 2, cp1, p=p, nmin=nmin, interval=interval, ...)
    if (missing(ylim))
        ylim <- range(x, na.rm=TRUE)
    plot(0, type="n", axes=FALSE, ann=FALSE, xlim=c(0.5, ncol(x)+0.5),
        ylim=ylim)
    axis(1, 1:ncol(x), colnames(x), lwd=0)
    axis(2)
    title(ylab=ylab, xlab=xlab, main=main)
    col <- rep(col, ncol(x))[1:ncol(x)]
    border <- rep(border, ncol(x))[1:ncol(x)]
    for (i in 1:ncol(x)) {
        if (!is.null(xx[[i]])) {
            polygon(0.45*c(-xx[[i]]$w, rev(xx[[i]]$w))+i,
                c(xx[[i]]$h, rev(xx[[i]]$h)), col=col[i], border=border[i])
            #lines(c(i,i), range(xx[[i]]$h), col=col[i], lwd=2)
        }
        lines(c(i-0.2, i+0.2), rep(median(x[,i], na.rm=TRUE), 2), lwd=3)
    }
    #points(1:ncol(x), colMeans(x, na.rm=TRUE), pch=21, cex=1.5)
    invisible(NULL)
}
fu <- function(x) sign(x) * plogis(log(abs(x/100)))
dfun <- function(pp) {
    ppp <- density(pp)
    list(x=ppp$xcol, y=ppp$yrow, z=t(ppp$v)/pp$n)
}
ord_fun <- function(TAX, north=TRUE, scaling=2, alpha=1, col.text=1, col.pts,
    cex.pts=1, cex.text=1, plot=TRUE, ...) {
    require(vegan)
    Excl <- c("WhiteSpruceCC20", "WhiteSpruceCC40", "WhiteSpruceCC60",
        "PineCC20", "PineCC40", "PineCC60",
        "DeciduousCC20",   "DeciduousCC40", "DeciduousCC60",
        "MixedwoodCC20", "MixedwoodCC40",   "MixedwoodCC60")
    if (TAX=="All") {
        TT <- if (north)
            "veg" else "soilnt"
        z <- do.call(rbind, lapply(AllIn, function(zz) zz[[TT]]))
    } else {
        z <- if (north)
            AllIn[[TAX]]$veg else AllIn[[TAX]]$soilnt
    }
    rownames(z) <- z[,1]
    z <- as.matrix(z[,-1])
    z <- z[,!grepl("\\.", colnames(z))]
    z <- z[,!(colnames(z) %in% Excl)]
    cn <- colnames(z)
    cn <- gsub("[[:digit:]]", "", cn)
    #cn <- gsub("CC", "", cn)
    zz <- groupMeans(z, 2, cn)
    #rownames(zz) <- make.cepnames(rownames(zz)) # species
    colnames(zz) <- make.cepnames(colnames(zz)) # habitats

    fs <- function(x) 0.5*(1+x/max(abs(x)))
    yy <- (t(zz))
    m <- cca(yy)
    s <- scores(m, 1:3)
    ColSi <- 2
    ColSp <- rgb(red=fs(s$species[,1]),
        green=fs(s$species[,2]), blue=fs(s$species[,3]), alpha=alpha)
    if (!missing(col.pts))
        ColSp <- col.pts
    if (plot) {
        plot(m, scaling=scaling, type="none", ...)
        #text(m, "species", col=ColSp, cex=0.5, scaling=scaling)
        points(m, "species", col=ColSp, cex=cex.pts, pch=19, scaling=scaling)
        text(m, "sites", col=col.text, cex=cex.text, scaling=scaling)
    }
    invisible(m)
}
