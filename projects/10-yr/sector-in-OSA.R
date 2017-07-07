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

if (FALSE) {
x <- rnorm(100)

op <- par(mfrow=c(2,2))
range(x)
d1 <- bounded_density(exp(x))
d2 <- bounded_density(exp(x), c(0, Inf))
plot(bounded_density(x, type="l"))

range(exp(x))
d1 <- bounded_density(exp(x))
d2 <- bounded_density(exp(x), c(0, Inf))
plot(d1, type="l", ylim=c(range(d1$y, d2$y)))
lines(d2, col=2)

range(plogis(x)*4-1)
d1 <- bounded_density(plogis(x)*4-1)
d2 <- bounded_density(plogis(x)*4-1, c(-1, 3))
plot(d1, type="l", ylim=c(range(d1$y, d2$y)))
lines(d2, col=2)

range(-exp(x)+3)
d1 <- bounded_density(-exp(x)+3)
d2 <- bounded_density(-exp(x)+3, c(-Inf, 3))
plot(d1, type="l", ylim=c(range(d1$y, d2$y)))
lines(d2, col=2)
par(op)
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
    for (i in 1:ncol(x)) {
        if (!is.null(xx[[i]])) {
            polygon(0.45*c(-xx[[i]]$w, rev(xx[[i]]$w))+i,
                c(xx[[i]]$h, rev(xx[[i]]$h)), col=col, border=border)
            #lines(c(i,i), range(xx[[i]]$h), col=col, lwd=2)
        }
        lines(c(i-0.2, i+0.2), rep(median(x[,i], na.rm=TRUE), 2), lwd=3)
    }
    #points(1:ncol(x), colMeans(x, na.rm=TRUE), pch=21, cex=1.5)
    invisible(NULL)
}
fu <- function(x) sign(x) * plogis(log(abs(x/100)))

se1 <- read.csv("c:/Users/Peter/sector-effects-birds-RefInRegion-inOSA.csv")
se2 <- read.csv("c:/Users/Peter/sector-effects-birds-RefUnderHF-inOSA.csv")
rownames(se1) <- se1[,1]
se1[,1] <- NULL
rownames(se2) <- se2[,1]
se2[,1] <- NULL

par(mfrow=c(1,2), las=1, yaxs="i")
p <- c(0,1)#p <- c(0.025, 0.975)
ylim <- c(-100, 200)
#ylim <- c(-50, 100)
vp(se1,
    main="% pop. change / ref. pop. in region (birds)", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(se2, #interval=c(-100,Inf),
    main="% pop. change / ref. pop. in HF (birds)", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)

