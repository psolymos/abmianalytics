plot_control_chart1 <-
function(x, y, type=c("qcc","cusum","ewma"), n=8, nsigmas=4.5,
main, use.date=TRUE, offset=0.2, lambda=0.2, sort=TRUE)
{
    require(qcc)
    if (missing(main))
        main <- deparse(substitute(y))
    variable <- main
    type <- match.arg(type)
    xx <- data.frame(Date=x, Value=y)
    if (sort)
        xx <- xx[order(xx$Date),]
    var <- xx$Value
    if (use.date) {
        tt <- xx$Date
        tt <- tt[!is.na(var)]
    }
    var <- var[!is.na(var)]
    if (length(var) <= n) {
        plot.new()
        return(xx)
    }

    isBase <- logical(length(var))
    isBase[1:n] <- TRUE

    if (type == "qcc") {
        ylab <- variable
        fit <- qcc(data.matrix(var[1:n]),
            newdata=data.matrix(var[(n+1):length(var)]),
            type="xbar.one",
            nsigmas=nsigmas,
            std.dev="SD",
            #labels=dat$t[1:n], newlabels=dat$t[(n+1):length(var)],
            xlab="Date", ylab=ylab, title="Shewhart Control Chart",
            plot=FALSE)
        xy <- cbind(seq_len(length(var)), var)
        zzz <- c(fit$center, fit$limits)
        ylim <- range(c(xy[,2], zzz))

    }
    if (type == "cusum") {
        ylab <- paste("Cusum", variable)
        fit <- cusum(data.matrix(var[1:n]),
            newdata=data.matrix(var[(n+1):length(var)]),
            type="xbar.one",
            decision.interval=nsigmas,
            std.dev="SD",
            #labels=dat$t[1:n], newlabels=dat$t[(n+1):length(var)],
            xlab="Date", ylab=ylab, title="CUSUM Control Chart",
            plot=FALSE)
        xy <- cbind(seq_len(length(var)), fit$pos, fit$neg)
        zzz <- c(0, -fit$decision.interval, fit$decision.interval)
        ylim <- range(c(xy[,2:3], zzz))
    }
    if (type == "ewma") {
        ylab <- paste("EWMA", variable)
        fit <- ewma(data.matrix(var[1:n]),
            newdata=data.matrix(var[(n+1):length(var)]),
            type="xbar.one",
            nsigmas=nsigmas,
            std.dev="SD",
            #labels=dat$t[1:n], newlabels=dat$t[(n+1):length(var)],
            xlab="Date", ylab=ylab, title="EWMA Control Chart",
            plot=FALSE,
            lambda=lambda)
        xy <- cbind(seq_len(length(var)), fit$y, var)
        zzz <- cbind(fit$center, fit$limits)
        ylim <- range(c(xy[,2:3], zzz))
    }
    if (use.date && type != "ewma") {
        xy <- data.frame(xy)
        xy[,1] <- tt
    }

    xlim <- range(xy[,1])
    ylim <- c(ylim[1]-diff(ylim)*offset, ylim[2]+diff(ylim)*offset)

    xoff <- diff(xlim)/10
    yoff <- diff(ylim)/10
    End <- mean(xy[c(max(which(isBase)), max(which(isBase))+1), 1])
    col1 <- switch(type,
        "qcc"="grey",
        "cusum"="tan",
        "ewma"="wheat3")
    col2 <- switch(type,
        "qcc"="gold",
        "cusum"="pink",
        "ewma"="palegreen2")
    plot(xy[,1:2], type="n", las=1, #xaxs = "i", yaxs = "i",
        xlim=xlim, ylim=ylim, main=main,
        xlab="Date", ylab=ylab)
    if (type == "ewma") {
        ## nicer curve is needed here
        ## or just dates are screwing it up???
        polygon(c(xlim[1]-xoff, xy[isBase,1], End,
            End, rev(xy[isBase,1]), xlim[1]-xoff),
            c(zzz[1,2], zzz[isBase,2], zzz[sum(isBase),2],
             zzz[sum(isBase),3], rev(zzz[isBase,3]), zzz[1,3]),
            col=col1, border=NA)

        polygon(c(xlim[2]+xoff, rev(xy[!isBase,1]), End,
            End, xy[!isBase,1], xlim[2]+xoff),
            c(zzz[nrow(zzz),2], rev(zzz[!isBase, 2]),zzz[!isBase, 2][1],
            rep(ylim[1]-yoff, sum(!isBase)+2)),
            col=col2, border=NA)
        polygon(c(xlim[2]+xoff, rev(xy[!isBase,1]), End,
            End, xy[!isBase,1], xlim[2]+xoff),
            c(zzz[nrow(zzz),3], rev(zzz[!isBase, 3]),zzz[!isBase, 3][1],
            rep(ylim[2]+yoff, sum(!isBase)+2)),
            col=col2, border=NA)
        abline(h=zzz[1,1], lty=1, col=1)
    } else {
        polygon(c(xlim[1]-xoff, End, End, xlim[1]-xoff),
            zzz[c(2,2,3,3)],
            col=col1, border=NA)
        polygon(c(xlim[2]+xoff, End, End, xlim[2]+xoff),
            c(zzz[c(2,2)], ylim[c(1,1)]-yoff),
            col=col2, border=NA)
        polygon(c(xlim[2]+xoff, End, End, xlim[2]+xoff),
            c(zzz[c(3,3)], ylim[c(2,2)]+yoff),
            col=col2, border=NA)
        abline(h=zzz[1], lty=1, col=1)
    }
    box()

    if (type=="qcc") {
        lines(xy, col=4)
        points(xy[isBase,,drop=FALSE], pch=19, col=4, cex=1.2)
        points(xy[!isBase,,drop=FALSE], pch=21, col=4, cex=1.2,
            bg="white")
        if (length(fit$violations$beyond.limits))
            points(xy[fit$violations$beyond.limits,,drop=FALSE],
                pch=21, col=4, cex=1.2, bg="red")
        if (length(fit$violations$violating.runs))
            points(xy[fit$violations$violating.runs,,drop=FALSE],
                pch=21, col=4, cex=1.2, bg="orange")
    }
    if (type=="cusum") {
        lines(xy[,1:2], col=4)
        lines(xy[,c(1,3)], col=4)
        points(xy[isBase,1:2,drop=FALSE], pch=19, col=4, cex=1.2)
        points(xy[isBase,c(1,3),drop=FALSE], pch=19, col=4, cex=1.2)
        points(xy[!isBase,1:2,drop=FALSE], pch=21, col=4, cex=1.2,
            bg="white")
        points(xy[!isBase,c(1,3),drop=FALSE], pch=21, col=4, cex=1.2,
            bg="white")
        if (length(fit$violations$lower))
            points(xy[fit$violations$lower,,drop=FALSE],
                pch=21, col=4, cex=1.2, bg="red")
        if (length(fit$violations$upper))
            points(xy[fit$violations$upper,,drop=FALSE],
                pch=21, col=4, cex=1.2, bg="red")
    }
    if (type == "ewma") {
        lines(xy[,1:2,drop=FALSE], col=4)
        points(xy[isBase,1:2,drop=FALSE], pch=19, col=4, cex=1.2)
        points(xy[!isBase,1:2,drop=FALSE], pch=21, col=4, cex=1.2,
            bg="white")
        points(xy[,c(1,3),drop=FALSE], pch=3, col=4, cex=1.2)

        if (length(fit$violations))
            points(xy[fit$violations,1:2,drop=FALSE],
                pch=21, col=4, cex=1.2, bg="red")
    }
    invisible(xx)
}

