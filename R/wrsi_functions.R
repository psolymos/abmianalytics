## wrsi: weighted relative suitability index

wrsi <-
function(y, x)
{
    Y <- as.integer(ifelse(y > 0, 1L, 0L))
    X <- data.matrix(x)
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
    #Var <- (1/colSums(Xu)) - (1/sum(Xu)) + (1/colSums(X)) - (1/sum(X))
    res <- data.frame(
        WRSI=WRSI,
        zWRSI=log(WRSI),
        rWRSI=(exp(2 * log(WRSI)) - 1)/(1 + exp(2 * log(WRSI))),
        Pused=Pused,
        Pavail=Pavail,
        Pw=Pw,
        u=colSums(Xu),
        a=colSums(X))
    rownames(res) <- colnames(X)
    class(res) <- c("wrsi", "data.frame")
    res
}

plot_wrsi <- 
function(Y, X, south=FALSE, ...)
{
    ## calculate the index to plot
    w <- wrsi(Y, X)
    rw <- w$rWRSI
    names(rw) <- rownames(w)
    ## set up a color palette
    if (south)
        col <- c(rep("brown",4), rep("grey30",2))
    if (!south)
        col <- brewer.pal(8, "Accent")[c(1,1,1,1, 2,2, 3,5, 7, 8,8, 8,8)]
    ## adjust margins and text direction
    op <- par(mar=c(8,4,2,2)+0.1, las=2)
    ## produce a barplot without axes
    tmp <- barplot(rw, horiz=FALSE, ylab="Affinity", 
        space=NULL, col=col, border=col, #width=sqrt(w$Pavail),
        ylim=c(-1,1), axes=FALSE, ...)
    ## add in the axis
    axis(2)
    ## a line
    abline(h=0, col="red4", lwd=2)
    ## reset par options
    par(op)
    ## return the index values invisibly
    invisible(w)
}

