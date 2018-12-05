#' ---
#' title: "Functions"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Nov 28, 2018"
#' output: pdf_document
#' ---
#'
#' SQL tables represent character not factor
make_char2fact <- function(x) {
    if (is.null(dim(x)))
        if (is.character(x))
            return(as.factor(x))
    for (i in seq_len(ncol(x)))
        if (is.character(x[,i]))
            x[,i] <- as.factor(x[,i])
        x
}
#' This turns (`1`, `11`) into (` 1`, `11`) [same char width]
fill_char <- function(x, fill=" ") {
    x <- as.character(x)
    n <- nchar(x)
    m <- max(n)
    sapply(seq_len(length(x)), function(i) {
        paste0(c(rep(fill, m-n[i]), x[i]), collapse="")
    })
}
#' cleans up AOU codes
normalize_sppcode <- function(x) {
    f <- function(x) {
        x <- toupper(x)
        x <- gsub(" ", "", x)
        x
    }
    if (is.factor(x)) {
        levels(x) <- f(levels(x))
        return(x)
    }
    if (is.character(x)) {
        return(f(x))
    }
    stop("x must be character or factor")
}
#' repeated calls to gsub
Gsub <- function(pattern, replacement, x, ...) {
    pattern <- as.list(pattern)
    replacement <- as.list(replacement)
    if (length(pattern) > length(replacement)) {
        j <- rep(seq_along(replacement), length(pattern))[seq_along(pattern)]
        replacement <- replacement[j]
    }
    if (length(pattern) < length(replacement)) {
        j <- rep(seq_along(pattern), length(replacement))[seq_along(replacement)]
        pattern <- pattern[j]
    }
    for (i in seq_along(pattern)) {
        x <- gsub(pattern[[i]], replacement[[i]], x, ...)
    }
    x
}
#' repeated calls to grepl, evaluating if ANY or all TRUE (NA treated as FALSE)
Grepl <- function(pattern, x, ..., ANY=TRUE) {
    if (length(pattern) == 1L)
        return(grepl(pattern, x, ...))
    z <- sapply(pattern, grepl, x=x, ...)
    z[is.na(z)] <- FALSE
    rs <- rowSums(z)
    if (ANY)
        rs > 0 else rs == ncol(z)
}
#' row standardization, avoids division by 0
row_std <- function(x) {
    rs <- rowSums(x)
    rs[rs == 0] <- 1
    #rs[is.na(rs)] <- 1
    x / rs
}
#' transforming space/climate variables
transform_clim <- function(x) {
    out <- data.frame(
        xX   = (x$X - (-113.3)) / 2.12,
        xY   = (x$Y - 54.1) / 2.43,
        xAHM = (x$AHM - 0) / 50,
        xPET = (x$PET - 0) / 800,
        xFFP = (x$FFP - 0) / 130,
        xMAP = (x$MAP - 0) / 2300,
        xMAT = (x$MAT - 0) / 6,
        xMCMT= (x$MCMT - 0) / 25,
        xMWMT= (x$MWMT - 0) / 20)
        #xASP= x$ASP
        #xSLP= log(x$SLP + 1)
        #xTRI= log(x$TRI / 5)
        #xCTI= log((x$CTI + 1) / 10)
        rownames(out) <- rownames(x)
        out$xX2 <- out$xX^2
        out$xY2 <- out$xY^2
        out
}
#' get terms from a list of formulas
get_terms <- function(mods, type=c("formula", "list"), intercept=TRUE) {
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
#' make formula term names sorted and predictable, i.e. always A:B instead of B:A
fix_names <- function(x, sep=":") {
    unlist(lapply(x, function(z) {
        paste(sort(strsplit(z, sep)[[1]]), collapse=sep)
    }))
}
#' logLik method for a failed model: make it the worst possible
'logLik.try-error' <- function (object, ...) {
    structure(-.Machine$double.xmax^(1/3), df = 1,
        nobs = 1, class = "logLik")
}
#' leave the baggage
glm_skeleton <- function(object, ..., CAICalpha=0.5) {
    ## strip parent env ref
    if (inherits(object, "try-error"))
        return(structure(as.character(object), class="try-error"))
    ## this is considered the skeleton
    out <- structure(list(
        call=object$call,
        formula=formula(object),
        coef=coef(object),
        converge=object$converge,
        logLik=as.numeric(logLik(object)),
        df=attr(logLik(object), "df"),
        nobs=nobs(object)), class="glm_skeleton")
    ## few more things that might come handy
    out$class0 <- class(object)[1L]
    out$aic <- -2*out$logLik + 2*out$df
    out$bic <- -2*out$logLik + log(out$nobs)*out$df
    out$caic <- CAICalpha * out$aic + (1-CAICalpha) * out$bic
    out
}
#' most useful stats
fstat <- function(x, level=0.95, ...) {
    .fstat <- function(x, level=0.95)
        c(Mean=mean(x, ...), Median=median(x, ...),
            quantile(x, c((1-level)/2, 1 - (1-level)/2), ...))
    if (is.null(dim(x)))
        .fstat(x, level, ...) else t(apply(x, 1, .fstat, level=level, ...))
}
#' workhorse function for branching model selection
#' - `j`: bootstrap run
#' - `i`: species
#' - `mods`: model formula list
#' - `CAICalpha`: defailts to AIC
#' - `wcol`: column ID in `DAT` for weights (weights=1 if NULL)
#' - `ssh_class`: NULL or column name with land cover classes
#'   (matching column names of SSH)
#' - `ssh_fit`: model stage to be used for getting fitted values for Lc cut
#' NULL or the stage of the model to be considered for thresholding for
#'   finding 1 km$^2$ level surrounding suitable habitat.
#'   If not NULL, it expectes a stage called SSH (and predictors SSH_KM and SSH05_KM).
#' need the following objects in parent env:
#' - `DAT`: data frame with predictors
#' - `OFF`: offsets matrix
#' - `BB`: bootstrap sample ID matrix
#' - `YY`: species cross tab
#' - `SSH`: optional matrix with surrounding compositional data
run_path1 <- function(j, i, mods, CAICalpha=1, wcol=NULL,
    ssh_class=NULL, ssh_fit=NULL)
{
    t0 <- proc.time()
    if (!is.null(ssh_class)) {
        if (!(ssh_class %in% colnames(DAT)))
            stop("ssh_class not found in DAT")
        if (is.null(ssh_fit))
            stop("specify model stage ssh_fit")
        if (!(ssh_fit %in% names(mods)))
            stop("ssh_fit not found in model stages")
        if (!("SSH" %in% names(mods)))
            stop("SSH not found in model stages")
        if (which(names(mods) == ssh_fit) >= which(names(mods) == "SSH"))
            stop("fitting stage must precede SSH stage in model list")
    }
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    w <- if (is.null(wcol))
        rep(1, nrow(x)) else x[,wcol]
    nmods <- length(mods)
    nnmods <- sapply(mods, length)
    mid <- numeric(nmods)
    bestList <- vector("list", nmods)
    caicList <- vector("list", nmods)
    names(mid) <- names(bestList) <- names(caicList) <- names(mods)
    null <- glm_skeleton(
        glm(
            formula=y ~ 1,
            data=x,
            family=poisson(),
            offset=off,
            weights=x$w,
            x=FALSE,
            y=FALSE,
            model=FALSE),
        CAICalpha=CAICalpha)
    best <- null
    res <- NULL
    for (l1 in 1:nmods) {
        if (nnmods[l1] > 0) {
            mlist <- vector("list", nnmods[l1])
            glist <- vector("list", nnmods[l1])
            for (l2 in 1:nnmods[l1]) {
                mlist[[l2]] <- glm_skeleton(try(update(object=best,
                    formula=mods[[l1]][[l2]]), silent=silent), CAICalpha=CAICalpha)
            }
            mcaic <- sapply(mlist, "[[", "caic")
            attr(mcaic, "StartCAIC") <- best$caic
            for (l2 in 1:length(mlist)) { # check convergence
                if (mlist[[l2]]$class != "try-error" && !mlist[[l2]]$converge)
                    mcaic[l2] <- 2*.Machine$double.xmax^(1/3)
            }
            dcaic <- mcaic - best$caic
            mmid <- which.min(dcaic)
            if (dcaic[mmid] < 0) {
                best <- mlist[[mmid]]
                mid[l1] <- mmid
                gofbest <- glist[[mmid]]
            }
            caicList[[l1]] <- mcaic
        }
        if (!is.null(ssh_class) && names(mods)[l1] == ssh_fit) {
            est <- best$coef
            Xn <- model.matrix(get_terms(mods, "formula"), x)
            colnames(Xn) <- fix_names(colnames(Xn))
            names(est) <- fix_names(names(est))
            Xn <- Xn[,names(est),drop=FALSE]
            pr <- exp(drop(Xn %*% est))

            g <- sum_by(pr, x[,ssh_class])
            l <- lorenz(g[,"x"]/g[,"by"], g[,"by"])
            s <- summary(l)
            lab <- rownames(l[l[,"x"] >= s["x[t]"],])
            x$SSH_KM <- rowSums(SSH[BB[,j],lab,drop=FALSE])
            x$SSH05_KM <- sqrt(x$SSH_KM)

            res <- list(lorenz=l, labels=lab)
        }
        bestList[[l1]] <- best
    }
    out <- list(species=i, iteration=j,
        null=null$coef,
        null_caic=null$caic,
        caic=caicList,
        coef=lapply(bestList, "[[", "coef"),
        mid=mid,
        alpha=CAICalpha,
        wcol=wcol,
        ssh=res,
        timer=proc.time()-t0)
    out
}

