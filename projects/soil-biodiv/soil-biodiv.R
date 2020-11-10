## I. non-negative matrix factorization using alternating least squares

# some libraries
library(ALS) # for ALS
library(mefa4) # for some data manipulation

## this is taken and modified from package ALS
silent_als <-
function (CList, PsiList, S = matrix(), WList = list(), thresh = 0.001,
    maxiter = 100, forcemaxiter = FALSE, optS1st = TRUE, x = 1:nrow(CList[[1]]),
    x2 = 1:nrow(S), baseline = FALSE, fixed = vector("list",
        length(PsiList)), uniC = FALSE, uniS = FALSE, nonnegC = TRUE,
    nonnegS = TRUE, normS = 0, closureC = list(), silent=FALSE) {
    RD <- 10^20
    PsiAll <- do.call("rbind", PsiList)
    resid <- vector("list", length(PsiList))
    if (length(WList) == 0) {
        WList <- vector("list", length(PsiList))
        for (i in 1:length(PsiList)) WList[[i]] <- matrix(1,
            nrow(PsiList[[1]]), ncol(PsiList[[1]]))
    }
    W <- do.call("rbind", WList)
    for (i in 1:length(PsiList)) resid[[i]] <- matrix(0, nrow(PsiList[[i]]),
        ncol(PsiList[[i]]))
    for (j in 1:length(PsiList)) {
        for (i in 1:nrow(PsiList[[j]])) {
            resid[[j]][i, ] <- PsiList[[j]][i, ] - CList[[j]][i,
                ] %*% t(S * WList[[j]][i, ])
        }
    }
    initialrss <- oldrss <- sum(unlist(resid)^2)
    if (silent) cat("Initial RSS", initialrss, "\n")
    iter <- 1
    b <- if (optS1st)
        1
    else 0
    oneMore <- FALSE
    while (((RD > thresh || forcemaxiter) && maxiter >= iter) ||
        oneMore) {
        if (iter%%2 == b)
            S <- getS(CList, PsiAll, S, W, baseline, uniS, nonnegS,
                normS, x2)
        else CList <- getCList(S, PsiList, CList, WList, resid,
            x, baseline, fixed, uniC, nonnegC, closureC)
        for (j in 1:length(PsiList)) {
            for (i in 1:nrow(PsiList[[j]])) {
                resid[[j]][i, ] <- PsiList[[j]][i, ] - CList[[j]][i,
                  ] %*% t(S * WList[[j]][i, ])
            }
        }
        rss <- sum(unlist(resid)^2)
        RD <- ((oldrss - rss)/oldrss)
        oldrss <- rss
        typ <- if (iter%%2 == b)
            "S"
        else "C"
        if (silent) cat("Iteration (opt. ", typ, "): ", iter, ", RSS: ",
            rss, ", RD: ", RD, "\n", sep = "")
        iter <- iter + 1
        oneMore <- (normS > S && (iter%%2 != b) && maxiter !=
            1) || (length(closureC) > 0 && (iter%%2 == b))
    }
    if (silent) cat("Initial RSS / Final RSS =", initialrss, "/", rss, "=",
        initialrss/rss, "\n")
    return(list(CList = CList, S = S, rss = rss, resid = resid,
        iter = iter))
}

## standardizes the output in a more meaningful format
run_als <- function(x, k, weights=NULL, ...) {
    als_getG <- function(object) {
       out <- m$S
       Fsum <- rowSums(t(m$CList[[1]]))
       out <- t(t(out) * Fsum)
       out
    }
    als_getF <- function(object) {
       out <- t(m$CList[[1]])
       Fsum <- rowSums(out)
       out <- out / Fsum
       out[is.na(out)] <- 0
       out
    }
    if (is.null(weights)) {
        weights <- x
        weights[] <- 1
    }
    t0 <- proc.time()
    m <- silent_als(
        CList=list(matrix(runif(ncol(x)*k),ncol(x), k)), # = t(F)
        S=matrix(1,nrow(x),k), # = G
        PsiList=list(t(x)), # = t(X)
        WList=list(t(weights)), # = t(W)
        ...)
    out <- list(
        results=m,
        x = x,
        weights=weights,
        sources=als_getF(m),
        mixing=als_getG(m),
        xhat = t(m$CList[[1]] %*% t(m$S)),
        timing=proc.time()-t0)
    class(out) <- "run_als"
    out
}

# load data
load("~/Desktop/Otu raw data and LDA out.RData")

# take OTU data, k = # of factors
k <- 5
a <- run_als(as.matrix(OtuData), k=k)

# it is not straightforward how to best select number of k
# sum of squares
sum(a$weights * (a$x - a$xhat)^2)
plot(a$x[1,], a$xhat[1,]);abline(0,1)

# input x: is n*m matrix
# sources: this gives the 'fingerprint' of the communities, k*m matrix (F)
# mixing: this gives the mixing proportions of the sources at each site, n*k (G)
# xhat: Xhat=FG+error, k=number of sources

# Classifying OTUs into sources
SRC <- t(a$sources)
colSums(SRC) # source proportions sum to 1
dimnames(SRC) <- list(colnames(OtuData), paste0("Source", 1:k))
b <- find_max(SRC)
# gives the source with the highest value
head(b)
table(b$index)

MIX <- a$mixing
dimnames(MIX) <- list(rownames(OtuData), paste0("Source", 1:k))
# these are NOT proportions, if you want that: divide by row total
head(MIX)
head(MIX / rowSums(MIX))

# silhouette based on G matrix
# find the cluster membership based on highest mixing value
# use that and Euclidean distance among samples to calculate
# 'silhouette width' see ?cluster::silhouette for explanation
# https://en.wikipedia.org/wiki/Silhouette_(clustering)
g2sil <- function(G, D) {
    cl <- apply(G, 1, which.max)
    si <- cluster::silhouette(cl, D)
    summary(si)$avg.width
}
g2sil(a$mixing, dist(a$x))

# fit different k values
K <- 2:10
fit <- lapply(K, function(k) run_als(as.matrix(OtuData), k=k))

# look for highest
sil <- sapply(fit, function(z) g2sil(z$mixing, dist(z$x)))
plot(K, sil, type="b")
# this is not really useful: overfitting gives lower SSQ
ssq <- sapply(fit, function(z) sum(z$weights * (z$x - z$xhat)^2))
plot(K, ssq, type="b")

# find best ALS based on maximum silhouette
BEST <- which.max(sil)
K[BEST]
m <- fit[[BEST]]
cl <- apply(m$mixing, 1, which.max)
si <- cluster::silhouette(cl, dist(m$x))
plot(si)
summary(si)


## II. composition based on LDA output

library(VGAM)

load("~/Desktop/Otu raw data and LDA out.RData")

ymat <- as.matrix(Theta)

vm <- vglm(formula=ymat ~ TrtGp.otu, family=multinomial)

# multinomial expects counts, so let's cheat

ymat <- round(as.matrix(Theta) * rowSums(OtuData))
vm <- vglm(formula=ymat ~ TrtGp.otu, family=multinomial)
summary(vm)
## this should give the predictions
fitted(vm)

# multinomial probabilities for
f <- fitted(vm)[!duplicated(TrtGp.otu),]
rownames(f) <- TrtGp.otu[!duplicated(TrtGp.otu)]

