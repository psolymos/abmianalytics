## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "ab-birds-mixedwood-2019-10-31.RData"
## project name for storing the output
PROJ <- "mixedwood"

## CAIC = alpha * AIC + (1 - alpha) * BIC, 1: AIC, 0: BIC
CAICalpha <- 1
## Number of bootstrap runs, 100 or 240
MaxB <- 240

## test suite uses limited sets
TEST <- FALSE

if (TEST) {
    MaxB <- 2
    cat("* Note: this is a test run!\n")
}

cat("* Loading packages and sourcing functions:")
library(parallel)
library(mefa4)
library(opticut)
source("~/repos/abmianalytics/birds/00-functions.R")

## this version of the function subsets the data to 2006-2015
## so that we can compare to PIF v3 estimates
.run_path2 <- function(j, i, mods, CAICalpha=1, wcol=NULL,
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
    x <- DAT
    x$order <- seq_len(nrow(x))
    x <- x[BB[,j],]
    x <- x[x$YEAR >= 2006 & x$YEAR <= 2015,]
    y <- as.numeric(YY[x$order, i])
    off <- OFF[x$order, i]
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
                    formula=mods[[l1]][[l2]]), silent=TRUE), CAICalpha=CAICalpha)
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
            x$SSH_KM <- rowSums(SSH[x$order,lab,drop=FALSE])
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
run_path2 <- function(...) try(.run_path2(...))

## Create an array from the NODESLIST environnement variable
if (interactive()) {
    nodeslist <- 2
    BBB <- 2
    setwd("d:/abmi/AB_data_v2018/data/analysis/birds")
} else {
    cat("OK\n* Getting nodes list ... ")
    nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
    BBB <- MaxB
    cat("OK\n  Nodes list:\n")
    print(nodeslist)
}

## Create the cluster with the nodes name.
## One process per count of node name.
## nodeslist = node1 node1 node2 node2, means we are starting 2 processes
## on node1, likewise on node2.
cat("* Spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")

cat("OK\n* Loading data on master ... ")
load(file.path("data", fn))

cat("OK\nload packages on workers .. .")
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, library(opticut))
tmpcl <- clusterEvalQ(cl, source("~/repos/abmianalytics/birds/00-functions.R"))

cat("OK\n* Exporting and data loading on workers ... ")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("d:/abmi/AB_data_v2018/data/analysis/birds"))
#tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))
clusterExport(cl, c("DAT", "YY", "OFF", "BB", ".run_path2"))

cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(YY)
DONE <- character(0)
if (interactive() | TEST)
    SPP <- SPP[1:2]

DONE <- substr(list.files(paste0("out/", PROJ)), 1, 4)
TOGO <- setdiff(SPP, DONE)

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
while (length(TOGO) > 0) {
    SPP1 <- sample(TOGO, 1)
    cat("\n  -", length(DONE), "done,", length(TOGO), "more to go, doing", SPP1, "on", date(), "... ")
    if (interactive())
        flush.console()
    t0 <- proc.time()
    #z <- run_path2(1, "AMRO", mods, CAICalpha=1, wcol="vegw", ssh_class="vegc", ssh_fit="Space")
    if (interactive()) {
        res <- pblapply(cl=cl, X=1:BBB, FUN=run_path2,
            i=SPP1, mods=mods, CAICalpha=CAICalpha)
    } else {
        res <- parLapply(cl, 1:BBB, run_path2,
            i=SPP1, mods=mods, CAICalpha=CAICalpha)
    }
    attr(res, "timing") <- proc.time() - t0
    attr(res, "proj") <- PROJ
    attr(res, "spp") <- SPP1
    attr(res, "CAICalpha") <- CAICalpha
    attr(res, "date") <- as.character(Sys.Date())
    attr(res, "ncl") <- length(cl)
    save(res,
        file=paste0("out/", PROJ, "/", SPP1, ".RData"))
    DONE <- substr(list.files(paste0("out/", PROJ)), 1, 4)
    TOGO <- setdiff(SPP, DONE)
    cat("OK")
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
