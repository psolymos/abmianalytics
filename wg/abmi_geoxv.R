library(MASS)
library(ResourceSelection)
library(mefa4)
source("~/repos/bragging/R/glm_skeleton.R")
#source("~/repos/abmianalytics/R/analysis_functions.R")

load("e:/peter/AB_data_v2016/out/birds/data/data-north-geoxv.Rdata")
OUTDIR <- "e:/peter/AB_data_v2016/out/birds/results/north-geoxv"

#load("e:/peter/AB_data_v2016/out/birds/data/data-south-geoxv.Rdata")
#OUTDIR <- "e:/peter/AB_data_v2016/out/birds/results/south-geoxv"

SPP <- colnames(YY)
SPP <- rev(SPP)

BBB <- max(BB[,2])

do_1spec1run_geoxv <- function(j, i, mods,
silent=FALSE, CAICalpha=1)
{
    jj <- if (j==0)
        BB[,1] else BB[BB[,2]!=j,1] # this is geoxv
    x <- DAT[jj,]
    y <- as.numeric(YY[jj, i])
    off <- if (i %in% colnames(OFF))
        OFF[jj, i] else OFFmean[jj]
    ## empty objects for storing results
    nmods <- length(mods)
    nnmods <- sapply(mods, length)
    mid <- numeric(nmods)
    bestList <- vector("list", nmods)
    caicList <- vector("list", nmods)
    ## Null
    null <- glm_skeleton(glm(y ~ 1,
        x,
        family=poisson(),
        offset=off,
        #weights=w,
        x=FALSE, y=FALSE, model=FALSE), CAICalpha=CAICalpha)
    best <- null
    ## looping through models list
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
        bestList[[l1]] <- best
    }
    ## final assembly
    out <- list(species=i, iteration=j,
        null=null$coef,
        null_caic=null$caic,
        caic=caicList,
        coef=lapply(bestList, "[[", "coef"),
        mid=mid,
        alpha=CAICalpha)
    out
}


## CAIC = alpha * AIC + (1 - alpha) * BIC
## 1: AIC, 0: BIC
CAICalpha <- 1

## catch errors that cannot be dealt with internally in glm_skeleton
wg_fun <- function(...) try(do_1spec1run_geoxv(...))

cat("OK\nstart running models:")
for (SPP1 in SPP) {
    cat("\t", SPP1, date(), "...")
    if (interactive())
        flush.console()
    res <- lapply(0:BBB, wg_fun, i=SPP1, mods=mods, CAICalpha=CAICalpha)
    save(res,
        file=paste0(OUTDIR, "/birds_abmi-north-geoxv_", SPP1, ".Rdata"))
    cat("OK\n")
}

