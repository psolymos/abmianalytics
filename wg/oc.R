#### preliminaries ####

if (!interactive()) {
.Last <- function() {
    if (getOption("CLUSTER_ACTIVE")) {
        stopCluster(cl)
        cat("active cluster stopped by .Last\n")
    } else {
        cat("no active cluster found\n")
    }
}
options("CLUSTER_ACTIVE" = FALSE)
}
cat("loading packages...")
library(snow)
if (!interactive())
    library(Rmpi)
library(MASS)
library(ResourceSelection)
library(mefa4)
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/analysis_functions.R")

cat("OK\ngetting args and setup...")
if (!interactive())
    (args <- commandArgs(trailingOnly = TRUE))

TEST <- interactive()

#### setup ####

nodes <- if (interactive())
    5 else as.numeric(args[1])
BBB <- if (TEST) 2 else 240 # = 4*5*12
ncl <- if (TEST) 2 else nodes*12

#### load all object on the master ####

if (interactive())
    setwd("e:/peter/AB_data_v2016/out/birds")

cat("OK\nload data on master...")
fn <- paste0("data-north.Rdata")
load(file.path("data", fn))

cut_1spec1run_noW <- function(j, i, z)
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    HABV <- x[,z]
    ## opticut based approach for core habitat delineation
    require(opticut)
    oc <- opticut(y ~ ROAD01, data=x, strata=HABV, dist="poisson",
        offset=off, comb="rank")
    part <- drop(bestpart(oc))
    habmod_oc <- glm_skeleton(bestmodel(oc)[[1]])
    ocres <- drop(as.matrix(summary(oc)$summary[,c("I","beta0","beta1","logLR","w")]))
    Prob <- table(HABV, part)[,"1"]
    ## missing/dropped levels are NaN=0/0
    Prob[is.na(Prob)] <- 0
    Hi_oc <- names(Prob)[Prob > 0]
    ## Lorenz-tangent approach for core habitat delineation
    habmod <- glm_skeleton(try(glm(y ~ HABV + ROAD01,
        x,
        family=poisson(), 
        offset=off, 
        #weights=w,
        x=FALSE, y=FALSE, model=FALSE)))
    ## need to correct for linear effects
    ## so that we estimate potential pop in habitats (and not realized)
    XHSH <- model.matrix(~ HABV + ROAD01, x)
    XHSH[,"ROAD01"] <- 0 # not predicting edge effects
    ## some levels might be dropped (e.g. Marsh)
    XHSH <- XHSH[,names(habmod$coef)]
    lam <- exp(drop(XHSH %*% habmod$coef))
    cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
    Freq <- table(hab=HABV, lc=ifelse(lam >= cv$lam, 1, 0))
    Prob <- Freq[,"1"] / rowSums(Freq)
    ## missing/dropped levels are NaN=0/0
    Prob[is.na(Prob)] <- 0
    Hi_lc <- names(Prob)[Prob > 0.5]
    ## final assembly
    out <- list(species=i, iteration=j,
        hi_oc=Hi_oc,
        hi_lc=Hi_lc,
        lc=cv,
        ocres=ocres,
        #nmax=nmax,
        #w_id=w_id,
        habmod_oc=habmod_oc$coef,
        habmod_lc=habmod$coef)
    out
}

#### spawning the slaves ####
cat("OK\nspawning slaves...")
cl <- if (interactive())
    makeCluster(ncl) else makeMPIcluster(ncl)
if (!interactive())
    options("CLUSTER_ACTIVE" = TRUE)

#### loading packages on slaves ####
cat("OK\nload packages on slaves...")
tmpcl <- clusterEvalQ(cl, library(ResourceSelection))
tmpcl <- clusterEvalQ(cl, library(MASS))
tmpcl <- clusterEvalQ(cl, library(mefa4))

library(opticut)
tmpcl <- clusterEvalQ(cl, library(opticut))

tmpcl <- clusterEvalQ(cl, source("~/repos/bragging/R/glm_skeleton.R"))
#tmpcl <- clusterEvalQ(cl, source("~/repos/abmianalytics/R/analysis_functions.R"))

#### load all the objects on the slaves ####

cat("OK\nexporting and data loading on slaves...")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("e:/peter/AB_data_v2016/out/birds"))
tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))

#### project identifier ####

PROJECT <- if (TEST)
    paste0("abmi-oc-test") else paste0("abmi-oc")

#### checkpoint ####
cat("OK\nsetting checkpoint...")
SPP <- colnames(YY)
done_fl <- list.files(paste0("results/oc"))
DONE <- substr(sapply(strsplit(done_fl, "_"), "[[", 3), 1, 4)
SPP <- setdiff(SPP, DONE)
if (TEST)
    SPP <- SPP[1:2]

cat("OK\nstart running models:")
for (SPP1 in SPP) {
    cat("\t", SPP1, date(), "...")
    if (interactive())
        flush.console()
    res <- parLapply(cl, 1:BBB, cut_1spec1run_noW, i=SPP1, z="hab1ec")
    save(res, file=paste0("results/oc/birds_", PROJECT, "_", SPP1, ".Rdata"))
    cat("OK\n")
}

#### shutting down ####
cat("shutting down\n")
stopCluster(cl)
if (!interactive()) {
    options("CLUSTER_ACTIVE" = FALSE)
    mpi.quit("no")
} else {
    quit("no")
}
