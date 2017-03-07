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

TEST <- FALSE#interactive()

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

f_opticut <- function(j, i, z, dist="poisson")
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    if (dist == "binomial") {
        y <- ifelse(y > 0, 1, 0)
        dist <- "binomial:cloglog"
    }
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    HABV <- x[,z]
    ## opticut based approach for core habitat delineation
    require(opticut)
    oc <- opticut(y ~ ROAD01, data=x, strata=HABV, dist=dist,
        offset=off, comb="rank")
    part <- drop(bestpart(oc))
    habmod_oc <- glm_skeleton(bestmodel(oc)[[1]])
    ocres <- drop(as.matrix(summary(oc)$summary[,c("I","beta0","beta1","logLR","w")]))
    Prob <- table(HABV, part)[,"1"]
    ## missing/dropped levels are NaN=0/0
    Prob[is.na(Prob)] <- 0
    Hi <- names(Prob)[Prob > 0]
    ## final assembly
    out <- list(species=i, iteration=j,
        oc1=oc$species[[1]],
        dist=dist,
        hi=Hi,
        ocres=ocres,
        habmod=habmod_oc$coef)
    out
}
f_lorenz <- function(j, i, z)
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    HABV <- x[,z]
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

    ## density based
    lam <- exp(drop(XHSH %*% habmod$coef))
    cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
    Freq <- table(hab=HABV, lc=ifelse(lam >= cv$lam, 1, 0))
    Prob <- Freq[,"1"] / rowSums(Freq)
    Prob[is.na(Prob)] <- 0
    Hi <- names(Prob)[Prob > 0.5]

    ## probability based
    p <- 1-exp(-lam)
    cv2 <- Lc_cut(lam, transform=FALSE) # $lam is threshold
    Freq2 <- table(hab=HABV, lc=ifelse(p >= cv2$lam, 1, 0))
    Prob2 <- Freq2[,"1"] / rowSums(Freq2)
    Prob2[is.na(Prob2)] <- 0
    Hi2 <- names(Prob2)[Prob2 > 0.5]

    ## final assembly
    out <- list(species=i, iteration=j, habmod=habmod$coef,
        lam=list(
            hi=Hi,
            freq=Freq,
            lc=cv),
        p=list(
            hi=Hi2,
            freq=Freq2,
            lc=cv2))
    out
}
f_optilevels <- function(j, i, z) {
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    require(opticut)
    Z <- model.matrix(~ ROAD01 - 1, x)
    Z[,1] <- Z[,1] - mean(Z[,1], na.rm=TRUE)
    optilevels(y=y, x=x[,z], z=Z, dist="poisson",
        offset=off)
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
SPP <- colnames(OFF)
done_fl <- list.files(paste0("results/oc"))
DONE <- substr(sapply(strsplit(done_fl, "_"), "[[", 3), 1, 4)
SPP <- setdiff(SPP, DONE)
if (TEST)
    SPP <- SPP[1:2]

if (FALSE) {
SPP1 <- "OVEN"
system.time(m1 <- f_opticut(1, i=SPP1, z="hab1ec", dist="binomial")) # 49 sec
gc()
system.time(m2 <- f_opticut(1, i=SPP1, z="hab1ec", dist="poisson")) # 42 sec
gc()
system.time(m3 <- f_lorenz(1, i=SPP1, z="hab1ec")) # 5 sec
gc()
system.time(m4 <- f_optilevels(1, i=SPP1, z="hab1ec")) #  sec
gc()

library(opticut)
SPP <- sort(colnames(OFF))
all_res <- list()
for (i in 1:length(SPP)) {
    SPP1 <- SPP[i]
    cat("\n", SPP1, date(), "(", i, "/", length(SPP), ") --- bin");flush.console()
    m1 <- f_opticut(1, i=SPP1, z="hab1ec", dist="binomial")
    cat("\tpois");flush.console()
    m2 <- f_opticut(1, i=SPP1, z="hab1ec", dist="poisson")
    cat("\tlc");flush.console()
    m3 <- f_lorenz(1, i=SPP1, z="hab1ec")
    res <- list(opticut_binomial=m1, opticut_poisson=m2, lorenz=m3)
    all_res[[SPP1]] <- res
}

save(all_res, file="~/Dropbox/josm/2016/oc/oc-local-all-spp-2016-11-04.Rdata")
}

cat("OK\nstart running models:")
for (SPP1 in SPP) {
    cat("\t", SPP1, date(), "...")
    if (interactive())
        flush.console()
    m1 <- parLapply(cl, 1:BBB, f_opticut, i=SPP1, z="hab1ec", dist="binomial")
    m2 <- parLapply(cl, 1:BBB, f_opticut, i=SPP1, z="hab1ec", dist="poisson")
    m3 <- parLapply(cl, 1:BBB, f_lorenz, i=SPP1, z="hab1ec")
    #m4 <- parLapply(cl, 1:BBB, f_optilevels, i=SPP1, z="hab1ec")
    res <- list(opticut_binomial=m1, opticut_poisson=m2,
        lorenz=m3)

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
