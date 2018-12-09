## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "ab-birds-validation-2018-12-07.RData"
## CAIC = alpha * AIC + (1 - alpha) * BIC, 1: AIC, 0: BIC
CAICalpha <- 1
## Number of bootstrap runs, 100 or 240
MaxB <- 100
## project name for storing the output
PROJ <- "validation"
## testing?

cat("loading packages etc...")
library(parallel)
library(mefa4)
library(opticut)
source("~/repos/abmianalytics/birds/00-functions.R")


## Create an array from the NODESLIST environnement variable
if (interactive()) {
    cat("OK\nNote: this is a test run...")
    nodeslist <- 2
    BBB <- 2
    setwd("e:/peter/AB_data_v2018/data/analytics/birds")
} else {
    cat("OK\ngetting nodes list...")
    nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
    BBB <- MaxB
    cat("OK\nNodes list:\n")
    print(nodeslist)
}

## Create the cluster with the nodes name.
## One process per count of node name.
## nodeslist = node1 node1 node2 node2, means we are starting 2 processes
## on node1, likewise on node2.
cat("spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")

cat("OK\nload data on master...")
load(file.path("data", fn))

cat("OK\nload packages on workers...")
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, library(opticut))
tmpcl <- clusterEvalQ(cl, source("~/repos/abmianalytics/birds/00-functions.R"))

cat("OK\nexporting and data loading on workers...")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("e:/peter/AB_data_v2018/data/analytics/birds"))
tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))

cat("OK\nsetting checkpoint...")
SPP <- colnames(YY)
#if (REVERSE)
#    SPP <- rev(SPP)
DONE <- substr(list.files(paste0("out/", PROJ)), 1, 4)
SPP <- setdiff(SPP, DONE)
if (interactive())
    SPP <- SPP[1:2]

cat("OK\nstart running models:")
for (SPP1 in SPP) {
    cat("\t", SPP1, date(), "...")
    if (interactive())
        flush.console()
    t0 <- proc.time()
    res <- parLapply(cl, 1:BBB, wg_fun, i=SPP1, mods=mods, CAICalpha=CAICalpha)
    attr(res, "timing") <- t0 - proc.time()
    attr(res, "proj") <- PROJ
    attr(res, "spp") <- SPP1
    attr(res, "CAICalpha") <- CAICalpha
    attr(res, "date") <- as.character(Sys.Date())
    attr(res, "ncl") <- length(cl)
    save(res,
        file=paste0("out/", PROJ, "/", SPP1, ".RData"))
    cat("OK\n")
}

## Releaseing resources.
cat("shutting down.\n\nDONE!\n")
stopCluster(cl)
q("no")
