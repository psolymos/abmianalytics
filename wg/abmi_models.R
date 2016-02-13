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
fid <- if (interactive())
    "north" else as.character(args[2])

subset <- if (interactive())
    "full" else as.character(args[3]) # full or useok

cat("OK\nload data on master...")
fn <- paste0("data-", subset, "-", fid, ".Rdata")
load(file.path("data", fn))

CAICalpha <- 1

USE_HSH <- if (interactive())
    TRUE else as.character(args[4]) == "do_hsh"
if (fid == "south")
    USE_HSH <- FALSE

if (USE_HSH) {
    hsh_name <- "hab1ec"
    mods <- mods2
    hshid <- "dohsh"
} else {
    hsh_name <- NA
    hshid <- "nohsh"
}

if (TEST)
    mods <- mods[1:2]

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
if (hshid == "dohsh") {
    library(opticut)
    tmpcl <- clusterEvalQ(cl, library(opticut))
}
tmpcl <- clusterEvalQ(cl, source("~/repos/bragging/R/glm_skeleton.R"))
tmpcl <- clusterEvalQ(cl, source("~/repos/abmianalytics/R/analysis_functions.R"))

#### load all the objects on the slaves ####

cat("OK\nexporting and data loading on slaves...")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("c:/p/AB_data_v2015/out/birds"))
tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))

#### project identifier ####

PROJECT <- if (TEST)
    paste0("abmi-", hshid, "-", fid, "-test") else paste0("abmi-", hshid, "-", fid)


#### checkpoint ####
cat("OK\nsetting checkpoint...")
SPP <- colnames(YY)
done_fl <- list.files("results")
done_fl <- done_fl[grepl(fid, done_fl)]
DONE <- substr(sapply(strsplit(done_fl, "_"), "[[", 3), 1, 4)
SPP <- setdiff(SPP, DONE)
if (TEST)
    SPP <- SPP[1:2]

## restrict to 2 species
#SPP <- c("WEWP","RWBL")
SPP <- c("OSFL","OVEN","CAWA","ALFL","BOCH")

#system.time(aaa <- do_1spec1run_noW(1, i="WEWP", mods=mods2, hsh_name="hab1ec", CAICalpha=1, method="oc"))
#system.time(aaa <- do_1spec1run_noW(1, i="RWBL", mods=mods,  hsh_name=NA, CAICalpha=1))
#system.time(aaa <- do_1spec1run_noW(1, i="WEWP", mods=mods,  hsh_name=NA, CAICalpha=1))
#system.time(aaa <- do_1spec1run_noW(1, i="ALFL", mods=mods,  hsh_name=NA, CAICalpha=1))
#system.time(aaa <- do_1spec1run_noW(1, i=SPP1, mods=mods,  hsh_name=NA, CAICalpha=1))
# j=2
# i="AMBI"

## catch errors that cannot be dealt with internally in glm_skeleton
wg_fun <- function(...) try(do_1spec1run_noW(...))

cat("OK\nstart running models:")
for (SPP1 in SPP) {
    cat("\t", SPP1, date(), "...")
    if (interactive())
        flush.console()
    t0 <- proc.time()
    res <- parLapply(cl, 1:BBB, wg_fun, i=SPP1, mods=mods,
        hsh_name=hsh_name, CAICalpha=CAICalpha)
#    res <- lapply(1:BBB, wg_fun, i=SPP1, mods=mods,
#        hsh_name=hsh_name, CAICalpha=CAICalpha)
    attr(res, "timing") <- t0 - proc.time()
    attr(res, "fid") <- fid
    attr(res, "spp") <- SPP1
    attr(res, "hsh_name") <- hsh_name
    attr(res, "CAICalpha") <- CAICalpha
    attr(res, "date") <- as.character(Sys.Date())
    attr(res, "ncl") <- ncl
    save(res, file=paste0("results/birds_", PROJECT, "_", SPP1, ".Rdata"))
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
