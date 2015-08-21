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
library(snow)
if (!interactive())
    library(Rmpi)
library(MASS)
library(ResourceSelection)
library(mefa4)
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/analysis_functions.R")
if (!interactive())
    (args <- commandArgs(trailingOnly = TRUE))

TEST <- interactive()

ROOT <- "~/abmi"

#### setup ####

nodes <- if (interactive())
    5 else as.numeric(args[1])
BBB <- if (TEST) 2 else 240 # = 4*5*12
ncl <- if (TEST) 2 else nodes*12

#### load all object on the master ####

if (interactive())
    setwd("c:/p/AB_data_v2015/out/birds")
fid <- if (interactive())
    "north" else as.character(args[2])
fn <- paste0("data-useok-", fid, ".Rdata")
load(file.path("data", fn))

if (TEST)
    mods <- mods[1:2]

#### spawning the slaves ####

cl <- if (interactive())
    makeCluster(ncl) else makeMPIcluster(ncl)
if (!interactive())
    options("CLUSTER_ACTIVE" = TRUE)

#### loading packages on slaves ####

tmpcl <- clusterEvalQ(cl, library(ResourceSelection))
tmpcl <- clusterEvalQ(cl, library(MASS))
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, source("~/repos/bragging/R/glm_skeleton.R"))
tmpcl <- clusterEvalQ(cl, source("~/repos/abmianalytics/R/analysis_functions.R"))

#### load all the objects on the slaves ####

tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("c:/p/AB_data_v2015/out/birds"))
tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))

#### project identifier ####

PROJECT <- if (TEST)
    paste0("abmi-", fid, "-test") else paste0("bam-", fid) 


#### checkpoint ####

SPP <- colnames(YY)
done_fl <- list.files("results")
done_fl <- done_fl[grepl(fid, done_fl)]
DONE <- substr(sapply(strsplit(done_fl, "_"), "[[", 3), 1, 4)
SPP <- setdiff(SPP, DONE)
SPP <- rev(SPP)
if (TEST)
    SPP <- SPP[1:2]

#system.time(aaa <- do_1spec1run_noW(1, i="ALFL", mods=mods,  hsh_name=NA, CAICalpha=1))
#system.time(aaa <- do_1spec1run_noW(1, i=SPP1, mods=mods,  hsh_name=NA, CAICalpha=1))
# j=2
# i="AMBI"
hsh_name <- NA
CAICalpha <- 1

## catch errors that cannot be dealt with internally in glm_skeleton
wg_fun <- function(...) try(do_1spec1run_noW(...))

for (SPP1 in SPP) {
    cat(SPP1, date(), "\n");flush.console()
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
    save(res, file=paste("results/birds_", PROJECT, "_", SPP1, ".Rdata", sep=""))
}

#### shutting down ####

stopCluster(cl)
if (!interactive()) {
    options("CLUSTER_ACTIVE" = FALSE)
    mpi.quit("no")
} else {
    quit("no")
}
