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
    "south" else as.character(args[2])
fn <- paste0("data-useok-", fid, ".Rdata")
load(file.path("data", fn))

if (TEST)
    mods <- mods[1:3]

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
DONE <- substr(sapply(strsplit(list.files("results"), "_"), "[[", 3), 1, 4)
SPP <- setdiff(SPP, DONE)
if (TEST)
    SPP <- SPP[1:2]

#system.time(aaa <- do_1spec1run_noW(1, i=SPP[1], mods=mods,  hsh_name=NA, CAICalpha=1))
# j=1
# i="ALFL"
hsh_name <- NA
CAICalpha <- 1

for (SPP1 in SPP) {
    cat(SPP1, date(), "\n");flush.console()
    res <- parLapply(cl, 1:BBB, do_1spec1run_noW, i=SPP1, mods=mods, 
        hsh_name=hsh_name, CAICalpha=CAICalpha)
    #res <- wg_fun2(1, i=SPP[1], mods=mods, 
    #    output="return", path="results", project=PROJECT, 
    #    use_wt=TRUE, ip_name=ip_name, CAICalpha=CAICalpha, nmax=nmax)
    attr(res, "fid") <- fid
    attr(res, "spp") <- SPP1
    attr(res, "hsh_name") <- hsh_name
    attr(res, "CAICalpha") <- CAICalpha
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
