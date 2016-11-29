## old-forest birds guild results

library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"

level <- 0.9

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
#    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

slt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(slt) <- slt$AOU
slt$comments <- NULL

## guild and members
gname <- "old-forest-birds"
NAM <- "Old Forest Birds"
gspp <- rownames(slt)[slt$oldforest == 1 & slt$map.pred]

## coefs

load(file.path(ROOT, "tables", "res_coef.Rdata"))

tp <- combine_spp_coefs(res_coef, gspp)
if (max(tp$max, na.rm=TRUE) > 3*min(tp$max, na.rm=TRUE)) {
    MAXn <- tp$max[1]
    MAXs <- tp$max[2]
} else {
    MAXn <- max(tp$max, na.rm=TRUE)
    MAXs <- max(tp$max, na.rm=TRUE)
}
if (is.na(MAXn))
    MAXn <- max(tp$max, na.rm=TRUE)
if (is.na(MAXs))
    MAXs <- max(tp$max, na.rm=TRUE)

prn <- tp$veg
NDAT <- length(tp$species$north)
## veghf
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-veghf-north.png"))
png(file=fname,width=1500,height=700)
fig_veghf(prn,
    paste0(NAM, " (n = ", NDAT, " species)"), ymax=MAXn)
dev.off()
## linear
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-linear-north.png"))
png(file=fname,width=350,height=400)
fig_linear(attr(prn, "linear"),
    paste0(NAM, "\nNorth (n = ", NDAT, " species)"))
dev.off()

prs <- tp$soil
NDAT <- length(tp$species$south)
## treed
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-soilhf-treed-south.png"))
png(file=fname,width=500,height=450)
fig_soilhf(prs$treed,
    paste0(NAM, ", South, Treed (n = ", NDAT, " species)"),
    ymax=MAXs)
dev.off()
## nontreed
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-soilhf-nontreed-south.png"))
png(file=fname,width=500,height=450)
fig_soilhf(prs$nontreed,
    paste0(NAM, ", South, Non-treed (n = ", NDAT, " species)"),
    ymax=MAXs)
dev.off()
## linear
fname <- file.path(ROOT, "guilds", gname,
    paste0(gname, "-linear-south.png"))
png(file=fname,width=350,height=400)
fig_linear(prs$linear,
    paste0(NAM, "\nSouth (n = ", NDAT, " species)"))
dev.off()

## sector effects

load(file.path(ROOT, "tables", "sector-effects.Rdata"))

gseff <- combine_spp_seff(seff_res, gspp)



    png(file.path(ROOT, "out", "birds", "figs", paste0("sector-", WHERE),
        paste0(as.character(tax[spp, "Spp"]), TAG, ".png")),
        width=600, height=600)
    dev.off()

