#### define inputs ----------------------------
## root folder for the files
#ROOT <- "w:/reports"
ROOT <- "http://ftp.public.abmi.ca/species.abmi.ca/reports"
## version
VER <- "2017"
## if bootstrap based confidence intervals are needed
CI <- TRUE
## confidence interval coverage
LEVEL <- 0.9
## species to be considered
## "all": use all species
## any of "birds" "lichens" "mammals" "mites" "mosses" "vplants":
##   all species in a given taxon
## otherwise the intersect of provided and available species IDs are used
#SPP <- "all"
#SPP <- "birds"
SPP <- "Ovenbird"
## pixel IDs for spatial subset
PIX <- c("182_362", "182_363", "182_364", "182_365", "182_366", "182_367",
    "182_368", "182_369", "182_370", "182_371", "182_372")
#### load common objects ----------------------------------
## KA_2012, KA_2014: sector areas by 1km unit
## KT: 10km unit mapping to 1km units
## XY: coordinates of 1km units
## SP: species lookup table
library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(raster)
LOCAL <- !startsWith(ROOT, "http://")
if (LOCAL) {
    load(file.path(ROOT, VER, "data", "kgrid_areas_by_sector.RData"))
} else {
    con <- url(file.path(ROOT, VER, "data", "kgrid_areas_by_sector.RData"))
    load(con)
    close(con)
}

#### define species and pixel lists ----------------------------------
PIX <- sort(intersect(PIX, rownames(KT)))
PIX10 <- sort(unique(as.character(KT[PIX, "Row10_Col10"])))
KTsub <- KT[PIX,,drop=FALSE]
A_2012 <- colSums(KA_2012[PIX,])
A_2014 <- colSums(KA_2014[PIX,])
if (length(SPP) == 1L && SPP %in% c("all","birds","lichens","mammals","mites","mosses","vplants")) {
    SPPfull <- if (SPP == "all")
        as.character(SP$SpeciesID) else as.character(SP[SP$taxon==SPP, "SpeciesID"])
} else {
    SPPfull <- intersect(SPP, as.character(SP$SpeciesID))
}
SPPlist <- list(
    "birds"=character(0),
    "lichens"=character(0),
    "mammals"=character(0),
    "mites"=character(0),
    "mosses"=character(0),
    "vplants"=character(0))
for (i in names(SPPlist)) {
    SPPlist[[i]] <- intersect(SPPfull, as.character(SP[SP$taxon==i, "SpeciesID"]))
}
a <- c(0.5*(1-LEVEL), 1-0.5*(1-LEVEL))
cn <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")
#### processing the data
#taxon <- "birds"
#spp <- "Ovenbird"
RESULTS <- list()
## loop for taxa
for (taxon in names(SPPlist)) {
## loop for species within taxon
RESULTS[[taxon]] <- list()
for (spp in SPPlist[[taxon]]) {

## need to check if file exists
t0 <- proc.time()
if (LOCAL) {
    load(file.path(ROOT, VER, "results", taxon, "sector", paste0(spp, ".RData")))
} else {
    con <- url(file.path(ROOT, VER, "results", taxon, "sector", paste0(spp, ".RData")))
    load(con)
    close(con)
}

MAX <- max(max(rowSums(SA.Curr)), max(rowSums(SA.Ref)))
SA.Curr <- SA.Curr[PIX,cn]
SA.Ref <- SA.Ref[PIX,cn]
MEAN <- max(mean(rowSums(SA.Curr)), mean(rowSums(SA.Ref)))
KEEP <- MEAN > 0.01*MAX
CS <- colSums(SA.Curr)
RS <- colSums(SA.Ref)
NC <- sum(CS)
NR <- sum(RS)
if (!KEEP) {
    CS[] <- NA
    RS[] <- NA
}
SI <- if (KEEP)
    100 * min(NC, NR) / max(NC, NR) else NA

if (CI && KEEP) {
    if (LOCAL) {
        load(file.path(ROOT, VER, "results", taxon, "boot", paste0(spp, ".RData")))
    } else {
        con <- url(file.path(ROOT, VER, "results", taxon, "boot", paste0(spp, ".RData")))
        load(con)
        close(con)
    }
    Curr.Boot <- Curr.Boot[PIX10,]
    Ref.Boot <- Ref.Boot[PIX10,]
    Curr.Boot <- Curr.Boot[match(KTsub$Row10_Col10, rownames(Curr.Boot)),]
    Ref.Boot <- Ref.Boot[match(KTsub$Row10_Col10, rownames(Ref.Boot)),]

    CB <- colSums(Curr.Boot)
    RB <- colSums(Ref.Boot)

    NC_CI <- quantile(CB, a)
    NR_CI <- quantile(RB, a)
    SI_CI <- quantile(100 * pmin(CB, RB) / pmax(CB, RB), a)
} else {
    NC_CI <- c(NA, NA)
    names(NC_CI) <- paste0(100*a, "%")
    NR_CI <- SI_CI <- NC_CI
}

Sector_Total <- (100 * (CS - RS) / NR)[-1]
Sector_UnderHF <- (100 * (CS - RS) / RS)[-1]
KA <- if (taxon == "birds") A_2012 else A_2014
Sector_Area <- (100 * KA / sum(KA))[names(Sector_Total)]
Sector_Unit <- 100 * Sector_Total / Sector_Area

OUT <- list(
    taxon=taxon,
    species=spp,
    max=MAX,
    mean=MEAN,
    keep=KEEP,
    ptime=proc.time()-t0,
    intactness=rbind(
        Current=c(Estimate=NC, NC_CI),
        Reference=c(Estimate=NR, NR_CI),
        Intactness=c(Estimate=SI, SI_CI)),
    sector=rbind(
        Area=Sector_Area,
        Total=Sector_Total,
        UnderHF=Sector_UnderHF,
        Unit=Sector_Unit))
RESULTS[[taxon]][[spp]] <- OUT

rm(SA.Curr, SA.Ref, OUT)
if (CI)
    rm(Curr.Boot, Ref.Boot)
gc()

} # end loop for species
} # end loop for taxa

