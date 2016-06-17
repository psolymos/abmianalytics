library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
source("~/repos/abmianalytics/R/maps_functions.R")

## mammals

dirin <- "e:\\peter\\sppweb2015-round2\\Mammals\\Km2 summaries\\"
files <- c("AllLeporids.csv",
    "MartenFisher.csv",
    "DomesticDog.csv",
    "Mink.csv",
    "Moose.csv",
    "Deer.csv",
    "ElkWapiti.csv",
    "Foxes.csv",
    "Coyote.csv",
    "GrayWolf.csv",
    "CanadaLynx.csv",
    "WeaselsAndErmine.csv",
    "RedSquirrel.csv")
SPP <- gsub(".csv", "", files)
fl <- paste0(dirin, files)
for (i in 1:length(fl)) {
    cat("Mammals:", SPP[i], "\n");flush.console()
    km <- read.csv(fl[i])
    km <- km[rowSums(is.na(km)) == 0,]
    #km <- km[match(kgrid$Row_Col, km$LinkID),]
    ndat <- normalize_data(rf=km$Ref, cr=km$Curr)
    ndat <- data.frame(LinkID=km$LinkID, ndat)
    write.csv(ndat, row.names=FALSE,
        file=paste0("w:\\normalized_data\\mammals\\maps\\", files[i]))
}

## birds

fln <- list.files(file.path(ROOT, "out", "birds", "results", "north"))
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- list.files(file.path(ROOT, "out", "birds", "results", "south"))
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)
SPP <- sort(union(fls, fln))

for (spp in SPP) {

    cat(spp, "\t");flush.console()
    load(file.path(ROOT, "out", "birds", "pred1cmb", paste0(spp, ".Rdata")))
    km <- data.frame(km)

    TYPE <- "C" # combo
    #if (!slt[spp, "veghf.north"])
    if (!(spp %in% fln))
        TYPE <- "S"
    #if (!slt[spp, "soilhf.south"])
    if (!(spp %in% fls))
        TYPE <- "N"

    wS <- 1-kgrid$pAspen
    if (TYPE == "S")
        wS[] <- 1
    if (TYPE == "N")
        wS[] <- 0
    wS[kgrid$useS] <- 1
    wS[kgrid$useN] <- 0

    cr <- wS * km$CurrS + (1-wS) * km$CurrN
    rf <- wS * km$RefS + (1-wS) * km$RefN
    ndat <- normalize_data(rf=rf, cr=cr)

}