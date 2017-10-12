library(mefa4)

ROOT <- "e:/peter/AB_data_v2017"

load(file.path(ROOT, "data", "analysis", "kgrid_table_km.Rdata"))
source("~/repos/abmianalytics/R/maps_functions.R")
Info0 <- readLines("~/repos/abmianalytics/mapping/README.md")

IC <- TRUE
WEB <- TRUE

#kgrid2 <- kgrid[,c("Row_Col", "POINT_X", "POINT_Y")]
#write.csv(kgrid2, row.names=FALSE, file="w:/gis/Grid1km_working.csv")

#taxon <- "lichens"
#taxon <- "vplants"
#taxon <- "mosses"
taxon <- "mites"
#taxon <- "mammals"
#taxon <- "birds"

if (taxon == "mammals") {
    ext <- ".csv" # csv or Rdata
    VER <- "5.0 (2017-07-13)" # version
    nam_col <- "species" # column used in README
    nam_in <- "sppid" # name used in input files
    nam_out <- "sppid" # name used in output files
}
if (taxon == "birds") {
    ext <- ".Rdata" # csv or Rdata
    VER <- "4.1 (2017-07-13)" # version
    nam_col <- "species" # column used in README
    nam_in <- "AOU"
    nam_out <- "sppid"
}
if (taxon %in% c("mites", "mosses", "lichens", "vplants")) {
    ext <- ".Rdata" # csv or Rdata
    VER <- "5.0 (2017-07-13)" # version
    nam_col <- "scinam" # column used in README
    nam_in <- "sppid" # name used in input files
    nam_out <- "sppid" # name used in output files
}

root_in <- "v:/contents/2017/species"
root_out <- "w:/species-2017"
spptab <- read.csv(file.path("~/repos/abmispecies/_data", paste0(taxon, ".csv")))
rownames(spptab) <- spptab[,nam_in]
sppid <- sort(rownames(spptab)[spptab$map.pred])

fexists <- sapply(file.path(root_in, taxon, "km2", paste0(sppid, ext)), file.exists)
table(fexists)
sppid[!fexists]
stopifnot(all(fexists))

for (i in 1:length(sppid)) {
    gc()
    spp <- sppid[i]
    sppout <- as.character(spptab[spp, nam_out])
    cat(taxon, spp, "\n");flush.console()
    if (ext == ".csv") {
        km <- read.csv(file.path(root_in, taxon, "km2", paste0(spp, ext)))
    } else {
        e <- new.env()
        load(file.path(root_in, taxon, "km2", paste0(spp, ext)), envir=e)
        if (taxon == "birds") {
            km <- data.frame(e$km2)
            km$LinkID <- rownames(km)
        } else {
            km <- e$RefCurr
        }
    }
    km <- km[rowSums(is.na(km)) == 0,]
    if (IC) {
        ndat <- normalize_data(rf=km$Ref, cr=km$Curr, q=0.99, normalize=TRUE)
        ndat <- data.frame(LinkID=km$LinkID, ndat)
        write.csv(ndat, row.names=FALSE,
            file=file.path(root_out, taxon, paste0(sppout, ".csv")))
    }
    if (WEB) {
        ndat <- normalize_data(rf=km$Ref, cr=km$Curr, q=0.99, normalize=FALSE)
        ndat <- data.frame(LinkID=km$LinkID, ndat)
        dir.create(file.path(root_out, taxon, sppout))
        write.csv(ndat, row.names=FALSE,
            file=file.path(root_out, taxon, sppout, paste0(sppout, ".csv")))
        Info <- gsub("&species&", spptab[spp, nam_col], Info0)
        Info <- gsub("&version&", VER, Info)
        writeLines(Info, file.path(root_out, taxon, sppout, "README.md"))

        setwd(file.path(root_out, taxon))
        zip(paste0(sppout, ".zip"),
            paste0("./", sppout, "/", c(paste0(sppout, ".csv"), "README.md")))
        unlink(file.path(root_out, taxon, sppout), recursive=TRUE, force=TRUE)
    }
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

## NN plants

f <- "e:/peter/sppweb2016/nn-plants/Combine regions/Km2 summaries/nNNSpp.csv"
km <- read.csv(f)

ndat <- normalize_data(rf=km$Ref, cr=km$Curr, q=0.99, normalize=FALSE)
ndat <- data.frame(LinkID=km$LinkID, ndat)
dir.create(file.path("e:/peter/sppweb2016/nn-plants", "nonnative-plants"))
write.csv(ndat, row.names=FALSE,
    file=file.path("e:/peter/sppweb2016/nn-plants", "nonnative-plants", "nonnative-plants.csv"))
Info <- gsub("&species&", "Non-native vascular plants", Info0)
Info <- gsub("&version&", "3.2", Info)
writeLines(Info, file.path("e:/peter/sppweb2016/nn-plants", "nonnative-plants", "README.md"))

setwd(file.path("e:/peter/sppweb2016/nn-plants"))
zip(paste0("nonnative-plants.zip"),
    paste0("./nonnative-plants/", c("nonnative-plants.csv", "README.md")))
unlink(file.path(root_out, taxon, sppout), recursive=TRUE, force=TRUE)
