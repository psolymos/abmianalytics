#remotes::install_github("ABbiodiversity/cure4insect", ref="v2018")

library(cure4insect)
library(sp)
library(rgdal)

opar <- set_options(path = "s:/reports") # this is optional if you have local copy
load_common_data()

## summarize spclim by NSR for DH

SPP <- get_all_species()
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))
xy <- kgrid
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, proj4string(rpa))

SOIL <- matrix(NA, length(SPP), nlevels(kgrid$NSRNAME))
dimnames(SOIL) <- list(SPP, levels(kgrid$NSRNAME))
VEG <- SOIL

for (spp in SPP) {
    #spp <- "AlderFlycatcher"
    cat(spp, which(SPP == spp), "/", length(SPP), "\n")
    flush.console()
    object <- load_spclim_data(spp)
    if (!is.null(object$rsoil)) {
        isoil <- extract(object$rsoil, xy, method="simple")
        ipa <- extract(rpa, xy, method="simple")
        vsoil <- object$caspen * ipa + isoil
        a <- aggregate(data.frame(soil=vsoil),
            list(NSR=kgrid$NSRNAME), mean, na.rm=TRUE)
        rownames(a) <- a[,1]
        SOIL[spp, rownames(a)] <- a$soil
    }
    if (!is.null(object$rveg)) {
        vveg <- extract(object$rveg, xy, method="simple")
        a <- aggregate(data.frame(veg=vveg),
            list(NSR=kgrid$NSRNAME), mean, na.rm=TRUE)
        rownames(a) <- a[,1]
        VEG[spp, rownames(a)] <- a$veg
    }
}
save(SOIL, VEG, file="~/Desktop/spclim-toDave-2019-09-17.RData", version=2)

## run species richness for website multi-spp

ID <- get_all_id()

for (tx in c("mammals", "birds", "mites", "mosses", "lichens", "vplants")) {
    subset_common_data(ID, tx)
    r <- make_multispecies_map("richness")
    writeRaster(r,
        file=paste0("d:/abmi/sppweb2018/multi-spp/richness-map-", tx, ".tif"))
}
