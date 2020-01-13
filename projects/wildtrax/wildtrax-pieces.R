## Rasters

library(mefa4)
library(cure4insect)
set_options(verbose=0, path = "d:/abmi/reports")
load_common_data()


spt <- droplevels(get_species_table("birds"))

rt <- .read_raster_template()
KT <- cure4insect:::.c4if$KT

rasterize_results_cr <- function (y, current=TRUE)
{
    NC <- if (current)
        rowSums(y$SA.Curr) else rowSums(y$SA.Ref)
    KT$NC <- NC[match(rownames(KT), names(NC))]
    KT$NC[is.na(KT$NC)] <- 0
    r <- .make_raster(KT$NC, rc = KT, rt = rt)
    r
}


Species <- as.character(spt$SpeciesID)
for (spp in Species) {
    cat(spp, "\n")
    y <- load_species_data(spp)
    r <- rasterize_results_cr(y)
    p <- p_bird(r, "ha", 7) # no pair adjustment, just area
    writeRaster(p, paste0("s:/WildTrax/raster/wt-raster-", spp, ".tif"), overwrite=TRUE)
}
write.csv(spt, row.names=FALSE, file="s:/WildTrax/raster/birds-lookup.csv")

## phi & tau
