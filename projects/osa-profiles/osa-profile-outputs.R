## Generate OSA profile figures

#devtools::install_github("ABbiodiversity/cure4insect")
library(cure4insect)
library(knitr)
library(rgdal)
set_options(verbose=0, path = "w:/reports")
load_common_data()

setwd("~")

Species <- c("Ovenbird",
    "BlackthroatedGreenWarbler",
    "CanadaWarbler",
    "BrownCreeper",
    "BaybreastedWarbler")

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837")))(100)
col3 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221"))(100)
#ply <- readOGR(dsn=system.file("extdata/OSA_bound.geojson", package="cure4insect"))
ply <- readOGR(dsn="e:/peter/AB_data_v2018/data/raw/xy/Oilsands-Boundaries.gdb",
    "OilsandRegionDissolve10TM")
ID <- overlay_polygon(ply)
## write IDs into a text file
#write.table(data.frame(SpatialID=ID), row.names=FALSE, file="SpatialID.txt")
AB <- readOGR(dsn=system.file("extdata/AB_bound.geojson", package="cure4insect"))
AB <- spTransform(AB, proj4string(r))
ply <- spTransform(ply, proj4string(r))

#species <- "Ovenbird"
#species <- "BlackthroatedGreenWarbler"
#species <- "CanadaWarbler"
#species <- "BrownCreeper"
#species <- "BaybreastedWarbler"
for (species in Species) {

    cat(species, "\n"); flush.console()
    info <- as.list(get_species_table()[species, c("SpeciesID", "CommonName", "ScientificName", "TSNID")])
    info <- lapply(info, as.character)
    y <- load_species_data(species)
    r <- rasterize_results(y)
    subset_common_data(id=NULL, species=species)
    x <- calculate_results(y)

    ## calculate regional stats
    subset_common_data(id=ID, species)
    xreg <- calculate_results(y)
    ## clip raster
    rreg <- crop(r, ply)
    rreg <- mask(rreg, ply)
    Max <- max(values(rreg[["NC"]]), values(rreg[["NR"]]), na.rm=TRUE)
    df <- (rreg[["NC"]] - rreg[["NR"]]) / Max
    df <- sign(df) * abs(df)^0.5
    df <- 100*df
    df[!is.na(values(df)) & values(df) > 100] <- 100
    Rng <- range(values(df), na.rm=TRUE)
    df[!is.na(values(df)) & values(df) == Rng[1]] <- -100
    df[!is.na(values(df)) & values(df) == Rng[2]] <- 100

    ## what we need for each of the species
    ## - AB detection map (from spp web)
    ## - veg-hf barplot
    ## - sector effects: total (in region)
    ## - sector effects: under hf (in region)
    ## - current abundance map, in region
    ## - difference map, in region


    png(paste0("osa-profiles-", species, "-habitat.png"),
        width=1200, height=600)
    plot_abundance(species, type="veg_coef")
    dev.off()

    png(paste0("osa-profiles-", species, "-sector-total.png"),
        width=500, height=500)
    plot_sector(xreg, type="regional", main=info$CommonName)
    dev.off()

    png(paste0("osa-profiles-", species, "-sector-underhf.png"),
        width=500, height=500)
    plot_sector(xreg, type="underhf", info$CommonName)
    dev.off()

    png(paste0("osa-profiles-", species, "-map-current-abundance.png"),
        width=500, height=500)
    plot(rreg[["NC"]], col=col1, axes=FALSE, box=FALSE, main="Current abundance in OSA")
    dev.off()

    png(paste0("osa-profiles-", species, "-map-difference.png"),
        width=500, height=500)
    plot(df, col=col3, axes=FALSE, box=FALSE, main="Difference in OSA")
    dev.off()
}


op <- par(mfrow=c(1,3), mar=c(1,1,1,1))
plot(r[["NC"]], col=col1, axes=FALSE, box=FALSE, main="Current abundance")
plot(r[["SE"]], col=col2, axes=FALSE, box=FALSE, main="Stadard error")
plot(r[["SI2"]], col=col3, axes=FALSE, box=FALSE, main="Species intactness")
par(op)

op <- par(mfrow=c(2,2), mar=2*c(1,1,1,1))
plot(rreg[["NC"]], col=col1, axes=FALSE, box=FALSE, main="Current abundance in OSA")
plot(rreg[["SE"]], col=col2, axes=FALSE, box=FALSE, main="Stadard error in OSA")
plot(rreg[["SI2"]], col=col3, axes=FALSE, box=FALSE, main="Species intactness in OSA")
plot(df, col=col3, axes=FALSE, box=FALSE, main="Difference in OSA")
par(op)


op <- par(mfrow=c(1,2))
par(op)

#df3 <- data.frame(xreg$intactness)
#colnames(df3)[2:3] <- c("Lower", "Upper")
#kable(df3, digits=2, caption="Table 1. Species intactness results in the OSA.")

#df4 <- data.frame(xreg$sector)[,-1]
#kable(df4, digits=2, caption="Table 2. Sector effects results in the OSA.")


