devtools::install_github("ABbiodiversity/cure4insect@v2018")

library(cure4insect)
library(mefa4)
set_options(path = "s:/reports")
load_common_data()

load("d:/abmi/sppweb2018/c4i/tables/sector-effects.RData")
ROOT <- "d:/abmi/sppweb2018/www/"

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

Rmaskn <- make_raster(as.integer(1-kgrid$useS), kgrid, rt)
values(Rmaskn)[values(Rmaskn) == 0] <- NA
Rmasks <- make_raster(as.integer(1-kgrid$useN), kgrid, rt)
values(Rmasks)[values(Rmasks) == 0] <- NA
#Rmaskm <- make_raster(as.integer(kgrid$NRNAME == "Rocky Mountain"), kgrid, rt)
#values(Rmaskm)[values(Rmaskm) == 0] <- NA
Rw <- make_raster(as.integer(kgrid$pWater > 0.99), kgrid, rt)
values(Rw)[values(Rw) == 0] <- NA

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
col3 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221"))(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude



gr <- "birds"
#gr <- "vplants"
#gr <- "lichens"
#gr <- "mosses"
#gr <- "mites"
SPP <- rownames(resn[resn$Taxon == gr,])
#spp <- "AlderFlycatcher"

for (spp in SPP) {

    cat(gr, spp, "\n");flush.console()

    TYPE <- "C"
    if (resn[spp, "model_north"] && !resn[spp, "model_south"])
        TYPE <- "N"
    if (!resn[spp, "model_north"] && resn[spp, "model_south"])
        TYPE <- "S"

    if (TYPE != "S") {
        ## hab-north
        png(paste0(ROOT, "/", gr, "/", spp, "-coef-north.png"),
            height=700, width=1500, res=150)
        layout(matrix(c(1,1,2), nrow=1))
        plot_abundance(spp, "veg_coef")
        par(mar=c(12,5,4,3))
        plot_abundance(spp, "veg_lin", main="")
        dev.off()

        png(paste0(ROOT, "/", gr, "/", spp, "-sector-north.png"),
            height=2*500, width=2*1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(resn[spp,], "unit")
        plot_sector(resn[spp,], "regional", main="")
        plot_sector(resn[spp,], "underhf", main="")
        par(op)
        dev.off()
    }
    if (TYPE != "N") {
        ## hab-south
        png(paste0(ROOT, "/", gr, "/", spp, "-coef-south.png"),
            height=700, width=1500, res=150)
        layout(matrix(c(1,2,3), nrow=1))
        p1 <- plot_abundance(spp, "soil_coef", paspen=0, plot=FALSE)
        p2 <- plot_abundance(spp, "soil_coef", paspen=1, plot=FALSE)
        plot_abundance(spp, "soil_coef", paspen=0, ylim=c(0, max(p1, p2)))
        plot_abundance(spp, "soil_coef", paspen=1, ylim=c(0, max(p1, p2)), main="")
        par(mar=c(11,5,4,3))
        plot_abundance(spp, "soil_lin", main="")
        dev.off()

        png(paste0(ROOT, "/", gr, "/", spp, "-sector-south.png"),
            height=2*500, width=2*1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(ress[spp,], "unit")
        plot_sector(ress[spp,], "regional", main="")
        plot_sector(ress[spp,], "underhf", main="")
        par(op)
        dev.off()
    }


    y <- load_species_data(spp)
    Curr <- y$SA.Curr[match(rownames(kgrid), rownames(y$SA.Curr)),]
    Ref <- y$SA.Ref[match(rownames(kgrid), rownames(y$SA.Ref)),]
    Dcr <- rowSums(Curr)
    q <- quantile(Dcr, 0.99)
    Dcr[Dcr > q] <- q
    Drf <- rowSums(Ref)
    q <- quantile(Drf, 0.99)
    Drf[Drf > q] <- q
    MAX <- max(Dcr, Drf)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)

    Rcr <- make_raster(cr, kgrid, rt)
    Rrf <- make_raster(rf, kgrid, rt)
    Rdf <- make_raster(df-100, kgrid, rt)
    if (TYPE == "S")
        Msk <- Rmasks
    if (TYPE == "N")
        Msk <- Rmaskn
    if (TYPE != "C") {
        Rcr <- mask(Rcr, Msk)
        Rrf <- mask(Rrf, Msk)
        Rdf <- mask(Rdf, Msk)
    }
    ## add here mask for Rockies if needed

    png(paste0(ROOT, "/", gr, "/", spp, "-map.png"),
        height=1500*1, width=1000*3, res=300)
    op <- par(mfrow=c(1,3), mar=c(2,1,2,3))
    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Reference", legend=FALSE)
    plot(Rrf, add=TRUE, col=col1[1:max(rf)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)
    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Current", legend=FALSE)
    plot(Rcr, add=TRUE, col=col1[1:max(cr)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)
    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Difference", legend=FALSE)
    plot(Rdf, add=TRUE, col=col3[min(df):max(df)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)
    par(op)
    dev.off()



}
