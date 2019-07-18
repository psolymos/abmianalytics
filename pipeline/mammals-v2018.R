#remotes::install_github("ABbiodiversity/cure4insect", ref="v2018")
library(cure4insect)
library(mefa4)
set_options(path = "s:/reports")
load_common_data()

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y

## ch2soil ch2veg trSoil trVeg
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]

compare_sets(ch2soil$cr, rownames(ts))
setdiff(ch2soil$cr, rownames(ts))
setdiff(rownames(ts), ch2soil$cr)

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2soil$rf2 <- ts$UseInAnalysisCoarse[match(ch2soil$rf, rownames(ts))]
ch2soil$cr2 <- ts$UseInAnalysisCoarse[match(ch2soil$cr, rownames(ts))]
ch2soil$sector <- ts$Sector61[match(ch2soil$cr, rownames(ts))]

ch2veg$rf2 <- tv$CoefTabs[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$CoefTabs[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

str(ch2soil)
str(ch2veg)

EXCL <- c("HWater", "SoilUnknown", "SoilWater", "Water", "HFor")
keeps <- rownames(ch2soil)[!(ch2soil$cr2 %in% EXCL) & !(ch2soil$rf2 %in% EXCL) ]
keepn <- rownames(ch2veg)[!(ch2veg$cr2 %in% EXCL) & !(ch2veg$rf2 %in% EXCL)]
trSoil <- trSoil[,keeps]
trVeg <- trVeg[,keepn]
ch2soil <- ch2soil[keeps,]
ch2veg <- ch2veg[keepn,]

rss <- rowSums(trSoil)
rss[rss==0] <- 1
trSoil <- trSoil / rss
rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn

summary(ch2soil)
summary(ch2veg)

stopifnot(all(rownames(trVeg) == rownames(kgrid)))
stopifnot(all(rownames(trSoil) == rownames(kgrid)))

## lat/long for centroids
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X+POINT_Y
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))

## pAspen: also part of kgrid
xy <- spTransform(xy, proj4string(rpa))
pAspen <- extract(rpa, xy, "bilinear")
summary(pAspen)
o <- 1:length(pAspen)
nas <- which(is.na(pAspen))
for (i in nas) {
    d <- colSums((drop(coordinates(xy[i,])) - t(coordinates(xy)))^2)
    d[i] <- Inf
    d[nas] <- Inf
    j <- which.min(d)
    o[i] <- o[j]
}
xy <- kgrid[o,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X+POINT_Y
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
xy <- spTransform(xy, proj4string(rpa))

pAspen <- extract(rpa, xy, "bilinear")
summary(pAspen)

kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

## checking levels
compare_sets(ch2veg$cr2, get_levels()$veg)
compare_sets(ch2soil$cr2, get_levels()$soil)
setdiff(ch2soil$cr2, get_levels()$soil)
setdiff(get_levels()$soil, ch2soil$cr2)

## species list
spt <- get_species_table("mammals")
(SPP <- get_all_species("mammals"))
SECT <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban",
    "Energy", "Transportation")

OUT <- list()

## species
for(spp in SPP) {
    cat("\n", spp, "\t")
    object <- load_spclim_data(spp)

    SA.Curr <- SA.Ref <- matrix(NA, nrow(kgrid), length(SECT), dimnames=list(rownames(kgrid), SECT))

    ## sector
    for (s in SECT) {
        cat(" ", s)
        flush.console()

        TYPE <- "C" # combo
        if (spt[spp, "model_south"] && !spt[spp, "model_north"])
            TYPE <- "S"
        if (!spt[spp, "model_south"] && spt[spp, "model_north"])
            TYPE <- "N"
        if (s == "Forestry") # no forestry in south
            TYPE <- "N"

        if (TYPE != "N") {
            soil <- trSoil[,ch2soil$sector == s]
            soilCr <- groupSums(soil, 2, as.character(ch2soil$cr2[ch2soil$sector == s]))
            soilRf <- groupSums(soil, 2, as.character(ch2soil$rf2[ch2soil$sector == s]))
        } else {
            soilCr <- soilRf <- NULL
        }
        if (TYPE != "S") {
            veg <- trVeg[,ch2veg$sector == s]
            vegCr <- groupSums(veg, 2, as.character(ch2veg$cr2[ch2veg$sector == s]))
            vegRf <- groupSums(veg, 2, as.character(ch2veg$rf2[ch2veg$sector == s]))
        } else {
            vegCr <- vegRf <- NULL
        }

        prCr <- predict_mat.c4ispclim(object, xy, vegCr, soilCr, method="bilinear")
        prRf <- predict_mat.c4ispclim(object, xy, vegRf, soilRf, method="bilinear")

        wS <- 1-kgrid$pAspen
        if (TYPE == "S")
            wS[] <- 1
        if (TYPE == "N")
            wS[] <- 0
        wS[kgrid$useS] <- 1
        wS[kgrid$useN] <- 0

        if (TYPE != "N") {
            NsoilCr <- rowSums(prCr$soil)
            NsoilRf <- rowSums(prRf$soil)
        } else {
            NsoilCr <- 0
            NsoilRf <- 0
        }
        if (TYPE != "S") {
            NvegCr <- rowSums(prCr$veg)
            NvegRf <- rowSums(prRf$veg)
        } else {
            NvegCr <- 0
            NvegRf <- 0
        }

        Curr <- wS * NsoilCr + (1-wS) * NvegCr
        Ref <- wS * NsoilRf + (1-wS) * NvegRf

        SA.Curr[,s] <- Curr
        SA.Ref[,s] <- Ref
    }
    save(SA.Curr, SA.Ref,
        file=paste0("d:/abmi/reports/2018/results/mammals/sector/", spp, ".RData"))
    #OUT[[spp]] <- list(cr=SA.Curr, rf=SA.Ref)
}

## comparing the 2 versions

library(cure4insect)
library(mefa4)

SPP <- c("AllLeporids", "CanadaLynx", "Coyote", "Deer", "ElkWapiti",
    "Foxes", "GrayWolf", "MartenFisher", "Mink", "Moose", "RedSquirrel",
    "WeaselsAndErmine")

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837")))(100)

#spp <- SPP[3]
for (spp in SPP) {

    clear_common_data()
    set_options(path = "s:/reports", version="2017")
    load_common_data()
    y <- load_species_data(spp)
    r17 <- rasterize_results(y)

    clear_common_data()
    set_options(path = "s:/reports", version="2018")
    load_common_data()
    y <- load_species_data(spp)
    r18 <- rasterize_results(y)

    df1 <- r18[["NC"]] - r17[["NC"]]
    df0 <- r18[["NR"]] - r17[["NR"]]

    png(paste0("d:/abmi/reports/", spp, "-1718.png"), width=1500,height=1000)
    op <- par(mfrow=c(2,3), mar=c(2,2,2,2))
    plot(r17, "NC", main="Curr v2017", col=col1, axes=FALSE, box=FALSE)
    plot(r18, "NC", main="Curr v2018", col=col1, axes=FALSE, box=FALSE)
    plot(df1, main="Curr diff", col=col2, axes=FALSE, box=FALSE)
    plot(r17, "NR", main="Ref v2017", col=col1, axes=FALSE, box=FALSE)
    plot(r18, "NR", main="Ref v2018", col=col1, axes=FALSE, box=FALSE)
    plot(df0, "NC", main="Ref diff", col=col2, axes=FALSE, box=FALSE)
    par(op)
    dev.off()
}


