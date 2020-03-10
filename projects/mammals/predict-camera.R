## summaries for camera mammals for data portal
library(mefa4)
library(raster)
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2019/misc/mammals"
load(file.path(ROOT,
    "Mammal coefficients North Feb 2019 Best model OFFICIAL coefficients ORIGINAL NAMES.Rdata"))

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y

load("d:/abmi/AB_data_v2019/misc/overlap/OverlapReg.RData")
rownames(OverlapReg) <- OverlapReg$LinkID
OverlapReg$pAspen <- kgrid[rownames(OverlapReg), "pAspen"]
OverlapReg$wN <- OverlapReg$pAspen / (OverlapReg$pAspen + (1-OverlapReg$pForest))
kgrid$wN <- ifelse(kgrid$NRNAME == "Grassland", 0, 1)
kgrid[rownames(OverlapReg), "wN"] <- OverlapReg$wN

xclim <- data.frame(
    transform_clim(kgrid),
    pAspen=kgrid$pAspen,
    pWater_KM=kgrid$pWater,
    pWater2_KM=kgrid$pWater^2)
xclim$NSR1Parkland <- as.integer(kgrid$NSRNAME=="Central Parkland" |
        kgrid$NSRNAME=="Foothills Parkland" | kgrid$NSRNAME=="Peace River Parkland")
xclim$NSR1DryMixedwood <- as.integer(kgrid$NSRNAME=="Dry Mixedwood")
xclim$NSR1CentralMixedwood <- as.integer(kgrid$NSRNAME=="Central Mixedwood")
xclim$NSR1Foothills <- as.integer(kgrid$NSRNAME=="Lower Foothills" |
        kgrid$NSRNAME=="Upper Foothills")
xclim$NSR1North <- as.integer(kgrid$NSRNAME=="Lower Boreal Highlands" |
        kgrid$NSRNAME=="Upper Boreal Highlands" | kgrid$NSRNAME=="Boreal Subarctic" |
        kgrid$NSRNAME=="Northern Mixedwood")
xclim$NSR1Shield <- as.integer(kgrid$NSRNAME=="Kazan Uplands" |
        kgrid$NSRNAME=="Peace-Athabasca Delta" | kgrid$NSRNAME=="Athabasca Plain")
xclim$NSR1Mountain <- as.integer(kgrid$NSRNAME=="Montane" |
        kgrid$NSRNAME=="Subalpine" | kgrid$NSRNAME=="Alpine")
xclim$Lat <- kgrid$POINT_Y
xclim$Long <- kgrid$POINT_X
xclim <- data.frame(xclim, kgrid[,c("AHM", "PET", "FFP", "MAP", "MAT", "MCMT", "MWMT")])
#xclim$Intercept<-1
xclim$LatLong<-xclim$Lat*xclim$Long
xclim$Lat2<-xclim$Lat*xclim$Lat
xclim$Lat3<-xclim$Lat*xclim$Lat*xclim$Lat
xclim$Long2<-xclim$Long*xclim$Long
xclim$Lat2Long2<-xclim$Lat*xclim$Lat*xclim$Long*xclim$Long
xclim$MAPFFP<-xclim$MAP*xclim$FFP
xclim$MAPPET<-xclim$MAP*xclim$PET
xclim$MATAHM<-xclim$MAT*xclim$AHM
xclim$MWMT2<-xclim$MWMT*xclim$MWMT
xclim$MAT2<-xclim$MAT*xclim$MAT
xclim$LongMAT<-xclim$Long*xclim$MAT
## this has pAspen for the south, otherwise all the same
Xclim <- model.matrix(~., xclim)
colnames(Xclim)[1] <- "Intercept"
Xclim <- Xclim[,colnames(Res.coef.official)]


## ch2soil ch2veg trSoil trVeg
#load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v61hf2016v3WildFireUpTo2016.Rdata")
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
ch2veg$rf2 <- tv$UseInAnalysisFineAge[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$UseInAnalysisFineAge[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

str(ch2soil)
str(ch2veg)

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn

CN <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")

SPP <- rownames(Coef.official)

for (spp in SPP) {
    gc()
    cat(spp)
    flush.console()

    estnclim <- Res.coef.official[spp,]
    munClim <- drop(Xclim[,names(estnclim)] %*% estnclim)
    munClim[munClim > 2] <- 2

    tmp <- Coef.official[spp,]
    names(tmp)[names(tmp) == "CultivationCrop"] <- "Crop"
    names(tmp)[names(tmp) == "CultivationTamePasture"] <- "TameP"
    names(tmp)[names(tmp) == "CultivationRoughPasture"] <- "RoughP"
    names(tmp)[names(tmp) == ""] <- "TreeShrubSwamp"
    #names(tmp)[names(tmp) == "RuralResidentialIndustrial"] <- "Urban"
    names(tmp)[names(tmp) == "WellSite"] <- "Well"
    names(tmp)[names(tmp) == "UrbanIndustrial"] <- "Industrial"
    names(tmp)[names(tmp) == "TreeShrubSwamp"] <- "Swamp"

    munHab <- structure(rep(NA, nlevels(ch2veg$cr2)), names=levels(ch2veg$cr2))
    munHab[intersect(names(tmp), names(munHab))] <- tmp[intersect(names(tmp), names(munHab))]

    names(tmp) <- gsub("TreedBog", "BSpr", names(tmp))
    munHab[names(tmp)[grep("BSpr", names(tmp))]] <- tmp[names(tmp)[grep("BSpr", names(tmp))]]
    munHab[grep("Larch", names(munHab))] <- tmp["TreedFen"]
    munHab[c("Seismic", "TrSoftLin", "EnSoftLin")] <- tmp["SoftLin"]
    munHab[c("GraminoidFen", "Marsh")] <- tmp["NonTreeFenMarsh"]
    munHab[c("Urban", "Rural")] <- tmp["RuralResidentialIndustrial"]

    munHab["BSpr9"] <- munHab["BSpr8"]
    munHab["Decid9"] <- munHab["Decid8"]
    munHab["Mixedwood9"] <- munHab["Mixedwood8"]
    munHab["Pine9"] <- munHab["Pine8"]
    munHab["Spruce9"] <- munHab["Spruce8"]

    munHab[is.na(munHab)] <- 0 # water, snow/ice, mine

    ## expand coefficients for north
    prnCr <- munHab[match(ch2veg$cr2, names(munHab))]
    prnRf <- munHab[match(ch2veg$rf2, names(munHab))]
    ## put pieces together for north
    ## multiplying with row normalized area gives the weighted averaging
    ADnCr <- t(prnCr * t(trVeg)) * exp(munClim)
    ADnRf <- t(prnRf * t(trVeg)) * exp(munClim)
    ## add up by sector for north
    ADnCrSect <- groupSums(ADnCr, 2, ch2veg$sector)
    ADnRfSect <- groupSums(ADnRf, 2, ch2veg$sector)
    Curr <- ADnCrSect[,CN]
    Ref <- ADnRfSect[,CN]

    save(Curr, Ref,
        file=paste0(ROOT, "/pred/", spp, ".RData"))

    cat(" - DONE\n")
}


## mapping

library(mefa4)
library(raster)

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

## checking results

cn <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")
rn <- rownames(kgrid)[kgrid$NRNAME != "Grassland"]
Aveg <- groupSums(trVeg, 2, ch2veg$sector)[,cn]
KA <- 100*colMeans(Aveg[rn,])

SEff <- list()

for (spp in SPP) {

    cat(spp, "\n")
    flush.console()

    load(paste0(ROOT, "/pred/", spp, ".RData"))

    CS <- colSums(Curr[rn,cn])
    RS <- colSums(Ref[rn,cn])
    NC <- sum(CS)
    NR <- sum(RS)
    Sector_Total <- (100 * (CS - RS) / NR)[-1]
    Sector_UnderHF <- (100 * (CS - RS) / RS)[-1]
    Sector_Area <- (100 * KA / sum(KA))[names(Sector_Total)]
    Sector_Unit <- 100 * Sector_Total / Sector_Area

#    round(cbind(CS=CS, RS=RS, Df=CS-RS), 2)
#    round(cbind(Total=Sector_Total,
#        Under=Sector_UnderHF,
#        Unit=Sector_Unit), 2)
    SEff[[spp]] <- cbind(Total=Sector_Total,
        Under=Sector_UnderHF,
        Unit=Sector_Unit)

    Dcr <- rowSums(Curr)
    q <- quantile(Dcr, 0.99)
    Dcr[Dcr > q] <- q
    summary(Dcr)
    Drf <- rowSums(Ref)
    q <- quantile(Drf, 0.99)
    Drf[Drf > q] <- q
    summary(Drf)
    MAX <- max(Dcr, Drf)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)

    if (FALSE) {
        tmp <- data.frame(
            ID=rownames(kgrid),
            Current=Dcr,
            Reference=Drf)
            #Color_Current=col1[cr],
            #Color_Reference=col1[rf],
            #Color_Difference=col3[df]
        tmp <- tmp[kgrid$useN,]
#        dbWriteTable(con, "test_num", tmp,
#            overwrite=TRUE, row.names=FALSE)
    }

    if (TRUE) {
        Rcr <- make_raster(cr, kgrid, rt)
        Rrf <- make_raster(rf, kgrid, rt)
        Rdf <- make_raster(df-100, kgrid, rt)
        Msk <- Rmaskn
        Rcr <- mask(Rcr, Msk)
        Rrf <- mask(Rrf, Msk)
        Rdf <- mask(Rdf, Msk)
        ## add here mask for Rockies if needed

        png(paste0(ROOT, "/maps/", spp, ".png"),
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
}


## compare

library(cure4insect)
set_options(path = "d:/abmi/reports")
load_common_data()

info <- droplevels(get_species_table("birds"))

for (spp in rownames(tax)) {
    if (tax[spp, "SpeciesID"] %in% info$SpeciesID) {
        cat(spp, "\n");flush.console()
        species <- as.character(tax[spp, "SpeciesID"])
        y <- load_species_data(species)
        TYPE <- "C" # combo
        if (info[species, "model_south"] && !info[species, "model_north"])
            TYPE <- "S"
        if (!info[species, "model_south"] && info[species, "model_north"])
            TYPE <- "N"

        Dcr <- rowSums(y$SA.Curr[match(rownames(kgrid), rownames(y$SA.Curr)),])
        Dcr[is.na(Dcr)] <- 0
        q <- quantile(Dcr, 0.99)
        Dcr[Dcr > q] <- q
        summary(Dcr)
        Drf <- rowSums(y$SA.Ref[match(rownames(kgrid), rownames(y$SA.Ref)),])
        Drf[is.na(Drf)] <- 0
        q <- quantile(Drf, 0.99)
        Drf[Drf > q] <- q
        summary(Drf)
        MAX <- max(Dcr, Drf)

        df <- (Dcr-Drf) / MAX
        df <- sign(df) * abs(df)^0.5
        df <- pmin(200, ceiling(99 * df)+100)
        df[df==0] <- 1
        cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
        rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
        #si <- 100 * pmin(Dcr, Drf)/pmax(Dcr, Drf)
        #si[is.na(si)] <- 100
        #si[si==0] <- 1

        Rcr <- make_raster(cr, kgrid, rt)
        Rrf <- make_raster(rf, kgrid, rt)
        Rdf <- make_raster(df-100, kgrid, rt)
        #Rsi <- make_raster(si, kgrid, rt)
        if (TYPE == "S")
            Msk <- Rmasks
        if (TYPE == "N")
            Msk <- Rmaskn
        if (TYPE != "C") {
            Rcr <- mask(Rcr, Msk)
            Rrf <- mask(Rrf, Msk)
            Rdf <- mask(Rdf, Msk)
            #Rsi <- mask(Rsi, Msk)
        }
        ## add here mask for Rockies if needed

        png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/figs/maps/00-compare-", spp, ".png"),
            height=1500*2, width=1000*3, res=300)
        op <- par(mfrow=c(2,3), mar=c(2,1,2,3))
        plot(Rcr, col=col1[1:max(cr)], axes=FALSE, box=FALSE, main="Current old", legend=FALSE)
        plot(Rrf, col=col1[1:max(rf)], axes=FALSE, box=FALSE, main="Reference old", legend=FALSE)
        plot(Rdf, col=col3[min(df):max(df)], axes=FALSE, box=FALSE, main="Difference old", legend=FALSE)




        load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-04-01/", spp, ".RData"))
        TYPE <- "C" # combo
        if (tax[spp, "ModelSouth"] && !tax[spp, "ModelNorth"])
            TYPE <- "S"
        if (!tax[spp, "ModelSouth"] && tax[spp, "ModelNorth"])
            TYPE <- "N"

        Dcr <- rowSums(Curr)
        q <- quantile(Dcr, 0.99)
        Dcr[Dcr > q] <- q
        summary(Dcr)
        Drf <- rowSums(Ref)
        q <- quantile(Drf, 0.99)
        Drf[Drf > q] <- q
        summary(Drf)
        MAX <- max(Dcr, Drf)

        df <- (Dcr-Drf) / MAX
        df <- sign(df) * abs(df)^0.5
        df <- pmin(200, ceiling(99 * df)+100)
        df[df==0] <- 1
        cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
        rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
        #si <- 100 * pmin(Dcr, Drf)/pmax(Dcr, Drf)
        #si[is.na(si)] <- 100
        #si[si==0] <- 1

        Rcr <- make_raster(cr, kgrid, rt)
        Rrf <- make_raster(rf, kgrid, rt)
        Rdf <- make_raster(df-100, kgrid, rt)
        #Rsi <- make_raster(si, kgrid, rt)
        if (TYPE == "S")
            Msk <- Rmasks
        if (TYPE == "N")
            Msk <- Rmaskn
        if (TYPE != "C") {
            Rcr <- mask(Rcr, Msk)
            Rrf <- mask(Rrf, Msk)
            Rdf <- mask(Rdf, Msk)
            #Rsi <- mask(Rsi, Msk)
        }
        ## add here mask for Rockies if needed

        plot(Rcr, col=col1[1:max(cr)], axes=FALSE, box=FALSE, main="Current new", legend=FALSE)
        plot(Rrf, col=col1[1:max(rf)], axes=FALSE, box=FALSE, main="Reference new", legend=FALSE)
        plot(Rdf, col=col3[min(df):max(df)], axes=FALSE, box=FALSE, main="Difference new", legend=FALSE)

        par(op)
        dev.off()


    }
}

## sector effects

library(mefa4)
library(raster)

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth | Lookup$ModelSouth,])
rownames(tax) <- tax$Code


## species tables for data portal

library(cure4insect)
set_options(path="d:/abmi/reports")
load_common_data()
SP <- get_species_table()
SPP <- rownames(SP)
KT <- get_id_table()

for (spp in SPP) {

    cat(spp, which(SPP==spp), "/", length(SPP), "\n")
    flush.console()

    tx <- as.character(SP[spp, "taxon"])

    TYPE <- "C" # combo
    if (SP[spp, "model_south"] && !SP[spp, "model_north"])
        TYPE <- "S"
    if (!SP[spp, "model_south"] && SP[spp, "model_north"])
        TYPE <- "N"

    y <- load_species_data(spp, boot=FALSE)
    Dcr <- rowSums(y$SA.Curr[match(rownames(KT), rownames(y$SA.Curr)),])
    q <- quantile(Dcr, 0.99, na.rm=TRUE)
    Dcr[!is.na(Dcr) & Dcr > q] <- q
    Drf <- rowSums(y$SA.Ref[match(rownames(KT), rownames(y$SA.Ref)),])
    q <- quantile(Drf, 0.99, na.rm=TRUE)
    Drf[!is.na(Drf) & Drf > q] <- q
    MAX <- max(Dcr, Drf, na.rm=TRUE)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
    d <- data.frame(Row_Col=rownames(KT), Curr=cr, Ref=rf, Diff=df)
    i <- !is.na(cr) & !is.na(rf)
    d <- d[i,]
    write.csv(d, row.names = FALSE,
        file=paste0("s:/BDQT/species-tables-1km/", tx, "/", spp, ".csv"))
}


