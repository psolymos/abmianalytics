#library(cure4insect)
library(mefa4)
library(intrval)

if (FALSE) {
## fix LUF regions
load("d:/abmi/reports/2017/data/kgrid_areas_by_sector.RData")

load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
all(rownames(kgrid) == rownames(KT))
all(kgrid$LUF_NAME == KT$reg_luf)
KT$reg_luf <- kgrid$LUF_NAME

VER$species <- as.numeric(table(SP$taxon)[rownames(VER)])

save(XY, KT, KA_2012, KA_2014, SP, QT2KT, VER, CF, CFbirds,
    file="d:/abmi/reports/2017/data/kgrid_areas_by_sector.RData")
save(XY, KT, KA_2012, KA_2014, SP, QT2KT, VER, CF, CFbirds,
    file="s:/reports/2017/data/kgrid_areas_by_sector.RData")

}

c4i0 <- new.env()
load("d:/abmi/reports/2017/data/kgrid_areas_by_sector.RData", envir=c4i0)

## same
XY <- c4i0$XY
str(XY)
QT2KT <- c4i0$QT2KT
str(QT2KT)

## lookup
str(c4i0$SP)
ROOT <- "~/GoogleWork/tmp/tables"
e1 <- new.env()
e2 <- new.env()
e3 <- new.env()
e4 <- new.env()
e5 <- new.env()
load(file.path(ROOT, paste0("StandardizedOutput-birds.RData")), envir=e1)
load(file.path(ROOT, paste0("StandardizedOutput-vplants.RData")), envir=e2)
load(file.path(ROOT, paste0("StandardizedOutput-mites.RData")), envir=e3)
load(file.path(ROOT, paste0("StandardizedOutput-mosses.RData")), envir=e4)
load(file.path(ROOT, paste0("StandardizedOutput-lichens.RData")), envir=e5)
tmp1 <- e1$Lookup
tmp1$taxon <- "birds"
tmp2 <- e2$Lookup
tmp2$taxon <- "vplants"
tmp3 <- e3$Lookup
tmp3$taxon <- "mites"
tmp4 <- e4$Lookup
tmp4$taxon <- "mosses"
tmp5 <- e5$Lookup
tmp5$taxon <- "lichens"

spt <- rbind(tmp1[,colnames(tmp2)], tmp2,
    tmp3[,colnames(tmp2)], tmp4[,colnames(tmp2)], tmp5[,colnames(tmp2)])
spt <- droplevels(spt[spt$ModelNorth | spt$ModelSouth, ])
rownames(spt) <- spt$SpeciesID
spt$native <- !spt$Nonnative
spt$model_north <- spt$ModelNorth
spt$model_south <- spt$ModelSouth
spt$Species <- ifelse(is.na(spt$CommonName), spt$ScientificName, spt$CommonName)
spt$habitat_assoc <- c4i0$SP$habitat_assoc[match(rownames(spt), rownames(c4i0$SP))]
spt$model_region <- factor("North and South", levels(c4i0$SP$model_region))
spt$model_region[!spt$model_north] <- "South"
spt$model_region[!spt$model_south] <- "North"
spt <- spt[,colnames(c4i0$SP)]
spt <- droplevels(rbind(spt, c4i0$SP[c4i0$SP$taxon == "mammals",]))

with(c4i0$SP, table(taxon, model_region))
with(spt, table(taxon, model_region))

with(spt, table(taxon, model_region)) - with(c4i0$SP, table(taxon, model_region))
SP <- spt

## version
VER <- c4i0$VER
VER$version <- 2018
VER$yr_last <- 2017
VER["mammals", "yr_last"] <- c4i0$VER["mammals", "yr_last"]
VER["mosses", "yr_last"] <- 2016
VER$hf <- "2016v3"
VER$veg <- "v6.1"
VER$species <- as.numeric(table(SP$taxon)[rownames(VER)])

c4i0$VER
VER

## regions
KT <- c4i0$KT

load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
all(rownames(kgrid) == rownames(KT))
all(kgrid$LUF_NAME == KT$reg_luf)
if (!all(kgrid$LUF_NAME == KT$reg_luf))
    KT$reg_luf <- kgrid$LUF_NAME

## HF by sector

load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(KT) == rownames(trVeg)))
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]

KA_2016 <- groupSums(trVeg[,rownames(ch2veg)], 2, ch2veg$sector)
#colnames(KA_2016)[colnames(KA_2016) == "Native"] <- "NATIVE"
KA_2016 <- KA_2016[,colnames(c4i0$KA_2014)]

KA_2016 <- KA_2016 / 10^6

rs <- rowSums(KA_2016)
rs[rs <= 1] <- 1
KA_2016 <- KA_2016 / rs

summary(rowSums(c4i0$KA_2012))
summary(rowSums(c4i0$KA_2014))
summary(rowSums(KA_2016))
sum(c4i0$KA_2012[,-1])/sum(c4i0$KA_2012)
sum(c4i0$KA_2014[,-1])/sum(c4i0$KA_2014)
sum(KA_2016[,-1])/sum(KA_2016)

## CF
load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")
colnames(c4i0$CF$coef$veg)
summary(c4i0$CF$coef$paspen)
pA <- as.matrix(OUT$pAspen)

Soil <- cbind(as.matrix(OUT$SoilhfSouthNontreed[,-(1:5)]),
    as.matrix(OUT$LinearSouth[,-(1:5)]))
Soil[is.na(Soil)] <- 0
SoilL <- Soil[,startsWith(colnames(Soil), "Lower_")]
SoilU <- Soil[,startsWith(colnames(Soil), "Upper_")]
Soil <- Soil[,!startsWith(colnames(Soil), "Lower_") & !startsWith(colnames(Soil), "Upper_")]
colnames(Soil)
colnames(c4i0$CF$coef$soil)
colnames(SoilL) <- gsub("Lower_", "", colnames(SoilL))
colnames(SoilL)
colnames(c4i0$CF$lower$soil)
colnames(SoilU) <- gsub("Upper_", "", colnames(SoilU))
colnames(SoilU)
colnames(c4i0$CF$higher$soil)

Veg <- cbind(as.matrix(OUT$VeghfNorth[,-(1:5)]),
    as.matrix(OUT$LinearNorth[,-(1:5)]))
Veg[is.na(Veg)] <- 0
VegL <- Veg[,startsWith(colnames(Veg), "Lower_")]
VegU <- Veg[,startsWith(colnames(Veg), "Upper_")]
Veg <- Veg[,!startsWith(colnames(Veg), "Lower_") & !startsWith(colnames(Veg), "Upper_")]
colnames(Veg)
colnames(c4i0$CF$coef$veg)
colnames(VegL) <- gsub("Lower_", "", colnames(VegL))
colnames(VegL)
colnames(c4i0$CF$lower$veg)
colnames(VegU) <- gsub("Upper_", "", colnames(VegU))
colnames(VegU)
colnames(c4i0$CF$higher$veg)


compare_sets(rownames(SP)[SP$model_south], rownames(Soil))
SP[rownames(SP) %ni% rownames(Soil) & SP$model_south,]

compare_sets(rownames(SP)[SP$model_north], rownames(Veg))
SP[rownames(SP) %ni% rownames(Veg) & SP$model_north,]

## no mammals included
CF <- list(
    coef=list(veg=Veg, soil=Soil, paspen=pA),
    lower=list(veg=VegL, soil=SoilL),
    higher=list(veg=VegU, soil=SoilU))

## save
write.csv(spt, row.names=FALSE, file="d:/abmi/reports/2018/data/species-info.csv")
save(XY, KT, KA_2016, SP, QT2KT, VER, CF, # CFbirds,
    file="d:/abmi/reports/2018/data/kgrid_areas_by_sector.RData")

write.csv(spt, row.names=FALSE, file="s:/reports/2018/data/species-info.csv")
save(XY, KT, KA_2016, SP, QT2KT, VER, CF, # CFbirds,
    file="s:/reports/2018/data/kgrid_areas_by_sector.RData")


## calculate spclim raster for /spclim

library(mefa4)
library(raster)
source("~/repos/abmianalytics/birds/00-functions.R")
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
#ROOT <- "~/GoogleWork/tmp"

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

cfs <- list(
    spclim=c("pWater_KM",
        "pWater2_KM", "xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
        "xMCMT", "xY", "xX", "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY"))

## kgrid
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y

xclim <- data.frame(
    transform_clim(kgrid),
    pAspen=kgrid$pAspen,
    pWater_KM=kgrid$pWater,
    pWater2_KM=kgrid$pWater^2)
## this has pAspen for the south, otherwise all the same
Xclim <- model.matrix(as.formula(paste0("~-1+", paste(cfs$spclim, collapse="+"))), xclim)
colnames(Xclim) <- fix_names(colnames(Xclim))

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

## birds
SPP <- rownames(SP)[SP$taxon == "birds"]
AOU <- as.character(e1$Lookup[SPP,"Code"])
for (i in 1:length(AOU)) {
    spp <- AOU[i]
    cat(spp, "\n");flush.console()

    TYPE <- "C" # combo
    if (SP[SPP[i], "model_south"] && !SP[SPP[i], "model_north"])
        TYPE <- "S"
    if (!SP[SPP[i], "model_south"] && SP[SPP[i], "model_north"])
        TYPE <- "N"

    if (TYPE != "N") {
        ests <- e1$CoefSouthBootlist[[spp]][1,cfs$spclim]
        musClim <- drop(Xclim[,cfs$spclim] %*% ests[cfs$spclim])
        rsoil <- make_raster(musClim, kgrid, rt)
    } else {
        rsoil <- NULL
    }

    if (TYPE != "S") {
        estn <- e1$CoefNorthBootlist[[spp]][1,cfs$spclim]
        munClim <- drop(Xclim[,cfs$spclim] %*% estn[cfn$spclim])
        rveg <- make_raster(munClim, kgrid, rt)
    } else {
        rveg <- NULL
    }

    save(rveg, rsoil, file=paste0("d:/abmi/reports/2018/results/birds/spclim/",
        SPP[i], ".RData"))

}


cn <- colnames(e2$SpclimNorth)
compare_sets(cn, colnames(kgrid))
kgrid$Intercept <- 1
kgrid$Lat <- kgrid$POINT_Y
kgrid$Long <- kgrid$POINT_X
kgrid$Lat2 <- kgrid$Lat^2
kgrid$Long2 <- kgrid$Long^2
kgrid$LatLong <- kgrid$Lat*kgrid$Long
kgrid$MAPPET <- kgrid$MAP*kgrid$PET
kgrid$MATAHM <- kgrid$MAT*kgrid$AHM
kgrid$MAPFFP <- kgrid$MAP*kgrid$FFP
kgrid$MAT2 <- kgrid$MAT^2
kgrid$MWMT2 <- kgrid$MWMT^2
setdiff(cn, colnames(kgrid))
Xclim <- as.matrix(kgrid[,cn])

tx <- "vplants"
for (tx in c("vplants", "mites", "mosses", "lichens")) {
    ee <- switch(tx,
        "vplants"=e2,
        "mites"=e3,
        "mosses"=e4,
        "lichens"=e5)
    SPP <- rownames(SP)[SP$taxon == tx]
    for (spp in SPP) {
        cat(tx, spp, "\n");flush.console()

        TYPE <- "C" # combo
        if (SP[spp, "model_south"] && !SP[spp, "model_north"])
            TYPE <- "S"
        if (!SP[spp, "model_south"] && SP[spp, "model_north"])
            TYPE <- "N"

        if (TYPE != "N") {
            ests <- ee$SpclimSouth[spp,,1]
            ests <- ests[names(ests) != "pAspen"]
            musClim <- drop(Xclim[,names(ests)] %*% ests)
            rsoil <- make_raster(musClim, kgrid, rt)
        } else {
            rsoil <- NULL
        }

        if (TYPE != "S") {
            estn <- ee$SpclimNorth[spp,,1]
            munClim <- drop(Xclim[,names(estn)] %*% estn)
            rveg <- make_raster(munClim, kgrid, rt)
        } else {
            rveg <- NULL
        }

        save(rveg, rsoil, file=paste0("d:/abmi/reports/2018/results/", tx, "/spclim/",
            spp, ".RData"))

    }
}

## /sector files

tx <- "vplants"
for (tx in c("vplants", "mites", "mosses", "lichens")) {

    dirin <- paste0("s:/Result from Ermias_2018/", tx, "/combined/Sector effects/Sector abundance summary/")

    SPP <- rownames(SP)[SP$taxon == tx]
    for (spp in SPP) {
        cat(tx, spp, "\n");flush.console()
        ee <- new.env()
        load(paste0(dirin, spp, ".RData"), envir=ee)
        SA.Curr <- as.matrix(ee$SA.curr[,colnames(KA_2016)])
        SA.Ref <- as.matrix(ee$SA.ref[,colnames(KA_2016)])

        save(SA.Curr, SA.Ref, file=paste0("d:/abmi/reports/2018/results/", tx, "/sector/",
            spp, ".RData"))
    }
}

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$ModelNorth | Lookup$ModelSouth,])
rownames(tax) <- tax$Code
SPP <- as.character(tax$SpeciesID)
AOU <- rownames(tax)
CN <- c("Native", "Misc", "Agriculture", "Forestry", "RuralUrban", "Energy", "Transportation")

for (i in 1:length(AOU)) {
    spp <- AOU[i]
    cat(spp, "\n");flush.console()

    ee <- new.env()
    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-04-01/", spp, ".RData"), envir=ee)
    SA.Curr <- as.matrix(ee$Curr[,CN])
    SA.Ref <- as.matrix(ee$Ref[,CN])

    save(SA.Curr, SA.Ref, file=paste0("d:/abmi/reports/2018/results/birds/sector/",
        SPP[i], ".RData"))

}

## compare species lists

library(cure4insect)
set_options(path = "s:/reports")
set_options(version = 2017)
load_common_data()
SP1 <- get_species_table()
clear_common_data()
set_options(version = 2018)
load_common_data()
SP2 <- get_species_table()


## use package to do sector effects updates

devtools::install_github("ABbiodiversity/cure4insect@v2018")
library(cure4insect)
set_options(path = "s:/reports")
set_options(version = 2018)
#clear_common_data()
load_common_data()


KT <- cure4insect:::.c4if$KT
KT$N <- KT$reg_nr != "Grassland" & KT$reg_nr != "Rocky Mountain" &
    KT$reg_nr != "Parkland" & KT$reg_nsr != "Dry Mixedwood"
KT$S <- KT$reg_nr == "Grassland" | KT$reg_nr == "Parkland" |
    KT$reg_nsr == "Dry Mixedwood"


ID <- rownames(KT)[KT$N]
Spp <- get_all_species()
subset_common_data(id=ID, species=Spp)
xxn <- report_all(cores=8)
resn <- do.call(rbind, lapply(xxn, flatten))
class(resn) <- c("data.frame", "c4idf")
#y <- load_species_data("AlderFlycatcher")
#x <- calculate_results(y)
#flatten(x)

ID <- rownames(KT)[KT$S]
subset_common_data(id=ID, species=Spp)
xxs <- report_all(cores=8)
ress <- do.call(rbind, lapply(xxs, flatten))
class(ress) <- c("data.frame", "c4idf")

save(resn, ress, file="d:/abmi/sppweb2018/c4i/tables/sector-effects.RData")


load("d:/abmi/sppweb2018/c4i/tables/sector-effects.RData")

SPP <- rownames(resn)

for (spp in SPP) {
    TYPE <- "C"
    if (resn[spp, "model_north"] && !resn[spp, "model_south"])
        TYPE <- "N"
    if (!resn[spp, "model_north"] && resn[spp, "model_south"])
        TYPE <- "S"
    if (TYPE != "S") {
        png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/figs/sector/", spp, "-north.png"),
            height=500, width=1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(resn[spp,], "unit")
        plot_sector(resn[spp,], "regional", main="")
        plot_sector(resn[spp,], "underhf", main="")
        par(op)
        dev.off()
    }
    if (TYPE != "S") {
        png(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/figs/sector/", spp, "-south.png"),
            height=500, width=1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(ress[spp,], "unit")
        plot_sector(ress[spp,], "regional", main="")
        plot_sector(ress[spp,], "underhf", main="")
        par(op)
        dev.off()

    }
}

