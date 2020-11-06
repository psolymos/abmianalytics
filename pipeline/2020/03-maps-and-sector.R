## mapping and sector effects
library(mefa4)

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")

load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]

## common functions and variables

## organize space-climate variables using kgrid table
## not sure if truncated latitude is needed: S study are excluded >56.7
make_clim <- function(x, birds=FALSE, truncate_latitude=FALSE) {
    if (birds) {
        z <- with(x, cbind(
            pWater_KM=pWater,
            pWater2_KM=pWater^2,
            xPET=(PET - 0) / 800,
            xMAT=(MAT - 0) / 6,
            xAHM=(AHM - 0) / 50,
            xFFP=(FFP - 0) / 130,
            xMAP=(MAP - 0) / 2300,
            xMWMT=(MWMT - 0) / 20,
            xMCMT=(MCMT - 0) / 25,
            xY=(POINT_Y - 54.1) / 2.43,
            xX=(POINT_X - (-113.3)) / 2.12))
        z <- cbind(z,
            xY2=z[,"xY"]^2,
            xX2=z[,"xX"]^2,
            `xFFP:xMAP`=z[,"xFFP"]*z[,"xMAP"],
            `xMAP:xPET`=z[,"xMAP"]*z[,"xPET"],
            `xAHM:xMAT`=z[,"xAHM"]*z[,"xMAT"],
            `xX:xY`=z[,"xX"]*z[,"xY"])
    } else {
        LAT <- pmin(x$POINT_Y, 56.5)
        z <- with(x, cbind(
            Intercept=1,
            Lat=LAT,
            Long=x$POINT_X,
            AHM=x$AHM,
            PET=x$PET,
            FFP=x$FFP,
            MAP=x$MAP,
            MAT=x$MAT,
            MCMT=x$MCMT,
            MWMT=x$MWMT,
            Lat2=LAT^2,
            Long2=x$POINT_X^2,
            LatLong=x$POINT_X*LAT,
            MAPPET=x$MAP*x$PET,
            MATAHM=x$MAT*x$AHM,
            MAPFFP=x$MAP*x$FFP,
            MAT2=x$MAT^2,
            MWMT2=x$MWMT^2))
    }
    rownames(z) <- rownames(x)
    z
}
## making space-climate model matrix
Xclim_bird <- make_clim(kgrid, birds=TRUE)
Xclim_nonb <- make_clim(kgrid, birds=FALSE)

kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

UseN <- rownames(kgrid)[!kgrid$useS]
UseS <- rownames(kgrid)[!kgrid$useN]

Xclim_bird_S <- Xclim_bird[UseS,]
Xclim_nonb_S <- Xclim_nonb[UseS,]
pA <- kgrid[UseS, "pAspen"]
Xclim_bird_N <- Xclim_bird[UseN,]
Xclim_nonb_N <- Xclim_nonb[UseN,]

## weights for overlap region
if (FALSE) {
load("d:/abmi/AB_data_v2019/misc/overlap/OverlapReg.RData")
rownames(OverlapReg) <- OverlapReg$LinkID
OverlapReg$pAspen <- kgrid[rownames(OverlapReg), "pAspen"]
OverlapReg$wN <- OverlapReg$pAspen / (OverlapReg$pAspen + (1-OverlapReg$pForest))

kgrid$wN <- ifelse(kgrid$NRNAME == "Grassland", 0, 1)
kgrid[rownames(OverlapReg), "wN"] <- OverlapReg$wN
kgrid$wN[kgrid$useN] <- 1
summary(kgrid$wN[kgrid$useS]) # S only
summary(kgrid$wN[kgrid$useN]) # N only
summary(kgrid$wN[!kgrid$useS & !kgrid$useN]) # overlap
}

FORE <- c("CCDecid1", "CCDecid2",
    "CCDecid3", "CCDecid4", "CCDecidR", "CCMixedwood1", "CCMixedwood2",
    "CCMixedwood3", "CCMixedwood4", "CCMixedwoodR", "CCPine1", "CCPine2",
    "CCPine3", "CCPine4", "CCPineR", "CCSpruce1", "CCSpruce2", "CCSpruce3",
    "CCSpruce4", "CCSpruceR",
    "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6", "Decid7", "Decid8",
    "DecidR", "Mixedwood1", "Mixedwood2",
    "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7",
    "Mixedwood8", "MixedwoodR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
    "Pine7", "Pine8", "PineR",
    "Spruce1", "Spruce2", "Spruce3", "Spruce4", "Spruce5", "Spruce6",
    "Spruce7", "Spruce8", "SpruceR", "TreedBog1",
    "TreedBog2", "TreedBog3", "TreedBog4", "TreedBog5", "TreedBog6",
    "TreedBog7", "TreedBog8", "TreedBogR", "TreedFen1", "TreedFen2",
    "TreedFen3", "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7",
    "TreedFen8", "TreedFenR", "TreedSwamp1", "TreedSwamp2", "TreedSwamp3",
    "TreedSwamp4", "TreedSwamp5", "TreedSwamp6", "TreedSwamp7", "TreedSwamp8",
    "TreedSwampR")
UNKN <- chSoil$rf=="UNK" | chSoil$cr=="UNK"
kgrid$pSoilUnk <- rowSums(trSoil[,UNKN]) / rowSums(trSoil)
kgrid$pFor <- rowSums(trVeg[,chVeg$cr %in% FORE]) / rowSums(trVeg)

kgrid$wS <- (1-kgrid$pAspen)*(1-kgrid$pSoilUnk)*(1-kgrid$pFor)
kgrid$wS[kgrid$useN] <- 0
kgrid$wS[kgrid$useS] <- 1
kgrid$wN <- kgrid$pAspen / (kgrid$pAspen + (1-kgrid$pFor))
kgrid$wN[kgrid$useN] <- 1
kgrid$wN[kgrid$useS] <- 0
wsm <- kgrid$wS + kgrid$wN
kgrid$wS <- kgrid$wS / wsm
kgrid$wN <- kgrid$wN / wsm
hist(kgrid$wN[!kgrid$useS & !kgrid$useN]) # overlap


## lookup tables
lt <- list(
    south=structure(list(Label = c("ClaySub", "Other", "Other", "Other",
        "Other", "Other", "RapidDrain", "Loamy", "Other", "Other", "Other",
        "Other", "Other", "Other", "ClaySub", "SandyLoam", "RapidDrain",
        "RapidDrain", "RapidDrain", "RapidDrain", "RapidDrain", "ThinBreak",
        "Blowout", "Other", "Other", "Blowout", "UNK", "Water", "Urban",
        "Urban", "Rural", "Industrial", "Industrial", "Rural", "Mine",
        "Mine", "Wellsites", "EnSoftLin", "EnSoftLin", "EnSeismic", "EnSeismic",
        "HardLin", "HardLin", "TrSoftLin", "TrSoftLin", "TrSoftLin",
        "Crop", "RoughP", "RoughP", "TameP", "Industrial", "Water", "Water",
        "Water", "Water", "UNK"), Sector = c("Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "RuralUrban", "RuralUrban",
        "RuralUrban", "RuralUrban", "Energy", "Misc", "Energy", "Misc",
        "Energy", "Energy", "Energy", "Energy", "Energy", "Transportation",
        "Transportation", "Transportation", "Transportation", "Transportation",
        "Agriculture", "Agriculture", "Agriculture", "Agriculture", "Agriculture",
        "Misc", "Misc", "Misc", "Misc", "Forestry")), row.names = c("Cy",
        "Len", "LenS", "LenSP", "LenT", "LenW", "Li", "Lo", "Ltc", "LtcC",
        "LtcD", "LtcH", "LtcS", "LtcR", "Sb", "Sy", "BdL", "CS", "Gr",
        "Sa", "SwG", "TB", "BlO", "LenA", "Ov", "SL", "UNK", "Water",
        "UrbanIndustrial", "UrbanResidence", "RuralResidentialIndustrial",
        "IndustrialSiteRural", "WindGenerationFacility", "OtherDisturbedVegetation",
        "MineSite", "PeatMine", "WellSite", "Pipeline", "TransmissionLine",
        "SeismicLineNarrow", "SeismicLineWide", "RoadHardSurface", "RailHardSurface",
        "RoadTrailVegetated", "RoadVegetatedVerge", "RailVegetatedVerge",
        "CultivationCrop", "CultivationAbandoned", "CultivationRoughPasture",
        "CultivationTamePasture", "HighDensityLivestockOperation", "BorrowpitsDugoutsSumps",
        "MunicipalWaterSewage", "Reservoirs", "Canals", "CutBlocks"), class = "data.frame"),
    north=structure(list(Label = c("DeciduousR", "Deciduous1", "Deciduous2",
        "Deciduous3", "Deciduous4", "Deciduous5", "Deciduous6", "Deciduous7",
        "Deciduous8", "MixedwoodR", "Mixedwood1", "Mixedwood2", "Mixedwood3",
        "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7", "Mixedwood8",
        "PineR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
        "Pine7", "Pine8", "WhiteSpruceR", "WhiteSpruce1", "WhiteSpruce2",
        "WhiteSpruce3", "WhiteSpruce4", "WhiteSpruce5", "WhiteSpruce6",
        "WhiteSpruce7", "WhiteSpruce8", "TreedBogR", "TreedBog1", "TreedBog2",
        "TreedBog3", "TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7",
        "TreedBog8", "TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3",
        "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8",
        "TreedSwamp", "TreedSwamp", "TreedSwamp", "TreedSwamp", "TreedSwamp",
        "TreedSwamp", "TreedSwamp", "TreedSwamp", "TreedSwamp", "GrassHerb",
        "Shrub", "GraminoidFen", "Marsh", "ShrubbyBog", "ShrubbyFen",
        "ShrubbySwamp", "Water", "Urban", "Urban", "Rural", "Industrial",
        "Industrial", "Rural", "Mine", "Mine", "Wellsites", "EnSoftLin",
        "EnSoftLin", "EnSeismic", "EnSeismic", "HardLin", "HardLin",
        "TrSoftLin", "TrSoftLin", "TrSoftLin", "Crop", "RoughP", "RoughP",
        "TameP", "TameP", "CCDeciduousR", "CCDeciduous1", "CCDeciduous2",
        "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR", "CCMixedwood1",
        "CCMixedwood2", "CCMixedwood3", "CCMixedwood4", "CCPineR", "CCPine1",
        "CCPine2", "CCPine3", "CCPine4", "CCWhiteSpruceR", "CCWhiteSpruce1",
        "CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "Bare",
        "Water", "Water", "Water", "Water", "SnowIce"), Sector = c("Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "Native", "Native", "Native", "Native", "Native", "Native", "Native",
        "RuralUrban", "RuralUrban", "RuralUrban", "RuralUrban", "Energy",
        "Misc", "Energy", "Misc", "Energy", "Energy", "Energy", "Energy",
        "Energy", "Transportation", "Transportation", "Transportation",
        "Transportation", "Transportation", "Agriculture", "Agriculture",
        "Agriculture", "Agriculture", "Agriculture", "Forestry", "Forestry",
        "Forestry", "Forestry", "Forestry", "Forestry", "Forestry", "Forestry",
        "Forestry", "Forestry", "Forestry", "Forestry", "Forestry", "Forestry",
        "Forestry", "Forestry", "Forestry", "Forestry", "Forestry", "Forestry",
        "Native", "Misc", "Misc", "Misc", "Misc", "Native")), row.names = c("DecidR",
        "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6", "Decid7",
        "Decid8", "MixedwoodR", "Mixedwood1", "Mixedwood2", "Mixedwood3",
        "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7", "Mixedwood8",
        "PineR", "Pine1", "Pine2", "Pine3", "Pine4", "Pine5", "Pine6",
        "Pine7", "Pine8", "SpruceR", "Spruce1", "Spruce2", "Spruce3",
        "Spruce4", "Spruce5", "Spruce6", "Spruce7", "Spruce8", "TreedBogR",
        "TreedBog1", "TreedBog2", "TreedBog3", "TreedBog4", "TreedBog5",
        "TreedBog6", "TreedBog7", "TreedBog8", "TreedFenR", "TreedFen1",
        "TreedFen2", "TreedFen3", "TreedFen4", "TreedFen5", "TreedFen6",
        "TreedFen7", "TreedFen8", "TreedSwampR", "TreedSwamp1", "TreedSwamp2",
        "TreedSwamp3", "TreedSwamp4", "TreedSwamp5", "TreedSwamp6", "TreedSwamp7",
        "TreedSwamp8", "GrassHerb", "Shrub", "GraminoidFen", "Marsh",
        "ShrubbyBog", "ShrubbyFen", "ShrubbySwamp", "Water", "UrbanIndustrial",
        "UrbanResidence", "RuralResidentialIndustrial", "IndustrialSiteRural",
        "WindGenerationFacility", "OtherDisturbedVegetation", "MineSite",
        "PeatMine", "WellSite", "Pipeline", "TransmissionLine", "SeismicLineNarrow",
        "SeismicLineWide", "RoadHardSurface", "RailHardSurface", "RoadTrailVegetated",
        "RoadVegetatedVerge", "RailVegetatedVerge", "CultivationCrop",
        "CultivationAbandoned", "CultivationRoughPasture", "CultivationTamePasture",
        "HighDensityLivestockOperation", "CCDecidR", "CCDecid1", "CCDecid2",
        "CCDecid3", "CCDecid4", "CCMixedwoodR", "CCMixedwood1", "CCMixedwood2",
        "CCMixedwood3", "CCMixedwood4", "CCPineR", "CCPine1", "CCPine2",
        "CCPine3", "CCPine4", "CCSpruceR", "CCSpruce1", "CCSpruce2",
        "CCSpruce3", "CCSpruce4", "Bare", "BorrowpitsDugoutsSumps", "Canals",
        "MunicipalWaterSewage", "Reservoirs", "SnowIce"), class = "data.frame"))

## define sectors to be used here
## i.e. change to attribution stuff or keep finer level transitions etc
chVeg$sector_use <- chVeg$sector
chSoil$sector_use <- chSoil$sector


## process south monster matrix and find some efficiencies

compare_sets(chSoil$cr, rownames(lt$south))
compare_sets(chSoil$rf, rownames(lt$south))
chSoil$cr2 <- lt$south$Label[match(chSoil$cr, rownames(lt$south))]
chSoil$rf2 <- lt$south$Label[match(chSoil$rf, rownames(lt$south))]
#chSoil$sector2 <- lt$south$Sector[match(chSoil$cr, rownames(lt$south))]
#with(chSoil, table(sector, sector2)) # sector definition is up to date

# aggregating the monster matrix
chSoil$tr2 <- paste0(chSoil$rf2, "->", chSoil$cr2)
chSoil$tr2[chSoil$rf2 == chSoil$cr2] <- chSoil$rf2[chSoil$rf2 == chSoil$cr2]
chSoil$tr2[chSoil$cr2=="UNK" | chSoil$rf2=="UNK"] <- "UNK"
length(unique(chSoil$tr2))/nrow(chSoil) # savings!

table(rn=grepl("UNK", rownames(chSoil)), tr2=chSoil$tr2=="UNK")

trSoil <- groupSums(trSoil, 2, chSoil$tr2)
chSoil <- nonDuplicated(chSoil, tr2, TRUE)[colnames(trSoil),c("sector_use","cr2", "rf2")]
## either soil is unknown or Sector is Forestry --> exclude these
s <- chSoil$cr2=="UNK" | chSoil$rf2=="UNK" | rownames(chSoil) == "UNK"
chSoil[s,]
trSoil <- trSoil[,!s]
chSoil <- chSoil[!s,]

## reference=water is not part of the landbase, so does not count for averaging
s <- chSoil$rf2=="Water" | chSoil$cr2=="Water"
chSoil[s,]
## but some of this is not Water->Water: we cannot attribute and it is implausible
sum(trSoil[,s])/sum(trSoil) - sum(trSoil[,chSoil$rf2=="Water" & chSoil$cr2=="Water"])/sum(trSoil)
## so we drop this ~3% together with open water
trSoil <- trSoil[,!s]
chSoil <- chSoil[!s,]

## now we make sure rows sum to 1 (or 0)
rs <- rowSums(trSoil)
Ps <- trSoil / ifelse(rs > 0, rs, 1)
## take subset that contains only the south study region
Ps <- Ps[UseS,]
summary(rowSums(Ps))

## trVeg processing

compare_sets(chVeg$cr, rownames(lt$north))
setdiff(chVeg$cr, rownames(lt$north))
setdiff(rownames(lt$north), chVeg$cr)
compare_sets(chVeg$rf, rownames(lt$north))
chVeg$cr2 <- lt$north$Label[match(chVeg$cr, rownames(lt$north))]
chVeg$rf2 <- lt$north$Label[match(chVeg$rf, rownames(lt$north))]
#chVeg$sector2 <- lt$north$Sector[match(chVeg$cr, rownames(lt$north))]
#with(chVeg, table(sector, sector2)) # sector definition is up to date

# aggregating the monster matrix
chVeg$tr2 <- paste0(chVeg$rf2, "->", chVeg$cr2)
chVeg$tr2[chVeg$rf2 == chVeg$cr2] <- chVeg$rf2[chVeg$rf2 == chVeg$cr2]
chVeg$tr2[chVeg$cr2=="UNK" | chVeg$rf2=="UNK"] <- "UNK"
length(unique(chVeg$tr2))/nrow(chVeg) # savings!

trVeg <- groupSums(trVeg, 2, chVeg$tr2)
chVeg <- nonDuplicated(chVeg, tr2, TRUE)[colnames(trVeg),c("sector_use","cr2", "rf2")]

## reference=water is not part of the landbase, so does not count for averaging
s <- chVeg$rf2=="Water" | chVeg$cr2=="Water"
chVeg[s,]
## but some of this is not Water->Water: we cannot attribute and it is implausible
sum(trVeg[,s])/sum(trVeg) - sum(trVeg[,chVeg$rf2=="Water" & chVeg$cr2=="Water"])/sum(trVeg)
## so we drop this ~3% together with open water
trVeg <- trVeg[,!s]
chVeg <- chVeg[!s,]

## now we make sure rows sum to 1 (or 0)
rn <- rowSums(trVeg)
Pn <- trVeg / ifelse(rn > 0, rn, 1)
## take subset that contains only the north study region
Pn <- Pn[UseN,]
summary(rowSums(Pn))


## need to make Pscr/Psrf and Pncr/Pnrf based on cr*sector and rf*sector
## this will reduce number of columns considerably
chSoil$cr2s <- paste0(chSoil$cr2, "/", chSoil$sector_use)
chSoil$rf2s <- paste0(chSoil$rf2, "/", chSoil$sector_use)
Psrf <- groupSums(Ps, 2, chSoil$rf2s)
chSrf <- nonDuplicated(chSoil, rf2s, TRUE)[colnames(Psrf), c("rf2", "sector_use")]
Pscr <- groupSums(Ps, 2, chSoil$cr2s)
chScr <- nonDuplicated(chSoil, cr2s, TRUE)[colnames(Pscr), c("cr2", "sector_use")]

chVeg$cr2s <- paste0(chVeg$cr2, "/", chVeg$sector_use)
chVeg$rf2s <- paste0(chVeg$rf2, "/", chVeg$sector_use)
Pnrf <- groupSums(Pn, 2, chVeg$rf2s)
chVrf <- nonDuplicated(chVeg, rf2s, TRUE)[colnames(Pnrf), c("rf2", "sector_use")]
Pncr <- groupSums(Pn, 2, chVeg$cr2s)
chVcr <- nonDuplicated(chVeg, cr2s, TRUE)[colnames(Pncr), c("cr2", "sector_use")]

sapply(chSoil,function(z) nlevels(droplevels(as.factor(z))))/ncol(Ps)
sapply(chVeg,function(z) nlevels(droplevels(as.factor(z))))/ncol(Pn)

## this is the current landscape only
Pncr2 <- groupSums(Pn, 2, chVeg$cr2)
Pscr2 <- groupSums(Ps, 2, chSoil$cr2)

## species specific part begins here

#spp <- "AlderFlycatcher"
#taxon <- "birds"

#spp <- "Actaea.rubra"
#taxon <- "vplants"

BOOT <- FALSE # do bootstrap?
SEFF <- TRUE  # non bootstrapped map/SE stuff
BMAX <- 100

ROOT <- "s:/AB_data_v2020/Results/pred"
ROOT2 <- "s:/AB_data_v2020/Results/boot"

COEFS$mammals <- COEFS2$mammals
COEFS$habitats <- COEFS2$habitats

TAXA <- names(COEFS)
for (taxon in TAXA) {

    if (SEFF && !dir.exists(file.path(ROOT, taxon)))
        dir.create(file.path(ROOT, taxon))
    if (BOOT && !dir.exists(file.path(ROOT2, taxon)))
        dir.create(file.path(ROOT2, taxon))

    SPPn <- if (taxon!="birds")
        dimnames(COEFS[[taxon]]$north)[[1]] else dimnames(COEFS[[taxon]]$north$joint)[[1]]
    SPPs <- if (taxon!="birds")
        dimnames(COEFS[[taxon]]$south)[[1]] else dimnames(COEFS[[taxon]]$south$joint)[[1]]
    SPP <- sort(union(SPPn, SPPs))

    for (spp in SPP) {

        cat(taxon, spp)
        flush.console()

        type <- "C" # combo species (N+S)
        M <- list(N=spp %in% SPPn, S=spp %in% SPPs)
        if (M$N & !M$S)
            type <- "N"
        if (!M$N & M$S)
            type <- "S"
        if (taxon == "birds") {
            cfn <- if (type == "S")
                NULL else COEFS[[taxon]]$north$joint[spp,,]
            cfs <- if (type == "N")
                NULL else COEFS[[taxon]]$south$joint[spp,,]
            XclimS <- Xclim_bird_S
            XclimN <- Xclim_bird_N
            #FUN <- poisson()$linkinv
            FUN <- function (eta)
                pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)
        } else {
            cfn <- if (type == "S")
                NULL else COEFS[[taxon]]$north[spp,,]
            cfs <- if (type == "N")
                NULL else COEFS[[taxon]]$south[spp,,]
            XclimS <- Xclim_nonb_S
            XclimN <- Xclim_nonb_N
            FUN <- binomial()$linkinv
        }
        if (taxon == "mammals")
            FUN <- function (eta)
                pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)
        if (taxon == "habitats")
            FUN <- binomial(COEFS[[taxon]]$species[spp, "link"])$linkinv


        ## bootstrap for current map
        ## note: it is still slow. wait and run final set on westgid
        if (BOOT) {
            for (i in 1:BMAX) {
                t0 <- proc.time()
                if (type != "N") {
                    gc()
                    ## south calculations for the i'th run
                    bscr <- cfs[colnames(Pscr2), i] # current land cover
                    ## space-climate coefs
                    bscl <- if (taxon == "birds")
                        cfs[colnames(Xclim_bird), i] else cfs[colnames(Xclim_nonb), i]
                    bspa <- cfs["pAspen", i]
                    ## additive components for south
                    muscl <- drop(XclimS %*% bscl) + pA * bspa
                    if (taxon == "birds") {
                        # exponential link is 3-4x faster without creating the matrix
                        muscr <- t(t(Pscr2) * FUN(bscr)) * FUN(muscl)
                        NScr <- rowSums(muscr)
                    } else {
                        muscr <- matrix(muscl, nrow=nrow(Pscr2), ncol=ncol(Pscr2))
                        muscr <- t(t(muscr) + bscr)
                        NScr <- rowSums(Pscr2 * FUN(muscr))
                        # using outer is surprisingly less efficient than transpose
                        #muscr <- outer(muscl, bscr, FUN="+")
                        #NScr <- rowSums(Pscr2 * FUN(muscr))

                    }
                } else {
                    NScr <- NULL
                }
                if (type != "S") {
                    gc()
                    ## north calculations for the i'th run
                    tmpn <- c(cfn[,i], Bare=-10^4, SnowIce= -10^4)
                    bncr <- tmpn[colnames(Pncr2)] # current land cover
                    ## space-climate coefs
                    bncl <- if (taxon == "birds")
                        cfn[colnames(Xclim_bird), i] else cfn[colnames(Xclim_nonb), i]
                    ## additive components for north
                    muncl <- drop(XclimN %*% bncl)
                    if (taxon == "birds") {
                        muncr <- t(t(Pncr2) * FUN(bncr)) * FUN(muncl)
                        NNcr <- rowSums(muncr)
                    } else {
                        muncr <- matrix(muncl, nrow=nrow(Pncr2), ncol=ncol(Pncr2))
                        muncr <- t(t(muncr) + bncr)
                        NNcr <- rowSums(Pncr2 * FUN(muncr))
                    }
                } else {
                    NNcr <- NULL
                }
                cat("\t", proc.time()[3]-t0[3], "\n")
            }
        }

        ## only doing 1 run for sector effects
        if (SEFF) {
            i <- 1
            t0 <- proc.time()
            if (type != "N") {
                gc()
                ## south calculations for the i'th run
                #compare_sets(rownames(cfs), chSoil$cr2)
                bscr <- cfs[chScr$cr2, i] # current land cover
                bsrf <- cfs[chSrf$rf2, i] # reference land cover
                ## space-climate coefs
                bscl <- if (taxon == "birds")
                    cfs[colnames(Xclim_bird), i] else cfs[colnames(Xclim_nonb), i]
                bscl[is.na(bscl)] <- 0 # this happens for habitat elements
                bspa <- cfs["pAspen", i]
                ## additive components for south
                if (taxon=="mammals") {
                    ## space/clim includes presence-absence piece
                    #Total.Abundance.approx<- exp(log(TA.vegHF) + SC.pOcc + pAspen.pa + pAspen.agp)
                    pApa <- COEFS[[taxon]]$pAspenPA[spp]
                    muscl <- drop(cbind(XclimS, pA) %*% c(bscl, pApa))
                    if (spp == "Pronghorn")
                        muscl <- (muscl - mean(muscl))/1000
                } else {
                    muscl <- drop(XclimS %*% bscl)
                }
                muspa <- pA * bspa

                muscr <- matrix(muscl + muspa, nrow=nrow(Pscr), ncol=ncol(Pscr))
                muscr <- t(t(muscr) + bscr)
                NScr <- as.matrix(groupSums(Pscr * FUN(muscr), 2, chScr$sector_use))
                NScr <- cbind(NScr, Forestry=0)

                musrf <- matrix(muscl + muspa, nrow=nrow(Psrf), ncol=ncol(Psrf))
                musrf <- t(t(musrf) + bsrf)
                NSrf <- as.matrix(groupSums(Psrf * FUN(musrf), 2, chSrf$sector_use))
                NSrf <- cbind(NSrf, Forestry=0)
            } else {
                NScr <- NULL
                NSrf <- NULL
            }
            if (type != "S") {
                gc()
                ## north calculations for the i'th run
                #compare_sets(rownames(cfn), chVeg$cr2)
                tmpn <- c(cfn[,i], Bare=-10^4, SnowIce= -10^4)
                bncr <- tmpn[chVcr$cr2] # current land cover
                bnrf <- tmpn[chVrf$rf2] # reference land cover
                ## space-climate coefs
                bncl <- if (taxon == "birds")
                    cfn[colnames(Xclim_bird), i] else cfn[colnames(Xclim_nonb), i]
                bncl[is.na(bncl)] <- 0 # this happens for habitat elements
                ## additive components for north
                muncl <- drop(XclimN %*% bncl)

                muncr <- matrix(muncl, nrow=nrow(Pncr), ncol=ncol(Pncr))
                muncr <- t(t(muncr) + bncr)
                NNcr <- as.matrix(groupSums(Pncr * FUN(muncr), 2, chVcr$sector_use))

                munrf <- matrix(muncl, nrow=nrow(Pnrf), ncol=ncol(Pnrf))
                munrf <- t(t(munrf) + bnrf)
                NNrf <- as.matrix(groupSums(Pnrf * FUN(munrf), 2, chVrf$sector_use))
            } else {
                NNcr <- NULL
                NNrf <- NULL
            }

            ## combine NS and NN together (weighted avg in overlap zone) for species with Combo (N+S)
            if (type == "C") {
                # averaging comes here
                NNcr <- NNcr[match(rownames(kgrid), rownames(NNcr)),]
                NNrf <- NNrf[match(rownames(kgrid), rownames(NNrf)),]
                NNcr[is.na(NNcr)] <- 0
                NNrf[is.na(NNrf)] <- 0
                rownames(NNcr) <- rownames(NNrf) <- rownames(kgrid)

                NScr <- NScr[match(rownames(kgrid), rownames(NScr)),]
                NSrf <- NSrf[match(rownames(kgrid), rownames(NSrf)),]
                NScr[is.na(NScr)] <- 0
                NSrf[is.na(NSrf)] <- 0
                rownames(NScr) <- rownames(NSrf) <- rownames(kgrid)

                Ncr <- kgrid$wN * NNcr + (1-kgrid$wN) * NScr
                Nrf <- kgrid$wN * NNrf + (1-kgrid$wN) * NSrf
            }
            if (type == "S") {
                Ncr <- NScr
                Nrf <- NSrf
            }
            if (type == "N") {
                Ncr <- NNcr
                Nrf <- NNrf
            }

            save(Ncr, Nrf,
                file=file.path(ROOT, taxon, paste0(spp, ".RData")))
            cat("\t", proc.time()[3]-t0[3], "\n")

        } # end SEFF

    }

}


## check correlations for use-avail

x1 <- c("WhiteSpruceR", "WhiteSpruce1", "WhiteSpruce2", "WhiteSpruce3",
"WhiteSpruce4", "WhiteSpruce5", "WhiteSpruce6", "WhiteSpruce7",
"WhiteSpruce8", "PineR", "Pine1", "Pine2", "Pine3", "Pine4",
"Pine5", "Pine6", "Pine7", "Pine8", "DeciduousR", "Deciduous1",
"Deciduous2", "Deciduous3", "Deciduous4", "Deciduous5", "Deciduous6",
"Deciduous7", "Deciduous8", "MixedwoodR", "Mixedwood1", "Mixedwood2",
"Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7",
"Mixedwood8", "TreedBogR", "TreedBog1", "TreedBog2", "TreedBog3",
"TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7", "TreedBog8",
"TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3", "TreedFen4",
"TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8", "GrassHerb",
"Shrub", "GraminoidFen", "Marsh", "RoughP", "TameP", "Crop",
"Industrial", "Rural", "Urban", "CCWhiteSpruceR", "CCWhiteSpruce1",
"CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "CCPineR",
"CCPine1", "CCPine2", "CCPine3", "CCPine4", "CCDeciduousR", "CCDeciduous1",
"CCDeciduous2", "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR",
"CCMixedwood1", "CCMixedwood2", "CCMixedwood3", "CCMixedwood4",
"ShrubbySwamp", "ShrubbyBog", "ShrubbyFen", "TreedSwamp", "Wellsites",
"EnSeismic", "EnSoftLin", "TrSoftLin", "HardLin")
x2 <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak",
"Other", "RoughP", "TameP", "Crop", "Urban", "Rural", "Industrial",
"EnSeismic", "EnSoftLin", "TrSoftLin", "HardLin", "Wellsites")

CFS <- CFN <- NULL
for (taxon in names(COEFS)) {
        if (taxon == "birds") {
            cfn <- COEFS[[taxon]]$north$marginal[,x1,1]
            cfs <- COEFS[[taxon]]$south$marginal[,x2,1]
            FUN <- function (eta) {
                lam <- pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)
                1-exp(-lam*7*2)
            }
        } else {
            cfn <- COEFS[[taxon]]$north[,x1,1]
            cfs <- COEFS[[taxon]]$south[,x2,1]
            FUN <- binomial()$linkinv
        }
    CFS <- rbind(CFS, FUN(cfs))
    CFN <- rbind(CFN, FUN(cfn))
}

max(CFN)
max(CFS)
u <- colnames(CFN)
u <- gsub("CC", "", u)
for (i in c("WhiteSpruce", "Pine", "Mixedwood", "Deciduous", "TreedBog", "TreedFen"))
    u[grep(i, u)] <- i
CFN2 <- mefa4::groupMeans(CFN, 2, u)

ds <- as.dist(1-abs(cor(CFS)))
dn <- as.dist(1-abs(cor(CFN2)))

plot(hclust(ds, "ward.D"))
plot(hclust(dn, "ward.D"))


## create sector effects plots

library(mefa4)
library(raster)

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
COEFS <- c(COEFS, COEFS2)
load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]

kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

rn <- rownames(kgrid)[!kgrid$useS & kgrid$NRNAME != "Rocky Mountain"]
rs <- rownames(kgrid)[!kgrid$useN & kgrid$NRNAME != "Rocky Mountain"]

AreaN <- colSums(groupSums(trVeg[rn,], 2, chVeg$sector))
AreaN <- 100 * AreaN/sum(AreaN)

AreaS <- colSums(groupSums(trVeg[rs,], 2, chVeg$sector))
AreaS <- 100 * AreaS/sum(AreaS)

ROOT <- "s:/AB_data_v2020/Results/pred"
ROOT2 <- "s:/AB_data_v2020/Results/web"
sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")

DEL_FOR <- FALSE
if (DEL_FOR)
    AreaS["Forestry"] <- 0

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt) {
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

Rmaskn <- make_raster(as.integer(1-kgrid$useS), kgrid, rt)
values(Rmaskn)[values(Rmaskn) == 0] <- NA
Rmasks <- make_raster(as.integer(1-kgrid$useN), kgrid, rt)
values(Rmasks)[values(Rmasks) == 0] <- NA
Rmaskm <- make_raster(as.integer(!kgrid$NRNAME == "Rocky Mountain"), kgrid, rt)
values(Rmaskm)[values(Rmaskm) == 0] <- NA
Rw <- make_raster(as.integer(kgrid$pWater > 0.99), kgrid, rt)
values(Rw)[values(Rw) == 0] <- NA

col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
col2 <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B",
    "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)
col3 <- colorRampPalette(c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221"))(200)
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude



if (FALSE) {
spp <- "AlderFlycatcher"
taxon <- "birds"

spp <- "Actaea.rubra"
taxon <- "vplants"

load(file.path(ROOT, taxon, paste0(spp, ".RData"))) # Ncr, Nrf

#cure4insect:::.plot_sector1
#function(Curr, Ref, Area, RefTotal, main, col, ylim, ylab, xlab, subset=NULL, ...)

RefN <- colSums(Nrf[rn,])
RefTotalN <- sum(RefN)
CurrN <- colSums(Ncr[rn,])

## delete forestry here???
RefS <- colSums(Nrf[rs,])
if (DEL_FOR)
    RefS["Forestry"] <- 0
RefTotalS <- sum(RefS)
CurrS <- colSums(Ncr[rs,])
if (DEL_FOR)
    CurrS["Forestry"] <- 0

SEn <- cure4insect:::.plot_sector1(CurrN[sectors], RefN[sectors], AreaN, RefTotalN, main=spp)
SEs <- cure4insect:::.plot_sector1(CurrS[sectors], RefS[sectors], AreaS, RefTotalS, main=spp)
}

SE <- list()
TAXA <- names(COEFS)
doSEFF <- TRUE
doMAPS <- FALSE
for (taxon in TAXA) {

    if (!dir.exists(file.path(ROOT, taxon)))
        dir.create(file.path(ROOT, taxon))

    SE[[taxon]] <- list()
    A <- if (taxon == "birds")
        COEFS[[taxon]]$north$marginal else COEFS[[taxon]]$north
    SPPn <- dimnames(A)[[1]]
    A <- if (taxon == "birds")
        COEFS[[taxon]]$south$marginal else COEFS[[taxon]]$south
    SPPs <- dimnames(A)[[1]]
    SPP <- sort(unique(c(SPPn, SPPs)))
    SPP <- SPP[SPP != "Do.not.analyze"]

    for (spp in SPP) {

        spp0 <- spp

        if (endsWith(spp, "."))
            spp <- substr(spp, 1, nchar(spp)-1)
        if (grepl("\\.\\.", spp))
            spp <- gsub("\\.\\.", "\\.", spp)

        cat(taxon, spp, "\n")
        flush.console()

        TYPE <- "C"
        if (spp %in% SPPs && !(spp %in% SPPn))
            TYPE <- "S"
        if (spp %in% SPPn && !(spp %in% SPPs))
            TYPE <- "N"

        if (!dir.exists(file.path(ROOT2, taxon, spp)))
            dir.create(file.path(ROOT2, taxon, spp))

        load(file.path(ROOT, taxon, paste0(spp0, ".RData"))) # Ncr, Nrf

        if (doSEFF) {
            if (spp %in% SPPn) {
                RefN <- colSums(Nrf[rn,])
                RefTotalN <- sum(RefN)
                CurrN <- colSums(Ncr[rn,])
                png(file.path(ROOT2, taxon, spp, "sector-north.png"), width=3*600, height=500)
                op <- par(mfrow=c(1,3))
                SEn1 <- cure4insect:::.plot_sector1(
                    Curr=CurrN[sectors], Ref=RefN[sectors],
                    Area=AreaN, RefTotal=RefTotalN, main=paste(spp, "Unit effect"))
                SEn2 <- cure4insect:::.plot_sector2(
                    Curr=CurrN[sectors], Ref=RefN[sectors],
                    RefTotal=RefTotalN, regional=TRUE, main="Regional effect")
                SEn3 <- cure4insect:::.plot_sector2(
                    Curr=CurrN[sectors], Ref=RefN[sectors],
                    RefTotal=RefTotalN, regional=FALSE, main="Under HF effect")
                par(op)
                dev.off()
                SEn <- list(unit=SEn1, regional=SEn2, underhf=SEn3)
            } else {
                SEn <- NULL
            }

            if (spp %in% SPPs) {
                ## delete forestry here???
                RefS <- colSums(Nrf[rs,])
                if (DEL_FOR)
                    RefS["Forestry"] <- 0
                RefTotalS <- sum(RefS)
                CurrS <- colSums(Ncr[rs,])
                if (DEL_FOR)
                    CurrS["Forestry"] <- 0
                png(file.path(ROOT2, taxon, spp, "sector-south.png"), width=3*600, height=500)
                op <- par(mfrow=c(1,3))
                SEs1 <- cure4insect:::.plot_sector1(
                    Curr=CurrS[sectors], Ref=RefS[sectors],
                    Area=AreaS, RefTotal=RefTotalS, main=paste(spp, "Unit effect"))
                SEs2 <- cure4insect:::.plot_sector2(
                    Curr=CurrS[sectors], Ref=RefS[sectors],
                    RefTotal=RefTotalS, regional=TRUE, main="Regional effect")
                SEs3 <- cure4insect:::.plot_sector2(
                    Curr=CurrS[sectors], Ref=RefS[sectors],
                    RefTotal=RefTotalS, regional=FALSE, main="Under HF effect")
                par(op)
                dev.off()
                SEs <- list(unit=SEs1, regional=SEs2, underhf=SEs3)
            } else {
                SEs <- NULL
            }

            SE[[taxon]][[spp]] <- list(north=SEn, south=SEs)
        }

        if (doMAPS) {
            Dcr <- rowSums(Ncr)[match(rownames(kgrid), rownames(Ncr))]
            Drf <- rowSums(Nrf)[match(rownames(kgrid), rownames(Nrf))]
            NA_VAL <- 0
            if (spp == "pH")
                NA_VAL <- 7

            Dcr[is.na(Dcr)] <- NA_VAL
            q <- quantile(Dcr, 0.99)
            Dcr[Dcr > q] <- q
            Drf[is.na(Drf)] <- NA_VAL
            q <- quantile(Drf, 0.99)
            Drf[Drf > q] <- q
            MAX <- max(Dcr, Drf)
            if (taxon == "habitats" && COEFS[[taxon]]$species[spp, "Unit"] == "%" && spp != "SoilCarbon")
                MAX <- 1
            if (spp == "pH") {
                Dcr[Dcr < 0] <- 0
                Drf[Drf < 0] <- 0
                MAX <- 14
            }

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
    #        if (taxon != "birds") {
    #            Rcr <- mask(Rcr, Rmaskm)
    #            Rrf <- mask(Rrf, Rmaskm)
    #            Rdf <- mask(Rdf, Rmaskm)
    #        }

            png(file.path(ROOT2, taxon, spp, "map.png"),
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

            if (!dir.exists(file.path("s:/AB_data_v2020/Results/normalized-maps", taxon)))
                dir.create(file.path("s:/AB_data_v2020/Results/normalized-maps", taxon))

            out <- data.frame(Current=cr, Reference=rf, Difference=df)
            rownames(out) <- rownames(kgrid)
            out <- out[rownames(Ncr),]
            save(out, file=file.path("s:/AB_data_v2020/Results/normalized-maps", taxon,
                paste0(spp, ".RData")))

            gc()
        }


    }

}

save(SE,  file="s:/AB_data_v2020/Results/SE-ESTIMATES.RData")




