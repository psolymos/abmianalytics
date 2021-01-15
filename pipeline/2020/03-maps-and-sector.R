## mapping and sector effects
library(mefa4)
library(qs)
library(raster)

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
load("s:/AB_data_v2020/Results/COEFS3-mammals.RData")

load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide_water.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]

if (FALSE) {
xvr <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "rf"])
xvc <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "cr"])
setdiff(colnames(xvr), colnames(xvc))
setdiff(colnames(xvc), colnames(xvr))
i <- intersect(colnames(xvc), colnames(xvr))
d <- (xvr[,i] - xvc[,i]) / 10^6
range(d)
}

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

## mammals south:
# "Lat3"       "Lat2Long2"  "LongMAT"    "PeaceRiver"
## mammals north:
# [1] "Lat3"                 "Lat2Long2"            "LongMAT"              "NSR1Parkland"
# [5] "NSR1DryMixedwood"     "NSR1CentralMixedwood" "NSR1Foothills"        "NSR1North"
# [9] "NSR1Shield"           "NSR1Mountain"
kgrid$NSR1<-c(rep("Parkland",3),"DryMixedwood","CentralMixedwood",rep("Foothills",2),
    rep("North",4),rep("Shield",3),rep("Mountain",3))[match(kgrid$NSRNAME,
        c("Central Parkland","Foothills Parkland","Peace River Parkland","Dry Mixedwood",
            "Central Mixedwood","Lower Foothills","Upper Foothills",
            "Lower Boreal Highlands","Upper Boreal Highlands","Boreal Subarctic",
            "Northern Mixedwood","Athabasca Plain","Kazan Uplands","Peace-Athabasca Delta",
            "Montane","Subalpine","Alpine"))]
kgrid$NSR1[is.na(kgrid$NSR1)] <- ""
kgrid$NSR1 <- as.factor(kgrid$NSR1)
kgrid$NSR1 <- relevel(kgrid$NSR1, "")
Xclim_mammal <- cbind(Xclim_nonb,
    Lat3=Xclim_nonb[,"Lat"]^3,
    Lat2Long2=Xclim_nonb[,"Lat"]^2 * Xclim_nonb[,"Long"]^2,
    LongMAT=Xclim_nonb[,"Long"] * Xclim_nonb[,"MAT"],
    PeaceRiver=ifelse(kgrid$NSRNAME=="Peace River Parkland", 1, 0),
    model.matrix(~NSR1,kgrid)[,-1])
Xclim_mammal[,"MAT2"] <- Xclim_mammal[,"MAT"]*(Xclim_mammal[,"MAT"]+10)


kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"

UseN <- rownames(kgrid)[!kgrid$useS]
UseS <- rownames(kgrid)[!kgrid$useN]

Xclim_bird_S <- Xclim_bird[UseS,]
Xclim_nonb_S <- Xclim_nonb[UseS,]
Xclim_mamm_S <- Xclim_mammal[UseS,]
pA <- kgrid[UseS, "pAspen"]
Xclim_bird_N <- Xclim_bird[UseN,]
Xclim_nonb_N <- Xclim_nonb[UseN,]
Xclim_mamm_N <- Xclim_mammal[UseN,]

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
kgrid$pSoilUnk <- rowSums(trSoil[,UNKN,drop=FALSE]) / rowSums(trSoil)
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
#hist(kgrid$wN[!kgrid$useS & !kgrid$useN]) # overlap


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
        "RuralUrban", "RuralUrban", "Energy", "RuralUrban", "Energy",
        "Misc", "Energy", "Energy", "Energy", "Energy", "Energy", "Transportation",
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
        "TameP", "Industrial", "CCDeciduousR", "CCDeciduous1", "CCDeciduous2",
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
        "RuralUrban", "Energy", "Misc", "Energy", "Energy", "Energy",
        "Energy", "Energy", "Transportation", "Transportation", "Transportation",
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
#lt$south["OtherDisturbedVegetation", "Sector"] <- "RuralUrban"
#lt$north["OtherDisturbedVegetation", "Sector"] <- "RuralUrban"
#lt$north["HighDensityLivestockOperation", "Label"] <- "Industrial"

## define sectors to be used here
## i.e. change to attribution stuff or keep finer level transitions etc
chVeg$sector_use <- chVeg$sector
chSoil$sector_use <- chSoil$sector

## fine sector attribution
if (FALSE) {
    tmp1 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v2020.csv")
    tmp2 <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
    chVeg$sector_fine <- tmp1$SectorFine[match(chVeg$cr, tmp1$ID)]
    chSoil$sector_fine <- tmp2$SectorFine[match(chSoil$cr, tmp2$ID)]
    chVeg$sector_use <- chVeg$sector_fine
    chSoil$sector_use <- chSoil$sector_fine
}

## process south monster matrix and find some efficiencies

compare_sets(chSoil$cr, rownames(lt$south))
compare_sets(chSoil$rf, rownames(lt$south))
chSoil$cr2 <- lt$south$Label[match(chSoil$cr, rownames(lt$south))]
chSoil$rf2 <- lt$south$Label[match(chSoil$rf, rownames(lt$south))]
#chSoil$sector2 <- lt$south$Sector[match(chSoil$cr, rownames(lt$south))]
#with(chSoil, table(sector, sector2)) # sector definition is up to date

chSoil[chSoil$cr=="WindGenerationFacility",]
rowSums(with(chSoil[chSoil$sector != "Native",], table(cr2, sector_use))>0)
with(chSoil[chSoil$sector != "Native",], table(cr2, sector_use))
chSoil[chSoil$cr2=="Mine",c("sector_use","cr2", "cr")]
chSoil[chSoil$cr2=="Industrial",c("sector_use","cr2", "cr")]

# aggregating the monster matrix
# need to use _sector b/c Industrial and Mine belongs to multiple sectors
chSoil$tr2 <- paste0(chSoil$rf2, "->", chSoil$cr2, "_", chSoil$sector_use)
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

chVeg[chVeg$cr=="WindGenerationFacility",]
rowSums(with(chVeg[chVeg$sector != "Native",], table(cr2, sector_use))>0)
with(chVeg[chVeg$sector != "Native",], table(cr2, sector_use))
chVeg[chVeg$cr2=="Mine",c("sector_use","cr2", "cr")]
chVeg[chVeg$cr2=="Industrial",c("sector_use","cr2", "cr")]

# aggregating the monster matrix
# need to use _sector b/c Industrial and Mine belongs to multiple sectors
chVeg$tr2 <- paste0(chVeg$rf2, "->", chVeg$cr2, "_", chVeg$sector_use)
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

## some checks
if (FALSE) {
xvr <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "rf2"])
xvc <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "cr2"])
setdiff(colnames(xvr), colnames(xvc))
setdiff(colnames(xvc), colnames(xvr))
i <- intersect(colnames(xvc), colnames(xvr))
d <- (xvr[,i] - xvc[,i]) / 10^6
range(d)
}

rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))
make_raster <- function(value, rc, rt) {
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}
f <- function(x) {
    if (is.null(x)) {
        v <- rep(0, nrow(kgrid))
    } else {
        v <- if (is.null(dim(x)))
            x else rowSums(x)
        q <- quantile(v, 0.99)
        v[v>q] <- q
    }
    make_raster(v, kgrid, rt)
}
# this helps quickly plot stuff on a map
p <- function(vcl, ...) {
    if (is.null(vcl))
        return(NULL)
    vcl <- if (is.null(dim(vcl))) vcl else rowSums(vcl)
    vcl <- vcl[match(rownames(kgrid),names(vcl))]
    print(summary(vcl))
    vcl[is.na(vcl)] <- -999
    vr <- f(vcl)
    values(vr)[values(vr) == -999] <- NA
    plot(vr, col=hcl.colors(100, "plasma"), ...)
    invisible(vr)
}

COEFS$mammals <- COEFS2$mammals
COEFS$habitats <- COEFS2$habitats
COEFS$nnplants <- COEFS2$nnplants

## species specific part begins here

#spp <- "AlderFlycatcher"
#spp <- "RedtailedHawk"
#taxon <- "birds"

#spp <- "Actaea.rubra"
#taxon <- "vplants"

BOOT <- TRUE # do bootstrap?
BMAX <- 100

ROOT <- "s:/AB_data_v2020/Results/pred"
#ROOT <- "s:/AB_data_v2020/Results/pred-fine" # fine sectors
ROOT2 <- "s:/AB_data_v2020/Results/pred-boot"
#ROOT2 <- "s:/AB_data_v2020/Results/pred-boot-fine"


TAXA <- names(COEFS)
#TAXA="birds"
#TAXA="mammals"
#spp="Moose"
for (taxon in TAXA) {

    if (!dir.exists(file.path(ROOT, taxon)))
        dir.create(file.path(ROOT, taxon))
    if (BOOT && !dir.exists(file.path(ROOT2, taxon)))
        dir.create(file.path(ROOT2, taxon))

    SPPn <- if (taxon!="birds")
        dimnames(COEFS[[taxon]]$north)[[1]] else dimnames(COEFS[[taxon]]$north$joint)[[1]]
    SPPs <- if (taxon!="birds")
        dimnames(COEFS[[taxon]]$south)[[1]] else dimnames(COEFS[[taxon]]$south$joint)[[1]]
    if (taxon == "mammals") {
        TAB <- COEFS2$mammals$species
        rownames(TAB) <- TAB$SpeciesID
        SPPn <- rownames(TAB)[TAB$ModelNorth]
        SPPs <- rownames(TAB)[TAB$ModelSouth]
#        SPPnx <- rownames(COEFS3$north$total)
#        SPPsx <- rownames(COEFS3$south$total)
    }
    SPP <- sort(union(SPPn, SPPs))
    #SPPn <- intersect(SPPn, SPPs)
    #SPP <- c("RuffedGrouse", "PileatedWoodpecker", "BrownheadedCowbird")#, "Ovenbird", "CanadaWarbler")

    for (spp in SPP) {

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
        if (taxon == "mammals") {
            # this is a robust exponential
            #FUN <- function (eta)
            #    pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)
            INFO <- TAB[spp,]
            II <- read.csv("s:/AB_data_v2020/Results/Camera mammal models revised June 2020/Mammal header table June 2020.csv")
            rownames(II) <- II$SpeciesID
            II <- II[spp,]
            Link <- list(
                N=as.character(II$LinkSpclimNorth),
                S=as.character(II$LinkSpclimSouth))
            if (is.na(Link$N))
                Link$N <- "Log" # option 3 - no climate model, use 0
            if (is.na(Link$S))
                Link$S <- "Log" # option 3 - no climate model, use 0
            # NOTE!!!!
            # mammal habitat coefs are NOT transfoed to the log/logit scale
            ## (other taxa have cfn and cfs as log(x)/pogit(x))
            cfn_tot <- if (type == "S") NULL else COEFS3$north$total[spp,]
            cfs_tot <- if (type == "N") NULL else COEFS3$south$total[spp,]
            cfn_pa <-  if (type == "S") NULL else COEFS3$north$pa[spp,]
            cfs_pa <-  if (type == "N") NULL else COEFS3$south$pa[spp,]
            cfn_agp <- if (type == "S") NULL else COEFS3$north$agp[spp,]
            cfs_agp <- if (type == "N") NULL else COEFS3$south$agp[spp,]
            if (!is.null(cfn_agp))
                cfn_agp[is.na(cfn_agp)] <- 0
            if (!is.null(cfs_agp))
                cfs_agp[is.na(cfs_agp)] <- 0
            cfn_cl <-  if (type == "S") NULL else COEFS3$north$clim[spp,]
            cfs_cl <-  if (type == "N") NULL else COEFS3$south$clim[spp,]
            XclimS <- Xclim_mamm_S[,names(cfs_cl)]
            XclimN <- Xclim_mamm_N[,names(cfn_cl)]
            cfs_asp <- if (type == "N") NULL else COEFS3$south$asp[spp,]

        }
        if (taxon == "habitats")
            FUN <- binomial(COEFS[[taxon]]$species[spp, "LinkHabitat"])$linkinv

        ## only doing 1 run for sector effects
        imax <- if (BOOT)
            BMAX else 1
        #i <- 1
        ## for loop for boot
        ivec <- seq_len(imax)
        #ivec <- 1:25
        #ivec <- 26:100
        for (i in ivec) {

            cat(taxon, spp, i)
            flush.console()
            t0 <- proc.time()
            if (taxon == "mammals") {
                if (type != "N") {
                    gc()
                    ## south calculations for the i'th run
                    ## space-climate coefs
                    bscl <- cfs_cl[colnames(XclimS)]
                    bscl[is.na(bscl)] <- 0 # this happens for habitat elements
                    muscl <- drop(XclimS %*% bscl)

                    muspaPA <-  pA * cfs_asp["pAspen.pa"]
                    muspaAGP <- pA * cfs_asp["pAspen.agp"]

                    ## habitat
                    if (Link$S == "Log") {
                        # pa
                        tmps <- c(cfs_pa, SnowIce=0)
                        tmps <- qlogis(tmps)
                        tmps[is.na(tmps)] <- -10^4
                        bscrp <- tmps[chScr$cr2] # current land cover
                        bsrfp <- tmps[chSrf$rf2] # reference land cover
                        # agp
                        tmps <- c(cfs_agp, SnowIce=0)
                        tmps <- log(tmps)
                        tmps[is.na(tmps)] <- -10^4
                        bscra <- tmps[chScr$cr2] # current land cover
                        bsrfa <- tmps[chSrf$rf2] # reference land cover

                        mpacr <- matrix(bscrp, nrow=nrow(Pscr), ncol=ncol(Pscr), byrow=TRUE) + muspaPA
                        magpcr <- matrix(bscra, nrow=nrow(Pscr), ncol=ncol(Pscr), byrow=TRUE) + muspaAGP
                        muscr <- matrix(muscl, nrow=nrow(Pscr), ncol=ncol(Pscr))
                        muscr <- exp(log(plogis(mpacr) * exp(magpcr)) + muscr)
                        NScr <- as.matrix(groupSums(Pscr * muscr, 2, chScr$sector_use))

                        mparf <- matrix(bsrfp, nrow=nrow(Psrf), ncol=ncol(Psrf), byrow=TRUE) + muspaPA
                        magprf <- matrix(bsrfa, nrow=nrow(Psrf), ncol=ncol(Psrf), byrow=TRUE) + muspaAGP
                        musrf <- matrix(muscl, nrow=nrow(Psrf), ncol=ncol(Psrf))
                        musrf <- exp(log(plogis(mparf) * exp(magprf)) + musrf)
                        NSrf <- as.matrix(groupSums(Psrf * musrf, 2, chSrf$sector_use))
                    } else {
                        # pa
                        tmpsp <- c(cfs_pa, SnowIce=0)
                        tmpsp <- qlogis(tmpsp)
                        tmpsp[is.na(tmpsp)] <- -10^4
                        bscrp <- tmpsp[chScr$cr2] # current land cover
                        bsrfp <- tmpsp[chSrf$rf2] # reference land cover
                        # agp
                        tmpsa <- c(cfs_agp, SnowIce=0)
                        tmpsa <- log(tmpsa)
                        tmpsa[is.na(tmpsa)] <- -10^4
                        bscra <- tmpsa[chScr$cr2] # current land cover
                        bsrfa <- tmpsa[chSrf$rf2] # reference land cover

                        mpacr <- matrix(bscrp, nrow=nrow(Pscr), ncol=ncol(Pscr), byrow=TRUE) + muspaPA
                        magpcr <- matrix(bscra, nrow=nrow(Pscr), ncol=ncol(Pscr), byrow=TRUE) + muspaAGP
                        #muscr0 <- matrix(muscl, nrow=nrow(Pscr), ncol=ncol(Pscr))
                        muscr <- plogis(mpacr + muscl) * exp(magpcr)
                        NScr <- as.matrix(groupSums(Pscr * muscr, 2, chScr$sector_use))

                        mparf <- matrix(bsrfp, nrow=nrow(Psrf), ncol=ncol(Psrf), byrow=TRUE) + muspaPA
                        magprf <- matrix(bsrfa, nrow=nrow(Psrf), ncol=ncol(Psrf), byrow=TRUE) + muspaAGP
                        #musrf0 <- matrix(muscl, nrow=nrow(Psrf), ncol=ncol(Psrf))
                        musrf <- plogis(mparf + muscl) * exp(magprf)
                        NSrf <- as.matrix(groupSums(Psrf * musrf, 2, chSrf$sector_use))

                    }
                    NScr <- cbind(NScr, Forestry=0)
                    NSrf <- cbind(NSrf, Forestry=0)
                } else {
                    NScr <- NULL
                    NSrf <- NULL
                }
                if (type != "S") {
                    gc()
                    ## north calculations for the i'th run
                    ## space-climate coefs
                    bncl <- cfn_cl[colnames(XclimN)]
                    bncl[is.na(bncl)] <- 0 # this happens for habitat elements
                    muncl <- drop(XclimN %*% bncl)
                    ## habitat
                    if (Link$N == "Log") {
                        # this is total, log transformed
                        tmpn <- c(cfn_tot, SnowIce=0)
                        tmpn <- log(tmpn)
                        tmpn[is.na(tmpn)] <- -10^4
                        bncr <- tmpn[chVcr$cr2] # current land cover
                        bnrf <- tmpn[chVrf$rf2] # reference land cover

                        muncr <- matrix(muncl, nrow=nrow(Pncr), ncol=ncol(Pncr))
                        muncr <- t(t(muncr) + bncr)
                        NNcr <- as.matrix(groupSums(Pncr * exp(muncr), 2, chVcr$sector_use))

                        munrf <- matrix(muncl, nrow=nrow(Pnrf), ncol=ncol(Pnrf))
                        munrf <- t(t(munrf) + bnrf)
                        NNrf <- as.matrix(groupSums(Pnrf * exp(munrf), 2, chVrf$sector_use))
                    } else {
                        # this is PA, qlogis transformed
                        tmpn <- c(cfn_pa, SnowIce=0)
                        tmpn <- qlogis(tmpn)
                        tmpn[is.na(tmpn)] <- -10^4
                        bncr <- tmpn[chVcr$cr2] # current land cover
                        bnrf <- tmpn[chVrf$rf2] # reference land cover
                        # this is the agp, untransformed
                        tmpna <- c(cfn_agp, SnowIce=0)
                        tmpna[is.na(tmpna)] <- 0
                        tmpna[is.infinite(tmpna)] <- max(tmpna[!is.infinite(tmpna)])
                        bncra <- tmpna[chVcr$cr2] # current land cover
                        bnrfa <- tmpna[chVrf$rf2] # reference land cover

                        muncr <- matrix(muncl, nrow=nrow(Pncr), ncol=ncol(Pncr))
                        muncr <- t(plogis(t(muncr) + bncr) * bncra)
                        NNcr <- as.matrix(groupSums(Pncr * muncr, 2, chVcr$sector_use))

                        munrf <- matrix(muncl, nrow=nrow(Pnrf), ncol=ncol(Pnrf))
                        munrf <- t(plogis(t(munrf) + bnrf) * bnrfa)
                        NNrf <- as.matrix(groupSums(Pnrf * munrf, 2, chVrf$sector_use))
                    }
                } else {
                    NNcr <- NULL
                    NNrf <- NULL
                }
            ## non-mammal taxa
            } else {
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
                    muscl <- drop(XclimS %*% bscl)
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

            ## mammals troubleshooting
            if (FALSE) {
                u <- read.csv(paste0(
                    "s:/AB_data_v2020/Results/Camera mammal models South revised Nov 2020/Combine regions/Km2 summaries Oct 2020/",
                    spp, ".csv"))
                NDcr <- u$Curr[match(rownames(kgrid), u$LinkID)]
                NDcr[is.na(NDcr)] <- 0
                Ncr <- Ncr[match(rownames(kgrid), rownames(Ncr)),]
                Ncr[is.na(Ncr)] <- 0
                if (!is.null(NNcr)) {
                    NNcr <- NNcr[match(rownames(kgrid), rownames(NNcr)),]
                    NNcr[is.na(NNcr)] <- 0
                }
                if (!is.null(NScr)) {
                    NScr <- NScr[match(rownames(kgrid), rownames(NScr)),]
                    NScr[is.na(NScr)] <- 0
                }
                l <- stack(list(south=f(NScr), north=f(NNcr), combo=f(Ncr), Dave=f(NDcr)))

                png(paste0("s:/AB_data_v2020/Results/web1/mammals/",spp,".png"))
                plot(l, col=hcl.colors(100)[30:100])
                dev.off()
            }

            if (BOOT) {
                if (!dir.exists(file.path(ROOT2, taxon, spp)))
                    dir.create(file.path(ROOT2, taxon, spp))
                # this is 5x faster and a littla bit smaller (load with qs::qloadm)
                qsavem(Ncr, Nrf,
                    file=file.path(ROOT2, taxon, spp, paste0(spp, "-", i, ".qrda")))
                #save(Ncr, Nrf,
                #    file=file.path(ROOT2, taxon, spp, paste0(spp, "-", i, ".RData")))
            } else {
                save(Ncr, Nrf,
                    file=file.path(ROOT, taxon, paste0(spp, ".RData")))
            }
            cat("\t", proc.time()[3]-t0[3], "\n")

        } # end of for loop for boot

    }

}


library(raster)
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))
make_raster <- function(value, rc, rt) {
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}
f <- function(x) {
    if (is.null(x)) {
        v <- rep(0, nrow(kgrid))
    } else {
        v <- if (is.null(dim(x)))
            x else rowSums(x)
        q <- quantile(v, 0.99)
        v[v>q] <- q
    }
    make_raster(v, kgrid, rt)
}
u <- read.csv(paste0(
    "s:/AB_data_v2020/Results/Camera mammal models South revised Nov 2020/Combine regions/Km2 summaries Oct 2020/",
    spp, ".csv"))
NDcr <- u$Curr[match(rownames(kgrid), u$LinkID)]
NDcr[is.na(NDcr)] <- 0
l <- stack(list(south=f(NScr), north=f(NNcr), combo=f(Ncr), Dave=f(NDcr)))

png(paste0("s:/AB_data_v2020/Results/web1/mammals/",spp,".png"))
plot(l, col=hcl.colors(100)[30:100])
dev.off()


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
load("s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
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
#ROOT <- "s:/AB_data_v2020/Results/pred1"
#ROOT2 <- "s:/AB_data_v2020/Results/web1"
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
doMAPS <- TRUE
#TAXA="mammals"
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
    if (taxon == "mammals") {
        TAB <- COEFS2$mammals$species
        rownames(TAB) <- TAB$SpeciesID
        SPPn <- rownames(TAB)[TAB$ModelNorth]
        SPPs <- rownames(TAB)[TAB$ModelSouth]
    }
    #SPPn <- intersect(SPPn, SPPs)
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

#            SE[[taxon]][[spp]] <- list(north=SEn, south=SEs)
            SE[[spp]] <- list(north=SEn, south=SEs)
        }

        if (doMAPS) {
            if (taxon=="mammals") {
                u <- read.csv(paste0(
                    "s:/AB_data_v2020/Results/Camera mammal models South revised Nov 2020/Combine regions/Km2 summaries Oct 2020/",
                    spp, ".csv"))
                NDcr <- u$Curr[match(rownames(kgrid), u$LinkID)]
                NDrf <- u$Ref[match(rownames(kgrid), u$LinkID)]
                names(NDcr) <- names(NDrf) <- rownames(kgrid)
#                Ncr <- Ncr[match(u$LinkID, rownames(Ncr)),]
#                Nrf <- Nrf[match(u$LinkID, rownames(Nrf)),]



            }

            Dcr <- rowSums(Ncr)[match(rownames(kgrid), rownames(Ncr))]
            Drf <- rowSums(Nrf)[match(rownames(kgrid), rownames(Nrf))]

            if (FALSE) {
                summary(NDcr)
                summary(Dcr)
                summary(NDrf)
                summary(Drf)
                cor(cbind(NDcr,Dcr,NDrf,Drf), use="c")
                library(ggplot2)
                ggplot(data.frame(NDcr,Dcr,NDrf,Drf), aes(x=Dcr, y=NDcr)) + geom_bin2d()
                ggplot(data.frame(NDcr,Dcr,NDrf,Drf), aes(x=Drf, y=NDrf)) + geom_bin2d()
                plot(NDcr,Dcr)
                plot(NDrf,Drf)

                op <- par(mfrow=c(2,2), mar=c(1,4,1,4))
                plot(make_raster(Dcr, kgrid, rt),col=hcl.colors(50))
                plot(make_raster(Drf, kgrid, rt),col=hcl.colors(50))
                plot(make_raster(NDcr, kgrid, rt),col=hcl.colors(50))
                plot(make_raster(NDrf, kgrid, rt),col=hcl.colors(50))
                par(op)
            }


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
            if (taxon == "habitats" && COEFS[[taxon]]$species[spp, "Comments"] == "%" && spp != "SoilCarbon")
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
    save(SE,
        file=paste0("s:/AB_data_v2020/Results/SEffect-", taxon, ".RData"))

}


