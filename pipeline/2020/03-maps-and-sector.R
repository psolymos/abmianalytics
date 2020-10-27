## mapping and sector effects
library(mefa4)

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]

## common functions and variables

## organize space-climate variables using kgrid table
make_clim <- function(x, birds=FALSE) {
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
        z <- with(x, cbind(
            Intercept=1,
            Lat=x$POINT_Y,
            Long=x$POINT_X,
            AHM=x$AHM,
            PET=x$PET,
            FFP=x$FFP,
            MAP=x$MAP,
            MAT=x$MAT,
            MCMT=x$MCMT,
            MWMT=x$MWMT,
            Lat2=x$POINT_Y^2,
            Long2=x$POINT_X^2,
            LatLong=x$POINT_X*x$POINT_Y,
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

TAXA <- names(COEFS)
for (taxon in TAXA) {

    if (!dir.exists(file.path(ROOT, taxon)))
        dir.create(file.path(ROOT, taxon))
    if (!dir.exists(file.path(ROOT2, taxon)))
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




## create sector effects plots

library(mefa4)
library(raster)

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
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
for (taxon in names(COEFS)) {
    if (!dir.exists(file.path(ROOT2, taxon)))
        dir.create(file.path(ROOT2, taxon))
    SE[[taxon]] <- list()
    A <- if (taxon == "birds")
        COEFS[[taxon]]$north$marginal else COEFS[[taxon]]$north
    SPPn <- dimnames(A)[[1]]
    A <- if (taxon == "birds")
        COEFS[[taxon]]$south$marginal else COEFS[[taxon]]$south
    SPPs <- dimnames(A)[[1]]
    SPP <- sort(unique(c(SPPn, SPPs)))

    for (spp in SPP) {
        cat(taxon, spp, "\n")
        flush.console()

        TYPE <- "C"
        if (spp %in% SPPs && !(spp %in% SPPn))
            TYPE <- "S"
        if (spp %in% SPPn && !(spp %in% SPPs))
            TYPE <- "N"

        if (!dir.exists(file.path(ROOT2, taxon, spp)))
            dir.create(file.path(ROOT2, taxon, spp))

        load(file.path(ROOT, taxon, paste0(spp, ".RData"))) # Ncr, Nrf

        if (spp %in% SPPn) {
            RefN <- colSums(Nrf[rn,])
            RefTotalN <- sum(RefN)
            CurrN <- colSums(Ncr[rn,])
            png(file.path(ROOT2, taxon, spp, "sector-north.png"), width=600, height=500)
            SEn <- cure4insect:::.plot_sector1(
                CurrN[sectors], RefN[sectors],
                AreaN, RefTotalN, main=spp)
            dev.off()
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
            png(file.path(ROOT2, taxon, spp, "sector-south.png"), width=600, height=500)
            SEs <- cure4insect:::.plot_sector1(
                CurrS[sectors], RefS[sectors],
                AreaS, RefTotalS, main=spp)
            dev.off()
        } else {
            SEs <- NULL
        }

        SE[[taxon]][[spp]] <- list(north=SEn, south=SEs)


        Dcr <- rowSums(Ncr)[match(rownames(kgrid), rownames(Ncr))]
        q <- quantile(Dcr, 0.99)
        Dcr[Dcr > q] <- q
        Drf <- rowSums(Nrf)[match(rownames(kgrid), rownames(Nrf))]
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
        if (taxon != "birds") {
            Rcr <- mask(Rcr, Rmaskm)
            Rrf <- mask(Rrf, Rmaskm)
            Rdf <- mask(Rdf, Rmaskm)
        }

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

save(SE,  file="s:/AB_data_v2020/Results/SE-ESTIMATES.RData")

## maps




spp <- "AlderFlycatcher"
taxon <- "birds"

spp <- "Actaea.rubra"
taxon <- "vplants"

load(file.path(ROOT, taxon, paste0(spp, ".RData"))) # Ncr, Nrf






# --


library(cure4insect)
library(mefa4)
#set_options(path = "s:/reports")
set_options(path = "d:/abmi/reports")
load_common_data()

load("d:/abmi/sppweb2018/c4i/tables/sector-effects.RData")
load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")
ROOT <- "d:/abmi/AB_data_v2018/www"
#ROOT <- "d:/abmi/sppweb2018/www/"

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


for (gr in c("birds","vplants","lichens","mosses","mites")) {
#gr <- "birds"
#gr <- "vplants"
#gr <- "lichens"
#gr <- "mosses"
#gr <- "mites"

#SPP <- rownames(resn[resn$Taxon == gr,])
SPP <- rownames(OUT$Species[(OUT$Species$ModelNorth |
        OUT$Species$ModelSouth) & OUT$Species$Group == gr,])
#spp <- "AlderFlycatcher"

for (spp in SPP) {

    cat(gr, spp, "\n");flush.console()

    TYPE <- "C"
    if (OUT$Species[spp, "ModelNorth"] && !OUT$Species[spp, "ModelSouth"])
        TYPE <- "N"
    if (!OUT$Species[spp, "ModelNorth"] && OUT$Species[spp, "ModelSouth"])
        TYPE <- "S"

    NAM <- OUT$Species[spp, "CommonName"]
    if (is.na(NAM))
        NAM <- OUT$Species[spp, "ScientificName"]
    NAM <- as.character(NAM)

    if (TYPE != "S") {
        ## hab-north
        png(paste0(ROOT, "/figs/", gr, "/", spp, "-coef-north.png"),
            height=800, width=2200, res=150)
        layout(matrix(c(1,1,2), nrow=1))
        plot_abundance(spp, "veg_coef")
        par(mar=c(12,5,4,3))
        plot_abundance(spp, "veg_lin", main="")
        dev.off()

        png(paste0(ROOT, "/figs/", gr, "/", spp, "-sector-north.png"),
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
        png(paste0(ROOT, "/figs/", gr, "/", spp, "-coef-south.png"),
            height=800, width=2000, res=150)
        layout(matrix(c(1,2,3), nrow=1))
        p1 <- plot_abundance(spp, "soil_coef", paspen=0, plot=FALSE)
        p2 <- plot_abundance(spp, "soil_coef", paspen=1, plot=FALSE)
        plot_abundance(spp, "soil_coef", paspen=0, ylim=c(0, max(p1, p2)), main=paste0(NAM, " - non treed"))
        plot_abundance(spp, "soil_coef", paspen=1, ylim=c(0, max(p1, p2)), main="treed")
        par(mar=c(11,5,4,3))
        plot_abundance(spp, "soil_lin", main="")
        dev.off()

        png(paste0(ROOT, "/figs/", gr, "/", spp, "-sector-south.png"),
            height=2*500, width=2*1500, res=150)
        op <- par(mfrow=c(1,3))
        plot_sector(ress[spp,], "unit")
        plot_sector(ress[spp,], "regional", main="")
        plot_sector(ress[spp,], "underhf", main="")
        par(op)
        dev.off()
    }

    y <- load_species_data(spp, boot=TRUE)
    ry <- rasterize_results(y)

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
    if (gr != "birds") {
        Rcr <- mask(Rcr, Rmaskm)
        Rrf <- mask(Rrf, Rmaskm)
        Rdf <- mask(Rdf, Rmaskm)
    }

if (TRUE) {
    writeRaster(Rcr, paste0(ROOT, "/normalized-maps/", gr, "/", spp, "-cr.tif"), overwrite=TRUE)
    writeRaster(Rrf, paste0(ROOT, "/normalized-maps/", gr, "/", spp, "-rf.tif"), overwrite=TRUE)
    writeRaster(Rdf, paste0(ROOT, "/normalized-maps/", gr, "/", spp, "-df.tif"), overwrite=TRUE)
}

    ## add here mask for Rockies if needed
if (TRUE) {
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-map.png"),
        height=1500*2, width=1000*2, res=300)
    op <- par(mfrow=c(2,2), mar=c(2,1,2,3))

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Reference", legend=FALSE)
    plot(Rrf, add=TRUE, col=col1[1:max(rf)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Current", legend=FALSE)
    plot(Rcr, add=TRUE, col=col1[1:max(cr)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Difference", legend=FALSE)
    plot(Rdf, add=TRUE, col=col3[min(df):max(df)])
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    plot(rt, col=CE, axes=FALSE, box=FALSE, main="Std. Error (Current)", legend=FALSE)
    plot(ry[["SE"]]/MAX, add=TRUE, col=rev(col2))
    plot(Rw, add=TRUE, col=CW, legend=FALSE)

    par(op)
    dev.off()

    gc()
}

}

}

## use avail figures

load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")

library(RColorBrewer)


tab <- OUT$Species
uan <- OUT$UseavailNorth
uas <- OUT$UseavailSouth

x <- as.matrix(uan[ ,c("Deciduous","Mixedwood","WhiteSpruce","Pine","BlackSpruce","TreedFen","Open","Wetland","HFor","Crop", "TameP", "RoughP","UrbInd","HardLin","SoftLin")])
HabLabel <- c("Deciduous","Mixedwood","Upland Spruce","Pine","Black Spruce","Treed Fen","Open Upland","Open Wetland","Forestry","Crop", "Tame Pasture", "Rough Pasture","Urban/Industry","Hard Linear","Soft Linear" )
col1<-brewer.pal(8, "Dark2")[c(1,1,1,1, 5,5, 6,7)]
col2<-brewer.pal(12, "Paired")[c(4,7,7,7,12,12,10)]
cols <- c(col1,col2)

#SPP <- rownames(tab)[tab$UseavailNorth]
SPP <- rownames(tab)[tab$UseavailNorth & tab$Group == "birds"]

for (spp in SPP) {
    gr <- tab[spp, "Group"]
    spnam <- if (is.na(tab[spp, "CommonName"])) {
        as.character(tab[spp, "ScientificName"])
    } else {
        paste0(as.character(tab[spp, "CommonName"]), " (", as.character(tab[spp, "ScientificName"]), ")")
    }
    cat(gr, spp, "\n");flush.console()
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-useavail-north.png"),
        height=480, width=600)
    op <- par(mar=c(6,4,2,2)+0.1, las=2)
    x1 <- barplot(as.vector(x [spp, ]), horiz=FALSE, ylab="Affinity",space=NULL, col=cols, border=cols, ylim=c(-1,1), axes=FALSE,axisnames=F )
    axis(side=2)
    abline(h=0, col="red4", lwd=2)
    mtext(side=3,at=x1[1],adj=0, spnam, cex=1.2,col="grey40",las=1)
    text(x=x1, y=par()$usr[3]-0.01,labels=HabLabel, srt=60, adj=1, col=cols, xpd=TRUE)
    par(op)
    dev.off()
}

x<- as.matrix(uas[ , c("Productive","Clay","Saline","RapidDrain","Crop","TameP","RoughP","UrbInd","HardLin","SoftLin")])
HabLabel <- c("Productive","Clay","Saline","Rapid Drain","Crop", "Tame Pasture", "Rough Pasture","Urban/Industry","Hard Linear","Soft Linear" )
col1<-brewer.pal(8, "Dark2")[c(7,7,7,7)]
col2<-brewer.pal(12, "Paired")[c(7,7,7,12,12,10)]
cols <- c(col1,col2)

#SPP <- rownames(tab)[tab$UseavailSouth]
SPP <- rownames(tab)[tab$UseavailSouth & tab$Group == "birds"]

for (spp in SPP) {
    gr <- tab[spp, "Group"]
    spnam <- if (is.na(tab[spp, "CommonName"])) {
        as.character(tab[spp, "ScientificName"])
    } else {
        paste0(as.character(tab[spp, "CommonName"]), " (", as.character(tab[spp, "ScientificName"]), ")")
    }
    cat(gr, spp, "\n");flush.console()
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-useavail-south.png"),
        height=480, width=600)
    op <- par(mar=c(6,4,2,2)+0.1, las=2)
    x1 <- barplot(as.vector(x[spp, ]), horiz=FALSE, ylab="Affinity",space=NULL, col=cols, border=cols, ylim=c(-1,1), axes=FALSE,axisnames=F )
    axis(side=2)
    abline(h=0, col="red4", lwd=2)
    mtext(side=3,at=x1[1],adj=0,spnam,cex=1.2,col="grey40",las=1)
    text(x=x1, y=par()$usr[3]-0.01,labels=HabLabel, srt=60, adj=1, col=cols, xpd=TRUE)
    par(op)
    dev.off()
}

## detection maps for non birds

m <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

rnr <- make_raster(as.integer(kgrid$NRNAME), kgrid, rt)
cnr <- c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae')
cnr <- cnr[c(5,6,1,2,4,3)]


ex <- new.env()
gr <- "mites"
load("s:/Result from Ermias_2018/mites/Species detection Mites 2018.RData", envir=ex)

ex <- new.env()
gr <- "lichens"
load("s:/Result from Ermias_2018/lichens/Species detection Lichens 2018.RData", envir=ex)

ex <- new.env()
gr <- "mosses"
load("s:/Result from Ermias_2018/mosses/Species detection Moss 2018.RData", envir=ex)

ex <- new.env()
gr <- "vplants"
load("s:/Result from Ermias_2018/vplants/Species detection Vascular plants 2018.RData", envir=ex)



site <- ex$dd$Site
og <- ex$dd$OnOffGrid == "OG"
ogs <- sapply(strsplit(site[og], "-"), "[[", 3)
site[og] <- ogs
site <- gsub("B", "", site)
siten <- as.integer(site)

xy <- data.frame(x=m$PUBLIC_LONGITUDE, y=m$PUBLIC_LATTITUDE)[match(site, m$SITE_ID),]
coordinates(xy) <- ~x+y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))

yy <- ex$dd[,6:ncol(ex$dd)]
compare_sets(colnames(yy), rownames(tab[tab$Group == gr,]))

SPP <- rownames(tab[tab$Group == gr,])

for (spp in SPP) {
    spnam <- if (is.na(tab[spp, "CommonName"])) {
        as.character(tab[spp, "ScientificName"])
    } else {
        paste0(as.character(tab[spp, "CommonName"]), " (", as.character(tab[spp, "ScientificName"]), ")")
    }
    cat(gr, spp, "\n");flush.console()
    xy0 <- xy[yy[,spp] == 0 & !duplicated(site),]
    xy1 <- xy[yy[,spp] > 0,]
    png(paste0(ROOT, "/figs/", gr, "/", spp, "-det.png"),
        height=1500*1.5, width=1000*1.5, res=300)
    op <- par(mar=c(1,1,1,1))
    plot(rnr,col=cnr, axes=FALSE, box=FALSE, main=spnam, legend=FALSE)
    plot(xy0, add=TRUE, pch=19, col="#aaaaaa88", legend=FALSE, cex=0.8)
    plot(xy1, add=TRUE, pch=19, col="red4", legend=FALSE, cex=0.8)
    par(op)
    dev.off()

}

## detection maps for birds

## detections
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)

ddd <- ee$dd
## subset based on analysis data -- same filtering applied
ddd <- ddd[unique(c(rownames(en$DAT), rownames(es$DAT))),]
ddd <- nonDuplicated(ddd, ddd$SS, TRUE)
yyy <- groupSums(ee$yy, 1, ee$dd$SS)[rownames(ddd),]
yyy[yyy > 0] <- 1
ss <- !is.na(ddd$X) & !is.na(ddd$NRNAME)
ddd <- ddd[ss,]
yyy <- yyy[ss,]

xy <- SpatialPoints(as.matrix(ddd[,c("X","Y")]))
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))
rt10 <- aggregate(rt, fact=10)
sam0 <- rasterize(xy, rt10, field=1, fun='sum')
values(sam0)[!is.na(values(sam0))] <- 1

load("d:/abmi/sppweb2018/c4i/tables/lookup-birds.RData")
tax <- droplevels(Lookup[Lookup$UseavailNorth | Lookup$UseavailSouth,])
rownames(tax) <- tax$Code

SPP <- rownames(OUT$Species[OUT$Species$Group == gr,])
tax <- tax[tax$SpeciesID %in% SPP,]

gr <- "birds"
for (spp in rownames(tax)) {
    cat(gr, spp, "\n");flush.console()
    xy1 <- SpatialPoints(as.matrix(ddd[yyy[,spp] > 0,c("X","Y")]))
    proj4string(xy1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    xy1 <- spTransform(xy1, proj4string(rt))
    sam1 <- rasterize(xy1, rt10, field=1, fun='last')
    png(paste0("d:/abmi/AB_data_v2018/www", "/figs/", gr, "/", as.character(tax[spp, "SpeciesID"]), "-det.png"),
        height=1500*1.5, width=1000*1.5, res=300)
    op <- par(mar=c(1,1,1,1))
    plot(rnr,col=cnr, axes=FALSE, box=FALSE, main=as.character(tax[spp, "CommonName"]), legend=FALSE)
    plot(sam0,add=TRUE, col="#aaaaaa88", legend=FALSE)
    plot(sam1,add=TRUE, col="red4", legend=FALSE)
    par(op)
    dev.off()
}




## mammal camera dump

fl <- list.files("s:/Camera mammals Mar 2019/Figures North/Best model")
x <- gsub("Veg+HF figure best model ", "", fl, fixed=TRUE)
#x[endsWith(x, "Winter.png")] <- gsub("Winter.png", "-winter.png", x[endsWith(x, "Winter.png")])
#x[endsWith(x, "Summer.png")] <- gsub("Summer.png", "-summer.png", x[endsWith(x, "Summer.png")])

i1 <- endsWith(x, "Winter.png")
i2 <- endsWith(x, "Summer.png")

## combined (fine coef, summer+winter)
file.copy(paste0("s:/Camera mammals Mar 2019/Figures North/Best model/", fl[!i1 & !i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-north-combo/", x[!i1 & !i2]))
## honest models, winter
file.copy(paste0("s:/Camera mammals Mar 2019/Figures North/Best model/", fl[i1]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-north-winter/",
    gsub("Winter.png", ".png", x[i1])))
## honest models, summer
file.copy(paste0("s:/Camera mammals Mar 2019/Figures North/Best model/", fl[i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-north-summer/",
    gsub("Summer.png", ".png", x[i2])))

fl1 <- list.files("s:/Camera mammals Mar 2019/Figures South/Best model/Treed")
fl2 <- list.files("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed")
x1 <- gsub("Soil+HF figure best model ", "", fl1, fixed=TRUE)
x2 <- gsub("Soil+HF figure best model ", "", fl2, fixed=TRUE)
all(x1 == x2)

i1 <- endsWith(x1, "Winter.png")
i2 <- endsWith(x1, "Summer.png")

## combined (fine coef, summer+winter)
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Treed/", fl1[!i1 & !i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-combo-treed/", x1[!i1 & !i2]))
## honest models, winter
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Treed/", fl1[i1]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-winter-treed/",
    gsub("Winter.png", ".png", x1[i1])))
## honest models, summer
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Treed/", fl1[i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-summer-treed/",
    gsub("Summer.png", ".png", x1[i2])))

file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed/", fl2[!i1 & !i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-combo-nontreed/", x2[!i1 & !i2]))
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed/", fl2[i1]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-winter-nontreed/",
    gsub("Winter.png", ".png", x2[i1])))
file.copy(paste0("s:/Camera mammals Mar 2019/Figures South/Best model/Non-treed/", fl2[i2]),
    paste0("d:/abmi/reports/2018/images/mammals-camera/coef-south-summer-nontreed/",
    gsub("Summer.png", ".png", x2[i2])))

## jpg to png

library(magick)
fl <- list.files("s:/Camera mammals Mar 2019/Maps/North Climate and spatial/")
x <- gsub(".jpg", ".png", fl, fixed=TRUE)
for (i in 1:length(fl)) {
    img <- image_read(paste0("s:/Camera mammals Mar 2019/Maps/North Climate and spatial/", fl[i]))
    image_write(img, paste0("s:/Camera mammals Mar 2019/Maps/North Climate and spatial/", x[i]), format="png")
}
fl <- list.files("s:/Camera mammals Mar 2019/Maps/South Climate and spatial/")
x <- gsub(".jpg", ".png", fl, fixed=TRUE)
for (i in 1:length(fl)) {
    img <- image_read(paste0("s:/Camera mammals Mar 2019/Maps/South Climate and spatial/", fl[i]))
    image_write(img, paste0("s:/Camera mammals Mar 2019/Maps/South Climate and spatial/", x[i]), format="png")
}


## lookup table

r <- "d:/abmi/reports/2018/images/mammals-camera/"
fl <- list.files("d:/abmi/reports/2018/images/mammals-camera/map-det")

d <- list.dirs(r, full.names=FALSE)[-1]

ls <- lapply(d, function(z) list.files(paste0(r, z)))
names(ls) <- d

spp <- sort(unique(unlist(ls)))

mefa4::compare_sets(fl, spp)
setdiff(fl, spp)
setdiff(spp, fl)
intersect(spp, fl)


SPP <- c(
    "Badger" = "Badger",
    "Beaver" = "Beaver",
    "Bighornsheep" = "Bighorn Sheep",
    "Bison" = "Bison",
    "BlackBear" = "Black Bear",
    "Bobcat" = "Bobcat",
    "CanadaLynx" = "Canada Lynx",
    "ColumbianGroundSquirrel" = "Columbian Ground Squirrel",
    "Cougar" = "Cougar",
    "Coyote" = "Coyote",
    "Deer" = "Deer",
    "Elk" = "Elk",
    "Fisher" = "Fisher",
    "Foxes" = "Foxes",
    "GoldenMantledGroundSquirrel" = "Golden Mantled Ground Squirrel",
    "GrayWolf" = "Gray Wolf",
    "Grizzlybear" = "Grizzly Bear",
    "Groundhog" = "Groundhog",
    "HoaryMarmot" = "Hoary Marmot",
    "LeastChipmunk" = "Least Chipmunk",
    "Marten" = "Marten",
    "Mink" = "Mink",
    "Moose" = "Moose",
    "Mountaingoat" = "Mountain Goat",
    "Muledeer" = "Muledeer",
    "Muskrat" = "Muskrat",
    "NorthernFlyingSquirrel" = "Northern Flying Squirrel",
    "Porcupine" = "Porcupine",
    "Pronghorn" = "Pronghorn",
    "Raccoon" = "Raccoon",
    "Redfox" = "Red Fox",
    "RedSquirrel" = "Red Squirrel",
    "RichardsonsGroundSquirrel" = "Richardson's Ground Squirrel",
    "RiverOtter" = "River Otter",
    "SnowshoeHare" = "Snowshoe Hare",
    "StripedSkunk" = "Striped Skunk",
    "VolesMiceandAllies" = "Voles, Mice and Allies",
    "WeaselsandErmine" = "Weasels and Ermine",
    "WhitetailedDeer" = "White-tailed Deer",
    "WhitetailedJackRabbit" = "Whitetailed Jack Rabbit",
    "Wolverine" = "Wolverine",
    "WolvesCoyotesandAllies" = "Wolves, Coyotes and Allies",
    "WoodlandCaribou" = "Woodland Caribou")

tab <- data.frame(SpeciesID=names(SPP), CommonName=SPP)
for (i in d) {
    tab[[i]] <- names(SPP) %in% gsub(".png", "", ls[[i]])
}


library(jsonlite)
toJSON(tab, rownames=FALSE,pretty=TRUE)




library(cure4insect)
library(mefa4)
set_options(path = "s:/reports")
load_common_data()

tab <- get_species_table()
tab$DisplayName <- paste0(tab$CommonName, " (", tab$ScientificName, ")")
toJSON(tab[c("AlderFlycatcher", "Ovenbird"),], rownames=FALSE,pretty=TRUE)


gr <- "birds"
SPP <- rownames(resn[resn$Taxon == gr,])


## API ----------------------------------
## all species: display name, group, speciesID

library(jsonlite)

load("d:/abmi/reports/2018/misc/DataPortalUpdate.RData")
tab <- OUT$Species

dot <- endsWith(rownames(tab), ".")
dotBad <- rownames(tab)[dot]
dotOK <- substr(dotBad, 1, nchar(dotBad)-1)
data.frame(Bad=dotBad,Good=dotOK, tab$Group[dot])

for (i in seq_along(dotBad)) {
    levels(tab$SpeciesID)[levels(tab$SpeciesID) == dotBad[i]] <- dotOK[i]
    rownames(tab)[rownames(tab) == dotBad[i]] <- dotOK[i]
}
any(endsWith(rownames(tab), "."))

tab$DisplayName <- ifelse(is.na(tab$CommonName), as.character(tab$ScientificName),
    paste0(tab$CommonName, " (", tab$ScientificName, ")"))

grs <- c("birds", "vplants", "mites", "lichens", "mosses", "mammals-camera")

tmp <- list()
for (g in grs[1:5]) {
    p <- file.path("d:/abmi/reports/2018", "images", g)
    f <- list.files(p)
    s <- unique(sapply(strsplit(f, "-"), "[[", 1))
    tmp[[g]] <- tab[s,c("SpeciesID","DisplayName", "Group")]
    tmp[[g]]$sppprevious <- c(rownames(tmp[[g]])[nrow(tmp[[g]])], rownames(tmp[[g]])[-nrow(tmp[[g]])])
    tmp[[g]]$sppnext <- c(rownames(tmp[[g]])[-1], rownames(tmp[[g]])[1])
}


SPP <- c(
    "Badger" = "Badger",
    "Beaver" = "Beaver",
    "Bighornsheep" = "Bighorn Sheep",
    "Bison" = "Bison",
    "BlackBear" = "Black Bear",
    "Bobcat" = "Bobcat",
    "CanadaLynx" = "Canada Lynx",
    "ColumbianGroundSquirrel" = "Columbian Ground Squirrel",
    "Cougar" = "Cougar",
    "Coyote" = "Coyote",
    "Deer" = "Deer",
    "Elk" = "Elk",
    "Fisher" = "Fisher",
    "Foxes" = "Foxes",
    "GoldenMantledGroundSquirrel" = "Golden Mantled Ground Squirrel",
    "GrayWolf" = "Gray Wolf",
    "Grizzlybear" = "Grizzly Bear",
    "Groundhog" = "Groundhog",
    "HoaryMarmot" = "Hoary Marmot",
    "LeastChipmunk" = "Least Chipmunk",
    "Marten" = "Marten",
    "Mink" = "Mink",
    "Moose" = "Moose",
    "Mountaingoat" = "Mountain Goat",
    "Muledeer" = "Muledeer",
    "Muskrat" = "Muskrat",
    "NorthernFlyingSquirrel" = "Northern Flying Squirrel",
    "Porcupine" = "Porcupine",
    "Pronghorn" = "Pronghorn",
    "Raccoon" = "Raccoon",
    "Redfox" = "Red Fox",
    "RedSquirrel" = "Red Squirrel",
    "RichardsonsGroundSquirrel" = "Richardson's Ground Squirrel",
    "RiverOtter" = "River Otter",
    "SnowshoeHare" = "Snowshoe Hare",
    "StripedSkunk" = "Striped Skunk",
    "VolesMiceandAllies" = "Voles, Mice and Allies",
    "WeaselsandErmine" = "Weasels and Ermine",
    "WhitetailedDeer" = "White-tailed Deer",
    "WhitetailedJackRabbit" = "Whitetailed Jack Rabbit",
    "Wolverine" = "Wolverine",
    "WolvesCoyotesandAllies" = "Wolves, Coyotes and Allies",
    "WoodlandCaribou" = "Woodland Caribou")


z <- read.csv("d:/abmi/sppweb2018/Camera mammals Mar 2019/Mammal header table May 2019.csv")
rownames(z) <- z[,1]
mefa4::compare_sets(names(SPP), rownames(z))
setdiff(names(SPP), rownames(z))
SPP <- SPP[names(SPP) %in% rownames(z)]
z <- z[names(SPP),]

g <- grs[6]
tmp[[g]] <- data.frame(SpeciesID=names(SPP), DisplayName=SPP, Group=g)
tmp[[g]]$sppprevious <- c(rownames(tmp[[g]])[nrow(tmp[[g]])], rownames(tmp[[g]])[-nrow(tmp[[g]])])
tmp[[g]]$sppnext <- c(rownames(tmp[[g]])[-1], rownames(tmp[[g]])[1])

ALL <- do.call(rbind, tmp)

writeLines(toJSON(ALL, rownames=FALSE), file.path("d:/abmi/reports/2018", "api", "index.json"))

for (g in grs) {
    writeLines(toJSON(tmp[[g]], rownames=FALSE), file.path("d:/abmi/reports/2018", "api", g, "index.json"))
}

figs <- c("det", "useavail-north", "useavail-south",
    "coef-north", "coef-south", "map", "sector-north", "sector-south")

for (g in grs[1:5]) {
    s <- as.character(tmp[[g]]$SpeciesID)
    f <- list.files(file.path("d:/abmi/reports/2018", "images", g))
    m <- as.data.frame(matrix(FALSE, length(s), length(figs)))
    dimnames(m) <- list(s, figs)

    for (i in s) {
#        ff <- f[startsWith(f, i)]
#        ff <- gsub(".png", "", gsub(paste0(i, "-"), "", ff))
        d <- OUT$Species[i,]
        ff <- c("det",
            if (d$ModelNorth) c("sector-north", "coef-north") else NULL,
            if (d$ModelSouth) c("sector-south", "coef-south") else NULL,
            if (d$ModelNorth || d$ModelSouth) c("map") else NULL,
            if (d$UseavailNorth && !d$ModelNorth) c("useavail-north") else NULL,
            if (d$UseavailSouth && !d$ModelSouth) c("useavail-south") else NULL
        )
        m[i,] <- figs %in% ff
        v <- cbind(tab[i,], tmp[[g]][i,c("sppprevious", "sppnext")], m[i,])
        #toJSON(as.list(v), rownames=FALSE,pretty=TRUE,auto_unbox=TRUE)
        if (!dir.exists(file.path("d:/abmi/reports/2018", "api", g, i)))
            dir.create(file.path("d:/abmi/reports/2018", "api", g, i))
        writeLines(toJSON(as.list(v), rownames=FALSE, auto_unbox=TRUE,pretty=TRUE),
            file.path("d:/abmi/reports/2018", "api", g, i, "index.json"))
    }
}

g <- grs[6]
figs <- list.dirs(file.path("d:/abmi/reports/2018", "images", g),full.names=FALSE)[-1]
s <- as.character(tmp[[g]]$SpeciesID)
f <- list.files(file.path("d:/abmi/reports/2018", "images", g))
m <- as.data.frame(matrix(FALSE, length(s), length(figs)))
dimnames(m) <- list(s, figs)
for (j in figs) {
    fff <- gsub(".png", "", list.files(file.path("d:/abmi/reports/2018", "images", g, j)))
    fff <- fff[fff %in% names(SPP)]
    m[fff,j] <- TRUE
}
sum(as.matrix(m))
for (i in s) {
    if (!z[i,"ModelNorth"]) {
        m[i,"coef-north-combo"] <- FALSE
        m[i,"map-spclim-north"] <- FALSE
    }
    if (!z[i,"ModelNorthSummer"])
        m[i,"coef-north-summer"] <- FALSE
    if (!z[i,"ModelNorthWinter"])
        m[i,"coef-north-winter"] <- FALSE
    if (!z[i,"UseavailNorth"])
        m[i,"useavail-north"] <- FALSE

    if (!z[i,"ModelSouth"]) {
        m[i,"coef-south-combo-treed"] <- FALSE
        m[i,"coef-south-combo-nontreed"] <- FALSE
        m[i,"map-spclim-south"] <- FALSE
    }
    if (!z[i,"ModelSouthSummer"]) {
        m[i,"coef-south-summer-treed"] <- FALSE
        m[i,"coef-south-summer-nontreed"] <- FALSE
    }
    if (!z[i,"ModelSouthWinter"]) {
        m[i,"coef-south-winter-treed"] <- FALSE
        m[i,"coef-south-winter-nontreed"] <- FALSE
    }
    if (!z[i,"UseavailSouth"])
        m[i,"useavail-south"] <- FALSE

    if (!z[i,"ModelNorth"] && !z[i,"ModelSouth"]) {
        m[i,"map-cr"] <- FALSE
        m[i,"map-rf"] <- FALSE
        m[i,"map-df"] <- FALSE
    }
}
sum(as.matrix(m))

for (i in s) {
    #ff <- f[startsWith(f, i)]
    #ff <- gsub(".png", "", gsub(paste0(i, "-"), "", ff))
    v <- cbind(tmp[[g]][i,], m[i,])
    #toJSON(as.list(v), rownames=FALSE,pretty=TRUE,auto_unbox=TRUE)
    if (!dir.exists(file.path("d:/abmi/reports/2018", "api", g, i)))
        dir.create(file.path("d:/abmi/reports/2018", "api", g, i))
    writeLines(toJSON(as.list(v), rownames=FALSE, auto_unbox=TRUE,pretty=TRUE),
        file.path("d:/abmi/reports/2018", "api", g, i, "index.json"))
}

## clean up dotted names
if (FALSE) {
fl <- list.files("d:/abmi/reports/2018", recursive = TRUE)
xx <- NULL
for (i in seq_along(dotBad)) {
    xx <- c(xx, grep(dotBad[i], fl))
}
xx <- sort(unique(xx))
str(xx)
fl[xx]
for (j in xx) {
    In <- Out <- fl[j]
    for (i in seq_along(dotBad)) {
        Out <- gsub(dotBad[i], dotOK[i], Out)
    }
    file.rename(file.path("d:/abmi/reports/2018", In),
        file.path("d:/abmi/reports/2018", Out))
}
}
