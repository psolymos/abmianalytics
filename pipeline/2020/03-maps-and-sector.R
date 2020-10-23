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

Xclim_bird_S <- Xclim_bird[!kgrid$useN,]
Xclim_nonb_S <- Xclim_nonb[!kgrid$useN,]
pA <- kgrid[!kgrid$useN, "pAspen"]
Xclim_bird_N <- Xclim_bird[!kgrid$useS,]
Xclim_nonb_N <- Xclim_nonb[!kgrid$useS,]

## weights for overlap region
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
        "CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "Water",
        "Water", "Water", "Water", "Water", "Water"), Sector = c("Native",
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
        "Misc", "Misc", "Misc", "Misc", "Misc", "Misc")), row.names = c("DecidR",
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
## either sol inf is unknown or Sector is Forestry --> exclude these
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
Ps <- Ps[!kgrid$useN,]
summary(rowSums(Ps))

## here comes trVeg processing

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
Pn <- Pn[!kgrid$useS,]
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


## species specific part begins here

spp <- "AlderFlycatcher"
taxon <- "birds"

spp <- "Actaea.rubra"
taxon <- "vplants"

i <- 1

type <- "C" # combo species (N+S)
if (taxon == "birds") {
    M <- list(
        N=spp %in% dimnames( COEFS[[taxon]]$north$joint)[[1]],
        S=spp %in% dimnames( COEFS[[taxon]]$south$joint)[[1]])
    if (M$N & !M$S)
        type <- "N"
    if (!M$N & M$S)
        type <- "S"
    cfn <- if (type == "S")
        NULL else COEFS[[taxon]]$north$joint[spp,,]
    cfs <- if (type == "N")
        NULL else COEFS[[taxon]]$south$joint[spp,,]
    XclimS <- Xclim_bird_S
    XclimN <- Xclim_bird_N
    FUN <- poisson()$linkinv
} else {
    M <- list(
        N=spp %in% dimnames( COEFS[[taxon]]$north)[[1]],
        S=spp %in% dimnames( COEFS[[taxon]]$south)[[1]])
    if (M$N & !M$S)
        type <- "N"
    if (!M$N & M$S)
        type <- "S"
    cfn <- if (type == "S")
        NULL else COEFS[[taxon]]$north[spp,,]
    cfs <- if (type == "N")
        NULL else COEFS[[taxon]]$south[spp,,]
    XclimS <- Xclim_nonb_S
    XclimN <- Xclim_nonb_N
    FUN <- binomial()$linkinv
}

## bootstrap specific part

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

    musrf <- matrix(muscl + muspa, nrow=nrow(Psrf), ncol=ncol(Psrf))
    musrf <- t(t(musrf) + bsrf)
    NSrf <- as.matrix(groupSums(Psrf * FUN(musrf), 2, chSrf$sector_use))
} else {
    NScr <- NULL
    NSrf <- NULL
}
if (type != "S") {
    gc()
    ## north calculations for the i'th run
    #compare_sets(rownames(cfn), chVeg$cr2)
    bncr <- cfn[chVcr$cr2, i] # current land cover
    bnrf <- cfn[chVrf$rf2, i] # reference land cover
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
    NScr <- NNcr[match(rownames(kgrid), rownames(NScr)),]
    NSrf <- NNrf[match(rownames(kgrid), rownames(NSrf)),]
    NScr[is.na(NScr)] <- 0
    NSrf[is.na(NSrf)] <- 0
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


ksc <- dd_2018$soil_current
ksc <- ksc / ifelse(rowSums(ksc) > 0, rowSums(ksc), 1)

cns <- rownames(COEFS$mites$south[1,,])
cns0 <- cns[1:(which(cns=="Intercept")-2)]
cns1 <- cns[which(cns=="Intercept"):length(cns)]
cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]
cnn1 <- cnn[which(cnn=="Intercept"):length(cnn)]
all(cns1==cnn1)

LTs <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(LTs) <- LTs[,1]
LTs$Map <- as.character(LTs$UseInAnalysis)
LTs$Map[LTs$Map == "Well"] <- "Wellsites"
LTs$Map[LTs$Map %in% c("SoilWater", "HWater")] <- "Water"
LTs$Map[LTs$Map %in% c("SoilUnknown", "HFor")] <- "UNK"
LTs <- LTs[,c("Map", "Sector")]



## -- this is not for sector effects!
##
## cfs: coefficients (soil type x boot)
## cl: climate variables (cell x var)
## pa: pAspen values (vector of length #cell)
## ps: soil class areas or proportions
## i: bootstrap iteration (1 by default)
## returns a #cell x sector matrix
pred_south_nonbird <- function(cfs, cl, pa, ps, i=1) {
    CL <- make_clim(x)
    LTs <- lt$south
    sclass <- rownames(LTs)[!(LTs$Map %in% c("UNK", "Water"))]
    sclass <- intersect(colnames(ps), sclass)
    LTsx <- LTs[sclass,]
    w <- groupSums(ps[,sclass], 2, LTs[sclass,"Map"])
    rs <- rowSums(w)
    w <- w / ifelse(rs > 0, rs, 1)

    mucl <- drop(CL %*% cfs[colnames(CL),i])
    mupa <- pa * cfs["pAspen",i]
    mu <- as.matrix(w)
    mu[] <- mucl + mupa
    mu <- t(t(mu) + cfs[colnames(mu),i])
    rabunds <- rowSums(w * plogis(mu))



}

ss <- which(kgrid$NRNAME=="Grassland")[1:100]
i <- 1 # bootstrap run
## spp coef table with boot
cfs <- COEFS$vplants$south[1,,1:10]
## climate info
cl <- kgrid[ss,]
## soil class proportions
ps <- ksc[ss,]
## pAspen
pa <- kgrid[ss,"pAspen"]




lt1 <- nonDuplicated(chVeg, cr, TRUE)[,c("cr", "sector")]
lt2 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v2020.csv")
rownames(lt2) <- lt2[,1]
lt2 <- lt2[!endsWith(rownames(lt2), "0"),]
lt2[,1] <- NULL
lt2$cr <- as.character(lt2$UseInAnalysis)
lt2$cr <- gsub("Decid", "Deciduous", lt2$cr)
lt2$cr <- gsub("Spruce", "WhiteSpruce", lt2$cr)
lt2$cr <- gsub("Well", "Wellsites", lt2$cr)
lt2$cr <- gsub("9", "8", lt2$cr)

cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]
cnn1 <- cnn[which(cnn=="Intercept"):length(cnn)]
compare_sets(cnn0, lt2$cr)
setdiff(cnn0, lt2$cr)
setdiff(lt2$cr, cnn0)

lt2$Label <- as.character(lt2$cr)
lt2$Sector <- as.character(lt2$Sector)
v=lt2[,c("Label", "Sector")]

compare_sets(rownames(lt1), rownames(v))
setdiff(rownames(lt1), rownames(v))
setdiff(rownames(lt2), rownames(v))
v <- v[v$Label %in% chVeg$cr,]
compare_sets(cnn0, v$Label)
setdiff(cnn0, v$Label)
setdiff(v$Label, cnn0)

v <- v[v$Label %in% cnn0,]
dput(v)
