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
    north=list())


## process south monster matrix and find some efficiencies

compare_sets(chSoil$cr, rownames(lt$south))
compare_sets(chSoil$rf, rownames(lt$south))
chSoil$cr2 <- lt$south$Label[match(chSoil$cr, rownames(lt$south))]
chSoil$rf2 <- lt$south$Label[match(chSoil$rf, rownames(lt$south))]
#chSoil$sector2 <- lt$south$Sector[match(chSoil$cr, rownames(lt$south))]
#with(chSoil, table(sector, sector2)) # sector definition is up to date

## define sectors to be used here
chSoil$sector_use <- chSoil$sector

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

## here comes the N/S weight variable

## figure out if a species is N/S/Combo

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

if (type != "N") {
    ## south calculations for the i'th run
    #compare_sets(rownames(cfs), chSoil$cr2)
    bscr <- cfs[chSoil$cr2, i] # current land cover
    bsrf <- cfs[chSoil$rf2, i] # reference land cover
    ## space-climate coefs
    bscl <- if (taxon == "birds")
        cfs[colnames(Xclim_bird), i] else cfs[colnames(Xclim_nonb), i]
    bspa <- cfs["pAspen", i]
    ## additive components for south
    muscl <- drop(XclimS %*% bscl)
    muspa <- pA * bspa
    mus <- as.matrix(0*Ps)
    mus[] <- muscl + muspa
    muscr <- t(t(mus) + bscr)
    musrf <- t(t(mus) + bsrf)
    NScr <- as.matrix(groupSums(Ps * FUN(muscr), 2, chSoil$sector_use))
    NSrf <- as.matrix(groupSums(Ps * FUN(musrf), 2, chSoil$sector_use))
} else {
    NScr <- NULL
    NSrf <- NULL
}
if (type != "S") {
    ## north calculations for the i'th run
    #compare_sets(rownames(cfn), chVeg$cr2)
    bncr <- cfn[chVeg$cr2, i] # current land cover
    bnrf <- cfn[chVeg$rf2, i] # reference land cover
    ## space-climate coefs
    bncl <- if (taxon == "birds")
        cfn[colnames(Xclim_bird), i] else cfn[colnames(Xclim_nonb), i]
    ## additive components for north
    muncl <- drop(XclimN %*% bncl)
    mun <- as.matrix(0*Pn)
    mun[] <- muscl
    muncr <- t(t(mun) + bncr)
    munrf <- t(t(mun) + nsrf)
    NNcr <- as.matrix(groupSums(Pn * FUN(muncr), 2, chVeg$sector_use))
    NNrf <- as.matrix(groupSums(Pn * FUN(munrf), 2, chVeg$sector_use))
} else {
    NNcr <- NULL
    NNrf <- NULL
}

## combine NS and NN together (weighted avg in overlap zone) for species with Combo (N+S)
if (type == "C") {
    # averaging comes here
    NNcr <- NNcr[rownames(kgrid),]
    NNrf <- NNrf[rownames(kgrid),]
    NNcr[is.na(NNcr)] <- 0
    NNrf[is.na(NNrf)] <- 0
    NScr <- NNcr[rownames(kgrid),]
    NSrf <- NNrf[rownames(kgrid),]
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




## plots


## south
cns <- c("Loamy", "SandyLoam", "RapidDrain", "ClaySub", "ThinBreak", "Blowout", "Other",
    "Crop", "TameP", "RoughP",
    "EnSeismic", "EnSoftLin", "TrSoftLin", "HardLin",
    "Wellsites", "Rural", "Urban", "Industrial",
    "Mine", "MineV", "Water")
cnn <- c(
    "WhiteSpruceR", "WhiteSpruce1", "WhiteSpruce2", "WhiteSpruce3",
    "WhiteSpruce4", "WhiteSpruce5", "WhiteSpruce6", "WhiteSpruce7",
    "WhiteSpruce8", "PineR", "Pine1", "Pine2", "Pine3", "Pine4",
    "Pine5", "Pine6", "Pine7", "Pine8", "DeciduousR", "Deciduous1",
    "Deciduous2", "Deciduous3", "Deciduous4", "Deciduous5", "Deciduous6",
    "Deciduous7", "Deciduous8", "MixedwoodR", "Mixedwood1", "Mixedwood2",
    "Mixedwood3", "Mixedwood4", "Mixedwood5", "Mixedwood6", "Mixedwood7",
    "Mixedwood8", "TreedBogR", "TreedBog1", "TreedBog2", "TreedBog3",
    "TreedBog4", "TreedBog5", "TreedBog6", "TreedBog7", "TreedBog8",
    "TreedSwamp", "ShrubbySwamp", "ShrubbyBog", "ShrubbyFen", "GraminoidFen",
    "Marsh", "Shrub", "CCWhiteSpruceR", "CCWhiteSpruce1", "CCWhiteSpruce2",
    "CCWhiteSpruce3", "CCWhiteSpruce4", "CCPineR", "CCPine1", "CCPine2",
    "CCPine3", "CCPine4", "CCDeciduousR", "CCDeciduous1", "CCDeciduous2",
    "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR", "CCMixedwood1",
    "CCMixedwood2", "CCMixedwood3", "CCMixedwood4", "Crop", "TameP",
    "RoughP", "Wellsites", "EnSeismic", "EnSoftLin", "TrSoftLin",
    "HardLin", "TreedFenR", "TreedFen1", "TreedFen2", "TreedFen3",
    "TreedFen4", "TreedFen5", "TreedFen6", "TreedFen7", "TreedFen8",
    "Rural", "Urban", "Industrial", "Mine", "MineV", "Water", "GrassHerb",
    "Intercept", "Lat", "Long", "AHM", "PET", "FFP", "MAP", "MAT",
    "MCMT", "MWMT", "Lat2", "Long2", "LatLong", "MAPPET", "MATAHM",
    "MAPFFP", "MAT2", "MWMT2")

taxon <- "vplants"
spp <- "Achillea.borealis"
cf <- COEFS[[taxon]]$south[spp,,]


FUN <- plogis

ra <- rbind(data.frame(pAspen="Non-treed", get_stats(FUN(cf[cns,]))),
    data.frame(pAspen="Treed", get_stats(FUN(t(t(cf[cns,]) + cf["pAspen",])))))

library(ggplot2)

p <- ggplot(ra, aes(x=Label, y=First, fill=Label, group=pAspen)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  width=0.2, color="#105A73") +
    ylab("Relative abundance") +
    xlab("Land cover class") +
    facet_wrap(~pAspen) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_minimal() +
    labs(title=paste0(spp, " (", taxon, ")")) +
    theme(legend.position="none")


