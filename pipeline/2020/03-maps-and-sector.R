
## mapping and sector effects

load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData")


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
    z
}

## south - non birds

lt <- list(
    south=structure(list(Map = c("ClaySub", "Other", "Other", "Other",
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


