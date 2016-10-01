make_vegHF_sector <- function(d, col.label, col.year=NULL, wide=TRUE, sparse=FALSE)
{
    dd <- make_vegHF_wide2(d, col.label, col.year=NULL, wide=FALSE)
    ## manipulate levels for veg
    hf <- dd$VEGHFAGEclass
    ## dont collapse because we need cutblock ages relative to recent year
    #levels(hf)[grepl("CC", levels(hf))] <- "CutBlocks"
    veglevs <- levels(dd$VEGAGEclass)
    hflevs <- setdiff(levels(hf), levels(dd$VEGAGEclass))
    dd$VEGchange <- interaction(dd$VEGAGEclass, hf, drop=FALSE, sep="->")
    levels(dd$VEGchange) <- sapply(strsplit(levels(dd$VEGchange), "->"),
        function(z) if (z[2] %in% hflevs) paste(z[1], z[2], sep="->") else z[1])
    ## manipulate levels for soil
    hf <- dd$SOILHFclass
    ## soil classes have CutBlocks without ages
    veglevs <- levels(dd$SOILclass)
    hflevs <- setdiff(levels(hf), levels(dd$SOILclass))
    dd$SOILchange <- interaction(dd$SOILclass, hf, drop=FALSE, sep="->")
    levels(dd$SOILchange) <- sapply(strsplit(levels(dd$SOILchange), "->"),
        function(z) if (z[2] %in% hflevs) paste(z[1], z[2], sep="->") else z[1])
    if (!wide)
        return(dd)
    ## do crosstab for veg and soil (no need to do cr and ref, it is combined)
    VegCh <- Xtab(Shape_Area ~ LABEL + VEGchange, dd)
    SoilCh <- Xtab(Shape_Area ~ LABEL + SOILchange, dd)
    out <- list(veg_change=VegCh,
        soil_change=SoilCh)
    ## have same checkpoints for wide and sparse
    if (!sparse)
        out <- lapply(out, as.matrix)
    out$sample_year <- THIS_YEAR
    out
}


make_vegHF_wide_v6 <-
function(d, col.label, col.year=NULL, col.HFyear=NULL, wide=TRUE, sparse=FALSE) {

    TreedClasses <- c("Conif", "Decid", "Fir", "Mixedwood", "Pine", "Spruce",
        "TreedBog-BSpr",
        "TreedFen-BSpr", "TreedFen-Decid", "TreedFen-Larch", "TreedFen-Mixedwood",
        "TreedSwamp-Conif", "TreedSwamp-Decid", "TreedSwamp-Fir",
        #"TreedSwamp-Forest",
        "TreedSwamp-Mixedwood", "TreedSwamp-Spruce",
        #"TreedWetland-Mixedwood",
        "AlpineLarch")
    NontreedClasses <- c("ShrubbyBog", "ShrubbyFen", "ShrubbySwamp",
        "Bog", "Marsh", "Swamp", "GraminoidFen",
        "GrassHerb", "Shrub", "Alkali",
        "Bare", "SnowIce",
        "Water")
    Fragment <- c("TreedWetland-Mixedwood", "TreedSwamp-Forest")
    HFLab <- c("BorrowpitsDugoutsSumps", "Canals", "CultivationCropPastureBareground",
        "HighDensityLivestockOperation", "IndustrialSiteRural", "MineSite",
        "MunicipalWaterSewage", "OtherDisturbedVegetation", "PeatMine",
        "Pipeline", "RailHardSurface", "RailVegetatedVerge", "Reservoirs",
        "RoadHardSurface", "RoadTrailVegetated", "RoadVegetatedVerge",
        "RuralResidentialIndustrial", "SeismicLine", "TransmissionLine",
        "Urban", "WellSite", "WindGenerationFacility")
    RfLab <- c(paste0(rep(TreedClasses, each=11),
        c("0","R","1","2","3","4","5","6","7","8","9")),
        Fragment, NontreedClasses)
    CrOnlyLab <- c(HFLab, paste0("CC", paste0(rep(TreedClasses, each=11),
        c("R","1","2","3","4"))))
    HLEVS <- c(TreedClasses, Fragment, NontreedClasses)

    ## designate a label column (there are different column names in use)
    d$LABEL <- d[,col.label]
    d$HF_Year <- d[,col.HFyear]
    if (any(is.na(d$LABEL)))
        stop("missing LABEL")
    #    d <- d[!is.na(d$LABEL),]
    ## designate a year column
    if (is.null(col.year)) {
        THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900
        d$SampleYear <- THIS_YEAR
    } else {
        if (is.numeric(col.year)) {
            if (length(col.year) > 1)
                stop("length of col.yeat > 1")
            THIS_YEAR <- col.year
            d$SampleYear <- THIS_YEAR
        } else {
            THIS_YEAR <- NA
            d$SampleYear <- d[,col.year]
        }
    }
    ## use upper-case labels for FEATURE_TY
    levels(d$FEATURE_TY) <- toupper(levels(d$FEATURE_TY))

    d$ORIGIN_YEAR <- d$Origin_Year
    #d$Origin_Year <- NULL
    d$HABIT <- d$c4
    #d$c4 <- NULL

    #### Footprint classes:
    ## check if we have all the feature types in the lookup table
    ## "" blank is for non-HF classes in current veg
    levels(d$FEATURE_TY)[levels(d$FEATURE_TY) == "''"] <- ""
    levels(d$FEATURE_TY)[levels(d$FEATURE_TY) == " "] <- ""
    if (!all(setdiff(levels(d$FEATURE_TY), rownames(hflt)) == ""))
        stop("HF diff:\n\t",
            dput(paste(setdiff(levels(d$FEATURE_TY), rownames(hflt)),
            collapse="\n\t", sep="")))
    ## classify feature types according to the chosen level of HF designation
    ## which comes from hf.level column of hflt (HF lookup table)
    d$HFclass <- hflt$HF_GROUP[match(d$FEATURE_TY, rownames(hflt))]
    ## HFclass inhgerits all levels from hflt[,hf.level]
    ## need to add in the blank for further parsing
    levels(d$HFclass) <- c(levels(d$HFclass), "")
    d$HFclass[is.na(d$HFclass)] <- ""

    ## slivers (tiny polys with no veg info):
    #stopifnot(max(d$Shape_Area[d$VEGclass == ""]) < 1)
    if (any(d$HABIT == ""))
        warning(paste("blank HABIT:", sum(d$Shape_Area[d$HABIT == ""]), "m^2"))
    d <- d[d$HABIT != "",]
    d$HABIT <- droplevels(d$HABIT)

    #### HABIT/EC classes:
    ## follow HABIT/EC classes, but there are few oddities when outside of AVI
    #d$VEGclass <- d$EC_Type
    d$VEGclass <- d$HABIT
    levels(d$VEGclass)[levels(d$VEGclass) == "Non-Veg"] <- "NonVeg"
    levels(d$VEGclass) <- gsub("/", "", levels(d$VEGclass))

    if (length(setdiff(d$VEGclass, HLEVS)) > 0)
        stop(paste("check HABIT classes", setdiff(d$VEGclass, HLEVS)))
    #tb <- cbind(c("", HLEVS), c("", HLEVS))
#return(setdiff(levels(d$VEGclass), tb[,1]))
    #levels(d$VEGclass) <- tb[match(levels(d$VEGclass), tb[,1]),2]

    #tmp <- aggregate(d$Shape_Area, list(lcc=d$EC_Type), sum)
    #tmp$p <- round(100*tmp$x/sum(tmp$x),2)
    #tmp2 <- aggregate(d$Shape_Area, list(lcc=d$VEGclass), sum)
    #tmp2$p <- round(100*tmp2$x/sum(tmp2$x),2)


#    NontreedClasses <- setdiff(levels(d$VEGclass), TreedClasses)
#    NontreedClasses <- NontreedClasses[NontreedClasses != ""]

    #### Age info for backfilled (Rf) and current (Cr)
    ## reference age class 0=no age (either not forest or no info)
    ## 1=0-19, 2=20-39, etc.
    d$AgeRf <- as.integer(sign(d$ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$ORIGIN_YEAR) / 20)))
    ## truncate reference age classes at 9 = 160+
    d$AgeRf[d$AgeRf > 9L] <- 9L
    ## placeholder for recent burn (0-9 years)
    tmp <- as.integer(sign(d$ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$ORIGIN_YEAR) / 10)))
    d$AgeRf[tmp == 1L] <- 999L
    ## set 0 year in treed habitats as max (assumed old forest)
#    d$AgeRf[d$AgeRf == 0L & d$VEGclass %in% TreedClasses] <- 9L
    ## unknown age is set to 0
    #table(d$AgeRf, d$VEGclass, useNA="a") # check NA ORIGIN_YEAR values
    #d$AgeRf[is.na(d$AgeRf)] <- 0L

    ## incorporate HF year for cutblocks
    d$CC_ORIGIN_YEAR <- d$ORIGIN_YEAR
    ii <- d$HFclass == "CutBlocks"
    ii[ii & !is.na(d$ORIGIN_YEAR) & d$HF_Year >= d$ORIGIN_YEAR] <- TRUE
    ii[ii & is.na(d$ORIGIN_YEAR)] <- TRUE
    d$CC_ORIGIN_YEAR[ii] <- d$HF_Year[ii]
    ## age for current with cutblock ages
    d$AgeCr <- as.integer(sign(d$CC_ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$CC_ORIGIN_YEAR) / 20)))
    ## truncate current age classes at 9
    d$AgeCr[d$AgeCr > 9L] <- 9L
    ## placeholder for recent CC (0-9 years)
    tmp <- as.integer(sign(d$CC_ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$CC_ORIGIN_YEAR) / 10)))
    d$AgeCr[tmp == 1L] <- 999L
    ## unknown age is set to 0
    #table(d$AgeCr, d$VEGclass, useNA="a") # check NA ORIGIN_YEAR values
    #d$AgeCr[is.na(d$AgeCr)] <- 0L
    #table(d$AgeCr,useNA="a")

    ## correcting reference age class based on cutblock info:
    ## these happened as a result of backfilling, so we accept HF age instead
    ## but this should be rare (ref age must be >= current)
    #ii <- !is.na(d$AgeCr) & d$AgeCr > d$AgeRf & d$AgeCr < 999L
    #if (sum(ii)>0)
    #    warning(paste("AgeCr > AgeRf for this many cases:", sum(ii)))
    #d$AgeRf[ii] <- d$AgeCr[ii]
    d$AgeRf[is.na(d$AgeRf)] <- 0L
    #table(rf=d$AgeRf,cr=d$AgeCr,useNA="a")
    ## turning age values into factor:
    ## 0=no age info,
    ## 1:9=valid age classes for treed veg types,
    ## ""=non-treed
    ## 999=placeholder for _R_ecent burn "R"
    d$AgeRf <- factor(d$AgeRf, levels=c(as.character(c(0:9, 999)), ""))
    ## NA --> "0" as unknown age class
    d$AgeRf[is.na(d$AgeRf)] <- "0"
    ## age is not relevant in non-treed veg types
    d$AgeRf[d$VEGclass %in% NontreedClasses] <- ""
    ## burn
    levels(d$AgeRf)[levels(d$AgeRf)=="999"] <- "R"

    ## making current age as factor
    d$AgeCr <- factor(d$AgeCr, levels=c(as.character(c(0:9, 999)), ""))
    ## NA --> "0" as unknown age class
    d$AgeCr[is.na(d$AgeCr)] <- "0"
    ## age is not relevant in non-treed veg types (no HF)
    d$AgeCr[d$VEGclass %in% NontreedClasses & d$HFclass == ""] <- ""
    ## age is not relevant outside of cutblocks
    d$AgeCr[!(d$HFclass %in% c("", "CutBlocks"))] <- ""
    ## recent CC
    levels(d$AgeCr)[levels(d$AgeCr)=="999"] <- "R"
    #table(current=d$AgeCr, reference=d$AgeRf)

    #### Combining VEG, HF and Age:
    ## reference VEG + Age labels:
    d$VEGAGEclass <- interaction(d$VEGclass, d$AgeRf, drop=TRUE, sep="", lex.order=TRUE)
    levels(d$VEGAGEclass) <- c(levels(d$VEGAGEclass),
        setdiff(RfLab, levels(d$VEGAGEclass)))

    ## manage CC labels
    ## current veg+hf
    d$VEGHFclass <- d$VEGclass
    #CClabels <- paste0("CC", levels(d$VEGclass)[levels(d$VEGclass) != ""])
    CClabels <- paste0("CC", levels(d$VEGclass))
    tmp <- setdiff(levels(d$HFclass), levels(d$VEGclass))
    tmp <- tmp[!(tmp %in% c("", "CutBlocks"))]
    levels(d$VEGHFclass) <- c(levels(d$VEGHFclass), tmp, CClabels)
    ## add non-CC HF types
    d$VEGHFclass[!(d$HFclass %in% c("", "CutBlocks"))] <- d$HFclass[!(d$HFclass %in% c("", "CutBlocks"))]
    ## should later the non-merchendisable forests with CC should be redistributed?
    ## e.g. after producing the wide format
    ## update CC labels obly for <= 80 yr CC (usually this does not happen
    ## just to make sure labels are OK)
    ## anything above age class >4 is turned into 4 to avoid labeling issues (shrubland)
    d$AgeCr[d$HFclass == "CutBlocks" & d$AgeCr %in% c("5","6","7","8","9")] <- "4"
    ii <- d$HFclass == "CutBlocks" & d$AgeCr %in% c("0","R","1","2","3","4")
    if (sum(ii) > 0)
        d$VEGHFclass[ii] <- paste0("CC", as.character(d$VEGclass[ii]))

    ## labels where backfilled cutblock label is not forested habitat
    ## right now I just collapse them to see % of the areas
    ## it is usually < 10% at this scale so it might be safe to ignore them
    ## usually young ages, but ranges R-1-2-3
    #ii <- unlist(lapply(paste0("CC", NontreedClasses), grep, x=levels(d$VEGHFclass)))
    #levels(d$VEGHFclass)[ii] <- "CCOpenTypes"
    ## unknown types under 'CC' considered as 'CCOpenTypes'
    #levels(d$VEGHFclass)[levels(d$VEGHFclass) == "CC"] <- "CCOpenTypes"
    ## treed wetlands
    #ii <- unlist(lapply(paste0("CC", TreedWetClasses), grep, x=levels(d$VEGHFclass)))
    #levels(d$VEGHFclass)[ii] <- "CCWetTypes"

    ## current VEG + HF + Age labels:
    d$VEGHFAGEclass <- interaction(d$VEGHFclass, d$AgeCr, drop=TRUE, sep="", lex.order=TRUE)
    ## labels with 0 age category are also to be fixed later ------> hard stuff
    #ii <- unlist(lapply(paste0(TreedClasses, 0), grep, x=levels(d$VEGHFAGEclass)))
    #levels(d$VEGHFAGEclass)[ii] <- "CCproblem"
    ## Labels for output columns
    AllLabels <- c(RfLab, CrOnlyLab)
    levels(d$VEGHFAGEclass) <- c(levels(d$VEGHFAGEclass), setdiff(AllLabels, levels(d$VEGHFAGEclass)))

    #### soils:
    SoilLab <- c("UNK", "Water", "BdL", "BlO", "CS", "Cy", "Gr", "LenA", "LenSP",
        "LenT", "LenS", "Li", "Lo", "LtcC", "LtcD", "LtcH", "LtcS", "Ov",
        "Sa", "Sb", "SL", "SwG", "Sy", "TB")

    #d$SOILclass <- d$SOIL_TYPE
    d$SOILclass <- d$Soil_Type_1
    ## need to have the UNKnown class to be able to deal with NAs
    if (!is.factor(d$SOILclass))
        d$SOILclass <- as.factor(d$SOILclass)
    if (!any(levels(d$SOILclass) == ""))
        levels(d$SOILclass) <- c(levels(d$SOILclass), "")
    ## dealing with NAs
    d$SOILclass[is.na(d$SOILclass)] <- ""
    ## unknown soil type outside of GVI and Dry Mixedwood
    levels(d$SOILclass)[levels(d$SOILclass) == ""] <- "UNK"
    levels(d$SOILclass)[levels(d$SOILclass) == " "] <- "UNK"
    ## get rid of modifiers
    levels(d$SOILclass) <- sapply(strsplit(levels(d$SOILclass), "-"), function(z) z[1L])
    ## add in Water label
    levels(d$SOILclass) <- c(levels(d$SOILclass), "Water")
    ## treat these as Water
    levels(d$SOILclass)[levels(d$SOILclass) %in% c("Len","LenW","Ltc","LtcR")] <- "Water"
    ## DEM/EC based Water class overrides soil
    d$SOILclass[d$VEGclass == "Water"] <- "Water"
    levels(d$SOILclass) <- c(levels(d$SOILclass), setdiff(SoilLab, levels(d$SOILclass)))
    d$SOILHFclass <- d$SOILclass
    levels(d$SOILHFclass) <- c(levels(d$SOILHFclass), levels(d$HFclass)[levels(d$HFclass) != ""])
    d$SOILHFclass[d$HFclass != ""] <- d$HFclass[d$HFclass != ""]
    SoilHFLab <- levels(d$SOILHFclass)
    ## NOTE: current UNK can be smaller than reference UNK, it can be turned into HF
    ## currently this is not tracked

    ## for point intersection or transition matrix processing, etc.
    if (!wide)
        return(d)

    #### crosstabs
    ## veg reference
    VegRf <- Xtab(Shape_Area ~ LABEL + VEGAGEclass, d)
    ## veg + HF current
    VegCr <- Xtab(Shape_Area ~ LABEL + VEGHFAGEclass, d)
    ## soils (`reference`)
    SoilRf <- Xtab(Shape_Area ~ LABEL + SOILclass, d)
    ## soils (`current`, soil + HF)
    SoilCr <- Xtab(Shape_Area ~ LABEL + SOILHFclass, d)

    rn <- rownames(VegRf) # make sure row labels are identical across tables
    VegRf <- VegRf[rn, RfLab, drop=FALSE]
    VegCr <- VegCr[rn, AllLabels, drop=FALSE]
    SoilRf <- SoilRf[rn, SoilLab, drop=FALSE]
    SoilCr <- SoilCr[rn, SoilHFLab, drop=FALSE]

    out <- list(veg_current=VegCr,
        veg_reference=VegRf,
        soil_current=SoilCr,
        soil_reference=SoilRf)
#        rs_veg_current=structure(as.integer(rsVegCr), names=names(rsVegCr)),
#        rs_veg_reference=structure(as.integer(rsVegRf), names=names(rsVegRf)),
#        rs_soil_current=structure(as.integer(rsSoilCr), names=names(rsSoilCr)),
#        rs_soil_reference=structure(as.integer(rsSoilRf), names=names(rsSoilRf)))
    if (!sparse)
        out <- lapply(out, as.matrix)

    ## year for each row
    tmp <- nonDuplicated(d, LABEL, TRUE)
    tmp <- tmp[rownames(VegCr),]
    #out$sample_year <- THIS_YEAR
    out$sample_year <- tmp$SampleYear

    out
}

make_vegHF_wide <-
function(d, col.label, col.year=NULL, col.HFyear=NULL, wide=TRUE, sparse=FALSE) {

    RfLab <- c("Conif0", "ConifR", "Conif1", "Conif2", "Conif3", "Conif4", "Conif5", "Conif6",
        "Conif7", "Conif8", "Conif9",
        "Decid0", "DecidR", "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6",
        "Decid7", "Decid8", "Decid9",
        "Mixwood0", "MixwoodR", "Mixwood1", "Mixwood2", "Mixwood3", "Mixwood4", "Mixwood5",
        "Mixwood6", "Mixwood7", "Mixwood8", "Mixwood9",
        "Pine0", "PineR", "Pine1", "Pine2", "Pine3", "Pine4",
        "Pine5", "Pine6", "Pine7", "Pine8", "Pine9",

        "GrassHerb", "NonVeg", "Shrub", "Water",

        "Swamp-Conif0", "Swamp-ConifR", "Swamp-Conif1", "Swamp-Conif2",
        "Swamp-Conif3", "Swamp-Conif4", "Swamp-Conif5", "Swamp-Conif6",
        "Swamp-Conif7", "Swamp-Conif8", "Swamp-Conif9",
        "Swamp-Decid0", "Swamp-DecidR", "Swamp-Decid1", "Swamp-Decid2",
        "Swamp-Decid3", "Swamp-Decid4", "Swamp-Decid5", "Swamp-Decid6",
        "Swamp-Decid7", "Swamp-Decid8", "Swamp-Decid9",
        "Swamp-Mixwood0", "Swamp-MixwoodR", "Swamp-Mixwood1", "Swamp-Mixwood2",
        "Swamp-Mixwood3", "Swamp-Mixwood4", "Swamp-Mixwood5",
        "Swamp-Mixwood6", "Swamp-Mixwood7", "Swamp-Mixwood8", "Swamp-Mixwood9",
        "Swamp-Pine0", "Swamp-PineR", "Swamp-Pine1", "Swamp-Pine2",
        "Swamp-Pine3", "Swamp-Pine4",
        "Swamp-Pine5", "Swamp-Pine6", "Swamp-Pine7", "Swamp-Pine8", "Swamp-Pine9",


        "Wetland-Bare", "Wetland-GrassHerb", "Wetland-Shrub",

        "Wetland-BSpr0", "Wetland-BSprR", "Wetland-BSpr1", "Wetland-BSpr2", "Wetland-BSpr3",
        "Wetland-BSpr4", "Wetland-BSpr5", "Wetland-BSpr6", "Wetland-BSpr7",
        "Wetland-BSpr8", "Wetland-BSpr9",
        "Wetland-Decid0", "Wetland-DecidR", "Wetland-Decid1", "Wetland-Decid2", "Wetland-Decid3",
        "Wetland-Decid4", "Wetland-Decid5", "Wetland-Decid6", "Wetland-Decid7",
        "Wetland-Decid8", "Wetland-Decid9",
        "Wetland-Larch0", "Wetland-LarchR", "Wetland-Larch1", "Wetland-Larch2", "Wetland-Larch3",
        "Wetland-Larch4", "Wetland-Larch5", "Wetland-Larch6", "Wetland-Larch7",
        "Wetland-Larch8", "Wetland-Larch9")

    CrOnlyLab <- c("BorrowpitsDugoutsSumps", "Canals", "CultivationCropPastureBareground",
        "HighDensityLivestockOperation", "IndustrialSiteRural", "MineSite",
        "MunicipalWaterSewage", "OtherDisturbedVegetation", "PeatMine",
        "Pipeline", "RailHardSurface", "RailVegetatedVerge", "Reservoirs",
        "RoadHardSurface", "RoadTrailVegetated", "RoadVegetatedVerge",
        "RuralResidentialIndustrial", "SeismicLine", "TransmissionLine",
        "Urban", "WellSite", "WindGenerationFacility",
        "CCDecid0", "CCDecidR", "CCDecid1", "CCDecid2",
        "CCDecid3", "CCDecid4",
        "CCMixwood0", "CCMixwoodR", "CCMixwood1", "CCMixwood2",
        "CCMixwood3", "CCMixwood4",
        "CCConif0", "CCConifR", "CCConif1", "CCConif2",
        "CCConif3", "CCConif4",
        "CCPine0", "CCPineR", "CCPine1", "CCPine2",
        "CCPine3", "CCPine4")

    ## designate a label column (there are different column names in use)
    d$LABEL <- d[,col.label]
    d$HF_Year <- d[,col.HFyear]
    if (any(is.na(d$LABEL)))
        stop("missing LABEL")
    #    d <- d[!is.na(d$LABEL),]
    ## designate a year column
    if (is.null(col.year)) {
        THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900
        d$SampleYear <- THIS_YEAR
    } else {
        if (is.numeric(col.year)) {
            if (length(col.year) > 1)
                stop("length of col.yeat > 1")
            THIS_YEAR <- col.year
            d$SampleYear <- THIS_YEAR
        } else {
            THIS_YEAR <- NA
            d$SampleYear <- d[,col.year]
        }
    }
    ## use upper-case labels for FEATURE_TY
    levels(d$FEATURE_TY) <- toupper(levels(d$FEATURE_TY))

    #### Footprint classes:
    ## check if we have all the feature types in the lookup table
    ## "" blank is for non-HF classes in current veg
    levels(d$FEATURE_TY)[levels(d$FEATURE_TY) == "''"] <- ""
    levels(d$FEATURE_TY)[levels(d$FEATURE_TY) == " "] <- ""
    if (!all(setdiff(levels(d$FEATURE_TY), rownames(hflt)) == ""))
        stop("HF diff:\n\t",
            dput(paste(setdiff(levels(d$FEATURE_TY), rownames(hflt)),
            collapse="\n\t", sep="")))
    ## classify feature types according to the chosen level of HF designation
    ## which comes from hf.level column of hflt (HF lookup table)
    d$HFclass <- hflt$HF_GROUP[match(d$FEATURE_TY, rownames(hflt))]
    ## HFclass inhgerits all levels from hflt[,hf.level]
    ## need to add in the blank for further parsing
    levels(d$HFclass) <- c(levels(d$HFclass), "")
    d$HFclass[is.na(d$HFclass)] <- ""

    ## slivers (tiny polys with no veg info):
    #stopifnot(max(d$Shape_Area[d$VEGclass == ""]) < 1)
    if (any(d$HABIT == ""))
        warning(paste("blank HABIT:", sum(d$Shape_Area[d$HABIT == ""]), "m^2"))
    d <- d[d$HABIT != "",]
    d$HABIT <- droplevels(d$HABIT)

    #### HABIT/EC classes:
    ## follow HABIT/EC classes, but there are few oddities when outside of AVI
    #d$VEGclass <- d$EC_Type
    d$VEGclass <- d$HABIT
    levels(d$VEGclass)[levels(d$VEGclass) == "Non-Veg"] <- "NonVeg"
    levels(d$VEGclass) <- gsub("/", "", levels(d$VEGclass))

    ## April 26 2015 update:
    ## complicated rule set with these levels:
    HLEVS <- c("Conif", "Decid", "GrassHerb", "Mixwood", "NonVeg", "Pine",
        "Shrub", "Swamp-Conif", "Swamp-Decid", "Swamp-Mixwood", "Swamp-Pine",
        "Water", "Wetland-Bare", "Wetland-BSpr", "Wetland-Decid", "Wetland-GrassHerb",
        "Wetland-Larch", "Wetland-Shrub")
    if (length(setdiff(d$VEGclass, HLEVS)) > 0)
        stop(paste("check HABIT classes", setdiff(d$VEGclass, HLEVS)))
    #tb <- cbind(c("", HLEVS), c("", HLEVS))
#return(setdiff(levels(d$VEGclass), tb[,1]))
    #levels(d$VEGclass) <- tb[match(levels(d$VEGclass), tb[,1]),2]

    #tmp <- aggregate(d$Shape_Area, list(lcc=d$EC_Type), sum)
    #tmp$p <- round(100*tmp$x/sum(tmp$x),2)
    #tmp2 <- aggregate(d$Shape_Area, list(lcc=d$VEGclass), sum)
    #tmp2$p <- round(100*tmp2$x/sum(tmp2$x),2)


    ## treed and non-treed veg classes where having age makes sense
#    TreedClasses <- c("BlackSpruce", "WhiteSpruce", "Deciduous",
#        "Mixedwood", "Pine", "LarchFen")
    TreedClasses <- c("Conif", "Decid", "Mixwood", "Pine",
        "Swamp-Conif", "Swamp-Decid", "Swamp-Mixwood", "Swamp-Pine",
        "Wetland-BSpr", "Wetland-Decid",
        "Wetland-Larch")
    ## this is needed to assign "CCWetTypes"
    TreedWetClasses <- c("Swamp-Conif", "Swamp-Decid", "Swamp-Mixwood", "Swamp-Pine",
        "Wetland-BSpr", "Wetland-Decid", "Wetland-Larch")
    NontreedClasses <- setdiff(levels(d$VEGclass), TreedClasses)
    NontreedClasses <- NontreedClasses[NontreedClasses != ""]

    #### Age info for backfilled (Rf) and current (Cr)
    ## reference age class 0=no age (either not forest or no info)
    ## 1=0-19, 2=20-39, etc.
    d$AgeRf <- as.integer(sign(d$ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$ORIGIN_YEAR) / 20)))
    ## truncate reference age classes at 9 = 160+
    d$AgeRf[d$AgeRf > 9L] <- 9L
    ## placeholder for recent burn (0-9 years)
    tmp <- as.integer(sign(d$ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$ORIGIN_YEAR) / 10)))
    d$AgeRf[tmp == 1L] <- 999L
    ## set 0 year in treed habitats as max (assumed old forest)
#    d$AgeRf[d$AgeRf == 0L & d$VEGclass %in% TreedClasses] <- 9L
    ## unknown age is set to 0
    #table(d$AgeRf, d$VEGclass, useNA="a") # check NA ORIGIN_YEAR values
    #d$AgeRf[is.na(d$AgeRf)] <- 0L

    ## incorporate HF year for cutblocks
    d$CC_ORIGIN_YEAR <- d$ORIGIN_YEAR
    ii <- d$HFclass == "CutBlocks"
    ii[ii & !is.na(d$ORIGIN_YEAR) & d$HF_Year >= d$ORIGIN_YEAR] <- TRUE
    ii[ii & is.na(d$ORIGIN_YEAR)] <- TRUE
    d$CC_ORIGIN_YEAR[ii] <- d$HF_Year[ii]
    ## age for current with cutblock ages
    d$AgeCr <- as.integer(sign(d$CC_ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$CC_ORIGIN_YEAR) / 20)))
    ## truncate current age classes at 9
    d$AgeCr[d$AgeCr > 9L] <- 9L
    ## placeholder for recent CC (0-9 years)
    tmp <- as.integer(sign(d$CC_ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$CC_ORIGIN_YEAR) / 10)))
    d$AgeCr[tmp == 1L] <- 999L
    ## unknown age is set to 0
    #table(d$AgeCr, d$VEGclass, useNA="a") # check NA ORIGIN_YEAR values
    #d$AgeCr[is.na(d$AgeCr)] <- 0L
    #table(d$AgeCr,useNA="a")

    ## correcting reference age class based on cutblock info:
    ## these happened as a result of backfilling, so we accept HF age instead
    ## but this should be rare (ref age must be >= current)
    #ii <- !is.na(d$AgeCr) & d$AgeCr > d$AgeRf & d$AgeCr < 999L
    #if (sum(ii)>0)
    #    warning(paste("AgeCr > AgeRf for this many cases:", sum(ii)))
    #d$AgeRf[ii] <- d$AgeCr[ii]
    d$AgeRf[is.na(d$AgeRf)] <- 0L
    #table(rf=d$AgeRf,cr=d$AgeCr,useNA="a")
    ## turning age values into factor:
    ## 0=no age info,
    ## 1:9=valid age classes for treed veg types,
    ## ""=non-treed
    ## 999=placeholder for _R_ecent burn "R"
    d$AgeRf <- factor(d$AgeRf, levels=c(as.character(c(0:9, 999)), ""))
    ## NA --> "0" as unknown age class
    d$AgeRf[is.na(d$AgeRf)] <- "0"
    ## age is not relevant in non-treed veg types
    d$AgeRf[d$VEGclass %in% NontreedClasses] <- ""
    ## burn
    levels(d$AgeRf)[levels(d$AgeRf)=="999"] <- "R"

    ## making current age as factor
    d$AgeCr <- factor(d$AgeCr, levels=c(as.character(c(0:9, 999)), ""))
    ## NA --> "0" as unknown age class
    d$AgeCr[is.na(d$AgeCr)] <- "0"
    ## age is not relevant in non-treed veg types (no HF)
    d$AgeCr[d$VEGclass %in% NontreedClasses & d$HFclass == ""] <- ""
    ## age is not relevant outside of cutblocks
    d$AgeCr[!(d$HFclass %in% c("", "CutBlocks"))] <- ""
    ## recent CC
    levels(d$AgeCr)[levels(d$AgeCr)=="999"] <- "R"
    #table(current=d$AgeCr, reference=d$AgeRf)

    #### Combining VEG, HF and Age:
    ## reference VEG + Age labels:
    d$VEGAGEclass <- interaction(d$VEGclass, d$AgeRf, drop=TRUE, sep="", lex.order=TRUE)
    levels(d$VEGAGEclass) <- c(levels(d$VEGAGEclass),
        setdiff(RfLab, levels(d$VEGAGEclass)))

    ## manage CC labels
    ## current veg+hf
    d$VEGHFclass <- d$VEGclass
    #CClabels <- paste0("CC", levels(d$VEGclass)[levels(d$VEGclass) != ""])
    CClabels <- paste0("CC", levels(d$VEGclass))
    tmp <- setdiff(levels(d$HFclass), levels(d$VEGclass))
    tmp <- tmp[!(tmp %in% c("", "CutBlocks"))]
    levels(d$VEGHFclass) <- c(levels(d$VEGHFclass), tmp, CClabels)
    ## add non-CC HF types
    d$VEGHFclass[!(d$HFclass %in% c("", "CutBlocks"))] <- d$HFclass[!(d$HFclass %in% c("", "CutBlocks"))]
    ## should later the non-merchendisable forests with CC should be redistributed?
    ## e.g. after producing the wide format
    ## update CC labels obly for <= 80 yr CC (usually this does not happen
    ## just to make sure labels are OK)
    ## anything above age class >4 is turned into 4 to avoid labeling issues (shrubland)
    d$AgeCr[d$HFclass == "CutBlocks" & d$AgeCr %in% c("5","6","7","8","9")] <- "4"
    ii <- d$HFclass == "CutBlocks" & d$AgeCr %in% c("0","R","1","2","3","4")
    if (sum(ii) > 0)
        d$VEGHFclass[ii] <- paste0("CC", as.character(d$VEGclass[ii]))

    ## labels where backfilled cutblock label is not forested habitat
    ## right now I just collapse them to see % of the areas
    ## it is usually < 10% at this scale so it might be safe to ignore them
    ## usually young ages, but ranges R-1-2-3
    #ii <- unlist(lapply(paste0("CC", NontreedClasses), grep, x=levels(d$VEGHFclass)))
    #levels(d$VEGHFclass)[ii] <- "CCOpenTypes"
    ## unknown types under 'CC' considered as 'CCOpenTypes'
    #levels(d$VEGHFclass)[levels(d$VEGHFclass) == "CC"] <- "CCOpenTypes"
    ## treed wetlands
    #ii <- unlist(lapply(paste0("CC", TreedWetClasses), grep, x=levels(d$VEGHFclass)))
    #levels(d$VEGHFclass)[ii] <- "CCWetTypes"

    ## current VEG + HF + Age labels:
    d$VEGHFAGEclass <- interaction(d$VEGHFclass, d$AgeCr, drop=TRUE, sep="", lex.order=TRUE)
    ## labels with 0 age category are also to be fixed later ------> hard stuff
    #ii <- unlist(lapply(paste0(TreedClasses, 0), grep, x=levels(d$VEGHFAGEclass)))
    #levels(d$VEGHFAGEclass)[ii] <- "CCproblem"
    ## Labels for output columns
    AllLabels <- c(RfLab, CrOnlyLab)
    levels(d$VEGHFAGEclass) <- c(levels(d$VEGHFAGEclass), setdiff(AllLabels, levels(d$VEGHFAGEclass)))

    #### soils:
    SoilLab <- c("UNK", "Water", "BdL", "BlO", "CS", "Cy", "Gr", "LenA", "LenSP",
        "LenT", "LenS", "Li", "Lo", "LtcC", "LtcD", "LtcH", "LtcS", "Ov",
        "Sa", "Sb", "SL", "SwG", "Sy", "TB")

    d$SOILclass <- d$SOIL_TYPE
    ## need to have the UNKnown class to be able to deal with NAs
    if (!is.factor(d$SOILclass))
        d$SOILclass <- as.factor(d$SOILclass)
    if (!any(levels(d$SOILclass) == ""))
        levels(d$SOILclass) <- c(levels(d$SOILclass), "")
    ## dealing with NAs
    d$SOILclass[is.na(d$SOILclass)] <- ""
    ## unknown soil type outside of GVI and Dry Mixedwood
    levels(d$SOILclass)[levels(d$SOILclass) == ""] <- "UNK"
    levels(d$SOILclass)[levels(d$SOILclass) == " "] <- "UNK"
    ## get rid of modifiers
    levels(d$SOILclass) <- sapply(strsplit(levels(d$SOILclass), "-"), function(z) z[1L])
    ## add in Water label
    levels(d$SOILclass) <- c(levels(d$SOILclass), "Water")
    ## treat these as Water
    levels(d$SOILclass)[levels(d$SOILclass) %in% c("Len","LenW","Ltc","LtcR")] <- "Water"
    ## DEM/EC based Water class overrides soil
    d$SOILclass[d$VEGclass == "Water"] <- "Water"
    levels(d$SOILclass) <- c(levels(d$SOILclass), setdiff(SoilLab, levels(d$SOILclass)))
    d$SOILHFclass <- d$SOILclass
    levels(d$SOILHFclass) <- c(levels(d$SOILHFclass), levels(d$HFclass)[levels(d$HFclass) != ""])
    d$SOILHFclass[d$HFclass != ""] <- d$HFclass[d$HFclass != ""]
    SoilHFLab <- levels(d$SOILHFclass)
    ## NOTE: current UNK can be smaller than reference UNK, it can be turned into HF
    ## currently this is not tracked

    ## for point intersection or transition matrix processing, etc.
    if (!wide)
        return(d)

    #### crosstabs
    ## veg reference
    VegRf <- Xtab(Shape_Area ~ LABEL + VEGAGEclass, d)
    ## veg + HF current
    VegCr <- Xtab(Shape_Area ~ LABEL + VEGHFAGEclass, d)
    ## soils (`reference`)
    SoilRf <- Xtab(Shape_Area ~ LABEL + SOILclass, d)
    ## soils (`current`, soil + HF)
    SoilCr <- Xtab(Shape_Area ~ LABEL + SOILHFclass, d)

    rn <- rownames(VegRf) # make sure row labels are identical across tables
    VegRf <- VegRf[rn, RfLab, drop=FALSE]
    VegCr <- VegCr[rn, AllLabels, drop=FALSE]
    SoilRf <- SoilRf[rn, SoilLab, drop=FALSE]
    SoilCr <- SoilCr[rn, SoilHFLab, drop=FALSE]

    out <- list(veg_current=VegCr,
        veg_reference=VegRf,
        soil_current=SoilCr,
        soil_reference=SoilRf)
#        rs_veg_current=structure(as.integer(rsVegCr), names=names(rsVegCr)),
#        rs_veg_reference=structure(as.integer(rsVegRf), names=names(rsVegRf)),
#        rs_soil_current=structure(as.integer(rsSoilCr), names=names(rsSoilCr)),
#        rs_soil_reference=structure(as.integer(rsSoilRf), names=names(rsSoilRf)))
    if (!sparse)
        out <- lapply(out, as.matrix)

    ## year for each row
    tmp <- nonDuplicated(d, LABEL, TRUE)
    tmp <- tmp[rownames(VegCr),]
    #out$sample_year <- THIS_YEAR
    out$sample_year <- tmp$SampleYear

    out
}


## this is used to bind matrices whose columns are the same
## multiple rows for same rownames are added up
## must have dimnames, but that is given from Xtab
## Note: melt-recast is more efficient than groupSums
## This also handles repeat rows that might happen at QS level
bind_fun <- function(x) {
    if (inherits(x[[1]], "data.frame")) {
        out <- do.call(rbind, x)
    }
    if (inherits(x[[1]], "Matrix")) {
        cn <- lapply(x, colnames)
        if (length(cn) > 1) {
            for (i in 2:length(cn))
                if (!all(cn[[i]] == cn[[i-1]]))
                    stop("colnames must be equal")
        }
        mel <- lapply(x, Melt)
        mel <- do.call(rbind, mel)
        out <- Xtab(value ~ rows + cols, mel)
        out <- out[,colnames(x[[1]])]
    }
    out
}
bind_fun2 <- function(x, y, check.col=TRUE) {
    if (check.col &&
        length(union(colnames(x), colnames(y))) != length(intersect(colnames(x), colnames(y))))
            stop("colnames must be same", dput(setdiff(union(colnames(x), colnames(y)),
                intersect(colnames(x), colnames(y)))))
    melx <- Melt(x)
    mely <- Melt(y)
    mel <- rbind(melx, mely)
    out <- Xtab(value ~ rows + cols, mel)
    out <- out[,colnames(x)]
    out
}

c_fun <- function(x, y) {
    df <- data.frame(SUM=c(x,y), nam=c(names(x), names(y)),
        w=c(rep(1L, length(x)), rep(2L, length(y))))
    xt <- Xtab(SUM ~ nam+w, df)
    structure(rowSums(xt), names=rownames(xt))
}

## redistributing unknown ages among known ages
fill_in_0ages <-
function(x, NSR)
{
    Target <- c("Conif", "Decid", "Mixwood", "Pine",
        "Swamp-Conif", "Swamp-Decid", "Swamp-Mixwood", "Swamp-Pine",
        "Wetland-BSpr", "Wetland-Decid", "Wetland-Larch")
    Ages <- c("0", "R", as.character(1:9))
    NSR <- droplevels(as.factor(NSR))
    NSRs <- levels(NSR)
    for (current in c(TRUE, FALSE)) {
        xx <- if (current)
            x$veg_current else x$veg_reference
        xx <- as.matrix(xx)
        ag <- if (current)
            AvgAges$current else AvgAges$reference
        for (nsr in NSRs) {
            cat(ifelse(current, "current:", "reference:"), nsr)
            flush.console()
            for (i in Target) {
                Cols <- paste0(i, Ages)
                j <- NSR == nsr
                if (any(j)) {
                    Mat <- xx[j, Cols, drop=FALSE]
                    ## multiply Mat[,1] (unknown age) with this matrix
                    Unk <- Mat[,1] * t(matrix(ag[i,,nsr], length(Ages), sum(j)))
                    Mat[,1] <- 0 # will be 0 and redistributed from Unk
                    Mat <- Mat + Unk
                    xx[j, Cols] <- Mat # ridiculously slow as sparse matrix
                }
            }
            cat(" --- OK\n")
        }
        if (current) {
            x$veg_current <- as(xx, "dgCMatrix")
        } else {
            x$veg_reference <- as(xx, "dgCMatrix")
        }
    }
    x
}
