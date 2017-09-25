## values: 2014_fine, 2014_coarse, 2012, 2010_coarse
HF_VERSION <- "2014v2_coarse"
PIXEL <- "km" # km or qs
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
SCALE <- paste0(PIXEL, HF_YEAR)

### Transition for 1K grid ------------------------------------------------

## this is based on the fix-fire fix-age0 version
## label collapsing as desired (swamp/wetland, ages?)

COL_LABEL <- if (PIXEL == "km")
    "Row_Col" else "LinkID"

source("~/repos/abmianalytics/R/veghf_functions.R")

load(file.path(ROOT, VER, "data", "analysis", paste0("kgrid_table_", PIXEL, ".Rdata")))
fl <- list.files(file.path(ROOT, VER, "data", "inter", "veghf", SCALE))

cc <- c("Row_Col","VEGAGEclass","VEGHFAGEclass","SOILclass","SOILHFclass","Shape_Area")

Start <- c(0:79*10+1, 802)

SCALE_NOTE <- if (PIXEL == "km")
    "1 km x 1 km prediction grid cells" else "QS prediction grid cells"
COL_LABEL <- if (PIXEL == "km")
    "Row_Col" else "LinkID"
if (HF_YEAR == 2014) {
COL_HABIT <- "Combined"
COL_SOIL <- "Soil_Type"
}
if (HF_YEAR == 2010) {
COL_HABIT <- "Combined_ChgByCWCS"
COL_SOIL <- "Soil_Type_1"
}

load(file.path(ROOT, VER, "data", "inter", "veghf", SCALE, fl[1]))
dd <- make_vegHF_wide_v6(d,
        col.label=COL_LABEL,
        col.year=HF_YEAR,
        col.HFyear="YEAR_MineCFOAgCutblock",
        col.HABIT=COL_HABIT,
        col.SOIL=COL_SOIL,
        sparse=TRUE,
        HF_fine=grepl("_fine", HF_VERSION),
        wide=FALSE)
ddd0 <- dd[character(0),cc]
xddd0 <- ddd0

for (s in 1:(length(Start)-1)) {
    cat("\n----------------------\n\nStarting block", s, "\n")
    for (i in Start[s]:(Start[s+1]-1)) {
        cat(i, "of", length(fl), "-", fl[i], "\t")
        flush.console()
        #d <- read.csv(file.path(ROOT, VER, "data", "kgrid", "tiles", fl[i]))
        load(file.path(ROOT, VER, "data", "inter", "veghf", SCALE, fl[i]))
        dd <- make_vegHF_wide_v6(d,
            col.label=COL_LABEL,
            col.year=HF_YEAR,
            col.HFyear="YEAR_MineCFOAgCutblock",
            col.HABIT=COL_HABIT,
            col.SOIL=COL_SOIL,
            sparse=TRUE,
            HF_fine=grepl("_fine", HF_VERSION),
            wide=FALSE)
        if (i == Start[s]) {
            dd0 <- dd[,cc]
        } else {
            dd0 <- rbind(dd0, dd[,cc])
        }
        cat("OK", nrow(dd0), "\n")
    }
    ddd0 <- rbind(ddd0, dd0)
    cat("\nFinished block", s, "dim:", nrow(ddd0), "\n")
    if (i %in% c(100, 200, 300, 400, 500, 600, 700, 801)) {
        save(ddd0, file=file.path(ROOT, VER,
            "data", "inter", "veghf", paste0("long_", SCALE), paste0("Long-part", i, ".Rdata")))
        ddd0 <- xddd0
        gc()
    }
}

## -- works on pre-saved chunks

load(file.path(ROOT, VER, "data", "analysis", paste0("kgrid_table_", PIXEL, ".Rdata")))
lu <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
su <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

lu$use_tr <- as.character(lu$ETA_UseInAnalysis_Sector)
allVegTr <- unique(c(lu$use_tr[!lu$IS_HF],
    paste0(rep(lu$use_tr[!lu$IS_HF], sum(lu$IS_HF)), "->",
    rep(lu$use_tr[lu$IS_HF], each=sum(!lu$IS_HF)))))
allVegTr <- allVegTr[!grepl("XXX", allVegTr)]
lu$use_tr <- as.factor(lu$use_tr)

su$use_tr <- as.character(su$Levels1)
su$use_tr[!is.na(su$HF)] <- as.character(su$SOILHF[!is.na(su$HF)])
allSoilTr <- unique(c(su$use_tr[is.na(su$HF)],
    paste0(rep(su$use_tr[is.na(su$HF)], sum(!is.na(su$HF))), "->",
    rep(su$use_tr[!is.na(su$HF)], each=sum(is.na(su$HF))))))
su$use_tr <- as.factor(su$use_tr)

#load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))

#Target0 <- c("WhiteSpruce0", "Decidious0", "Mixedwood0", "Pine0", "BlackSpruce")

fl3 <- list.files(file.path(ROOT, VER, "data", "inter", "veghf", paste0("long_", SCALE)))


## do one LUFxNSR class at a time and save it
#i <- "LowerAthabasca_CentralMixedwood"
for (ii in 1:nlevels(kgrid$LUFxNSR)) {
    i <- levels(kgrid$LUFxNSR)[ii]
    cat("\n---------\nStarting", i, "-", ii, "of", nlevels(kgrid$LUFxNSR))
    #j <- 4
    units <- list()
    sunits <- list()
    for (j in 1:length(fl3)) {
        cat("\n", j);flush.console()
        load(file.path(ROOT, VER, "data", "inter", "veghf", paste0("long_", SCALE), fl3[j]))
        flush.console()
        ddd0$LUFxNSR <- kgrid$LUFxNSR[match(ddd0$Row_Col, kgrid$Row_Col)]
        #nsr <- as.character(kgrid[which(kgrid$LUFxNSR==i)[1], "NSRNAME"])

        if (any(ddd0$LUFxNSR == i)) {
            cat(" processing ... ")

            xx <- ddd0[ddd0$LUFxNSR == i,,drop=FALSE]
            xx$Row_Col <- droplevels(xx$Row_Col)
            xx$LUFxNSR <- NULL

            xx$soil <- su$use_tr[match(xx$SOILclass, rownames(su))]
            xx$shf <- su$use_tr[match(xx$SOILHFclass, rownames(su))]

            xx$veg <- lu$use_tr[match(xx$VEGAGEclass, rownames(lu))]
            xx$vhf <- lu$use_tr[match(xx$VEGHFAGEclass, rownames(lu))]

            xx$soilTr <- ifelse(as.character(xx$soil) == as.character(xx$shf),
                as.character(xx$soil), paste0(as.character(xx$soil),
                "->", as.character(xx$shf)))

            xx$vegTr <- ifelse(as.character(xx$veg) == as.character(xx$vhf),
                as.character(xx$veg), paste0(as.character(xx$veg),
                "->", as.character(xx$vhf)))

            sxt <- Xtab(Shape_Area ~ Row_Col + soilTr, xx)
            sxxx <- Melt(sxt)
            colnames(sxxx) <- c("Row_Col", "soilTr", "Shape_Area")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")

            ch2veg <- data.frame(t(sapply(strsplit(as.character(xxx$vegTr), "->"),
                function(z) if (length(z)==1) z[c(1,1)] else z[1:2])))
            colnames(ch2veg) <- c("rf","cr")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xxx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")
            units[[j]] <- xxx
            sunits[[j]] <- sxxx
        } else cat(" onto the next chunk")
    }
    units <- do.call(rbind, units)
    levels(units$vegTr) <- c(levels(units$vegTr),
        setdiff(allVegTr, levels(units$vegTr)))
    sunits <- do.call(rbind, sunits)
    levels(sunits$soilTr) <- c(levels(sunits$soilTr),
        setdiff(allSoilTr, levels(sunits$soilTr)))

    trVeg <- Xtab(Shape_Area ~ Row_Col + vegTr, units)
    trVeg <- trVeg[,allVegTr]
    trSoil <- Xtab(Shape_Area ~ Row_Col + soilTr, sunits)
    trSoil <- trSoil[rownames(trVeg),allSoilTr]
    range(rowSums(trVeg)/10^6)
    range(rowSums(trSoil)/10^6)

    if (sum(trVeg[,grep("0", colnames(trVeg)),]) > 0)
        stop("Reference age 0 issue (3)")

    cat("\nSaving", i, "\n\n")
    flush.console()
    save(trVeg, trSoil,
        file=file.path(ROOT, VER, "data", "analysis", paste0("transitions_", SCALE),
        paste0(i, ".Rdata")))
}

## OSA to include
osagrid <- read.csv(file.path(ROOT, VER, "data", "raw", "xy", "Grid1km_in3OilSandRegions.csv"))
kgrid$OSANAME <- osagrid$SHORTNAME[match(kgrid$Row_Col, osagrid$Row_Col)]
table(kgrid$OSANAME)


ivals <- c("Athabasca Oilsand Area", "Cold Lake Oilsand Area", "Peace River Oilsand Area")
for (i in ivals) {
    cat("\n---------", i)
    #j <- 4
    units <- list()
    sunits <- list()
    for (j in 1:length(fl3)) {
        cat("\n", j);flush.console()
        load(file.path(ROOT, VER, "data", "inter", "veghf", paste0("long_", SCALE), fl3[j]))
        flush.console()
        ddd0$OSANAME <- kgrid$OSANAME[match(ddd0$Row_Col, kgrid$Row_Col)]
        #nsr <- as.character(kgrid[which(kgrid$OSANAME==i)[1], "NSRNAME"])

        if (any(ddd0$OSANAME == i)) {
            cat(" processing ... ")

            xx <- ddd0[ddd0$OSANAME == i,,drop=FALSE]
            xx$Row_Col <- droplevels(xx$Row_Col)
            xx$OSANAME <- NULL

            xx$soil <- su$use_tr[match(xx$SOILclass, rownames(su))]
            xx$shf <- su$use_tr[match(xx$SOILHFclass, rownames(su))]

            xx$veg <- lu$use_tr[match(xx$VEGAGEclass, rownames(lu))]
            xx$vhf <- lu$use_tr[match(xx$VEGHFAGEclass, rownames(lu))]

            xx$soilTr <- ifelse(as.character(xx$soil) == as.character(xx$shf),
                as.character(xx$soil), paste0(as.character(xx$soil),
                "->", as.character(xx$shf)))

            xx$vegTr <- ifelse(as.character(xx$veg) == as.character(xx$vhf),
                as.character(xx$veg), paste0(as.character(xx$veg),
                "->", as.character(xx$vhf)))

            sxt <- Xtab(Shape_Area ~ Row_Col + soilTr, xx)
            sxxx <- Melt(sxt)
            colnames(sxxx) <- c("Row_Col", "soilTr", "Shape_Area")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")

            ch2veg <- data.frame(t(sapply(strsplit(as.character(xxx$vegTr), "->"),
                function(z) if (length(z)==1) z[c(1,1)] else z[1:2])))
            colnames(ch2veg) <- c("rf","cr")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xxx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")
            units[[j]] <- xxx
            sunits[[j]] <- sxxx
        } else cat(" onto the next chunk")
    }
    units <- do.call(rbind, units)
    levels(units$vegTr) <- c(levels(units$vegTr),
        setdiff(allVegTr, levels(units$vegTr)))
    sunits <- do.call(rbind, sunits)
    levels(sunits$soilTr) <- c(levels(sunits$soilTr),
        setdiff(allSoilTr, levels(sunits$soilTr)))

    trVeg <- Xtab(Shape_Area ~ Row_Col + vegTr, units)
    trVeg <- trVeg[,allVegTr]
    trSoil <- Xtab(Shape_Area ~ Row_Col + soilTr, sunits)
    trSoil <- trSoil[rownames(trVeg),allSoilTr]
    range(rowSums(trVeg)/10^6)
    range(rowSums(trSoil)/10^6)

    if (sum(trVeg[,grep("0", colnames(trVeg)),]) > 0)
        stop("Reference age 0 issue (3)")

    cat("\nSaving", i, "\n\n")
    flush.console()
    save(trVeg, trSoil,
        file=file.path(ROOT, VER, "data", "analysis", paste0("transitions_", SCALE),
        paste0("00OSA ", i, ".Rdata")))
}

    units <- list()
    sunits <- list()
    for (j in 1:length(fl3)) {
        cat("\n", j);flush.console()
        load(file.path(ROOT, VER, "data", "inter", "veghf", paste0("long_", SCALE), fl3[j]))
        flush.console()
        ddd0$OSANAME <- kgrid$OSANAME[match(ddd0$Row_Col, kgrid$Row_Col)]
        #nsr <- as.character(kgrid[which(kgrid$OSANAME==i)[1], "NSRNAME"])

        if (any(ddd0$OSANAME %in% ivals)) {
            cat(" processing ... ")

            xx <- ddd0[ddd0$OSANAME %in% ivals,,drop=FALSE]
            xx$Row_Col <- droplevels(xx$Row_Col)
            xx$OSANAME <- NULL

            xx$soil <- su$use_tr[match(xx$SOILclass, rownames(su))]
            xx$shf <- su$use_tr[match(xx$SOILHFclass, rownames(su))]

            xx$veg <- lu$use_tr[match(xx$VEGAGEclass, rownames(lu))]
            xx$vhf <- lu$use_tr[match(xx$VEGHFAGEclass, rownames(lu))]

            xx$soilTr <- ifelse(as.character(xx$soil) == as.character(xx$shf),
                as.character(xx$soil), paste0(as.character(xx$soil),
                "->", as.character(xx$shf)))

            xx$vegTr <- ifelse(as.character(xx$veg) == as.character(xx$vhf),
                as.character(xx$veg), paste0(as.character(xx$veg),
                "->", as.character(xx$vhf)))

            sxt <- Xtab(Shape_Area ~ Row_Col + soilTr, xx)
            sxxx <- Melt(sxt)
            colnames(sxxx) <- c("Row_Col", "soilTr", "Shape_Area")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")

            ch2veg <- data.frame(t(sapply(strsplit(as.character(xxx$vegTr), "->"),
                function(z) if (length(z)==1) z[c(1,1)] else z[1:2])))
            colnames(ch2veg) <- c("rf","cr")

            xt <- Xtab(Shape_Area ~ Row_Col + vegTr, xxx)
            xxx <- Melt(xt)
            colnames(xxx) <- c("Row_Col", "vegTr", "Shape_Area")
            units[[j]] <- xxx
            sunits[[j]] <- sxxx
        } else cat(" onto the next chunk")
    }
    units <- do.call(rbind, units)
    levels(units$vegTr) <- c(levels(units$vegTr),
        setdiff(allVegTr, levels(units$vegTr)))
    sunits <- do.call(rbind, sunits)
    levels(sunits$soilTr) <- c(levels(sunits$soilTr),
        setdiff(allSoilTr, levels(sunits$soilTr)))

    trVeg <- Xtab(Shape_Area ~ Row_Col + vegTr, units)
    trVeg <- trVeg[,allVegTr]
    trSoil <- Xtab(Shape_Area ~ Row_Col + soilTr, sunits)
    trSoil <- trSoil[rownames(trVeg),allSoilTr]
    range(rowSums(trVeg)/10^6)
    range(rowSums(trSoil)/10^6)

    if (sum(trVeg[,grep("0", colnames(trVeg)),]) > 0)
        stop("Reference age 0 issue (3)")

    save(trVeg, trSoil,
        file=file.path(ROOT, VER, "data", "analysis", paste0("transitions_", SCALE),
        paste0("00OSA All3.Rdata")))
