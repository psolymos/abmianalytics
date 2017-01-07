source("~/repos/abmianalytics/veghf/veghf-setup.R")

### 1K grid --------------------------------------------------------

## Sample year is current year, so that forest ages are relative to present
## and not relative to HF or veg inventory year.

#fl <- list.files(file.path(ROOT, VER, "data", "kgrid-V6", "tiles"))
fl <- list.files(file.path(ROOT, VER, "data", "kgrid-V6dec", "tiles-rdata"))

## test feature types and save in Rdata format
if (FALSE) {

library(data.table)
getOption("datatable.fread.datatable")
options(datatable.fread.datatable=FALSE)
getOption("datatable.fread.datatable")

fl0 <- list.files(file.path(ROOT, VER, "data", "kgrid-V6dec", "tiles"))
HF <- character(0)
VEG <- character(0)
A1 <- numeric(length(fl0))
A2 <- numeric(length(fl0))
for (i in seq_len(length(fl0))) {
    fn <- fl0[i]
    cat("\nchecking", i, "/", length(fl0));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
    #d <- fread(f)
    d <- read.csv(f)
    d2 <- c4_fun(d)
    HF <- union(HF, levels(d$FEATURE_TY))
    VEG <- union(VEG, levels(d2$c4))
    A1[i] <- sum(d$Shape_Area)
    A2[i] <- sum(d2$Shape_Area)
    save(d, file=file.path(ROOT, VER, "data", "kgrid-V6dec", "tiles-rdata",
        gsub(".csv", ".Rdata", fn, fixed = TRUE)))
}

setdiff(HF, c("", rownames(hflt)))
summary(A1-A2)
## todo
## - check c4 lookup table and rules
## - make sure that Areas are conserved
## - do crosstab
## - prepare maps

recl <- read.csv("c:/Users/Peter/Dropbox/abmi/V6/veg-V6-combined3.csv")

HF <- character(0)
C4 <- character(0)
for (i in seq_len(length(fl))) {
    fn <- fl[i]
    cat("\nchecking", i, "/", length(fl));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid-V6dec", "tiles-rdata", fn)
    e <- new.env()
    load(f, envir=e)
    d <- e$d
    HF <- union(HF, levels(d$FEATURE_TY))
    VEG <- union(VEG, levels(d$c4))
}

NEW <- character(0)
HF <- character(0)
VEG <- character(0)
VEG3 <- character(0)
SOIL <- character(0)
LARCH <- numeric(0)
CWCS <- character(0)
SRC <- character(0)
blank_n <- numeric(length(fl)) # no. of cases
blank_a <- numeric(length(fl)) # total area
nro <- 0
for (i in seq_len(length(fl))) {
    fn <- fl[i]
    cat("\nchecking", i, "/", length(fl));flush.console()
#    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles-rdata", fn)
#    d <- read.csv(f)
    e <- new.env()
    load(f, envir=e)
    d <- e$d
    diff <- setdiff(levels(d$FEATURE_TY), c("", rownames(hflt)))
    if (length(diff)) {
        NEW <- union(NEW, diff)
        cat("\n\t", length(NEW), "new types found\n")
    }
    d$VEG3 <- interaction(d$Veg_Type, d$Moisture_Reg, d$PreBackfill_Source,
        drop=TRUE, sep="::")
    d$c3 <- recl$Combined[match(d$VEG3, recl$Veg_3)]
    d$needCWCS <- recl$need_CWCS[match(d$VEG3, recl$Veg_3)] == "x"
    tmp <- data.frame(cwcs=d$CWCS_Class[d$needCWCS], c4=d$c3[d$needCWCS],
        src=d$PostBackfill_Source[d$needCWCS])
    tmp$c4x <- as.character(tmp$c4)

    ## grass/GraminoidWetland
    tmp$c4x[tmp$cwcs == "Marsh" & tmp$c4 == "GraminoidWetland"] <- "Marsh"
    tmp$c4x[tmp$cwcs == "Bog" & tmp$c4 == "GraminoidWetland"] <- "GraminoidBog"
    tmp$c4x[tmp$cwcs == "Fen" & tmp$c4 == "GraminoidWetland"] <- "GraminoidFen"
    #tmp$c4x[!(tmp$cwcs %in% c("Marsh","Fen","Bog")) & tmp$c4 == "GraminoidWetland"] <-
    #    "Marsh"
    ## nothing else remains, so no need for setting unknowns to Marsh

    ## shrub/ShrubbyWetland
    tmp$c4x[tmp$cwcs %in% c("Fen","Bog") & tmp$c4 == "ShrubbyWetland"] <-
        paste0("Shrubby",
        as.character(tmp$cwcs[tmp$cwcs %in% c("Fen","Bog") & tmp$c4 == "ShrubbyWetland"]))
    #tmp$c4x[!(tmp$cwcs %in% c("Fen","Bog")) & tmp$c4 == "ShrubbyWetland"] <-
    #    "ShrubbySwamp"
    ## nothing else remains, so no need for setting unknowns to ShrubbySwamp

    ## muskeg
    #tmp$c4x[tmp$cwcs == "Fen" & tmp$c4 == "Muskeg"] <-
    #    "GraminoidFen"
    #tmp$c4x[tmp$cwcs != "Fen" & tmp$c4 == "Muskeg"] <-
    #    "GraminoidBog" # can be: "None"   "Soil"   "ABMILC" "AVIE" "PLVI"  "Phase1" "MTNP"  "EINP"
    #tmp$c4x[!(tmp$cwcs %in% c("Swamp","Fen","Bog")) & tmp$c4 == "Muskeg" &
    #    tmp$src != "WBNP"] <- "GraminoidFen"

    d$c4 <- as.character(d$c3)
    d$c4[d$needCWCS] <- tmp$c4x

    ## treedwetland-mixedwood
    d$c4[d$c4 == "TreedWetland-Mixedwood" & d$cwcs == "Fen"] <- "TreedFen-Mixedwood"
    d$c4[d$c4 == "TreedWetland-Mixedwood" & d$cwcs == "Bog"] <- "TreedBog-BSpr"
    d$c4[d$c4 == "TreedWetland-Mixedwood" & d$cwcs == "Swamp"] <- "TreedSwamp-Mixedwood"
    ## call the rest (not fen/bog/swamp) as TreedSwamp-Mixedwood ?

    ## larch
    d$c4[d$Pct_of_Larch >= 5] <- "TreedFen-Larch"
    d$c4[d$Pct_of_Larch > 0 & d$c4 == "TreedBog-BSpr"] <- "TreedFen-BSpr"
    d$c4[d$Pct_of_Larch > 0 & d$Pct_of_Larch < 5 & d$c4 %in% c("Decid",
        "Fir","Mixedwood","Pine","Spruce")] <- paste0("TreedSwamp-", d$c4[d$Pct_of_Larch > 0 &
        d$Pct_of_Larch < 5 & d$c4 %in% c("Decid", "Fir","Mixedwood","Pine","Spruce")])

    #d$c4 <- factor(d$c4, c(levels(d$c3)[!(levels(d$c3) %in%
    #    c("GraminoidWetland","ShrubbyWetland","Muskeg"))], "GraminoidBog"))
    d$c4 <- factor(d$c4, c(levels(d$c3), "GraminoidBog", "TreedSwamp-Pine"))
    if (any(is.na(d$c4))) break

    #d$VEG5 <- interaction(d$Veg_Type, d$Moisture_Reg, d$PreBackfill_Source,
    #    d$CWCS_Class, ifelse(d$Pct_of_Larch>0, "Larch1", "Larch0"), drop=TRUE, sep="::")
    HF <- union(HF, levels(d$FEATURE_TY))
    VEG <- union(VEG, levels(d$Veg_Type))
    #VEG3 <- union(VEG3, levels(d$VEG3))
    VEG3 <- union(VEG3, levels(droplevels(d$c4[d$CutYear != 0])))
    SOIL <- union(SOIL, levels(d$Soil_Type_1))
    CWCS <- union(CWCS, levels(d$CWCS_Class))
    #SRC <- union(SRC, levels(d$PostBackfill_Source))
    SRC <- union(SRC, levels(droplevels(tmp$src)))
    tmp <- d[d$HABIT == "",,drop=FALSE]
    if (nrow(tmp)) {
        blank_n[i] <- nrow(tmp)
        blank_a[i] <- sum(tmp$Shape_Area)
        cat("\n\tfound blanks:", nrow(tmp), "\n")
    }
#    LARCH <- sort(unique(c(LARCH, d$Pct_of_Larch)))
    cat("\n\tveg:", length(VEG), #"- veg^3:", length(VEG3), "- larch:", length(LARCH),
        "- soil:", length(SOIL), "- hf:", length(HF), "\n")
    #print(summary(d$Pct_of_Larch))

    nro <- nro + nrow(d)
#    cc <- factor(ifelse(d$CutYear != 0, "CC", "F"), c("CC","F"))
    cc <- factor(ifelse(d$FEATURE_TY == "CUTBLOCK", "CC", "F"), c("CC","F"))
    #lr <- factor(ifelse(d$Pct_of_Larch != 0, "Larch", "Not"), c("Larch","Not"))
    lr <- cut(d$Pct_of_Larch, c(-1, 0, 4, 10))
    levels(lr) <- c("PctL0","PctL1-4","PctL5-10")
    if (i == 1) {
        nro2 <- Xtab(Shape_Area ~ c4 + cc, d)
        nro3 <- Xtab(Shape_Area ~ c3 + cc, d)
    } else {
        nro2 <- nro2 + Xtab(Shape_Area ~ c4 + cc, d)
        nro3 <- nro3 + Xtab(Shape_Area ~ c3 + cc, d)
    }
    if (i == 1) {
        al <- Xtab(Shape_Area ~ c4 + lr, d)
        al3 <- Xtab(Shape_Area ~ c3 + lr, d)
    } else {
        al <- al + Xtab(Shape_Area ~ c4 + lr, d)
        al3 <- al3 + Xtab(Shape_Area ~ c3 + lr, d)
    }

    d$VEG3 <- d$c3 <- d$needCWCS <- NULL
    save(d, file=file.path(ROOT, VER, "data", "kgrid-V6", "tiles-rdata",
        gsub(".csv", ".Rdata", fn, fixed = TRUE)))
}

write.csv(as.matrix(nro2/10^6), file=file.path(ROOT, VER, "data", "kgrid-V6", "veg-V6-cutblocks.csv"))
write.csv(as.matrix(nro3/10^6), file=file.path(ROOT, VER, "data", "kgrid-V6", "veg-V6-cutblocks3cc.csv"))
write.csv(as.matrix(al3/10^6), file=file.path(ROOT, VER, "data", "kgrid-V6", "veg-V6-larch3.csv"))
write.csv(as.matrix(al/10^6), file=file.path(ROOT, VER, "data", "kgrid-V6", "veg-V6-larch.csv"))

## no blanks
(blanks <- which(blank_n > 0))
sum(blank_a)
blank_a[blanks]
## no NEW HF labels
NEW
#x <- data.frame(Veg_Type=sort(VEG))
#write.csv(x, row.names=FALSE, file=file.path(ROOT, VER, "data", "kgrid-V6", "veg-V6.csv"))

MAP <- read.csv("c:/Users/Peter/Dropbox/abmi/V6/veg-V6-reclass.csv")
x <- data.frame(Veg_3=sort(VEG3))
tmp <- strsplit(sort(VEG3), "::")
x$Veg_Type <- sapply(tmp, "[[", 1)
x$Moisture_Reg <- sapply(tmp, "[[", 2)
x$PreBackfill_Source <- sapply(tmp, "[[", 3)
x$Veg_V5 <- MAP$Combined[match(x$Veg_Type, MAP$Veg_Type)]
x <- droplevels(x[!(x$Veg_Type %in% c("Black spruce", "Lentic shrub")),])
write.csv(x, row.names=FALSE, file="c:/Users/Peter/Dropbox/abmi/V6/veg-V6-combined3.csv")

## reclass (v6->v5)
MAP <- read.csv("c:/Users/Peter/Dropbox/abmi/V6/veg-V6-reclass.csv")
for (fn in fl) {
    cat("checking", which(fl == fn), "/", length(fl), "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
    d <- read.csv(f)
    d$HABIT <- reclass(d$Veg_Type, map=MAP, all=TRUE)
    dd <- make_vegHF_wide_v6(d, col.label="Row_Col",
        col.year=HF_YEAR, col.HFyear="CutYear", wide=FALSE)
    tmp <- colSums(is.na(dd[,c("VEGAGEclass",
        "VEGHFAGEclass","SOILclass","SOILHFclass")]))
    natrack[[fn]] <- tmp
    if (tmp[1] > 0)
        break
}

## tracking the strange area mismatch
natrack <- list()
for (fn in fl) {
    cat("checking", which(fl == fn), "/", length(fl), "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
    d <- read.csv(f)
    dd <- make_vegHF_wide(d, col.label="Row_Col",
        col.year=NULL, col.HFyear="CutYear", wide=FALSE)
    tmp <- colSums(is.na(dd[,c("VEGAGEclass",
        "VEGHFAGEclass","SOILclass","SOILHFclass")]))
    natrack[[fn]] <- tmp
    if (tmp[1] > 0)
        break
}

}

## processing csv files in batches of 50

Start <- c(1, 51, 101, 151, 201, 251, 301, 351, 401, 451,
    501, 551, 601, 651, 701, 751, 802)
tmplist <- list()

for (s in 1:(length(Start)-1)) {

    gc()
    fn <- fl[Start[s]]
    cat("\n\n------------- batch", s, "----------------\n")
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
#    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
#    d <- read.csv(f)
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles-rdata", fn)
    e <- new.env()
    load(f, envir=e)
    d <- e$d
    ## HF year is used as base year for prediction purposes
    Basic <- Xtab(Shape_Area ~ Row_Col + c4, d)
    dd <- make_vegHF_wide_v6(d, col.label="Row_Col",
        col.year=HF_YEAR, col.HFyear="CutYear", sparse=TRUE)
    veg_current <- dd$veg_current
    veg_reference <- dd$veg_reference
    soil_current <- dd$soil_current
    soil_reference <- dd$soil_reference
    sample_year <- dd$sample_year[1]

#lapply(dd[1:4], sum)

    for (i in (Start[s]+1):(Start[s+1]-1)) {

        fn <- fl[i]
        cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
#        f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
#        d <- read.csv(f)
        f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles-rdata", fn)
        e <- new.env()
        load(f, envir=e)
        d <- e$d
        Basic2 <- Xtab(Shape_Area ~ Row_Col + c4, d)
        Basic <- bind_fun2(Basic, Basic2)
        dd <- make_vegHF_wide_v6(d, col.label="Row_Col",
            col.year=HF_YEAR, col.HFyear="CutYear", sparse=TRUE)
        veg_current <- bind_fun2(veg_current, dd$veg_current)
        veg_reference <- bind_fun2(veg_reference, dd$veg_reference)
        soil_current <- bind_fun2(soil_current, dd$soil_current)
        soil_reference <- bind_fun2(soil_reference, dd$soil_reference)

    }
    tmplist[[s]] <- list(
        veg_current = veg_current,
        veg_reference = veg_reference,
        soil_current = soil_current,
        soil_reference = soil_reference,
        sample_year = sample_year,
        scale = "1 km x 1 km prediction grid cells",
        c4 = Basic)
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
Basic <- tmplist[[1]]$c4
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
    Basic <- bind_fun2(Basic, tmplist[[j]]$c4)
}

## assembling return object
dd1km_pred <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year = tmplist[[1]]$sample_year,
    scale = "1 km x 1 km prediction grid cells",
    v6veg = Basic)

## this has the climate stuff
kgrid <- read.csv(
    file.path(ROOT, VER, "data", "kgrid-V6",
    "Grid1km_template_final_clippedBy_ABBound_with_atts_to_Peter.csv"))
rownames(kgrid) <- kgrid$Row_Col
## this is the correct lat/long (i.e. more decimal places)
kgrid2 <- read.csv(
    file.path(ROOT, VER, "data", "kgrid-V6",
    "Grid1km_template_latLong.csv"))
rownames(kgrid2) <- kgrid2$Row_Col
kgrid2 <- kgrid2[rownames(kgrid),]
stopifnot(all(rownames(kgrid) == rownames(kgrid2)))
kgrid$POINT_X <- kgrid2$POINT_X
kgrid$POINT_Y <- kgrid2$POINT_Y
rm(kgrid2)

## this is has Brandt boreal and BCR
kgrid2 <- read.csv(
    file.path(ROOT, VER, "data", "kgrid-V6",
    "Grid1km_template_final_clippedBy_ABBound_with_atts_to_Peter_BCR_BRANDT_Done.csv"))
kgrid2 <- nonDuplicated(kgrid2, Row_Col, TRUE)
#rownames(kgrid2) <- kgrid2$Row_Col
kgrid2 <- kgrid2[rownames(kgrid),]
stopifnot(all(rownames(kgrid) == rownames(kgrid2)))
kgrid$TYPE_BRANDT <- kgrid2$TYPE_BRANDT
#kgrid$BCR <- kgrid2$BCR
#kgrid$BCR_NAME <- kgrid2$BCR_NAME
kgrid$BCRCODE <- kgrid2$BCRCODE
rm(kgrid2)

compare.sets(rownames(dd1km_pred$veg_current), rownames(kgrid))
#stopifnot(all(rownames(dd1km_pred$veg_current) == rownames(kgrid)))

## NSR x LUF regions used as prediction regions in sector effects etc.
kgrid$nsr_luf <- with(kgrid, paste(as.integer(NSRNAME), as.integer(LUF_NAME), sep="_"))
colnames(kgrid)[colnames(kgrid) == "col"] <- "Col"
colnames(kgrid)[colnames(kgrid) == "Eref"] <- "PET"
colnames(kgrid)[colnames(kgrid) == "Populus_tremuloides_brtpred_nofp"] <- "pAspen"

## 10 x 10 km grid
kgrid$Row10 <- 1 + kgrid$Row %/% 10
kgrid$Col10 <- 1 + kgrid$Col %/% 10
kgrid$Row10_Col10 <- interaction(kgrid$Row10, kgrid$Col10, sep="_", drop=TRUE)

## random pick from 10K grid
tmp <- as.integer(kgrid$Row10_Col10)
kgrid$Rnd10 <- integer(length(tmp))
set.seed(1234)
for (i in seq_len(max(tmp))) {
    lg <- tmp == i
    kgrid$Rnd10[lg] <- sample.int(sum(lg))
}

#load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
dd1km_pred$veg_current <- dd1km_pred$veg_current[rownames(kgrid),]
dd1km_pred$veg_reference <- dd1km_pred$veg_reference[rownames(kgrid),]
dd1km_pred$soil_current <- dd1km_pred$soil_current[rownames(kgrid),]
dd1km_pred$soil_reference <- dd1km_pred$soil_reference[rownames(kgrid),]
dd1km_pred$v6veg <- dd1km_pred$v6veg[rownames(kgrid),]

## check area diff
range(sapply(dd1km_pred[1:4], sum) / 10^6)

## proportion of water -- for mapping purposes
kgrid$pWater <- dd1km_pred$veg_current[,"Water"] / rowSums(dd1km_pred$veg_current)

kgrid$pSoil <- 1 - rowSums(dd1km_pred$soil_reference[,c("UNK","Water")]) / rowSums(dd1km_pred$soil_reference)

## veg based area < soil based area, thus using the max
kgrid$Area_km2 <- rowSums(dd1km_pred$veg_current) / 10^6

kgrid$NEAR_DIST <- NULL

## UTM projection for fake maps
library(raster)
library(sp)
library(rgdal)
XYlatlon <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(XYlatlon) <- ~ POINT_X + POINT_Y
proj4string(XYlatlon) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
XY <- as.data.frame(spTransform(XYlatlon, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))
kgrid$X <- XY$POINT_X
kgrid$Y <- XY$POINT_Y


## fill-in NA values with nearest

lnas <- is.na(kgrid[,"pAspen"])
wnas <- which(!lnas)
for (i in which(lnas)) {
    j <- wnas[which.min(sqrt((kgrid$X[!lnas] - kgrid$X[i])^2 +
        (kgrid$Y[!lnas] - kgrid$Y[i])^2))]
    kgrid[i,"pAspen"] <- kgrid[j,"pAspen"]
}

cvs <- c("AHM", "PET", "FFP", "MAP", "MAT", "MCMT", "MWMT")
lnas <- is.na(kgrid[,cvs[1]])
wnas <- which(!lnas)
for (i in which(lnas)) {
    j <- wnas[which.min(sqrt((kgrid$X[!lnas] - kgrid$X[i])^2 +
        (kgrid$Y[!lnas] - kgrid$Y[i])^2))]
    kgrid[i,cvs] <- kgrid[j,cvs]
}

sum(is.na(kgrid))

kgrid$LUFxNSR <- interaction(kgrid$LUF_NAME, kgrid$NSRNAME, drop=TRUE, sep="_")
levels(kgrid$LUFxNSR) <- gsub(" ", "", levels(kgrid$LUFxNSR))

## topo variables
kgrid2 <- read.csv(
    file.path(ROOT, VER, "data", "kgrid-V6",
    "Grid1kmCenter_topo.csv"))
rownames(kgrid2) <- kgrid2$Row_Col
kgrid2 <- kgrid2[rownames(kgrid),]
cvs <- c("slope","slpasp","tri","cti")
lnas <- is.na(kgrid2[,cvs[3]])
wnas <- which(!lnas)
for (i in which(lnas)) {
    j <- wnas[which.min(sqrt((kgrid$X[!lnas] - kgrid$X[i])^2 +
        (kgrid$Y[!lnas] - kgrid$Y[i])^2))]
    kgrid2[i,cvs] <- kgrid2[j,cvs]
}
sum(is.na(kgrid2))

kgrid$SLP <- kgrid2$slope
kgrid$ASP <- kgrid2$slpasp
kgrid$TRI <- kgrid2$tri
kgrid$CTI <- kgrid2$cti

dd1km_nsr <- dd1km_pred
dd1km_nsr$veg_current <- groupSums(dd1km_pred$veg_current, 1, kgrid$NSRNAME)
dd1km_nsr$veg_reference <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NSRNAME)
dd1km_nsr$soil_current <- groupSums(dd1km_pred$soil_current, 1, kgrid$NSRNAME)
dd1km_nsr$soil_reference <- groupSums(dd1km_pred$soil_reference, 1, kgrid$NSRNAME)
dd1km_nsr$v6veg <- groupSums(dd1km_pred$v6veg, 1, kgrid$NSRNAME)

dd1km_nr <- dd1km_pred
dd1km_nr$veg_current <- groupSums(dd1km_pred$veg_current, 1, kgrid$NRNAME)
dd1km_nr$veg_reference <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NRNAME)
dd1km_nr$soil_current <- groupSums(dd1km_pred$soil_current, 1, kgrid$NRNAME)
dd1km_nr$soil_reference <- groupSums(dd1km_pred$soil_reference, 1, kgrid$NRNAME)
dd1km_nr$v6veg <- groupSums(dd1km_pred$v6veg, 1, kgrid$NRNAME)

if (SAVE) { ## needed for recalculating average ages
    save(dd1km_pred,
        file=file.path(ROOT, VER, "data/kgrid-V6", "veg-hf_1kmgrid_v6.Rdata"))
    save(dd1km_nsr, dd1km_nr,
        file=file.path(ROOT, VER, "data/kgrid-V6", "veg-hf_nsr_v6.Rdata"))
    #save(kgrid,
    #    file=file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
}

## fix age 0 in saved files -----------------------------

source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

## 1 km grid
load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

sum(dd1km_pred[[1]][,Target0])
sum(dd1km_pred[[2]][,Target0])
sum(dd1km_pred[[1]])
sum(dd1km_pred[[2]])
dd1km_pred <- fill_in_0ages(dd1km_pred, kgrid$NSRNAME)
sum(dd1km_pred[[1]][,Target0])
sum(dd1km_pred[[2]][,Target0])
sum(dd1km_pred[[1]])
sum(dd1km_pred[[2]])


if (SAVE) {
save(dd1km_pred,
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))
}

## summaries

if (FALSE) {
    x <- as.matrix(groupSums(dd1km_pred[[1]]/10^6, 1, kgrid$NRNAME))
    xx <- x[,"WindGenerationFacility"]
    xxx <- rowSums(x)
    data.frame(Wind=xx,Area=xxx,Perc=100*xx/xxx)

    x <- as.matrix(groupSums(dd1km_pred[[1]]/10^6, 1, kgrid$LUFxNSR))
    write.csv(x, file=file.path(ROOT, VER, "out/kgrid", "current-veghf-area-km2.csv"))
}


## check areas for V6 rule set

## check all fragments, nut just these

source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data/kgrid-V6", "veg-hf_nsr_v6.Rdata"))
Fragment <- c("TreedWetland-Mixedwood", "TreedSwamp-Forest",
    "GraminoidWetland", "ShrubbyWetland", "Muskeg")
aa <- as.matrix(dd1km_nr$v6veg[,Fragment])/10^6

## GraminoidWetland
round(x <- as.matrix(dd1km_nr$v6veg)[,c("GraminoidWetland",
    "Marsh", "GraminoidBog", "GraminoidFen")] / 10^6, 2)
round(100*x[,-1]/rowSums(x[,-1]),2)
find_max(100*x[,-1]/rowSums(x[,-1]))
r1 <- data.frame(round(x[,1,drop=FALSE],2),
    round(100*x[,-1]/rowSums(x[,-1]),2),
    find_max(100*x[,-1]/rowSums(x[,-1])))
print(r1, digits=2)
## regional level solution:
## Marsh in Grassland, GraminoidFen otherwise (one outlier is close in %)

## ShrubbyWetland
round(x <- as.matrix(dd1km_nr$v6veg)[,c("ShrubbyWetland",
    "ShrubbyBog", "ShrubbyFen", "ShrubbySwamp")] / 10^6, 2)
round(100*x[,-1]/rowSums(x[,-1]),2)
find_max(100*x[,-1]/rowSums(x[,-1]))
r2 <- data.frame(round(x[,1,drop=FALSE],2),
    round(100*x[,-1]/rowSums(x[,-1]),2),
    find_max(100*x[,-1]/rowSums(x[,-1])))
print(r2, digits=2)
## table level solution:
## call it ShrubbyFen (not present in Grassland, Alpine is ShrubbySwamp but very little area)

## TreedWetland-Mixedwood
round(x <- as.matrix(dd1km_nr$v6veg)[,c("TreedWetland-Mixedwood",
    "TreedFen-Mixedwood", "TreedBog-BSpr","TreedSwamp-Mixedwood")] / 10^6, 2)
round(100*x[,-1]/rowSums(x[,-1]),2)
find_max(100*x[,-1]/rowSums(x[,-1]))
## call it ShrubbyFen (not present in Grassland, Alpine is ShrubbySwamp but very little area)
r4 <- data.frame(round(x[,1,drop=FALSE],2),
    round(100*x[,-1]/rowSums(x[,-1]),2),
    find_max(100*x[,-1]/rowSums(x[,-1])))
print(r4, digits=2)
## table level solution:
## call it TreedBog-BSpr (Grassland is TreedSwamp-Mix but very little area)

## Muskeg
round(x <- as.matrix(dd1km_nr$v6veg)[,c("Muskeg",
    "GraminoidBog", "GraminoidFen",
    "ShrubbyBog", "ShrubbyFen",
    "TreedBog-BSpr", "TreedFen-BSpr", "TreedFen-Decid", "TreedFen-Larch",
    "TreedFen-Mixedwood")] / 10^6, 2)
round(100*x[,-1]/rowSums(x[,-1]),2)
find_max(100*x[,-1]/rowSums(x[,-1]))
## call it ShrubbyFen (not present in Grassland, Alpine is ShrubbySwamp but very little area)
r3 <- data.frame(x[,1,drop=FALSE], 100*x[,-1]/rowSums(x[,-1]), find_max(100*x[,-1]/rowSums(x[,-1])))
print(r3, digits=2)
aa <- data.frame(round(x[,1,drop=FALSE], 2), find_max(100*x[,-1]/rowSums(x[,-1])))
aa[order(aa[,1]),]
## ??? level solution:
## Shield:ShrubbyFen, Parkland:GraminoidFen, else: TreedBog-BSpr

cn <- c("GraminoidBog", "ShrubbyBog", "TreedBog-BSpr")
round(100*x[,cn]/rowSums(x[,cn]),2)
cn <- c("GraminoidFen", "ShrubbyFen", "TreedFen-BSpr",
    "TreedFen-Decid", "TreedFen-Larch","TreedFen-Mixedwood")
round(100*x[,cn]/rowSums(x[,cn]),2)

