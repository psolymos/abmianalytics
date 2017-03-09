HF_VERSION <- "2014"

source("~/repos/abmianalytics/veghf/veghf-setup.R")

SCALE <- "km" # km or qs

### Provincial grid (km or QS scale) ---------------------------------------

## Sample year is current year, so that forest ages are relative to present
## and not relative to HF or veg inventory year.

## save in Rdata format
fl0 <- list.files(file.path(ROOT, VER, "data", "raw", "veghf", SCALE))
for (i in seq_len(length(fl0))) {
    fn <- fl0[i]
    cat("\nchecking", i, "/", length(fl0));flush.console()
    f <- file.path(ROOT, VER, "data", "raw", "veghf", SCALE, fn)
    #d <- fread(f)
    d <- read.csv(f)
    save(d, file=file.path(ROOT, VER, "data", "inter", "veghf", SCALE,
        gsub(".csv", ".Rdata", fn, fixed = TRUE)))
}


## check Rdata files

fl0 <- list.files(file.path(ROOT, VER, "data", "inter", "veghf", SCALE))
SLIVER <- numeric(length(fl0))
HF <- character(0)
VEG <- character(0)
COMB <- character(0)
VEGissue <- character(0)
c5 <- character(0)
A1 <- numeric(length(fl0))
A2 <- numeric(length(fl0))
Acomb <- numeric(length(fl0))
Awetmix <- Awetmix0 <- numeric(length(fl0))
Awm <- matrix(NA, length(fl0), 3)
for (i in seq_len(length(fl0))) {
    fn <- fl0[i]
    cat("checking", i, "/", length(fl0), "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "inter", "veghf", SCALE, fn)
    load(f)

    d$c6 <- factor(as.character(d$Combined_ChgByCWCS), VEG_LEVS)

#    Awetmix0[i] <- sum(d$Shape_Area[d$c6 == "TreedWetland-Mixedwood"])
#    d$c6[d$c6 == "TreedWetland-Mixedwood" & d$CWCS == "Fen"] <- "TreedFen-Mixedwood"
#    d$c6[d$c6 == "TreedWetland-Mixedwood" & d$CWCS == "Bog"] <- "TreedBog-BSpr"
#    d$c6[d$c6 == "TreedWetland-Mixedwood" & d$CWCS == "Swamp"] <- "TreedSwamp-Mixedwood"
#    Awetmix[i] <- sum(d$Shape_Area[d$c6 == "TreedWetland-Mixedwood"])
    Awm[i,1] <- sum(d$Shape_Area[d$c6 == "TreedFen-Mixedwood"])
    Awm[i,2] <- sum(d$Shape_Area[d$c6 == "TreedBog-Mixedwood"])
    Awm[i,3] <- sum(d$Shape_Area[d$c6 == "TreedSwamp-Mixedwood"])

    d2 <- c4_fun(d)
    HF <- sort(union(HF, levels(d$FEATURE_TY)))
    VEG <- sort(union(VEG, levels(d2$c4)))
    COMB <- sort(union(COMB, as.character(d$Combined_ChgByCWCS)))
    A1[i] <- sum(d$Shape_Area)/10^6
    A2[i] <- sum(d2$Shape_Area)/10^6
    SLIVER[i] <- sum(d2$Shape_Area[d2$issue > 0])
    Acomb[i] <- sum(d$Shape_Area[d$Combined_ChgByCWCS == ""])
    VEGissue <- sort(union(VEGissue, as.character(d2$VEG3)[d2$issue > 0]))
    d2$c5 <- interaction(d$Veg_Type, d$Moisture_Reg, d$PreBackfill_Source,
        d$CWCS_Class, drop=TRUE, sep="::")
    c5 <- sort(union(c5, levels(d2$c5)))
    #save(d2, file=file.path(ROOT, VER, "data", "kgrid-V6dec", "tiles-rdata-x", fn))
}

setdiff(HF, c("", rownames(hflt)))
setdiff(VEG, levels(recl$Combined))
setdiff(COMB, levels(recl$Combined))
setdiff(levels(recl$Combined), COMB)
summary(A1-A2)
sum(SLIVER,na.rm=TRUE)/10^6
VEGissue
sum(Awetmix0,na.rm=TRUE)/10^6
sum(Awetmix,na.rm=TRUE)/10^6


tmp <- strsplit(c5, "::")
tmp <- do.call(rbind, lapply(tmp, function(z) if (length(z) < 4) c(z, "") else z))
colnames(tmp) <- c("Veg_Type", "Moisture_Reg", "PreBackfill_Source", "CWCS_Class")
d5 <- data.frame(Combine4=c5, tmp)
d5$Combine3 <- interaction(d5$Veg_Type, d5$Moisture_Reg, d5$PreBackfill_Source,
    drop=TRUE, sep="::")
d5$Final <- recl$Combined[match(d5$Combine3, recl$Veg_Moist_preBkfSrs)]
d5$needCWCS <- recl$need_CWCS[match(d5$Combine3, recl$Veg_Moist_preBkfSrs)]

## processing csv files in batches of 50

Start <- c(1, 51, 101, 151, 201, 251, 301, 351, 401, 451,
    501, 551, 601, 651, 701, 751, 802)
fl <- list.files(file.path(ROOT, VER, "data", "inter", "veghf", SCALE))

SCALE_NOTE <- if (SCALE == "km")
    "1 km x 1 km prediction grid cells" else "QS prediction grid cells"
tmplist <- list()
for (s in 1:(length(Start)-1)) {

    gc()
    fn <- fl[Start[s]]
    cat("\n\n------------- batch", s, "----------------\n")
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "inter", "veghf", SCALE, fn)
    e <- new.env()
    load(f, envir=e)
    d <- e$d
    ## HF year is used as base year for prediction purposes
    dd <- make_vegHF_wide_v6(d, col.label="Row_Col",
        col.year=HF_YEAR, col.HFyear="YEAR_MineCFOAgCutblock", sparse=TRUE)
    stopifnot(sum(duplicated(colnames(dd$veg_current))) < 1)
    veg_current <- dd$veg_current
    veg_reference <- dd$veg_reference
    soil_current <- dd$soil_current
    soil_reference <- dd$soil_reference
    sample_year <- dd$sample_year[1]

#lapply(dd[1:4], sum)

    for (i in (Start[s]+1):(Start[s+1]-1)) {

        fn <- fl[i]
        cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
        f <- file.path(ROOT, VER, "data", "inter", "veghf", SCALE, fn)
        e <- new.env()
        load(f, envir=e)
        d <- e$d
        dd <- make_vegHF_wide_v6(d, col.label="Row_Col",
            col.year=HF_YEAR, col.HFyear="YEAR_MineCFOAgCutblock", sparse=TRUE)
        veg_current <- bind_fun2(veg_current, dd$veg_current)
        veg_reference <- bind_fun2(veg_reference, dd$veg_reference)
        soil_current <- bind_fun2(soil_current, dd$soil_current)
        soil_reference <- bind_fun2(soil_reference, dd$soil_reference)

        cs <- colSums(veg_current)
        print(round(sum(cs[grep("Cultivation", names(cs))])/10^6,2))

    }
    tmplist[[s]] <- list(
        veg_current = veg_current,
        veg_reference = veg_reference,
        soil_current = soil_current,
        soil_reference = soil_reference,
        sample_year = sample_year,
        scale = SCALE_NOTE)
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
}

## assembling return object
dd1km_pred <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year = tmplist[[1]]$sample_year,
    scale = SCALE_NOTE)

df <- data.frame(h=colnames(dd1km_pred$veg_current), A=colSums(dd1km_pred$veg_current)/10^6)
df$P <- round(100*df$A / sum(df$A), 2)
rownames(df) <- NULL
df[,-2]
write.csv(df[,c(1,3,2)], file="veg6.csv")


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


load(file.path(ROOT, VER, "out", "kgrid", "kgrid_table.Rdata")) # kgrid
dd1km_pred$veg_current <- dd1km_pred$veg_current[rownames(kgrid),]
dd1km_pred$veg_reference <- dd1km_pred$veg_reference[rownames(kgrid),]
dd1km_pred$soil_current <- dd1km_pred$soil_current[rownames(kgrid),]
dd1km_pred$soil_reference <- dd1km_pred$soil_reference[rownames(kgrid),]
dd1km_pred$v6veg <- dd1km_pred$v6veg[rownames(kgrid),]

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
        file=file.path(ROOT, VER, "data/kgrid-V6dec", "veg-hf_1kmgrid_v6.Rdata"))
    save(dd1km_nsr, dd1km_nr,
        file=file.path(ROOT, VER, "data/kgrid-V6dec", "veg-hf_nsr_v6.Rdata"))
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

