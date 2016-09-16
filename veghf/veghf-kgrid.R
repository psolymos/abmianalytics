source("~/repos/abmianalytics/veghf/veghf-setup.R")

### 1K grid --------------------------------------------------------

## Sample year is current year, so that forest ages are relative to present
## and not relative to HF or veg inventory year.

fl <- list.files(file.path(ROOT, VER, "data", "kgrid-V6", "tiles"))

## test feature types
if (FALSE) {
NEW <- character(0)
HF <- character(0)
VEG <- character(0)
SOIL <- character(0)
blank_n <- numeric(length(fl)) # no. of cases
blank_a <- numeric(length(fl)) # total area
for (i in seq_len(length(fl))) {
    fn <- fl[i]
    cat("\nchecking", i, "/", length(fl));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
    d <- read.csv(f)
    diff <- setdiff(levels(d$FEATURE_TY), c("", rownames(hflt)))
    if (length(diff)) {
        NEW <- union(NEW, diff)
        cat("\n\t", length(NEW), "new types found\n")
    }
    HF <- union(HF, levels(d$FEATURE_TY))
    VEG <- union(VEG, levels(d$Veg_Type))
    SOIL <- union(SOIL, levels(d$Soil_Type_1))
    tmp <- d[d$HABIT == "",,drop=FALSE]
    if (nrow(tmp)) {
        blank_n[i] <- nrow(tmp)
        blank_a[i] <- sum(tmp$Shape_Area)
        cat("\n\tfound blanks:", nrow(tmp), "\n")
    }
    cat("\n\tveg:", length(VEG), "- soil:", length(SOIL), "- hf:", length(HF), "\n")
    save(d, file=file.path(ROOT, VER, "data", "kgrid-V6", "tiles-rdata",
        gsub(".csv", ".Rdata", fn, fixed = TRUE)))
}

## no blanks
(blanks <- which(blank_n > 0))
sum(blank_a)
blank_a[blanks]
## no NEW HF labels
NEW
x <- data.frame(Veg_Type=sort(VEG))
write.csv(x, row.names=FALSE, file=file.path(ROOT, VER, "data", "kgrid-V6", "veg-V6.csv"))

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
    f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
    d <- read.csv(f)
    ## HF year is used as base year for prediction purposes
    dd <- make_vegHF_wide(d, col.label="Row_Col",
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
        f <- file.path(ROOT, VER, "data", "kgrid-V6", "tiles", fn)
        d <- read.csv(f)
        dd <- make_vegHF_wide(d, col.label="Row_Col",
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
        scale = "1 km x 1 km prediction grid cells")
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
    scale = "1 km x 1 km prediction grid cells")

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

dd1km_pred$veg_current <- dd1km_pred$veg_current[rownames(kgrid),]
dd1km_pred$veg_reference <- dd1km_pred$veg_reference[rownames(kgrid),]
dd1km_pred$soil_current <- dd1km_pred$soil_current[rownames(kgrid),]
dd1km_pred$soil_reference <- dd1km_pred$soil_reference[rownames(kgrid),]

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

if (SAVE) { ## needed for recalculating average ages
    save(dd1km_pred,
        file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
    save(kgrid,
        file=file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
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
