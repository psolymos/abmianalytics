source("~/repos/abmianalytics/veghf/veghf-setup.R")

### 1K grid --------------------------------------------------------

## Sample year is current year, so that forest ages are relative to present
## and not relative to HF or veg inventory year.

fl <- list.files(file.path(ROOT, VER, "data", "kgrid", "tiles"))

## test feature types
if (FALSE) {
NEW <- character(0)
for (fn in fl) {
    cat("checking", which(fl == fn), "/", length(fl));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
    d <- read.csv(f)
    diff <- setdiff(levels(d$FEATURE_TY), rownames(hflt))
    if (length(diff))
        NEW <- union(NEW, diff)
    cat("\t", length(NEW), "new types found\n")
}

## tracking the strange area mismatch
natrack <- list()
for (fn in fl) {
    cat("checking", which(fl == fn), "/", length(fl), "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
    d <- read.csv(f)
    dd <- make_vegHF_wide(d, col.label="Row_Col", 
        col.year=NULL, col.HFyear="CutYear", wide=FALSE)
    tmp <- colSums(is.na(dd[,c("VEGAGEclass",
        "VEGHFAGEclass","SOILclass","SOILHFclass")]))
    natrack[[fn]] <- tmp
    if (tmp[1] > 0)
        break
}

## check blank HABIT cases
blank_n <- numeric(length(fl)) # no. of cases
blank_a <- numeric(length(fl)) # total area
for (i in 1:length(fl)) {
    cat("\nchecking", i, "/", length(fl));flush.console()
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fl[i])
    d <- read.csv(f)
    tmp <- d[d$HABIT == "",,drop=FALSE]
    if (nrow(tmp)) {
        blank_n[i] <- nrow(tmp)
        blank_a[i] <- sum(tmp$Shape_Area)
        cat("\tfound:", nrow(tmp))
    }
}

(blanks <- which(blank_n > 0))
sum(blank_a) # 1.262801 < 2 m^2
blank_a[blanks]

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
    f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
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
        f <- file.path(ROOT, VER, "data", "kgrid", "tiles", fn)
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
    file.path(ROOT, VER, "data", "kgrid", 
    "Grid1km_template_final_clippedBy_ABBound_with_atts_to_Peter.csv"))
rownames(kgrid) <- kgrid$Row_Col
## this is the correct lat/long (i.e. more decimal places)
kgrid2 <- read.csv(
    file.path(ROOT, VER, "data", "kgrid", 
    "Grid1km_template_latLong.csv"))
rownames(kgrid2) <- kgrid2$Row_Col
kgrid2 <- kgrid2[rownames(kgrid),]
stopifnot(all(rownames(kgrid) == rownames(kgrid2)))
kgrid$POINT_X <- kgrid2$POINT_X
kgrid$POINT_Y <- kgrid2$POINT_Y
rm(kgrid2)

## this is has Brandt boreal and BCR
kgrid2 <- read.csv(
    file.path(ROOT, VER, "data", "kgrid", 
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
    file.path(ROOT, VER, "data", "kgrid", 
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


if (SAVE) {
save(dd1km_pred, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
save(kgrid,
    file=file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
}

## fix age 0 in saved files -----------------------------

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

save(dd1km_pred, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))


### Transition for 1K grid ------------------------------------------------

## this is based on the fix-fire fix-age0 version
## label collapsing as desired (swamp/wetland, ages?)

source("~/repos/abmianalytics/R/veghf_functions.R")

load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
fl <- list.files(file.path(ROOT, VER, "data", "kgrid", "tiles"))

cc <- c("Row_Col","VEGAGEclass","VEGHFAGEclass","SOILclass","SOILHFclass","Shape_Area")

Start <- c(0:79*10+1, 802)


d <- read.csv(file.path(ROOT, VER, "data", "kgrid", "tiles", fl[1]))
dd <- make_vegHF_wide(d, col.label="Row_Col", 
    col.year=NULL, col.HFyear="CutYear", wide=FALSE)
ddd0 <- dd[character(0),cc]
xddd0 <- ddd0

for (s in 1:(length(Start)-1)) {
    cat("----------------------\nStarting block", s, "\n")
    for (i in Start[s]:(Start[s+1]-1)) {
        cat(i, "of", length(fl), "-", fl[i], "\t")
        flush.console()
        d <- read.csv(file.path(ROOT, VER, "data", "kgrid", "tiles", fl[i]))
        dd <- make_vegHF_wide(d, col.label="Row_Col", 
            col.year=NULL, col.HFyear="CutYear", wide=FALSE)
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
            "data", "kgrid", "long", paste0("Long-part", i, ".Rdata")))
        ddd0 <- xddd0
        gc()
    }
}

## -- works on pre-saved chunks

load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))
lu <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
su <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

lu$use_tr <- as.character(lu$VEGAGE_use)
lu$use_tr[!is.na(lu$HF)] <- as.character(lu$VEGHFAGE[!is.na(lu$HF)])
lu$use_tr[lu$use_tr == "WetBare"] <- "NonVeg"
allVegTr <- unique(c(lu$use_tr[is.na(lu$HF)], 
    paste0(rep(lu$use_tr[is.na(lu$HF)], sum(!is.na(lu$HF))), "->",
    rep(lu$use_tr[!is.na(lu$HF)], each=sum(is.na(lu$HF))))))
lu$use_tr <- as.factor(lu$use_tr)

su$use_tr <- as.character(su$Levels1)
su$use_tr[!is.na(su$HF)] <- as.character(su$SOILHF[!is.na(su$HF)])
allSoilTr <- unique(c(su$use_tr[is.na(su$HF)], 
    paste0(rep(su$use_tr[is.na(su$HF)], sum(!is.na(su$HF))), "->",
    rep(su$use_tr[!is.na(su$HF)], each=sum(is.na(su$HF))))))
su$use_tr <- as.factor(su$use_tr)

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))
Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0", 
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0", 
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")

recl <- list(
    bf=c("Conif", "Decid", "Mixwood", "Pine", "Swamp-Conif", "Swamp-Decid", 
        "Swamp-Mixwood", "Swamp-Pine", "Wetland-BSpr", "Wetland-Decid", "Wetland-Larch"),
    target=c("Conif0", "Decid0", "Mixwood0", "Pine0", "Swamp-Conif0", "Swamp-Decid0", 
        "Swamp-Mixwood0", "Swamp-Pine0", "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0"),
    reclass=c("Conif0", "Decid0", "Mixwood0", "Pine0", "Conif0", "Decid0", 
        "Mixwood0", "Pine0", "BSpr0", "Decid0", "Larch0"))
    
fl3 <- list.files(file.path(ROOT, VER, "data", "kgrid", "long"))


## do one LUFxNSR class at a time and save it
#i <- "UpperAthabasca_CentralMixedwood"
for (ii in 1:nlevels(kgrid$LUFxNSR)) {
    i <- levels(kgrid$LUFxNSR)[ii]
    cat("\n---------", i)
    #j <- 4
    units <- list()
    sunits <- list()
    for (j in 1:length(fl3)) {
        cat("\n", j);flush.console()
        load(file.path(ROOT, VER, "data", "kgrid", "long", fl3[j]))
        flush.console()
        ddd0$LUFxNSR <- kgrid$LUFxNSR[match(ddd0$Row_Col, kgrid$Row_Col)]
        nsr <- as.character(kgrid[which(kgrid$LUFxNSR==i)[1], "NSRNAME"])
        
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
            xxx0 <- xxx[xxx$vegTr %in% Target0,,drop=FALSE]
            if (nrow(xxx0)>0) {
                cat("age0")
                xxx1 <- xxx[!(xxx$vegTr %in% Target0),,drop=FALSE]
                xxx0$vegTr <- as.character(xxx0$vegTr)
                xxx0$veg <- sapply(strsplit(as.character(xxx0$vegTr), "->"), "[[", 1)
                xxx0$vhf <- sapply(strsplit(as.character(xxx0$vegTr), "->"), 
                    function(z) z[length(z)])
                xxx0$vhf[xxx0$vhf == xxx0$veg] <- ""

                ## needs to sum to 1, include availability
                ages <- AvgAges$reference[,,nsr]
                areas <- AvgAges$area_rf[nsr,]
                bf0 <- groupMeans(ages * areas, 1, recl$reclass)[,-1]
                bf0 <- bf0 / rowSums(bf0)

                tmp <- list()
                for (k in 1:10) {
                    tmpv <- xxx0
                    target <- substr(tmpv$veg, 1, nchar(tmpv$veg)-1)
                    tmpv$Shape_Area <- tmpv$Shape_Area * bf0[match(tmpv$veg, rownames(bf0)),k]
                    tmpv$veg <- paste0(target, colnames(bf0)[k])
                    tmpv$vegTr <- ifelse(tmpv$vhf == "", tmpv$veg,
                        paste0(tmpv$veg, "->", tmpv$vhf))
                    tmp[[k]] <- tmpv[,colnames(xxx1)]
                }
                xxx0v <- do.call(rbind, tmp)
                xxx <- rbind(xxx1, xxx0v)
            }
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
    
    save(trVeg, trSoil, file=file.path(ROOT, VER, "out", "transitions", 
        paste0(i, ".Rdata")))
}

