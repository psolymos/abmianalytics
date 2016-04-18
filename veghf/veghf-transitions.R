source("~/repos/abmianalytics/veghf/veghf-setup.R")

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
    col.year=HF_YEAR, col.HFyear="CutYear", wide=FALSE)
ddd0 <- dd[character(0),cc]
xddd0 <- ddd0

for (s in 1:(length(Start)-1)) {
    cat("\n----------------------\n\nStarting block", s, "\n")
    for (i in Start[s]:(Start[s+1]-1)) {
        cat(i, "of", length(fl), "-", fl[i], "\t")
        flush.console()
        d <- read.csv(file.path(ROOT, VER, "data", "kgrid", "tiles", fl[i]))
        dd <- make_vegHF_wide(d, col.label="Row_Col", 
            col.year=HF_YEAR, col.HFyear="CutYear", wide=FALSE)
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
    reclass=c("Conif0", "Decid0", "Mixwood0", "Pine0", "Swamp", "Swamp", 
        "Swamp", "Swamp", "BSpr0", "Larch0", "Larch0"))
    
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

    cat("\nSaving", i, "\n\n")
    flush.console()
    save(trVeg, trSoil, file=file.path(ROOT, VER, "out", "transitions", 
        paste0(i, ".Rdata")))
}
