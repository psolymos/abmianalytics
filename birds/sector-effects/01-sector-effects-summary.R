## Sector effects

## check partial backill classes to match sectors used here!!!


library(mefa4)
library(intrval)
library(raster)
library(rgdal)
library(sp)

INDIR <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/sector")

sectors_all <- c("Agr", "Energy", "For", "Urban","Transp")

load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y


## ch2soil ch2veg trSoil trVeg
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
tv$ao <- as.factor(paste0(as.character(tv[, "UseInAnalysisFine"]), ifelse(tv[, "MatureOld"], "O", "")))

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2veg$rf2 <- tv$UseInAnalysisFineAge[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$UseInAnalysisFineAge[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]
ch2veg$rf3 <- tv$ao[match(ch2veg$rf, rownames(tv))]
ch2veg$cr3 <- tv$ao[match(ch2veg$cr, rownames(tv))]

with(ch2veg, table(cr3, sector))

EXCL <- c("HWater", "SoilUnknown", "SoilWater", "Water")
MODIF <- c("SoftLin", "Well", "EnSoftLin", "TrSoftLin", "Seismic")

keepn <- rownames(ch2veg)[!(ch2veg$cr2 %in% EXCL) & !(ch2veg$rf2 %in% EXCL)]
trVeg <- trVeg[,keepn]
ch2veg <- ch2veg[keepn,]
ch2veg$modif <- ch2veg$cr2 %in% MODIF
rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn


## define OSR
ogrListLayers("d:/spatial/Oilsands-Boundaries.gdb")
pl <- readOGR("d:/spatial/Oilsands-Boundaries.gdb", "OilsandRegionDissolve10TM")
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(pl))
o <- over(xy, pl)
#plot(xy, pch=".", col=ifelse(is.na(o$FIELDCODE), 1, 4))

ss <- !is.na(o$FIELDCODE)
trVegSS <- trVeg[ss,]
AVegSS <- colSums(trVegSS)

levs <- c(
    "MineSite",
    "Pipeline",
    "TransmissionLine",
    "RailHardSurface",
    "RailVegetatedVerge",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge",
    "SeismicLineNarrow",
    "SeismicLineWide",
    "IndustrialSiteRural",
    "UrbanIndustrial",
    "WellSite")

sect <- as.character(ch2veg$sector)
sect[ch2veg$cr %in% levs] <- as.character(ch2veg$cr)[ch2veg$cr %in% levs]
sect[sect %in% c("SeismicLineNarrow", "SeismicLineWide")] <- "SeismicLine"
sect[sect %in% c("Pipeline", "TransmissionLine")] <- "PipeTransLine"
sect[sect %in% c("RailHardSurface",
    "RailVegetatedVerge",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge")] <- "RoadRailVerge"
sect[sect %in% c("IndustrialSiteRural", "UrbanIndustrial")] <- "Industrial"
sect[sect == "Energy"] <- "OtherEnergy"
data.frame(x=table(sect))


SPP <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATTW", "BAOR", "BARS",
    "BAWW", "BBMA", "BBWA", "BBWO", "BCCH", "BHCO", "BHVI", "BLJA",
    "BLPW", "BOCH", "BRBL", "BRCR", "BTNW", "CAWA", "CCSP", "CEDW",
    "CHSP", "CMWA", "CONW", "CORA", "COYE", "DEJU", "DOWO", "EAKI",
    "EVGR", "FOSP", "GCKI", "GRAJ", "GRCA", "HAWO", "HETH", "HOLA",
    "HOWR", "LCSP", "LEFL", "LISP", "MAWA", "MODO", "MOWA", "NOFL",
    "NOWA", "OCWA", "OSFL", "OVEN", "PAWA", "PHVI", "PISI", "PIWO",
    "PUFI", "RBGR", "RBNU", "RCKI", "REVI", "RUBL", "RWBL", "SAVS",
    "SOSP", "SWSP", "SWTH", "TEWA", "TRES", "VATH", "VEER", "VESP",
    "WAVI", "WBNU", "WCSP", "WETA", "WEWP", "WIWA", "WIWR", "WTSP",
    "WWCR", "YBFL", "YBSA", "YEWA", "YRWA")


sector_res <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()

    ## !!! evaluate $range !!!

    Nsect <- list()
    for (SEC in c("All", sectors_all)) {
        e <- new.env()
        load(paste0(INDIR, "/", spp, "/", spp, "-HF-", toupper(SEC), ".RData"), envir=e)

        HCR <- groupSums(e$HABCR, 1, ch2veg$sector)
        HRF <- groupSums(e$HABRF, 1, ch2veg$sector)
        A <- groupSums(matrix(AVegSS, ncol=1), 1, ch2veg$sector)

        Nsect[[SEC]] <- cbind(cr=apply(HCR, 1, median), rf=apply(HRF, 1, median), a=A[,1])
    }
#    e <- new.env()
#    load(paste0(INDIR, "/", spp, "/", spp, "-SPACE-ALL.RData"), envir=e)
#    HCR <- groupSums(e$HABCR, 1, ch2veg$sector)
#    HRF <- groupSums(e$HABRF, 1, ch2veg$sector)
#    A <- groupSums(matrix(AVegSS, ncol=1), 1, ch2veg$sector)
#    NsectSp <- cbind(cr=apply(HCR, 1, median), rf=apply(HRF, 1, median), a=A[,1])

    ## reference abundance should not differ
    #round(100*sapply(Nsect, function(z) (z[,2]-Nsect[[1]][,2])/sum(Nsect[[1]][,2])), 3)
    #max(abs(round(100*sapply(Nsect, function(z) (z[,2]-Nsect[[1]][,2])/sum(Nsect[[1]][,2])), 3)))

    Nref <- Nsect[[1]][,2]
    Diffs <- sapply(Nsect, function(z) z[,1]-Nref)
    Diffs <- Diffs[,c("All", "For", "Transp", "Energy", "Urban", "Agr")]
    Ahf <- Nsect[[1]][,3]
    Df2 <- cbind(Native=0,Diffs[,-1])
    Df2 <- Df2[,c("Native", "For", "Transp", "Energy", "Urban", "Agr")]

    Tots <- cbind(#All=Diffs[,1],
        Joint=Diffs[rownames(Diffs) != "Misc","All"], # joint
        Direct=diag(Df2), # pbf direct (pbf=partial backfilled)
        Indirect=colSums(Df2)-diag(Df2)) # indirect effects: pbf all - pbf direct
    Tots <- 100*Tots/sum(Nref)
    U <- (Tots/100) / (Ahf[rownames(Tots)]/sum(Ahf[rownames(Tots)]))
    attr(Tots, "Nref") <- sum(Nref)
    attr(Tots, "synergy") <- 100*(sum(Df2)-sum(Diffs[,1]))/sum(Nref) # (isolated-joint)/Nref

#    round(Tots, 2)
#    100*Ahf/sum(Ahf)

    sector_res[[spp]] <- list(diff=Diffs, total=Tots, unit=U, ref=Nref, ahf=Ahf)

}

library(intrval)

v <- lapply(1:6, function(i) t(sapply(sector_res, function(z) z$total[i,1:2])))
names(v) <- rownames(sector_res[[1]]$total)

op <- par(mfrow=c(2,3))
hist(v[[1]][,1], main="Native", sub="", breaks=20, xlim=c(-50,250),
    xlab="Indirect", col="#00000066")
abline(v=0, col=3, lty=1, lwd=2)
for (i in 2:6) {
    vv <- v[[i]]
    vv <- vv[order(abs(vv[,1] - vv[,2])),]
    lim <- c(-10, 25)
    if (i == 6)
        lim <- c(-50, 300)
    col <- ifelse(vv[,1] > vv[,2], "#ff000066", "#0000ff66")
    col[abs(vv[,1] - vv[,2]) %[]% c(0, 1)] <- "#00000066"
    plot(vv, col=col, main=names(v)[i],
        xlim=lim, ylim=lim, pch=19, cex=1.5)
    abline(0,1, col="grey")
    abline(h=0,v=0, col=3, lty=1)
}
par(op)

vvv <- data.frame(sapply(v, function(z) z[,1] - z[,2]))
plot(vvv)

op <- par(mfrow=c(2,3))
for (i in 1:6) {
hist(vvv[,i], main=colnames(vvv)[i], sub="", col="#00000066",
    xlab="Indirect", breaks=20)
abline(v=0, col=3, lty=1, lwd=2)
}
par(op)
