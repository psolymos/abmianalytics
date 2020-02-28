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

save(SPP, AVegSS, ch2veg,
    file=paste0(INDIR, "/sector-tools.RData"))

## summarize species

library(mefa4)
library(intrval)

#INDIR <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/sector")
INDIR <- paste0("~/GoogleWork/tmp/sector")
load(paste0(INDIR, "/sector-tools.RData"))

A <- groupSums(matrix(AVegSS, ncol=1), 1, ch2veg$sector)

load_spp <- function(spp) {
    out <- list()
    for (SEC in c("All", "Agr", "For", "Energy", "Urban","Transp")) {
        e <- new.env()
        load(paste0(INDIR, "/", spp, "/", spp, "-HF-", toupper(SEC), ".RData"), envir=e)
        out[[SEC]] <- list(
            HCR = groupSums(e$HABCR, 1, ch2veg$sector),
            HRF = groupSums(e$HABRF, 1, ch2veg$sector))
    }
    out
}

se_fun <- function(i, H) {
    Nsect <- list()
    for (SEC in names(H)) {
        Nsect[[SEC]] <- cbind(
            cr=H[[SEC]]$HCR[,i],
            rf=H[[SEC]]$HRF[,i])
    }
    Nref <- Nsect[[1]][,2]
    Diffs <- sapply(Nsect, function(z) z[,1]-Nref)
    Diffs <- Diffs[,c("All", "For", "Transp", "Energy", "Urban", "Agr")]
    Df2 <- cbind(Native=0,Diffs[,-1])
    Df2 <- Df2[,c("Native", "For", "Transp", "Energy", "Urban", "Agr")]
    Tots <- data.frame(#All=Diffs[,1],
        Joint=Diffs[rownames(Diffs) != "Misc","All"], # joint
        Direct=diag(Df2), # pbf direct (pbf=partial backfilled)
        Spillover=colSums(Df2)-diag(Df2)) # indirect effects: pbf all - pbf direct
    Tots <- 100*Tots/sum(Nref)
    Tots$Indirect <- Tots[,"Joint"]-Tots[,"Direct"]
    as.matrix(Tots)
}

sector <- function(spp, level=0.95) {
    a <- c((1-level)/2, 1-((1-level)/2))
    H <- load_spp(spp)
    B <- ncol(H[[1]][[1]])
    tmp <- se_fun(1, H)
    S <- array(0, c(dim(tmp), B))
    dimnames(S) <- list(rownames(tmp), colnames(tmp), NULL)
    for (i in seq_len(B))
        S[,,i] <- se_fun(i, H)
    OUT <- list(estimate=tmp, lower=tmp, upper=tmp)
    for (j in seq_len(nrow(tmp))) {
        for (k in seq_len(ncol(tmp))) {
            z <- quantile(S[j,k,], c(0.5, a))
            OUT$estimate[j,k] <- z[1]
            OUT$lower[j,k] <- z[2]
            OUT$upper[j,k] <- z[3]
        }
    }
    OUT
}

sector_res <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    sector_res[[spp]] <- sector(spp)
}
save(sector_res, file=paste0(INDIR, "/sector-res.RData"))

load(paste0(INDIR, "/sector-res.RData"))
library(intrval)

v <- lapply(1:6, function(i)
    t(sapply(sector_res, function(z) z$estimate[i,1:2])))
names(v) <- rownames(sector_res[[1]]$estimate)

pdf(paste0(INDIR, "/sector-all.pdf"), width=12, height=8)
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
dev.off()

A <- groupSums(matrix(AVegSS, ncol=1), 1, ch2veg$sector)
A <- A[,1]/sum(A)

tab_fun <- function(spp, cn, unit=FALSE) {
    z <- sector_res[[spp]]
    if (unit) {
        z[[1]] <- z[[1]]/A[rownames(z[[1]])]
        z[[2]] <- z[[2]]/A[rownames(z[[1]])]
        z[[3]] <- z[[3]]/A[rownames(z[[1]])]
    }
    out <- data.matrix(c(z[[1]][,cn], lower=z[[2]][,cn], upper=z[[3]][,cn]))
    out <- data.frame(spp=spp, cn=cn, unit=unit, as.data.frame(t(out)))
}
Res <- rbind(do.call(rbind, lapply(SPP, tab_fun, "Joint")),
    do.call(rbind, lapply(SPP, tab_fun, "Joint", TRUE)),
    do.call(rbind, lapply(SPP, tab_fun, "Direct")),
    do.call(rbind, lapply(SPP, tab_fun, "Direct", TRUE)))
write.csv(Res, row.names = FALSE, file=paste0(INDIR, "/sector-results.csv"))

plot_sect <- function(spp, unit=FALSE) {

    z <- sector_res[[spp]]
    if (unit) {
        z[[1]] <- z[[1]]/A[rownames(z[[1]])]
        z[[2]] <- z[[2]]/A[rownames(z[[1]])]
        z[[3]] <- z[[3]]/A[rownames(z[[1]])]
    }
    rr1 <- range(z[[1]][,1:2])
    rr2 <- range(z[[1]][,1:2], z[[2]][,1:2], z[[3]][,1:2])
    rr <- rr2
    rr[1] <- if (abs(rr2[1] - rr1[1])/diff(rr1) > 0.5)
        rr1[1]-0.5*diff(rr1) else rr2[1]
    rr[2] <- if (abs(rr2[2] - rr1[2])/diff(rr1) > 0.5)
        rr1[2]+0.5*diff(rr1) else rr2[2]
    #rr <- rr - (rr %% 10)
    #rr <- rr + c(-10, 10)

    YLAB <- if (unit)
        "Unit effect" else "Regional effect"
    plot(0, type="n", xlim=c(0.5, 6.5), ylim=rr, ann=FALSE, axes=FALSE)
    w <- 0.4
    Cols <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')
    for (i in 1:6) {
        COL <- paste0(Cols[i], '66')
        polygon(
            i + c(0, -w, -w, 0),
            c(0, 0, z$estimate[i, "Direct"], z$estimate[i, "Direct"]),
            col=paste0(Cols[i], '44'), border=paste0(Cols[i], 'aa'))
        polygon(
            i + c(0, -w, -w, 0) + w,
            c(z$estimate[i, "Direct"], z$estimate[i, "Direct"],
              z$estimate[i, "Joint"], z$estimate[i, "Joint"]),
            col=paste0(Cols[i], '88'), border=paste0(Cols[i], 'aa'))
        NONSIG <- 0 %[]% c(z$lower[i, "Joint"], z$upper[i, "Joint"])
        lines(
            i + c(w, w)/2,
            c(z$lower[i, "Joint"], z$upper[i, "Joint"]),
            lwd=if (NONSIG) 1.5 else 3,lend=1,
            col=if (NONSIG) paste0(Cols[i], 'aa') else paste0(Cols[i], 'ff'))
        lines(
            i + c(0, w),
            c(z$estimate[i, "Joint"], z$estimate[i, "Joint"]),
            lwd=3, lend=1,
            col=if (NONSIG) paste0(Cols[i], 'aa') else paste0(Cols[i], 'ff'))
    }
    title(ylab=YLAB)
    axis(2, col="grey")
    abline(h=0, col="grey")
    if (unit)
        abline(h=c(-100, 100), lty=2)
    axis(1, 1:6, substr(rownames(z[[1]]),1,1), tick=FALSE)
    invisible(z)
}

pdf(paste0(INDIR, "/sector-results.pdf"), onefile=TRUE)
for (spp in SPP) {
    op <- par(mfrow=c(2,1), las=1, mar=c(4,4,4,1))
    plot_sect(spp)
    title(main=spp)
    plot_sect(spp, TRUE)
    par(op)
}
dev.off()

vvv <- data.frame(sapply(v, function(z) z[,1] - z[,2]))
plot(vvv)

op <- par(mfrow=c(2,3))
for (i in 1:6) {
hist(vvv[,i], main=colnames(vvv)[i], sub="", col="#00000066",
    xlab="Indirect", breaks=20)
abline(v=0, col=3, lty=1, lwd=2)
}
par(op)
