# processing bootstrap results

## ----------------------------------------------
# Bayne et al. 2016
# we don't have this for PIWO
## ----------------------------------------------

library(intrval)
x <- read.csv("~/repos/abmianalytics/projects/osm-oven/data/bayne-2016/condor-suppl.csv", stringsAsFactors=FALSE)
x1 <- x[,1:4]
x1[] <- lapply(x1, as.factor)
rownames(x1) <- x1$AOUcode
x2 <- x[,5:13]
x3 <- x[,14:22]

f1 <- function(i) {
    z <- t(sapply(strsplit(x2[,i], "/"), as.integer))
    colnames(z) <- paste0(colnames(x2)[i], c("_Impact", "_Control"))
    z
}
x2 <- do.call(cbind, lapply(1:ncol(x2), f1))
rownames(x2) <- x1$AOUcode

f21 <- function(z) {
    if (z == "")
        return(c(NA, NA, NA))
    zz <- strsplit(z, " ")
    zzz <- strsplit(zz[[1]][2], "-")
    as.numeric(gsub("[^0-9\\.]", "", c(zz[[1]][1], zzz[[1]])))
}
f2 <- function(i) {
    z <- t(sapply(x3[,i], f21))
    colnames(z) <- paste0(colnames(x3)[i], c("", "_Lo", "_Hi"))
    z
}
x3 <- do.call(cbind, lapply(1:ncol(x3), f2))
rownames(x3) <- x1$AOUcode
xx <- data.frame(x1, x2, x3)

spp <- "BHCO"

m <- matrix(x3[spp,], nrow=3)
m0 <- matrix(m[1,], 3, 3)
mL <- matrix(m[2,], 3, 3)
mU <- matrix(m[3,], 3, 3)
main <- paste(x1[spp, "CommonName"], x1[spp, "Habitat"])
colnames(m0) <- c("Seismic", "Pipeline", "Wellpad")
rownames(m0) <- c("50m", "100m", "Unlimited")
dimnames(mL) <- dimnames(mU) <- dimnames(m0)
m0[is.na(m0)] <- 0
mL[is.na(mL)] <- 0
mU[is.na(mU)] <- 0
M1 <- list(est=m0, lwr=mL, upr=mU)


spp <- "RUGR"

m <- matrix(x3[spp,], nrow=3)
m0 <- matrix(m[1,], 3, 3)
mL <- matrix(m[2,], 3, 3)
mU <- matrix(m[3,], 3, 3)
main <- paste(x1[spp, "CommonName"], x1[spp, "Habitat"])
colnames(m0) <- c("Seismic", "Pipeline", "Wellpad")
rownames(m0) <- c("50m", "100m", "Unlimited")
dimnames(mL) <- dimnames(mU) <- dimnames(m0)
m0[is.na(m0)] <- 0
mL[is.na(mL)] <- 0
mU[is.na(mU)] <- 0
M2 <- list(est=m0, lwr=mL, upr=mU)
dput(M)

MM <- list(BHCO=M1, RUGR=M2)
dput(MM)

v <- barplot(m0,beside=TRUE, main=main, col=col, ylim=c(0, min(10, max(mU))))
abline(h=1, col='#8da0cb')
for (k in 1:9) {
    int <- c(mL[k], mU[k])
    lwd <- if (1 %[]% int) 1 else 3
    lines(c(v[k], v[k]), int, lwd=lwd)
}

## ----------------------------------------------
## Calculate bootstrap for maps and sector effects
## ----------------------------------------------

## define OSR
library(rgdal)
library(sp)
library(mefa4)
library(raster)
library(qs)

load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]

ogrListLayers("d:/spatial/Oilsands-Boundaries.gdb")
pl <- readOGR("d:/spatial/Oilsands-Boundaries.gdb", "OilsandRegionDissolve10TM")
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(pl))
o <- over(xy, pl)
#plot(xy, pch=".", col=ifelse(is.na(o$FIELDCODE), 1, 4))

ss <- !is.na(o$FIELDCODE)
pveg_mine <- 0.44 # proportion of vegetated mines in OSR from NDVI


# get uncertainty maps

SPP <- c("RuffedGrouse", "BrownheadedCowbird", "PileatedWoodpecker")


spp <- SPP[3]

rn <- rownames(kgrid)[ss]

## pop sizes
ROOT <- "s:/AB_data_v2020/Results/pred-boot"
NNN <- list()
for (spp in SPP) {
    gc()
    NNi <- matrix(NA, 2, 100)
    rownames(NNi) <- c("Ref", "Curr")
    for (i in 1:100) {
        cat(spp, i, "\n")
        flush.console()
        qreadm(file.path(ROOT, "birds", spp, paste0(spp, "-", i, ".qrda"))) # Ncr, Nrf
        NNi["Ref", i] <- sum(Nrf[rn,])
        NNi["Curr", i] <- sum(Ncr[rn,])
    }
    NNN[[spp]] <- NNi
}
## these are summaed average densities
## x100 to get number of males, x200 to get pop size
save(NNN, file="s:/AB_data_v2020/Results/osm/NNN.RData")
save(NNN, file="~/GoogleWork/abmi/osm-2021/NNN.RData")


## scaling

x <- 1:1000
y <- ifelse(x <= 100, x, plogis((x-100)/100)*200)
plot(x, y, type="l")

## usual sectors
ROOT <- "s:/AB_data_v2020/Results/pred-boot"

AreaN <- colSums(groupSums(trVeg[rn,], 2, chVeg$sector))
AreaN <- 100 * AreaN/sum(AreaN)
sectors <- names(AreaN)

RES1 <- list()
for (spp in SPP) {
    gc()
    SE <- NULL # sector effects
    #DDcr <- DDrf <- matrix(0, length(rn), 100)
    #rownames(DDcr) <- rownames(DDrf) <- rn
    DDcr <- DDrf <- matrix(NA, nrow(kgrid), 100)
    rownames(DDcr) <- rownames(DDrf) <- rownames(kgrid)
    for (i in 1:100) {
        cat(spp, i, "\n")
        flush.console()

        qreadm(file.path(ROOT, "birds", spp, paste0(spp, "-", i, ".qrda"))) # Ncr, Nrf

        Ref <- colSums(Nrf[rn,])
        RefTotal <- sum(Ref)
        Curr <- colSums(Ncr[rn,])

        Curr[is.na(Curr)] <- 0
        Ref[is.na(Ref)] <- 0
        total.effect <- (100 * (Curr - Ref)/RefTotal)[sectors]
        unit.effect <- 100 * total.effect/AreaN[sectors]
        total.effect[is.na(total.effect)] <- 0
        unit.effect[is.na(unit.effect)] <- 0
        AREAS <- AreaN[sectors]
        underhf.effect <- (100 * (Curr - Ref)/Ref)[sectors]
        underhf.effect[is.na(underhf.effect)] <- 0
        SE <- rbind(SE,
            c(Area=AREAS, Total=total.effect, Unit=unit.effect, UnderHF=underhf.effect))

        #DDcr[,i] <- rowSums(Ncr[rn,])
        #DDrf[,i] <- rowSums(Nrf[rn,])
        DDcr[rownames(Ncr),i] <- rowSums(Ncr)
        DDrf[rownames(Nrf),i] <- rowSums(Nrf)


    }
    RES1[[spp]] <- list(se=SE, cr=DDcr, rf=DDrf)
}


## fine sectors
ROOT <- "s:/AB_data_v2020/Results/pred-boot-fine"

tmp1 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v2020.csv")
chVeg$sector_fine <- tmp1$SectorFine[match(chVeg$cr, tmp1$ID)]

AreaN <- colSums(groupSums(trVeg[rn,], 2, chVeg$sector_fine))
AreaN <- 100 * AreaN/sum(AreaN)
sectors <- names(AreaN)

for (spp in SPP) {
    SE <- NULL # sector effects
    for (i in 1:100) {
        cat(spp, i, "\n")
        flush.console()

        qreadm(file.path(ROOT, "birds", spp, paste0(spp, "-", i, ".qrda"))) # Ncr, Nrf

        Ref <- colSums(Nrf[rn,])
        RefTotal <- sum(Ref)
        Curr <- colSums(Ncr[rn,])

        Curr[is.na(Curr)] <- 0
        Ref[is.na(Ref)] <- 0
        total.effect <- (100 * (Curr - Ref)/RefTotal)[sectors]
        unit.effect <- 100 * total.effect/AreaN[sectors]
        total.effect[is.na(total.effect)] <- 0
        unit.effect[is.na(unit.effect)] <- 0
        AREAS <- AreaN[sectors]
        underhf.effect <- (100 * (Curr - Ref)/Ref)[sectors]
        underhf.effect[is.na(underhf.effect)] <- 0
        SE <- rbind(SE,
            c(Area=AREAS, Total=total.effect, Unit=unit.effect, UnderHF=underhf.effect))

    }
    RES1[[spp]]$fine <- SE
}

fstat <- function(x, level=0.95, ...) {
    .fstat <- function(x, level=0.95, ...)
        c(Mean=mean(x, ...),
            Median=median(x, ...),
            quantile(x, c((1-level)/2, 1 - (1-level)/2), ...),
            First=x[1],
            SD=sd(x, ...),
            IQR=IQR(x, ...))
    if (is.null(dim(x)))
        .fstat(x, level, ...) else t(apply(x, 1, .fstat, level=level, ...))
}

## stats
for (spp in SPP) {
    RES1[[spp]]$Dcr <- fstat(RES1[[spp]]$cr, na.rm=TRUE)
    RES1[[spp]]$Drf <- fstat(RES1[[spp]]$rf, na.rm=TRUE)
}

## intactness
for (spp in SPP) {
    is1 <- 100 * ifelse(RES1[[spp]]$cr > RES1[[spp]]$rf,
        RES1[[spp]]$rf / RES1[[spp]]$cr,
        RES1[[spp]]$cr / RES1[[spp]]$rf)
    is2 <- is1
    z <- RES1[[spp]]$cr > RES1[[spp]]$rf
    z[is.na(z)] <- FALSE
    is2[z] <- 200 - is1[z]

    RES1[[spp]]$is1 <- fstat(is1, na.rm=TRUE)
    RES1[[spp]]$is2 <- fstat(is2, na.rm=TRUE)
}

## delete big tables
for (spp in SPP) {
    RES1[[spp]]$cr <- NULL
    RES1[[spp]]$rf <- NULL
}

## save
save(RES1, file="s:/AB_data_v2020/Results/osm/RES1.RData")
save(RES1, file="~/GoogleWork/abmi/osm-2021/RES1.RData")

