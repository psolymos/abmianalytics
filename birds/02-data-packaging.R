#' ---
#' title: "Bird data packaging"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Nov 28, 2018"
#' output: pdf_document
#' ---
#'
#' # Preamble
#'
#' In this script we filter, transform, mutate, and get some offsets calculated.
#'
library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")
knitr::opts_chunk$set(eval=FALSE)
load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData")


## dealing with sine NAs
## checking veg+soil+hf summaries
dd$vshf <- ifelse(rowSums(is.na(vc1)) > 0, NA, 0)
(aa <- data.frame(n_of_NAs=colSums(is.na(dd))))
## ROAD: needs HF info
## DATE/DATI: use constant sra
## XY: should probably be dropped too
dd <- droplevels(dd[!is.na(dd$X) & !is.na(dd$vshf),])
dd$MAXDUR[is.na(dd$MAXDUR)] <- 1 # ABMISM bits
## climate: probably outside of AB bound: can use nearest if distance is small
dd <- droplevels(dd[!is.na(dd$pAspen) & !is.na(dd$NRNAME),])

## shrink yy
yy <- yy[rownames(dd), colSums(yy > 0) >= 20]

dd <- data.frame(dd, vs0[rownames(dd),])

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
stopifnot(all(colnames(vc1) == rownames(tv)))
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]
stopifnot(all(colnames(sc1) == rownames(ts)))

dd$pWater <- rowSums(row_std(vc1[rownames(dd),])[,tv$is_water])
dd$pRoad <- rowSums(row_std(vc1[rownames(dd),])[,tv$is_road])
dd$pWell <- row_std(vc1[rownames(dd),])[,"WellSite"]
dd$pRoadVeg <- rowSums(row_std(vc1[rownames(dd),])[,tv$is_road_veg])

lcc4 <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"LCC4"]))
tmp <- find_max(lcc4)
dd$LCC4 <- tmp$index
dd$LCC2 <- dd$LCC4
levels(dd$LCC2) <- c("OpenWet", "Forest", "Forest", "OpenWet")
dd$pClosed <- rowSums(row_std(vc1[rownames(dd),])[,tv$is_closed])
dd$pHarest <- rowSums(row_std(vc1[rownames(dd),])[,tv$is_harvest])
dd$pHF <- rowSums(row_std(vc1[rownames(dd),])[,tv$is_HF])


dd$pWater_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv$is_water])
dd$pWet_KM <- rowSums(row_std(vc2[rownames(dd),])[,tv$is_wet])
dd$pWetWater_KM <-dd$pWater_KM + dd$pWet_KM


vc1r <- row_std(groupSums(vc1[rownames(dd),], 2, tv[colnames(vc1),"UseInAnalysisFine"]))
vr1r <- row_std(groupSums(vr1[rownames(dd),], 2, tv[colnames(vr1),"UseInAnalysisFine"]))
sc1r <- row_std(groupSums(sc1[rownames(dd),], 2, ts[colnames(sc1),"UseInAnalysisCoarse"]))
sr1r <- row_std(groupSums(sr1[rownames(dd),], 2, ts[colnames(sr1),"UseInAnalysisCoarse"]))


dd$useSouth <- FALSE
dd$useSouth[dd$NRNAME %in% c("Grassland", "Parkland")] <- TRUE
dd$useSouth[dd$NSRNAME %in% c("Dry Mixedwood")] <- TRUE
dd$useSouth[dd$useSouth & dd$Y > 56.7] <- FALSE # ~80 points
dd$useSouth[dd$useSouth & dd$pWater > 0.5] <- FALSE
dd$useSouth[dd$useSouth & sr1r[,"SoilUnknown"] > 0] <- FALSE
dd$useSouth[dd$useSouth & sc1r[,"HFor"] > 0] <- FALSE

## use in north
dd$useNorth <- TRUE
dd$useNorth[dd$NRNAME == "Grassland"] <- FALSE
dd$useNorth[dd$useNorth & dd$pWater > 0.5] <- FALSE


## road
dd$ROAD[is.na(dd$ROAD) & dd$pRoad > 0.04] <- 1
dd$ROAD[is.na(dd$ROAD)] <- 0
table(dd$ROAD, dd$pRoad > 0.04, useNA="a")

## wells
dd$WELL <- ifelse(dd$pWell > 0.1, 1, 0) # >3/4 of a 1 ha well site is inside)

## soil classification

## reclass and find max
tmp <- find_max(sc1r[,colnames(sc1r) %ni% c("SoilWater", "SoilUnknown", "HWater",
    "SoftLin", "HardLin", "HFor", "Well")])
dd$soilc <- tmp$index
dd$soilv <- tmp$value

tmp <- find_max(vc1r[,colnames(vc1r) %ni% c("Water","HWater", "OtherDisturbedVegetation",
    "Seismic", "HardLin", "TrSoftLin", "SnowIce", "Bare", "Well")])
dd$vegc <- tmp$index
dd$vegv <- tmp$value

#dd$WELL[dd$soilc == "UrbInd"] <- 0
with(dd[dd$useSouth,], table(soilc, WELL))
with(dd[dd$useNorth,], table(vegc, WELL))



ys <- groupSums(yy[dd$useSouth,], 1, dd$SS[dd$useSouth])
ys[ys > 0] <- 1
ys <- ys[,colSums(ys) >= 20]
ds <- nonDuplicated(dd, SS, TRUE)[rownames(ys),]

table(ds$soilc, ds$WELL)
hist(ds$soilv)

library(opticut)
library(parallel)
cl <- makeCluster(8)

o <- opticut(as.matrix(ys) ~ ROAD, data=ds, strata=ds$soilc, dist="binomial:cloglog",
    weights=ds$soilv, cl=cl)

stopCluster(cl)

## derive all the variables and run non BB model
## check if habCl is needed
## get offsets (lcc and tree!)


yn <- groupSums(yy[dd$useNorth,], 1, dd$SS[dd$useNorth])
yn[yn > 0] <- 1
yn <- yn[,colSums(yn) >= 100]
dn <- nonDuplicated(dd, SS, TRUE)[rownames(yn),]

data.frame(table(dn$vegc))
hist(dn$vegv)

hist(vc1r[vc1r[,"Well"]>0.1,"Well"]*7)

library(opticut)
library(parallel)
cl <- makeCluster(8)

o <- opticut(as.matrix(yn) ~ ROAD, data=dn, strata=dn$vegc, dist="binomial:cloglog",
    weights=dn$vegv, cl=cl)

stopCluster(cl)


## JDAY
dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
dd$JDAY <- dd$JULIAN / 365
## TSSR
library(maptools)
Coor <- as.matrix(dd[,c("X", "Y")])
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dd$srise <- NA
dd$srise[subset] <- sr
dd$start <- as.POSIXlt(dd$DATE)$hour + as.POSIXlt(dd$DATE)$min / 60
dd$TSSR <- (dd$start - dd$srise) / 24


offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
offdat$DSLS2 <- offdat$DSLS^2
offdat$LCC4 <- as.character(offdat$HAB_NALC2)
offdat$LCC4[offdat$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
offdat$LCC4[offdat$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
offdat$LCC4 <- factor(offdat$LCC4,
    c("DecidMixed", "Conif", "Open", "Wet"))
offdat$LCC2 <- as.character(offdat$LCC4)
offdat$LCC2[offdat$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
offdat$LCC2[offdat$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
offdat$LCC2 <- factor(offdat$LCC2, c("Forest", "OpenWet"))
table(offdat$LCC4, offdat$LCC2)
offdat$MAXDIS <- offdat$MAXDIS / 100

Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","DSLS","TSSR2","JDAY2","DSLS2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
    LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))
#offdat$OKp <- rowSums(is.na(offdat[,c("TSSR","JDAY","DSLS")])) == 0
#offdat$OKq <- rowSums(is.na(offdat[,c("TREE","LCC4")])) == 0
#Xp <- model.matrix(~TSSR+TSSR2+JDAY+JDAY2+DSLS+DSLS2, offdat[offdat$OKp,])
#Xq <- model.matrix(~LCC2+LCC4+TREE, offdat[offdat$OKq,])

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp

#spp <- "OVEN"
for (spp in sppp) {
    cat(spp, "\n");flush.console()
    p <- rep(NA, nrow(offdat))
    A <- q <- p

    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model
    mi <- bestmodelBAMspecies(spp, type="BIC")
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    #vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0

    p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])
    unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
    A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * offdat$MAXDIS[!OKq]^2)
    q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))

    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
    unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
    A[OKq] <- ifelse(unlim, pi * tau1^2, pi * offdat$MAXDIS[OKq]^2)
    q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))

    ii <- which(p == 0)
    p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])

    OFF[,spp] <- log(p) + log(A) + log(q)

}

(Ra <- apply(OFF, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,])) # BARS GCSP
which(!is.finite(Ra[2,]))

SPP <- sppp
save(OFF, SPP,
    file=file.path(ROOT, "out", "offsets-v3_2017-04-19.Rdata"))
offdat <- offdat[,c("PKEY","TSSR","JDAY","DSLS","TREE","LCC4","MAXDUR","MAXDIS")]
save(offdat,
    file=file.path(ROOT, "out", "offsets-v3data_2016-12-01.Rdata"))


