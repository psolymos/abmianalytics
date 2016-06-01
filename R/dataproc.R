##% Processing data for intactness, all-in-one, website things
##% P Solymos
##% May 31, 2016

library(RODBC)
library(mefa4)
library(raster)
library(sp)
library(rgdal)
library(maptools)

ROOT <- "e:/peter/bam/Apr2016"
ROOT2 <- "e:/peter/AB_data_v2016"

## proces ABMI data

e1 <- new.env()
e2 <- new.env()
load(file.path(ROOT2, "data", "species", "OUT_birdsrf_2016-05-27.Rdata"), envir=e1)
load(file.path(ROOT2, "data", "species", "OUT_birdssm_2016-05-30.Rdata"), envir=e2)

(mrf <- e1$m)
(msm <- Mefa(e2$xt, e2$x))
compare_sets(colnames(mrf), colnames(msm))
setdiff(colnames(mrf), colnames(msm))
setdiff(colnames(msm), colnames(mrf))
compare_sets(rownames(mrf), rownames(msm))
yy_abmi <- mbind(xtab(mrf), xtab(msm), fill=0)

pcabmi <- samp(mrf)
ssabmi <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

tmp <- do.call(rbind, sapply(levels(pcabmi$Label), strsplit, "_"))
colnames(tmp) <- c("Protocol", "OnOffGrid", "DataProvider", "SiteLabel", "YYYY", "Visit", "SubType", "BPC")
tmp2 <- sapply(tmp[,"SiteLabel"], strsplit, "-")
tmp3 <- sapply(tmp2, function(z) if (length(z)==1) "ABMI" else z[2])
tmp4 <- sapply(tmp2, function(z) if (length(z)==1) z[1] else z[3])
tmp <- data.frame(tmp, ClosestABMISite=tmp4)
tmp$DataProvider <- as.factor(tmp3)
tmp$Label <- with(tmp, paste(OnOffGrid, DataProvider, SiteLabel, YYYY, Visit, "PC", BPC, sep="_"))
tmp$Label2 <- with(tmp, paste(OnOffGrid, DataProvider, SiteLabel, YYYY, Visit, sep="_"))
tmp$ClosestABMISite <- as.integer(as.character(tmp$ClosestABMISite))
tmp$lat <- ssabmi$PUBLIC_LATTITUDE[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
tmp$long <- ssabmi$PUBLIC_LONGITUDE[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
#tmp$NatReg <- ssabmi$NATURAL_REGIONS[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
#tmp$boreal <- tmp$NatReg %in% c(c("Boreal", "Canadian Shield", "Foothills", "Rocky Mountain"))

pcabmi <- data.frame(pcabmi, tmp[match(pcabmi$Label, rownames(tmp)),])

## PKEY table and proper date format
PKEY_abmi <- nonDuplicated(pcabmi, pcabmi$Label, TRUE)
tmp <- PKEY_abmi$ADATE
tmp <- sapply(as.character(tmp), strsplit, split="-")
for (i in 1:length(tmp)) {
    if (length(tmp[[i]])<3) {
        tmp[[i]] <- rep("99", 3)
    }
}
table(sapply(tmp, "[[", 2))
for (i in 1:length(tmp)) {
    tmp[[i]][2] <- switch(tmp[[i]][2],
        "May"=5, "Jun"=6, "Jul"=7, "Aug"=8, "99"=99)
}
tmp <- sapply(tmp, function(z) paste("20",z[3],"-",z[2],"-",z[1], sep=""))
tmp[tmp=="2099-99-99"] <- NA
PKEY_abmi$Date <- as.POSIXct(tmp, tz="America/Edmonton")

## TSSR
Coor <- as.matrix(cbind(as.numeric(PKEY_abmi$long),as.numeric(PKEY_abmi$lat)))
JL <- as.POSIXct(PKEY_abmi$Date, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
PKEY_abmi$srise_MDT <- NA
PKEY_abmi$srise_MDT[subset] <- sr

tmp <- strsplit(as.character(PKEY_abmi$TBB_START_TIME), ":")
id <- sapply(tmp,length)==2
tmp <- tmp[id]
tmp <- as.integer(sapply(tmp,"[[",1)) + as.integer(sapply(tmp,"[[",2))/60
PKEY_abmi$start_time <- NA
PKEY_abmi$start_time[id] <- tmp
PKEY_abmi$srise <- PKEY_abmi$srise_MDT
PKEY_abmi$TSSR <- (PKEY_abmi$start_time - PKEY_abmi$srise) / 24 # MDT offset is 0

## Julian day
PKEY_abmi$jan1 <- as.Date(paste(PKEY_abmi$YEAR, "-01-01", sep=""))
PKEY_abmi$JULIAN <- as.numeric(as.Date(PKEY_abmi$Date)) - as.numeric(PKEY_abmi$jan1) + 1
PKEY_abmi$JULIAN[PKEY_abmi$JULIAN > 365] <- NA
PKEY_abmi$JDAY <- PKEY_abmi$JULIAN / 365

pcsm <- samp(msm)

tmp <- do.call(rbind, sapply(levels(pcsm$SITE_LABEL), strsplit, "_"))
colnames(tmp) <- c("Protocol", "OnOffGrid", "DataProvider", "SiteLabel", "YYYY", "Visit", "SubType", "BPC")
tmp2 <- sapply(tmp[,"SiteLabel"], strsplit, "-")
tmp3 <- sapply(tmp2, function(z) if (length(z)==1) "ABMI" else z[2])
tmp4 <- sapply(tmp2, function(z) if (length(z)==1) z[1] else z[3])
tmp <- data.frame(tmp, ClosestABMISite=tmp4)
tmp$DataProvider <- as.factor(tmp3)
tmp$Label <- with(tmp, paste(OnOffGrid, DataProvider, SiteLabel, YYYY, Visit, 
    "STATION", BPC, sep="_"))
tmp$Label2 <- with(tmp, paste(OnOffGrid, DataProvider, SiteLabel, YYYY, Visit, sep="_"))
tmp$ClosestABMISite <- as.integer(as.character(tmp$ClosestABMISite))
tmp$lat <- ssabmi$PUBLIC_LATTITUDE[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
tmp$long <- ssabmi$PUBLIC_LONGITUDE[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]

pcsm <- data.frame(pcsm, tmp[match(pcsm$SITE_LABEL, rownames(tmp)),])

## TSSR
Coor <- as.matrix(cbind(as.numeric(pcsm$long),as.numeric(pcsm$lat)))
JL <- as.POSIXct(pcsm$Start, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
pcsm$srise_MDT <- NA
pcsm$srise_MDT[subset] <- sr

tmp <- strsplit(as.character(pcsm$TBB_START_TIME), ":")
id <- sapply(tmp,length)==2
tmp <- tmp[id]
tmp <- as.integer(sapply(tmp,"[[",1)) + as.integer(sapply(tmp,"[[",2))/60
pcsm$start_time <- NA
pcsm$start_time[id] <- tmp
pcsm$srise <- pcsm$srise_MDT
pcsm$TSSR <- (pcsm$start_time - pcsm$srise) / 24 # MDT offset is 0

xx_rf <- with(PKEY_abmi, data.frame(
    PCODE="ABMI",
    PKEY=as.factor(paste0(Label, ":1")),
    SS=as.factor(Label),
    SITE=as.factor(Label2),
    YEAR=YEAR,
    TSSR=TSSR,
    JDAY=JDAY,
    MAXDUR=10,
    MAXDIS=Inf,
    TREE=NA
))
rownames(xx_rf) <- xx_rf$SS

xx_sm <- with(pcsm, data.frame(
    PCODE="ABMI",
    PKEY=as.factor(PKEY),
    SS=as.factor(SITE_LABEL),
    SITE=as.factor(Label2),
    YEAR=YEAR,
    TSSR=TSSR,
    JDAY=ToY/365,
    MAXDUR=3,
    MAXDIS=Inf,
    TREE=NA
))
rownames(xx_sm) <- xx_sm$PKEY

xx_abmi <- rbind(xx_rf, xx_sm)
#rownames(xx_abmi) <- xx_abmi$PKEY

compare_sets(rownames(yy_abmi), rownames(xx_abmi))
xx_abmi <- xx_abmi[rownames(yy_abmi),]

## use preprocessed national data set for BAM & BBS

load(file.path(ROOT, "out", "data_package_2016-04-18.Rdata"))
#load(file.path(ROOT, "out", "offsets-v3_2016-04-18.Rdata"))
TAX <- nonDuplicated(TAX, Species_ID, TRUE)
DAT <- data.frame(PKEY, SS[match(PKEY$SS, SS$SS),])
DAT$SS.1 <- NULL
DAT$PCODE.1 <- NULL
#rm(PCTBL, PKEY, SS)
table(duplicated(DAT$PKEY))
table(DAT$JURS, useNA="a")
DAT <- DAT[!is.na(DAT$JURS) & DAT$JURS == "AB",]
table(duplicated(DAT$PKEY))
rownames(DAT) <- DAT$PKEY

YY <- Xtab(ABUND ~ PKEY + SPECIES_ALL, PCTBL, cdrop="NONE")
ii <- sort(intersect(rownames(DAT), rownames(YY)))
DAT <- DAT[ii,]
YY <- YY[ii,]
#OFF <- OFF[ii,]
#YY <- YY[,colSums(YY) > 0]
## insert 0s into BAM+BBS:
#    BarnOwl
#    RedPhalarope
YY <- cBind(YY, BANO=0, REPH=0)
tax <- droplevels(TAX[colnames(YY),])
tax$Spp <- tax$English_Name
levels(tax$Spp) <- nameAlnum(levels(tax$Spp), capitalize="mixed", collapse="")
## fix labels:
#    MacgillivrayWarbler
#    BlackAndWhiteWarbler
levels(tax$Spp)[levels(tax$Spp)=="BlackandwhiteWarbler"] <- "BlackAndWhiteWarbler"
levels(tax$Spp)[levels(tax$Spp)=="MacGillivraysWarbler"] <- "MacgillivrayWarbler"
compare_sets(colnames(yy_abmi), levels(tax$Spp))
sort(setdiff(colnames(yy_abmi), levels(tax$Spp)))
#sort(setdiff(levels(tax$Spp), colnames(yy_abmi)))

## join the 2 tables by intersecting names
SPPx <- sort(intersect(colnames(yy_abmi), levels(tax$Spp)))
SPPx <- SPPx[SPPx != "UnidentifiedTern"]
TAX$Spp <- TAX$English_Name
levels(TAX$Spp) <- nameAlnum(levels(TAX$Spp), capitalize="mixed", collapse="")
levels(TAX$Spp)[levels(TAX$Spp)=="BlackandwhiteWarbler"] <- "BlackAndWhiteWarbler"
levels(TAX$Spp)[levels(TAX$Spp)=="MacGillivraysWarbler"] <- "MacgillivrayWarbler"
TAX <- TAX[TAX$Spp %in% SPPx,]
setdiff(SPPx, TAX$Spp)
TAX <- TAX[order(rownames(TAX)),]

#OFF <- OFF[,colnames(OFF) %in% rownames(TAX)]
YY <- YY[,rownames(TAX)]
yy_abmi <- yy_abmi[,as.character(TAX$Spp)]
colnames(yy_abmi) <- rownames(TAX)
YYY <- rBind(YY, yy_abmi)

## exclude 0 sum columns?
table(colSums(YYY))

xx_abmi$HAB_NALC1 <- NA
xx_abmi$HAB_NALC2 <- NA
compare_sets(colnames(DAT), colnames(xx_abmi))
setdiff(colnames(DAT), colnames(xx_abmi))
setdiff(colnames(xx_abmi), colnames(DAT))
cn <- intersect(colnames(DAT), colnames(xx_abmi))
DDAT <- rbind(DAT[,cn], xx_abmi[,cn])
DDAT <- DDAT[rownames(YYY),]

## habitat info

load(file.path(ROOT2, "out", "abmi_onoff", 
    "veg-hf-clim-reg_abmi-onoff_Birds-RF-SM_incl2015.Rdata"))
load(file.path(ROOT2, "out", "bambbs", 
    "veg-hf_bambbs_fix-fire_fix-age0.Rdata"))

compare_sets(rownames(YY), rownames(climPoint_bambbs))
compare_sets(rownames(yy_abmi), c(rownames(climRF), rownames(climSM)))
compare_sets(rownames(xx_rf), rownames(climRF))
compare_sets(rownames(xx_sm), rownames(climSM))

## RF: clean up issues
compare_sets(xx_rf$SS, climRF$Label2)
setdiff(xx_rf$SS, climRF$Label2)
setdiff(climRF$Label2, xx_rf$SS)
grep("-11_", rownames(climRF)) # this is hard stuff...

setdiff(rownames(climRF), rownames(xx_rf))

## SM: need to map SS to PKEY for rownames
## Centre is the right label, crude fix follows
climSM$Label <- gsub("_Center", "_1", climSM$Label)
climSM$Label0 <- gsub("_Center", "_1", climSM$Label0)
climSM$Label2 <- gsub("_Center", "_1", climSM$Label2)
for (i in 1:4) {
    rownames(dd150m_SM[[i]]) <- climSM$Label
    rownames(dd1km_SM[[i]]) <- climSM$Label
}


## we have DDAT and YYY
## now let us make CLIM
compare_sets(xx_sm$SS, climSM$Label)

climRF$YEAR <- climRF$Year
climRF$Year <- NULL
climSM$YEAR <- climSM$Year
climSM$Year <- NULL

climRF$PKEY <- climRF$Label

climSM <- climSM[match(xx_sm$SS, climSM$Label),]
climSM$PKEY <- rownames(xx_sm)
rownames(climSM) <- rownames(xx_sm)
climSM <- climSM[!is.na(climSM$Label),]

climPoint_bambbs$POINT_X <- DAT$X[match(climPoint_bambbs$PKEY, DAT$PKEY)]
climPoint_bambbs$POINT_Y <- DAT$Y[match(climPoint_bambbs$PKEY, DAT$PKEY)]

setdiff(colnames(climPoint_bambbs), colnames(climRF))
setdiff(colnames(climRF), colnames(climPoint_bambbs))

climPoint_bambbs$Part <- "BAMBBS"
climRF$Part <- "ABMIRF"
climSM$Part <- "ABMISM"
cn <- intersect(colnames(climPoint_bambbs), colnames(climRF))
CLIM <- rbind(climPoint_bambbs[,cn], climRF[,cn], climSM[,cn])
#CLIM <- CLIM[!is.na(CLIM$CTI),]
compare_sets(rownames(YYY), rownames(CLIM))
rrn <- intersect(rownames(CLIM), rownames(DDAT))

## now unify all habitat stuff into DD150m and DD1km
dd150m <- dd1km <- list()
for (i in 1:4) {
    tmp <- dd150m_SM[[i]]
    ii <- match(xx_sm$SS, rownames(tmp))
    ii <- ii[!is.na(ii)]
    tmp <- tmp[ii,]
    rownames(tmp) <- rownames(climSM)
    dd150m[[i]] <- rBind(dd150m_bambbs[[i]], dd150m_RF[[i]], tmp)
    dd150m[[i]] <- dd150m[[i]][rrn,]

    tmp <- dd1km_SM[[i]]
    ii <- match(xx_sm$SS, rownames(tmp))
    ii <- ii[!is.na(ii)]
    tmp <- tmp[ii,]
    rownames(tmp) <- rownames(climSM)
    dd1km[[i]] <- rBind(dd1km_bambbs[[i]], dd1km_RF[[i]], tmp)
    dd1km[[i]] <- dd1km[[i]][rrn,]
}

CLIM <- CLIM[rrn,]
YYY <- YYY[rrn,]
DDAT <- DDAT[rrn,]

data.frame(x=colSums(is.na(DDAT)))
data.frame(x=colSums(is.na(CLIM)))
sum(colSums(YYY)==0)

DAT <- cbind(droplevels(DDAT), droplevels(CLIM))
DAT$Part <- as.factor(DAT$Part)
DAT$MAP <- gsub(",", "", DAT$MAP)
DAT$MAP <- as.numeric(DAT$MAP)
YY <- YYY
rm(DDAT, CLIM, YYY)

dd150m[[5]] <- dd150m_RF[[5]]
dd1km[[5]] <- dd1km_RF[[5]]
names(dd150m) <- names(dd150m_bambbs)[c(1,2,3,4,6)]
names(dd1km) <- names(dd1km_bambbs)[c(1,2,3,4,6)]

## lookup tables & reclassed tables

ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

SoilKmCr <- groupSums(dd1km$soil_current, 2, ts[colnames(dd1km$soil_current), "UseInAnalysis"])
SoilKmRf <- groupSums(dd1km$soil_reference, 2, ts[colnames(dd1km$soil_reference), "UseInAnalysis"])
SoilPcCr <- groupSums(dd150m$soil_current, 2, ts[colnames(dd150m$soil_current), "UseInAnalysis"])
SoilPcRf <- groupSums(dd150m$soil_reference, 2, ts[colnames(dd150m$soil_reference), "UseInAnalysis"])

psoilhf <- groupSums(dd150m$soil_current, 2, ts[colnames(dd150m$soil_current), "UseInAnalysis"])
psoilhf <- psoilhf[,!(colnames(psoilhf) %in% c("SoilUnknown", "SoilWater",
    "HWater", "HFor", "SoftLin", "HardLin"))]
psoilhf <- as.matrix(psoilhf / rowSums(psoilhf))

#SoilKmCr <- SoilKmCr[,colnames(SoilKmCr) != "HFor"] # exclude forestry
SoilKmCr <- as.matrix(SoilKmCr / rowSums(SoilKmCr))
SoilKmRf <- as.matrix(SoilKmRf / rowSums(SoilKmRf))
#SoilPcCr <- SoilPcCr[,colnames(SoilPcCr) != "HFor"] # exclude forestry
SoilPcCr <- as.matrix(SoilPcCr / rowSums(SoilPcCr))
SoilPcRf <- as.matrix(SoilPcRf / rowSums(SoilPcRf))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")

if (FALSE) {
tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "Type"])
#tmp <- tmp[,!(colnames(tmp) %in% c("NonVeg", "Water", "HWater"))]
tmp <- as.matrix(tmp / rowSums(tmp))
OK <- tmp[,"XXX"] <= 0.8
table(OK)
tmp <- tmp[OK,colnames(tmp) != "XXX"]
any(is.na(tmp))
#tmp[is.na(tmp)] <- 0
iv <- find_max(tmp)
h <- iv$index
v <- iv$value
data.frame(table(h))
print(aggregate(data.frame(v=v), list(h=h), quantile, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)), 
    digits=3)
}

VegKmCr <- groupSums(dd1km$veg_current, 2, tv[colnames(dd1km$veg_current), "UseInAnalysis"])
VegKmRf <- groupSums(dd1km$veg_reference, 2, tv[colnames(dd1km$veg_reference), "UseInAnalysis"])
VegPcCr <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "UseInAnalysis"])
VegPcRf <- groupSums(dd150m$veg_reference, 2, tv[colnames(dd150m$veg_reference), "UseInAnalysis"])

VegKmCr <- as.matrix(VegKmCr / rowSums(VegKmCr))
VegKmRf <- as.matrix(VegKmRf / rowSums(VegKmRf))
VegPcCr <- as.matrix(VegPcCr / rowSums(VegPcCr))
VegPcRf <- as.matrix(VegPcRf / rowSums(VegPcRf))


## hab & CC in north

## 150m scale reference dominant habitat
tmp <- groupSums(dd150m$veg_reference, 2, tv[colnames(dd150m$veg_reference), "Type"])
tmp <- tmp[,!(colnames(tmp) %in% c("XXX"))]
tmp <- as.matrix(tmp / rowSums(tmp))
iv <- find_max(tmp)
DAT$hab0 <- iv$index
DAT$hab0_value <- iv$value

tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "Type"])
pveghf <- tmp[,!(colnames(tmp) %in% c("XXX"))]
pveghf <- as.matrix(pveghf / rowSums(pveghf))

## 150m scale current dominant habitat with CC
tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "TypeCC"])
tmp <- tmp[,!(colnames(tmp) %in% c("XXX","SoftLin","HardLin"))]
tmp <- as.matrix(tmp / rowSums(tmp))
iv <- find_max(tmp)
DAT$hab1cc <- iv$index
DAT$hab1cc_value <- iv$value

DAT$hab1 <- as.character(DAT$hab1cc)
ii <- !is.na(DAT$hab1cc) & !is.na(DAT$hab0) &
    DAT$hab1cc == "CC" & DAT$hab0 %in% c("Conif","Decid","Mixwood","Pine")
DAT$hab1[ii] <- as.character(DAT$hab0[ii])
DAT$isCC <- ifelse(ii, 1L, 0L)
## reset non-merchendizable CC classes to backfilled (isCC is 0)
ii <- !is.na(DAT$hab1cc) & !is.na(DAT$hab0) & DAT$hab1 == "CC"
DAT$hab1[ii] <- as.character(DAT$hab0[ii])

DAT$hab0 <- relevel(DAT$hab0, "Decid")
DAT$hab1 <- relevel(as.factor(DAT$hab1), "Decid")
DAT$hab1cc <- relevel(DAT$hab1cc, "Decid")

table(DAT$hab1, DAT$hab1cc, useNA="a")
table(DAT$hab1,DAT$isCC)


## age 

AgePtCr <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "AGE"])
AgePtCr <- as.matrix(AgePtCr / rowSums(AgePtCr))
## exclude unknown (0) and non-forest (blank)
AgePtCr <- t(AgePtCr[,c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9")])
AgeMin <- structure(c(0,10,20,40,60,80,100,120,140,160)/200,
    names=c("R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
DAT$wtAge <- colSums(AgePtCr * AgeMin) / colSums(AgePtCr)
DAT$wtAge[is.na(DAT$wtAge)] <- 0
FORclasses <- c("Conif","Decid","Mixwood","Pine","BSpr", "Larch")
DAT$wtAge[!(DAT$hab1 %in% FORclasses)] <- 0
## recent fire (non CC)
DAT$isRR <- ifelse(DAT$wtAge < 0.1 & !DAT$isCC, 1L, 0L)
DAT$isFor <- DAT$hab1 %in% FORclasses
DAT$isConif <- DAT$hab1 %in% c("Conif","Pine","BSpr", "Larch")
DAT$isRR[!DAT$isFor] <- 0
DAT$hab1b <- DAT$hab1
levels(DAT$hab1b) <- c(levels(DAT$hab1b), "Burn")
DAT$hab1b[DAT$isRR == 1] <- "Burn"

table(burn=DAT$isRR, cc=DAT$isCC)
table(DAT$hab1b,DAT$hab1)

## forestry convergence
## fCC1: linear
MAXFOR <- 50/200
DAT$fCC1 <- 0
DAT$fCC1[DAT$isCC==1] <- pmax(0, 1 - (DAT$isCC * DAT$wtAge/MAXFOR)[DAT$isCC==1])
plot(fCC1 ~ wtAge, DAT[DAT$isCC==1,])
## Dave`s recovery trajectories
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3, 
    82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85, 
    88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100
DAT$fCC2 <- 0
tmp1 <- approxfun(age, decid)(DAT$wtAge)
tmp1[is.na(tmp1)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- DAT$isFor & !DAT$isConif & DAT$isCC==1
DAT$fCC2[ii] <- tmp1[ii]
tmp2 <- approxfun(age, conif)(DAT$wtAge)
tmp2[is.na(tmp2)] <- 0 # this happens when wtAge > 0.4 (out of range of approx)
ii <- DAT$isConif & DAT$isCC==1
DAT$fCC2[ii] <- tmp2[ii]
plot(DAT$wtAge, DAT$fCC1, col=2, pch=".")
points(DAT$wtAge, DAT$fCC2, col=4, pch=".")
sum(is.na(DAT$fCC2))

if (FALSE) {
## Dave`s recovery trajectories
age <- c(0, 1:20*4)/200
conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3, 
    82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85, 
    88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100
Age <- seq(0, 100, by=0.5) / 200
plot(Age*200, approxfun(age, decid)(Age), type="l", col=2, lwd=2, 
    xlab="Years since disturbance", ylab="Transformed value")
lines(Age*200, approxfun(age, conif)(Age), col=4, lwd=2)
lines(Age*200, pmax(0, 1 - (Age * 200 / 50)), col=3, lwd=2)
lines(Age*200, 1/((Age*200+1)^(1/3)), col=1, lwd=2)
legend("topright", col=1:4, lty=1, lwd=2, 
    legend=c("cubic-root","expert decid","linear","expert conif"))
}

table(DAT$hab1,round(DAT$wtAge,1),useNA="a")
table(DAT$hab1,DAT$isCC,useNA="a")
table(DAT$hab1,DAT$isRR,useNA="a")
table(DAT$isRR,DAT$isCC,useNA="a")

DAT$ECage <- factor("", levels=c("", "R","A", "B", "C", "D"))
DAT$ECage[DAT$isFor & DAT$wtAge >= 0] <- "R"
DAT$ECage[DAT$wtAge >= 10/200] <- "A"
DAT$ECage[DAT$wtAge >= 20/200] <- "B"
DAT$ECage[DAT$wtAge >= 40/200] <- "C"
DAT$ECage[DAT$isFor & !DAT$isConif & DAT$wtAge >= 60/200] <- "D" # decid
DAT$ECage[DAT$isConif & DAT$wtAge >= 80/200] <- "D" # conif
DAT$hab1ec <- interaction(DAT$hab1, DAT$ECage, drop=TRUE, sep="")
DAT$hab1ec <- relevel(DAT$hab1ec, "DecidD")

## LCC_combo & TREE: use hab info

tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "EC_AGE"])
tmp <- as.matrix(tmp / rowSums(tmp))
DAT$ClosedCanopy <- tmp[,"C"] + tmp[,"D"]

DAT$TREE[is.na(DAT$TREE)] <- DAT$ClosedCanopy[is.na(DAT$TREE)]
DAT$TREE[is.na(DAT$TREE)] <- DAT$ClosedCanopy[is.na(DAT$TREE)]
DAT$TREE3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
DAT$TREE3[DAT$TREE < 0.25] <- "Open"
DAT$TREE3[DAT$TREE >= 0.25 & DAT$TREE < 0.60] <- "Sparse"
DAT$TREE3[DAT$TREE >= 0.60] <- "Dense"

## ROAD x habitat interaction based on forest cover
DAT$habCl <- ifelse(DAT$ClosedCanopy > 0.5, 1L, 0L) # EC age is C or D

## land cover for offsets
if (FALSE) {
tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "LCC5"])
tmp <- as.matrix(tmp / rowSums(tmp))
DAT$hab_lcc <- find_max(tmp)$index
DAT$hab_lcc <- as.integer(as.character(DAT$hab_lcc))
DAT$hab_lcc <- factor(DAT$hab_lcc, 1:5)
}

#DAT$LCC_combo[is.na(DAT$LCC_combo)] <- DAT$hab_lcc[is.na(DAT$LCC_combo)] 

## soil in south

round(100*colMeans(SoilPcRf[SoilPcRf[,"SoilUnknown"]==0,]),2)
round(100*colMeans(SoilPcCr[SoilPcRf[,"SoilUnknown"]==0,]),2)

tmp <- SoilPcRf[,!(colnames(SoilPcRf) %in% c("SoilWater","SoilWetland"))]
iv <- find_max(tmp / rowSums(tmp))
DAT$soil0 <- iv$index
DAT$soil0_value <- iv$value
DAT$soil0[DAT$soil0 == "SoilUnknown"] <- NA
DAT$soil0 <- droplevels(DAT$soil0)
DAT$soil0 <- relevel(DAT$soil0, "Productive")
table(DAT$soil0, useNA="a")

tmp <- SoilPcCr[,!(colnames(SoilPcCr) %in% c("SoilWater","SoilWetland",
    "HWater","SoftLin", "HardLin", "HFor"))]
iv <- find_max(tmp / rowSums(tmp))
DAT$soil1 <- iv$index
DAT$soil1_value <- iv$value
DAT$soil1[DAT$soil1 == "SoilUnknown"] <- NA
DAT$soil1 <- droplevels(DAT$soil1)
DAT$soil1 <- relevel(DAT$soil1, "Productive")
table(DAT$soil1, useNA="a")
ii <- !is.na(DAT$soil1) & !is.na(DAT$soil0) &
    !(DAT$soil1 %in% c("Cult","UrbInd"))
DAT$soil1[ii] <- as.character(DAT$soil0[ii])

table(DAT$soil1, DAT$soil0, useNA="a")

DAT$soil1v <- DAT$soil1
levels(DAT$soil1v)[levels(DAT$soil1) %in% c("Saline", "Clay")] <- "SalineAndClay"
table(DAT$soil1, DAT$soil1v, useNA="a")

DAT$soil1vv <- DAT$soil1v
levels(DAT$soil1vv)[levels(DAT$soil1v) %in% c("RapidDrain","SalineAndClay")] <- "NonProductive"
table(DAT$soil1, DAT$soil1vv, useNA="a")

## Water

## need to exclude bare/OW sites based on v=0
## use WET and WETWATER columns to calculate % wet and wetwater in buffers
DAT$pWater <- rowSums(dd150m$veg_current[,c("Water",
    "BorrowpitsDugoutsSumps","MunicipalWaterSewage","Reservoirs","Canals")]) /
    rowSums(dd150m$veg_current)

## ROAD and sof lin PC

DAT$Road_PC <- rowSums(dd150m$veg_current[,c("RailHardSurface",
    "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / rowSums(dd150m$veg_current)
DAT$ROAD01 <- ifelse(DAT$Road_PC >= 0.1, 1L, 0L)
DAT$ROAD01[DAT$PCODE == "BBSAB"] <- 1L
DAT$SoftLin_PC <- rowSums(dd150m$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge")]) /
    rowSums(dd150m$veg_current)

## ARU

DAT$ARU2 <- factor("TRAD", c("TRAD","ARU"))
DAT$ARU2[DAT$PCODE %in% c("ABMI","EMCLA","EMCLA2014")] <- "ARU"
DAT$ARU2[DAT$PCODE %in% c()] <- "RF"
table(DAT$PCODE, DAT$ARU2)

DAT$ARU3 <- factor("TRAD", c("TRAD","SM","RF"))
DAT$ARU3[DAT$PCODE %in% c("EMCLA","EMCLA2014")] <- "SM"
DAT$ARU3[DAT$PCODE %in% c("ABMI")] <- "RF"
DAT$ARU3[DAT$Part %in% c("ABMISM")] <- "SM"
table(DAT$PCODE, DAT$ARU3)

## transformations

## mapping projection
XYlatlon <- DAT[,c("POINT_X", "POINT_Y")]
coordinates(XYlatlon) <- ~ POINT_X + POINT_Y
proj4string(XYlatlon) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
XY <- as.data.frame(spTransform(XYlatlon, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))
DAT$X <- XY[,"POINT_X"]
DAT$Y <- XY[,"POINT_Y"]

#tmp <- DAT[,c("SLP"  ,        "ASP"     ,     "TRI"     ,     "CTI" )]

transform_CLIM <- function(x, ID="PKEY") {
    z <- x[,ID,drop=FALSE]
    z$xlong <- (x$POINT_X - (-113.7)) / 2.15
    z$xlat <- (x$POINT_Y - 53.8) / 2.28
    z$xAHM <- (x$AHM - 0) / 50
    z$xPET <- (x$PET - 0) / 800
    z$xFFP <- (x$FFP - 0) / 130
    z$xMAP <- (x$MAP - 0) / 2200
    z$xMAT <- (x$MAT - 0) / 6
    z$xMCMT <- (x$MCMT - 0) / 25
    z$xMWMT <- (x$MWMT - 0) / 20
    z$xASP <- x$ASP
    z$xSLP <- log(x$SLP + 1)
    z$xTRI <- log(x$TRI / 5)
    z$xCTI <- log((x$CTI + 1) / 10)
    z
}
#DAT$MAP <- gsub(",", "", DAT$MAP)
#DAT$MAP <- as.numeric(DAT$MAP)
tc <- transform_CLIM(DAT)
DAT <- data.frame(DAT, tc[,-1])
#DAT$xSLP[is.na(DAT$xSLP)] <- mean(DAT$xSLP, na.rm=TRUE)
#DAT$xASP[is.na(DAT$xASP)] <- mean(DAT$xASP, na.rm=TRUE)
#DAT$xTRI[is.na(DAT$xTRI)] <- mean(DAT$xTRI, na.rm=TRUE)
#DAT$xCTI[is.na(DAT$xCTI)] <- mean(DAT$xCTI, na.rm=TRUE)
#DAT$PKEY.1 <- NULL


## surrounding HF at 1km scale

#DAT$THF_KM <- rowSums(VegKmCr[,setdiff(colnames(VegKmCr), colnames(VegKmRf))])
DAT$THF_KM <- rowSums(dd1km$veg_current[,setdiff(colnames(dd1km$veg_current), colnames(dd1km$veg_reference))]) / rowSums(dd1km$veg_current)
DAT$Lin_KM <- rowSums(dd1km$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailHardSurface", "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / rowSums(dd1km$veg_current)
DAT$Nonlin_KM <- DAT$THF_KM - DAT$Lin_KM
DAT$Cult_KM <- rowSums(dd1km$veg_current[,c("CultivationCropPastureBareground",
    "HighDensityLivestockOperation")]) / rowSums(dd1km$veg_current)
DAT$Noncult_KM <- DAT$THF_KM - DAT$Cult_KM
CClabs <- colnames(dd1km$veg_current)[grep("CC", colnames(dd1km$veg_current))]
DAT$Succ_KM <- rowSums(dd1km$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge", 
    CClabs)]) / rowSums(dd1km$veg_current)
DAT$Alien_KM <- DAT$THF_KM - DAT$Succ_KM

DAT$THF2_KM <- DAT$THF_KM^2
DAT$Succ2_KM <- DAT$Succ_KM^2
DAT$Alien2_KM <- DAT$Alien_KM^2
DAT$Noncult2_KM <- DAT$Noncult_KM^2
DAT$Nonlin2_KM <- DAT$Nonlin_KM^2

if (FALSE) {

## PC level HF summaries following that of 1km above
#DAT$THF_PC <- rowSums(VegPcCr[,setdiff(colnames(VegPcCr), colnames(VegPcRf))])
DAT$THF_PC <- rowSums(dd150m$veg_current[,setdiff(colnames(dd150m$veg_current), colnames(dd150m$veg_reference))]) / rowSums(dd150m$veg_current)
DAT$Lin_PC <- rowSums(dd150m$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailHardSurface", "RailVegetatedVerge","RoadHardSurface","RoadTrailVegetated",
    "RoadVegetatedVerge")]) / rowSums(dd150m$veg_current)
DAT$Nonlin_PC <- DAT$THF_PC - DAT$Lin_PC
DAT$Cult_PC <- rowSums(dd150m$veg_current[,c("CultivationCropPastureBareground",
    "HighDensityLivestockOperation")]) / rowSums(dd150m$veg_current)
DAT$Noncult_PC <- DAT$THF_PC - DAT$Cult_PC
CClabs <- colnames(dd150m$veg_current)[grep("CC", colnames(dd1km$veg_current))]
DAT$Succ_PC <- rowSums(dd150m$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge", 
    CClabs)]) / rowSums(dd150m$veg_current)
DAT$Alien_PC <- DAT$THF_PC - DAT$Succ_PC

DAT$THF2_PC <- DAT$THF_PC^2
DAT$Succ2_PC <- DAT$Succ_PC^2
DAT$Alien2_PC <- DAT$Alien_PC^2
DAT$Noncult2_PC <- DAT$Noncult_PC^2
DAT$Nonlin2_PC <- DAT$Nonlin_PC^2


}

## 1km scale wet/water
DAT$WetKM <- rowSums(dd1km$veg_current[,tv[colnames(dd1km$veg_current), "WET"]==1]) / rowSums(dd1km$veg_current)
DAT$WaterKM <- rowSums(dd1km$veg_current[,tv[colnames(dd1km$veg_current), "WATER"]==1]) / rowSums(dd1km$veg_current)
DAT$WetWaterKM <- rowSums(dd1km$veg_current[,tv[colnames(dd1km$veg_current), "WETWATER"]==1]) / rowSums(dd1km$veg_current)

## HSH matrix using hab1ec (EC classes)
cn <- colnames(dd1km$veg_current)
HSH <- groupSums(dd1km$veg_current, 2, paste0(tv[cn,"Type"], tv[cn,"EC_AGE"]))
HSH <- as.matrix(HSH / rowSums(HSH))
HSH <- HSH[,levels(DAT$hab1ec)]
DAT$HSH <- 0
DAT$HSH2 <- 0

## restricted data use

PCODE_useOK <- c("ABCAWAWEST", "MGLE", "FTL", "THIN", 
    "LMWELL", "EMB-ASP", "EMB-BS", "EMB-NOISE", 
    "ECJOSM", "ECJOSM_JRB", 
    "EMCLA", "EMCLA2014", "ABMI", "BBSAB")
DAT$useOK <- DAT$PCODE %in% PCODE_useOK
DAT$useOK[DAT$YEAR > 2007 & DAT$PCODE == "CL"] <- TRUE # Calling Lake


## YR

DAT$YR <- (DAT$YEAR - 1997) / 10
DAT$YR5 <- cut(DAT$YEAR, c(1996, 2001, 2006, 2010, 2015))
table(DAT$YEAR, DAT$YR5, useNA="a")
table(DAT$YR5, useNA="a")

## Bootstrap blocks

## NR and lat/long
DAT$bootid <- as.character(DAT$NRNAME)
tmp <- as.character(interaction(ifelse(DAT$POINT_Y < (56.5), "S", "N"),
    ifelse(DAT$POINT_X < (-115.5), "W", "E"), sep=""))
DAT$bootid[DAT$bootid %in% c("Boreal","Canadian Shield")] <-
    tmp[DAT$bootid %in% c("Boreal","Canadian Shield")]
#with(DAT, plot(X, Y, col=as.factor(bootid)))
table(DAT$bootid, DAT$YR5)
DAT$bootid <- interaction(DAT$bootid, as.integer(DAT$YR5))

## use in south
DAT$useSouth <- FALSE
DAT$useSouth[DAT$NRNAME %in% c("Grassland", "Parkland")] <- TRUE
DAT$useSouth[DAT$NSRNAME %in% c("Dry Mixedwood")] <- TRUE
DAT$useSouth[DAT$useSouth & DAT$POINT_Y > 56.7] <- FALSE # ~80 points
DAT$useSouth[DAT$useSouth & DAT$pWater > 0.5] <- FALSE
DAT$useSouth[DAT$useSouth & SoilPcRf[,"SoilUnknown"] > 0] <- FALSE
DAT$useSouth[DAT$useSouth & is.na(DAT$soil1)] <- FALSE

## use in north
DAT$useNorth <- !is.na(DAT$hab1) # veg info available
DAT$useNorth[DAT$NRNAME == "Grassland"] <- FALSE
DAT$useNorth[DAT$useNorth & DAT$pWater > 0.5] <- FALSE
## this is not necessary, plenty of data
#DAT$useNorth[DAT$useNorth & DAT$POINT_Y < 52.8] <- FALSE

## use JOSM
DAT$useJosm <- !is.na(DAT$hab1) # veg info available
DAT$useJosm[DAT$useJosm & DAT$pWater > 0.5] <- FALSE

## within year visits

DAT <- DAT[sample.int(nrow(DAT), nrow(DAT)),]
## ABMI sites are used as SS, but only points are revisited 
## (i.e. spatial replication should not count)
#DAT$SS0 <- as.character(DAT$SS)
#tmp <- DAT$SS0[DAT$PCODE=="ABMI"]

DAT$SS_YR <- interaction(DAT$SS, DAT$YEAR, drop=TRUE)
table(table(DAT$SS_YR))
dup <- unique(DAT$SS_YR[duplicated(DAT$SS_YR)])
tmp1 <- DAT[DAT$SS_YR %in% dup, c("PKEY","SS_YR")]
rownames(tmp1) <- tmp1$PKEY
dim(tmp1)
tmp2 <- nonDuplicated(tmp1[sample.int(nrow(tmp1)),], SS_YR)
rownames(tmp2) <- tmp2$PKEY
dim(tmp2)
dim(tmp1) - dim(tmp2)
tmp1 <- tmp1[setdiff(rownames(tmp1), rownames(tmp2)),]
dim(tmp1)
dim(tmp1) + dim(tmp2)
DAT$Revisit <- DAT$PKEY %in% rownames(tmp1)
table(DAT$Revisit)
table(DAT$PCODE,DAT$Revisit)
DAT$SS_YR <- NULL

DAT$SITE <- as.character(DAT$SITE)
DAT$SITE[is.na(DAT$SITE)] <- as.character(DAT$SS)[is.na(DAT$SITE)]
DAT$SITE <- as.factor(DAT$SITE)

## derived variables

DAT$xlat2 <- DAT$xlat^2
DAT$xlong2 <- DAT$xlong^2
DAT$wtAge2 <- DAT$wtAge^2
DAT$wtAge05 <- DAT$wtAge^0.5

#DAT$hab_lcc3 <- DAT$hab_lcc
#levels(DAT$hab_lcc3) <- c("1", "1", "2", "2", "3")
#DAT$hab_lcc2 <- DAT$hab_lcc
#levels(DAT$hab_lcc2) <- c("1", "1", "1", "1", "2")

## simply treat Mixed as intercept with Decid as reference
DAT$isMix <- ifelse(DAT$hab1 == "Mixwood", 1L, 0L)
DAT$isWSpruce <- ifelse(DAT$hab1 == "Conif", 1L, 0L)
DAT$isPine <- ifelse(DAT$hab1 == "Pine", 1L, 0L)
DAT$isBSpruce <- ifelse(DAT$hab1 == "BSpr", 1L, 0L)
DAT$isLarch <- ifelse(DAT$hab1 == "Larch", 1L, 0L)
DAT$isBSLarch <- ifelse(DAT$hab1 %in% c("BSpr", "Larch"), 1L, 0L)
DAT$isUpCon <- ifelse(DAT$hab1 %in% c("Conif", "Pine"), 1L, 0L)
DAT$isCon <- ifelse(DAT$hab1 %in% c("BSpr", "Larch",
    "Conif", "Pine"), 1L, 0L)

source("~/repos/abmianalytics/R/analysis_models.R")
source("~/repos/bragging/R/glm_skeleton.R")
compare_sets(getTerms(modsSoil, "list"), colnames(DAT))
setdiff(getTerms(modsSoil, "list"), colnames(DAT))
stopifnot(length(setdiff(getTerms(modsSoil, "list"), colnames(DAT)))==0)
compare_sets(getTerms(modsVeg, "list"), colnames(DAT))
setdiff(getTerms(modsVeg, "list"), colnames(DAT))
stopifnot(length(setdiff(getTerms(modsVeg, "list"), colnames(DAT)))==0)

## screen for CTI if that is part of the model!
keep <- DAT$YEAR >= 1997 & !is.na(DAT$hab1) & DAT$pWater <= 0.5
#keep <- DAT$YEAR >= 1997 & !is.na(DAT$hab1) & !DAT$Revisit & DAT$pWater <= 0.5
#keep <- !is.na(DAT$hab1) & !DAT$Revisit & DAT$pWater <= 0.5

DAT$keep <- keep
plot(DAT$X, DAT$Y, col=ifelse(DAT$keep, 1, 2), pch=19, cex=0.2)
data.frame(x=colSums(is.na(DAT)))
DAT$PKEY.1 <- NULL
DAT$YEAR.1 <- NULL

## offsets
## QPADv3

library(mefa4)
library(QPAD)
#source("~/repos/bamanalytics/R/dataprocessing_functions.R")
#load(file.path(ROOT, "out", paste0("data_package_2016-04-18.Rdata")))

load_BAM_QPAD(3)
getBAMversion()
sppp <- getBAMspecieslist()
compare_sets(sppp, colnames(YY))
sppp <- intersect(sppp, colnames(YY))

offdat <- DAT[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS",
    "TREE","TREE3","HAB_NALC1","HAB_NALC2")]
summary(offdat)

offdat$TSSR[is.na(offdat$TSSR)] <- mean(offdat$TSSR, na.rm=TRUE)
offdat$JDAY[is.na(offdat$JDAY)] <- mean(offdat$JDAY, na.rm=TRUE)
offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
#offdat$DSLS2 <- offdat$DSLS^2

#table(DAT$hab1, DAT$HAB_NALC2)
tmp <- groupSums(dd150m$veg_current, 2, tv[colnames(dd150m$veg_current), "LCC4"])
tmp <- as.matrix(tmp / rowSums(tmp))
tmp2 <- find_max(tmp[is.na(offdat$HAB_NALC2),])

offdat$LCC4 <- as.character(offdat$HAB_NALC2)
offdat$LCC4[offdat$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
offdat$LCC4[offdat$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
offdat$LCC4 <- factor(offdat$LCC4,
        c("DecidMixed", "Conif", "Open", "Wet"))
offdat$LCC4[is.na(offdat$LCC4)] <- tmp2$index
offdat$LCC2 <- as.character(offdat$LCC4)
offdat$LCC2[offdat$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
offdat$LCC2[offdat$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
offdat$LCC2 <- factor(offdat$LCC2, c("Forest", "OpenWet"))
table(offdat$LCC4, offdat$LCC2, useNA="a")
offdat$MAXDIS <- offdat$MAXDIS / 100

Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
    LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))
OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- rownames(offdat)
colnames(OFF) <- sppp

#spp <- "OVEN"
for (spp in sppp) {
p <- rep(NA, nrow(offdat))
A <- q <- p

## constant for NA cases
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
## best model
mi <- bestmodelBAMspecies(spp, type="BIC", model.sra=0:8)
cat(spp, unlist(mi), "\n");flush.console()
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
A[OKq] <- ifelse(unlim, pi * tau1, pi * offdat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))

ii <- which(p == 0)
p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])

OFF[,spp] <- log(p) + log(A) + log(q)

}

(Ra <- apply(OFF, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))

OFFmean <- log(rowMeans(exp(OFF)))

#YY <- YY[rownames(DAT),]
#HSH <- HSH[rownames(DAT),]
#pveghf <- pveghf[rownames(DAT),]
#psoilhf <- psoilhf[rownames(DAT),]
save(DAT, YY, OFF, OFFmean, TAX, HSH, # pveghf, psoilhf, 
    file=file.path(ROOT2, "out", "birds", "data", "data-full-withrevisit.Rdata"))

## subsets -------------------------------------------------------------

library(mefa4)
#ROOT <- "e:/peter/bam/Apr2016"
ROOT2 <- "e:/peter/AB_data_v2016"

load(file.path(ROOT2, "out", "birds", "data", "data-full-withrevisit.Rdata"))
source("~/repos/abmianalytics/R/analysis_models.R")


DAT <- droplevels(DAT[DAT$keep,])
YY <- YY[rownames(DAT),]
OFF <- OFF[rownames(DAT),]
OFFmean <- OFFmean[rownames(DAT)]

plot(DAT$X, DAT$Y, col=ifelse(DAT$useOK, 1, 2), pch=19, cex=0.2)

DATs <- droplevels(DAT[DAT$useSouth & DAT$useOK,])
DATn <- droplevels(DAT[DAT$useNorth & DAT$useOK & DAT$POINT_Y>50,])
DATj <- droplevels(DAT)

aas <- colSums(is.na(DATs))
aan <- colSums(is.na(DATn))
aaj <- colSums(is.na(DATj))
aas[aas>0]
aan[aan>0]
aaj[aaj>0]

sapply(list(full=DAT, southOK=DATs, northOK=DATn, JOSM=DATj), nrow)

## bootids

library(detect)
B <- 239

bbfun <- function(DAT1, B, out=0.1, seed=1234) {
    set.seed(seed)
    ## make sure that time intervals are considered as blocks
    ## keep out 10% of the data for validation
    id2 <- list()
    for (l in levels(DAT1$bootid)) {
        sset0 <- which(DAT1$bootid == l)
        ## resample revisits
        if (length(sset0) > 1)
            sset0 <- sample(sset0)
        dpl <- DAT1$Revisit[sset0]
        ## combine dpl=F with !duplicated dpl=T elements
        tmp <- sset0[dpl]
        sset <- c(sset0[!dpl], tmp[!duplicated(DAT1[tmp, "SS"])])
        id2[[l]] <- sample(sset, floor(length(sset) * 1-out), FALSE)
    }
    KEEP_ID <- unname(unlist(id2))
    HOLDOUT_ID <- setdiff(seq_len(nrow(DAT1)), KEEP_ID)

    DAT1k <- DAT1[KEEP_ID,]
    DAT1 <- DAT1[c(KEEP_ID, HOLDOUT_ID),]
    DAT1k$SITE <- droplevels(DAT1k$SITE)
    BB1 <- hbootindex(DAT1k$SITE, DAT1k$bootid, B=B)
    BB1
}
BBs <- bbfun(DATs, B)
BBn <- bbfun(DATn, B)
BBj <- bbfun(DATj, B)

## figure out sets of species to analyze

## pa maps and habitat suitability: DAT[useOK,]
## S/N: nmin=25
## look at taxonomy???

DAT0 <- DAT
YY0 <- YY
OFF0 <- OFF
OFFmean0 <- OFFmean
#TAX0 <- TAX
#HSH0 <- HSH
nmin <- 25

## south
DAT <- DATs
YY <- YY0[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= nmin]
BB <- BBs
OFF <- OFF0[rownames(DAT),colnames(OFF0) %in% colnames(YY)]
OFFmean <- OFFmean0[rownames(DAT)]
mods <- modsSoil
mods$Topo <- NULL
names(mods)
dim(YY)
save(DAT, YY, OFF, OFFmean, mods, BB, # no HSH in the south
    file=file.path(ROOT2, "out", "birds", "data", "data-south.Rdata"))

DAT <- DATn
YY <- YY0[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= nmin]
BB <- BBn
#HSH <- HSH0[rownames(DAT),]
OFF <- OFF0[rownames(DAT),colnames(OFF0) %in% colnames(YY)]
OFFmean <- OFFmean0[rownames(DAT)]
mods <- modsVeg
mods$Topo <- NULL
names(mods)
dim(YY)
save(DAT, YY, OFF, OFFmean, mods, BB, # HSH,
    file=file.path(ROOT2, "out", "birds", "data", "data-north.Rdata"))

DAT <- DATj
YY <- YY0[rownames(DAT),]
YY <- YY[,colSums(YY>0) >= nmin]
YY <- YY[,colnames(YY) %in% colnames(OFF0)]
BB <- BBj
#HSH <- HSH0[rownames(DAT),]
OFF <- OFF0[rownames(DAT),colnames(OFF0) %in% colnames(YY)]
OFFmean <- OFFmean0[rownames(DAT)]
mods <- modsVeg
mods$Topo <- NULL
names(mods)
dim(YY)
save(DAT, YY, OFF, mods, BB, # HSH, OFFmean, 
    file=file.path(ROOT2, "out", "birds", "data", "data-josm.Rdata"))

## TODO
## make data with all detections per location for useOK surveys
## take the max and binarize for each SS (goupSums and nonDuplicated will do it)
## use this for detection maps
## use this for wrsi calculations
## create all the relevant buffer percentage summaries (i.e. not only hab1!)

## opticut stuff
if (FALSE) {

dat <- droplevels(DAT[DAT$PCODE == "ABMI",])
yy <- YY[rownames(dat),]
yy <- yy[,colSums(yy) > 0]
tax <- TAX[colnames(yy),]

ClosedCalopyForest <- c("Conif1", "Conif2", "Conif3", "Conif4", 
    "Conif5", "Conif6", "Conif7", "Conif8", "Conif9", 
    "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6", "Decid7", 
    "Decid8", "Decid9", "Mixwood1", "Mixwood2", 
    "Mixwood3", "Mixwood4", "Mixwood5", "Mixwood6", "Mixwood7", "Mixwood8", 
    "Mixwood9", "Pine1", "Pine2", "Pine3", "Pine4", 
    "Pine5", "Pine6", "Pine7", "Pine8", "Pine9", "Swamp-Conif1", 
    "Swamp-Conif2", "Swamp-Conif3", "Swamp-Conif4", "Swamp-Conif5", 
    "Swamp-Conif6", "Swamp-Conif7", "Swamp-Conif8", "Swamp-Conif9", 
    "Swamp-Decid1", "Swamp-Decid2", 
    "Swamp-Decid3", "Swamp-Decid4", "Swamp-Decid5", "Swamp-Decid6", 
    "Swamp-Decid7", "Swamp-Decid8", "Swamp-Decid9", 
    "Swamp-Mixwood1", "Swamp-Mixwood2", "Swamp-Mixwood3", 
    "Swamp-Mixwood4", "Swamp-Mixwood5", "Swamp-Mixwood6", "Swamp-Mixwood7", 
    "Swamp-Mixwood8", "Swamp-Mixwood9", 
    "Swamp-Pine1", "Swamp-Pine2", "Swamp-Pine3", "Swamp-Pine4", "Swamp-Pine5", 
    "Swamp-Pine6", "Swamp-Pine7", "Swamp-Pine8", "Swamp-Pine9", 
    "CCDecid1", "CCDecid2", "CCDecid3", "CCDecid4", 
    "CCMixwood1", "CCMixwood2", "CCMixwood3", "CCMixwood4", 
    "CCConif1", "CCConif2", "CCConif3", "CCConif4", 
    "CCPine1", "CCPine2", "CCPine3", "CCPine4")
MixedDecidForest <- c(
    "Decid1", "Decid2", "Decid3", "Decid4", "Decid5", "Decid6", "Decid7", 
    "Decid8", "Decid9", "Mixwood1", "Mixwood2", 
    "Mixwood3", "Mixwood4", "Mixwood5", "Mixwood6", "Mixwood7", "Mixwood8", 
    "Mixwood9", 
    "Swamp-Decid1", "Swamp-Decid2", 
    "Swamp-Decid3", "Swamp-Decid4", "Swamp-Decid5", "Swamp-Decid6", 
    "Swamp-Decid7", "Swamp-Decid8", "Swamp-Decid9", 
    "Swamp-Mixwood1", "Swamp-Mixwood2", "Swamp-Mixwood3", 
    "Swamp-Mixwood4", "Swamp-Mixwood5", "Swamp-Mixwood6", "Swamp-Mixwood7", 
    "Swamp-Mixwood8", "Swamp-Mixwood9", 
    "CCDecid1", "CCDecid2", "CCDecid3", "CCDecid4", 
    "CCMixwood1", "CCMixwood2", "CCMixwood3", "CCMixwood4")
Wetland <- c("Swamp-Conif0", "Swamp-ConifR", "Swamp-Conif1", 
    "Swamp-Conif2", "Swamp-Conif3", "Swamp-Conif4", "Swamp-Conif5", 
    "Swamp-Conif6", "Swamp-Conif7", "Swamp-Conif8", "Swamp-Conif9", 
    "Swamp-Decid0", "Swamp-DecidR", "Swamp-Decid1", "Swamp-Decid2", 
    "Swamp-Decid3", "Swamp-Decid4", "Swamp-Decid5", "Swamp-Decid6", 
    "Swamp-Decid7", "Swamp-Decid8", "Swamp-Decid9", "Swamp-Mixwood0", 
    "Swamp-MixwoodR", "Swamp-Mixwood1", "Swamp-Mixwood2", "Swamp-Mixwood3", 
    "Swamp-Mixwood4", "Swamp-Mixwood5", "Swamp-Mixwood6", "Swamp-Mixwood7", 
    "Swamp-Mixwood8", "Swamp-Mixwood9", "Swamp-Pine0", "Swamp-PineR", 
    "Swamp-Pine1", "Swamp-Pine2", "Swamp-Pine3", "Swamp-Pine4", "Swamp-Pine5", 
    "Swamp-Pine6", "Swamp-Pine7", "Swamp-Pine8", "Swamp-Pine9", "Wetland-Bare", 
    "Wetland-GrassHerb", "Wetland-Shrub", "Wetland-BSpr0", "Wetland-BSprR", 
    "Wetland-BSpr1", "Wetland-BSpr2", "Wetland-BSpr3", "Wetland-BSpr4", 
    "Wetland-BSpr5", "Wetland-BSpr6", "Wetland-BSpr7", "Wetland-BSpr8", 
    "Wetland-BSpr9", "Wetland-Decid0", "Wetland-DecidR", "Wetland-Decid1", 
    "Wetland-Decid2", "Wetland-Decid3", "Wetland-Decid4", "Wetland-Decid5", 
    "Wetland-Decid6", "Wetland-Decid7", "Wetland-Decid8", "Wetland-Decid9", 
    "Wetland-Larch0", "Wetland-LarchR", "Wetland-Larch1", "Wetland-Larch2", 
    "Wetland-Larch3", "Wetland-Larch4", "Wetland-Larch5", "Wetland-Larch6", 
    "Wetland-Larch7", "Wetland-Larch8", "Wetland-Larch9")


hf <- dat[,c("PKEY","SS","YEAR","NSRNAME","NRNAME","POINT_X","POINT_Y",
    "xAHM","xPET","xFFP","xMAP","xMAT","xMCMT","xMWMT","xlong","xlat",
    "THF_KM", "Succ_KM", "Alien_KM", "THF_PC", "Succ_PC", "Alien_PC")]
hf$pDecMix <- rowSums(dd150m$veg_current[rownames(hf),MixedDecidForest])
hf$pForest <- rowSums(dd150m$veg_current[rownames(hf),ClosedCalopyForest])
hf$pWet <- rowSums(dd150m$veg_current[rownames(hf),Wetland])
hf$HFpc <- factor("Nv", c("Nv", "Am", "Ah", "Sm", "Sh"))
hf$HFpc[hf$Alien_PC >= 0.5 & hf$Succ_PC < 0.5] <- "Am"
hf$HFpc[hf$Succ_PC >= 0.5 & hf$Alien_PC < 0.5] <- "Sm"
hf$HFpc[hf$Alien_PC >= 0.75 & hf$Succ_PC < 0.5] <- "Ah"
hf$HFpc[hf$Succ_PC >= 0.75 & hf$Alien_PC < 0.5] <- "Sh"
table(succ=cut(hf$Succ_PC, c(-1, 0.5, 0.75, 2)), alien=cut(hf$Alien_PC, c(-1, 0.5, 0.75, 2)))
table(hf$HFpc)

save(hf, yy, tax, file="~/Dropbox/collaborations/opticut/R/abmi-birds-for-opticut.Rdata")

library(opticut)
yyy <- as.matrix(yy[,colSums(yy) > 10])
oc0p <- opticut(yyy ~ 1, strata=hf$HFpc, dist="poisson")
oc0b <- opticut(ifelse(yyy>0,1,0) ~ 1, strata=hf$HFpc, dist="binomial")
oc1p <- opticut(yyy ~ pForest + pWet + xMAT, hf, strata=hf$HFpc, dist="poisson")
oc1b <- opticut(ifelse(yyy>0,1,0) ~ pForest + pWet + xMAT, hf, strata=hf$HFpc, dist="binomial")

os0p <- summary(oc0p)$summary
os0b <- summary(oc0b)$summary
os1p <- summary(oc1p)$summary
os1b <- summary(oc1b)$summary

data.frame(spp=rownames(os0b), o0=os0b$split, o1=os1b$split)
data.frame(spp=rownames(os0b), b=os1b$split, p=os1p$split)
}
