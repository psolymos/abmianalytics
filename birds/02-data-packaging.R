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
source("~/repos/abmianalytics/birds/00-functions.R")
knitr::opts_chunk$set(eval=FALSE)
load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData")


## dealing with sine NAs

## XY: should probably drop
dd <- droplevels(dd[!is.na(dd$X),])
dd$MAXDUR[is.na(dd$MAXDUR)] <- 1 # ABMISM bits
## ROAD: needs HF info
## DATE/DATI: use constant sra
## climate: probably outside of AB bound: can use nearest if distance is small
## region: can intersect and use nearest -- it is not crucial
## check veg/hf/soil
##
aa <- data.frame(colSums(is.na(dd)))
dd <- droplevels(dd[!is.na(dd$pAspen),])

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


