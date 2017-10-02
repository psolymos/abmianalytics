library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
source("~/repos/abmianalytics/R/maps_functions.R")
regs <- levels(kgrid$LUFxNSR)

load(file.path(ROOT, "out", "transitions", paste0(regs[1], ".Rdata")))
Aveg <- rbind(colSums(trVeg))
rownames(Aveg) <- regs[1]
colnames(Aveg) <- colnames(trVeg)
Asoil <- rbind(colSums(trSoil))
rownames(Asoil) <- regs[1]
colnames(Asoil) <- colnames(trSoil)

for (i in 2:length(regs)) {
    cat(regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "transitions", paste0(regs[i], ".Rdata")))
    Aveg <- rbind(Aveg, colSums(trVeg))
    rownames(Aveg) <- regs[1:i]
    Asoil <- rbind(Asoil, colSums(trSoil))
    rownames(Asoil) <- regs[1:i]
}
Aveg <- Aveg / 10^4
Asoil <- Asoil / 10^4

## sector effect

tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv0$Sector2 <- factor(ifelse(is.na(tv0$Sector), "NATIVE", as.character(tv0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
tv <- droplevels(tv0[!is.na(tv0$Sector),])

ts0 <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")
ts0$All <- as.character(ts0$HF)
ts0$All[is.na(ts0$All)] <- as.character(ts0$Levels1)[is.na(ts0$All)]
ts0$Sector2 <- factor(ifelse(is.na(ts0$Sector), "NATIVE", as.character(ts0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
ts <- droplevels(ts0[!is.na(ts0$Sector),])

ch2veg <- t(sapply(strsplit(colnames(trVeg), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2veg <- data.frame(ch2veg)
colnames(ch2veg) <- c("rf","cr")
rownames(ch2veg) <- colnames(Aveg)
ch2veg$isHF <- ch2veg$cr %in% c("BorrowpitsDugoutsSumps",
    "Canals", "CCConif0", "CCConif1", "CCConif2", "CCConif3", "CCConif4",
    "CCConifR", "CCDecid0", "CCDecid1", "CCDecid2", "CCDecid3", "CCDecid4",
    "CCDecidR", "CCMixwood0", "CCMixwood1", "CCMixwood2", "CCMixwood3",
    "CCMixwood4", "CCMixwoodR", "CCPine0", "CCPine1", "CCPine2",
    "CCPine3", "CCPine4", "CCPineR",
    "CultivationCropPastureBareground", "HighDensityLivestockOperation",
    "IndustrialSiteRural",
    "MineSite",
    "MunicipalWaterSewage", "OtherDisturbedVegetation",
    "PeatMine", "Pipeline", "RailHardSurface",
    "RailVegetatedVerge", "Reservoirs", "RoadHardSurface", "RoadTrailVegetated",
    "RoadVegetatedVerge", "RuralResidentialIndustrial", "SeismicLine",
    "TransmissionLine", "Urban", "WellSite",
    "WindGenerationFacility")
compare_sets(ch2veg$cr, tv0$Combined)
ch2veg$Sector <- tv0$Sector2[match(ch2veg$cr, tv0$Combined)]
ch2veg$Sector[is.na(ch2veg$Sector)] <- "NATIVE" # strage Swamp thing
table(ch2veg$Sector,useNA="a")

ch2soil <- t(sapply(strsplit(colnames(trSoil), "->"),
    function(z) if (length(z)==1) z[c(1,1)] else z[1:2]))
ch2soil <- data.frame(ch2soil)
colnames(ch2soil) <- c("rf","cr")
rownames(ch2soil) <- colnames(Asoil)
ch2soil$isHF <- ch2soil$cr %in% c("BorrowpitsDugoutsSumps", "Canals",
    "CultivationCropPastureBareground",
    "CutBlocks", "HighDensityLivestockOperation", "IndustrialSiteRural",
    "MineSite", "MunicipalWaterSewage", "OtherDisturbedVegetation",
    "PeatMine", "Pipeline", "RailHardSurface", "RailVegetatedVerge",
    "Reservoirs", "RoadHardSurface", "RoadTrailVegetated",
    "RoadVegetatedVerge", "RuralResidentialIndustrial",
    "SeismicLine", "TransmissionLine",
    "Urban", "WellSite", "WindGenerationFacility")
compare_sets(ch2soil$cr, ts0$All)
ch2soil$Sector <- ts0$Sector2[match(ch2soil$cr, ts0$All)]
table(ch2soil$Sector,useNA="a")

lxn <- nonDuplicated(kgrid[,c("LUF_NAME","NRNAME","NSRNAME")], kgrid$LUFxNSR, TRUE)
lxn$N <- lxn$NRNAME != "Grassland" & lxn$NRNAME != "Rocky Mountain" &
    lxn$NRNAME != "Parkland" & lxn$NSRNAME != "Dry Mixedwood"
lxn$S <- lxn$NRNAME == "Grassland" | lxn$NRNAME == "Parkland" |
    lxn$NSRNAME == "Dry Mixedwood"
table(lxn$NRNAME, lxn$N)
table(lxn$NRNAME, lxn$S)
lxn <- lxn[regs,]
all(rownames(Aveg) == regs)
all(rownames(Asoil) == regs)

CN <- c("NATIVE", "Misc", "Agriculture", "RuralUrban", "Energy", "Transportation", "Forestry")

## check Area when defining a subregion !!! ------ !!!

AvegN <- colSums(Aveg[lxn$N,])
AvegN <- AvegN / sum(AvegN)
AsoilS <- colSums(Asoil[lxn$S,])
AsoilS <- AsoilS / sum(AsoilS)

TAX <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(TAX) <- TAX$AOU
SPP <- as.character(TAX$AOU)[TAX$map.pred]
EST <- TAX[,c("veghf.north", "soilhf.south")]
dimnames(EST) <- list(TAX$AOU, c("N", "S"))

fun_se <- function(POPS) {
    ## (cr-rf)/Ntot-in-region
    dN1n <- 100*(POPS[,"CRn"] - POPS[,"RFn"]) / sum(POPS[,"RFn"])
    dN1s <- 100*(POPS[,"CRs"] - POPS[,"RFs"]) / sum(POPS[,"RFs"])
    ## (cr-rf)/Ntot-in-HF
    dN2n <- 100*(POPS[,"CRn"] - POPS[,"RFn"]) / sum(POPS[rownames(POPS)!="NATIVE","RFn"])
    dN2s <- 100*(POPS[,"CRs"] - POPS[,"RFs"]) / sum(POPS[rownames(POPS)!="NATIVE","RFs"])
    ## (cr-rf)/Ntot-in-sectorHF
    dN3n <- 100*(POPS[,"CRn"] - POPS[,"RFn"]) / POPS[,"RFn"]
    dN3s <- 100*(POPS[,"CRs"] - POPS[,"RFs"]) / POPS[,"RFs"]
    f <- function(x) ifelse(is.na(x), 0, x)
    P1n <- 100*AAn/sum(AAn)
    P1s <- 100*AAs/sum(AAs)
    P2n <- 100*AAn/sum(AAn[names(AAn)!="NATIVE"])
    P2n["NATIVE"] <- 0
    P2s <- 100*AAs/sum(AAs[names(AAs)!="NATIVE"])
    P2s["NATIVE"] <- 0

    out <- list(
        north=data.frame(CR=POPS[,"CRn"], RF=POPS[,"RFn"],
            dN1=f(dN1n), dN2=f(dN2n), dN3=f(dN3n),
            P1=P1n, P2=P2n, U1=100*f(dN1n/P1n), U2=100*f(dN2n/P2n)),
        south=data.frame(CR=POPS[,"CRs"], RF=POPS[,"RFs"],
            dN1=f(dN1s), dN2=f(dN2s), dN3=f(dN3s),
            P1=P1s, P2=P2s, U1=100*f(dN1s/P1s), U2=100*f(dN2s/P2s)))
    rownames(out$north) <- rownames(out$south) <- rownames(POPS)
    out
}
ff <- function(a, x, y="north") {
    out <- t(sapply(a, function(z) z[[y]][c("Agriculture", "RuralUrban", "Energy",
        "Transportation", "Forestry"), x]))
    colnames(out) <- c("Agriculture", "RuralUrban", "Energy", "Transportation", "Forestry")
    out
}

## sector effects

PRED_DIR_IN <- "pred1-seismic-as-ES" # "pred1" # "pred1-shf"

AAn <- colSums(groupSums(Aveg, 2, ch2veg$Sector)[lxn$N,,drop=FALSE]) # in ha
AAs <- colSums(groupSums(Asoil, 2, ch2soil$Sector)[lxn$S,,drop=FALSE]) # in ha

#spp <- "BRCR"
Totals <- list()
for (spp in SPP) {
    fl <- list.files(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp))
    hbNcr <- hbNrf <- hbScr <- hbSrf <- Cells <- NULL
    for (i in 1:length(fl)) {
        cat(spp, i, "/", length(fl), "\n");flush.console()
        e <- new.env()
        load(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp, fl[i]), envir=e)
        hbNcr <- rbind(hbNcr, e$hbNcr1[,1])
        hbNrf <- rbind(hbNrf, e$hbNrf1[,1])
        hbScr <- rbind(hbScr, e$hbScr1[,1])
        hbSrf <- rbind(hbSrf, e$hbSrf1[,1])
        Cells <- c(Cells, e$Cells) # names(Cells) gives subset IDs
    }
    if (!EST[spp,"N"]) {
        hbNcr[] <- 0
        hbNrf[] <- 0
    }
    if (!EST[spp,"S"]) {
        hbScr[] <- 0
        hbSrf[] <- 0
    }
    regs2 <- gsub("\\.Rdata", "", fl)
    stopifnot(all(regs %in% regs2))
    ## NOTE: this Area definition only applies to the N/S version
    ## for subsets: you need subset specific areas
    ## Assemble North
    dimnames(hbNcr) <- dimnames(hbNrf) <- list(regs, colnames(Aveg))
    hbNcr[is.na(hbNcr)] <- 0
    hbNrf[is.na(hbNrf)] <- 0
    hbNcr <- hbNcr * Aveg[regs2,]
    hbNrf <- hbNrf * Aveg[regs2,]
    ## Assemble South
    dimnames(hbScr) <- dimnames(hbSrf) <- list(regs, colnames(Asoil))
    hbScr[is.na(hbScr)] <- 0
    hbSrf[is.na(hbSrf)] <- 0
    hbScr <- hbScr * Asoil[regs2,]
    hbSrf <- hbSrf * Asoil[regs2,]

    CRn <- colSums(groupSums(hbNcr, 2, ch2veg$Sector)[lxn$N,,drop=FALSE])
    RFn <- colSums(groupSums(hbNrf, 2, ch2veg$Sector)[lxn$N,,drop=FALSE])
    CRs <- colSums(groupSums(hbScr, 2, ch2soil$Sector)[lxn$S,,drop=FALSE])
    RFs <- colSums(groupSums(hbSrf, 2, ch2soil$Sector)[lxn$S,,drop=FALSE])
    POPS <- data.frame(
        CRn = CRn[CN],
        RFn = RFn[CN],
        CRs = CRs[CN],
        RFs = RFs[CN])
    rownames(POPS) <- CN
    Totals[[spp]] <- fun_se(POPS)
}

save(Totals, file=file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-NandS.Rdata"))

tab1 <- data.frame(Species=TAX[SPP, "species"],
    PopEffect=ff(Totals, "dN1", "north"),
    UnitEffect=ff(Totals, "U1", "north"),
    SingleSectorEffect=ff(Totals, "dN3", "north"))
tab1 <- tab1[EST[SPP,"N"],]
tab1 <- tab1[order(tab1$Species),]
write.csv(tab1, file=file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-North.csv"))

tab2 <- data.frame(Species=TAX[SPP, "species"],
    PopEffect=ff(Totals, "dN1", "south"),
    UnitEffect=ff(Totals, "U1", "south"),
    SingleSectorEffect=ff(Totals, "dN3", "south"))
tab2 <- tab2[EST[SPP,"S"],]
tab2 <- tab2[order(tab2$Species),]
write.csv(tab2, file=file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-South.csv"))

## same for OSA

ee <- new.env()
load("e:/peter/AB_data_v2016/out/transitions/00OSA All3.Rdata", envir=ee)
#e:/peter/AB_data_v2016/out/transitions/00OSA Athabasca Oilsand Area.Rdata
#e:/peter/AB_data_v2016/out/transitions/00OSA Cold Lake Oilsand Area.Rdata
#e:/peter/AB_data_v2016/out/transitions/00OSA Peace River Oilsand Area.Rdata

AAn <- colSums(groupSums(ee$trVeg, 2, ch2veg$Sector)) / 10^4 # in ha
AAs <- colSums(groupSums(ee$trSoil, 2, ch2soil$Sector)) / 10^4 # in ha
ATOT <- sum(AAn)
PRED_DIR_IN <- "pred1-seismic-as-ES-OSA" # "pred1" # "pred1-shf"

TotalsOSA <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    e <- new.env()
    load(file.path(ROOT, "out", "birds", PRED_DIR_IN, spp, "00OSA All3.Rdata"), envir=e)
    hbNcr <- e$hbNcr1[,1]
    hbNrf <- e$hbNrf1[,1]
    if (!EST[spp, "N"]) {
        hbNcr[] <- 0
        hbNrf[] <- 0
    }
    ## NOTE: this Area definition only applies to the N/S version
    ## for subsets: you need subset specific areas
    ## Assemble North
    hbNcr[is.na(hbNcr)] <- 0
    hbNrf[is.na(hbNrf)] <- 0
    hbNcr <- t(data.matrix(hbNcr * ATOT))
    hbNrf <- t(data.matrix(hbNrf * ATOT))

    CRn <- colSums(groupSums(hbNcr, 2, ch2veg$Sector))
    RFn <- colSums(groupSums(hbNrf, 2, ch2veg$Sector))
    POPS <- data.frame(
        CRn = CRn[CN],
        RFn = RFn[CN])
    rownames(POPS) <- CN
    POPS$CRs <- POPS$RFs <- 0
    TotalsOSA[[spp]] <- fun_se(POPS)[1]
}
save(TotalsOSA, file=file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-OSA.Rdata"))

tab3 <- data.frame(Species=TAX[SPP, "species"],
    PopEffect=ff(TotalsOSA, "dN1"),
    UnitEffect=ff(TotalsOSA, "U1"),
    SingleSectorEffect=ff(TotalsOSA, "dN3"))
tab3 <- tab3[EST[SPP,"N"],]
tab3 <- tab3[order(tab3$Species),]
write.csv(tab3, file=file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-OSA.csv"))

## making carrot plots

ROOT <- "e:/peter/AB_data_v2016"

x1 <- read.csv(file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-North.csv"))
rownames(x1) <- x1$Species
x2 <- read.csv(file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-South.csv"))
rownames(x2) <- x2$Species
x3 <- read.csv(file.path(ROOT, "out", "birds", "tables",
    "sector-effects-3kinds-seismic-as-ES-OSA.csv"))
rownames(x3) <- x3$Species


source("~/repos/abmianalytics/projects/10-yr/ch4-functions.R")

ylim=c(-100,100)
pdf(file.path(ROOT, "out", "birds", "tables", "Birds_SectorEffects.pdf"), height=18, width=12)
par(las=1, mfrow=c(3,2), yaxs="i")

tmp <- x1
tmp <- tmp[,grep("PopEffect", colnames(x1))]
colnames(tmp) <- gsub("PopEffect\\.", "", colnames(tmp))
vp(tmp, main="Sector Effects, Birds - North", ylim=ylim)
abline(h=0, lty=2)

tmp <- x1
tmp <- tmp[,grep("SingleSectorEffect", colnames(x1))]
colnames(tmp) <- gsub("SingleSectorEffect\\.", "", colnames(tmp))
vp(tmp, main="Single Sector Effects, Birds - North", ylim=ylim)
abline(h=0, lty=2)

tmp <- x2
tmp <- tmp[,grep("PopEffect", colnames(x2))]
colnames(tmp) <- gsub("PopEffect\\.", "", colnames(tmp))
tmp <- tmp[,colnames(tmp) != "Forestry"]
vp(tmp, main="Sector Effects, Birds - South", ylim=ylim)
abline(h=0, lty=2)

tmp <- x2
tmp <- tmp[,grep("SingleSectorEffect", colnames(x2))]
colnames(tmp) <- gsub("SingleSectorEffect\\.", "", colnames(tmp))
tmp <- tmp[,colnames(tmp) != "Forestry"]
vp(tmp, main="Single Sector Effects, Birds - South", ylim=ylim)
abline(h=0, lty=2)

tmp <- x3
tmp <- tmp[,grep("PopEffect", colnames(x3))]
colnames(tmp) <- gsub("PopEffect\\.", "", colnames(tmp))
vp(tmp, main="Sector Effects, Birds - OSA", ylim=ylim)
abline(h=0, lty=2)

tmp <- x3
tmp <- tmp[,grep("SingleSectorEffect", colnames(x3))]
colnames(tmp) <- gsub("SingleSectorEffect\\.", "", colnames(tmp))
vp(tmp, main="Single Sector Effects, Birds - OSA", ylim=ylim)
abline(h=0, lty=2)

dev.off()

