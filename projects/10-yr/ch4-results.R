library(mefa4)
ROOT <- "v:/contents/2017/species"
TAXA <- c("vplants", "mosses", "lichens", "mites")
AllIn <- list()
vegcols <- read.csv("c:/Users/Peter/repos/abmianalytics/projects/10-yr/veg-cols.csv")

TAX <- "mammals"

lt0 <- read.csv(file.path(ROOT, TAX, "lookup.csv"))
linn0 <- read.csv(file.path(ROOT, TAX, "lin10-north.csv"))
lins0 <- read.csv(file.path(ROOT, TAX, "lin10-south.csv"))
sectn0 <- read.csv(file.path(ROOT, TAX, "sector-north.csv"))
sects0 <- read.csv(file.path(ROOT, TAX, "sector-south.csv"))
soilnt0 <- read.csv(file.path(ROOT, TAX, "soilhf-nontreed-south.csv"))
soiltr0 <- read.csv(file.path(ROOT, TAX, "soilhf-treed-south.csv"))
usen0 <- read.csv(file.path(ROOT, TAX, "useavail-north.csv"))
uses0 <- read.csv(file.path(ROOT, TAX, "useavail-south.csv"))
veg0 <- read.csv(file.path(ROOT, TAX, "veghf-north.csv"))
sects0$PopEffect.Forestry <- sects0$UnitEffect.Forestry <- NULL
AllIn[["mammals"]] <- list(lt=lt0, usen=usen0, uses=uses0,
    veg=veg0, linn=linn0, soilnt=soilnt0, soiltr=soiltr0,
    lins=lins0, sectn=sectn0, sects=sects0)

TAX <- "birds"

lt0 <- read.csv(file.path(ROOT, TAX, "lookup.csv"))
linn0 <- read.csv(file.path(ROOT, TAX, "lin10-north.csv"))
lins0 <- read.csv(file.path(ROOT, TAX, "lin10-south.csv"))
sectn0 <- read.csv(file.path(ROOT, TAX, "sector-north.csv"))
sects0 <- read.csv(file.path(ROOT, TAX, "sector-south.csv"))
soilnt0 <- read.csv(file.path(ROOT, TAX, "soilhf-nontreed-south.csv"))
soiltr0 <- read.csv(file.path(ROOT, TAX, "soilhf-treed-south.csv"))
usen0 <- read.csv(file.path(ROOT, TAX, "useavail-north.csv"))
uses0 <- read.csv(file.path(ROOT, TAX, "useavail-south.csv"))
veg0 <- read.csv(file.path(ROOT, TAX, "veghf-north.csv"))
sects0$PopEffect.Forestry <- sects0$UnitEffect.Forestry <- NULL
AllIn[["birds"]] <- list(lt=lt0, usen=usen0, uses=uses0,
    veg=veg0, linn=linn0, soilnt=soilnt0, soiltr=soiltr0,
    lins=lins0, sectn=sectn0, sects=sects0)

#TAX <- "mites"
for (TAX in TAXA) {

lt <- read.csv(file.path(ROOT, TAX, "lookup.csv"))
linn <- read.csv(file.path(ROOT, TAX, "lin10-north.csv"))
lins <- read.csv(file.path(ROOT, TAX, "lin10-south.csv"))
paspen <- read.csv(file.path(ROOT, TAX, "paspen.csv"))
sectn <- read.csv(file.path(ROOT, TAX, "sector-north.csv"))
sects <- read.csv(file.path(ROOT, TAX, "sector-south.csv"))
veg <- read.csv(file.path(ROOT, TAX, "veghf-north.csv"))
soil <- read.csv(file.path(ROOT, TAX, "soilhf-south.csv"))
usen <- read.csv(file.path(ROOT, TAX, "useavail-north.csv"))
uses <- read.csv(file.path(ROOT, TAX, "useavail-south.csv"))

## lookup
SPP <- as.character(lt$Analysis_Name)
N <- length(SPP)
fl <- gsub(".png", "", list.files(file.path(ROOT, TAX, "map-det")))
compare_sets(SPP, fl)
lt1 <- with(lt, data.frame(
    SpeciesID=Analysis_Name,
    Species=Scientific_Name_Analysis,
    CommonName=NA,
    ScientificName=Scientific_Name_Analysis,
    SppCode=SppCode,
    TSNID=TSNID,
    nSites=nSites,
    nQuadrant=Occurrences,
    map.det=fl %in% SPP,
    veghf.north=Analysis.North,
    soilhf.south=Analysis.South,
    map.pred=Maps,
    useavail.north=UseAvailability.North,
    useavail.south=UseAvailability.South,
    origin=NA,
    comments=NA
))
if (TAX=="vplants")
    lt1$origin <- lt$Origin
lt1$useavail.north[lt1$veghf.north] <- FALSE
lt1$useavail.south[lt1$soilhf.south] <- FALSE
table(mod=lt1$veghf.north, use=lt1$useavail.north)
table(mod=lt1$soilhf.south, use=lt1$useavail.south)
stopifnot(sum(diag(table(modAll=lt1$veghf.north | lt1$soilhf.south, map=lt1$map.pred))) == N)

## veghf
veg1 <- veg
veg1 <- veg1[,colnames(veg1) %in% vegcols$In]
#colnames(veg1) <- vegcols$Out[match(colnames(veg1), vegcols$In)]
veg1 <- veg1[,match(vegcols$In, colnames(veg1))]
colnames(veg1) <- vegcols$Out
compare_sets(lt1$Species, veg1$Species)

## lin N
linn1 <- data.frame(Species=veg[,"Species"], linn,
    veg[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
colnames(linn1) <- c("Species", "AverageCoef", "SoftLin10", "HardLin10", "SoftLin",
    "SoftLin.LCL", "SoftLin.UCL", "HardLin", "HardLin.LCL", "HardLin.UCL")
compare_sets(lt1$Species, linn1$Species)

## soil, nontreed and treed
soilnt1 <- soil
soilnt1$X <- NULL
colnames(soilnt1) <- gsub(".LCI", ".LCL", colnames(soilnt1))
colnames(soilnt1) <- gsub(".UCI", ".UCL", colnames(soilnt1))
soilnt1 <- soilnt1[,colnames(soilnt0)]
soiltr1 <- data.frame(Species=soilnt1$Species,
    plogis(qlogis(as.matrix(soilnt1[,-1])) + paspen$pAspen))
compare_sets(lt1$Species, soilnt1$Species)

## lin S
lins1 <- data.frame(Species=soil[,"Species"], lins,
    soil[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
colnames(lins1) <- c("Species", "AverageCoef", "SoftLin10", "HardLin10", "SoftLin",
    "SoftLin.LCL", "SoftLin.UCL", "HardLin", "HardLin.LCL", "HardLin.UCL")
compare_sets(lt1$Species, lins1$Species)

## use-avail
usen1 <- usen
usen1$Species <- usen1$SpLabel
usen1$SpLabel <- NULL
usen1 <- usen1[,!grepl("\\.WRSI", colnames(usen1))]
colnames(usen1) <- gsub("\\.rWRSI", "", colnames(usen1))
#compare_sets(lt1$Species, usen1$Species)

uses1 <- uses
uses1$Species <- uses1$SpLabel
uses1$SpLabel <- NULL
uses1 <- uses1[,!grepl("\\.WRSI", colnames(uses1))]
colnames(uses1) <- gsub("\\.rWRSI", "", colnames(uses1))

## sector
sectn1 <- sectn
sectn1 <- data.frame(Species=sectn1$Sp, sectn1)
sectn1$Sp <- sectn1$pcTotalProvAbundance <- sectn1$Reg_Model <- NULL
colnames(sectn1) <- gsub("pcTotalEffect", "PopEffect", colnames(sectn1))
colnames(sectn1) <- gsub("pcUnitEffect", "UnitEffect", colnames(sectn1))
sectn1 <- sectn1[,colnames(sectn0)]
compare_sets(lt1$Species, sectn1$Species)

## note: this does not have Forestry (Grassland!)
sects1 <- sects
sects1 <- data.frame(Species=sects1$Sp, sects1)
sects1$Sp <- sects1$pcTotalProvAbundance <- sects1$Reg_Model <- NULL
#sects1$pcTotalEffect.Forestry <- sects1$pcUnitEffect.Forestry <- 0
#sects1 <- sects1[,colnames(sectn1)]
colnames(sects1) <- gsub("pcTotalEffect", "PopEffect", colnames(sects1))
colnames(sects1) <- gsub("pcUnitEffect", "UnitEffect", colnames(sects1))
sects1 <- sects1[,colnames(sects0)]
compare_sets(lt1$Species, sects1$Species)

AllIn[[TAX]] <- list(lt=lt1, usen=usen1, uses=uses1,
    veg=veg1, linn=linn1, soilnt=soilnt1, soiltr=soiltr1,
    lins=lins1, sectn=sectn1, sects=sects1)
}
save(AllIn, file="~/Dropbox/abmi/10yr/R/AllTables.Rdata")

## summarizing intactness

library(mefa4)

load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))

load("~/Dropbox/abmi/10yr/R/AllTables.Rdata")
ROOT <- "v:/contents/2017/species"
TAXA <- c("mammals", "birds", "vplants", "mosses", "lichens", "mites")

AllSI <- list()
MatSI <- matrix(NA, nrow(kgrid), length(TAXA))
dimnames(MatSI) <- list(rownames(kgrid), TAXA)

TAX <- "mammals"
fl <- list.files(file.path(ROOT, TAX, "km2"))
nam <- sapply(strsplit(fl, "\\."), "[[", 1)
ext <- strsplit(fl, "\\.")[[1]][2]
lt <- AllIn[[TAX]]$lt
SPP <-as.character(lt$SpeciesID[lt$map.pred])
compare_sets(nam, SPP)
sppSI <- list()
simat <- matrix(NA, nrow(kgrid), length(nam))
dimnames(simat) <- list(rownames(kgrid), nam)
for (i in nam) {
    cat(TAX, i, "\n");flush.console()
    km <- read.csv(file.path(ROOT, TAX, "km2", paste0(i, ".", ext)))
    km <- km[match(rownames(kgrid), km$LinkID),c("Ref", "Curr")]
    SI <- 100 * pmin(km$Curr, km$Ref) / pmax(km$Curr, km$Ref)
    kmnr <- groupSums(as.matrix(km), 1, kgrid$NRNAME, na.rm=TRUE)
    kmnr <- rbind(kmnr, Alberta=colSums(kmnr))
    kmnr <- cbind(kmnr, SI=100*pmin(kmnr[,1], kmnr[,2])/pmax(kmnr[,1], kmnr[,2]))
    sppSI[[i]] <- kmnr
    simat[,i] <- SI
}
SIavg <- rowMeans(simat, na.rm=TRUE)
#SIavg[is.na(SIavg)] <- 100
AllSI[[TAX]] <- sppSI
MatSI[,TAX] <- SIavg

TAX <- "birds"
fl <- list.files(file.path(ROOT, TAX, "km2"))
nam <- sapply(strsplit(fl, "\\."), "[[", 1)
ext <- strsplit(fl, "\\.")[[1]][2]
lt <- AllIn[[TAX]]$lt
SPP <-as.character(lt$AOU[lt$map.pred])
compare_sets(nam, SPP)
sppSI <- list()
simat <- matrix(NA, nrow(kgrid), length(nam))
dimnames(simat) <- list(rownames(kgrid), nam)
for (i in nam) {
    cat(TAX, i, "\n");flush.console()
    e <- new.env()
    load(file.path(ROOT, TAX, "km2", paste0(i, ".", ext)), envir=e)
    km <- as.data.frame(e$km2)
    SI <- 100 * pmin(km$Curr, km$Ref) / pmax(km$Curr, km$Ref)
    kmnr <- groupSums(as.matrix(km), 1, kgrid$NRNAME, na.rm=TRUE)
    kmnr <- rbind(kmnr, Alberta=colSums(kmnr))
    kmnr <- cbind(kmnr, SI=100*pmin(kmnr[,1], kmnr[,2])/pmax(kmnr[,1], kmnr[,2]))
    sppSI[[i]] <- kmnr
    simat[,i] <- SI
}
SIavg <- rowMeans(simat, na.rm=TRUE)
#SIavg[is.na(SIavg)] <- 100
AllSI[[TAX]] <- sppSI
MatSI[,TAX] <- SIavg

if (FALSE) {
ByNSR <- list()
for (i in nam) {
    cat(TAX, i, "\n");flush.console()
    e <- new.env()
    load(file.path(ROOT, TAX, "km2", paste0(i, ".", ext)), envir=e)
    km <- as.data.frame(e$km2)
    kmnr <- groupSums(as.matrix(km), 1, kgrid$NSRNAME, na.rm=TRUE)
    ByNSR[[i]] <- 100*kmnr # (D per ha * 100 = inds/km^2)
}
ToSave <- do.call(rbind, lapply(1:length(ByNSR), function(z) {
    data.frame(Species=names(ByNSR)[z], NSR=rownames(ByNSR[[z]]), N_cr=ByNSR[[z]][,"Curr"],
    N_rf=ByNSR[[z]][,"Ref"])
}))
write.csv(ToSave,row.names=FALSE,file="Alberta-PopSizes-by-NSR.csv")

en <- new.env()
load(file.path("e:/peter/AB_data_v2016/out/birds", "data", "data-north.Rdata"), envir=en)

e <- new.env()
load(file.path("e:/peter/AB_data_v2016/out/birds", "data", "data-wrsi.Rdata"), envir=e)
TAX <- droplevels(e$TAX)

zz <- read.csv("Alberta-PopSizes-by-NSR.csv")
zz$HasOffset <- zz$Species %in% colnames(en$OFF)
zz <- data.frame(zz,
    TAX[match(zz$Species, rownames(TAX)),c("English_Name", "Scientific_Name",
    "Family_Sci", "Order")])
write.csv(zz,row.names=FALSE,file="Alberta-PopSizes-by-NSR_wStuff.csv")

}

#TAX <- "mites"
for (TAX in TAXA[-(1:2)]) {
    fl <- list.files(file.path(ROOT, TAX, "km2"))
    nam <- substr(fl, 1, nchar(fl)-6)
    ext <- substr(fl[1], nchar(fl)-4, nchar(fl))
    lt <- AllIn[[TAX]]$lt
    SPP <-as.character(lt$SpeciesID[lt$map.pred])
    if (TAX=="vplants") {
        Native <- !is.na(lt$origin) & lt$origin != "Exotic"
        SPP <-as.character(lt$SpeciesID[lt$map.pred & Native])
    }
    compare_sets(nam, SPP)
    sppSI <- list()
    simat <- matrix(NA, nrow(kgrid), length(nam))
    dimnames(simat) <- list(rownames(kgrid), nam)
    for (i in nam) {
        cat(TAX, i, "\n");flush.console()
        e <- new.env()
        load(file.path(ROOT, TAX, "km2", paste0(i, ".", ext)), envir=e)
        km <- as.data.frame(e$RefCurr)
        km <- km[match(rownames(kgrid), km$LinkID),c("Ref", "Curr")]
        SI <- 100 * pmin(km$Curr, km$Ref) / pmax(km$Curr, km$Ref)
        kmnr <- groupSums(as.matrix(km), 1, kgrid$NRNAME, na.rm=TRUE)
        kmnr <- rbind(kmnr, Alberta=colSums(kmnr))
        kmnr <- cbind(kmnr, SI=100*pmin(kmnr[,1], kmnr[,2])/pmax(kmnr[,1], kmnr[,2]))
        sppSI[[i]] <- kmnr
        simat[,i] <- SI
    }
    SIavg <- rowMeans(simat, na.rm=TRUE)
    #SIavg[is.na(SIavg)] <- 100
    AllSI[[TAX]] <- sppSI
    MatSI[,TAX] <- SIavg
}

save(AllSI, MatSI, file="e:/peter/sppweb2017/all-intactness.Rdata")

## summarizing GoF measures

coefbsn <- coefbss <- list()

load("e:/peter/sppweb2017/birds-r2-auc-occc-boot-subset.Rdata")
cfn <- array(NA, c(length(boot_res_N), nrow(boot_res_N[[1]]), 100))
cfs <- array(NA, c(length(boot_res_S), nrow(boot_res_S[[1]]), 100))
dimnames(cfn) <- list(names(boot_res_N), rownames(boot_res_N[[1]]), 1:100)
dimnames(cfs) <- list(names(boot_res_S), rownames(boot_res_S[[1]]), 1:100)
for (i in 1:length(boot_res_N))
    cfn[i,,] <- boot_res_N[[i]][,1:100]
for (i in 1:length(boot_res_S))
    cfs[i,,] <- boot_res_S[[i]][,1:100]
coefbsn$birds <- cfn
coefbss$birds <- cfs

e <- new.env()
load("w:/All-In_One_Ver2017/Mites/Analysis north/R object Mite intactness North by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbsn$mites <- e$Coef.bs
e <- new.env()
load("w:/All-In_One_Ver2017/Mites/Analysis south/R object Mite intactness South by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbss$mites <- e$Coef.bs

e <- new.env()
load("w:/All-In_One_Ver2017/Lichen/Analysis north/R object Lichen intactness North by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbsn$lichens <- e$Coef.bs
e <- new.env()
load("w:/All-In_One_Ver2017/Lichen/Analysis south/R object Lichen intactness South by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbss$lichens <- e$Coef.bs

e <- new.env()
load("w:/All-In_One_Ver2017/Moss/Analysis north/R object Moss intactness North by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbsn$mosses <- e$Coef.bs
e <- new.env()
load("w:/All-In_One_Ver2017/Moss/Analysis south/R object Moss intactness South by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbss$mosses <- e$Coef.bs

e <- new.env()
load("w:/All-In_One_Ver2017/VPlant/Analysis north/R object Plant intactness North by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbsn$vplants <- e$Coef.bs
e <- new.env()
load("w:/All-In_One_Ver2017/VPlant/Analysis south/R object Plant intactness South by 10x10km unit BS 2017 10Km2.Rdata", envir=e)
coefbss$vplants <- e$Coef.bs

rm(e)
sapply(coefbsn, dim)
sapply(coefbss, dim)

x <- coefbsn$birds[1,,]

occc_plot <- function(x, ...) {
    require(epiR)
    x <- x[!(rownames(x) %in% c("SoftLin","HardLin" )),,drop=FALSE]
    z <- t(x)
    oc <- epi.occc(x)
    cm <- apply(z, 2, median)
    z <- z[,order(cm)]
    z <- z / max(cm)
    q <- apply(z, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))
    matplot(t(q), type="l", lty=1,
        col=c("lightblue", "blue", "black", "blue", "lightblue"), axes=FALSE,
        xlab=paste0("OCCC=", round(oc$occc,2), ", Location=", round(oc$oprec, 2),
        ", Scale=", round(oc$oaccu, 2)),
        ylab="Relative abundance", ...)
    axis(2)
    abline(h=1, col="grey")
    invisible(NULL)
}

oc <- t(sapply(1:dim(coefbsn$vplants)[1], function(i) unlist(epi.occc(coefbsn$vplants[i,,])[1:3])))
rownames(oc) <- dimnames(coefbsn$vplants)[[1]]
head(oc[order(oc[,1]),])
head(oc[order(oc[,2]),])
head(oc[order(oc[,3]),])
tail(oc[order(oc[,1]),])

oc <- t(sapply(1:dim(coefbsn$birds)[1], function(i) unlist(epi.occc(coefbsn$birds[i,,])[1:3])))
rownames(oc) <- dimnames(coefbsn$birds)[[1]]
head(oc[order(oc[,1]),])
head(oc[order(oc[,2]),])
head(oc[order(oc[,3]),])
tail(oc[order(oc[,1]),])

pdf("e:/peter/sppweb2017/occc_fig_birds_north.pdf", onefile=TRUE)
for (i in 1:dim(coefbsn$birds)[1])
    occc_plot(coefbsn$birds[i,,], main=dimnames(coefbsn$birds)[[1]][i])
dev.off()

pdf("e:/peter/sppweb2017/occc_fig_birds_south.pdf", onefile=TRUE)
for (i in 1:dim(coefbss$birds)[1])
    occc_plot(coefbss$birds[i,,], main=dimnames(coefbss$birds)[[1]][i])
dev.off()

pdf("e:/peter/sppweb2017/occc_fig_vplants_north.pdf", onefile=TRUE)
for (i in 1:dim(coefbsn$vplants)[1])
    occc_plot(coefbsn$vplants[i,,], main=dimnames(coefbsn$vplants)[[1]][i])
dev.off()

pdf("e:/peter/sppweb2017/occc_fig_vplants_south.pdf", onefile=TRUE)
for (i in 1:dim(coefbss$vplants)[1])
    occc_plot(coefbss$vplants[i,,], main=dimnames(coefbss$vplants)[[1]][i])
dev.off()

gofs <- gofn <- list()
gofn$birds <- cbind(all_acc_N, t(sapply(occc_res_N, function(z) unlist(z[1:3]))), Wdist=NA)
gofs$birds <- cbind(all_acc_S, t(sapply(occc_res_S, function(z) unlist(z[1:3]))), Wdist=NA)

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Mites/North/AUC, pseudoR2, and coeffcient of discrimination_Mites_North.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Mites/North/Repeatability measures for habitat coefficients_Mites_North.csv")
gofn$mites <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Mites/South/AUC, pseudoR2, and coeffcient of discrimination_Mites_South.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Mites/South/Repeatability measures for habitat coefficients_Mites_South.csv")
gofs$mites <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Lichen/North/AUC, pseudoR2, and coeffcient of discrimination_Lichen_North.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Lichen/North/Repeatability measures for habitat coefficients_Lichen_North.csv")
gofn$lichens <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Lichen/South/AUC, pseudoR2, and coeffcient of discrimination_Lichen_South.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Lichen/South/Repeatability measures for habitat coefficients_Lichen_South.csv")
gofs$lichens <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Moss/North/AUC, pseudoR2, and coeffcient of discrimination_Moss_North.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Moss/North/Repeatability measures for habitat coefficients_Moss_North.csv")
gofn$mosses <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Moss/South/AUC, pseudoR2, and coeffcient of discrimination_Moss_South.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/Moss/South/Repeatability measures for habitat coefficients_Moss_South.csv")
gofs$mosses <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/VPlant/North/AUC, pseudoR2, and coeffcient of discrimination_VPlant_North.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/VPlant/North/Repeatability measures of habitat coefficients_VPlant_North.csv")
gofn$vplants <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

t1 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/VPlant/South/AUC, pseudoR2, and coeffcient of discrimination_VPlant_South.csv")
t2 <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables for 10 Year Review/VPlant/South/Repeatability measures for habitat coefficients_VPlant_South.csv")
gofs$vplants <- as.matrix(cbind(t1[,c("PseudoR2_VegHF", "PseudoR2_VegHF.ClimSpace", "AUC_Veg", "AUC_ALL")],
    t2[,c("occc", "oprec", "oaccu", "Within.Agreement")]))

save(gofn, gofs, coefbsn, coefbss, file="~/Dropbox/abmi/10yr/R/AllBoot.Rdata")

load("~/Dropbox/abmi/10yr/R/AllBoot.Rdata")
if (FALSE) {
summary(gofn$birds)
summary(gofs$birds)
wan <- read.csv("~/Dropbox/abmi/10yr/R/Repeatability measures_Birds_North.csv")
was <- read.csv("~/Dropbox/abmi/10yr/R/Repeatability measures_Birds_South.csv")
rownames(wan) <- wan$X
rownames(was) <- was$X
gofn$birds[intersect(wan$X,rownames(gofn$birds)),"Wdist"] <- wan[
    intersect(wan$X,rownames(gofn$birds)), "Within.Agreement"]
gofs$birds[intersect(was$X,rownames(gofs$birds)),"Wdist"] <- was[
    intersect(was$X,rownames(gofs$birds)), "Within.Agreement"]
summary(gofn$birds)
summary(gofs$birds)
save(gofn, gofs, coefbsn, coefbss, file="~/Dropbox/abmi/10yr/R/AllBoot.Rdata")
}

gof_plot <- function(x, type=c("r2", "auc", "oc", "d"), alpha=1, ...) {
    f <- function(x) {
        x <- x-min(x)
        x <- x/max(x)
        x
    }
    x <- x[rowSums(is.na(x[,1:7]))==0,]
    if (type=="r2") {
        x[x[,1] < 0,1] <- 0
        x[x[,2] < 0,2] <- 0
        plot(x[,1:2], col=rgb(f(x[,2]-x[,1]), 1-f(x[,1]), 1-f(x[,2]), alpha=alpha), pch=19,
        ylim=c(0,1), xlim=c(0,1), xlab=expression(R^2*(landcover)),
        ylab=expression(R^2*(landcover+climate)), ...)
        abline(0,1,col="grey")
    }
    if (type=="auc") {
        plot(x[,3:4], col=rgb(f(x[,4]-x[,3]), 1-f(x[,3]), 1-f(x[,4]), alpha=alpha), pch=19,
        ylim=c(0,1), xlim=c(0,1), xlab="AUC (landcover)",
        ylab="AUC (landcover+climate)", ...)
        abline(0,1,col="grey")
        abline(h=0.5, v=0.5,col="grey")
    }
    if (type=="oc") {
        plot(x[,6:7], col=rgb(f(x[,5]), 1-f(x[,6]), 1-f(x[,7]), alpha=alpha), pch=19,
        ylim=c(0,1), xlim=c(0,1), xlab="Location",
        ylab="Scale", ...)
        abline(0,1,col="grey")
    }
    if (type=="d") {
        plot(density(x[,8]), xlim=c(0,1), ...)
    }
    invisible(NULL)
}

gof_fig_all <- function(TAX, ...) {
    if (TAX=="All") {
        xn <- do.call(rbind, gofn)
        xs <- do.call(rbind, gofs)
    } else {
        xn <- gofn[[TAX]]
        xs <- gofs[[TAX]]
    }
    op <- par(mfrow=c(2,3), las=1, mar=c(5,5,2,2)+0.1)
    gof_plot(xn, "r2", main=paste(TAX, "north"), ...)
    gof_plot(xn, "auc", ...)
    gof_plot(xn, "oc", ...)
    gof_plot(xs, "r2", main=paste(TAX, "south"), ...)
    gof_plot(xs, "auc", ...)
    gof_plot(xs, "oc", ...)
    par(op)
    invisible(NULL)
}
pdf("e:/peter/sppweb2017/r2-auc-oc.pdf", onefile=TRUE, height=8, width=12)
for (i in c("All", names(gofn)))
    gof_fig_all(i, alpha=0.5)
dev.off()

## making sense

library(mefa4)
library(intrval)
load("~/Dropbox/abmi/10yr/R/AllTables.Rdata")
AllIn$vplants$lt$pOcc <- AllIn$vplants$lt$nSites / 1598

get_stuff0 <- function(TAX, TABLE, COL, north=TRUE) {
    SPP <- AllIn[[TAX]]$lt$Species
    i <- if (north)
        AllIn[[TAX]]$lt[,"veghf.north"] else AllIn[[TAX]]$lt[,"soilhf.south"]
    SPP <- as.character(SPP[i])
    z <- AllIn[[TAX]][[TABLE]]
    rownames(z) <- z$Species
    #compare_sets(rownames(z), SPP)
    z[SPP, COL]
}
rugged_mat <- function(x) {
    n <- sapply(x, length)
    out <- matrix(NA, max(n), length(x))
    colnames(out) <- names(x)
    for (i in 1:length(x))
        out[1:length(x[[i]]),i] <- x[[i]]
    out
}
get_stuff <- function(TABLE, COL, north=TRUE)
    rugged_mat(lapply(structure(names(AllIn), names=names(AllIn)),
        get_stuff0, TABLE=TABLE, COL=COL, north=north))
bounded_density <- function(x, interval=c(-Inf, Inf), ...) {
    require(intrval)
    interval <- sort(interval[1L:2L])
    if (!all(x %()% interval))
        stop("found values of x outside of open interval")
    a <- interval[1L]
    b <- interval[2L]
    bounds  <- is.finite(interval)
    if (!bounds[1L] && !bounds[2L]) { # (Inf, Inf)
        f <- finv <- function(x) x
    }
    if (bounds[1L] && !bounds[2L]) { # (a, Inf)
        f <- function(x)
            log(x-a)
        finv <- function(z)
            exp(z) + a
    }
    if (!bounds[1L] && bounds[2L]) { # (Inf, b)
        f <- function(x)
            log(b-x)
        finv <- function(z)
            b-exp(z)
    }
    if (bounds[1L] && bounds[2L]) { # (a, b)
        f <- function(x)
            qlogis((x-a) / (b-a))
        finv <- function(z)
            plogis(z) * (b-a) + a
    }
    fx <- f(x)
    d <- density(fx, ...)
    v <- d$x
    dv <- diff(v)
    h <- d$y
    n <- length(h)
    dh <- rowMeans(cbind(h[-n], h[-1L]))
    A <- dv * dh
    vinv <- finv(v)
    hinv <- A / abs(diff(vinv))
    hinv <- c(hinv[1L], hinv)
    data.frame(x=vinv, y=hinv)
}

if (FALSE) {
x <- rnorm(100)

op <- par(mfrow=c(2,2))
range(x)
d1 <- bounded_density(exp(x))
d2 <- bounded_density(exp(x), c(0, Inf))
plot(bounded_density(x, type="l"))

range(exp(x))
d1 <- bounded_density(exp(x))
d2 <- bounded_density(exp(x), c(0, Inf))
plot(d1, type="l", ylim=c(range(d1$y, d2$y)))
lines(d2, col=2)

range(plogis(x)*4-1)
d1 <- bounded_density(plogis(x)*4-1)
d2 <- bounded_density(plogis(x)*4-1, c(-1, 3))
plot(d1, type="l", ylim=c(range(d1$y, d2$y)))
lines(d2, col=2)

range(-exp(x)+3)
d1 <- bounded_density(-exp(x)+3)
d2 <- bounded_density(-exp(x)+3, c(-Inf, 3))
plot(d1, type="l", ylim=c(range(d1$y, d2$y)))
lines(d2, col=2)
par(op)
}

cp1 <- function(x, p=c(0, 1), nmin=5, interval=c(-Inf, Inf), ...) {
    x2 <- x[!is.na(x)]
    q <- quantile(x2, p)
    x2 <- x2[x2 > q[1] & x2 < q[2]]
    if (length(unique(x2)) < nmin) {
        m <- mean(x2)
        s <- sd(x2)
        v <- seq(min(min(x2), m-4*s), max(max(x2), m+4*s), length.out = 512)
        v <- v %()% interval
        d <- dnorm(v, m, s)
        v <- c(v[1L], v, v[length(v)])
        d <- c(0, d, 0)
        out <- data.frame(h=v, w=d/max(d))
    } else {
        d <- bounded_density(x2, interval=interval, ...)
        out <- data.frame(h=d$x, w=d$y/max(d$y))
    }
    out
}
vp <- function(x, p=c(0,1),
col="#FF664D80", border="#FF664D",
#col="#40E0D180", border="#40E0D1",
ylim, ylab="", xlab="", main="", nmin=5, interval=c(-Inf, Inf), ...) {
    xx <- apply(x, 2, cp1, p=p, nmin=nmin, interval=interval, ...)
    if (missing(ylim))
        ylim <- range(x, na.rm=TRUE)
    plot(0, type="n", axes=FALSE, ann=FALSE, xlim=c(0.5, ncol(x)+0.5),
        ylim=ylim)
    axis(1, 1:ncol(x), colnames(x), lwd=0)
    axis(2)
    title(ylab=ylab, xlab=xlab, main=main)
    col <- rep(col, ncol(x))[1:ncol(x)]
    border <- rep(border, ncol(x))[1:ncol(x)]
    for (i in 1:ncol(x)) {
        if (!is.null(xx[[i]])) {
            polygon(0.45*c(-xx[[i]]$w, rev(xx[[i]]$w))+i,
                c(xx[[i]]$h, rev(xx[[i]]$h)), col=col[i], border=border[i])
            #lines(c(i,i), range(xx[[i]]$h), col=col[i], lwd=2)
        }
        lines(c(i-0.2, i+0.2), rep(median(x[,i], na.rm=TRUE), 2), lwd=3)
    }
    #points(1:ncol(x), colMeans(x, na.rm=TRUE), pch=21, cex=1.5)
    invisible(NULL)
}
fu <- function(x) sign(x) * plogis(log(abs(x/100)))

colTd <- RColorBrewer::brewer.pal(6, "Dark2")
colTl <- paste0(colTd, "80")

## --- linear

softn <- get_stuff("linn", "SoftLin10", TRUE) /
    get_stuff("linn", "AverageCoef", TRUE)
hardn <- get_stuff("linn", "HardLin10", TRUE) /
    get_stuff("linn", "AverageCoef", TRUE)
softs <- get_stuff("lins", "SoftLin10", FALSE) /
    get_stuff("lins", "AverageCoef", FALSE)
hards <- get_stuff("lins", "HardLin10", FALSE) /
    get_stuff("lins", "AverageCoef", FALSE)
softn[softn > 10] <- NA
softs[softs > 10] <- NA
hardn[hardn > 10] <- NA
hards[hards > 10] <- NA

## might need to calculate mean density for birds instead of current AvgCoefficient
pdf("~/Dropbox/abmi/10yr/ch4/figs/lin.pdf", height=10, width=10)
par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(0, 3)
p <- c(0.025, 0.975)
vp(softn, main="Soft Linear, North", ylim=ylim, p=p, ylab="Std. Effect", col=colTl, border=colTd)
abline(h=1, lty=2)
vp(hardn, main="Hard Linear, North", ylim=ylim, p=p, ylab="Std. Effect", col=colTl, border=colTd)
abline(h=1, lty=2)
vp(softs, main="Soft Linear, South", ylim=ylim, p=p, ylab="Std. Effect", col=colTl, border=colTd)
abline(h=1, lty=2)
vp(hards, main="Hard Linear, South", ylim=ylim, p=p, ylab="Std. Effect", col=colTl, border=colTd)
abline(h=1, lty=2)
dev.off()

## --- sector

tnagr <- get_stuff("sectn", "PopEffect.Agriculture", TRUE)
tnfor <- get_stuff("sectn", "PopEffect.Forestry", TRUE)
tneng <- get_stuff("sectn", "PopEffect.Energy", TRUE)
tnurb <- get_stuff("sectn", "PopEffect.RuralUrban", TRUE)
tntra <- get_stuff("sectn", "PopEffect.Transportation", TRUE)

tsagr <- get_stuff("sects", "PopEffect.Agriculture", FALSE)
tsfor <- get_stuff("sects", "PopEffect.Forestry", FALSE)
tseng <- get_stuff("sects", "PopEffect.Energy", FALSE)
tsurb <- get_stuff("sects", "PopEffect.RuralUrban", FALSE)
tstra <- get_stuff("sects", "PopEffect.Transportation", FALSE)

unagr <- get_stuff("sectn", "UnitEffect.Agriculture", TRUE)
unfor <- get_stuff("sectn", "UnitEffect.Forestry", TRUE)
uneng <- get_stuff("sectn", "UnitEffect.Energy", TRUE)
unurb <- get_stuff("sectn", "UnitEffect.RuralUrban", TRUE)
untra <- get_stuff("sectn", "UnitEffect.Transportation", TRUE)

usagr <- get_stuff("sects", "UnitEffect.Agriculture", FALSE)
usfor <- get_stuff("sects", "UnitEffect.Forestry", FALSE)
useng <- get_stuff("sects", "UnitEffect.Energy", FALSE)
usurb <- get_stuff("sects", "UnitEffect.RuralUrban", FALSE)
ustra <- get_stuff("sects", "UnitEffect.Transportation", FALSE)

pdf("~/Dropbox/abmi/10yr/ch4/figs/sect1.pdf", height=10, width=10)
par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(-1, 1)
p <- c(0, 1)
INT <- c(-1,1)
vp(fu(uneng), main="Energy, North", ylim=ylim, p=p, ylab="Unit Effect",
    interval=INT, col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unagr), main="Agriculture, North", ylim=ylim, p=p, ylab="Unit Effect",
    interval=INT, col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unurb), main="Rural/Urban, North", ylim=ylim, p=p, ylab="Unit Effect",
    interval=INT, col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(untra), main="Transportation, North", ylim=ylim, p=p, ylab="Unit Effect",
    interval=INT, col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unfor), main="Forestry, North", ylim=ylim, p=p, ylab="Unit Effect",
    interval=INT, col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
plot.new()
dev.off()

pdf("~/Dropbox/abmi/10yr/ch4/figs/sect2.pdf", height=10, width=10)
par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(-1, 1)
p <- c(0, 1)
vp(fu(useng), main="Energy, South", ylim=ylim, p=p, ylab="Unit Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(usagr), main="Agriculture, South", ylim=ylim, p=p, ylab="Unit Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(usurb), main="Rural/Urban, South", ylim=ylim, p=p, ylab="Unit Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(ustra), main="Transportation, South", ylim=ylim, p=p, ylab="Unit Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
dev.off()


pdf("~/Dropbox/abmi/10yr/ch4/figs/sect3.pdf", height=10, width=10)
par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(-10, 10)
p <- c(0.025, 0.975)
vp(tneng, main="Energy, North", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tnagr, main="Agriculture, North", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tnurb, main="Rural/Urban, North", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tntra, main="Transportation, North", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tnfor, main="Forestry, North", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
plot.new()
dev.off()

pdf("~/Dropbox/abmi/10yr/ch4/figs/sect4.pdf", height=10, width=10)
par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(-20, 20)
p <- c(0.025, 0.975)
vp(tseng, main="Energy, South", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tsagr, main="Agriculture, South", ylim=c(-100, 100), p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tsurb, main="Rural/Urban, South", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
vp(tstra, main="Transportation, South", ylim=ylim, p=p, ylab="Total Effect", col=colTl, border=colTd)
abline(h=0, lty=2)
dev.off()

## --- ordination

ord_fun <- function(TAX, north=TRUE, scaling=2, alpha=1, col.text=1, col.pts,
    cex.pts=1, cex.text=1, ...) {
    require(vegan)
    Excl <- c("WhiteSpruceCC20", "WhiteSpruceCC40", "WhiteSpruceCC60",
        "PineCC20", "PineCC40", "PineCC60",
        "DeciduousCC20",   "DeciduousCC40", "DeciduousCC60",
        "MixedwoodCC20", "MixedwoodCC40",   "MixedwoodCC60")
    if (TAX=="All") {
        TT <- if (north)
            "veg" else "soilnt"
        z <- do.call(rbind, lapply(AllIn, function(zz) zz[[TT]]))
    } else {
        z <- if (north)
            AllIn[[TAX]]$veg else AllIn[[TAX]]$soilnt
    }
    rownames(z) <- z[,1]
    z <- as.matrix(z[,-1])
    z <- z[,!grepl("\\.", colnames(z))]
    z <- z[,!(colnames(z) %in% Excl)]
    cn <- colnames(z)
    cn <- gsub("[[:digit:]]", "", cn)
    #cn <- gsub("CC", "", cn)
    zz <- groupMeans(z, 2, cn)
    #rownames(zz) <- make.cepnames(rownames(zz)) # species
    colnames(zz) <- make.cepnames(colnames(zz)) # habitats

    fs <- function(x) 0.5*(1+x/max(abs(x)))
    yy <- (t(zz))
    m <- cca(yy)
    s <- scores(m, 1:3)
    ColSi <- 2
    ColSp <- rgb(red=fs(s$species[,1]),
        green=fs(s$species[,2]), blue=fs(s$species[,3]), alpha=alpha)
    if (!missing(col.pts))
        ColSp <- col.pts
    plot(m, scaling=scaling, type="none", ...)
    #text(m, "species", col=ColSp, cex=0.5, scaling=scaling)
    points(m, "species", col=ColSp, cex=cex.pts, pch=19, scaling=scaling)
    text(m, "sites", col=col.text, cex=cex.text, scaling=scaling)
    invisible(m)
}

par(mfrow=c(2,2), las=1)
ord_fun("birds", TRUE, main="Birds, North")
ord_fun("birds", FALSE, main="Birds, South")
ord_fun("mites", TRUE, main="Mites, North")
ord_fun("mites", FALSE, main="Mites, South")

par(mfrow=c(2,2), las=1)
ord_fun("lichens", TRUE, main="Lichens, North")
ord_fun("lichens", FALSE, main="Lichens, South")
ord_fun("mosses", TRUE, main="Bryophytes, North")
ord_fun("mosses", FALSE, main="Bryophytes, South")

par(mfrow=c(1,2), las=1)
ord_fun("vplants", TRUE, main="Vascular Plants, North")
ord_fun("vplants", FALSE, main="Vascular Plants, South")

par(mfrow=c(1,2), las=1)
mn <- ord_fun("All", TRUE, main="All Species, North", alpha=0.5, scaling=0, col.text=1)
ms <- ord_fun("All", FALSE, main="All Species, South", alpha=0.5, scaling=0, col.text=1)

library(spatstat)
dfun <- function(pp) {
    ppp <- density(pp)
    list(x=ppp$xcol, y=ppp$yrow, z=t(ppp$v)/pp$n)
}

mn <- ord_fun("All", TRUE, main="All Species, North", alpha=0.5, scaling=2, col.text=1)
s <- scores(mn, 1:3)
tax <- do.call(c, lapply(names(AllIn), function(zz) rep(zz, nrow(AllIn[[zz]][["veg"]]))))
names(tax) <- do.call(c, lapply(AllIn, function(zz) as.character(zz[["veg"]][,1])))
tax <- tax[rownames(s$species)]
ppn <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    as.ppp(s$species[j,1:2], c(-2,4,-2,2)) #c(range(s$species[,1]), range(s$species[,2])))
})
mn <- ord_fun("All", TRUE, main="All Species, North", alpha=0.5, scaling=2, col.text=1)
for (i in 1:length(ppn))
    contour(dfun(ppn[[i]]), levels=0.05, labels=names(ppn)[i], col=i, add=TRUE)
legend("topright", lty=1, col=1:length(ppn), legend=names(ppn), bty="n")

ms <- ord_fun("All", FALSE, main="All Species, South", alpha=0.5, scaling=2, col.text=1)
s <- scores(ms, 1:3)
tax <- do.call(c, lapply(names(AllIn), function(zz) rep(zz, nrow(AllIn[[zz]][["soilnt"]]))))
names(tax) <- do.call(c, lapply(AllIn, function(zz) as.character(zz[["soilnt"]][,1])))
tax <- tax[rownames(s$species)]
pps <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    as.ppp(s$species[j,1:2], c(-2,3,-3,2)) #c(range(s$species[,1]), range(s$species[,2])))
})

ms <- ord_fun("All", FALSE, main="All Species, South", alpha=0.5, scaling=2, col.text=1)
for (i in 1:length(pps))
    contour(dfun(pps[[i]]), levels=0.05, labels=names(pps)[i], col=i, add=TRUE)
legend("topright", lty=1, col=1:length(pps), legend=names(pps), bty="n")

## chull or pp density

pdf("~/Dropbox/abmi/10yr/ch4/figs/ord-ell.pdf", height=7, width=14)
par(mfrow=c(1,2), las=1)
mn <- ord_fun("All", TRUE, main="All Species, North", alpha=0.5, scaling=2, col.text=1,
    col.pts=NA, cex.pts=0.5, cex.text=1)
s <- scores(mn, 1:3)
tax <- do.call(c, lapply(names(AllIn), function(zz) rep(zz, nrow(AllIn[[zz]][["veg"]]))))
names(tax) <- do.call(c, lapply(AllIn, function(zz) as.character(zz[["veg"]][,1])))
tax <- tax[rownames(s$species)]
xy <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    s$species[j,1:2]
})
chn <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    xx <- s$species[j,1:2]
    xx[chull(xx),]
})
ppn <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    as.ppp(s$species[j,1:2], c(-2,4,-2,2)) #c(range(s$species[,1]), range(s$species[,2])))
})
ell <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    X <- s$species[j,1:2]
    W <- weights(mn, display="species")[j]
    mat <- cov.wt(X, W)
    t <- sqrt(qchisq(0.95, 2))
    xy <- vegan:::veganCovEllipse(mat$cov, mat$center, t)
})
for (i in 1:length(xy))
    points(xy[[i]], col=paste0(colTd[i], "80"), pch=19, cex=0.5)
#for (i in 1:length(chn))
#    polygon(chn[[i]], border=colTd[i], col=paste0(colTd[i], "10"))
#for (i in 1:length(ppn))
#    contour(dfun(ppn[[i]]), levels=0.05, labels=names(ppn)[i], col=colTd[i], add=TRUE)
for (i in 1:length(ell))
    lines(ell[[i]], col=colTd[i])
#legend("topleft", lty=1, lwd=2, col=colTd, legend=names(chn), bty="n")

ms <- ord_fun("All", FALSE, main="All Species, South", alpha=0.5, scaling=2, col.text=1,
    col.pts=NA, cex.pts=0.5, cex.text=1)
s <- scores(ms, 1:3)
tax <- do.call(c, lapply(names(AllIn), function(zz) rep(zz, nrow(AllIn[[zz]][["soilnt"]]))))
names(tax) <- do.call(c, lapply(AllIn, function(zz) as.character(zz[["soilnt"]][,1])))
tax <- tax[rownames(s$species)]
xy <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    s$species[j,1:2]
})
chs <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    xx <- s$species[j,1:2]
    xx[chull(xx),]
})
ppn <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    as.ppp(s$species[j,1:2], c(-2,4,-2,2)) #c(range(s$species[,1]), range(s$species[,2])))
})
ell <- lapply(structure(unique(tax), names=unique(tax)), function(i) {
    j <- tax == i
    X <- s$species[j,1:2]
    W <- weights(ms, display="species")[j]
    mat <- cov.wt(X, W)
    t <- sqrt(qchisq(0.95, 2))
    xy <- vegan:::veganCovEllipse(mat$cov, mat$center, t)
})
for (i in 1:length(xy))
    points(xy[[i]], col=paste0(colTd[i], "80"), pch=19, cex=0.5)
#for (i in 1:length(chs))
#    polygon(chs[[i]], border=colTd[i], col=paste0(colTd[i], "10"))
#for (i in 1:length(ppn))
#    contour(dfun(ppn[[i]]), levels=0.05, labels=names(ppn)[i], col=colTd[i], add=TRUE)
for (i in 1:length(ell))
    lines(ell[[i]], col=colTd[i])
legend("bottomright", lty=1, lwd=2, col=colTd, legend=names(chn), bty="n")
dev.off()

## --- SI maps

## plot SI

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(gstat)
library(raster)
#library(viridis)
load("e:/peter/sppweb2017/all-intactness.Rdata")
load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))

xy <- kgrid
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

rt <- raster(file.path("e:/peter/AB_data_v2016", "data", "kgrid", "AHM1k.asc"))
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
projection(rt) <- crs
xy <- spTransform(xy, crs)
mat0 <- as.matrix(rt)
matR <- as.matrix(Xtab(ifelse(kgrid$NRNAME=="Rocky Mountain", 1, 0) ~ Row + Col, kgrid))


setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(rt))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province

col <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B",
    "#FFFFBF","#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)

pdf("e:/peter/sppweb2017/SI-maps.pdf", onefile=TRUE, height=8, width=5)
for (i in colnames(MatSI)) {
    cat(i, "\n");flush.console()
    SI <- MatSI[,i]
    SI[is.na(SI)] <- 100 # water/ice
    mat <- as.matrix(Xtab(SI ~ Row + Col, kgrid))
    mat[is.na(mat0)] <- NA
    rout <- raster(x=mat, template=rt)
    op <- par(mar=rep(1, 4))
    plot(rout, col=col, axes=FALSE, box=FALSE, main=i)
    plot(ABnr, add=TRUE)
    par(op)
}
dev.off()

pdf("~/Dropbox/abmi/10yr/ch4/figs/si-map.pdf", height=16, width=15)
op <- par(mar=c(2,1,2,6), mfrow=c(2,3))
for (i in colnames(MatSI)) {
    cat(i, "\n");flush.console()
    SI <- MatSI[,i]
    SI[is.na(SI)] <- 100 # water/ice
    mat <- as.matrix(Xtab(SI ~ Row + Col, kgrid))
    mat[is.na(mat0)] <- NA
#    mat[matR==1] <- NA
    rout <- raster(x=mat, template=rt)
    plot(rout, col=col, axes=FALSE, box=FALSE,
        main=if(i=="vplants") paste(i, "(native)") else i, legend=i=="mites")
    plot(ABnr, add=TRUE, col=c(NA,NA,NA,NA,NA,"grey"))
}
par(op)
dev.off()

## --- SI distributions: violplots, Prov/N/S  by taxon

get_stuff0si <- function(TAX, REG) {
    SPP <- if (TAX=="birds")
        AllIn[[TAX]]$lt$AOU else AllIn[[TAX]]$lt$SpeciesID
    i <- AllIn[[TAX]]$lt[,"map.pred"]
    if (TAX=="vplants") {
        nnn <- !is.na(AllIn[[TAX]]$lt$origin) & AllIn[[TAX]]$lt$origin == "Exotic"
        i[nnn] <- FALSE
    }
    SPP <- as.character(SPP[i])
    z <- AllSI[[TAX]][SPP]
    SI <- sapply(z, function(zz) zz[REG,"SI"])
    if (TAX=="vplants") {
        NN <- 100*(1-AllIn[[TAX]]$lt[nnn,"pOcc"])
        names(NN) <- AllIn[[TAX]]$lt$SpeciesID[nnn]
#        SI <- c(SI, NN)
    }
    SI
}
get_stuffsi <- function(REG)
    rugged_mat(lapply(structure(names(AllSI), names=names(AllSI)),
        get_stuff0si, REG=REG))
get_stuffsi2 <- function(TAX)
    sapply(structure(rownames(AllSI$birds[[1]]), names=rownames(AllSI$birds[[1]])),
        function(z) get_stuff0si(REG=z, TAX=TAX))

## knocking out regions where regional reference is <2% of AB total reference
#for (i in names(AllSI)) {
#    for (j in names(AllSI[[i]])) {
#        tt <- AllSI[[i]][[j]]
#        Less <- tt[,"Ref"] < 0.02* tt["Alberta","Ref"]
#        AllSI[[i]][[j]][Less, "SI"] <- NA
#    }
#}

matSI <- lapply(structure(rownames(AllSI$birds[[1]]), names=rownames(AllSI$birds[[1]])),
    get_stuffsi)

SItoSave <- lapply(structure(names(AllSI), names=names(AllSI)), get_stuffsi2)
for (i in names(SItoSave))
    write.csv(SItoSave[[i]], file=paste0("e:/peter/sppweb2017/SI-", i, ".csv"))

ylim <- c(0, 100)
p <- c(0, 1)
pdf("~/Dropbox/abmi/10yr/ch4/figs/si.pdf", height=10, width=10)
par(mfrow=c(3,2), las=1, yaxs="i")
for (i in c(1,3,4,6,7))
    vp(matSI[[i]], main=names(matSI)[i], ylim=ylim, p=p, ylab="Species Intactness",
        col=colTl, border=colTd)
dev.off()

par(mfrow=c(1,1), las=1, yaxs="i")
vp(matSI[["Alberta"]], main="Alberta", ylim=ylim, p=p, ylab="Species Intactness")

## --- SI vs HF: gam/loess splines Prov/N/S (thf, alien, succ by taxon)

## this is 2012 HF -- might want to update to 2014?
load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))
load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata"))

kgrid$THF <- rowSums(dd1km_pred$veg_current[,setdiff(colnames(dd1km_pred$veg_current),
    colnames(dd1km_pred$veg_reference))]) / rowSums(dd1km_pred$veg_current)
CClabs <- colnames(dd1km_pred$veg_current)[grep("CC", colnames(dd1km_pred$veg_current))]
kgrid$Succ <- rowSums(dd1km_pred$veg_current[,c("SeismicLine","TransmissionLine","Pipeline",
    "RailVegetatedVerge","RoadTrailVegetated","RoadVegetatedVerge",
    CClabs)]) / rowSums(dd1km_pred$veg_current)
kgrid$Alien <- kgrid$THF - kgrid$Succ

si_cf <- list()
pdf("~/Dropbox/abmi/10yr/ch4/figs/si-tri.pdf", height=8, width=12)
par(mfrow=c(2,3), las=1)
for (i in colnames(MatSI)) {

SItmp <- MatSI[,i]
SItmp[SItmp<1] <- 1
kgrid$iSI <- log(SItmp/100)

#EXCL <- if (i %in% c("lichens", "mosses"))
#     c("Rocky Mountain", "Grassland") else "Rocky Mountain"
EXCL <- "Rocky Mountain"
kgrid2 <- kgrid[!(kgrid$NRNAME %in% EXCL),]
pol <- poly(kgrid2$Succ, kgrid2$Alien, degree=3, raw=TRUE)
m <- lm(iSI ~ pol-1, kgrid2)
#summary(m)
si_cf[[i]] <- coef(m)

BY <- 0.005
Succ <- seq(0, 1, by=BY)
Alien <- seq(0, 1, by=BY)
nd <- expand.grid(Succ=Succ, Alien=Alien)
nd$sum <- nd$Succ+nd$Alien
npol <- poly(nd$Succ, nd$Alien, degree=3, raw=TRUE)
pr <- drop(exp(npol %*% coef(m)))

tmp <- matrix(nd$sum, length(Succ), length(Alien))
tmp <- ifelse(tmp>1, NA, 1)
l <- list(x=Succ*100, y=Alien*100, z=100*tmp*matrix(pr, length(Succ), length(Alien)))
rr <- round(range(l$z, na.rm=TRUE))
rr[1] <- max(rr[1], 1)
rr[2] <- min(rr[2], 100)
image(l, col=col[rr[1]:rr[2]], main=i, xlab="Successional %", ylab="Alienating %", axes=FALSE)
contour(l, add=TRUE)
axis(1);axis(2)
abline(100,-1)
}
dev.off()

dd <- round(100*kgrid[,c("Succ", "Alien")])+1
dd$Succ <- factor(dd$Succ, 1:101)
dd$Alien <- factor(dd$Alien, 1:101)
xt <- as.matrix(Xtab(~Succ+Alien,dd))
BY <- 0.01
Succ <- seq(0, 1, by=BY)
Alien <- seq(0, 1, by=BY)
nd <- expand.grid(Succ=Succ, Alien=Alien)
nd$sum <- nd$Succ+nd$Alien
tmp <- matrix(nd$sum, length(Succ), length(Alien))
tmp <- ifelse(tmp>1, NA, 1)
l <- list(x=Succ*100, y=Alien*100, z=100*tmp*log(xt+1))
image(l, main="Representation",
    xlab="Successional %", ylab="Alienating %", axes=FALSE,
    col=grey(100:0/100))
axis(1);axis(2)

## 1D plots
do.call(cbind, si_cf)

BY <- 0.005
v <- seq(0, 1, by=BY)

Succ <- v
Alien <- 0
nd <- expand.grid(Succ=Succ, Alien=Alien)
nd$sum <- nd$Succ+nd$Alien
npol <- poly(nd$Succ, nd$Alien, degree=3, raw=TRUE)
lpr1 <- sapply(si_cf, function(z) drop(exp(npol %*% z)))

Succ <- 0
Alien <- v
nd <- expand.grid(Succ=Succ, Alien=Alien)
nd$sum <- nd$Succ+nd$Alien
npol <- poly(nd$Succ, nd$Alien, degree=3, raw=TRUE)
lpr2 <- sapply(si_cf, function(z) drop(exp(npol %*% z)))

pdf("~/Dropbox/abmi/10yr/ch4/figs/si-bi.pdf", height=7, width=12)
par(mfrow=c(1,2), las=1, yaxs="i", xaxs="i")
col1 <- colTd
matplot(v*100, lpr1*100, lty=1, type="l", ylim=c(0,100), lwd=2, col=col1,
    ylab="Intactness", xlab="Successional Footprint % (Alienating=0)")
abline(h=c(20, 40, 60, 80), v=c(20, 40, 60, 80), col="lightgrey")
matlines(v*100, lpr1*100, lty=1, lwd=3, col=col1)
legend("bottomleft", lty=1, lwd=3, col=col1, legend=colnames(lpr1), bty="n")
matplot(v*100, lpr2*100, lty=1, type="l", ylim=c(0,100), lwd=2, col=col1,
    ylab="Intactness", xlab="Alienating Footprint % (Successional=0)")
abline(h=c(20, 40, 60, 80), v=c(20, 40, 60, 80), col="lightgrey")
matlines(v*100, lpr2*100, lty=1, lwd=3, col=col1)
dev.off()

## carrot graphs for AUC and other validation metrics
load("~/Dropbox/abmi/10yr/R/AllBoot.Rdata")

cn <- c("birds", "vplants", "lichens", "mosses", "mites")
r2n <- rugged_mat(lapply(gofn, function(z) z[,2]))[,cn]
r2s <- rugged_mat(lapply(gofs, function(z) z[,2]))[,cn]
aucn <- rugged_mat(lapply(gofn, function(z) z[,4]))[,cn]
aucs <- rugged_mat(lapply(gofs, function(z) z[,4]))[,cn]
dwn <- rugged_mat(lapply(gofn, function(z) z[,8]))[,cn]
dws <- rugged_mat(lapply(gofs, function(z) z[,8]))[,cn]

pdf("~/Dropbox/abmi/10yr/ch4/figs/valid.pdf", height=5, width=5, onefile=TRUE)
par(las=1, yaxs="i")
ylim <- c(0, 1)
vp(r2n, main="R^2, North", ylim=ylim, p=p, ylab="R^2", col=colTl[-1], border=colTd[-1])
vp(r2s, main="R^2, South", ylim=ylim, p=p, ylab="R^2", col=colTl[-1], border=colTd[-1])
vp(aucn, main="AUC, North", ylim=ylim, p=p, ylab="AUC", col=colTl[-1], border=colTd[-1])
vp(aucs, main="AUC, South", ylim=ylim, p=p, ylab="AUC", col=colTl[-1], border=colTd[-1])
vp(dwn, main="Consistency, North", ylim=ylim, p=p, ylab="Consistency", col=colTl[-1], border=colTd[-1])
vp(dws, main="Consistency, South", ylim=ylim, p=p, ylab="Consistency", col=colTl[-1], border=colTd[-1])
dev.off()


## geoxv results

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(gstat)
library(raster)
#library(viridis)
load("e:/peter/sppweb2017/all-intactness.Rdata")
load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata"))

xy <- kgrid
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

rt <- raster(file.path("e:/peter/AB_data_v2016", "data", "kgrid", "AHM1k.asc"))
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
projection(rt) <- crs
xy <- spTransform(xy, crs)
mat0 <- as.matrix(rt)

setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(rt))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province

load("e:/peter/sppweb2017/birds-auc-geoxv.Rdata")

cln <- read.csv("e:/peter/AB_data_v2017/data/analysis/validation-clusters/Northern-site-clusters-for-geographical-model-validation.csv")
cls <- read.csv("e:/peter/AB_data_v2017/data/analysis/validation-clusters/Southern-site-clusters-for-geographical-model-validation.csv")

sxy <- read.csv("c:/Users/Peter/repos/abmianalytics/lookup/sitemetadata.csv")
sxy$Ngeo <- cln$Geo_Group[match(sxy$SITE_ID, cln$Site)]
sxy$Sgeo <- cls$Geo_Group[match(sxy$SITE_ID, cls$Site)]

coordinates(sxy) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(sxy) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sxy <- spTransform(sxy, crs)

COL <- c('#e6f5c9','#f4cae4','#b3e2cd','#fff2ae','#fdcdac','#cbd5e8')
library(RColorBrewer)
geocol <- brewer.pal(8, "Dark2")

geolabelN <- c("1 North-East","2 South-West","3 South","4 South-East","5 North-West",
    "6 Central-East","7 Central-West","8 Centre")
geolabelS <- c("1 North","2 East","3 Centre","4 West", "5 South-East","6 South-West",
    "7 North-West")

cn <- groupMeans(coordinates(sxy)[!is.na(sxy@data$Ngeo),], 1,
    as.character(sxy@data$Ngeo)[!is.na(sxy@data$Ngeo)])
cn <- cn[order(rownames(cn)),]
rownames(cn) <- geolabelN
cs <- groupMeans(coordinates(sxy)[!is.na(sxy@data$Sgeo),], 1,
    as.character(sxy@data$Sgeo)[!is.na(sxy@data$Sgeo)])
cs <- cs[order(rownames(cs)),]
rownames(cs) <- geolabelS

pdf("e:/peter/sppweb2017/geoxv-maps.pdf", onefile=TRUE, height=8, width=10)
op <- par(mar=rep(1,4), mfrow=c(1,2))
plot(ABnr, col=paste0(COL, "40"), border=paste0(COL, "80"), main="Geographic clusters, North")
points(sxy, col=paste0(geocol, "80")[sxy@data$Ngeo], pch=19)
plot(ABpr, border="grey", add=TRUE)
text(cn[,1], cn[,2], rownames(cn))

plot(ABnr, col=paste0(COL, "40"), border=paste0(COL, "80"), main="Geographic clusters, South")
points(sxy, col=paste0(geocol, "80")[sxy@data$Sgeo], pch=19)
plot(ABpr, border="grey", add=TRUE)
text(cs[,1], cs[,2], rownames(cs))
par(op)
dev.off()

## load plants results
pn <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables and Figures for 10 Year Review/VPlant/North/GeoValidation/VPlants_North_AUC All validation and modelling regions.csv")
ps <- read.csv("w:/ABMI 10 Year Review Tables to Peter/Result tables and Figures for 10 Year Review/VPlant/South/GeoValidation/VPlants_South_AUC All validation and modelling regions.csv")

pxn <- data.frame(AUC=c(pn$ROC.mod, pn$ROC.val),
    Cross=factor(c(paste0("In.", pn$Reg), paste0("Out.", pn$Reg)), levels(gxn$Cross)))
pxs <- data.frame(AUC=c(ps$ROC.mod, ps$ROC.val),
    Cross=factor(c(paste0("In.", ps$Reg), paste0("Out.", ps$Reg)), levels(gxs$Cross)))


par(mfrow=c(2,2),las=1)

boxplot(AUC ~ Cross, gxn, range=0, ylab="AUC", ylim=c(0,1),
    main="Birds, Geographical Cross-validation, North")
matlines(1:nlevels(gxn$Cross), t(mln), col=rgb(0,0,0,0.25), lty=1)
boxplot(AUC ~ Cross, gxn, range=0, add=TRUE, col=rgb(1,0.5,0.5,0.4))
abline(h=0.5, lty=2, col="grey")

boxplot(AUC ~ Cross, pxn, range=0, ylab="AUC", ylim=c(0,1),
    main="Vplants, Geographical Cross-validation, North")
matlines(1:nlevels(gxn$Cross), t(mln), col=rgb(0,0,0,0.25), lty=1)
boxplot(AUC ~ Cross, pxn, range=0, add=TRUE, col=rgb(1,0.5,0.5,0.4))
abline(h=0.5, lty=2, col="grey")

boxplot(AUC ~ Cross, gxs, range=0, ylab="AUC", ylim=c(0,1),
    main="Birds, Geographical Cross-validation, South")
matlines(1:nlevels(gxs$Cross), t(mls), col=rgb(0,0,0,0.25), lty=1)
boxplot(AUC ~ Cross, gxs, range=0, add=TRUE, col=rgb(1,0.5,0.5,0.4))
abline(h=0.5, lty=2, col="grey")

boxplot(AUC ~ Cross, pxs, range=0, ylab="AUC", ylim=c(0,1),
    main="Vplants, Geographical Cross-validation, South")
matlines(1:nlevels(gxs$Cross), t(mls), col=rgb(0,0,0,0.25), lty=1)
boxplot(AUC ~ Cross, pxs, range=0, add=TRUE, col=rgb(1,0.5,0.5,0.4))
abline(h=0.5, lty=2, col="grey")

opn<- par(mai=c(1.4,1 ,0.6,0.3), mfrow=c(2,2), las=1)
YLIM <- c(0.4, 1)
nn <- c(1,1.8, 3,3.8, 5,5.8, 7, 7.8, 9,9.8, 11, 11.8, 13,13.8, 15,15.8)
boxplot(AUC ~ Cross, gxn, xaxt='n', at=nn, ylab="AUC",
    ylim = YLIM, col = paste0(rep(geocol[1:8], each=2), rep(c("60","BB"), 8)),
    notch = FALSE,   outline = FALSE )
axis(1, at=seq(1.5,15.5,2), labels= FALSE)
text(x= seq(1.5,15.5,2),par("usr")[3],labels =geolabelN, srt =60, adj = c(1.2,1.3),
    xpd = TRUE, cex=1)
title('Birds, North', cex.main=1)
legend("bottomright", legend=c("Modelling data", "Validation data"), col="black",pt.bg=c("#1B9E7760","#1B9E77BB")  ,pch= 22 , cex=1,pt.cex=1.2 ,bty="n")

nn <- c(1,1.8, 3,3.8, 5,5.8, 7, 7.8, 9,9.8, 11, 11.8, 13,13.8, 15,15.8)
boxplot(AUC ~ Cross, pxn, xaxt='n', at=nn, ylab="AUC",
    ylim = YLIM, col = paste0(rep(geocol[1:8], each=2), rep(c("60","BB"), 8)),
    notch = FALSE,   outline = FALSE )
axis(1, at=seq(1.5,15.5,2), labels= FALSE)
text(x= seq(1.5,15.5,2),par("usr")[3],labels =geolabelN, srt =60, adj = c(1.2,1.3),
    xpd = TRUE, cex=1)
title('Vascular Plants, North', cex.main=1)
legend("bottomright", legend=c("Modelling data", "Validation data"), col="black",pt.bg=c("#1B9E7760","#1B9E77BB")  ,pch= 22 , cex=1,pt.cex=1.2 ,bty="n")

nn <- c(1,1.8, 3,3.8, 5,5.8, 7, 7.8, 9,9.8, 11, 11.8, 13,13.8 )
boxplot(AUC ~ Cross, gxs, xaxt='n',at=nn, ylab="AUC",
    ylim = YLIM, col = paste0(rep(geocol[1:7], each=2), rep(c("60","BB"), 7)),
    notch = FALSE, outline = FALSE)
axis(1, at=seq(1.5,13.5,2), labels= FALSE)
text(x= seq(1.5,13.5,2),par("usr")[3],labels =geolabelS, srt =60, adj = c(1.2,1.3),
    xpd = TRUE, cex=1)
title('Birds, South', cex.main=1)
legend("bottomright", legend=c("Modelling data", "Validation data"), col="black",pt.bg=c("#1B9E7760","#1B9E77BB")  ,pch= 22 , cex=1,pt.cex=1.2 ,bty="n")

nn <- c(1,1.8, 3,3.8, 5,5.8, 7, 7.8, 9,9.8, 11, 11.8, 13,13.8 )
boxplot(AUC ~ Cross, pxs, xaxt='n',at=nn, ylab="AUC",
    ylim = YLIM, col = paste0(rep(geocol[1:7], each=2), rep(c("60","BB"), 7)),
    notch = FALSE, outline = FALSE)
axis(1, at=seq(1.5,13.5,2), labels= FALSE)
text(x= seq(1.5,13.5,2),par("usr")[3],labels =geolabelS, srt =60, adj = c(1.2,1.3),
    xpd = TRUE, cex=1)
title('Vascular Plants, South', cex.main=1)
legend("bottomright", legend=c("Modelling data", "Validation data"), col="black",pt.bg=c("#1B9E7760","#1B9E77BB")  ,pch= 22 , cex=1,pt.cex=1.2 ,bty="n")

par(opn)


## saving Excel files

library(xlsx)
load("~/Dropbox/abmi/10yr/R/AllTables.Rdata")
TAXA <- c("mammals", "birds", "mites", "mosses", "lichens", "vplants")
FTAXA <- c(
    mammals="Mammals (Winter Active)",
    birds="Birds",
    mites="Soil Mites",
    mosses="Bryophytes",
    lichens="Lichens",
    vplants="Vascular Plants")
PROT <- c(
    mammals="Terrestrial (snow tracking transects)",
    birds="Terrestrial (point counts)",
    mites="Terrestrial (soil cores)",
    mosses="Terrestrial (1 ha plots)",
    lichens="Terrestrial (1 ha plots)",
    vplants="Terrestrial (1 ha plots)")
TABS <- names(AllIn[[1]])
SHEETS <- c(lt="TaxonomicInfo",
    usen="UseavailNorth",
    uses="UseavailSouth",
    veg="VegetationNorth",
    linn="LinearNorth",
    soilnt="SoilNontreedSouth",
    soiltr="SoilTreedSouth",
    lins="LinearSouth",
    sectn="SectorNorth",
    sects="SectorSouth")
SHEETS <- SHEETS[TABS]
MET <- read.csv("e:/peter/sppweb2017/Metadata-template.csv")
VER <- "5.0"

for (i in TAXA) {
    fn <- paste0("ABMI-species-v", VER, "_", i, "_", Sys.Date(), ".xlsx")
    fn <- file.path("e:/peter/sppweb2017/tables", fn)
    Descr <- data.frame(Description=c(Source=paste0("Alberta Biodiversity Monitoring Institute (2016). ABMI Species Website, version 5.0 (", Sys.Date(), "). URL: http://species.abmi.ca"),
        Taxon=unname(FTAXA[i]),
        Protocol=unname(PROT[i]),
        Contact="solymos@ualberta.ca",
        License="CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/"))
    cat(i)
    write.xlsx(Descr, fn, sheetName="Description", row.names=TRUE)
    for (j in TABS) {
        cat("\t", j);flush.console()
        write.xlsx(AllIn[[i]][[j]], fn, sheetName=unname(SHEETS[j]), append=TRUE, row.names=FALSE)
    }
    write.xlsx(MET, fn, sheetName="Metadata", append=TRUE, row.names=FALSE)
    cat("\n")
}

