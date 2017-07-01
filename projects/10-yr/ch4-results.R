library(mefa4)
ROOT <- "v:/contents/2017/species"
TAXA <- c("mites", "mosses", "lichens", "vplants")
AllIn <- list()
vegcols <- read.csv("c:/Users/Peter/repos/abmianalytics/projects/10-yr/veg-cols.csv")

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
TAXA <- c("birds", "mites", "mosses", "lichens", "vplants")

AllSI <- list()
MatSI <- matrix(NA, nrow(kgrid), length(TAXA))
dimnames(MatSI) <- list(rownames(kgrid), TAXA)

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

#TAX <- "mites"
for (TAX in TAXA[-1]) {
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
load("~/Dropbox/abmi/10yr/R/AllTables.Rdata")

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
vp1 <- function(x, p=c(0, 1), ...) {
    x2 <- x[!is.na(x)]
    if (length(unique(x2)) < 5)
        return(NULL)
    q <- quantile(x2, p)
    x2 <- x2[x2 > q[1] & x2 < q[2]]
    d <- density(x2, ...)
    data.frame(h=d$x, w=d$y/max(d$y))
}
vp <- function(x, p=c(0,1),
col="#FF664D80", border="#FF664D",
#col="#40E0D180", border="#40E0D1",
ylim, ylab="", xlab="", main="", ...) {
    xx <- apply(x, 2, vp1, p=p, ...)
    if (missing(ylim))
        ylim <- range(x, na.rm=TRUE)
    plot(0, type="n", axes=FALSE, ann=FALSE, xlim=c(0.5, ncol(x)+0.5),
        ylim=ylim)
    axis(1, 1:ncol(x), colnames(x), lwd=0)
    axis(2)
    title(ylab=ylab, xlab=xlab, main=main)
    for (i in 1:ncol(x)) {
        if (!is.null(xx[[i]])) {
            polygon(0.45*c(-xx[[i]]$w, rev(xx[[i]]$w))+i,
                c(xx[[i]]$h, rev(xx[[i]]$h)), col=col, border=border)
            #lines(c(i,i), range(xx[[i]]$h), col=col, lwd=2)
            lines(c(i-0.2, i+0.2), rep(median(x[,i], na.rm=TRUE), 2), lwd=3)
        }
    }
    #points(1:ncol(x), colMeans(x, na.rm=TRUE), pch=21, cex=1.5)
    invisible(NULL)
}
fu <- function(x) sign(x) * plogis(log(abs(x/100)))

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
par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(0, 3)
p <- c(0.025, 0.975)
vp(softn, main="Soft Linear, North", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)
vp(hardn, main="Hard Linear, North", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)
vp(softs, main="Soft Linear, South", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)
vp(hards, main="Hard Linear, South", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)

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

par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(-1, 1)
p <- c(0, 1)
vp(fu(uneng), main="Energy, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unagr), main="Algiculture, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unurb), main="Rural/Urban, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(untra), main="Transportation, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unfor), main="Forestry, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
plot.new()

par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(-1, 1)
p <- c(0, 1)
vp(fu(useng), main="Energy, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(usagr), main="Algiculture, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(usurb), main="Rural/Urban, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(ustra), main="Transportation, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")


par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(-10, 10)
p <- c(0.025, 0.975)
vp(tneng, main="Energy, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tnagr, main="Algiculture, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tnurb, main="Rural/Urban, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tntra, main="Transportation, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tnfor, main="Forestry, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
plot.new()

par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(-20, 20)
p <- c(0.025, 0.975)
vp(tseng, main="Energy, South", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tsagr, main="Algiculture, South", ylim=c(-100, 100), p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tsurb, main="Rural/Urban, South", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tstra, main="Transportation, South", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)

## --- ordination

ord_fun <- function(TAX, north=TRUE, scaling=2, alpha=1, col.text=1, ...) {
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
    rownames(zz) <- make.cepnames(rownames(zz))
    colnames(zz) <- make.cepnames(colnames(zz))

    library(vegan)
    fs <- function(x) 0.5*(1+x/max(abs(x)))
    yy <- (t(zz))
    m <- cca(yy)
    s <- scores(m, 1:3)
    ColSi <- 2
    ColSp <- rgb(red=fs(s$species[,1]),
        green=fs(s$species[,2]), blue=fs(s$species[,3]), alpha=alpha)
    plot(m, scaling=scaling, type="none", ...)
    #text(m, "species", col=ColSp, cex=0.5, scaling=scaling)
    points(m, "species", col=ColSp, cex=1, pch=19, scaling=scaling)
    text(m, "sites", col=col.text, cex=1, scaling=scaling)
    invisible(NULL)
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
ord_fun("All", TRUE, main="All Species, North", alpha=0.5, scaling=0, col.text=1)
ord_fun("All", FALSE, main="All Species, South", alpha=0.5, scaling=0, col.text=1)

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

setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, proj4string(rt))
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province

col <- colorRampPalette(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B",
    "#FFFFBF","#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837"))(100)

## lichens & mosses are not OK
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

## --- SI distributions: violplots, Prov/N/S  by taxon

get_stuff0si <- function(TAX, REG) {
    SPP <- if (TAX=="birds")
        AllIn[[TAX]]$lt$AOU else AllIn[[TAX]]$lt$SpeciesID
    i <- AllIn[[TAX]]$lt[,"map.pred"]
    if (TAX=="vplants")
        i[!is.na(AllIn[[TAX]]$lt$origin) & AllIn[[TAX]]$lt$origin == "Exotic"] <- FALSE
    SPP <- as.character(SPP[i])
    z <- AllSI[[TAX]][SPP]
    sapply(z, function(zz) zz[REG,"SI"])
}
get_stuffsi <- function(REG)
    rugged_mat(lapply(structure(names(AllSI), names=names(AllSI)),
        get_stuff0si, REG=REG))

matSI <- lapply(structure(rownames(AllSI$birds[[1]]), names=rownames(AllSI$birds[[1]])),
    get_stuffsi)

par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(0, 100)
p <- c(0, 1)
for (i in 1:6)
    vp(matSI[[i]], main=names(matSI)[i], ylim=ylim, p=p, ylab="Species Intactness")

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

par(mfrow=c(3,2), las=1)
for (i in colnames(MatSI)) {

SItmp <- MatSI[,i]
SItmp[SItmp<1] <- 1
kgrid$iSI <- log(SItmp/100)

EXCL <- if (i %in% c("lichens", "mosses"))
     c("Rocky Mountain", "Grassland") else "Rocky Mountain"
kgrid2 <- kgrid[!(kgrid$NRNAME %in% EXCL),]
pol <- poly(kgrid2$Succ, kgrid2$Alien, degree=3, raw=TRUE)
m <- lm(iSI ~ pol-1, kgrid2)
summary(m)

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
image(l, col=col, main=i, xlab="Successional %", ylab="Alienating %", axes=FALSE)
contour(l, add=TRUE)
axis(1);axis(2)
}

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
