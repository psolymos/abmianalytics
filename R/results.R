library(mefa4)
library(pbapply)
library(RColorBrewer)

ROOT <- "c:/p/AB_data_v2015/out/birds"

level <- 0.9

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
    source("~/repos/abmianalytics/R/results_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

e <- new.env()
load(file.path(ROOT, "data", "data-full.Rdata"), envir=e)
dat <- e$DAT
dat <- dat[dat$useOK,]
yy <- e$YY[rownames(dat),]
tax <- droplevels(e$TAX[colnames(yy),])

en <- new.env()
load(file.path(ROOT, "data", "data-useok-north.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "data", "data-useok-south.Rdata"), envir=es)
xns <- es$DAT
modss <- es$mods
yys <- es$YY
rm(e, en, es)

## model for species
fl <- list.files(file.path(ROOT, "results"))
fln <- fl[grep("-north_", fl)]
fln <- sub("birds_bam-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- fl[grep("-south_", fl)]
fls <- sub("birds_bam-south_", "", fls)
fls <- sub(".Rdata", "", fls)

tax$ndet <- colSums(yy>0)
tax$modelN <- rownames(tax) %in% fln
tax$modelS <- rownames(tax) %in% fls
tax$ndet_n <- colSums(yyn>0)[match(colnames(yy), colnames(yyn))]
tax$ndet_s <- colSums(yyn>0)[match(colnames(yy), colnames(yys))]

yy <- yy[,tax$ndet > 0]
tax <- droplevels(tax[colnames(yy),])

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

stage_hab_n <- 5
stage_hab_s <- 2


## spp specific output

spp <- "CAWA"


## map-det
## useavail-north
## useavail-south
## table: useavail-north
## table: useavail-south

## FIXME produce plots / save tables --------------------------------------------- FIXME

## veghf-north
## soilhf-treed-south
## soilhf-nontreed-south
## linear-north
## linear-south
## table: veghf-north
## table: soilhf-north

resn <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-north_", spp, ".Rdata")))
ress <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-south_", spp, ".Rdata")))
estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
ests_hab <- getEst(ress, stage=stage_hab_s, na.out=FALSE, Xns)
prn <- pred_veghf(estn_hab, Xnn)
prs <- pred_soilhf(ests_hab, Xns)

## FIXME produce plots / save tables--------------------------------------------- FIXME

## surroundinghf-north
## surroundinghf-south
## table: residual climate coefs north (with surrounding hf)
## table: residual climate coefs south (with surrounding hf)

## climate & surrounding hf
resn <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-north_", spp, ".Rdata")))
ress <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-south_", spp, ".Rdata")))
cn <- c("xPET", "xMAT", "xAHM", "xFFP", 
    "xMAP", "xMWMT", "xMCMT", "xlat", "xlong", "xlat2", "xlong2", 
    "THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", 
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", 
    "Noncult2_KM")
estn_sp <- getEst(resn, stage=stage_hab_n + 2, na.out=FALSE, Xnn)
ests_sp <- getEst(ress, stage=stage_hab_s + 2, na.out=FALSE, Xns)
sp_n <- colMeans(estn_sp[,cn])
sp_s <- colMeans(ests_sp[,cn])

## FIXME surrounding plot--------------------------------------------- FIXME
## FIXME save tables--------------------------------------------- FIXME

## trend-north
## trend-south
## table: north and south trend estimates

resn <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-north_", spp, ".Rdata")))
ress <- loadSPP(file.path(ROOT, "results", paste0("birds_bam-south_", spp, ".Rdata")))
estn_yr <- getEst(resn, stage=stage_hab_n + 3, na.out=FALSE, Xnn)
ests_yr <- getEst(ress, stage=stage_hab_s + 3, na.out=FALSE, Xns)
yr_n <- 100 * (exp(estn_yr[,"YR"]) - 1)
yr_s <- 100 * (exp(ests_yr[,"YR"]) - 1)

yr <- unlist(c(North=fstat(yr_n), 
    South=fstat(yr_s), 
    tax[spp,c("ndet", "ndet_n", "ndet_s")]))








## --
setwd("~/Dropbox/josm/dataproc3/")
RESDIR <- "~/Dropbox/josm/dataproc3/wg/results/"

use_restricted <- FALSE#TRUE
do_soils <- F#TRUE

level <- 0.9

library(mefa4)
library(pbapply)

source("~/Dropbox/josm/dataproc3/AB_04_making_sense_functions.R")
source("~/Dropbox/josm/dataproc3/AB_04_making_sense_functions2.R")

load("~/Dropbox/josm/dataproc3/wg/WGData-2015-03-11.Rdata")

ip_name <- "hab1ec"

## covariates

Terms <- getTerms(mods, "list")
#Terms <- c("MAXDISTANCE", "MAXDURATION", "JDAY", "TSSR", "TREEcombo", "LCC", 
#    "wgrid", "HabitatA3", Terms)

if (FALSE) { # Anj

table(samp(mm)$YR, samp(mm)$YR5)

load("c:/p/AB_data_v2014/josm2015/HabDataForAnjolene.Rdata")
setdiff(Terms, colnames(samp(mm)))
samp(mm)$YR <- 17
samp(mm)$YR5 <- factor(2, levels=0:2)
setdiff(Terms, colnames(samp(mm)))

XHF <- as.matrix(samp(mm)[,grep("IP_", colnames(samp(mm)))])
colnames(XHF) <- sub("IP_", "", colnames(XHF))
cn <- levels(samp(mm)$hab1ec)
XHF <- XHF[,cn]

}

## left-out portion is included here
setdiff(Terms, colnames(samp(mm)))
xn <- model.frame(getTerms(mods, "formula"), samp(mm)[,Terms], 
    na.action=na.pass)
#xn$Remn_QS <- as.numeric(0)
#xn$Remn2_QS <- as.numeric(0)
iout <- setdiff(1:nrow(xn), BB[,1])
dim(BB)
range(iout)
iout <- 1:nrow(xn) %in% iout

Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
#Xn <- cbind(Xn, "Wzone:YR52"=0)

#rn <- sample(rownames(xn), nrow(xn))
rn <- sample(rownames(xn), 5000)
xnqs <- xn[rn,]
Xnqs <- Xn[rn,]

fl <- list.files("~/Dropbox/josm/dataproc3/wg/results/")
fl <- sub("birds_josm_", "", fl)
fl <- sub(".Rdata", "", fl)

e <- new.env()
load("c:/p/AB_data_v2014/josm2015/JOSM-Data-2015-03-10.Rdata", envir=e)
samp(mm)$MAXDURATION <- samp(e$mm)$MAXDURATION
samp(mm)$MAXDISTANCE <- samp(e$mm)$MAXDISTANCE

stopifnot(all(rownames(XHF) == rownames(Xn)))

## -------------- plots etc.-----------------------------------------------------

if (FALSE) { # Anj
spp <- "CAWA"
res <- loadSPP(spp)
est <- getEst(res, stage=6)
pr4 <- getDataPred(res, stage=4, X=Xn, remn=FALSE)
pr6 <- getDataPred(res, stage=6, X=Xn, remn=TRUE)
#pr7 <- getDataPred(res, stage=7, X=Xn, remn=TRUE)

v4 <- t(apply(exp(pr4), 1, quantile, c(0.5, 0.05, 0.95)))
v6 <- t(apply(exp(pr7), 1, quantile, c(0.5, 0.05, 0.95)))
#v7 <- t(apply(exp(pr7), 1, quantile, c(0.5, 0.05, 0.95)))

out <- data.frame(pkey=rownames(mm), local=v4, local_qs=v6)
write.csv(out, file="CAWA-density-per-ha.csv", row.names=FALSE)

}

DIR <- "c:/p/AB_data_v2014/josm2015/josm-report-output/"

spp <- "OVEN"

## habitat
habres <- list()
for (spp in fl) {
NAM <- as.character(taxa(mm)[spp, "COMMON_NAME"])
cat(NAM, "\n");flush.console()
png(paste0(DIR, "hab/", spp, ".png"), width=1200, height=600)
#res <- loadSPP(spp, do_soils)
habres[[spp]] <- try(fig_hab(spp, burn=FALSE))
#fig_hab(spp, burn=FALSE)
dev.off()
}
save(habres, file="c:/p/AB_data_v2014/josm2015/habres.Rdata")

## surrounding hf
for (spp in fl) {
NAM <- as.character(taxa(mm)[spp, "COMMON_NAME"])
cat(NAM, "\n");flush.console()
png(paste0(DIR, "hf/", spp, ".png"), width=1200, height=400)
#res <- loadSPP(spp, do_soils)
try(fig_hf2(spp))
dev.off()
}


## mid figures
for (spp in fl) {
NAM <- as.character(taxa(mm)[spp, "COMMON_NAME"])
cat(NAM, "\n");flush.console()
res <- loadSPP(spp)
png(paste0(DIR, "mid/", spp, ".png"), width=800, height=800)
try(plotMid(res, mods))
dev.off()
}

## Lc figures
lcres <- list()
for (spp in fl) {
NAM <- as.character(taxa(mm)[spp, "COMMON_NAME"])
cat(NAM, "\n");flush.console()
png(paste0(DIR, "lc/", spp, ".png"), width=600, height=600)
lcres[[spp]] <- try(plotLc(spp))
dev.off()
}
save(lcres, file="c:/p/AB_data_v2014/josm2015/lcres.Rdata")

## gof figures
xtab(mm)[47530,"HOLA"] <- 1
gofres <- list()
for (spp in fl) {
NAM <- as.character(taxa(mm)[spp, "COMMON_NAME"])
cat(NAM, "\n");flush.console()
res <- loadSPP(spp)
png(paste0(DIR, "gof/", spp, ".png"), width=2400, height=1600, res=120)
gofres[[spp]] <- try(plotGof(spp))
dev.off()
}
save(gofres, file="c:/p/AB_data_v2014/josm2015/gofres.Rdata")

## maps & sector in separate file

