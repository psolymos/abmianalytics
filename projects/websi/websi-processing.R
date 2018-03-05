## processing common data (kgrid, areas, species lookups)

library(mefa4)
library(rgdal)
library(rgeos)
library(sp)
library(raster)

## HF 2012
load("e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv0$Sector2 <- factor(ifelse(is.na(tv0$Sector), "NATIVE", as.character(tv0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
VHF <- dd1km_pred[[1]]
tv0 <- tv0[colnames(VHF),]
dd1km_pred$sample_year
VHF <- VHF/10^6
summary(rowSums(VHF))
#RS <- rowSums(VHF)
#RS[RS < 1] <- 1
#VHF <- VHF / RS
#summary(rowSums(VHF))
KA_2012 <- groupSums(VHF, 2, tv0$Sector2) # in km^2

## HF 2014
load("e:/peter/AB_data_v2017/data/inter/veghf/veg-hf_km2014-grid_v6hf2014v2_coarse-fixage0.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
VHF <- dd1km_pred[[1]]
VHF <- VHF[,colnames(VHF) != "CutBlocks"]
cn <- colnames(VHF)
#cn[cn %in% c("UrbanIndustrial", "UrbanResidence")] <- "Urban"
#cn[cn %in% c("SeismicLineNarrow", "SeismicLineWide")] <- "SeismicLine"
#cn[cn %in% c("CultivationAbandoned", "CultivationRoughPasture",
#    "CultivationTamePasture", "CultivationCrop")] <- "CultivationCropPastureBareground"
#VHF <- groupSums(VHF, 2, cn)
tv0 <- tv0[colnames(VHF),]
dd1km_pred$sample_year
VHF <- VHF/10^6
summary(rowSums(VHF))
RS <- rowSums(VHF)
RS[RS < 1] <- 1
VHF <- VHF / RS
summary(rowSums(VHF))
KA_2014 <- groupSums(VHF, 2, tv0$ETA_UseInAnalysis_Sector) # in km^2

load("e:/peter/AB_data_v2017/data/analysis/kgrid_table_km.Rdata")
stopifnot(all(rownames(kgrid) == rownames(KA_2012)))
KA_2014 <- KA_2014[rownames(kgrid),]
stopifnot(all(rownames(kgrid) == rownames(KA_2014)))

KT <- kgrid[,c("Row", "Col", "Row10_Col10")]
XY <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(XY) <- ~ POINT_X + POINT_Y
proj4string(XY) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

SP <- read.csv("e:/peter/sppweb2017/tables/AllTaxaCombined-TaxonomicInfo.csv")
SP <- droplevels(SP[SP$map.pred,])
SP$infoOK <- TRUE
#SP$infoOKboot <- TRUE
for (i in 1:nrow(SP)) {
    fn1 <- file.path("w:/reports/2017/results", as.character(SP[i,"taxon"]), "sector",
        paste0(as.character(SP[i,"SpeciesID"]), ".RData"))
    #fn2 <- file.path("w:/reports/2017/results", as.character(SP[i,"taxon"]), "sector",
    #    paste0(as.character(SP[i,"SpeciesID"]), ".RData"))
    if (!file.exists(fn1))
        SP$infoOK[i] <- FALSE
    #if (!file.exists(fn2))
    #    SP$infoOKboot[i] <- FALSE
}
#table(SP$infoOK, SP$infoOKboot)
table(SP$infoOK, SP$taxon)

SP <- SP[SP$infoOK,]
SP$infoOK <- NULL
rownames(SP) <- SP$SpeciesID

SP$model_region <- ""
SP$model_region[SP$veghf.north] <- "North"
SP$model_region[SP$soilhf.south] <- "South"
SP$model_region[SP$veghf.north & SP$soilhf.south] <- "North and South"
SP$model_region <- as.factor(SP$model_region)
table(SP$taxon, SP$model_region, useNA="a")
SP <- SP[order(SP$taxon, SP$SpeciesID),]

SP$model_north <- SP$veghf.north
SP$model_south <- SP$soilhf.south
SP <- SP[,c("SpeciesID", "taxon",
    "Species", "CommonName", "ScientificName", "TSNID",
    "model_north", "model_south", "model_region")]

## species coefs

cv <- read.csv("e:/peter/sppweb2017/tables/AllTaxaCombined-VegetationNorth.csv")
cs <- read.csv("e:/peter/sppweb2017/tables/AllTaxaCombined-SoilNontreedSouth.csv")
cst <- read.csv("e:/peter/sppweb2017/tables/AllTaxaCombined-SoilTreedSouth.csv")
compare_sets(cv$Species, SP$Species[SP$veghf.north])
compare_sets(cs$Species, SP$Species[SP$soilhf.south])
lv <- read.csv("e:/peter/sppweb2017/tables/AllTaxaCombined-LinearNorth.csv")
ls <- read.csv("e:/peter/sppweb2017/tables/AllTaxaCombined-LinearSouth.csv")

cv <- cv[match(SP$Species[SP$model_north], cv$Species),]
cs <- cs[match(SP$Species[SP$model_south], cs$Species),]
cst <- cst[match(SP$Species[SP$model_south], cst$Species),]
lv <- lv[match(SP$Species[SP$model_north], lv$Species),]
ls <- ls[match(SP$Species[SP$model_south], ls$Species),]
rownames(lv) <- rownames(cv) <- rownames(SP)[SP$model_north]
rownames(ls) <- rownames(cs) <- rownames(cst) <- rownames(SP)[SP$model_south]

cv <- as.matrix(cv[,-(1:3)])
cs <- as.matrix(cs[,-(1:3)])
cst <- as.matrix(cst[,-(1:3)])
lv <- as.matrix(lv[,-(1:3)])
ls <- as.matrix(ls[,-(1:3)])

## insert updated mammal results

MP <- "e:/peter/sppweb2017/tables/mammal_latest"
msp <- rownames(SP)[SP$taxon=="mammals"]

m_s <- read.csv(file.path(MP, "soilhf.csv"))
m_ls <- read.csv(file.path(MP, "slin.csv"))
m_pa <- read.csv(file.path(MP, "pA.csv"))
compare_sets(msp, rownames(m_s))
compare_sets(msp, rownames(m_ls))
compare_sets(msp, rownames(m_pa))

m_km <- read.csv(file.path(MP, "kmdays.csv"))
m_km <- m_km[rownames(m_km) != "DomesticDog",]
m_v <- read.csv(file.path(MP, "veghf.csv"))
m_v <- m_v[rownames(m_v) != "DomesticDog",]
m_ln <- read.csv(file.path(MP, "nlin.csv"))
m_ln <- m_ln[rownames(m_ln) != "DomesticDog",]
compare_sets(msp, rownames(m_km))
compare_sets(msp, rownames(m_v))
compare_sets(msp, rownames(m_ln))

## kmdays
kmd <- structure(c(m_km$kmd*5 + m_km$kmd2*5^2), names=rownames(m_km))
m_s <- plogis(qlogis(as.matrix(m_s)) + kmd[rownames(m_s)])
m_v <- plogis(qlogis(as.matrix(m_v)) + kmd[rownames(m_v)])
m_ls <- plogis(qlogis(as.matrix(m_ls)) + kmd[rownames(m_ls)])
m_ln <- plogis(qlogis(as.matrix(m_ln)) + kmd[rownames(m_ln)])
## treed south
m_st <- plogis(qlogis(as.matrix(m_s)) + m_pa[rownames(m_s),])

c(sum(rownames(cv) %in% msp), nrow(m_v))
c(sum(rownames(cs) %in% msp), nrow(m_s))

cv[rownames(m_v),] <- m_v[,colnames(cv)]
cs[rownames(m_s),] <- m_s[,colnames(cs)]
cst[rownames(m_st),] <- m_st[,colnames(cst)]
lv[rownames(m_ln),] <- m_ln[,colnames(lv)]
ls[rownames(m_ls),] <- m_ls[,colnames(ls)]

## restructure

f <- function(x) {
    cn <- colnames(x)
    lo <- grepl("\\.LCL", cn)
    hi <- grepl("\\.UCL", cn)
    x0 <- x[,!hi & !lo]
    colnames(x0) <- cn[!hi & !lo]
    x1 <- x[,lo]
    colnames(x1) <- gsub("\\.LCL", "", cn[lo])
    x2 <- x[,hi]
    colnames(x2) <- gsub("\\.UCL", "", cn[hi])
    out <- list(coef=x0, lower=x1, higher=x2)
    out
}

pa <- rowMeans(qlogis(f(cst)[[1]]) - qlogis(f(cs)[[1]]))
ii <- SP[rownames(cs), "taxon"] == "birds"
tst <- f(cst)[[1]][ii,]
tst[tst == 0] <- 10^-6
ts <- f(cs)[[1]][ii,]
ts[ts == 0] <- 10^-6
pa2 <- rowMeans(log(tst) - log(ts))
pa[ii] <- pa2
pa[rownames(m_pa)] <- m_pa[,1] # not really necessary...
PA <- matrix(pa, ncol=1)
colnames(PA) <- "pAspen"
rownames(PA) <- names(pa)

## joint model results for birds to be used with climate
if (FALSE) {
library(mefa4)
ROOT <- "e:/peter/AB_data_v2016/out/birds"
level <- 0.9
up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
    invisible(NULL)
}
up()

en <- new.env()
load(file.path(ROOT, "data", "data-north.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "data", "data-south.Rdata"), envir=es)
xns <- es$DAT
modss <- es$mods
yys <- es$YY

rm(en, es)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

stage_hab_n <- 5
stage_hab_s <- 3

TAX <- read.csv("c:/Users/Peter/repos/abmispecies/_data/birds.csv")
rownames(TAX) <- TAX$AOU

## spp specific output
## raw=TRUE returns mu (lam=exp(mu))
res_coef <- list()
for (spp in rownames(TAX)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(TAX[spp, "sppid"])
    res_coef[[NAM]] <- list(
        veg=NULL, soil=NULL, paspen=NULL,
        vegj=NULL, soilj=NULL, paspenj=NULL)
    if (TAX[spp, "veghf.north"]) {
        resn <- loadSPP(file.path(ROOT, "results", "north",
            paste0("birds_abmi-north_", spp, ".Rdata")))
        estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
        estn_habj <- getEst(resn, stage=stage_hab_n+1, na.out=FALSE, Xnn)
        tmp <- pred_veghf(estn_hab, Xnn, burn_included=FALSE, raw=TRUE)[,1]
        tmpj <- pred_veghf(estn_habj, Xnn, burn_included=FALSE, raw=TRUE)[,1]
        nm <- names(tmp)
        nm <- gsub("Conif", "WhiteSpruce", nm)
        nm <- gsub("Decid", "Deciduous", nm)
        nm <- gsub("Mixwood", "Mixedwood", nm)
        nm <- gsub("BSpr", "BlackSpruce", nm)
        names(tmp) <- names(tmpj) <- nm
        names(tmp) <- names(tmpj) <- gsub(" ", "", names(tmp))
        res_coef[[NAM]]$veg <- tmp
        res_coef[[NAM]]$vegj <- tmpj
    }
    if (TAX[spp, "soilhf.south"]) {
        ress <- loadSPP(file.path(ROOT, "results", "south",
            paste0("birds_abmi-south_", spp, ".Rdata")))
        ests_hab <- getEst(ress, stage=stage_hab_s, na.out=FALSE, Xns)
        ests_habj <- getEst(ress, stage=stage_hab_s+1, na.out=FALSE, Xns)
        res_coef[[NAM]]$soil <- pred_soilhf(ests_hab, Xns, raw=TRUE)[,1]
        res_coef[[NAM]]$soilj <- pred_soilhf(ests_habj, Xns, raw=TRUE)[,1]
        res_coef[[NAM]]$paspen <- ests_hab[1,"pAspen"]
        res_coef[[NAM]]$paspenj <- ests_habj[1,"pAspen"]
    }
}
save(res_coef, file="w:/reports/2017/data/birds-coef-temp.RData")
#compare_sets(nm, colnames(CF$marginal$veg$coef))
#setdiff(nm, colnames(CF$marginal$veg$coef))
#setdiff(colnames(CF$marginal$veg$coef), nm)
}

load("w:/reports/2017/data/birds-coef-temp.RData")

CFall <- list(
    veg=f(cbind(cv, lv)),
    soil=f(cbind(cs, ls)),
    paspen=PA)

SPPv <- rownames(cv)[SP[rownames(cv), "taxon"] == "birds"]
SPPs <- rownames(cs)[SP[rownames(cs), "taxon"] == "birds"]

cf_veg_m <- t(sapply(SPPv, function(i) {
    z <- res_coef[[i]]$veg
    z <- z[match(colnames(CFall$veg$coef), names(z))]
    names(z) <- colnames(CFall$veg$coef)
    z
}))
cf_veg_j <- t(sapply(SPPv, function(i) {
    z <- res_coef[[i]]$vegj
    z <- z[match(colnames(CFall$veg$coef), names(z))]
    names(z) <- colnames(CFall$veg$coef)
    z
}))
cf_soil_m <- t(sapply(SPPs, function(i) {
    z <- res_coef[[i]]$soil
    z <- z[match(colnames(CFall$soil$coef), names(z))]
    names(z) <- colnames(CFall$soil$coef)
    z
}))
cf_soil_j <- t(sapply(SPPs, function(i) {
    z <- res_coef[[i]]$soilj
    z <- z[match(colnames(CFall$soil$coef), names(z))]
    names(z) <- colnames(CFall$soil$coef)
    z
}))
cf_pa_m <- t(sapply(SPPs, function(i) {
    res_coef[[i]]$paspen
}))
cf_pa_m <- array(cf_pa_m, c(length(SPPs), 1), list(SPPs, "pAspen"))
cf_pa_j <- t(sapply(SPPs, function(i) {
    res_coef[[i]]$paspenj
}))
cf_pa_j <- array(cf_pa_j, c(length(SPPs), 1), list(SPPs, "pAspen"))

CF <- list(coef=list(veg=CFall$veg$coef, soil=CFall$soil$coef, paspen=CFall$paspen),
    lower=list(veg=CFall$veg$lower, soil=CFall$soil$lower),
    higher=list(veg=CFall$veg$higher, soil=CFall$soil$higher))

CFbirds <- list(
    marginal=list(veg=cf_veg_m, soil=cf_soil_m, paspen=cf_pa_m),
    joint=list(veg=cf_veg_j, soil=cf_soil_j, paspen=cf_pa_j))

VER <- data.frame(
    taxon=  c("mammals", "birds", "mites", "mosses", "lichens", "vplants"),
    version=2017,
    yr_first=c(2001,     1997,    2007,    2003,     2003,      2003),
    yr_last= c(2013,     2015,    2016,    2015,     2016,      2016),
    method=c("snow_tracking","point_count","soil_core","centre_plot","centre_plot","centre_plot"),
    hf=     c("2014v2", "2012",   "2014v2", "2014v2", "2014v2", "2014v2"),
    veg=    c("v6",     "v5",     "v6",     "v6",     "v6",     "v6"),
    model=c("binomial_logit", "poisson_log", "binomial_logit", "binomial_logit", "binomial_logit", "binomial_logit"))

## save common data
save(KA_2012, KA_2014, KT, XY, SP, CF, CFbirds, VER,
    file="w:/reports/2017/data/kgrid_areas_by_sector.RData")
write.csv(SP, row.names=FALSE, file="w:/reports/2017/data/species-info.csv")


## no need to normalize files
#"Curr.Boot" "Ref.Boot" -- > CB/RB
#w:/reports/2017/Files from Ermias/Mites/Combine regions/Boot data/
#"SA.Curr" "SA.Ref" -- > CS/RS
#w:/reports/2017/Files from Ermias/Mites/Combine regions/Sector effects/Sector abundance summary/
#Native Misc Agriculture RuralUrban Energy Transportation Forestry
#Sectors <- c("Native", "Misc", "Agriculture", "RuralUrban", "Energy", "Transportation", "Forestry")

## process the rasters for the climate piece for poly level prediction

library(mefa4)
library(cure4insect)

ROOT <- "e:/peter/AB_data_v2016"
OUTDIR <- "e:/peter/AB_data_v2016/out/birds/web"

## set here if shf is wanted
load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata")) # kgrid
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/abmianalytics/R/results_functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

## climate
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
kgrid <- data.frame(kgrid, transform_CLIM(kgrid, "Row_Col"))
kgrid$xlong2 <- kgrid$xlong^2
kgrid$xlat2 <- kgrid$xlat^2

kgrid$Row_Col.1 <- NULL
kgrid$OBJECTID <- NULL
#kgrid$Row <- NULL
#kgrid$Col <- NULL
kgrid$AHM <- NULL
kgrid$PET <- NULL
kgrid$FFP <- NULL
kgrid$MAP <- NULL
kgrid$MAT <- NULL
kgrid$MCMT <- NULL
kgrid$MWMT <- NULL
kgrid$Row10 <- NULL
kgrid$Col10 <- NULL
kgrid$ASP <- NULL
kgrid$TRI <- NULL
kgrid$SLP <- NULL
kgrid$CTI <- NULL

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
kgrid$WetKM <- rowSums(dd1km_pred$veg_current[,tv[colnames(dd1km_pred$veg_current), "WET"]==1]) / rowSums(dd1km_pred$veg_current)
#kgrid$WaterKM <- rowSums(dd1km_pred$veg_current[,tv[colnames(dd1km_pred$veg_current), "WATER"]==1]) / rowSums(dd1km_pred$veg_current)
kgrid$WetWaterKM <- rowSums(dd1km_pred$veg_current[,tv[colnames(dd1km_pred$veg_current), "WETWATER"]==1]) / rowSums(dd1km_pred$veg_current)

## climate (North & South)
cnClim <- c("xPET", "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT",
    "xMCMT", "xlat", "xlong", "xlat2", "xlong2", "xFFP:xMAP",
    "xMAP:xPET", "xAHM:xMAT", "xlat:xlong", "WetKM", "WetWaterKM")
## model matrix for Clim
fclim <- as.formula(paste("~ - 1 +", paste(cnClim, collapse=" + ")))

en <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-north.Rdata"), envir=en)
xnn <- en$DAT[1:500,]
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-south.Rdata"), envir=es)
xns <- es$DAT[1:500,]
modss <- es$mods
yys <- es$YY
rm(en, es)

## model for species
fln <- list.files(file.path(ROOT, "out", "birds", "results", "north"))
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- list.files(file.path(ROOT, "out", "birds", "results", "south"))
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

SPP0 <- union(fln, fls)
TAX <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(TAX) <- TAX$AOU
SPP <- sort(as.character(TAX$AOU)[TAX$map.pred])
compare_sets(SPP,SPP0)

STAGE <- list(
    veg =which(names(modsn) == "Space"),
    soil=which(names(modss) == "Space"))

Xclim <- model.matrix(fclim, kgrid[,,drop=FALSE])
colnames(Xclim) <- fixNames(colnames(Xclim))
r0 <- .read_raster_template()
rpa <- .make_raster(kgrid$pAspen, kgrid, r0)

for (spp in SPP) { # species START

    ## raster is on linear predictor scale
    cat(spp, which(spp==SPP), "/", length(SPP))
    if (TAX[spp,"veghf.north"]) {
        cat(" - veg")
        flush.console()
        fn <- file.path(ROOT, "out", "birds", "results", "north",
            paste0("birds_abmi-north_", spp, ".Rdata"))
        resn <- suppressWarnings(loadSPP(fn))
        estn <- suppressWarnings(getEst(resn, stage=STAGE$veg, na.out=FALSE, Xnn))
        ## north - current
        estnClim <- estn[,colnames(Xclim),drop=FALSE]
        logPNclim1 <- Xclim %*% estnClim[1,]
        rveg <- .make_raster(logPNclim1, kgrid, rpa)
        #writeRaster(rn, paste0("w:/reports/2017/results/birds/clim-veg/",
        #    TAX[spp, "sppid"], ".tif"), overwrite=TRUE)
    } else {
        rveg <- NULL
    }
    if (TAX[spp,"soilhf.south"]) {
        cat(" - soil")
        flush.console()
        fs <- file.path(ROOT, "out", "birds", "results", "south",
            paste0("birds_abmi-south_", spp, ".Rdata"))
        ress <- suppressWarnings(loadSPP(fs))
        ests <- suppressWarnings(getEst(ress, stage=STAGE$soil, na.out=FALSE, Xns))
        ## south - current
        pA <- ests[1,"pAspen"]
        estsClim <- ests[,colnames(Xclim),drop=FALSE]
        logPSclim1 <- Xclim %*% estsClim[1,]
        rs <- .make_raster(logPSclim1, kgrid, rpa)
        rsoil <- pA * rpa + rs
        #writeRaster(rs, paste0("w:/reports/2017/results/birds/clim-soil/",
        #    TAX[spp, "sppid"], ".tif"), overwrite=TRUE)
    } else {
        rsoil <- NULL
    }
    save(rveg, rsoil, file=paste0("w:/reports/2017/results/birds/spclim/",
        TAX[spp, "sppid"], ".RData"))
    cat("\n")
}

## pAspen raster
writeRaster(rpa, paste0("w:/reports/2017/data/pAspen.tif"), overwrite=TRUE)

load_common_data()
SP <- cure4insect:::.c4if$SP
PA <- cure4insect:::.c4if$CF$coef$paspen
SPP <- rownames(SP)[SP$taxon != "birds"]
SPPv <- rownames(SP)[SP$veghf.north & SP$taxon != "birds"]
SPPs <- rownames(SP)[SP$soilhf.south & SP$taxon != "birds"]
rpa <- raster("w:/reports/2017/data/pAspen.tif")

a <- "w:/Files from Ermias/Climate and Space prediction/"
for (i in SPP) {
    cat(i, which(i==SPP), "/", length(SPP))
    tx <- as.character(SP[i,"taxon"])
    if (i %in% SPPv) {
        cat(" - veg")
        flush.console()
        rveg <- raster(paste0(a, tx, "/clim-veg/", i, ".tif")) + 1 - 1
    } else {
        rveg <- NULL
    }
    if (i %in% SPPs) {
        cat(" - soil")
        flush.console()
        rs <- raster(paste0(a, tx, "/clim-soil/", i, ".tif")) + 1 - 1
        rsoil <- PA[i,"pAspen"] * rpa + rs
    } else {
        rsoil <- NULL
    }
    save(rveg, rsoil, file=paste0("w:/reports/2017/results/", tx, "/spclim/",
        i, ".RData"))
    cat("\n")
}







## ---



load("e:/peter/AB_data_v2017/data/analysis/site/veg-hf_siteCenter_v6-fixage0.Rdata")
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
vhf <- as.matrix(dd_1ha$veg_current)
vhf <- vhf[,colnames(vhf) != "CutBlocks"]
a <- vhf / rowSums(vhf)
a <- a[,rownames(tv)]
a <- groupSums(a, 2, tv$ETA_UseInAnalysis_Sector)
summary(a)


## processing vplant data for sample based non-native sector effects

library(mefa4)

load("e:/peter/AB_data_v2017/data/analysis/site/veg-hf_siteCenter_v6-fixage0.Rdata")
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")
vhf <- as.matrix(dd_1ha$veg_current)
vhf <- vhf[,colnames(vhf) != "CutBlocks"]
a <- vhf / rowSums(vhf)
a <- a[,rownames(tv)]
a <- groupSums(a, 2, tv$ETA_UseInAnalysis_Sector)
summary(a)

load("e:/peter/AB_data_v2017/data/misc/VPlant.RData")

str(dd)
Y <- as.matrix(dd[,6:ncol(dd)])
Y[Y>0] <- 1

X <- dd[,1:6]
compare_sets(rownames(X), rownames(a))
Sites <- intersect(rownames(X), rownames(a))

X <- X[Sites,]
Y <- Y[Sites,]
A <- a[Sites,]

sp <- read.csv("e:/peter/AB_data_v2017/data/misc/VPlant Species analysis name_April 2017.csv")
compare_sets(sp$Analysis_Name, colnames(Y))
Z <- sp[match(colnames(Y), sp$Analysis_Name),]
rownames(Z) <- Z$Analysis_Name
Z <- droplevels(Z[Z$RANK_NAME == "Species",])
Z$NN <- Z$Origin_Analysis != "Native"
Y <- Y[,rownames(Z)]

#save(X, Y, A, Z, file="e:/peter/AB_data_v2017/data/reporting/vpalnts-2017-11-16.Rdata")

hf <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-class.csv")
knitr::kable(hf[,c("HF_GROUP", "Sector")])

library(mefa4)
load("e:/peter/AB_data_v2017/data/reporting/Summary of VPlant abundance by sector.RData")
load("e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv0$Sector2 <- factor(ifelse(is.na(tv0$Sector), "NATIVE", as.character(tv0$Sector)),
    c("NATIVE", "Agriculture", "Energy", "Forestry", "Misc", "RuralUrban", "Transportation"))
VHF <- dd1km_pred[[1]]
tv0 <- tv0[colnames(VHF),]

CS <- colSums(groupSums(VHF/10^6, 2, tv0$Sector2))
CS <- CS / sum(CS)

Curr <- VP_SectorAbundance.Curr["Achillea.millefolium",]
Ref <- VP_SectorAbundance.Ref["Achillea.millefolium",]
Area <- 100*CS

## Sector effect plot from Dave
plot_sector_1 <- function(Curr, Ref, Area, main="") {
    sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
    sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")
    c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")

    total.effect <- (100 * (Curr - Ref) / sum(Ref))[sectors]
    unit.effect <- 100 * total.effect / Area[sectors]
    ymax <- ifelse(max(abs(unit.effect))<20,20,
        ifelse(max(abs(unit.effect))<50,50,round(max(abs(unit.effect))+50,-2)))
    ymin <- ifelse(ymax>50,min(-100,round(min(unit.effect)-50,-2)),-ymax)
    ymax <- max(ymax,max(unit.effect)+0.08*(max(unit.effect)-min(unit.effect,0)))
    ymin <- min(ymin,min(unit.effect)-0.08*(max(unit.effect,0)-min(unit.effect)))
    q <- barplot(unit.effect,
        width=Area[sectors],
        space=0,col=c1,border=c1,ylim=c(ymin,ymax),
        ylab="Unit effect (%)",xlab="Area (% of region)",
        xaxt="n",cex.lab=1.3,cex.axis=1.2,tcl=0.3,
        xlim=c(0,round(sum(Area[sectors])+1,0)),
        bty="n",col.axis="grey40",col.lab="grey40",las=2)

    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray88",border="gray88")
    x.at<-pretty(c(0,sum(Area[sectors])))
    axis(side=1,tck=1,at=x.at,lab=rep("",length(x.at)),col="grey95")
    y.at<-pretty(c(ymin,ymax),n=6)
    axis(side=2,tck=1,at=y.at,lab=rep("",length(y.at)),col="grey95")
    q <- barplot(unit.effect,
        width=Area[sectors],
        space=0,col=c1,border=c1,ylim=c(ymin,ymax),
        ylab="Unit effect (%)",xlab="Area (% of region)",
        xaxt="n",cex.lab=1.3,cex.axis=1.2,tcl=0.3,
        xlim=c(0,round(sum(Area[sectors])+1,0)),
        bty="n",col.axis="grey40",col.lab="grey40",las=2,add=TRUE)
    box(bty="l",col="grey40")
    #mtext(side=1,line=2,at=x.at,x.at,col="grey40",cex=1.2)
    axis(side=1,at=x.at,tcl=0.3,lab=rep("",length(x.at)),col="grey40",
        col.axis="grey40",cex.axis=1.2,las=1)
    abline(h=0,lwd=2,col="grey40")
    mtext(side=1,at=q+c(0,0,-1,0,+1),sector.names,col=c1,cex=1.3,
        adj=0.5,line=c(0.1,0.1,1.1,0.1,1.1))
    y <- unit.effect+0.025*(ymax-ymin)*sign(unit.effect)
    if (abs(y[3]-y[4])<0.05*(ymax-ymin))
        y[3:4]<-mean(y[3:4])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[3:4])]
    if (abs(y[4]-y[5])<0.05*(ymax-ymin))
        y[4:5]<-mean(y[4:5])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[4:5])]
    text(q,y,paste(ifelse(total.effect>0,"+",""),
        sprintf("%.1f",total.effect),"%",sep=""),col="darkblue",cex=1.4)
    mtext(side=3,line=1,at=0,adj=0, main, cex=1.4,col="grey40")
    invisible(rbind(total=total.effect, unit=unit.effect, area=Area[sectors]))
}
plot_sector_2 <- function(Curr, Ref, regional=TRUE, main="") {
    sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
    sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")
    c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")
    total.effect <- if (regional)
        100 * (Curr - Ref)/sum(Ref) else 100 * (Curr - Ref)/Ref
    total.effect <- total.effect[sectors]
    off <- 0.25
    a <- 1-0.5-off
    b <- 5+0.5+off
    ymax <- ifelse(max(abs(total.effect))<20,20,
        ifelse(max(abs(total.effect))<50,50,round(max(abs(total.effect))+50,-2)))
    ymin <- ifelse(ymax>50,min(-100,round(min(total.effect)-50,-2)),-ymax)
    ymax <- max(ymax,max(total.effect)+0.08*(max(total.effect)-min(total.effect,0)))
    ymin <- min(ymin,min(total.effect)-0.08*(max(total.effect,0)-min(total.effect)))
    yax <- pretty(c(ymin,ymax))
    op <- par(las=1, xpd = TRUE)
    on.exit(par(op))
    plot(0, type="n", xaxs="i", yaxs = "i", ylim=c(ymin,ymax), xlim=c(a, b),
        axes=FALSE, ann=FALSE)
    polygon(c(a,a,b,b), c(ymin, ymax, ymax, ymin), col="grey88", border="grey88")
    segments(x0=rep(a, length(yax)), x1=rep(b,length(yax)),y0=yax, col="white")
    axis(2, yax, paste0(ifelse(yax>0, "+", ""), yax), tick=FALSE)
    rug(yax, side=2, ticksize=0.01, col="grey40", quiet=TRUE)
    lines(c(a,a), c(ymin, ymax), col="grey40", lwd=1)
    for (i in 1:5) {
        h <- total.effect[i]
        polygon(c(i-0.5, i-0.5, i+0.5, i+0.5), c(0,h,h,0), col=c1[i], border=NA)
    }
    lines(c(a,b), c(0, 0), col="grey40", lwd=2)
    title(ylab=if (regional) "Regional sector effects (%)" else "Local sector effects (%)",
        cex=1.3, col="grey40")
    mtext(side=1,at=1:5,sector.names,col=c1,cex=1.3,adj=0.5,line=0.5)

    y <- total.effect+0.025*(ymax-ymin)*sign(total.effect)
    if (abs(y[3]-y[4])<0.05*(ymax-ymin))
        y[3:4]<-mean(y[3:4])+(c(-0.015,0.015)*(ymax-ymin))[rank(y[3:4])]
    text(1:5,y,paste(sprintf("%.1f",total.effect),"%",sep=""),col="darkblue",cex=1.2)
    mtext(side=3,line=1,at=0,adj=0, main, cex=1.4,col="grey40")
    invisible(total.effect)
}

plot_sector_1(Curr, Ref, Area, main="Common Yarrow")

op <- par(mfrow=c(1,2))
plot_sector_2(Curr, Ref, regional=TRUE, main="Common Yarrow")
plot_sector_2(Curr, Ref, regional=FALSE, main="Common Yarrow")
par(op)

png(paste0("combinedSE/", species_name, ".png")), width=600, height=600)
plot_sector_1(Curr, Ref, Area, main=species_name)
dev.off()

png(paste0("regionalSE/", species_name, ".png")), width=600, height=600)
plot_sector_2(Curr, Ref, regional=TRUE, main=species_name)
dev.off()

png(paste0("localSE/", species_name, ".png")), width=600, height=600)
plot_sector_2(Curr, Ref, regional=FALSE, main=species_name)
dev.off()

sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
SEreg <- 100 * (VP_SectorAbundance.Curr - VP_SectorAbundance.Ref) / rowSums(VP_SectorAbundance.Ref)
#A <- 100*CS
#SEunit <- t(100 * t(SEreg) / A)
SEloc <- 100 * (VP_SectorAbundance.Curr - VP_SectorAbundance.Ref) / VP_SectorAbundance.Ref
SEreg <- SEreg[,sectors]
#SEunit <- SEunit[,sectors]
SEloc <- SEloc[,sectors]

plot_sector_3 <- function(x, ylab="Sector effects (%)", type="kde",
breaks = "Sturges", ...) {
    type <- match.arg(type, c("kde", "fft", "hist"))
    if (!is.list(x))
        x <- as.data.frame(x)
    sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")
    sector.names <- c("Agriculture","Forestry","Energy","RuralUrban","Transport")
    c1 <- c("tan3","palegreen4","indianred3","skyblue3","slateblue2")
    ymin <- -100
    ymax <- 100
    off <- 0.25
    a <- 1-0.5-off
    b <- 5+0.5+off
    v <- 0.1
    yax <- pretty(c(ymin,ymax))
    op <- par(las=1)
    on.exit(par(op))
    plot(0, type="n", xaxs="i", yaxs = "i", ylim=c(ymin,ymax), xlim=c(a, b),
        axes=FALSE, ann=FALSE)
    polygon(c(a,a,b,b), c(ymin, ymax, ymax, ymin), col="grey88", border="grey88")
    segments(x0=rep(a, length(yax)), x1=rep(b,length(yax)),y0=yax, col="white")
    axis(2, yax, paste0(ifelse(yax>0, "+", ""), yax), tick=FALSE)
    rug(yax, side=2, ticksize=0.01, col="grey40", quiet=TRUE)
    lines(c(a,a), c(ymin, ymax), col="grey40", lwd=1)
    lines(c(a,b), c(0, 0), col="grey40", lwd=2)
    out <- list()
    for (i in 1:5) {
        xx <- sort(x[[i]])
        k <- xx <= ymax
        out[[i]] <- sum(!k)
        st <- boxplot.stats(xx)
        s <- st$stats
        k[which(!k)[1]] <- TRUE
        if (type == "kde")
            d <- KernSmooth::bkde(xx[k]) # uses Normal kernel
        if (type == "fft")
            d <- density(xx[k]) # uses FFT
        if (type == "hist") {
            h <- hist(xx[k], plot=FALSE, breaks=breaks)
            xv <- rep(h$breaks, each=2)
            yv <- c(0, rep(h$density, each=2), 0)
        } else {
            xv <- d$x
            yv <- d$y
            j <- xv >= min(xx) & xv <= max(xx)
            xv <- xv[j]
            yv <- yv[j]
        }
        yv <- 0.4 * yv / max(yv)
        polygon(c(-yv, rev(yv))+i, c(xv, rev(xv)), col=c1[i], border=c1[i])
        polygon(c(-v,-v,v,v)+i, s[c(2,4,4,2)], col="#40404080", border=NA)
        lines(c(-v,v)+i, s[c(3,3)], lwd=2, col="grey30")
    }
    title(ylab=ylab, cex=1.3, col="grey40")
    mtext(side=1,at=1:5,sector.names,col=c1,cex=1.3,adj=0.5,line=0.5)
    op <- par(xpd = TRUE)
    on.exit(par(op), add=TRUE)
    out <- unlist(out)
    points(1:5, rep(105, 5), pch=19,
        cex=ifelse(out==0, 0, 0.5+2*out/max(out)), col=c1)
    invisible(x)
}
op <- par(mfrow=c(3,2))
plot_sector_3(SEreg, ylab="Regional (total) sector effects (%)", type="kde")
title(main="Kernel density")
plot_sector_3(SEloc, ylab="Local (within footprint) sector effects (%)", type="kde")

plot_sector_3(SEreg, ylab="Regional (total) sector effects (%)", type="fft")
title(main="Fast Fourier transform based density")
plot_sector_3(SEloc, ylab="Local (within footprint) sector effects (%)", type="fft")

plot_sector_3(SEreg, ylab="Regional (total) sector effects (%)", type="hist",
    breaks=20)
title(main="Binning")
plot_sector_3(SEloc, ylab="Local (within footprint) sector effects (%)", type="hist",
    breaks=20)
par(op)

## load
library(mefa4)
load("e:/peter/AB_data_v2017/data/reporting/vpalnts-2017-11-16.Rdata")

## subset here to the region in question

P <- t(sapply(colnames(Y), function(i) colMeans(A * Y[,i])))
SE <- 100 * P[,colnames(P) != "NATIVE"] / P[,"NATIVE"]
SI <- 100 * (1 - colMeans(Y))

P1 <- rowSums(P[,colnames(P) != "NATIVE"])
P0 <- P[,"NATIVE"]
SIv <- 100 * pmin(P1, P0) / pmax(P1, P0)

#plot(SI[Z$NN], SIv[Z$NN])

tmp <- rbind(data.frame(SI=SI[Z$NN], SE[Z$NN,], NN=1),
    data.frame(SI=SI[!Z$NN], SE[!Z$NN,], NN=0))
tmp <- tmp[,colnames(tmp) != "Misc"]
#write.csv(tmp, file="e:/peter/AB_data_v2017/data/reporting/nn-plant-sector-effects.csv")

NN <- VP_SectorAbundance.Curr / rowSums(VP_SectorAbundance.Curr)
DD <- t(t(NN) / CS)

SE2 <- 100 * DD[,colnames(DD) != "Native"] / DD[,"Native"]
SE2 <- SE2[rownames(SE),]

op <- par(mfrow=c(2,3))
for (i in 2:6) {
plot(SE[Z$NN,i], SE2[Z$NN,i], xlab="Data based", ylab="Prediction based",
     ylim=c(0, 4000), xlim=c(0, 500), main=colnames(SE)[i])
abline(0,1,col="grey")
abline(h=100,v=100,col="grey", lty=2)
}
par(op)

op <- par(mfrow=c(2,1))
boxplot(SE[Z$NN,], ylim=c(0,2000))
abline(h=100)
boxplot(SE2[Z$NN,], ylim=c(0,10000))
abline(h=100)
par(op)





## producing some pdf's

Area <- 100*CS
pdf("sector-effect_local-regional.pdf", onefile=TRUE, width=21, height=7)
for (spp in rownames(Z)[!Z$NN]) {
    cat(spp, "\n")
    flush.console()
    Curr <- VP_SectorAbundance.Curr[spp,]
    Ref <- VP_SectorAbundance.Ref[spp,]
    op <- par(mfrow=c(1,3))
    plot_sector_1(Curr, Ref, Area, main=paste("Combined -", spp))
    plot_sector_2(Curr, Ref, regional=TRUE, main=paste("Regional (Total) -", spp))
    plot_sector_2(Curr, Ref, regional=FALSE, main=paste("Local -", spp))
    par(op)
}
dev.off()


## making some polys:

library(rgdal)
library(rgeos)
library(sp)

## shape files for boundary
crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/nsr")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
setwd(od)
AB <- spTransform(AB, crs)
AB0 <- gUnaryUnion(AB, rep(1, nrow(AB@data)))
AB0 <- as(AB0, "SpatialPolygonsDataFrame")
AB0@data <- data.frame(Province="Alberta")
100*object.size(AB0)/object.size(AB)
writeOGR(AB0, "AB_bound.geojson", layer="Alberta", driver="GeoJSON")

crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/osa")
OSA <- readOGR(".", "JOSM_Boundary_Sept2012") # rgdal
setwd(od)
OSA <- spTransform(OSA, crs)
OSA0 <- gUnaryUnion(OSA, rep(1, nrow(OSA@data)))
OSA0 <- as(OSA0, "SpatialPolygonsDataFrame")
OSA0@data <- data.frame(Region="OSA")
writeOGR(OSA0, "OSA_bound.geojson", layer="OSA", driver="GeoJSON")


input.gd <- "e:/peter/AB_data_v2017/data/raw/xy/OilSandsBoundaries.gdb"
## List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(input.gd)
fc_list # Determine what layers you wish to load
## Read the feature class
os.peaceriver <- readOGR(dsn = input.gd, layer = "Oilsand3RegionDissolve10TM")
os.area <- readOGR(dsn = input.gd, layer = "OilsandRegionDissolve10TM")
os.minable <- readOGR(dsn = input.gd, layer = "Oilsand_Mineable10TM")
plot(os.area)

## checking bootstrap distr

library(cure4insect)
library(pbapply)
library(mefa4)
opar <- set_options(path = "w:/reports")
getOption("cure4insect")
load_common_data()
subset_common_data(id=get_all_id(),
    species=get_all_species())

.c4if=cure4insect:::.c4if
.c4is=cure4insect:::.c4is
KT <- .c4if$KT

spp <- "Achillea.millefolium"
y <- load_species_data(spp)

fun <- function(y, what="cr", plot=TRUE) {
    if (what == "cr") {
        x1 <- y$SA.Curr
        x10 <- y$Curr.Boot
    }
    if (what == "rf") {
        x1 <- y$SA.Ref
        x10 <- y$Ref.Boot
    }
    x1 <- x1[match(rownames(KT), rownames(x1)),]
    x1[is.na(x1)] <- 0
    x1 <- groupMeans(x1, 1, KT$Row10_Col10)
    x1 <- rowSums(x1[rownames(x10),])
    cm <- colMeans(cbind(x1, x10))
    names(cm) <- c("km", paste0("B", 1:100))
    if (plot) {
        hist(cm[-1], col="lightgrey", border="darkgrey", xlab="Mean abundance",
            xlim=range(cm), main=paste(spp, what))
        abline(v=mean(cm[-1]), col=1, lty=2, lwd=2)
        abline(v=cm[1], col=2, lty=1, lwd=2)
        legend("topright", bty="n", col=1:2,lty=2:1, lwd=2,
            legend=c("Mean 10km B", "Mean 1km 1st"))
    }
    cm
}

SPP <- get_all_species()
set_options(verbose="0")
CR <- pbsapply(SPP, function(spp) {
    y <- load_species_data(spp)
    cr <- fun(y, what="cr", plot=FALSE)
    rf <- fun(y, what="rf", plot=FALSE)
    c(cr=cr, rf=rf)
})
CR <- t(CR)
write.csv(CR, file="abundance-comparison-1vs10km_updated.csv")

tb <- get_species_table()
cr1 <- CR[,1]
cr10 <- rowMeans(CR[,2:101])
rf1 <- CR[,102]
rf10 <- rowMeans(CR[,103:202])

par(mfrow=c(1,2))
plot(cr1, cr10, col=tb$taxon, xlim=c(0,1), ylim=c(0,1), main="cr")
abline(0,1,col="grey",lty=2)
legend("bottomright", col=1:6, pch=19, legend=levels(tb$taxon))
plot(rf1, rf10, col=tb$taxon, xlim=c(0,1), ylim=c(0,1), main="cr")
abline(0,1,col="grey",lty=2)
legend("bottomright", col=1:6, pch=19, legend=levels(tb$taxon))


f <- function(x1, x10) {
    x <- x[match(rownames(KT), rownames(x)),]
    x[is.na(x)] <- 0

}
Nc1 <- y$SA.Curr[match() rowSums(


## see how these compare
system.time(res <- report_all(cores=1))
#system.time(res <- report_all(cores=2))
#system.time(res <- report_all(cores=4))
## this is for testing only
#system.time(res <- .report_all_by1())
(set_options(opar)) # reset options

