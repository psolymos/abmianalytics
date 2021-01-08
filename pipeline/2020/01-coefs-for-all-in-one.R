## species lookup table to begin with
## include model performance matrics AUC R2

ROOT <- "s:/AB_data_v2020/Results/Results from Ermias"



sp <- list()
sp$lichens <- read.csv(file.path(ROOT, "Species lookup for Lichens 2020 - with decisions.csv"))
sp$mites <- read.csv(file.path(ROOT, "Species lookup for Mites 2020 - with decisions.csv"))
sp$mosses <- read.csv(file.path(ROOT, "Species lookup for Moss 2020 - with decisions.csv"))
sp$vplants <- read.csv(file.path(ROOT, "Species lookup for Plants 2020 - with decision.csv"))

ms <- list()
mn <- list()
f <- function(x) {
    e <- new.env()
    load(x, envir=e)
    e$ModValid.df
}
ms$lichens <- f(file.path(ROOT, "Model performance tables","Lichen_South_Model performance table.RData"))
mn$lichens <- f(file.path(ROOT, "Model performance tables","Lichen_North_Model performance table.RData"))
ms$mites <- f(file.path(ROOT, "Model performance tables","Mites_South_Model performance table.RData"))
mn$mites <- f(file.path(ROOT, "Model performance tables","Mites_North_Model performance table.RData"))
ms$mosses <- f(file.path(ROOT, "Model performance tables","Moss_South_Model performance table.RData"))
mn$mosses <- f(file.path(ROOT, "Model performance tables","Moss_North_Model performance table.RData"))
ms$vplants <- f(file.path(ROOT, "Model performance tables","VPlants_South_Model performance table.RData"))
mn$vplants <- f(file.path(ROOT, "Model performance tables","VPlants_North_Model performance table.RData"))

f <- function(x) {
    e <- new.env()
    load(x, envir=e)
    e$Lookup
}
st <- list()
st$lichens <- f("s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for Lichens 2020.RData")
st$mites <- f("s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for Mites 2020.RData")
st$mosses <- f("s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for Moss 2020.RData")
st$vplants <- f("s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for VPlants 2020.RData")


taxon <- "mites"
SPTAB <- list()
for (taxon in names(sp)) {

    St <- st[[taxon]]
    Sp <- sp[[taxon]]
    rownames(Sp) <- Sp$Analysis_Name
    Mn <- mn[[taxon]]
    Ms <- ms[[taxon]]
    Sp$ExcludeN <- grepl("NO", as.character(Sp$North_Model_use))
    Sp$ExcludeS <- grepl("NO", as.character(Sp$South_Model_use))

    SPP <- as.character(Sp$Analysis_Name)
    Sp <- Sp[SPP,]
    St <- St[SPP,]
    Mn <- Mn[SPP,]
    Ms <- Ms[SPP,]

    xx <- data.frame(
        SpeciesID=SPP,
        ScientificName=Sp$Scientific_Name_Analysis,
        TSNID=Sp$TSNID,
        CommonName=NA,
        ModelNorth=Sp$Analysis.North & !Sp$ExcludeN,
        ModelSouth=Sp$Analysis.South & !Sp$ExcludeS,
        UseavailNorth=Sp$UseAvailability.North,
        UseavailSouth=Sp$UseAvailability.South,
        Occurrences=Sp$Occurrences,
        nSites=Sp$nSites,
        SizeNorth=NA,
        SizeSouth=NA,
        Nonnative=FALSE,
        LinkHabitat="logit",
        LinkSpclim="logit",
        AUCNorth=Mn$AUC.All,
        AUCSouth=Ms$AUC.All,
        R2North=Mn$PseudoR2.All,
        R2South=Ms$PseudoR2.All,
        Comments=NA,
        Group=taxon)
    if (taxon == "vplants") {
        xx$CommonName <- Sp$Common.Name
        xx$Nonnative <- Sp$Origin != "Native"
    }
    rownames(xx) <- xx$SpeciesID
    xx$UseavailNorth[xx$ModelNorth] <- FALSE
    xx$UseavailSouth[xx$ModelSouth] <- FALSE


    SPTAB[[taxon]] <- xx
}

for (i in 1:length(SPTAB))
    for (j in 1:ncol(SPTAB[[i]]))
        if (is.factor(SPTAB[[i]][,j]))
            SPTAB[[i]][,j] <- as.character(SPTAB[[i]][,j])

## organizing results for mites, lichens, mosses, vplants
if (FALSE) {
ROOT <- "s:/AB_data_v2020/Results/Results from Ermias/Preliminary/"

TAXA <- c("lichens", "mites", "mosses", "vplants")
#taxon <- "mites"

COEFS <- list()
B <- 100

for (taxon in TAXA) {

    if (taxon == "lichens") {
        fn <- "Lichen_North Coefficient tables 2020.Rdata"
        fs <- "Lichen_South Coefficient tables 2020.Rdata"
    }
    if (taxon == "mites") {
        fn <- "Mite_North Coefficient tables 2020.Rdata"
        fs <- "Mite_South Coefficient tables 2020.Rdata"
    }
    if (taxon == "mosses") {
        fn <- "Moss_North Coefficient tables 2020.Rdata"
        fs <- "Moss_South Coefficient tables 2020.Rdata"
    }
    if (taxon == "vplants") {
        fn <- "VPlant_North Coefficient tables 2020.Rdata"
        fs <- "VPlant_South Coefficient tables 2020.Rdata"
    }

    #' load North and south coefs into their respective environments
    en <- new.env()
    load(paste0(ROOT, fn), envir=en)
    es <- new.env()
    load(paste0(ROOT, fs), envir=es)
    et <- new.env()

    #' # North
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    cfn <- en$Coef
    colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
    cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=cfn[,"Grass"])
    cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass"))]
    xCFn <- cbind(qlogis(cfn), en$Sc.coef)
    #' array: bootstrap x species x coefs
    CFn <- array(0, c(nrow(xCFn), ncol(xCFn), B),
        dimnames=list(rownames(xCFn), colnames(xCFn), NULL))
    for (i in 1:B) {
        cfn <- en$Coef
        colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
        cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=cfn[,"Grass"])
        cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass"))]
        ses <- en$Coef.se
        colnames(ses) <- gsub("BlackSpruce", "TreedBog", colnames(ses))
        ses <- cbind(ses,
            TreedFenR=ses[,"TreedFen"],
            TreedFen1=ses[,"TreedFen"],
            TreedFen2=ses[,"TreedFen"],
            TreedFen3=ses[,"TreedFen"],
            TreedFen4=ses[,"TreedFen"],
            TreedFen5=ses[,"TreedFen"],
            TreedFen6=ses[,"TreedFen"],
            TreedFen7=ses[,"TreedFen"],
            TreedFen8=ses[,"TreedFen"],
            Rural=ses[,"UrbInd"],
            Urban=ses[,"UrbInd"],
            Industrial=ses[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=ses[,"Grass"])
        ses <- ses[,!(colnames(ses) %in% c("TreedFen", "UrbInd", "Grass"))]
        if (i > 1) {
            se <- ses
            se[se > 5] <- 5
            eps <- rnorm(length(se), 0, se)
        } else {
            eps <- 0
        }
        CFn[,,i] <- cbind(qlogis(cfn) + eps, en$Sc.coef)
    }


    #' # South
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    cfs <- es$Coef
    cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
    cfs <- cfs[,colnames(cfs) != "UrbInd"]
    xCFs <- cbind(qlogis(cfs), pAspen=es$Coef.pAspen, es$Sc.coef)
    #' array: bootstrap x species x coefs
    CFs <- array(0, c(nrow(xCFs), ncol(xCFs), B),
        dimnames=list(rownames(xCFs), colnames(xCFs), NULL))
    for (i in 1:B) {
        cfs <- es$Coef
        cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        cfs <- cfs[,colnames(cfs) != "UrbInd"]
        ses <- es$Coef.se
        ses <- cbind(ses,
            Rural=ses[,"UrbInd"],
            Urban=ses[,"UrbInd"],
            Industrial=ses[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        ses <- ses[,colnames(ses) != "UrbInd"]
        if (i > 1) {
            se <- ses
            se[se > 5] <- 5
            eps <- rnorm(length(se), 0, se)
        } else {
            eps <- 0
        }
        CFs[,,i] <- cbind(qlogis(cfs) + eps, pAspen=es$Coef.pAspen, es$Sc.coef)
    }
    #range(CFn)
    #range(CFs)
    CFn[CFn < -10^4] <- -10^4
    CFn[CFn > 10^4] <- 10^4
    CFs[CFs < -10^4] <- -10^4
    CFs[CFs > 10^4] <- 10^4


    COEFS[[taxon]] <- list(north=CFn, south=CFs)

}

save(COEFS, file="s:/AB_data_v2020/Results/COEFS-EA.RData")
}

## bootstrap based coef tables (i.e. not normal SEs)

ROOT <- "s:/AB_data_v2020/Results/Results from Ermias/Boot coef/"

TAXA <- c("lichens", "mites", "mosses", "vplants")
#taxon <- "mites"

COEFS <- list()
B <- 250

for (taxon in TAXA) {

    if (taxon == "lichens") {
        fn <- "Lichen_North Bootstrap coefficents 2020.Rdata"
        fs <- "Lichen_South Bootstrap coefficents 2020.Rdata"
        #ft <- "s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for Lichens 2020.RData"
    }
    if (taxon == "mites") {
        fn <- "Mite_North Bootstrap coefficents 2020.Rdata"
        fs <- "Mite_South Bootstrap coefficents 2020.Rdata"
        #ft <- "s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for Mites 2020.RData"
    }
    if (taxon == "mosses") {
        fn <- "Moss_North Bootstrap coefficents 2020.Rdata"
        fs <- "Moss_South Bootstrap coefficents 2020.Rdata"
        #ft <- "s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for Moss 2020.RData"
    }
    if (taxon == "vplants") {
        fn <- "VPlant_North Bootstrap coefficents 2020.Rdata"
        fs <- "VPlant_South Bootstrap coefficents 2020.Rdata"
        #ft <- "s:/AB_data_v2020/Results/Results from Ermias/Species look up tables/Species lookup for VPlants 2020.RData"
    }

    #' load North and south coefs into their respective environments
    en <- new.env()
    load(paste0(ROOT, fn), envir=en)
    es <- new.env()
    load(paste0(ROOT, fs), envir=es)
    #load(ft, envir=et)
    #st <- et$Lookup
    st <- SPTAB[[taxon]]
    #if (taxon=="vplants")
    #    st$SpeciesID <- st$Analysis_Name

    #' # North
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    i <- 1
    cfn <- en$Coef.bs[,,i]
    colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
    cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0, Bare=0, SnowIce=0,
            GrassHerb=cfn[,"Grass"])
    cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass"))]
    xCFn <- cbind(qlogis(cfn), en$Sc.coef.bs[,,i])
    #' array: bootstrap x species x coefs
    CFn <- array(0, c(nrow(xCFn), ncol(xCFn), B),
        dimnames=list(rownames(xCFn), colnames(xCFn), NULL))
    for (i in 1:B) {
        cfn <- en$Coef.bs[,,i]
        colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
        cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0, Bare=0, SnowIce=0,
            GrassHerb=cfn[,"Grass"])
        cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass"))]
        CFn[,,i] <- cbind(qlogis(cfn), en$Sc.coef.bs[,,i])
    }


    #' # South
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    i <- 1
    cfs <- es$Coef.bs[,,i]
    cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            #HFor=0,
            Mine=0, MineV=0, Water=0)
    cfs <- cfs[,colnames(cfs) != "UrbInd"]
    xCFs <- cbind(qlogis(cfs), pAspen=es$Coef.pAspen.bs[[i]], es$Sc.coef.bs[,,i])
    #' array: bootstrap x species x coefs
    CFs <- array(0, c(nrow(xCFs), ncol(xCFs), B),
        dimnames=list(rownames(xCFs), colnames(xCFs), NULL))
    for (i in 1:B) {
        cfs <- es$Coef.bs[,,i]
        cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            #HFor=0,
            Mine=0, MineV=0, Water=0)
        cfs <- cfs[,colnames(cfs) != "UrbInd"]
        CFs[,,i] <- cbind(qlogis(cfs), pAspen=es$Coef.pAspen.bs[[i]], es$Sc.coef.bs[,,i])
    }
    CFn[CFn < -10^4] <- -10^4
    CFn[CFn > 10^4] <- 10^4
    CFs[CFs < -10^4] <- -10^4
    CFs[CFs > 10^4] <- 10^4

    cat(taxon, "--------------------\n")
    SPPn <- st$SpeciesID[st$ModelNorth]
    SPPs <- st$SpeciesID[st$ModelSouth]
    print(mefa4::compare_sets(dimnames(CFs)[[1]], SPPs))
    print(mefa4::compare_sets(dimnames(CFn)[[1]], SPPn))
    #print(head(st))
    cat("\n\n\n")


    COEFS[[taxon]] <- list(north=CFn[SPPn,,], south=CFs[SPPs,,], species=st)

}

save(COEFS, file="s:/AB_data_v2020/Results/COEFS-EAboot.RData")


#' Birds
#'
#' Organize coefs: compare land cover classes N/S
#' Organize spclim standardizations
#' Organize marginal & joint versions

library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")
source("~/repos/abmianalytics/pipeline/2020/00-functions.R")

ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds"

ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)
bt <- ee$tax
rm(ee)

en <- new.env()
load(file.path(ROOT, "data", "ab-birds-north-2020-09-23.RData"), envir=en)
es <- new.env()
load(file.path(ROOT, "data", "ab-birds-south-2020-12-04.RData"), envir=es)
Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v2020.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

load("s:/AB_data_v2020/Results/COEFS-EAboot.RData")

## all species
SPPn <- substr(list.files(file.path(ROOT, "out", "north")), 1, 4)
names(SPPn) <- bt[SPPn, "sppid"]
SPPs <- substr(list.files(file.path(ROOT, "out", "south")), 1, 4)
names(SPPs) <- bt[SPPs, "sppid"]
## same B as for other taxa
B <- dim(COEFS[[1]]$north)[3]

## take info from revisions:
## c=N+S combo
## n/s=N/S model only
## u=use avail (no model)
## o=exclude (passing through, extinct, bogus)
blist <- read.csv("~/repos/abmianalytics/pipeline/2020/birds-v2020.csv")
compare_sets(blist$common, bt$species)
blist$id <- bt$sppid[match(blist$common, bt$species)]
rownames(blist) <- blist$id
SPPn <- SPPn[names(SPPn) %in% rownames(blist)[blist$show %in% c("c", "n")]]
SPPs <- SPPs[names(SPPs) %in% rownames(blist)[blist$show %in% c("c", "s")]]


AUCNorth <- list()
for (spp in SPPn) {
    cat(spp, "N\n");flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    yn <- as.numeric(en$YY[,spp])
    off <- if (spp %in% colnames(en$OFF))
        en$OFF[,spp] else en$OFFmean
    lamn <- exp(predict_with_SSH(resn, Xn, en$SSH, stage="Space") + off)
    rocn <- simple_roc(ifelse(yn > 0, 1, 0), rowMeans(lamn))
    aucn <- simple_auc(rocn)
    AUCNorth[[spp]] <- aucn
}
AUCNorth <- unlist(AUCNorth)
save(AUCNorth, file="s:/AB_data_v2020/Results/BIRDS-North-AUC.RData")

AUCSouth <- list()
for (spp in SPPs) {
    cat(spp, "S\n");flush.console()
    ress <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    ys <- as.numeric(es$YY[,spp])
    off <- if (spp %in% colnames(es$OFF))
        es$OFF[,spp] else es$OFFmean
    lams <- exp(predict_with_SSH(ress, Xs, es$SSH, stage="Space") + off)
    rocs <- simple_roc(ifelse(ys > 0, 1, 0), rowMeans(lams))
    aucs <- simple_auc(rocs)
    AUCSouth[[spp]] <- aucs
}
AUCSouth <- unlist(AUCSouth)
save(AUCSouth, file="s:/AB_data_v2020/Results/BIRDS-South-AUC.RData")

cfn <- list()
for (spp in SPPn) {
    cat(spp, "\n")
    flush.console()

    res <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))

    est1 <- suppressWarnings(get_coef(res, Xn, stage="ARU", na.out=FALSE))
    est2 <- suppressWarnings(get_coef(res, Xn, stage="Space", na.out=FALSE))

    BB <- min(B, nrow(est1))
    cf1 <- sapply(1:BB, function(i) get_coef_north(est1, subset=i))
    cf2 <- sapply(1:BB, function(i) get_coef_north(est2, subset=i))

    cfn[[spp]] <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}
save(cfn, file="s:/AB_data_v2020/Results/BIRDS-North-rev1.RData")

cfs <- list()
for (spp in SPPs) {
    cat(spp, "\n")
    flush.console()

    res <- load_species(file.path(ROOT, "out", "south", paste0(spp, ".RData")))
    est1 <- suppressWarnings(get_coef(res, Xs, stage="ARU", na.out=FALSE))
    est2 <- suppressWarnings(get_coef(res, Xs, stage="Space", na.out=FALSE))

    BB <- min(B, nrow(est1))
    cf1 <- sapply(1:BB, function(i) get_coef_south(est1, subset=i))
    cf2 <- sapply(1:BB, function(i) get_coef_south(est2, subset=i))

    cfs[[spp]] <- list(estARU=est1, estSpace=est2, coefARU=cf1, coefSpace=cf2)
}

save(cfs, file="s:/AB_data_v2020/Results/BIRDS-South-rev2.RData")

if (FALSE) {
#for (spp in union(names(SPPn), names(SPPs))) {
for (spp in union(SPPn, SPPs)) {
    file.copy(paste0("s:/AB_data_v2020/Results/web1/birds/", spp, "/map.png"),
        paste0("s:/AB_data_v2020/Results/web1/tmp/", spp, ".png"),
        overwrite=TRUE)
}
#for (spp in names(SPPs)) {
for (spp in SPPs) {
    file.copy(paste0("s:/AB_data_v2020/Results/web1/birds/", spp, "/soilhf.png"),
        paste0("s:/AB_data_v2020/Results/web1/soil/", spp, ".png"),
        overwrite=TRUE)
}
#for (spp in names(SPPs)) {
for (spp in SPPs) {
    file.copy(paste0("s:/AB_data_v2020/Results/web1/birds/", spp, "/veghf.png"),
        paste0("s:/AB_data_v2020/Results/web1/veg/", spp, ".png"),
        overwrite=TRUE)
}
}

## combine birds and other taxa

library(mefa4)
source("~/repos/abmianalytics/pipeline/2020/00-functions.R")

## species names etc
ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds"
ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)
TAX <- ee$tax
#rm(ee)

#load("s:/AB_data_v2020/Results/COEFS-EA.RData")
load("s:/AB_data_v2020/Results/COEFS-EAboot.RData")
load("s:/AB_data_v2020/Results/BIRDS-North-rev1.RData")
load("s:/AB_data_v2020/Results/BIRDS-South-rev2.RData")
load("s:/AB_data_v2020/Results/BIRDS-North-AUC.RData")
load("s:/AB_data_v2020/Results/BIRDS-South-AUC.RData")

if (FALSE) { # testing CC effects
spp <- "OSFL"
h <- c("CCWhiteSpruceR", "CCWhiteSpruce1",
"CCWhiteSpruce2", "CCWhiteSpruce3", "CCWhiteSpruce4", "CCPineR",
"CCPine1", "CCPine2", "CCPine3", "CCPine4", "CCDeciduousR", "CCDeciduous1",
"CCDeciduous2", "CCDeciduous3", "CCDeciduous4", "CCMixedwoodR",
"CCMixedwood1", "CCMixedwood2", "CCMixedwood3", "CCMixedwood4")
h1 <- gsub("CC", "", h)
m  <- matrix(rowMeans(cfn[[spp]]$coefARU)[h], 5, 4)
m1 <- matrix(rowMeans(cfn[[spp]]$coefARU)[h1], 5, 4)
plot(0, type="n", xlim=c(0,4), ylim=c(0,max(m,m1)))
for (i in 1:4) {
    lines(0:4, m[,i], col=i, lty=2)
    lines(0:4, m1[,i], col=i, lty=1)
}
}

blist <- read.csv("~/repos/abmianalytics/pipeline/2020/birds-v2020.csv")
compare_sets(blist$common, TAX$species)
TAX$type <- blist$show[match(TAX$species, blist$common)]

#SPPn <- as.character(TAX$sppid[match(names(cfn), TAX$code)])
#SPPs <- as.character(TAX$sppid[match(names(cfs), TAX$code)])
SPPn <- names(ALLBIRDSPP$north[match(names(cfn), ALLBIRDSPP$north)])
SPPs <- names(ALLBIRDSPP$south[match(names(cfs), ALLBIRDSPP$south)])
rownames(TAX) <- TAX$sppid
TAX <- TAX[sort(union(SPPn,SPPs)),]
names(cfs) <- SPPs
names(cfn) <- SPPn


## organizing bird coefs

cns <- rownames(COEFS$mites$south[1,,])
cns0 <- cns[1:(which(cns=="Intercept")-2)]
cns1 <- cns[(which(cns=="Intercept")-1):length(cns)]
cnn <- rownames(COEFS$mites$north[1,,])
cnn0 <- cnn[1:(which(cnn=="Intercept")-1)]
cnn1 <- cnn[which(cnn=="Intercept"):length(cnn)]

cnb <- c("pWater_KM", "pWater2_KM", "xPET",
    "xMAT", "xAHM", "xFFP", "xMAP", "xMWMT", "xMCMT", "xY", "xX",
    "xY2", "xX2", "xFFP:xMAP", "xMAP:xPET", "xAHM:xMAT", "xX:xY")

BMAX <- 250

# south
spp <- 1
B <- min(BMAX, nrow(cfs[[spp]]$estSpace))
cfsm <- rbind(log(cfs[[spp]]$coefARU), pAspen=cfs[[spp]]$estARU[1:B,"pAspen"])
cfsj <- rbind(log(cfs[[spp]]$coefSpace), t(cfs[[spp]]$estSpace[1:B,c("pAspen", cnb)]))
cfsm[cfsm > 10^4] <- 10^4
cfsm[cfsm < -10^4] <- -10^4
cfsj[cfsj > 10^4] <- 10^4
cfsj[cfsj < -10^4] <- -10^4
if (ncol(cfsm) < BMAX) {
    b <- sample.int(ncol(cfsm), BMAX-ncol(cfsm), replace=TRUE)
    cfsm <- cbind(cfsm, cfsm[,b])
    cfsj <- cbind(cfsj, cfsj[,b])
}
compare_sets(rownames(cfsm),dimnames(COEFS$lichens$south)[[2]])

CFsm <- array(0, c(length(cfs), nrow(cfsm), BMAX))
dimnames(CFsm) <- list(names(cfs), rownames(cfsm), NULL)
CFsj <- array(0, c(length(cfs), nrow(cfsj), BMAX))
dimnames(CFsj) <- list(names(cfs), rownames(cfsj), NULL)

for (spp in names(cfs)) {
    B <- min(BMAX, nrow(cfs[[spp]]$estSpace))
    cfsm <- rbind(log(cfs[[spp]]$coefARU), pAspen=cfs[[spp]]$estARU[1:B,"pAspen"])
    cfsj <- rbind(log(cfs[[spp]]$coefSpace), t(cfs[[spp]]$estSpace[1:B,c("pAspen", cnb)]))
    cfsm[cfsm > 10^4] <- 10^4
    cfsm[cfsm < -10^4] <- -10^4
    cfsj[cfsj > 10^4] <- 10^4
    cfsj[cfsj < -10^4] <- -10^4
    if (ncol(cfsm) < BMAX) {
        b <- sample.int(ncol(cfsm), BMAX-ncol(cfsm), replace=TRUE)
        cfsm <- cbind(cfsm, cfsm[,b])
        cfsj <- cbind(cfsj, cfsj[,b])
    }
    CFsm[spp,,] <- cfsm
    CFsj[spp,,] <- cfsj
}

# north
spp <- 1
B <- min(BMAX, nrow(cfn[[spp]]$estSpace))
cfnm <- rbind(log(cfn[[spp]]$coefARU))
cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[1:B,cnb]))
cfnm[cfnm > 10^4] <- 10^4
cfnm[cfnm < -10^4] <- -10^4
cfnj[cfnj > 10^4] <- 10^4
cfnj[cfnj < -10^4] <- -10^4
if (ncol(cfnm) < BMAX) {
    b <- sample.int(ncol(cfnm), BMAX-ncol(cfnm), replace=TRUE)
    cfnm <- cbind(cfnm, cfnm[,b])
    cfnj <- cbind(cfnj, cfnj[,b])
}
compare_sets(rownames(cfnm),dimnames(COEFS$lichens$north)[[2]])

CFnm <- array(0, c(length(cfn), nrow(cfnm), BMAX))
dimnames(CFnm) <- list(names(cfn), rownames(cfnm), NULL)
CFnj <- array(0, c(length(cfn), nrow(cfnj), BMAX))
dimnames(CFnj) <- list(names(cfn), rownames(cfnj), NULL)

for (spp in names(cfn)) {
    B <- min(BMAX, nrow(cfn[[spp]]$estSpace))
    cfnm <- rbind(log(cfn[[spp]]$coefARU))
    cfnj <- rbind(log(cfn[[spp]]$coefSpace), t(cfn[[spp]]$estSpace[1:B,cnb]))
    cfnm[cfnm > 10^4] <- 10^4
    cfnm[cfnm < -10^4] <- -10^4
    cfnj[cfnj > 10^4] <- 10^4
    cfnj[cfnj < -10^4] <- -10^4
    if (ncol(cfnm) < BMAX) {
        b <- sample.int(ncol(cfnm), BMAX-ncol(cfnm), replace=TRUE)
        cfnm <- cbind(cfnm, cfnm[,b])
        cfnj <- cbind(cfnj, cfnj[,b])
    }
    CFnm[spp,,] <- cfnm
    CFnj[spp,,] <- cfnj
}


#COEFS$lichens$north <- COEFS$lichens$north[dimnames(COEFS$lichens$north)[[1]] != "Do.not.analyze",,]
#COEFS$lichens$south <- COEFS$lichens$south[dimnames(COEFS$lichens$south)[[1]] != "Do.not.analyze",,]

## Bird lookup table here

uu <- ee$tax
uu$keep <- as.character(uu$species) %in%
    as.character(blist[blist$show != "o","common"])
uu <- uu[uu$keep,]
uu$show <- blist$show[match(uu$species, blist$common)]
uu <- uu[!is.na(uu$show),]
table(uu$show)

uu$AUCn <- AUCNorth[match(uu$code, names(AUCNorth))]
uu$AUCs <- AUCSouth[match(uu$code, names(AUCSouth))]
uu$pocc <- colSums(ee$yy[,rownames(uu)] > 0)

spb <- data.frame(
        SpeciesID=uu$sppid,
        ScientificName=uu$scinam,
        TSNID=NA,
        CommonName=uu$species,
        ModelNorth=uu$show %in% c("c", "n"),
        ModelSouth=uu$show %in% c("c", "s"),
        UseavailNorth=NA,
        UseavailSouth=NA,
        Occurrences=uu$pocc,
        nSites=NA,
        SizeNorth=NA,
        SizeSouth=NA,
        Nonnative=FALSE,
        LinkHabitat="log",
        LinkSpclim="log",
        AUCNorth=uu$AUCn,
        AUCSouth=uu$AUCs,
        R2North=NA,
        R2South=NA,
        Comments=uu$code,
        Group="birds")
for (j in 1:ncol(spb))
    if (is.factor(spb[,j]))
        spb[,j] <- as.character(spb[,j])
rownames(spb) <- spb$SpeciesID

COEFS$birds <- list(
    north=list(marginal=CFnm, joint=CFnj),
    south=list(marginal=CFsm, joint=CFsj),
    species=spb)

save(COEFS,
    file="s:/AB_data_v2020/Results/COEFS-ALL.RData")



## habitat elements

library(mefa4)
ROOT <- "s:/AB_data_v2020/Results/Habitat elements May 2020"
load("s:/AB_data_v2020/Results/COEFS-EAboot.RData")
cfn <- dimnames(COEFS$mites$north)[[2]]
cfn0 <- cfn[1:(which(cfn=="Intercept")-1)]
cfn1 <- cfn[!(cfn %in% cfn0)]
cfs <- dimnames(COEFS$mites$south)[[2]]
cfs0 <- cfs[1:(which(cfs=="Intercept")-1)]
cfs1 <- cfs[!(cfs %in% cfs0)]
rm(COEFS)

t1 <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(t1) <- t1[,1]
t2 <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(t2) <- t2[,1]
compare_sets(rownames(t1), rownames(t2))
setdiff(rownames(t2), rownames(t1))
tt <- data.frame(New=t2[,2], Old=t1[match(rownames(t2), rownames(t1)),c("UseInAnalysisFine")])
rownames(tt) <- rownames(t2)

#CS <- c(
#    Loamy="Productive",
#    SandyLoam="Productive",
#    RapidDrain="RapidDrain",
#    ClaySub="Clay",
#    ThinBreak="RapidDrain",
#    Blowout="Saline",
#    Other="SoilWetland", # "SoilWetland" "Lotic"
#    Crop="Crop",
#    TameP="TameP",
#    RoughP="RoughP",
#    Wellsites="Wells",
#    EnSeismic="SoftLin",
#    EnSoftLin="SoftLin",
#    TrSoftLin="SoftLin",
#    HardLin="HardLin",
#    Rural="RuralResInd",
#    Urban="UrbInd",
#    Industrial="RuralResInd")
CN <- c(
    WhiteSpruceR="WhiteSpruce_0-10",
    WhiteSpruce1="WhiteSpruce_10-20",
    WhiteSpruce2="WhiteSpruce_20-40",
    WhiteSpruce3="WhiteSpruce_40-60",
    WhiteSpruce4="WhiteSpruce_60-80",
    WhiteSpruce5="WhiteSpruce_80-100",
    WhiteSpruce6="WhiteSpruce_100-120",
    WhiteSpruce7="WhiteSpruce_120-140",
    WhiteSpruce8="WhiteSpruce_140+",
    PineR="Pine_0-10",
    Pine1="Pine_10-20",
    Pine2="Pine_20-40",
    Pine3="Pine_40-60",
    Pine4="Pine_60-80",
    Pine5="Pine_80-100",
    Pine6="Pine_100-120",
    Pine7="Pine_120-140",
    Pine8="Pine_140+",
    DeciduousR="Deciduous_0-10",
    Deciduous1="Deciduous_10-20",
    Deciduous2="Deciduous_20-40",
    Deciduous3="Deciduous_40-60",
    Deciduous4="Deciduous_60-80",
    Deciduous5="Deciduous_80-100",
    Deciduous6="Deciduous_100-120",
    Deciduous7="Deciduous_120-140",
    Deciduous8="Deciduous_140+",
    MixedwoodR="Mixedwood_0-10",
    Mixedwood1="Mixedwood_10-20",
    Mixedwood2="Mixedwood_20-40",
    Mixedwood3="Mixedwood_40-60",
    Mixedwood4="Mixedwood_60-80",
    Mixedwood5="Mixedwood_80-100",
    Mixedwood6="Mixedwood_100-120",
    Mixedwood7="Mixedwood_120-140",
    Mixedwood8="Mixedwood_140+",
    TreedBogR="BlackSpruce_0-10",
    TreedBog1="BlackSpruce_10-20",
    TreedBog2="BlackSpruce_20-40",
    TreedBog3="BlackSpruce_40-60",
    TreedBog4="BlackSpruce_60-80",
    TreedBog5="BlackSpruce_80-100",
    TreedBog6="BlackSpruce_100-120",
    TreedBog7="BlackSpruce_120-140",
    TreedBog8="BlackSpruce_140+",
    TreedFenR="TreedFen",
    TreedFen1="TreedFen",
    TreedFen2="TreedFen",
    TreedFen3="TreedFen",
    TreedFen4="TreedFen",
    TreedFen5="TreedFen",
    TreedFen6="TreedFen",
    TreedFen7="TreedFen",
    TreedFen8="TreedFen",
    CCWhiteSpruceR="CCWhiteSpruce_0-10",
    CCWhiteSpruce1="CCWhiteSpruce_10-20",
    CCWhiteSpruce2="CCWhiteSpruce_20-40",
    CCWhiteSpruce3="CCWhiteSpruce_40-60",
    CCWhiteSpruce4="CCWhiteSpruce_60-80",
    CCPineR="CCPine_0-10",
    CCPine1="CCPine_10-20",
    CCPine2="CCPine_20-40",
    CCPine3="CCPine_40-60",
    CCPine4="CCPine_60-80",
    CCDeciduousR="CCDeciduous_0-10",
    CCDeciduous1="CCDeciduous_10-20",
    CCDeciduous2="CCDeciduous_20-40",
    CCDeciduous3="CCDeciduous_40-60",
    CCDeciduous4="CCDeciduous_60-80",
    CCMixedwoodR="CCMixedwood_0-10",
    CCMixedwood1="CCMixedwood_10-20",
    CCMixedwood2="CCMixedwood_20-40",
    CCMixedwood3="CCMixedwood_40-60",
    CCMixedwood4="CCMixedwood_60-80",
    TreedSwamp="TreeShrubSwamp",
    ShrubbySwamp="TreeShrubSwamp",
    ShrubbyBog="TreeShrubSwamp",
    ShrubbyFen="TreeShrubSwamp",
    GraminoidFen="NonTreeFenMarsh",
    Marsh="NonTreeFenMarsh",
    Shrub="Shrub",
    GrassHerb="Grass",
    Crop="Crop",
    TameP="TameP",
    RoughP="RoughP",
    Wellsites="Wells",
    EnSeismic="SoftLin",
    EnSoftLin="SoftLin",
    TrSoftLin="SoftLin",
    HardLin="HardLin",
    Rural="RuralResInd",
    Urban="UrbInd",
    Industrial="RuralResInd",
    Water="Water",
    Bare="Bare")

st <- read.csv("~/repos/abmianalytics/lookup/habitatelements-2020-lookup.csv")
rownames(st) <- st[,1]
st$link <- ifelse(st$LinkHabitat == "Binom", "logit", "log")
st$link[st$LinkHabitat == "Normal"] <- "identity"
colnames(st)[1] <- "sppid"

## south
es <- new.env()
#load(file.path(ROOT,
#    "Habitat coefficients South May 2020 OFFICIAL coefficients.Rdata"), envir=es)
load(file.path("s:/AB_data_v2020/Results/Habitat elements Oct 2020",
    "Habitat coefficients South Oct 2020 OFFICIAL coefficients.Rdata"), envir=es)


names(es)
#l10s <- read.csv(file.path(ROOT,
#    "Linear 10pc figure values for habitat South.csv"))
pA <- read.csv(file.path(ROOT,
    "Habitat coefficients South May 2020 pApen coefficients.csv"))
pA <- structure(pA$pAspen, names=as.character(pA$Sp))
xs0 <- es$Coef.official
colnames(xs0)[colnames(xs0) == "Well"] <- "Wellsites"
xs0 <- cbind(xs0, Mine=0, MineV=0, Water=0)
xs0 <- xs0[,cfs0[cfs0 != "pAspen"]]


xsl0 <- es$Coef.official.lci
colnames(xsl0)[colnames(xsl0) == "Well"] <- "Wellsites"
xsl0 <- cbind(xsl0, Mine=0, MineV=0, Water=0)
xsl0 <- xsl0[,cfs0[cfs0 != "pAspen"]]

xsu0 <- es$Coef.official.uci
colnames(xsu0)[colnames(xsu0) == "Well"] <- "Wellsites"
xsu0 <- cbind(xsu0, Mine=0, MineV=0, Water=0)
xsu0 <- xsu0[,cfs0[cfs0 != "pAspen"]]

rownames(xs0)[rownames(xs0) == "CanopyCover"] <- "CanopyClosure"
rownames(xsl0)[rownames(xsl0) == "CanopyCover"] <- "CanopyClosure"
rownames(xsu0)[rownames(xsu0) == "CanopyCover"] <- "CanopyClosure"

## link function comes here
for (i in rownames(xs0)) {
    FUN <- binomial(st[i,"link"])$linkfun
    xs0[i,] <- FUN(xs0[i,])
    xsl0[i,] <- FUN(xsl0[i,])
    xsu0[i,] <- FUN(xsu0[i,])
}

names(pA)[names(pA) == "CanopyCover"] <- "CanopyClosure"
pA <- pA[rownames(xs0)]
all(names(pA)==rownames(xs0))

xscl <- es$Res.coef.official[,cfs1]
rownames(xscl)[rownames(xscl) == "CanopyCover"] <- "CanopyClosure"
xscl <- xscl[rownames(xs0),]

CFs <- array(NA, c(nrow(xs0), length(cfs), 3))
dimnames(CFs) <- list(rownames(xs0), cfs, c("Estimate", "LCL", "UCL"))
CFs[,cfs,"Estimate"] <- cbind(xs0, pAspen=pA, xscl)[,cfs]
CFs[,cfs0[cfs0 != "pAspen"],"LCL"] <- xsl0[,cfs0[cfs0 != "pAspen"]]
CFs[,cfs0[cfs0 != "pAspen"],"UCL"] <- xsu0[,cfs0[cfs0 != "pAspen"]]
CFs[!is.na(CFs) & CFs < -10^4] <- -10^4
CFs[!is.na(CFs) & CFs > 10^4] <- 10^4

## north

en <- new.env()
load(file.path(ROOT,
    "Habitat coefficients North May 2020 OFFICIAL coefficients ALL.Rdata"), envir=en)
names(en)

xn <- en$Coef.official.ALL
xnl <- en$Coef.official.ALL.lci
xnu <- en$Coef.official.ALL.uci

rownames(xn)[rownames(xn) == "CanopyCover"] <- "CanopyClosure"
rownames(xnl)[rownames(xnl) == "CanopyCover"] <- "CanopyClosure"
rownames(xnu)[rownames(xnu) == "CanopyCover"] <- "CanopyClosure"

rownames(xn)[rownames(xn) == "Live Decid BA"] <- "Live Deciduous BA"
rownames(xnl)[rownames(xnl) == "Live Decid BA"] <- "Live Deciduous BA"
rownames(xnu)[rownames(xnu) == "Live Decid BA"] <- "Live Deciduous BA"

xn0 <- xn[,CN]
colnames(xn0) <- names(CN)
xn0 <- cbind(xn0, Mine=0, MineV=0, SnowIce=0)

xnl0 <- xnl[,CN]
colnames(xnl0) <- names(CN)
xnl0 <- cbind(xnl0, Mine=0, MineV=0, SnowIce=0)

xnu0 <- xnu[,CN]
colnames(xnu0) <- names(CN)
xnu0 <- cbind(xnu0, Mine=0, MineV=0, SnowIce=0)

## link function comes here
for (i in rownames(xn0)) {
    FUN <- binomial(st[i,"link"])$linkfun
    xn0[i,] <- FUN(xn0[i,])
    xnl0[i,] <- FUN(xnl0[i,])
    xnu0[i,] <- FUN(xnu0[i,])
}

xncl <- en$Res.coef.official[,cfn1]
xncl <- xncl[match(rownames(xn0), rownames(xncl)),]
rownames(xncl) <- rownames(xn0)

CFn <- array(NA, c(nrow(xn), length(cfn), 3))
dimnames(CFn) <- list(rownames(xn), cfn, c("Estimate", "LCL", "UCL"))
CFn[,,"Estimate"] <- cbind(xn0, xncl)[,cfn]
CFn[,cfn0[cfn0 != "pAspen"],"LCL"] <- xnl0[,cfn0[cfn0 != "pAspen"]]
CFn[,cfn0[cfn0 != "pAspen"],"UCL"] <- xnu0[,cfn0[cfn0 != "pAspen"]]
CFn[!is.na(CFn) & CFn < -10^4] <- -10^4
CFn[!is.na(CFn) & CFn > 10^4] <- 10^4

dimnames(CFs)[[1]] <- gsub(" ", "", dimnames(CFs)[[1]])
dimnames(CFn)[[1]] <- gsub(" ", "", dimnames(CFn)[[1]])
rownames(st) <- gsub(" ", "",rownames(st))
dimnames(CFs)[[1]] <- gsub("_", ".", dimnames(CFs)[[1]])
dimnames(CFn)[[1]] <- gsub("_", ".", dimnames(CFn)[[1]])
rownames(st) <- gsub("_", ".",rownames(st))
st$sppid <- rownames(st)

CFs[,"pAspen",2] <- CFs[,"pAspen",3] <- CFs[,"pAspen",1]

COEFS2 <- list(habitats = list(north=CFn, south=CFs, species=st))

## mammals (camera)

library(mefa4)
load("s:/AB_data_v2020/Results/COEFS-EAboot.RData")
cfn <- dimnames(COEFS$mites$north)[[2]]
cfn0 <- cfn[1:(which(cfn=="Intercept")-1)]
cfn1 <- cfn[!(cfn %in% cfn0)]
cfs <- dimnames(COEFS$mites$south)[[2]]
cfs0 <- cfs[1:(which(cfs=="Intercept")-1)]
cfs1 <- cfs[!(cfs %in% cfs0)]
rm(COEFS)

ROOT <- "s:/AB_data_v2020/Results"
st <- read.csv(file.path(ROOT,
    "Camera mammal models revised June 2020", "Mammal header table June 2020.csv"))
pA <- read.csv(file.path(ROOT,
    "Camera mammal models revised June 2020", "South", "Mammal coefficients South June 2020 Best model pApen coefficients.csv"))
rownames(pA) <- pA[,1]
en <- new.env()
es <- new.env()
load(file.path(ROOT,
    "Camera mammal models revised June 2020", "North",
    "Mammal coefficients North June 2020 Best model OFFICIAL coefficients.Rdata"
), envir=en)
load(file.path(ROOT,
    "Camera mammal models revised June 2020", "South",
    "Mammal coefficients South June 2020 Best model OFFICIAL coefficients.Rdata"
), envir=es)
names(en)

# south
CS <- c(
    Loamy="Productive",
    SandyLoam="Productive",
    RapidDrain="RapidDrain",
    ClaySub="Clay",
    ThinBreak="RapidDrain",
    Blowout="Saline",
    Other="SoilWetland", # "SoilWetland" "Lotic"
    Crop="Crop",
    TameP="TameP",
    RoughP="RoughP",
    Wellsites="Well",
    EnSeismic="EnSeismic",
    EnSoftLin="EnSoftLin",
    TrSoftLin="TrSoftLin",
    HardLin="HardLin",
    Rural="Rural",
    Urban="Urban",
    Industrial="Industrial")

names(es)
xs <- es$Coef.official
xs0 <- xs[,CS]
colnames(xs0) <- names(CS)
xs0 <- cbind(xs0, Mine=0, MineV=0, Water=0)
xs0[,"Other"] <- 0.5 * (xs[,"SoilWetland"] + xs[,"Lotic"])


xsl <- es$Coef.official.lci
xsl0 <- xsl[,CS]
colnames(xsl0) <- names(CS)
xsl0 <- cbind(xsl0, Mine=0, MineV=0, Water=0)
xsl0[,"Other"] <- 0.5 * (xsl[,"SoilWetland"] + xsl[,"Lotic"])

xsu <- es$Coef.official.uci
xsu0 <- xsu[,CS]
colnames(xsu0) <- names(CS)
xsu0 <- cbind(xsu0, Mine=0, MineV=0, Water=0)
xsu0[,"Other"] <- 0.5 * (xsu[,"SoilWetland"] + xsu[,"Lotic"])

FUN <- log
xs0 <- FUN(xs0)
xsl0 <- FUN(xsl0)
xsu0 <- FUN(xsu0)

## pa and agp ???
pA <- pA[rownames(xs), ]

xscl <- es$Res.coef.official[rownames(xs),cfs1]

CFs <- array(NA, c(nrow(xs), length(cfs), 3))
dimnames(CFs) <- list(rownames(xs), cfs, c("Estimate", "LCL", "UCL"))
CFs[,,"Estimate"] <- cbind(xs0, pAspen=pA[,"pAspen.agp"], xscl)[,cfs]
CFs[,cfs0[cfs0 != "pAspen"],"LCL"] <- xsl0[,cfs0[cfs0 != "pAspen"]]
CFs[,cfs0[cfs0 != "pAspen"],"UCL"] <- xsu0[,cfs0[cfs0 != "pAspen"]]
CFs[!is.na(CFs) & CFs < -10^4] <- -10^4
CFs[!is.na(CFs) & CFs > 10^4] <- 10^4

CFs[,"pAspen",2] <- CFs[,"pAspen",3] <- CFs[,"pAspen",1]

pA2 <- pA[,"pAspen.pa"]
names(pA2) <- rownames(pA)


# north

names(en)
xn <- en$Coef.official
xnl <- en$Coef.official.lci
xnu <- en$Coef.official.uci

compare_sets(colnames(xn0), cfn0)
setdiff(colnames(xn0), cfn0)
setdiff(cfn0, colnames(xn0))

colnames(xn) <- gsub("Spruce", "WhiteSpruce", colnames(xn))
colnames(xn) <- gsub("Decid", "Deciduous", colnames(xn))
colnames(xn) <- gsub("BSpr", "TreedBog", colnames(xn))
colnames(xn) <- gsub("Larch", "TreedFen", colnames(xn))
colnames(xn) <- gsub("Well", "Wellsites", colnames(xn))
colnames(xnl) <- colnames(xnu) <- colnames(xn)

xn0 <- cbind(xn, ShrubbyFen=xn[,"ShrubbyBog"], Water=0, MineV=0)[,cfn0]
xnl0 <- cbind(xnl, ShrubbyFen=xnl[,"ShrubbyBog"], Water=0, MineV=0)[,cfn0]
xnu0 <- cbind(xnu, ShrubbyFen=xnu[,"ShrubbyBog"], Water=0, MineV=0)[,cfn0]

FUN <- log
xn0 <- FUN(xn0)
xnl0 <- FUN(xnl0)
xnu0 <- FUN(xnu0)

xncl <- en$Res.coef.official[rownames(xn0),cfn1]

CFn <- array(NA, c(nrow(xn), length(cfn), 3))
dimnames(CFn) <- list(rownames(xn), cfn, c("Estimate", "LCL", "UCL"))
CFn[,,"Estimate"] <- cbind(xn0, xncl)[,cfn]
CFn[,cfn0[cfn0 != "pAspen"],"LCL"] <- xnl0[,cfn0[cfn0 != "pAspen"]]
CFn[,cfn0[cfn0 != "pAspen"],"UCL"] <- xnu0[,cfn0[cfn0 != "pAspen"]]
CFn[!is.na(CFn) & CFn < -10^4] <- -10^4
CFn[!is.na(CFn) & CFn > 10^4] <- 10^4

WM <- data.frame(Estimate=xn[,"WetlandMargin"],
    xnl[,"WetlandMargin"], xnu[,"WetlandMargin"])[rownames(xn0),]

COEFS2$mammals <- list(north=CFn, south=CFs, species=st, pAspenPA=pA2, WetlandMargin=WM)

save(COEFS2,
    file="s:/AB_data_v2020/Results/COEFS-ALL2.RData")

## mammal results - latest: Dec 16, 2020
library(mefa4)

### organizing all the coefs

ROOT <- "s:/AB_data_v2020/Results/mammal-coefs-2020-12-15"
x0 <- read.csv(file.path(ROOT,
    "Mammal coefficients South June 2020 Best model pApen coefficients.csv"))
x2 <- new.env()
load(file.path(ROOT,
    "Mammal coefficients North June 2020 Best model OFFICIAL coefficients.Rdata"),
    envir=x2)
x3 <- new.env()
load(file.path(ROOT,
    "Mammal coefficients South Oct 2020 Best model OFFICIAL coefficients.Rdata"),
    envir=x3)
names(x0)
x1 <- as.matrix(x0[,-1])
rownames(x1) <- x0$Sp
x1 <- x1[rownames(x3$Coef.official),]

names(x2)
identical(rownames(x2$Coef.official), rownames(x2$Res.coef.official))
identical(dimnames(x2$Coef.official), dimnames(x2$Coef.pa.all))
identical(dimnames(x2$Coef.official), dimnames(x2$Coef.official.lci))
identical(dimnames(x2$Coef.official), dimnames(x2$Coef.official.uci))

names(x3)
identical(rownames(x3$Coef.official), rownames(x1))
identical(rownames(x3$Coef.official), rownames(x3$Res.coef.official))
identical(dimnames(x3$Coef.official), dimnames(x3$Coef.pa.all))
identical(dimnames(x3$Coef.official), dimnames(x3$Coef.official.lci))
identical(dimnames(x3$Coef.official), dimnames(x3$Coef.official.uci))

CF <- list(
    north=list(
        total=x2$Coef.official,
        pa=x2$Coef.pa.all,
        agp=x2$Coef.official/x2$Coef.pa.all,
        lwr=x2$Coef.official.lci,
        upr=x2$Coef.official.uci,
        clim=x2$Res.coef.official
    ),
    south=list(
        total=x3$Coef.official,
        pa=x3$Coef.pa.all,
        agp=x3$Coef.official/x3$Coef.pa.all,
        lwr=x3$Coef.official.lci,
        upr=x3$Coef.official.uci,
        clim=x3$Res.coef.official,
        asp=x1
    )
)
for (i in 1:length(CF$south))
    colnames(CF$south[[i]])[colnames(CF$south[[i]]) == "Well"] <- "Wellsites"
cn <- colnames(CF$north$total)
cn[cn == "Well"] <- "Wellsites"
cn <- gsub("Decid", "Deciduous", cn)
cn <- gsub("Spruce", "WhiteSpruce", cn)
cn <- gsub("BSpr", "TreedBog", cn)
cn <- gsub("Larch", "TreedFen", cn)
for (i in c("total", "pa", "agp", "lwr", "upr")) {
    colnames(CF$north[[i]]) <- cn
    CF$north[[i]] <- CF$north[[i]][,!grepl("9", cn)]
    CF$north[[i]] <- cbind(CF$north[[i]], Water=0, ShrubbyFen=CF$north[[i]][,"ShrubbyBog"])
}



## mapping soil types to soilhf categories
s <- colnames(CF$south$total)
compare_sets(lt$south$Label, s)
setdiff(lt$south$Label, s)
setdiff(s, lt$south$Label)

## mapping veg types to veghf categories
s <- colnames(CF$north$total)
compare_sets(lt$north$Label, s)
setdiff(lt$north$Label, s)
setdiff(s, lt$north$Label)

COEFS3 <- CF

save(COEFS3,
    file="s:/AB_data_v2020/Results/COEFS3-mammals.RData")


## updating habitat element and mammal lookup tables

load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")

uu <- COEFS2$habitats$species

xx <- data.frame(
        SpeciesID=uu$sppid,
        ScientificName=NA,
        TSNID=NA,
        CommonName=uu$CommonName,
        ModelNorth=uu$ModelNorth,
        ModelSouth=uu$ModelSouth,
        UseavailNorth=FALSE,
        UseavailSouth=FALSE,
        Occurrences=NA,
        nSites=NA,
        SizeNorth=NA,
        SizeSouth=NA,
        Nonnative=FALSE,
        LinkHabitat=uu$link,
        LinkSpclim=uu$link,
        AUCNorth=uu$AUCNorth,
        AUCSouth=uu$AUCNorth,
        R2North=NA,
        R2South=NA,
        Comments=uu$Unit,
        Group="habitats")
for (j in 1:ncol(xx))
    if (is.factor(xx[,j]))
        xx[,j] <- as.character(xx[,j])
COEFS2$habitats$species <- xx

uu <- COEFS2$mammals$species

xx <- data.frame(
        SpeciesID=uu$SpeciesID,
        ScientificName=uu$ScientificName,
        TSNID=NA,
        CommonName=uu$CommonName,
        ModelNorth=uu$ModelNorth,
        ModelSouth=uu$ModelSouth,
        UseavailNorth=uu$UseavailNorth,
        UseavailSouth=uu$UseavailSouth,
        Occurrences=NA,
        nSites=NA,
        SizeNorth=NA,
        SizeSouth=NA,
        Nonnative=FALSE,
        LinkHabitat=paste0("N:", uu$LinkHabitatNorth, "-S:", uu$LinkHabitatSouth),
        LinkSpclim=paste0("N:", uu$LinkSpclimNorth, "-S:", uu$LinkSpclimSouth),
        AUCNorth=uu$AUCNorth,
        AUCSouth=uu$AUCNorth,
        R2North=NA,
        R2South=NA,
        Comments=uu$Notes,
        Group="mammals")
for (j in 1:ncol(xx))
    if (is.factor(xx[,j]))
        xx[,j] <- as.character(xx[,j])
COEFS2$mammals$species <- xx

load("s:/AB_data_v2020/Results/COEFS3-mammals.RData")

COEFS2$mammals$fullcf <- COEFS3

save(COEFS2,
    file="s:/AB_data_v2020/Results/COEFS-ALL2.RData")

## adding NN richness

taxon <- "nnplants"

xx <- data.frame(
        SpeciesID="NNplantRich",
        ScientificName=NA,
        TSNID=NA,
        CommonName="Non-native Plant Richness",
        ModelNorth=TRUE,
        ModelSouth=TRUE,
        UseavailNorth=FALSE,
        UseavailSouth=FALSE,
        Occurrences=NA,
        nSites=NA,
        SizeNorth=NA,
        SizeSouth=NA,
        Nonnative=TRUE,
        LinkHabitat="log",
        LinkSpclim="log",
        AUCNorth=NA,
        AUCSouth=NA,
        R2North=NA,
        R2South=NA,
        Comments=NA,
        Group=taxon)
for (j in 1:ncol(xx))
    if (is.factor(xx[,j]))
        xx[,j] <- as.character(xx[,j])

fn <- "s:/AB_data_v2020/Results/Results from Ermias/Non-native vascular plants/R objects North NNplants Coefficient tables 2020_Quad.Rdata"
fs <- "s:/AB_data_v2020/Results/Results from Ermias/Non-native vascular plants/R objects South NNplants Coefficient tables 2020_Quad.Rdata"

en <- new.env()
load(fn, envir=en)
es <- new.env()
load(fs, envir=es)

B <- 100

    #' # North
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    cfn <- en$Coef
    colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
    cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=cfn[,"Grass"])
    cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass")),drop=FALSE]
    xCFn <- cbind(log(cfn), en$Sc.coef)
    #' array: bootstrap x species x coefs
    CFn <- array(0, c(nrow(xCFn), ncol(xCFn), B),
        dimnames=list(rownames(xCFn), colnames(xCFn), NULL))
    for (i in 1:B) {
        cfn <- en$Coef
        colnames(cfn) <- gsub("BlackSpruce", "TreedBog", colnames(cfn))
        cfn <- cbind(cfn,
            TreedFenR=cfn[,"TreedFen"],
            TreedFen1=cfn[,"TreedFen"],
            TreedFen2=cfn[,"TreedFen"],
            TreedFen3=cfn[,"TreedFen"],
            TreedFen4=cfn[,"TreedFen"],
            TreedFen5=cfn[,"TreedFen"],
            TreedFen6=cfn[,"TreedFen"],
            TreedFen7=cfn[,"TreedFen"],
            TreedFen8=cfn[,"TreedFen"],
            Rural=cfn[,"UrbInd"],
            Urban=cfn[,"UrbInd"],
            Industrial=cfn[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=cfn[,"Grass"])
        cfn <- cfn[,!(colnames(cfn) %in% c("TreedFen", "UrbInd", "Grass")),drop=FALSE]
        ses <- en$Coef.se
        colnames(ses) <- gsub("BlackSpruce", "TreedBog", colnames(ses))
        ses <- cbind(ses,
            TreedFenR=ses[,"TreedFen"],
            TreedFen1=ses[,"TreedFen"],
            TreedFen2=ses[,"TreedFen"],
            TreedFen3=ses[,"TreedFen"],
            TreedFen4=ses[,"TreedFen"],
            TreedFen5=ses[,"TreedFen"],
            TreedFen6=ses[,"TreedFen"],
            TreedFen7=ses[,"TreedFen"],
            TreedFen8=ses[,"TreedFen"],
            Rural=ses[,"UrbInd"],
            Urban=ses[,"UrbInd"],
            Industrial=ses[,"UrbInd"],
            Mine=0, MineV=0, Water=0,
            GrassHerb=ses[,"Grass"])
        ses <- ses[,!(colnames(ses) %in% c("TreedFen", "UrbInd", "Grass")),drop=FALSE]
        if (i > 1) {
            se <- ses
            se[se > 5] <- 5
            eps <- rnorm(length(se), 0, se)
        } else {
            eps <- 0
        }
        CFn[,,i] <- cbind(log(cfn) + eps, en$Sc.coef)
    }


    #' # South
    #' Combine together the pieces into a logit scaled matrix (no surrounding effects)
    cfs <- es$Coef
    cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
    cfs <- cfs[,colnames(cfs) != "UrbInd",drop=FALSE]
    xCFs <- cbind(log(cfs), pAspen=es$Coef.pAspen, es$Sc.coef)
    #' array: bootstrap x species x coefs
    CFs <- array(0, c(nrow(xCFs), ncol(xCFs), B),
        dimnames=list(rownames(xCFs), colnames(xCFs), NULL))
    for (i in 1:B) {
        cfs <- es$Coef
        cfs <- cbind(cfs,
            Rural=cfs[,"UrbInd"],
            Urban=cfs[,"UrbInd"],
            Industrial=cfs[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        cfs <- cfs[,colnames(cfs) != "UrbInd",drop=FALSE]
        ses <- es$Coef.se
        ses <- cbind(ses,
            Rural=ses[,"UrbInd"],
            Urban=ses[,"UrbInd"],
            Industrial=ses[,"UrbInd"],
            Mine=0, MineV=0, Water=0)
        ses <- ses[,colnames(ses) != "UrbInd",drop=FALSE]
        if (i > 1) {
            se <- ses
            se[se > 5] <- 5
            eps <- rnorm(length(se), 0, se)
        } else {
            eps <- 0
        }
        CFs[,,i] <- cbind(log(cfs) + eps, pAspen=es$Coef.pAspen, es$Sc.coef)
    }
    #range(CFn)
    #range(CFs)
    CFn[CFn < -10^4] <- -10^4
    CFn[CFn > 10^4] <- 10^4
    CFs[CFs < -10^4] <- -10^4
    CFs[CFs > 10^4] <- 10^4


    COEFS2[[taxon]] <- list(north=CFn, south=CFs, species=xx)

rownames(COEFS2$mammals$species) <- COEFS2$mammals$species$SpeciesID
rownames(COEFS2$habitats$species) <- COEFS2$habitats$species$SpeciesID
rownames(COEFS2$nnplants$species) <- COEFS2$nnplants$species$SpeciesID

save(COEFS2,
    file="s:/AB_data_v2020/Results/COEFS-ALL2.RData")
