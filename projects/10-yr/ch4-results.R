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
    comments=NA
))
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

## lin N
linn1 <- data.frame(Species=veg[,"Species"], linn,
    veg[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
colnames(linn1) <- c("Species", "AverageCoef", "SoftLin10", "HardLin10", "SoftLin",
    "SoftLin.LCL", "SoftLin.UCL", "HardLin", "HardLin.LCL", "HardLin.UCL")

## soil, nontreed and treed
soilnt1 <- soil
soilnt1$X <- NULL
colnames(soilnt1) <- gsub(".LCI", ".LCL", colnames(soilnt1))
colnames(soilnt1) <- gsub(".UCI", ".UCL", colnames(soilnt1))
soilnt1 <- soilnt1[,colnames(soilnt0)]
soiltr1 <- data.frame(Species=soilnt1$Species,
    plogis(qlogis(as.matrix(soilnt1[,-1])) + paspen$pAspen))

## lin S
lins1 <- data.frame(Species=soil[,"Species"], lins,
    soil[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
colnames(lins1) <- c("Species", "AverageCoef", "SoftLin10", "HardLin10", "SoftLin",
    "SoftLin.LCL", "SoftLin.UCL", "HardLin", "HardLin.LCL", "HardLin.UCL")

## use-avail
usen1 <- usen
usen1$Species <- usen1$SpLabel
usen1$SpLabel <- NULL
usen1 <- usen1[,!grepl("\\.WRSI", colnames(usen1))]
colnames(usen1) <- gsub("\\.rWRSI", "", colnames(usen1))

uses1 <- uses
uses1$Species <- uses1$SpLabel
uses1$SpLabel <- NULL
uses1 <- uses1[,!grepl("\\.WRSI", colnames(uses1))]
colnames(uses1) <- gsub("\\.rWRSI", "", colnames(uses1))

## sector
sectn1 <- sectn
sectn1 <- data.frame(Species=sectn1$Sp, sectn1)
sectn1$Sp <- sectn1$pcTotalProvAbundance <- sectn1$Reg_Model <- NULL

## note: this does not have Forestry (Grassland!)
sects1 <- sects
sects1 <- data.frame(Species=sects1$Sp, sects1)
sects1$Sp <- sects1$pcTotalProvAbundance <- sects1$Reg_Model <- NULL
sects1$pcTotalEffect.Forestry <- sects1$pcUnitEffect.Forestry <- 0
sects1 <- sects1[,colnames(sectn1)]

AllIn[[TAX]] <- list(lt=lt1, usen=usen1, uses=uses1,
    veg=veg1, linn=linn1, soilnt=soilnt1, soiltr=soiltr1,
    lins=lins1, sectn=sectn1, sects=sects1)
}


