#' ---
#' title: "Bird data processing"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Nov 28, 2018"
#' output: pdf_document
#' ---
#'
#' # Preamble
#'
library(mefa4)
library(odbc)
library(DBI)
source("~/.ssh/boreal")
source("~/repos/abmianalytics/birds/00-functions.R")
knitr::opts_chunk$set(eval=FALSE)
#'
#' # Species data
#'
#' Auxiliary data fields:
#'
#' - `PCODE`: project code (methodology is usually constant within `PCODE`)
#' - `SS`: location
#' - `SSYR`: location and year (for year varying predictors)
#' - `PKEY`: survey (visit) to a location, sampling event
#' - `YEAR`, `DATE`, `DATI`: year, date and date/time of the sampling event
#' - `MAXDUR`, `MAXDIS`: maximum duration and distance for sampling event
#' - `CMETHOD`: counting method (device etc)
#'
#' Data subsets:
#'
#' - `*_sm`: ABMI SM units
#' - `*_rf`: ABMI river forks units
#' - `*_bu`: Bioacoustic Unit
#' - `*_bb`: BAM+BBS human PCs
#'
#' ## ABMI RF and SM
#'
ROOT <- "d:/abmi/AB_data_v2018"
e1 <- new.env()
load(file.path(ROOT, "data", "inter", "species", "birds-revisit.Rdata"),
    envir=e1)

resRF <- e1$resRF
resRF$DATE <- strptime(as.character(resRF$FieldDate), "%d-%b-%y")
tmp <- paste0(as.character(resRF$DATE), " ",
    fill_char(as.character(resRF[["StartofPointCount(24hourclock)"]]), "0"), ":00")
tmp[grepl("NA", tmp)] <- NA
resRF$DATI <- as.POSIXlt(tmp)
resRF$MAXDUR <- 10
resRF$MAXDIS <- Inf
resRF$PCODE <- "ABMIRF"
resRF$SS <- paste0(resRF$ABMISite, "_", resRF$subunit)
resRF$SSYR <- resRF$site_year_sub
resRF$PKEY <- resRF$site_year_sub
resRF$YEAR <- resRF$Year
resRF$CMETHOD <- "RF" # River Fork
y_rf <- Xtab(~ PKEY + SpeciesID, resRF, cdrop = c("NONE","SNI", "VNA", "DNC", "PNA"))
d_rf <- droplevels(nonDuplicated(resRF, PKEY, TRUE))

resSM <- e1$resSM
resSM$DATE <- strptime(as.character(resSM$RecordingDate), "%d-%b-%y")
resSM$DATI <- resSM$Start
resSM$MAXDUR <- resSM$Duration
resSM$MAXDIS <- Inf
resSM$PCODE <- "ABMISM"
resSM$SS <- paste0(resSM$ABMISite, "_", resSM$subunit)
resSM$SSYR <- resSM$site_year_sub
tmp <- gsub(" ", "-", as.character(resSM$dati))
resSM$PKEY <- paste0(resSM$site_year_sub, "_", tmp)
resSM$YEAR <- resSM$Year
resRF$CMETHOD <- "SM" # Song Meter
y_sm <- Xtab(~ PKEY + SpeciesID, resSM, cdrop = c("NONE","SNI", "VNA", "DNC", "PNA"))
d_sm <- droplevels(nonDuplicated(resSM, PKEY, TRUE))
#'
#' ## BAM and BBS
#'
e2 <- new.env()
load(file.path("d:/bam/Apr2016/out", "data_package_2016-04-18.Rdata"),
    envir=e2)

TAX <- nonDuplicated(e2$TAX, Species_ID, TRUE)
DAT <- data.frame(e2$PKEY, e2$SS[match(e2$PKEY$SS, e2$SS$SS),])
DAT$SS.1 <- NULL
DAT$PCODE.1 <- NULL
DAT <- droplevels(DAT[!is.na(DAT$JURS) & DAT$JURS == "AB",])
rownames(DAT) <- DAT$PKEY
DAT$CMETHOD <- "HS" # Homo sapiens is the device
DAT$SSYR <- paste0(DAT$SS, "_", DAT$YEAR)
DAT$DATI <- DAT$DATE
DAT$DATE <- as.Date(DAT$DATE)

y_bb <- Xtab(ABUND ~ PKEY + SPECIES_ALL, e2$PCTBL, cdrop="NONE")
d_bb <- droplevels(DAT)
rm(DAT)
#'
#' ## BU
#'
con <- dbConnect(
    odbc::odbc(),
    dsn = "BOREAL",
    database = "EMCLA_SQL_Database",
    driver = "SQL Server",
    uid = ..boreal_db_access$uid,
    pwd = ..boreal_db_access$pwd)
#dbl <- dbListTables(con)
d0 <- dbReadTable(con, "ViewSpecies_LongFormCount")
s0 <- dbReadTable(con, "list Species Code")
dbDisconnect(con)

d0 <- make_char2fact(d0)
s0 <- make_char2fact(s0)

d <- droplevels(d0[d0$Method %in% as.character(c(0, 8, 11:14)) & d0$Replicate == 1,])
d$MAXDUR <- 3
d$MAXDUR[d$Method %in% c("12", "13")] <- 1
d$MAXDUR[d$Method == "0"] <- 10
d$MAXDIS <- Inf
d$SS <- with(d, interaction(
    ProjectID,
    Cluster,
    SITE,
    STATION,
    sep="::", drop=TRUE))
tmp1 <- as.character(d$RECORDING_DATE)
tmp2 <- sapply(strsplit(as.character(d$RECORD_TIME), " "), "[[", 2)
d$YEAR <- as.numeric(substr(tmp1, 1, 4))
d$DATI <- as.POSIXlt(paste(tmp1, tmp2))
d$PKEY <- as.factor(paste0(
    as.character(d$SS),
    "_",
    tmp1,
    "-",
    tmp2))

y_bu <- Xtab(Abundance ~ PKEY + SPECIES, d)
d_bu <- droplevels(nonDuplicated(d, PKEY, TRUE))
#'
#' # Combining species data
#'
if (FALSE) {
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

#d$SppName <- s0$ENGLISH.NAME[match(d$SPECIES, s0$CODE)]
#d$SppeciesID <- d$SppName
#levels(d$SpeciesID) <- nameAlnum(levels(d$SpeciesID), capitalize="mixed", collapse="")


## VEG SOIL HF DATA =============================================

file.path("d:/abmi/AB_data_v2018", "data/analysis/site", "veg-hf_BAM-BBS-BU_v6verified.Rdata")
file.path("d:/abmi/AB_data_v2018", "data/analysis/site", "veg-hf_CameraARU_v6verified.Rdata")

load(file.path("d:/abmi/AB_data_v2018", "data", "analysis", "site",
               "veg-hf_BAM-BBS-BU_v6verified.Rdata"))

mefa4::compare_sets(d$SS, dd_point$SS)

## CLIMATE AND TERRAIN DATA =============================================

file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "cam-aru-bird-2003-2016.csv")
file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "site-center-2003-2016.csv")
file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "1_CamARU2017_v2_Summary_Climate_data.csv")
file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "1_SiteCentre2017_Summary_Climate_data.csv")


}
