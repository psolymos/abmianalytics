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
resSM$CMETHOD <- "SM" # Song Meter
y_sm <- Xtab(~ PKEY + SpeciesID, resSM, cdrop = c("NONE","SNI", "VNA", "DNC", "PNA"))
d_sm <- droplevels(nonDuplicated(resSM, PKEY, TRUE))
#'
#' ## BAM and BBS
#'
e2 <- new.env()
load(file.path("d:/bam/Apr2016/out", "data_package_2016-04-18.Rdata"),
    envir=e2)

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
colnames(y_bb) <- normalize_sppcode(colnames(y_bb))
d_bb <- droplevels(DAT)
y_bb <- y_bb[rownames(y_bb) %in% rownames(d_bb),]
## Pacific-slope Flycatcher: Western Flycatcher, Cordilleran Flycatcher
y_bb[,"PSFL"] <- y_bb[,"PSFL"] + y_bb[,"WEFL"] + y_bb[,"COFL"]
y_bb[,"WEFL"] <- 0
y_bb[,"COFL"] <- 0
y_bb <- y_bb[,colSums(y_bb)>0]
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
d$PCODE <- "BU"
d$CMETHOD <- "SM"
d$DATE <- d$RECORDING_DATE
d$SSYR <- paste0(d$SS, "_", d$YEAR)
d$SPECIES <- normalize_sppcode(d$SPECIES)

y_bu <- Xtab(Abundance ~ PKEY + SPECIES, d)
d_bu <- droplevels(nonDuplicated(d, PKEY, TRUE))
#'
#' # Combining species data
#'
#' - `code`: AOU code
#' - `species`: common name
#' - `scinam`: scientific name
#' - `family`, `order`: family and order
#' - `sppid`: CamelCase name w/o punctuation
#'
## this is lookup table from previous years
tax0 <- read.csv("~/repos/abmispecies/_data/birds.csv")
tax0 <- tax0[,c("AOU","species","scinam","family","sppid")]
colnames(tax0)[1] <- "code"
tax0$code <- normalize_sppcode(tax0$code)
rownames(tax0) <- tax0$code
## this is BAM lookup table
TAX <- nonDuplicated(e2$TAX, Species_ID, TRUE)
TAX <- TAX[,c("Species_ID", "English_Name", "Scientific_Name", "Family_Sci")]
TAX$sppid <- TAX$English_Name
levels(TAX$sppid) <- nameAlnum(levels(TAX$sppid), capitalize="mixed", collapse="")
colnames(TAX) <- colnames(tax0)
TAX$code <- normalize_sppcode(TAX$code)
rownames(TAX) <- TAX$code
## this is BU lookup table
st <- s0[,c("CODE", "ENGLISH.NAME")]
st$scinam <- paste(s0$Genus, s0$Species)
st$family <- nameAlnum(tolower(s0$Family), "first")
st$sppid <- st$ENGLISH.NAME
levels(st$sppid) <- nameAlnum(levels(st$sppid), capitalize="mixed", collapse="")
colnames(st) <- colnames(tax0)
st$code <- normalize_sppcode(st$code)
st <- nonDuplicated(st, code, TRUE)
## merging
tax <- rbind(tax0, TAX[setdiff(rownames(TAX), rownames(tax0)),])
tax <- rbind(tax, st[setdiff(rownames(st), rownames(tax)),])
tax <- droplevels(tax)

compare_sets(colnames(y_bb), rownames(tax))
compare_sets(colnames(y_bu), rownames(tax))
compare_sets(colnames(y_rf), tax$sppid)
compare_sets(colnames(y_sm), tax$sppid)

setdiff(colnames(y_bb), rownames(tax)) # nothing here
setdiff(colnames(y_bu), rownames(tax)) # probably OK?
setdiff(colnames(y_rf), tax$sppid)
setdiff(colnames(y_sm), tax$sppid)
cn <- c("AmericanGoldenplover" = "AmericanGoldenPlover",
    "EasternWoodpewee" = "EasternWoodPewee",
    "EurasianCollareddove" = "EurasianCollaredDove",
    "GraycrownedRosyfinch" = "GraycrownedRosyFinch",
    "IndianPeafowl" = "CommonPeafowl",
    "MccownsLongspur" = "McCownsLongspur",
    "MyrtlesWarbler" = "YellowrumpedWarbler", # MYWA > YRWA
    "NorthernPygmyowl" = "NorthernPygmyOwl",
    "RossGoose" = "RosssGoose",
    "WesternWoodpewee" = "WesternWoodPewee")
for (i in seq_along(cn)) {
    colnames(y_rf)[colnames(y_rf) == names(cn)[i]] <- cn[i]
    colnames(y_sm)[colnames(y_sm) == names(cn)[i]] <- cn[i]
}

colnames(y_sm) <- rownames(tax)[match(colnames(y_sm), tax$sppid)]
y_sm <- y_sm[,!is.na(colnames(y_sm))]
colnames(y_rf) <- rownames(tax)[match(colnames(y_rf), tax$sppid)]
y_rf <- y_rf[,!is.na(colnames(y_rf))]

SPP <- unique(c(colnames(y_bb[,colSums(y_bb)>0]), colnames(y_bu[,colSums(y_bu)>0]),
    colnames(y_sm[,colSums(y_sm)>0]), colnames(y_rf[,colSums(y_rf)>0])))
SPP <- intersect(SPP, rownames(tax))
SPP <- SPP[!startsWith(SPP, "UN")]
SPP <- SPP[!endsWith(SPP, "_UNI")]
SPP <- SPP[SPP != "NONE"]
SPP <- sort(SPP)
tax <- tax[SPP,]
tax <- tax[!grepl("hybr", tolower(as.character(tax$species))),]
tax <- tax[!grepl("uniden", tolower(as.character(tax$species))),]

tax$order <- as.character(e2$TAX$Order[match(rownames(tax), e2$TAX$Species_ID)])
tax$order2 <- as.character(s0$ORDER[match(rownames(tax), s0$CODE)])
tax$order[is.na(tax$order)] <- tax$order2[is.na(tax$order)]
tax$order <- as.factor(nameAlnum(tolower(tax$order), "first"))
tax <- tax[!(tax$order %in% c("Abiotic", "Anura", "Artiodactyla",
    "Carnivora", "Chiroptera", "Coleoptera", "Lagomorpha", "Perissodactyla", "Rodentia")),]

tax <- tax[!is.na(tax$scinam),]
tax <- tax[!is.na(tax$code),]
SPP <- sort(rownames(tax))
tax <- droplevels(tax[SPP,c("code", "sppid", "species", "scinam", "order", "family")])

#tax$bb <- colSums(y_bb>0)[match(rownames(tax), colnames(y_bb))]
#tax$bu <- colSums(y_bu>0)[match(rownames(tax), colnames(y_bu))]
#tax$rf <- colSums(y_rf>0)[match(rownames(tax), colnames(y_rf))]
#tax$sm <- colSums(y_sm>0)[match(rownames(tax), colnames(y_sm))]
#tax[,c(3,7:10)]

bump <- function(y, ref) {
    m <- Melt(y)
    levels(m$cols) <- c(levels(m$cols), setdiff(ref, levels(m$cols)))
    out <- Xtab(value ~ rows + cols, m)
    out[,ref]
}
y_bb <- bump(y_bb, rownames(tax))
y_bu <- bump(y_bb, rownames(tax))
y_rf <- bump(y_bb, rownames(tax))
y_sm <- bump(y_bb, rownames(tax))

yy <- rbind(y_bb, y_bu, y_rf, y_sm)
cn <- c("PCODE",
    "SS",
    "SSYR",
    "PKEY",
    "YEAR",
    "DATE",
    "DATI",
    "MAXDUR",
    "MAXDIS",
    "CMETHOD")
dd <- make_char2fact(rbind(d_bb[,cn], d_bu[,cn], d_rf[,cn], d_sm[,cn]))
compare_sets(rownames(yy), rownames(dd))
ii <- intersect(rownames(yy), rownames(dd))
yy <- yy[ii,]
dd <- droplevels(dd[ii,])

#'
#' # Veg/soil/HF data
#'

file.path("d:/abmi/AB_data_v2018", "data/analysis/site", "veg-hf_BAM-BBS-BU_v6verified.Rdata")
file.path("d:/abmi/AB_data_v2018", "data/analysis/site", "veg-hf_CameraARU_v6verified.Rdata")

load(file.path("d:/abmi/AB_data_v2018", "data", "analysis", "site",
               "veg-hf_BAM-BBS-BU_v6verified.Rdata"))

mefa4::compare_sets(d$SS, dd_point$SS)

#'
#' # Climate data
#'

file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "cam-aru-bird-2003-2016.csv")
file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "site-center-2003-2016.csv")
file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "1_CamARU2017_v2_Summary_Climate_data.csv")
file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "1_SiteCentre2017_Summary_Climate_data.csv")


