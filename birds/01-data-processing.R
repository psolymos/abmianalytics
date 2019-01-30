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
#' - `VKEY`: key to match with veg/HF info
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
resRF$YEAR <- resRF$Year
resRF$MAXDUR <- 10
resRF$MAXDIS <- Inf
resRF$PCODE <- "ABMIRF"
resRF$SS <- paste0(resRF$ABMISite, "_", resRF$subunit)
resRF$SSYR <- paste0(resRF$ABMISite, "_", resRF$subunit, "_", resRF$YEAR)
resRF$PKEY <- resRF$SSYR
resRF$CMETHOD <- "RF" # River Fork
y_rf <- Xtab(~ PKEY + SpeciesID, resRF, cdrop = c("NONE","SNI", "VNA", "DNC", "PNA"))
d_rf <- make_char2fact(droplevels(nonDuplicated(resRF, PKEY, TRUE)))
d_rf$VKEY <- d_rf$PKEY

resSM <- e1$resSM
resSM$DATE <- strptime(as.character(resSM$RecordingDate), "%d-%b-%y")
resSM$DATI <- resSM$Start
resSM$YEAR <- resSM$Year
resSM$MAXDUR <- resSM$Duration
resSM$MAXDIS <- Inf
resSM$PCODE <- "ABMISM"
resSM$SS <- paste0(resSM$ABMISite, "_", resSM$subunit)
resSM$SSYR <- paste0(resSM$ABMISite, "_", resSM$subunit, "_", resSM$YEAR)
tmp <- gsub(" ", "-", as.character(resSM$DATI))
resSM$PKEY <- paste0(resSM$SSYR, "_", tmp)
resSM$CMETHOD <- "SM" # Song Meter
y_sm <- Xtab(~ PKEY + SpeciesID, resSM, cdrop = c("NONE","SNI", "VNA", "DNC", "PNA"))
d_sm <- make_char2fact(droplevels(nonDuplicated(resSM, PKEY, TRUE)))
d_sm$VKEY <- d_sm$SSYR
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
d_bb <- make_char2fact(droplevels(DAT))
d_bb$VKEY <- d_bb$PKEY
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
xybu <- dbReadTable(con, "viewCoord_LatLong")
d0 <- dbReadTable(con, "ViewSpecies_LongFormCount")
s0 <- dbReadTable(con, "list Species Code")

dbDisconnect(con)

d0 <- make_char2fact(d0)
s0 <- make_char2fact(s0)
xybu <- make_char2fact(xybu)

xybu$SS <- with(xybu, interaction(
    ProjectID,
    Cluster,
    SITE,
    STATION,
    sep="::", drop=TRUE))

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
d$PCODE <- paste0("BU_", d$ProjectID)
d$CMETHOD <- "SM"
d$DATE <- d$RECORDING_DATE
d$SSYR <- paste0(d$SS, "_", d$YEAR)
d$SPECIES <- normalize_sppcode(d$SPECIES)
d$X <- xybu$Longitude[match(d$SS, xybu$SS)]
d$Y <- xybu$Latitude[match(d$SS, xybu$SS)]

y_bu <- Xtab(Abundance ~ PKEY + SPECIES, d)
d_bu <- make_char2fact(droplevels(nonDuplicated(d, PKEY, TRUE)))
d_bu$VKEY <- d_bu$SS
#'
#' # Veg/soil/HF data
#'
#' - `vs0`: veg+soil current+reference at point level (data frame)
#' - `vc1`: veg current 150m buffer (matrix)
#' - `vc2`: veg current 564m buffer (matrix)
#' - `vr1`: veg reference 150m buffer (matrix)
#' - `vr2`: veg reference 564m buffer (matrix)
#' - `sc1`: soil current 150m buffer (matrix)
#' - `sc2`: soil current 564m buffer (matrix)
#' - `sr1`: soil reference 150m buffer (matrix)
#' - `sr2`: soil reference 564m buffer (matrix)
#'
## checking BAM+BBS+BU
e3 <- new.env()
load(file.path("d:/abmi/AB_data_v2018", "data/analysis/site",
    "veg-hf_BAM-BBS-BU_v6verified.Rdata"), envir=e3)

head(rownames(e3$dd_150m[[1]]))
head(d_bb$PKEY)

compare_sets(e3$dd_point$UID, d_bb$PKEY)
compare_sets(rownames(e3$dd_150m[[1]]), d_bb$PKEY)
head(setdiff(rownames(e3$dd_150m[[1]]), d_bb$PKEY))

head(rownames(e3$dd_150m[[1]]))
head(d_bu$SS)

compare_sets(rownames(e3$dd_150m[[1]]), d_bu$SS)
compare_sets(rownames(e3$dd_150m[[1]]), c(as.character(d_bb$PKEY), as.character(d_bu$SS)))
## BAM+BBS
vs0_bb <- e3$dd_point[match(d_bb$VKEY, e3$dd_point$UID),]
vc1_bb <- as.matrix(e3$dd_150m$veg_current)[match(d_bb$VKEY, rownames(e3$dd_150m$veg_current)),]
vr1_bb <- as.matrix(e3$dd_150m$veg_reference)[match(d_bb$VKEY, rownames(e3$dd_150m$veg_reference)),]
sc1_bb <- as.matrix(e3$dd_150m$soil_current)[match(d_bb$VKEY, rownames(e3$dd_150m$soil_current)),]
sr1_bb <- as.matrix(e3$dd_150m$soil_reference)[match(d_bb$VKEY, rownames(e3$dd_150m$soil_reference)),]
vc2_bb <- as.matrix(e3$dd_564m$veg_current)[match(d_bb$VKEY, rownames(e3$dd_564m$veg_current)),]
vr2_bb <- as.matrix(e3$dd_564m$veg_reference)[match(d_bb$VKEY, rownames(e3$dd_564m$veg_reference)),]
sc2_bb <- as.matrix(e3$dd_564m$soil_current)[match(d_bb$VKEY, rownames(e3$dd_564m$soil_current)),]
sr2_bb <- as.matrix(e3$dd_564m$soil_reference)[match(d_bb$VKEY, rownames(e3$dd_564m$soil_reference)),]
rownames(vs0_bb) <-
    rownames(vc1_bb) <- rownames(vc2_bb) <-
    rownames(vr1_bb) <- rownames(vr2_bb) <-
    rownames(sc1_bb) <- rownames(sc2_bb) <-
    rownames(sr1_bb) <- rownames(sr2_bb) <- rownames(d_bb)
## BU
vs0_bu <- e3$dd_point[match(d_bu$VKEY, e3$dd_point$UID),]
vc1_bu <- as.matrix(e3$dd_150m$veg_current)[match(d_bu$VKEY, rownames(e3$dd_150m$veg_current)),]
vr1_bu <- as.matrix(e3$dd_150m$veg_reference)[match(d_bu$VKEY, rownames(e3$dd_150m$veg_reference)),]
sc1_bu <- as.matrix(e3$dd_150m$soil_current)[match(d_bu$VKEY, rownames(e3$dd_150m$soil_current)),]
sr1_bu <- as.matrix(e3$dd_150m$soil_reference)[match(d_bu$VKEY, rownames(e3$dd_150m$soil_reference)),]
vc2_bu <- as.matrix(e3$dd_564m$veg_current)[match(d_bu$VKEY, rownames(e3$dd_564m$veg_current)),]
vr2_bu <- as.matrix(e3$dd_564m$veg_reference)[match(d_bu$VKEY, rownames(e3$dd_564m$veg_reference)),]
sc2_bu <- as.matrix(e3$dd_564m$soil_current)[match(d_bu$VKEY, rownames(e3$dd_564m$soil_current)),]
sr2_bu <- as.matrix(e3$dd_564m$soil_reference)[match(d_bu$VKEY, rownames(e3$dd_564m$soil_reference)),]
rownames(vs0_bu) <-
    rownames(vc1_bu) <- rownames(vc2_bu) <-
    rownames(vr1_bu) <- rownames(vr2_bu) <-
    rownames(sc1_bu) <- rownames(sc2_bu) <-
    rownames(sr1_bu) <- rownames(sr2_bu) <- rownames(d_bu)
## checking ABMI RF+SM
e4 <- new.env()
load(file.path("d:/abmi/AB_data_v2018", "data/analysis/site", "veg-hf_CameraARU_v6verified.Rdata"),
    envir=e4)

head(rownames(e4$dd_150m[[1]]))
head(d_rf$PKEY)
head(d_sm$SSYR)

compare_sets(rownames(e4$dd_150m[[1]]), d_rf$PKEY)
compare_sets(setdiff(rownames(e4$dd_150m[[1]]), d_rf$PKEY), d_sm$SSYR)
compare_sets(rownames(e4$dd_150m[[1]]), c(as.character(d_sm$SSYR), as.character(d_rf$PKEY)))

as.character(setdiff(rownames(e4$dd_150m[[1]]), c(as.character(d_sm$SSYR), as.character(d_rf$PKEY))))
as.character(setdiff(c(as.character(d_sm$SSYR), as.character(d_rf$PKEY)), rownames(e4$dd_150m[[1]])))
## ABMI RF
vs0_rf <- e4$dd_point[match(d_rf$VKEY, e4$dd_point$Site_bird_year),]
vc1_rf <- as.matrix(e4$dd_150m$veg_current)[match(d_rf$VKEY, rownames(e4$dd_150m$veg_current)),]
vr1_rf <- as.matrix(e4$dd_150m$veg_reference)[match(d_rf$VKEY, rownames(e4$dd_150m$veg_reference)),]
sc1_rf <- as.matrix(e4$dd_150m$soil_current)[match(d_rf$VKEY, rownames(e4$dd_150m$soil_current)),]
sr1_rf <- as.matrix(e4$dd_150m$soil_reference)[match(d_rf$VKEY, rownames(e4$dd_150m$soil_reference)),]
vc2_rf <- as.matrix(e4$dd_564m$veg_current)[match(d_rf$VKEY, rownames(e4$dd_564m$veg_current)),]
vr2_rf <- as.matrix(e4$dd_564m$veg_reference)[match(d_rf$VKEY, rownames(e4$dd_564m$veg_reference)),]
sc2_rf <- as.matrix(e4$dd_564m$soil_current)[match(d_rf$VKEY, rownames(e4$dd_564m$soil_current)),]
sr2_rf <- as.matrix(e4$dd_564m$soil_reference)[match(d_rf$VKEY, rownames(e4$dd_564m$soil_reference)),]
rownames(vs0_rf) <-
    rownames(vc1_rf) <- rownames(vc2_rf) <-
    rownames(vr1_rf) <- rownames(vr2_rf) <-
    rownames(sc1_rf) <- rownames(sc2_rf) <-
    rownames(sr1_rf) <- rownames(sr2_rf) <- rownames(d_rf)
## ABMI SM
vs0_sm <- e4$dd_point[match(d_sm$VKEY, e4$dd_point$Site_bird_year),]
vc1_sm <- as.matrix(e4$dd_150m$veg_current)[match(d_sm$VKEY, rownames(e4$dd_150m$veg_current)),]
vr1_sm <- as.matrix(e4$dd_150m$veg_reference)[match(d_sm$VKEY, rownames(e4$dd_150m$veg_reference)),]
sc1_sm <- as.matrix(e4$dd_150m$soil_current)[match(d_sm$VKEY, rownames(e4$dd_150m$soil_current)),]
sr1_sm <- as.matrix(e4$dd_150m$soil_reference)[match(d_sm$VKEY, rownames(e4$dd_150m$soil_reference)),]
vc2_sm <- as.matrix(e4$dd_564m$veg_current)[match(d_sm$VKEY, rownames(e4$dd_564m$veg_current)),]
vr2_sm <- as.matrix(e4$dd_564m$veg_reference)[match(d_sm$VKEY, rownames(e4$dd_564m$veg_reference)),]
sc2_sm <- as.matrix(e4$dd_564m$soil_current)[match(d_sm$VKEY, rownames(e4$dd_564m$soil_current)),]
sr2_sm <- as.matrix(e4$dd_564m$soil_reference)[match(d_sm$VKEY, rownames(e4$dd_564m$soil_reference)),]
rownames(vs0_sm) <-
    rownames(vc1_sm) <- rownames(vc2_sm) <-
    rownames(vr1_sm) <- rownames(vr2_sm) <-
    rownames(sc1_sm) <- rownames(sc2_sm) <-
    rownames(sr1_sm) <- rownames(sr2_sm) <- rownames(d_sm)
#'
#' # Climate data
#'
## ABMI RF+SM
pat <- c(paste0("-", c("NW", "NE", "SW", "SE"), "-BOTH"),
    paste0("-", c("NW", "NE", "SW", "SE"), "-CAM"),
    paste0("-", c("NW", "NE", "SW", "SE"), "-ARU"), "-BOTH")
c1 <- read.csv(file.path("d:/abmi/AB_data_v2018", "data/raw/clim", "cam-aru-bird-2003-2016.csv"))
c1$SITE <- c1$Site_ID
a <- levels(c1$SITE)
a[endsWith(a, "B")] <- substr(a[endsWith(a, "B")], 1, nchar(a[endsWith(a, "B")])-1)
a[endsWith(a, "-")] <- substr(a[endsWith(a, "-")], 1, nchar(a[endsWith(a, "-")])-1)
a <- Gsub("-b", "", a)
a <- Gsub(pat, "", a)
levels(c1$SITE) <- a
c1$SS <- as.factor(paste0(c1$SITE, "_", c1$Cam_ARU_Bird_Location))
c1 <- nonDuplicated(c1, SS, TRUE)

c2 <- read.csv(file.path("d:/abmi/AB_data_v2018", "data/raw/clim",
    "1_CamARU2017_v2_Summary_Climate_data.csv"))
c2$pAspen <- c2$Populus_tremuloides_brtpred_nofp
c2$PET <- c2$Eref
levels(c2$Cam_ARU_Bird_Location) <- gsub("-b", "", levels(c2$Cam_ARU_Bird_Location))
c2$SITE <- c2$Site_ID
a <- levels(c2$SITE)
a[endsWith(a, "B")] <- substr(a[endsWith(a, "B")], 1, nchar(a[endsWith(a, "B")])-1)
a[endsWith(a, "-")] <- substr(a[endsWith(a, "-")], 1, nchar(a[endsWith(a, "-")])-1)
a <- Gsub("-b", "", a)
a <- Gsub(pat, "", a)
levels(c2$SITE) <- a
c2 <- droplevels(c2[!endsWith(as.character(c2$SITE), ")"),])
c2$SS <- as.factor(paste0(c2$SITE, "_", c2$Cam_ARU_Bird_Location))
c2 <- nonDuplicated(c2, SS, TRUE)

d2 <- e4$dd_point
d2$SS <- as.factor(paste0(
    as.character(d2$SITE), "_",
    ifelse(is.na(d2$Cam_ARU_Bird_Location), "NA", as.character(d2$Cam_ARU_Bird_Location))))

compare_sets(d2$SS, c2$SS)
setdiff(c2$SS, d2$SS)

cn <- c("SS", "SITE", "AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen")
ca <- rbind(c1[,cn], c2[,cn])
compare_sets(c(as.character(d_rf$SS),as.character(d_sm$SS)), ca$SS)

ca$nearest <- as.numeric(as.character(ca$SITE))
tmp <- as.character(ca$SITE)[is.na(ca$nearest)]
out <- tmp
out[] <- NA
stmp <- strsplit(tmp, "-")
out[startsWith(tmp, "OGW-")] <- sapply(stmp[startsWith(tmp, "OGW-")], "[[", 3)
out[startsWith(tmp, "OG-")] <- sapply(stmp[startsWith(tmp, "OG-")], "[[", 3)
ca$nearest[is.na(ca$nearest)] <- as.integer(out)
xy0 <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
ca$X <- xy0$PUBLIC_LONGITUDE[match(ca$nearest, xy0$SITE_ID)]
ca$Y <- xy0$PUBLIC_LATTITUDE[match(ca$nearest, xy0$SITE_ID)]
## BBS+BAM+BU
library(sp)
library(raster)
library(rgdal)

xy1 <- rbind(nonDuplicated(d_bb[,c("X", "Y")], d_bb$SS, TRUE),
    nonDuplicated(d_bu[,c("X", "Y")], d_bu$SS, TRUE))
xy1 <- xy1[!is.na(xy1[,1]),]
coordinates(xy1) <- ~ X + Y
proj4string(xy1) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

z <- c("AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "Eref",
    "Populus_tremuloides_brtpred_nofp")
r <- list()
for (i in z)
    r[[i]] <- raster(paste0("d:/spatial/ab-climate/", i, ".asc"))
r <- stack(r)
proj4string(r) <- proj4string(cure4insect::.read_raster_template())

tmp <- extract(r, xy1)
cb <- data.frame(coordinates(xy1), tmp)
cb$SS <- rownames(cb)
cb$PET <- cb$Eref
cb$pAspen <- cb$Populus_tremuloides_brtpred_nofp
## Regions
d_bb$NRNAME <- e3$dd_point$NRNAME[match(d_bb$VKEY, e3$dd_point$UID)]
d_bb$NSRNAME <- e3$dd_point$NSRNAME[match(d_bb$VKEY, e3$dd_point$UID)]
d_bb$LUF_NAME <- e3$dd_point$LUF_NAME[match(d_bb$VKEY, e3$dd_point$UID)]

d_bu$NRNAME <- e3$dd_point$NRNAME[match(d_bu$VKEY, e3$dd_point$UID)]
d_bu$NSRNAME <- e3$dd_point$NSRNAME[match(d_bu$VKEY, e3$dd_point$UID)]
d_bu$LUF_NAME <- e3$dd_point$LUF_NAME[match(d_bu$VKEY, e3$dd_point$UID)]

d_rf$NRNAME <- e4$dd_point$NRNAME[match(d_rf$VKEY, e4$dd_point$Site_bird_year)]
d_rf$NSRNAME <- e4$dd_point$NSRNAME[match(d_rf$VKEY, e4$dd_point$Site_bird_year)]
d_rf$LUF_NAME <- e4$dd_point$LUF_NAME[match(d_rf$VKEY, e4$dd_point$Site_bird_year)]

d_sm$NRNAME <- e4$dd_point$NRNAME[match(d_sm$VKEY, e4$dd_point$Site_bird_year)]
d_sm$NSRNAME <- e4$dd_point$NSRNAME[match(d_sm$VKEY, e4$dd_point$Site_bird_year)]
d_sm$LUF_NAME <- e4$dd_point$LUF_NAME[match(d_sm$VKEY, e4$dd_point$Site_bird_year)]
## adding climate and xy
cn <- c("AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen", "X", "Y")
d_rf <- data.frame(d_rf, ca[match(d_rf$SS, ca$SS), cn])
d_sm <- data.frame(d_sm, ca[match(d_sm$SS, ca$SS), cn])
d_bb <- data.frame(d_bb, cb[match(d_bb$SS, cb$SS), cn])
d_bu <- data.frame(d_bu, cb[match(d_bu$SS, cb$SS), cn])
#'
#' # Combining species data and other tables
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
    "Carnivora", "Chiroptera", "Coleoptera", "Lagomorpha",
    "Perissodactyla", "Rodentia")),]

tax <- tax[!is.na(tax$scinam),]
tax <- tax[!is.na(tax$code),]
SPP <- sort(rownames(tax))
tax <- droplevels(tax[SPP,c("code", "sppid", "species", "scinam",
    "order", "family")])

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
yy <- rbind(
    bump(y_bb, rownames(tax)),
    bump(y_bu, rownames(tax)),
    bump(y_rf, rownames(tax)),
    bump(y_sm, rownames(tax)))
cn <- c("PCODE",
    "SS",
    "SSYR",
    "PKEY",
    "YEAR",
    "DATE",
    "DATI",
    "MAXDUR",
    "MAXDIS",
    "CMETHOD",
    "ROAD",
    "AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", "PET", "pAspen",
    "X", "Y", "NRNAME", "NSRNAME", "LUF_NAME")
d_bu$ROAD <- NA
d_sm$ROAD <- NA
d_rf$ROAD <- NA
dd <- rbind(d_bb[,cn], d_bu[,cn], d_rf[,cn], d_sm[,cn])
compare_sets(rownames(yy), rownames(dd))
ii <- intersect(rownames(yy), rownames(dd))
yy <- yy[ii,]
dd <- droplevels(dd[ii,])
cn <- c("HFclass", "VEGclass", "AgeRf", "AgeCr",
    "VEGAGEclass", "VEGHFclass", "VEGHFAGEclass", "SOILclass", "SOILHFclass")
vs0 <- rbind(vs0_bb[,cn], vs0_bu[,cn], vs0_rf[,cn], vs0_sm[,cn])[ii,]
## don't store as sparse matrix with NAs (slow and big)
vc1 <- rbind(vc1_bb, vc1_bu, vc1_rf, vc1_sm)[ii,]
vr1 <- rbind(vr1_bb, vr1_bu, vr1_rf, vr1_sm)[ii,]
sc1 <- rbind(sc1_bb, sc1_bu, sc1_rf, sc1_sm)[ii,]
sr1 <- rbind(sr1_bb, sr1_bu, sr1_rf, sr1_sm)[ii,]
vc2 <- rbind(vc2_bb, vc2_bu, vc2_rf, vc2_sm)[ii,]
vr2 <- rbind(vr2_bb, vr2_bu, vr2_rf, vr2_sm)[ii,]
sc2 <- rbind(sc2_bb, sc2_bu, sc2_rf, sc2_sm)[ii,]
sr2 <- rbind(sr2_bb, sr2_bu, sr2_rf, sr2_sm)[ii,]
#'
#' # Save
#'
save(tax, yy, dd, vs0, vc1, vr1, sc1, sr1, vc2, vr2, sc2, sr2,
    file="d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData")
#' Save pieces for EMB
if (FALSE) {
    out <- data.frame(dd, as.matrix(yy), as.matrix(vc2))
    write.csv(out, row.names=FALSE, file="bird-data-AB-all-2019-01-22.csv")
}
#' The End
