od <- setwd("~/repos/recurring/veghf")

## 20200130_SC_veghf2017_NWSAR ---------------------------------

FILE       = "s:/GC_eric/FromEric/IC_Reporting_Recovery/NWSAR/20200130_SC_veghf2017_NWSAR.sqlite"
TABLE      = "SC_NWSAR_project_Summary_Vegetation_HFI_rawdata"
SUB_COL    = NULL
SUB_VAL    = ""
UID_COL    = "Region"
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "NWSAR Reporting - veg+hf veg+vhf3by7 // 2020-01-30"
TOL        = 0
SAVE       = NULL

source("function.R")

## 20200130_SC_veghf3by7_2017_NWSAR -------------------------

FILE       = "s:/GC_eric/FromEric/IC_Reporting_Recovery/NWSAR/20200130_SC_veghf3by7_2017_NWSAR.sqlite"
TABLE      = "SC_NWSAR_project_Summary_Vegetation_VHF3by7_rawdata"
SUB_COL    = NULL
SUB_VAL    = ""
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Area_m2"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "NWSAR Reporting - veg+hf veg+vhf3by7 // 2020-01-30"
TOL        = 0
SAVE       = "lookup"

source("00-setup.R")
source("01-data.R")

table(d$year_3by7)
d$UID <- paste0(d$Region, "_", d$ABMI_ID, "_", d$year_3by7)
UID_COL <- "UID"

source("02-long.R")
source("03-wide.R")

#summary(rowSums(d_wide[[1]])/(21*10^6))
#sum(rowSums(d_wide[[1]])/(21*10^6) < 0.9999)/nrow(d_wide[[1]])
tmp <- strsplit(rownames(d_wide[[1]]), "_")
lookup <- data.frame(
  region=sapply(tmp, "[[", 1),
  siteid=as.integer(sapply(tmp, "[[", 2)),
  year=as.integer(sapply(tmp, "[[", 3)))
rownames(lookup) <- rownames(d_wide[[1]])

source("04-save.R")


## 20200131_SC_veghf2017_AB_LUF_NR_OSR_Caribou_7Gen ----------------------


FILE       = "s:/GC_eric/FromEric/IC_Reporting_Recovery/AB_LUF_NR_OSR_Caribou_7Gen/20200131_SC_veghf2017_AB_LUF_NR_OSR_Caribou_7Gen.sqlite"
SUB_COL    = NULL
SUB_VAL    = ""
UID_COL    = "Region"
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "AB_LUF_NR_OSR_Caribou_7Gen // 2020-01-31"
TOL        = 0
SAVE       = NULL

TABLES <- c("SC_Caribou_Range_Vegetation_HFI2017_rawdata",
  "SC_LAND_USE_FRAMEWORK_Vegetation_HFI2017_rawdata",
  "SC_NaturalSubRegion_Vegetation_HFI2017_rawdata",
  "SC_Oilsand3Region_Vegetation_HFI2017_rawdata",
  "SC_Oilsand_Mineable_Vegetation_HFI2017_rawdata",
  "SC_SevenGeneration_Vegetation_HFI2017_rawdata")
for (i in TABLES) {

  TABLE      = i
  OUTPUT     = paste0("s:/GC_eric/FromEric/IC_Reporting_Recovery/AB_LUF_NR_OSR_Caribou_7Gen/", TABLE, ".RData")

  source("function.R")
}

## 20200131_SC_veghf3by7_2017_AB_LUF_NR_OSR_Caribou_7Gen ----------------------

FILE       = "s:/GC_eric/FromEric/IC_Reporting_Recovery/AB_LUF_NR_OSR_Caribou_7Gen/20200131_SC_veghf3by7_2017_AB_LUF_NR_OSR_Caribou_7Gen.sqlite"
SUB_COL    = NULL
SUB_VAL    = ""
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Area_m2"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "AB_LUF_NR_OSR_Caribou_7Gen 3x7 // 2020-01-31"
TOL        = 0.001
SAVE       = "lookup"

source("00-setup.R")

TABLES <- c(
  "SC_Caribou_Range_Vegetation_VHF3by7_rawdata",
  "SC_LAND_USE_FRAMEWORK_Vegetation_VHF3by7_rawdata",
  "SC_NaturalSubRegion_Vegetation_VHF3by7_rawdata",
  "SC_Oilsand3Region_Vegetation_VHF3by7_rawdata",
  "SC_Oilsand_Mineable_Vegetation_VHF3by7_rawdata",
  "SC_SevenGeneration_Vegetation_VHF3by7_rawdata")
for (i in 1:length(TABLES)) {

  cat("--------------------", i, "---------------------\n")
  TABLE <- TABLES[i]
  OUTPUT <- paste0("s:/GC_eric/FromEric/IC_Reporting_Recovery/AB_LUF_NR_OSR_Caribou_7Gen/", TABLE,
    "_", Sys.Date(), ".RData")

  source("01-data.R")

  d$UID <- paste0(d$Region, "_", d$ABMI_ID, "_", d$year_3by7)
  UID_COL <- "UID"

  source("02-long.R")
  source("03-wide.R")

  tmp <- strsplit(rownames(d_wide[[1]]), "_")
  lookup <- data.frame(
    region=sapply(tmp, "[[", 1),
    siteid=as.integer(sapply(tmp, "[[", 2)),
    year=as.integer(sapply(tmp, "[[", 3)))
  rownames(lookup) <- rownames(d_wide[[1]])

  source("04-save.R")
}

setwd(od)
