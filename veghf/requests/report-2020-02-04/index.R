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

od <- setwd("~/repos/recurring/veghf")
source("function.R")
setwd(od)

## 20200130_SC_veghf3by7_2017_NWSAR -------------------------

FILE       = "s:/GC_eric/FromEric/IC_Reporting_Recovery/NWSAR/20200130_SC_veghf3by7_2017_NWSAR.sqlite"
TABLE      = "SC_NWSAR_project_Summary_Vegetation_VHF3by7_rawdata"
SUB_COL    = NULL
SUB_VAL    = ""
UID_COL    = "Region"
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Area_m2"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "NWSAR Reporting - veg+hf veg+vhf3by7 // 2020-01-30"
TOL        = 0

od <- setwd("~/repos/recurring/veghf")
source("function.R")
setwd(od)

## 20200131_SC_veghf2017_AB_LUF_NR_OSR_Caribou_7Gen ----------------------

od <- setwd("~/repos/recurring/veghf")

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
setwd(od)

## 20200131_SC_veghf3by7_2017_AB_LUF_NR_OSR_Caribou_7Gen ----------------------

od <- setwd("~/repos/recurring/veghf")

FILE       = "s:/GC_eric/FromEric/IC_Reporting_Recovery/AB_LUF_NR_OSR_Caribou_7Gen/20200131_SC_veghf3by7_2017_AB_LUF_NR_OSR_Caribou_7Gen.sqlite"
SUB_COL    = NULL
SUB_VAL    = ""
UID_COL    = "Region"
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Area_m2"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "AB_LUF_NR_OSR_Caribou_7Gen 3x7 // 2020-01-31"
TOL        = 0.001

TABLES <- c("SC_Caribou_Range_Vegetation_VHF3by7_rawdata",
  "SC_LAND_USE_FRAMEWORK_Vegetation_VHF3by7_rawdata", # CCDecid0 CCPine0
  "SC_NaturalSubRegion_Vegetation_VHF3by7_rawdata", # CCDecid0 CCPine0
  "SC_Oilsand3Region_Vegetation_VHF3by7_rawdata", # done
  "SC_Oilsand_Mineable_Vegetation_VHF3by7_rawdata", # done
  "SC_SevenGeneration_Vegetation_VHF3by7_rawdata") # CCPine0
for (i in TABLES) {

  TABLE      = i
  OUTPUT     = paste0("s:/GC_eric/FromEric/IC_Reporting_Recovery/AB_LUF_NR_OSR_Caribou_7Gen/", TABLE, ".RData")

  source("function.R")
}
setwd(od)

h <- c("HARVEST-AREA", "HARVEST_AREA", "CUTBLOCK")
j <- d_long$FEATURE_TY %in% h & is.na(d_long$YEAR)
sum(j)
sum(d_long$Shape_Area[j])/10^4 # ha
