od <- setwd("~/repos/recurring/veghf")

##  ---------------------------------

FILE       = "d:/abmi/AB_data_v2019/data/raw/veghf/20200318_SC_veghf2017.sqlite"
TABLE      = "SC_LandscapeUnit_Vegetation_HFI2017_rawdata"
SUB_COL    = NULL
SUB_VAL    = ""
UID_COL    = "Region"
VEG_COL    = "Vegetation"
BASE_YR    = 2017
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = NULL
COMMENTS   = "Veg+HF summary for OSM LU's - veg+hf2017 // 2020-03-18"
TOL        = 0
SAVE       = NULL
UNROUND    = FALSE

source("00-setup.R")
source("01-data.R")
source("02-long.R")
source("03-wide.R")
source("04-save.R")

