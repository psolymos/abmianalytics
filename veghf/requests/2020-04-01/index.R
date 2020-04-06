od <- setwd("~/repos/recurring/veghf")

##  ---------------------------------

FILE       = "S:/AB_data_v2018/data/raw/veghf/site_all/20190129_SummaryTables_CAMARU_2017_2018_Veg61_vHFSPOT2017.sqlite"
TABLE      = "Summary_Buffers"
SUB_COL    = "Section"
SUB_VAL    = c("0-56m", "56-100m", "100-150m") # 0-56m 100-150m 150-564m  56-100m
UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = "s:/AB_data_v2019/data/analysis/site/veg-hf_CAMARU-2017-2018_Veg61-vHF.Rdata"
COMMENTS   = "Veg+HF summary for cam/aru - veg61+vhf2017 // 2020-04-01"
TOL        = 0
SAVE       = NULL
UNROUND    = FALSE

source("00-setup.R")
source("01-data.R")
source("02-long.R")
source("03-wide.R")
source("04-save.R")

