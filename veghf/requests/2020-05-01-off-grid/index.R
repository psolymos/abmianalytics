od <- setwd("~/repos/recurring/veghf")

## CamARU off grid sites ---------------------------------

rm(list=ls())
source("00-setup.R")

FILE       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200429_SC_Sites_SummaryTables_Round_2020_batch02.sqlite")
#20200429_OffGrid_Sites_SurveyYear_2019_2020_Buffers_batch02
#20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_ClimateData_batch02
#20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_batch02

## climate raster values

db <- dbConnect(RSQLite::SQLite(), FILE)
#tb <- dbListTables(db)
d <- dbReadTable(db, "20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_ClimateData_batch02")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
summary(clim)
dbDisconnect(db)

## points

TABLE      = "20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_batch02"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = FALSE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
d_long_pt <- d_long

##  large buffers

TABLE      = "20200429_OffGrid_Sites_SurveyYear_2019_2020_Buffers_batch02"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1km <- d_long
d_wide_1km <- d_wide

## small buffers

TABLE      = "20200429_OffGrid_Sites_SurveyYear_2019_2020_Buffers_batch02"
SUB_COL    = "Section"
SUB_VAL    = c("0-56m", "56-100m", "100-150m")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = "s:/AB_data_v2020/data/analysis/site/veg-hf_CAMARU-2019-2020-offgrid_Veg61-vHF.Rdata"
COMMENTS   = "Veg+HF summary for offgrid sites - veg61+vhf // 2020-05-01"
TOL        = 0
SAVE       = c("clim", "d_long_1km", "d_wide_1km", "d_long_pt")
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
source("04-save.R")
