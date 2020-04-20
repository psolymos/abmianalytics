od <- setwd("~/repos/recurring/veghf")

##  1ha site centres---------------------------------

rm(list=ls())
source("00-setup.R")

FILE       = "s:/GC_eric/FromEric/Sites_summaries/Round2020/20200224_SC_Sites_SummaryTables_Round_2020.sqlite"
db <- dbConnect(RSQLite::SQLite(), FILE)
#tb <- dbListTables(db)
d <- dbReadTable(db, "20200224_TERRESTRIAL_SurveyYear_2019_Points_ClimateData_batch01")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
dbDisconnect(db)


##  large buffer around site centres

TABLE      = "20200224_TERRESTRIAL_SurveyYear_2019_Buffers_batch01"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
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

## 1ha level
TABLE      = "20200224_TERRESTRIAL_SurveyYear_2019_Buffers_batch01"
SUB_COL    = "Section"
SUB_VAL    = c("NE", "NW", "SE", "SW")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = "s:/AB_data_v2020/data/analysis/site/veg-hf_SITE1HA-2019_Veg61-vHF.Rdata"
COMMENTS   = "Veg+HF summary for site centre 1ha - veg61+vhf // 2020-04-20"
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1ha <- d_long
d_wide_1ha <- d_wide

## quadrant level
TABLE      = "20200224_TERRESTRIAL_SurveyYear_2019_Buffers_batch01"
SUB_COL    = "Section"
SUB_VAL    = c("NE", "NW", "SE", "SW")
UID_COL    = "UID" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = "s:/AB_data_v2020/data/analysis/site/veg-hf_SITE1HA-2019_Veg61-vHF.Rdata"
COMMENTS   = "Veg+HF summary for site centre 1ha - veg61+vhf // 2020-04-20"
TOL        = 0
SAVE       = c("clim", "d_long_1km", "d_wide_1km",  "d_long_1ha", "d_wide_1ha")
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))

source("04-save.R")



##  150m Cam/ARU buffers ---------------------------------

rm(list=ls())
source("00-setup.R")

FILE       = "s:/GC_eric/FromEric/Sites_summaries/Round2020/20200224_SC_Sites_SummaryTables_Round_2020.sqlite"
db <- dbConnect(RSQLite::SQLite(), FILE)
#tb <- dbListTables(db)
d <- dbReadTable(db, "20200224_CAMARU_SurveyYear_2019_Points_ClimateData_batch01")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
dbDisconnect(db)

## points

TABLE      = "20200224_CAMARU_SurveyYear_2019_Points_batch01"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = FALSE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
d_long_pt <- d_long

##  large buffer around site centres

TABLE      = "20200224_CAMARU_SurveyYear_2019_Buffers_batch01"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
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


TABLE      = "20200224_CAMARU_SurveyYear_2019_Buffers_batch01"
SUB_COL    = "Section"
SUB_VAL    = c("0-56m", "56-100m", "100-150m")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
OUTPUT     = "s:/AB_data_v2020/data/analysis/site/veg-hf_CAMARU-2019_Veg61-vHF.Rdata"
COMMENTS   = "Veg+HF summary for site centre 150m buffer - veg61+vhf // 2020-04-20"
TOL        = 0
SAVE       = c("clim", "d_long_1km", "d_wide_1km", "d_long_pt")
UNROUND    = FALSE

source("01-data.R")
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
source("04-save.R")
