rm(list=ls())
od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")

SCALE <- 5
SCALE <- 10
SCALE <- 20

for (SCALE in c(5,10,20)) {
FILE <- sprintf(
    "s:\\AB_data_v2020\\data\\raw\\veghf\\Summary_buf%sm_VegHfSoil_AllYears.csv",
    SCALE)
d <- read.csv(FILE)

UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE


source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]])/(SCALE^2*pi))

save(d_long, d_wide,
    file=sprintf(
        "s:\\AB_data_v2020\\data\\analysis\\veghf\\Summary_buf%sm_VegHfSoil_AllYears.RData",
        SCALE))
}
