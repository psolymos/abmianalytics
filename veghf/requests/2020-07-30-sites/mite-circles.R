rm(list=ls())
od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")

FILE <- "s:\\AB_data_v2020\\data\\raw\\veghf\\Summary_buf5m_VegHfSoil_AllYears.csv"
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
summary(rowSums(d_wide[[1]])/(5^2*pi))

save(d_long, d_wide,
    file="s:\\AB_data_v2020\\data\\analysis\\veghf\\Summary_buf5m_VegHfSoil_AllYears.RData")
