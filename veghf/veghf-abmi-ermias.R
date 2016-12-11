source("~/repos/abmianalytics/veghf/veghf-setup.R")

### Ermias wellsite stuff

## ABMI sites (on+off) cetre 1 ha
f <- "e:/peter/AB_data_v2016/raw_new/wellsite_ermias/Boreal_sde_sites_jc_20161208.csv"
d <- read.csv(f)
## summarize for Site_Id and Quadrant_I
## sampling year is 2014
d$survey_year <- 2014
d$Site_YEAR <- with(d, interaction(ABMI_Assigned_Site_ID, survey_year, sep="_", drop=TRUE))
head(d)

dd <- make_vegHF_wide(d, col.label = "Site_YEAR",
    col.year="survey_year", col.HFyear="year_")
dd$scale <- "1 ha square around site centre"

