source("~/repos/abmianalytics/veghf/veghf-setup.R")

### Ermias wellsite stuff

## ABMI sites (on+off) cetre 1 ha
f <- "e:/peter/AB_data_v2016/raw_new/wellsite_ermias/BorealSites_sde_20161212.csv"
d <- read.csv(f)
## summarize for Site_Id and Quadrant_I
## sampling year is 2014
with(d, table(Site_ID, Quadrant_I))
d$survey_year <- 2014
d$SiteQ <- with(d, interaction(Site_ID, Quadrant_I, sep="_", drop=TRUE))
head(d)

dd <- make_vegHF_wide(d, col.label = "SiteQ",
    col.year="survey_year", col.HFyear="CutYear")
dd$scale <- "1 ha and 1/4 ha square at well sites"

save(dd, file="e:/peter/AB_data_v2016/raw_new/wellsite_ermias/BorealSites_sde_20161212.Rdata")