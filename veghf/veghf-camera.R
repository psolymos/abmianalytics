HF_VERSION <- "2014_coarse"

source("~/repos/abmianalytics/veghf/veghf-setup.R")

d <- read.csv(file.path(ROOT, VER, "data", "raw", "veghf", "site",
    "cam_aru_4birds_on_W2WBkfVegHFV6.csv"))
d <- droplevels(d[d$Year != 0,])
d$Label <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, Year, sep="::", drop=TRUE))
dd <- make_vegHF_wide_v6(d, col.label="Label",
    col.year="Year", col.HFyear="YEAR_MineCFOAgCutblock", wide=FALSE,
    HF_fine=HF_VERSION=="2014")

write.csv(dd, row.names=FALSE, file.path(ROOT, VER, "data", "inter", "veghf", "site",
    "cam_aru_W2WBkfVegHFV6-2015-2016.csv"))


## 2016 code

source("~/repos/abmianalytics/veghf/veghf-setup.R")
recl <- read.csv("c:/Users/Peter/Dropbox/abmi/V6/veg-V6-combined3.csv")

fn <- "camARU_pts_vs_bkfV6.csv"
f <- file.path(ROOT, VER, "data", "kgrid-V6", fn)
d <- read.csv(f)
d <- c4_fun(d)
d$Year_[is.na(d$Year_)] <- 2017
d$Site_YEAR <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, Year_, sep="::", drop=TRUE))
save(d, file=file.path(ROOT, VER, "data", "kgrid-V6", "camARU_pts_vs_bkfV6"))

fn <- "camARU_buf150m_vs_bkfV6.csv"
f <- file.path(ROOT, VER, "data", "kgrid-V6", fn)
d <- read.csv(f)
d <- c4_fun(d)
d$Year_[is.na(d$Year_)] <- 2017
d$Site_YEAR <- with(d, interaction(Site_ID, deployment, Cam_ARU_Bird_Location, Year_, sep="::", drop=TRUE))

dd <- make_vegHF_wide_v6(d, col.label = "Site_YEAR",
    col.year="Year_", col.HFyear="CutYear")
dd$scale <- "150 m radius circle around site centre"

save(dd, file=file.path(ROOT, VER, "data", "kgrid-V6", "camARU_buf150m_vs_bkfV6"))


