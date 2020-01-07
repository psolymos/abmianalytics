HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

ROOT <- "d:/abmi/AB_data_v2019/data/raw/veghf/3x7"

f <- file.path(ROOT, "20191219_SC_veghf3by7_2017.sqlite")
db <- dbConnect(RSQLite::SQLite(), f)
lt <- dbListTables(db)

d <- dbReadTable(db, "SC_FMA_Current_ALPAC_Vegetation_VHF3by7_rawdata")
yr <- sort(unique(d$year_3by7))
d$Shape_Area <- d$Area_m2
str(d)

yr <- as.integer(substr(lt, nchar(lt)-3, nchar(lt)))

for (yri in yr) {
    cat(yri, "\n")
    flush.console()
    dd <- make_vegHF_wide_v6(d[d$year_3by7 == yri,],
        col.label="ABMI_ID",
        col.year=yri,
        col.HFyear="YEAR",
        col.HABIT="Vegetation",
        col.SOIL="Soil_Type_1",
        sparse=TRUE, HF_fine=TRUE) # use refined classes
    dd$scale <- "Backfill61 + VerifiedHF in 3x7km"
    dx <- nonDuplicated(d, ABMI_ID, TRUE)[rownames(dd[[1]]),]
    dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
    assign(paste0("dd_", yri), dd)
}


save(list=paste0("dd_", yr),
    file="d:/abmi/AB_data_v2019/data/analysis/alpac/SC_FMA_Current_ALPAC_Vegetation_VHF3by7_veghfSummary.Rdata")
