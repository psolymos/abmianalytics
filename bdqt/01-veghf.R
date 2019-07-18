#e:/peter/AB_data_v2018/data/raw/veghf/abmi

HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")




f <- "s:/GC_eric/FromEric/BDQT_2019_veg61RefCond2017_HFI2017/20190705_BDQT_Veg61RefCond2017HFI2017.sqlite"

db <- dbConnect(RSQLite::SQLite(), f)

dbListTables(db)
dbListFields(db, "Veg61HF2017BDQT")
#d <- dbGetQuery(db, "SELECT * FROM `Veg61HF2017BDQT` LIMIT 100000")

d <- dbGetQuery(db,
    "SELECT
      UID, Easting, Northing,
      Origin_Year,
      Pct_of_Larch,
      NSRNAME,
      Soil_Type_1,
      FEATURE_TY,
      Combined_ChgByCWCS,
      YEAR,
      Shape_Area
    FROM
      `Veg61HF2017BDQT`;")
dbDisconnect(db)

d <- make_char2fact(d)
gc()

d2 <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE)

## deal with age 0
levs <- levels(d2$VEGHFAGEclass)[endsWith(levels(d2$VEGHFAGEclass), "0")]
levs <- substr(levs, 1, nchar(levs)-1)
for (i in seq_along(levs)) {
    ii <- paste0(levs[i], c("0", "9"))
    oo <- paste0(levs[i], c("8"))
    d2$VEGAGEclass[d2$VEGAGEclass %in% ii] <- oo
    d2$VEGHFAGEclass[d2$VEGHFAGEclass %in% ii] <- oo
}

data.frame(table(d2$VEGAGEclass))

d2 <- d2[,c("UID", "Easting", "Northing", "NSRNAME",
    "VEGAGEclass", "VEGHFAGEclass",
    "SOILclass", "SOILHFclass")]
gc()


mydb <- dbConnect(RSQLite::SQLite(), "s:/AB_data_v2019/bdqt/bdqt-labels-2017hfi_2019-07-18.sqlite")

dbWriteTable(mydb, "veghf", d2, overwrite = TRUE)

dbDisconnect(mydb)

