HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

ROOT <- "d:/abmi/AB_data_v2019/data/raw/veghf/3x7"

#f <- file.path(ROOT, "20191219_SC_veghf3by7_2017.sqlite")
f <- "s:/GC_eric/FromEric/IC_Reporting_Recovery/20200102_SC_veghf3by7_2017.sqlite"
db <- dbConnect(RSQLite::SQLite(), f)
lt <- dbListTables(db)
lt

## FMA


d0 <- dbReadTable(db, "SC_FMA_Current_ALPAC_Vegetation_HFI2017_in3by7_rawdata")
d0$Shape_Area <- d0$Area_m2
#d1 <- dbReadTable(db, "SC_FMA_Current_ALPAC_Summary_HF_By_Category_unit_hectares")

d <- dbReadTable(db, "SC_FMA_Current_ALPAC_Vegetation_VHF3by7_rawdata")
yr <- sort(unique(d$year_3by7))
d$Shape_Area <- d$Area_m2
str(d)

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

dd <- make_vegHF_wide_v6(d0,
    col.label="ABMI_ID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "Backfill61 + 2017HFI in 3x7km"
dx <- nonDuplicated(d, ABMI_ID, TRUE)[rownames(dd[[1]]),]
dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
dd_hfi2017 <- dd

dbDisconnect(db)

save(list=c("dd_hfi2017", paste0("dd_", yr)),
    file="d:/abmi/AB_data_v2019/data/analysis/alpac/SC_FMA_ALPAC_3by7_veghfSummary_20200102.RData")


## AEI

d0 <- dbReadTable(db, "SC_AEI_Alpac_Mistik_AB_only_Vegetation_HFI2017_in3by7_rawdata")
d0$Shape_Area <- d0$Area_m2
#d1 <- dbReadTable(db, "SC_AEI_Alpac_Mistik_AB_only_Summary_HF_By_Category_unit_hectares")

d <- dbReadTable(db, "SC_AEI_Alpac_Mistik_AB_only_Vegetation_VHF3by7_rawdata")
yr <- sort(unique(d$year_3by7))
d$Shape_Area <- d$Area_m2
str(d)

for (yri in yr) {
    cat(yri)
    flush.console()
    x <- d[d$year_3by7 == yri,]
    issue <- which(x$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & is.na(x$YEAR))
    if (length(issue)) {
        cat("issue", yri)
        x[issue, "YEAR"] <- round(mean(x$YEAR[x$FEATURE_TY %in% c("CUTBLOCK", "HARVEST-AREA") & !is.na(x$YEAR)]))
    }
    cat("\n")
    dd <- make_vegHF_wide_v6(x,
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

dd <- make_vegHF_wide_v6(d0,
    col.label="ABMI_ID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "Backfill61 + 2017HFI in 3x7km"
dx <- nonDuplicated(d, ABMI_ID, TRUE)[rownames(dd[[1]]),]
dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
dd_hfi2017 <- dd

dbDisconnect(db)

save(list=c("dd_hfi2017", paste0("dd_", yr)),
    file="d:/abmi/AB_data_v2019/data/analysis/alpac/SC_AEI_Alpac_Mistik_AB_only_3by7_veghfSummary_20200102.RData")


## Total area

f <- "s:/GC_eric/FromEric/IC_Reporting_Recovery/20200110_SC_veghf2017.sqlite"
db <- dbConnect(RSQLite::SQLite(), f)
(lt <- dbListTables(db))

d1 <- dbReadTable(db, "SC_FMA_Current_ALPAC_Vegetation_HFI2017_rawdata")
d2 <- dbReadTable(db, "SC_AEI_Alpac_Mistik_AB_only_Vegetation_HFI2017_rawdata")

d1$ID <- "ALL"
d2$ID <- "ALL"
dbDisconnect(db)

dd <- make_vegHF_wide_v6(d1,
    col.label="ID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "Backfill61 + 2017HFI in 3x7km"
dd <- fill_in_0ages_v6(dd, "Central Mixedwood", ages_list)
dd_fma <- dd

dd <- make_vegHF_wide_v6(d2,
    col.label="ID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "Backfill61 + 2017HFI in 3x7km"
dd <- fill_in_0ages_v6(dd, "Central Mixedwood", ages_list)
dd_aei <- dd


save(dd_aei, dd_fma,
    file="d:/abmi/AB_data_v2019/data/analysis/alpac/SC_AEIandFMA_veghfSummary_20200110.RData")




d <- read.csv("~/Desktop/Sites_id_VegVhf3by7_2017.csv")

dd <- make_vegHF_wide_v6(d,
    col.label="OBJECTID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    wide=FALSE, HF_fine=TRUE) # use refined classes
str(dd)

d2 <- read.csv("~/Desktop/Sites_id_VegHf_2017.csv")

dd2 <- make_vegHF_wide_v6(d2,
    col.label="OBJECTID",
    col.year=2017,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    wide=FALSE, HF_fine=TRUE) # use refined classes
str(dd2)

write.csv(dd, row.names=FALSE, file="~/Desktop/Sites_id_VegVhf3by7_2017-output.csv")
write.csv(dd2, row.names=FALSE, file="~/Desktop/Sites_id_VegHf_2017-output.csv")
