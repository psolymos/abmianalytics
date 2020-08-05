od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)


f <- "s:/GC_eric/FromEric/Veghf/20200618_SC_veghf2018_2010_grid1skqm.sqlite"


db <- dbConnect(RSQLite::SQLite(), f)

dbListTables(db)
dbListFields(db, "veghf2018_2010_grid1sqkm")
di <- dbGetQuery(db, "SELECT * FROM `veghf2018_2010_grid1sqkm` LIMIT 1000;")

d <- dbGetQuery(db,
    "SELECT
      GRID_LABEL,
      Origin_Year,
      Origin_Year_cor,
      Origin_Year_2010cond,
      Origin_Year_cor_2010,
      Pct_of_Larch,
      NSRNAME,
      Soil_Type_1,
      FEATURE_TY,
      FEATURE_TY_2010cond,
      Combined_ChgByCWCS,
      YEAR,
      YEAR_2010cond,
      SHAPE_Area
    FROM
      `veghf2018_2010_grid1sqkm`;")


dbDisconnect(db)

d <- make_char2fact(d)
save(d, file="s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm.RData")


cn <- colnames(d)[colnames(d) != "SHAPE_Area"]
tmp <- NULL
for (i in cn) {
    tmp <- paste0(tmp, "_", d[[i]])
}
table(duplicated(tmp))
tmp <- as.factor(tmp)
tmp <- as.integer(tmp)
d$uid <- tmp
rm(tmp)
gc()

s <- mefa4::sum_by(d$SHAPE_Area, d$uid)
head(s)
ds <- d[!duplicated(d$uid),]
rownames(ds) <- ds$uid

rm(d)
gc()

s <- s[rownames(ds),]

ds$SHAPE_Area <- s[, "x"]
ds$n_poly <- s[, "by"]

rm(s)
gc()

save(ds, file="s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged.RData")

ds$oyear18 <- ds$Origin_Year
ds$oyear18[ds$Origin_Year_cor != 9999] <- ds$Origin_Year[ds$Origin_Year_cor != 9999] +
    ds$Origin_Year_cor[ds$Origin_Year_cor != 9999]

ds$oyear10 <- ds$Origin_Year_2010cond
ds$oyear10[ds$Origin_Year_cor_2010 != 9999] <- ds$Origin_Year_2010cond[ds$Origin_Year_cor_2010 != 9999] +
    ds$Origin_Year_cor_2010[ds$Origin_Year_cor_2010 != 9999]

ds$Origin_Year_2018cond <- ds$Origin_Year
ds$Origin_Year_cor_2018 <- ds$Origin_Year_cor
ds$Origin_Year_cor <- NULL

ds$FEATURE_TY_2018cond <- ds$FEATURE_TY

ds$Shape_Area <- ds$SHAPE_Area

ds$FEATURE_TY_2018cond[which(ds$FEATURE_TY_2018cond == "WELL-UNKNOWN")] <- "WELL-OTHER"

ds$Origin_Year <- NULL
ds$FEATURE_TY <- NULL
ds$uid <- NULL
ds$n_poly <- NULL
ds$FEATURE_TY <- NULL
ds$SHAPE_Area <- NULL

ds$Origin_Year_2010cond <- NULL
ds$Origin_Year_2018cond <- NULL
ds$Origin_Year_cor_2010 <- NULL
ds$Origin_Year_cor_2018 <- NULL
ds$YEAR_2018cond <- ds$YEAR
ds$YEAR <- NULL
ds <- ds[,sort(colnames(ds))]

save(ds, file="s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged_clean.RData")

## 2018 w2w

ds$Origin_Year <- ds$oyear18
ds$FEATURE_TY <- ds$FEATURE_TY_2018cond
dd18 <- make_vegHF_wide_v6(ds,
    col.label="GRID_LABEL",
    col.year=2018,
    col.HFyear="YEAR_2018cond",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)

x <- nonDuplicated(ds, GRID_LABEL, TRUE)[rownames(dd18)[[1]],]
dd18 <- fill_in_0ages_v6(dd18, x$NSRNAME, ages_list)

save(dd18, file="s:/AB_data_v2020/data/analysis/veghf/w2w_veghf2018_grid1sqkm.RData")


## transitions


ds$Origin_Year <- ds$oyear18
ds$FEATURE_TY <- ds$FEATURE_TY_2018cond
d18 <- make_vegHF_wide_v6(ds,
    col.label="GRID_LABEL",
    col.year=2018,
    col.HFyear="YEAR_2018cond",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE)


ds$Origin_Year <- ds$oyear10
ds$FEATURE_TY <- ds$FEATURE_TY_2010cond
d10 <- make_vegHF_wide_v6(ds,
    col.label="GRID_LABEL",
    col.year=2010,
    col.HFyear="YEAR_2010cond",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=FALSE)




