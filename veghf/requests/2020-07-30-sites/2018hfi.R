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
#load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged.RData")


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

load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged_clean.RData")
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
ds$Origin_Year <- ds$FEATURE_TY <- NULL

save(dd18, file="s:/AB_data_v2020/data/analysis/veghf/w2w_veghf2018_grid1sqkm.RData")


## transitions
od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)

load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged_clean.RData")


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

dd <- data.frame(d10[,c("GRID_LABEL", "NSRNAME", "Shape_Area")])
dd$rf_veg <- d10$VEGAGEclass
dd$rf_age <- d10$AgeRf
dd$rf_soil <- d10$SOILclass
dd$cr_veg_10 <- d10$VEGHFAGEclass
dd$cr_age_10 <- d10$AgeCr
dd$cr_soil_10 <- d10$SOILHFclass
dd$cr_veg_18 <- d18$VEGHFAGEclass
dd$cr_age_18 <- d18$AgeCr
dd$cr_soil_18 <- d18$SOILHFclass

levels(dd$rf_veg) <- gsub("0", "8", levels(dd$rf_veg))
levels(dd$rf_veg) <- gsub("9", "8", levels(dd$rf_veg))

levels(dd$cr_veg_10) <- gsub("0", "8", levels(dd$cr_veg_10))
levels(dd$cr_veg_10) <- gsub("9", "8", levels(dd$cr_veg_10))

levels(dd$cr_veg_18) <- gsub("0", "8", levels(dd$cr_veg_18))
levels(dd$cr_veg_18) <- gsub("9", "8", levels(dd$cr_veg_18))

levels(dd$rf_age) <- gsub("0", "8", levels(dd$rf_age))
levels(dd$rf_age) <- gsub("9", "8", levels(dd$rf_age))

levels(dd$cr_age_10) <- gsub("0", "8", levels(dd$cr_age_10))
levels(dd$cr_age_10) <- gsub("9", "8", levels(dd$cr_age_10))

levels(dd$cr_age_18) <- gsub("0", "8", levels(dd$cr_age_18))
levels(dd$cr_age_18) <- gsub("9", "8", levels(dd$cr_age_18))

rm(ds, d10, d18)
gc()
dd$u <- with(dd, paste(GRID_LABEL, rf_veg, rf_soil, cr_veg_10, cr_veg_18, cr_soil_10, cr_soil_18))
#dd$u <- as.integer(as.factor(dd$u))
table(duplicated(dd$u))

dd$rf_age <- factor(as.character(dd$rf_age), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8"))
dd$cr_age_10 <- factor(as.character(dd$cr_age_10), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8"))
dd$cr_age_18 <- factor(as.character(dd$cr_age_18), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8"))



s <- mefa4::sum_by(dd$Shape_Area, dd$u)
head(s)
dds <- dd[!duplicated(dd$u),]
rownames(dds) <- dds$u

s <- s[rownames(dds),]
dds$Shape_Area <- s[, "x"]

sum(dds$Shape_Area)
sum(dd$Shape_Area)

dd <- dds
dd$u <- NULL
rm(dds, s)
gc()


sum(dd$Shape_Area[dd$rf_age=="0"])/sum(dd$Shape_Area)
sum(dd$Shape_Area[dd$cr_age_10=="0"])/sum(dd$Shape_Area)
sum(dd$Shape_Area[dd$cr_age_18=="0"])/sum(dd$Shape_Area)


save(dd, file="s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged_clean_long.RData")

## transitions

load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged_clean_long.RData")

veg10 <- as.character(dd$cr_veg_10)
veg18 <- as.character(dd$cr_veg_18)

dd$tr <- paste0(veg10, "->", veg18)
dd$tr[veg10 == veg18] <- veg10[veg10 == veg18]

keep <- c("tr", "cr_veg_10", "cr_age_10", "cr_veg_18", "cr_age_18")
dt <- nonDuplicated(dd[,keep], tr, TRUE)
colnames(dt) <- c("tr", "v10", "a10", "v18", "a18")
for (i in 1:ncol(dt))
  dt[[i]] <- as.character(dt[[i]])

trVeg <- Xtab(Shape_Area ~ GRID_LABEL + tr, dd)

rownames(dt) <- dt$tr
dt <- dt[colnames(trVeg),]
dt$area <- colSums(trVeg)

dt <- dt[order(rownames(dt)),]

save(dt, file="s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_transitions.RData")

## cleaning up tr
library(mefa4)

load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_transitions.RData")

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
hf <- rownames(tv)[tv$is_HF]
hfw <- rownames(tv)[tv$is_HF & tv$is_water]
hfcc <- rownames(tv)[tv$is_harvest]
fores <- rownames(tv)[tv$is_forest]
cult <- hf[grep("Cultivation", hf)]

up <- function(dt) {

  dt$changed <- dt$v10 != dt$v18
  dt$washf <- dt$v10 %in% hf
  dt$ishf <- dt$v18 %in% hf
  dt$washw <- dt$v10 %in% hfw
  dt$ishw <- dt$v18 %in% hfw
  dt$wascc <-  dt$v10 %in% hfcc
  dt$iscc <-  dt$v18 %in% hfcc
  dt$wasforest <-  dt$v10 %in% fores
  dt$isforest <-  dt$v18 %in% fores
  dt$wascc2 <- startsWith(dt$v10, "CC")
  dt$iscc2 <- startsWith(dt$v18, "CC")
  dt$wascult <- dt$v10 %in% cult
  dt$iscult <- dt$v10 %in% cult

  dt$an10 <- factor(dt$a10, c("", "R", "1", "2", "3", "4", "5", "6", "7", "8"))
  levels(dt$an10) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140")
  dt$an10 <- as.integer(dt$an10) - 1
  dt$an18 <- factor(dt$a18, c("", "R", "1", "2", "3", "4", "5", "6", "7", "8"))
  levels(dt$an18) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140")
  dt$an18 <- as.integer(dt$an18) - 1
  dt$diff <- dt$an18 - dt$an10

  dt$tr2 <-  ifelse(dt$changed, paste0(dt$v10, "->", dt$v18), dt$v10)

  dt$type <- factor("Other", c("NewFor", "OldFor", "NewHF", "OldHF", "Fire", "Aging0", "Aging1", "Other"))
  dt$type[dt$isforest & dt$wasforest] <- "Aging0"
  dt$type[dt$ishf & !dt$washf] <- "NewHF"
  dt$type[dt$ishf & dt$washf] <- "OldHF"
  dt$type[dt$iscc & !dt$wascc] <- "NewFor"
  dt$type[dt$iscc & dt$wascc] <- "OldFor"
  dt$type[!dt$ishf & !dt$washf & dt$diff < 0] <- "Fire"
  dt$type[!dt$ishf & !dt$washf & dt$diff > 0] <- "Aging1"

  dt

}

dt <- up(dt)
summary(dt)
plot(table(dt$diff))


## HF -> forest with >1 age diff (0.11%)
ss <- dt$washf & !dt$ishf & dt$diff > 1
100*sum(dt$area[ss])/sum(dt$area)
dt$tr2[ss]
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)
plot(table(dt$diff))

## assume 2018 veghf is better when age diff >1 (0.01%)
ss <- dt$diff > 1
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)
plot(table(dt$diff))

## HF water -> HF !water (0.002%)
ss <- dt$washw & dt$ishf & !dt$ishw
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## Non merchantable harvest (0.001%)
ss <- dt$wascc2 & !dt$wascc
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

ss <- dt$iscc2 & !dt$iscc
rownames(dt)[ss]
100*sum(dt$area[ss])/sum(dt$area)

## 2010 pipeline fixed in 2018 (0.34%)
ss <- dt$v10 == "Pipeline"
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)


## forest -> harvest where age is increasing (0.004%)
ss <- !dt$wascc & dt$iscc & dt$diff > 0
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## open -> forest (<0.0001%)
ss <- !dt$wasforest & dt$isforest & dt$diff > 0
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## open -> harvest irrespective of age diff (0.01%)
ss <- !dt$wasforest & dt$iscc
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## forest -> harvest where age is NOT decreasing (0.01%)
table(dt$diff[!dt$wascc & dt$iscc])
ss <- !dt$wascc & dt$iscc & dt$diff == 0
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## young forest (should not be harvested) -> harvest (0.008%)
ss <- !dt$wascc & dt$iscc & dt$an10 < 4
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## recently harvested forest --> harvest where CC is not R (0.03%)
ss <- !dt$wascc & dt$iscc & dt$an18 > 0
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## hatvest -> harvest that got younger (harvested before regen) (0.006%)
ss <- dt$wascc & dt$iscc & dt$diff < 0
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## harvest -> forest (<0.0001%)
ss <- dt$wascc & dt$isforest & !dt$iscc
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## HF water -> !HF (0.001%)
ss <- dt$washw & !dt$ishf
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## !CC HF -> !HF (0.008%)
ss <- dt$washf & !dt$ishf
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)

## !CC forest -> !CC forest and age does not go to R (<0.00001%)
ss <- dt$wasforest & !dt$wascc & dt$isforest & !dt$iscc & dt$diff < 0 & dt$a18 != "R"
dt$tr2[ss]
100*sum(dt$area[ss])/sum(dt$area)
dt$v10[ss] <- dt$v18[ss]
dt$a10[ss] <- dt$a18[ss]
dt <- up(dt)



## overall changes: lots of labels, 0.314% of area
100*sum(dt$tr != dt$tr2)/nrow(dt)
100*sum(dt$area[dt$tr != dt$tr2])/sum(dt$area)
table(dt$type)

save(dt, file="s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_transitions.RData")


s <- mefa4::sum_by(dt$area, dt$tr2)
dt2 <- mefa4::nonDuplicated(dt, tr2, TRUE)
s <- s[rownames(dt2),]
dt2$area <- s[, "x"]

dt2 <- dt2[order(dt2$tr2),]
dt2 <- dt2[order(dt2$type),]
dt2 <- dt2[,c("v10", "a10",  "an10", "v18", "a18", "an18",  "area", "diff", "tr2", "type")]

write.csv(dt2, row.names = FALSE, file="~/GoogleWork/abmi/tr1018/transitions-w2w-2010-2018.csv")



## use clean transitions to process dd -----------------------------------------------------



od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)

checks <- function(xvr, xvc) {
  i <- intersect(colnames(xvc), colnames(xvr))
  d <- (xvr[,i] - xvc[,i]) / 10^6
  list(rf_only=setdiff(colnames(xvr), colnames(xvc)),
    cr_only=setdiff(colnames(xvc), colnames(xvr)),
    rf_minus_cr_range=range(d))
}

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
hf <- rownames(tv)[tv$is_HF]
hf <- c(hf, "CutBlocks")

load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_grid1sqkm_merged_clean_long.RData")
load("s:/AB_data_v2020/data/inter/veghf/veghf2018_2010_transitions.RData")

veg10 <- as.character(dd$cr_veg_10)
veg18 <- as.character(dd$cr_veg_18)

dd$rf_age <- dd$cr_age_10 <- dd$cr_age_18 <- NULL
dd$NSRNAME <- NULL
for (i in c("rf_veg", "rf_soil", "cr_veg_10", "cr_soil_10", "cr_veg_18", "cr_soil_18"))
  dd[[i]] <- as.character(dd[[i]])

dd$tr <- paste0(veg10, "->", veg18)
dd$tr[veg10 == veg18] <- veg10[veg10 == veg18]

dd$tr2 <- dt$tr2[match(dd$tr, dt$tr)]
dd$cr_veg_10 <- dt$v10[match(dd$tr, dt$tr)]
dd$cr_veg_18 <- dt$v18[match(dd$tr, dt$tr)]

ss <- dd$cr_soil_10 %in% hf & !(dd$cr_veg_10 %in% hf)
100*sum(dd$Shape_Area[ss])/sum(dd$Shape_Area)
dd$cr_soil_10[ss] <- dd$rf_soil[ss]

ss <- dd$cr_soil_18 %in% hf & !(dd$cr_veg_18 %in% hf)
100*sum(dd$Shape_Area[ss])/sum(dd$Shape_Area)

## year specific reference
dd$rf_veg_10 <- dd$cr_veg_10
dd$rf_veg_10[dd$cr_veg_10 %in% hf] <- dd$rf_veg[dd$cr_veg_10 %in% hf]
dd$rf_veg_18 <- dd$cr_veg_18
dd$rf_veg_18[dd$cr_veg_18 %in% hf] <- dd$rf_veg[dd$cr_veg_18 %in% hf]
dd$rf_soil_10 <- dd$cr_soil_10
dd$rf_soil_10[dd$cr_soil_10 %in% hf] <- dd$rf_soil[dd$cr_soil_10 %in% hf]
dd$rf_soil_18 <- dd$cr_soil_18
dd$rf_soil_18[dd$cr_soil_18 %in% hf] <- dd$rf_soil[dd$cr_soil_18 %in% hf]

dd$tr <- dd$tr2 <- NULL
dd$rf_veg <- dd$rf_soil <- NULL

with(dd, table(isHF=cr_veg_10 %in% hf, cr_rf_same=cr_veg_10 == rf_veg_10))
with(dd, table(isHF=cr_veg_18 %in% hf, cr_rf_same=cr_veg_18 == rf_veg_18))
with(dd, table(isHF=cr_soil_10 %in% hf, cr_rf_same=cr_soil_10 == rf_soil_10))
with(dd, table(isHF=cr_soil_18 %in% hf, cr_rf_same=cr_soil_18 == rf_soil_18))


## 2010 and 2018 w2w
dd_2010 <- list(
  veg_current=Xtab(Shape_Area ~ GRID_LABEL + cr_veg_10, dd),
  veg_reference=Xtab(Shape_Area ~ GRID_LABEL + rf_veg_10, dd),
  soil_current=Xtab(Shape_Area ~ GRID_LABEL + cr_soil_10, dd),
  soil_reference=Xtab(Shape_Area ~ GRID_LABEL + rf_soil_10, dd))
dd_2018 <- list(
  veg_current=Xtab(Shape_Area ~ GRID_LABEL + cr_veg_18, dd),
  veg_reference=Xtab(Shape_Area ~ GRID_LABEL + rf_veg_18, dd),
  soil_current=Xtab(Shape_Area ~ GRID_LABEL + cr_soil_18, dd),
  soil_reference=Xtab(Shape_Area ~ GRID_LABEL + rf_soil_18, dd))

checks(dd_2010$veg_reference, dd_2010$veg_current)
checks(dd_2010$soil_reference, dd_2010$soil_current)

checks(dd_2018$veg_reference, dd_2018$veg_current)
checks(dd_2018$soil_reference, dd_2018$soil_current)

## ref->2018 transitions
dd$trv_1018 <- ifelse(dd$cr_veg_10 == dd$cr_veg_18, dd$cr_veg_10, paste0(dd$cr_veg_10, "->", dd$cr_veg_18))
dd$trv_rf18 <- ifelse(dd$rf_veg_18 == dd$cr_veg_18, dd$rf_veg_18, paste0(dd$rf_veg_18, "->", dd$cr_veg_18))
dd$trs_rf18 <- ifelse(dd$rf_soil_18 == dd$cr_soil_18, dd$rf_soil_18, paste0(dd$rf_soil_18, "->", dd$cr_soil_18))

trVeg_1018 <- Xtab(Shape_Area ~ GRID_LABEL + trv_1018, dd)
trVeg <- Xtab(Shape_Area ~ GRID_LABEL + trv_rf18, dd)
trSoil <- Xtab(Shape_Area ~ GRID_LABEL + trs_rf18, dd)

lts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v2020.csv")
rownames(lts) <- lts[,1]
colnames(lts)[1] <- "ID"
lts <- lts[colnames(dd_2018$soil_current),]
ltv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v2020.csv")
rownames(ltv) <- ltv[,1]
colnames(ltv)[1] <- "ID"
ltv <- ltv[colnames(dd_2018$veg_current),]

cn <- colnames(trVeg_1018)
tmp <- strsplit(cn, "->")
chVeg_1018 <- data.frame(
  label=cn,
  cr10=sapply(tmp, "[[", 1),
  cr18=sapply(tmp, function(z) if (length(z) > 1) z[2] else z[1]))
chVeg_1018$sector10 <- ltv$Sector[match(chVeg_1018$cr10, ltv$ID)]
chVeg_1018$sector18 <- ltv$Sector[match(chVeg_1018$cr18, ltv$ID)]
chVeg_1018$type <- dt$type[match(chVeg_1018$label, dt$tr2)]
rownames(chVeg_1018) <- cn

cn <- colnames(trVeg)
tmp <- strsplit(cn, "->")
chVeg <- data.frame(
  label=cn,
  rf=sapply(tmp, "[[", 1),
  cr=sapply(tmp, function(z) if (length(z) > 1) z[2] else z[1]))
chVeg$sector <- ltv$Sector[match(chVeg$cr, ltv$ID)]
rownames(chVeg) <- cn

cn <- colnames(trSoil)
tmp <- strsplit(cn, "->")
chSoil <- data.frame(
  label=cn,
  rf=sapply(tmp, "[[", 1),
  cr=sapply(tmp, function(z) if (length(z) > 1) z[2] else z[1]))
chSoil$sector <- lts$Sector[match(chSoil$cr, lts$ID)]
rownames(chSoil) <- cn

save(
  ltv, lts,
  dd_2010,
  file="s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_2010_wide.RData"
)
save(
  ltv, lts,
  dd_2018,
  file="s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_2018_wide.RData"
)

save(
  trVeg,  trSoil,
  chVeg, chSoil,
  file="s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData"
)

save(trVeg_1018, chVeg_1018,
  file="s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_2010_2018_transitions_wide.RData"
)


## checking HFI compared to 2016


library(mefa4)
library(raster)

load("s:/AB_data_v2018/data/analysis/grid/veg-hf_grid_v61hf2016v3WildFireUpTo2016.Rdata")
load("s:/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
load("s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_2010_wide.RData")
load("s:/AB_data_v2020/data/analysis/veghf/veghf_w2w_2018_wide.RData")


rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

make_raster <- function(value, rc, rt)
{
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}

cv <- rownames(ltv)
cs <- rownames(lts)

dd_2016 <- dd_kgrid
a9 <- c("Decid9", "Mixedwood9", "Pine9", "Spruce9", "TreedBog9", "TreedFen9", "TreedSwamp9")
a8 <- c("Decid8", "Mixedwood8", "Pine8", "Spruce8", "TreedBog8", "TreedFen8", "TreedSwamp8")
cn1 <- colnames(dd_kgrid[[1]])
names(cn1) <- cn1
cn1[a9] <- a8
dd_2016[[1]] <- groupSums(dd_kgrid[[1]], 2, cn1)[,cv]

cn2 <- colnames(dd_kgrid[[2]])
names(cn2) <- cn2
cn2[a9] <- a8
dd_2016[[2]] <- groupSums(dd_kgrid[[2]], 2, cn2)[,colnames(dd_2010[[2]])]

compare_sets(cv, colnames(dd_2016[[1]]))
compare_sets(cv, colnames(dd_2010[[1]]))
compare_sets(cv, colnames(dd_2018[[1]]))

compare_sets(cs, colnames(dd_2016[[3]]))
compare_sets(cs, colnames(dd_2010[[3]]))
compare_sets(cs, colnames(dd_2018[[3]]))


v2010 <- dd_2010[[1]][,cv]
v2016 <- dd_2016[[1]][rownames(v2010),cv]
v2018 <- dd_2018[[1]][rownames(v2010),cv]
v2010 <- v2010/rowSums(v2010)
v2016 <- v2016/rowSums(v2016)
v2018 <- v2018/rowSums(v2018)

s2010 <- dd_2010[[3]][,cs]
s2016 <- dd_2016[[3]][rownames(s2010),cs]
s2018 <- dd_2018[[3]][rownames(s2010),cs]
s2010 <- s2010/rowSums(s2010)
s2016 <- s2016/rowSums(s2016)
s2018 <- s2018/rowSums(s2018)

dv10 <- v2016-v2010
dv18 <- v2018-v2016
ds10 <- s2016-s2010
ds18 <- s2018-s2016

dv <- data.frame(
  #max10=apply(abs(dv10), 2, max),
  avg10=apply(abs(dv10), 2, mean),
  #max18=apply(abs(dv18), 2, max),
  avg18=apply(abs(dv18), 2, mean))

ds <- data.frame(
  #max10=apply(abs(dv10), 2, max),
  avg10=apply(abs(ds10), 2, mean),
  #max18=apply(abs(dv18), 2, max),
  avg18=apply(abs(ds18), 2, mean))

round(data.frame(avg=colMeans(abs(v2018-v2010))), 4)
round(data.frame(avg=colMeans(abs(s2018-s2010))), 4)

v2010x <- groupSums(v2010, 2, gsub("[:0-9:]", "", cv))
v2016x <- groupSums(v2016, 2, gsub("[:0-9:]", "", cv))
v2018x <- groupSums(v2018, 2, gsub("[:0-9:]", "", cv))

ddv <- data.frame(
  tv10 = colSums(v2010)/sum(v2010),
  tv16 = colSums(v2016)/sum(v2016),
  tv18 = colSums(v2018)/sum(v2018))
ddvx <- data.frame(
  tv10 = colSums(v2010x)/sum(v2010x),
  tv16 = colSums(v2016x)/sum(v2016x),
  tv18 = colSums(v2018x)/sum(v2018x))

dds <- data.frame(
  ts10 = colSums(s2010)/sum(s2010),
  ts16 = colSums(s2016)/sum(s2016),
  ts18 = colSums(s2018)/sum(s2018))

plot(ts16 ~ ts10, dds);abline(0,1)
plot(ts18 ~ ts10, dds);abline(0,1)
plot(ts18 ~ ts16, dds);abline(0,1)

plot(tv16 ~ tv10, ddv);abline(0,1)
plot(tv18 ~ tv10, ddv);abline(0,1)
plot(tv18 ~ tv16, ddv);abline(0,1)

plot(tv16 ~ tv10, ddvx);abline(0,1)
plot(tv18 ~ tv10, ddvx);abline(0,1)
plot(tv18 ~ tv16, ddvx);abline(0,1)



