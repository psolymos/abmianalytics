
# b1-3: batches
FILE1       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200224_SC_Sites_SummaryTables_Round_2020.sqlite")
FILE2       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200429_SC_Sites_SummaryTables_Round_2020_batch02.sqlite")
FILE3       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200506_SC_Sites_SummaryTables_Round_2020_batch03.sqlite")

db <- dbConnect(RSQLite::SQLite(), FILE1)
tt <- dbListTables(db)
for (i in tt) {
    cat("\n---------", i, "\n")
    print(head(dbGetQuery(db, paste0("SELECT * FROM '", i, "' LIMIT 5"))$UID))
}
dbDisconnect(db)

db <- dbConnect(RSQLite::SQLite(), FILE2)
tt <- dbListTables(db)
for (i in tt) {
    cat("\n---------", i, "\n")
    print(head(dbGetQuery(db, paste0("SELECT * FROM '", i, "' LIMIT 5"))$UID))
}
dbDisconnect(db)

db <- dbConnect(RSQLite::SQLite(), FILE3)
tt <- dbListTables(db)
for (i in tt) {
    cat("\n---------", i, "\n")
    print(head(dbGetQuery(db, paste0("SELECT * FROM '", i, "' LIMIT 5"))$UID))
}
dbDisconnect(db)


##  1ha site centres from batch 1 ---------------------------------

rm(list=ls())
od <- setwd("~/repos/recurring/veghf")
FILE1       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200224_SC_Sites_SummaryTables_Round_2020.sqlite")

source("00-setup.R")

FILE       = FILE1
db <- dbConnect(RSQLite::SQLite(), FILE)
(tb <- dbListTables(db))

d <- dbReadTable(db, "20200224_TERRESTRIAL_SurveyYear_2019_Points_ClimateData_batch01")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
dbDisconnect(db)

clim <- clim

##  large buffer around site centres

TABLE      = "20200224_TERRESTRIAL_SurveyYear_2019_Buffers_batch01"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1km <- d_long
d_wide_1km <- d_wide

## 1ha level
TABLE      = "20200224_TERRESTRIAL_SurveyYear_2019_Buffers_batch01"
SUB_COL    = "Section"
SUB_VAL    = c("NE", "NW", "SE", "SW")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1ha <- d_long
d_wide_1ha <- d_wide

## quadrant level
TABLE      = "20200224_TERRESTRIAL_SurveyYear_2019_Buffers_batch01"
SUB_COL    = "Section"
SUB_VAL    = c("NE", "NW", "SE", "SW")
UID_COL    = "UID" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))


SAVE       = c("clim", "d_long_1km", "d_wide_1km",  "d_long_1ha", "d_wide_1ha")
OUTPUT     = "s:/AB_data_v2020/data/analysis/site/veg-hf_SITE1HA-2019_Veg61-vHF.Rdata"
COMMENTS   = "Veg+HF summary for site centre 1ha - veg61+vhf // 2020-04-20"

source("04-save.R")



## CamARU off grid sites from batches 2-3 ---------------------------------

rm(list=ls())
setwd("~/repos/recurring/veghf")
FILE2       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200429_SC_Sites_SummaryTables_Round_2020_batch02.sqlite")
FILE3       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200506_SC_Sites_SummaryTables_Round_2020_batch03.sqlite")

source("00-setup.R")

## climate raster values

FILE       = FILE2
db <- dbConnect(RSQLite::SQLite(), FILE)
#tb <- dbListTables(db)
d <- dbReadTable(db, "20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_ClimateData_batch02")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
summary(clim)
dbDisconnect(db)
clim2 <- clim

FILE       = FILE3
db <- dbConnect(RSQLite::SQLite(), FILE)
#tb <- dbListTables(db)
d <- dbReadTable(db, "20200506_OnGrid_CAMARU_SurveyYear_2019_2020_Points_ClimateData_batch03")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
summary(clim)
dbDisconnect(db)

clim_off <- rbind(clim2, clim)

## points

FILE       = FILE2
TABLE      = "20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_batch02"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = FALSE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]

source("02-long.R")
d_long_pt2 <- d_long

FILE       = FILE3
TABLE      = "20200506_OnGrid_CAMARU_SurveyYear_2019_2020_Points_batch03"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = FALSE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]

source("02-long.R")
d_long_pt3 <- d_long

compare_sets(colnames(d_long_pt2), colnames(d_long_pt3))
setdiff(colnames(d_long_pt3), colnames(d_long_pt2))
isect <- intersect(colnames(d_long_pt3), colnames(d_long_pt2))
d_long_pt_off <- rbind(d_long_pt2[,isect], d_long_pt3[,isect])

##  large buffers

FILE       = FILE2
TABLE      = "20200429_OffGrid_Sites_SurveyYear_2019_2020_Buffers_batch02"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
d_long_1km2 <- d_long
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_wide_1km2 <- d_wide

FILE       = FILE3
TABLE      = "20200506_OnGrid_CAMARU_SurveyYear_2019_2020_Buffers_batch03"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
d_long_1km3 <- d_long
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_wide_1km3 <- d_wide

compare_sets(colnames(d_long_1km2), colnames(d_long_1km3))
setdiff(colnames(d_long_1km3), colnames(d_long_1km2))
isect <- intersect(colnames(d_long_1km2), colnames(d_long_1km3))
d_long_1km_off <- rbind(d_long_1km2[,isect], d_long_1km3[,isect])

d_wide_1km_off <- d_wide_1km2
d_wide_1km_off[[1]] <- rbind(d_wide_1km2[[1]], d_wide_1km3[[1]])
d_wide_1km_off[[2]] <- rbind(d_wide_1km2[[2]], d_wide_1km3[[2]])
d_wide_1km_off[[3]] <- rbind(d_wide_1km2[[3]], d_wide_1km3[[3]])
d_wide_1km_off[[4]] <- rbind(d_wide_1km2[[4]], d_wide_1km3[[4]])
d_wide_1km_off[[5]] <- c(d_wide_1km2[[5]], d_wide_1km3[[5]])




## small buffers

FILE       = FILE2
TABLE      = "20200429_OffGrid_Sites_SurveyYear_2019_2020_Buffers_batch02"
SUB_COL    = "Section"
SUB_VAL    = c("0-56m", "56-100m", "100-150m")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
table(fire_year_update)
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long2 <- d_long
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_wide2 <- d_wide

FILE       = FILE3
TABLE      = "20200506_OnGrid_CAMARU_SurveyYear_2019_2020_Buffers_batch03"
SUB_COL    = "Section"
SUB_VAL    = c("0-56m", "56-100m", "100-150m")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "Survey_Year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
d_long3 <- d_long
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_wide3 <- d_wide

compare_sets(colnames(d_long2), colnames(d_long3))
setdiff(colnames(d_long3), colnames(d_long2))
isect <- intersect(colnames(d_long2), colnames(d_long3))
d_long_off <- rbind(d_long2[,isect], d_long3[,isect])

d_wide_off <- d_wide2
d_wide_off[[1]] <- rbind(d_wide2[[1]], d_wide3[[1]])
d_wide_off[[2]] <- rbind(d_wide2[[2]], d_wide3[[2]])
d_wide_off[[3]] <- rbind(d_wide2[[3]], d_wide3[[3]])
d_wide_off[[4]] <- rbind(d_wide2[[4]], d_wide3[[4]])
d_wide_off[[5]] <- c(d_wide2[[5]], d_wide3[[5]])


##  150m Cam/ARU buffers on grid from batch 1 ---------------------------------

#rm(list=ls())
#source("00-setup.R")

FILE       = "s:/GC_eric/FromEric/Sites_summaries/Round2020/20200224_SC_Sites_SummaryTables_Round_2020.sqlite"
db <- dbConnect(RSQLite::SQLite(), FILE)
(tb <- dbListTables(db))

d <- dbReadTable(db, "20200224_CAMARU_SurveyYear_2019_Points_ClimateData_batch01")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
dbDisconnect(db)
clim_on <- clim

## points

TABLE      = "20200224_CAMARU_SurveyYear_2019_Points_batch01"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID"
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = FALSE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
d_long_pt_on <- d_long

##  large buffer around site centres

TABLE      = "20200224_CAMARU_SurveyYear_2019_Buffers_batch01"
SUB_COL    = NULL
SUB_VAL    = NULL
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1km_on <- d_long
d_wide_1km_on <- d_wide


TABLE      = "20200224_CAMARU_SurveyYear_2019_Buffers_batch01"
SUB_COL    = "Section"
SUB_VAL    = c("0-56m", "56-100m", "100-150m")
UID_COL    = "UID_old" # UID is site_year_section, UID_old is site_year
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
SAVE       = c("clim", "d_long_1km", "d_wide_1km", "d_long_pt")
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_on <- d_long
d_wide_on <- d_wide

d_long_on$Survey_Year <- d_long_on$survey_year
compare_sets(colnames(d_long_on), colnames(d_long_off))
setdiff(colnames(d_long_on), colnames(d_long_off))
setdiff(colnames(d_long_off), colnames(d_long_on))

d_long_pt_on$Survey_Year <- d_long_pt_on$survey_year
compare_sets(colnames(d_long_pt_on), colnames(d_long_pt_off))
setdiff(colnames(d_long_pt_on), colnames(d_long_pt_off))
setdiff(colnames(d_long_pt_off), colnames(d_long_pt_on))

d_long_1km_on$Survey_Year <- d_long_1km_on$survey_year
compare_sets(colnames(d_long_1km_on), colnames(d_long_1km_off))
setdiff(colnames(d_long_1km_on), colnames(d_long_1km_off))
setdiff(colnames(d_long_1km_off), colnames(d_long_1km_on))


clim <- rbind(clim_on, clim_off)
ii <- intersect(colnames(d_long_off), colnames(d_long_on))
d_long <- rbind(d_long_on[,ii], d_long_off[,ii])
ii <- intersect(colnames(d_long_pt_off), colnames(d_long_pt_on))
d_long_pt <- rbind(d_long_pt_on[,ii], d_long_pt_off[,ii])
ii <- intersect(colnames(d_long_1km_off), colnames(d_long_1km_on))
d_long_1km <- rbind(d_long_1km_on[,ii], d_long_1km_off[,ii])

d_wide <- d_wide_on
d_wide[[1]] <- rbind(d_wide_on[[1]], d_wide_off[[1]])
d_wide[[2]] <- rbind(d_wide_on[[2]], d_wide_off[[2]])
d_wide[[3]] <- rbind(d_wide_on[[3]], d_wide_off[[3]])
d_wide[[4]] <- rbind(d_wide_on[[4]], d_wide_off[[4]])
d_wide[[5]] <- c(d_wide_on[[5]], d_wide_off[[5]])

d_wide_1km <- d_wide_1km_on
d_wide_1km[[1]] <- rbind(d_wide_1km_on[[1]], d_wide_1km_off[[1]])
d_wide_1km[[2]] <- rbind(d_wide_1km_on[[2]], d_wide_1km_off[[2]])
d_wide_1km[[3]] <- rbind(d_wide_1km_on[[3]], d_wide_1km_off[[3]])
d_wide_1km[[4]] <- rbind(d_wide_1km_on[[4]], d_wide_1km_off[[4]])
d_wide_1km[[5]] <- c(d_wide_1km_on[[5]], d_wide_1km_off[[5]])

SAVE       = c("clim", "d_long_1km", "d_wide_1km", "d_long_pt")
COMMENTS   = "Veg+HF summary for on and offgrid sites - veg61+vhf // 2020-05-01"
OUTPUT     = "s:/AB_data_v2020/data/analysis/site/veg-hf_CAMARU-2019_Veg61-vHF.Rdata"

source("04-save.R")

## ------------------ batch 5: 2017-2018 site updates

rm(list=objects())
od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
FILE       = paste0("s:/GC_eric/FromEric/Sites_summaries/Round2020/",
                    "20200826_SC_Sites_SummaryTables_Round_2020_batch05.sqlite")


db <- dbConnect(RSQLite::SQLite(), FILE)
(tb <- dbListTables(db))

d <- dbReadTable(db, "20200826_SummaryTable_terr_sites_2017_points_Climate_batch05")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
clim17 <- clim

d <- dbReadTable(db, "20200826_SummaryTable_terr_sites_2018_points_Climate_batch05")
clim <- data.frame(as.matrix(Xtab(RASTERVALU ~ UID + param, d)))
clim$PET <- clim$Eref
clim$Eref <- NULL
clim$pAspen <- clim$Populus_tremuloides_brtpred_nofp
clim$Populus_tremuloides_brtpred_nofp <- NULL
clim18 <- clim

dbDisconnect(db)

## 2017

TABLE      = "20200826_SummaryTable_terr_sites_2017_batch05_Buffer"
SUB_COL    = NULL
SUB_VAL    = NULL
VEG_COL    = "Combined_ChgByCWCS"
BASE_YR    = "survey_year"
AREA_COL   = "Shape_Area"
AREA       = TRUE
TOL        = 0
UNROUND    = FALSE

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
tmp <- strsplit(as.character(d$UID), "_")
table(sapply(tmp, "[[", 3))
d$SUB <- sapply(tmp, "[[", 3)
d0 <- d
x17 <- nonDuplicated(d0, UID_old, TRUE)

## 564m buffer
UID_COL    = "UID_old"
d[[UID_COL]] <- droplevels(d[[UID_COL]])
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1km_17 <- d_long
d_wide_1km_17 <- d_wide

## 1ha
UID_COL    = "UID_old"
d <- d0[d0$SUB %in% c("NE", "NW", "SE", "SW"),]
d[[UID_COL]] <- droplevels(d[[UID_COL]])
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1ha_17 <- d_long
d_wide_1ha_17 <- d_wide

## quadrants (1/4 ha)
UID_COL    = "UID"
d <- d0[d0$SUB %in% c("NE", "NW", "SE", "SW"),]
d[[UID_COL]] <- droplevels(d[[UID_COL]])
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_qha_17 <- d_long
d_wide_qha_17 <- d_wide

## 2018

TABLE      = "20200826_SummaryTable_terr_sites_2018_batch05_Buffer"

source("01-data.R")
fire_year_update <- !is.na(d$Origin_Year) & d$WILDFIRE_YEAR > d$Origin_Year & d$WILDFIRE_YEAR < d[[BASE_YR]]
d$Origin_Year[fire_year_update] <- d$WILDFIRE_YEAR[fire_year_update]
tmp <- strsplit(as.character(d$UID), "_")
table(sapply(tmp, "[[", 3))
d$SUB <- sapply(tmp, "[[", 3)
d0 <- d
x18 <- nonDuplicated(d0, UID_old, TRUE)

## 564m buffer
UID_COL    = "UID_old"
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1km_18 <- d_long
d_wide_1km_18 <- d_wide

## 1ha
UID_COL    = "UID_old"
d <- d0[d0$SUB %in% c("NE", "NW", "SE", "SW"),]
d[[UID_COL]] <- droplevels(d[[UID_COL]])
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_1ha_18 <- d_long
d_wide_1ha_18 <- d_wide

## quadrants (1/4 ha)
UID_COL    = "UID"
d <- d0[d0$SUB %in% c("NE", "NW", "SE", "SW"),]
d[[UID_COL]] <- droplevels(d[[UID_COL]])
source("02-long.R")
source("03-wide.R")
summary(rowSums(d_wide[[1]]))
d_long_qha_18 <- d_long
d_wide_qha_18 <- d_wide

## combine the 2 years

x17 <- x17[rownames(clim17),]
x18 <- x18[rownames(clim18),]

clim17 <- data.frame(clim17, x17[,c("survey_year", "NSRNAME", "NRNAME", "LUF_NAME")])
clim18 <- data.frame(clim18, x18[,c("survey_year", "NSRNAME", "NRNAME", "LUF_NAME")])

clim <- rbind(clim17, clim18)
table(clim$survey_year)

d_wide_1km <- d_wide_1km_17
for (i in 1:4)
    d_wide_1km[[i]] <- rbind(d_wide_1km_17[[i]], d_wide_1km_18[[i]])
d_wide_1km[[5]] <- c(d_wide_1km_17[[5]], d_wide_1km_18[[5]])

d_wide_1ha <- d_wide_1ha_17
for (i in 1:4)
    d_wide_1ha[[i]] <- rbind(d_wide_1ha_17[[i]], d_wide_1ha_18[[i]])
d_wide_1ha[[5]] <- c(d_wide_1ha_17[[5]], d_wide_1ha_18[[5]])

d_wide_qha <- d_wide_qha_17
for (i in 1:4)
    d_wide_qha[[i]] <- rbind(d_wide_qha_17[[i]], d_wide_qha_18[[i]])
d_wide_qha[[5]] <- c(d_wide_qha_17[[5]], d_wide_qha_18[[5]])

compare_sets(rownames(clim), rownames(d_wide_1km[[1]]))
compare_sets(rownames(clim), rownames(d_wide_1ha[[1]]))
compare_sets(rownames(clim), rownames(d_wide_qha[[1]]))


save(clim, d_wide_1km, d_wide_1ha, d_wide_qha,
    file="s:/AB_data_v2020/data/analysis/site/veg-hf_SITES-2017-2018_Veg61-vHF.Rdata")
save(clim, d_wide_1km, d_wide_1ha, d_wide_qha,
    file="d:/abmi/AB_data_v2020/data/analysis/site/veg-hf_SITES-2017-2018_Veg61-vHF.Rdata")



## 2018 CAMARU data

rm(list=ls())
od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")

FILE = "d:/abmi/AB_data_v2020/data/raw/birds/20190129_SummaryTables_CAMARU_2017_2018_Veg61_vHFSPOT2017.sqlite"

db <- dbConnect(RSQLite::SQLite(), FILE)
dbListTables(db)
d <- dbReadTable(db, "Summary_Buffers")
dc <- dbReadTable(db, "Summary_Points")

dbDisconnect(db)

table(d$Survey_Year)
table(d$Section)

