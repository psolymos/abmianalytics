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
      SHAPE_Area
    FROM
      `Veg61HF2017BDQT`;")

Area <- dbGetQuery(db,
    "SELECT
      UID,
      SHAPE_Area
    FROM
      `Veg61HF2017BDQT`;")

Nr <- dbGetQuery(db,
    "SELECT
      UID, Easting, Northing,
      NRNAME,
      NSRNAME
    FROM
      `Veg61HF2017BDQT`;")
rownames(Nr) <- as.character(Nr$UID)
Nr$NRNAME <- as.factor(Nr$NRNAME)
Nr$NSRNAME <- as.factor(Nr$NSRNAME)
str(Nr)

dbDisconnect(db)

#save(Area, file="s:/AB_data_v2019/bdqt/bdqt-poly-area_2020-05-12.RData")
#save(Nr, file="s:/AB_data_v2019/bdqt/bdqt-poly-nr-nsr_2020-09-16.RData")


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
    "SOILclass", "SOILHFclass", "SHAPE_Area")]
gc()


mydb <- dbConnect(RSQLite::SQLite(), "s:/AB_data_v2019/bdqt/bdqt-labels-2017hfi_2019-07-18.sqlite")

dbWriteTable(mydb, "veghf", d2, overwrite = TRUE)

dbDisconnect(mydb)


## saving intermediate files

library(cure4insect)
library(mefa4)
library(DBI)
set_options(path = "d:/abmi/reports")
load_common_data()

## load pre-processed poly data
fdb <- "s:/AB_data_v2019/bdqt/bdqt-labels-2017hfi_2019-07-18.sqlite"
con <- dbConnect(RSQLite::SQLite(), fdb)
dbListTables(con)
d <- dbReadTable(con, "veghf")
dbDisconnect(con)

rownames(d) <- d$UID

## make sp points object
xy <- d[,c("Easting", "Northing")]
coordinates(xy) <- ~ Easting + Northing
proj4string(xy) <- proj4string(.read_raster_template())

## make veg/soil/hf table (current only)
x <- d[,c("UID", "VEGHFAGEclass", "SOILHFclass")]
x$VEGHFAGEclass <- as.factor(x$VEGHFAGEclass)
x$SOILHFclass <- as.factor(x$SOILHFclass)

rm(d)
gc()

compare_sets(get_levels()$soil, levels(x$SOILHFclass))
ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]
cbind(class=levels(x$SOILHFclass), reclass=as.character(ts[levels(x$SOILHFclass), "UseInAnalysisCoarse"]))

levels(x$SOILHFclass) <- as.character(ts[levels(x$SOILHFclass), "UseInAnalysisCoarse"])
## NA will be treated as 0 (water and HFor)
x$SOILHFclass[x$SOILHFclass %in% setdiff(levels(x$SOILHFclass), get_levels()$soil)] <- NA
x$SOILHFclass <- droplevels(x$SOILHFclass)

compare_sets(get_levels()$soil, levels(x$SOILHFclass))
setdiff(get_levels()$soil, levels(x$SOILHFclass))
setdiff(levels(x$SOILHFclass), get_levels()$soil)

compare_sets(get_levels()$veg, levels(x$VEGHFAGEclass))
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
cbind(class=levels(x$VEGHFAGEclass), reclass=as.character(tv[levels(x$VEGHFAGEclass), "CoefTabsBDQT"]))

levels(x$VEGHFAGEclass) <- as.character(tv[levels(x$VEGHFAGEclass), "CoefTabsBDQT"])
## NA will be treated as 0 (SnowIce)
x$VEGHFAGEclass[x$VEGHFAGEclass %in% setdiff(levels(x$VEGHFAGEclass), get_levels()$veg)] <- NA
x$VEGHFAGEclass <- droplevels(x$VEGHFAGEclass)

compare_sets(get_levels()$veg, levels(x$VEGHFAGEclass))

x <- droplevels(x)
str(x)

## extract wNorth values
rw <- raster(system.file("extdata/wNorth.tif", package="cure4insect"))
wNorth <- extract(rw, xy)
wNorth2 <- extract(rw, xy[is.na(wNorth),], method="bilinear")
summary(wNorth2)
wNorth[is.na(wNorth)] <- wNorth2

#plot(rw)
#plot(xy[is.na(wNorth),], add=TRUE)
wNorth[is.na(wNorth)] <- 1
summary(wNorth)

x$wNorth <- wNorth
rm(wNorth, wNorth2)
gc()

## define chunks
c(N=sum(x$wNorth == 1)/nrow(x),
    S=sum(x$wNorth == 0)/nrow(x),
    O=sum(x$wNorth < 1 & x$wNorth > 0)/nrow(x))
x$chunk <- factor(NA, c("S1", "O1", "O2", paste0("N", 1:11)))
x$chunk[x$wNorth == 0] <- "S1"
x$chunk[x$wNorth < 1 & x$wNorth > 0] <- c(rep("O1", 5894543), rep("O2", 5894542))
x$chunk[is.na(x$chunk)] <- sample(paste0("N", 1:11), sum(is.na(x$chunk)), replace=TRUE)
table(x$chunk, useNA="a")/10^6

#load("s:/AB_data_v2019/bdqt/bdqt-poly-xy_2019-12-10.RData")
#load("s:/AB_data_v2019/bdqt/bdqt-poly-hab_2019-12-10.RData")
#x$chunko=x$chunk
levels(x$chunk) <- c(levels(x$chunk), "O3")
x$chunk[x$chunk %in% c("O1", "O2")] <- sample(c("O1", "O2", "O3"),
  sum(x$chunk %in% c("O1", "O2")), replace=TRUE)
#table(x$chunk,x$chunko)

all(rownames(x) == as.character(x$UID))
all(rownames(xy) == as.character(x$UID))
#save(xy, file="s:/AB_data_v2019/bdqt/bdqt-poly-xy_2019-12-10.RData")
#save(x, file="s:/AB_data_v2019/bdqt/bdqt-poly-hab_2019-12-10.RData")

for (i in c("O1", "O2", "O3")) {#levels(x$chunk)) {
    cat(i, "\n")
    xi <- x[x$chunk == i,]
    xyi <- xy[x$chunk == i,]
    xi$chunk <- NULL
    if (!(i %in% c("O1", "O2", "O3")))
        xi$wNorth <- NULL
    save(xi, xyi, file=paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-", i, "_2019-12-10.RData"))
    gc()
}


## checking overlap zone and assigning weights as distance to grassland NR

library(sp)
library(rgdal)
library(rgeos)

p <- readOGR("s:/Base_shapefiles/NSR/Natural_Regions_Subregions_of_Alberta.shp")
pp <- p[p@data$NRNAME == "Grassland",]
pp <- gUnaryUnion(pp, id = rep(1, nrow(pp@data)))
p <- gUnaryUnion(p, id=p@data$NRNAM)

PART <- "O1"
f1 <- paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-",
  PART, "_2019-12-10.RData")
f2 <- paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-",
  PART, "_2020-09-16.RData")
load(f1)
pp <- spTransform(pp, proj4string(xyi))
p <- spTransform(p, proj4string(xyi))

q <- quantile(1:nrow(xi), seq(0, 1, 1/4000))
k <- cut(1:nrow(xi), q, include.lowest=TRUE, labels=FALSE)
xi$wDistance <- 0

for (i in 1:4000) {
  if (i %% 10 == 0)
    cat(".")
  if (i %% 100 == 0)
    cat(i, "\n")
  flush.console()
  xy <- xyi[k==i,]
  d <- gDistance(xy, pp, byid=TRUE) / 10^3
  d <- d/50
  d[d > 1] <- 1
  xi$wDistance[k==i] <- d
}

save(xi, xyi, file=f2)

