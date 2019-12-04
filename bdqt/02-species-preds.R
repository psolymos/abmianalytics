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

## make sp points object
xy <- d[,c("Easting", "Northing")]
coordinates(xy) <- ~ Easting + Northing
proj4string(xy) <- proj4string(.read_raster_template())
#save(xy, file="s:/AB_data_v2019/bdqt/bdqt-poly-xy_2019-12-04.RData")

## make veg/soil/hf table (current only)
x <- d[,c("UID", "VEGHFAGEclass", "SOILHFclass")]
x$VEGHFAGEclass <- as.factor(x$VEGHFAGEclass)
x$SOILHFclass <- as.factor(x$SOILHFclass)
rm(d)
gc()

## manage labels

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
cbind(class=levels(x$VEGHFAGEclass), reclass=as.character(tv[levels(x$VEGHFAGEclass), "CoefTabs"]))
levels(x$VEGHFAGEclass) <- as.character(tv[levels(x$VEGHFAGEclass), "CoefTabs"])
compare_sets(get_levels()$veg, levels(x$VEGHFAGEclass))

x <- droplevels(x)
str(x)

## extract wNorth values
rw <- raster(system.file("extdata/wNorth.tif", package="cure4insect"))
wNorth <- extract(rw, xy)
wNorth2 <- extract(rw, xy[is.na(wNorth),], method="bilinear")
summary(wNorth2)
wNorth[is.na(wNorth)] <- wNorth2

plot(rw)
plot(xy[is.na(wNorth),], add=TRUE)
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

#save(x, file="s:/AB_data_v2019/bdqt/bdqt-poly-hab_2019-12-04.RData")

for (i in levels(x$chunk)) {
    cat(i, "\n")
    xi <- x[x$chunk == i,]
    xyi <- xy[x$chunk == i,]
    xi$chunk <- NULL
    if (!(i %in% c("O1", "O2")))
        xi$wNorth <- NULL
    save(xi, xyi, file=paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-", i, "_2019-12-04.RData"))
    gc()
}

## read in species
## O - N+S averaged
## N & S: no averaging just veg or soil
## birds: p-transformation
## store results in matrix that is saved


library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()

ch <- "O1"
type <- substr(ch, 1, 1)
load(paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-", ch, "_2019-12-04.RData"))

tax <- "birds"
spp <- "Ovenbird"


cat(spp, "in", ch, "\n")
flush.console()

object <- load_spclim_data(spp)

if (type == "S") {
    pr <- predict(object, xyi, soil=xi$SOILHFclass)$soil
}
    veg <- xi$SOILHFclass
if (type == "N") {
    pr <- predict(object, xyi, veg=xi$VEGHFAGEclass)$veg
}
if (type == "O") {
    #pr <- predict(object, xyi, veg=xi$VEGHFAGEclass, soil=xi$SOILHFclass)$comb # 31sec
    prs <- predict(object, xyi, soil=xi$SOILHFclass)$soil
    prv <- predict(object, xyi, veg=xi$VEGHFAGEclass)$veg
    pr <- xi$wNorth * prv + (1 - xi$wNorth) * prs # 28sec
}
pr[is.na(pr)] <- 0
summary(pr)

if (tax == "birds")
    pr <- p_bird(pr, area="ha", pair_adj=2)


