#e:/peter/AB_data_v2018/data/raw/veghf/abmi

HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))
meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

## 2019

## 1 ha in 4 x 0.25ha quadrants
f <- "d:/abmi/AB_data_v2019/data/analysis/alpac/20191210_SummaryTables_Alpac_veghf_2010_2016.sqlite"
db <- dbConnect(RSQLite::SQLite(), f)
dbListTables(db)
d <- dbReadTable(db, "Summary")
dbDisconnect(db)

str(d)

d2010 <- d
d2010$FEATURE_TY <- d2010$FEATURE_TY_2010
d2010$Origin_Year <- d2010$Origin_Year_2010cond
d2010$YEAR <- d2010$YEAR_2010

table(d2010$Origin_Year == d$Origin_Year)
ii <- d2010$Origin_Year == d$Origin_Year &
    d2010$Origin_Year < 2000 &
    d$Origin_Year < 2000 &
    d2010$Origin_Year %% 10 == 0
table(ii)
set.seed(1)
rnd <- sample.int(10, sum(ii), replace=TRUE)-5
d2010$Origin_Year[ii] <- d2010$Origin_Year[ii] + rnd
d$Origin_Year[ii] <- d$Origin_Year[ii] + rnd

dd2010 <- make_vegHF_wide_v6(d2010,
    col.label="GRID_LABEL",
    col.year=2010,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    wide=FALSE, sparse=TRUE, HF_fine=TRUE) # use refined classes
dd2016 <- make_vegHF_wide_v6(d,
    col.label="GRID_LABEL",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    wide=FALSE, sparse=TRUE, HF_fine=TRUE) # use refined classes

dd <- d[,c("GRID_LABEL", "Shape_Area")]
dd$veghf0 <- as.character(dd2010$VEGAGEclass)
dd$veghf10 <- as.character(dd2010$VEGHFAGEclass)
dd$veghf16 <- as.character(dd2016$VEGHFAGEclass)
#dd$veghfch <- dd$veghf10
#i <- dd$veghf10 != dd$veghf16
#table(i) / 10^6
#dd$veghfch[i] <- paste0(dd$veghf10[i], "->", dd$veghf16[i])
dd$veghfch <- paste0(dd$veghf10, "->", dd$veghf16)

if (FALSE) {
## troubleshooting the ages
table(was=dd2010$AgeCr, is=dd2016$AgeCr)
lev1 <- c("was_R", "was_1", "was_2", "was_3", "was_4", "was_5",
    "was_6", "was_7", "was_8", "was_9", "was_0", "was_")
lev2 <- c("is_R", "is_1", "is_2", "is_3", "is_4", "is_5",
    "is_6", "is_7", "is_8", "is_9", "is_0", "is_")
dd$age10 <- factor(paste0("was_", dd2010$AgeCr), lev1)
dd$age16 <- factor(paste0("is_", dd2016$AgeCr), lev2)
x <- as.matrix(Xtab(Shape_Area ~ age10 + age16, dd)) / 10^6
round(x, 0)

aa=sum_by(d$Shape_Area/10^6, d$Origin_Year)
aa=data.frame(aa[order(aa[,1], decreasing = TRUE),])
aa$perc <- 100*aa[,1]/sum(aa[,1])
#data.frame(table(d$Origin_Year_2010))

dd$yr <- d$Origin_Year
ii <- dd$yr %% 10 == 0 & dd$yr < 2000
dd$yr2 <- dd$yr
dd$yr2[ii] <- dd$yr[ii] + sample.int(10, sum(ii), replace=TRUE)-5
plot(table(dd$yr2[dd$yr < 9999]))
}

## using AlPac FMA only
#trVeg <- Xtab(Shape_Area ~ GRID_LABEL + veghfch, dd[d$FMA_NAME=="Alberta-Pacific Forest Industries Inc.",])
## using AlPac FMA + other bits
trVeg <- Xtab(Shape_Area ~ GRID_LABEL + veghfch, dd)
dim(trVeg)


tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
hf <- rownames(tv)[tv$is_HF]

Tot <- data.frame(veghfch=colnames(trVeg), area=colSums(trVeg)/10^6, prop=colSums(trVeg)/sum(trVeg))
tmp <- strsplit(colnames(trVeg), "->")
Tot$veghf10 <- sapply(tmp, "[[", 1)
Tot$veghf16 <- sapply(tmp, "[[", 2)
Tot$changed <- Tot$veghf10 != Tot$veghf16
Tot$washf <- Tot$veghf10 %in% hf
Tot$ishf <- Tot$veghf16 %in% hf

## ~5km^2 was HF and now is not --> stamp old HF over new
table(Tot[Tot$washf & ! Tot$ishf,"veghf10"])
sum(Tot[Tot$washf & ! Tot$ishf,"area"])
#Tot[Tot$washf & ! Tot$ishf,"veghf16"] <- Tot[Tot$washf & ! Tot$ishf,"veghf10"]
#Tot$changed <- Tot$veghf10 != Tot$veghf16
#Tot$washf <- Tot$veghf10 %in% hf
#Tot$ishf <- Tot$veghf16 %in% hf

## CC-nt forest ~ 1km^2 --> override by 2016 condition
i <- Tot$veghf10 %in% setdiff(Tot$veghf10, rownames(tv))
table(Tot[i,"veghf10"])
sum(Tot[i,"area"])
Tot[i,"veghf10"] <- substr(Tot[i,"veghf10"], 3, nchar(Tot[i,"veghf10"]))

i <- Tot$veghf10 %in% setdiff(Tot$veghf10, rownames(tv))
table(Tot[i,"veghf10"])
sum(Tot[i,"area"])
Tot[i,"veghf10"] <- substr(Tot[i,"veghf10"], 1, nchar(Tot[i,"veghf10"])-1)

## looking at ages & sectors
Tot$sect10 <- tv$Sector61[match(Tot$veghf10, rownames(tv))]
Tot$sect16 <- tv$Sector61[match(Tot$veghf16, rownames(tv))]
Tot$age10 <- tv$AGE[match(Tot$veghf10, rownames(tv))]
Tot$age16 <- tv$AGE[match(Tot$veghf16, rownames(tv))]

summary(Tot)
## unknown ages
sum(Tot[endsWith(Tot$veghf10, "0"), "area"])
sum(Tot[endsWith(Tot$veghf16, "0"), "area"])
## unknown ages in CC
sum(Tot[startsWith(Tot$veghf10, "CC") & endsWith(Tot$veghf10, "0"), "area"])
sum(Tot[startsWith(Tot$veghf10, "CC") & endsWith(Tot$veghf16, "0"), "area"])
Tot$veghf16[startsWith(Tot$veghf10, "CC") & endsWith(Tot$veghf16, "0")] <-
    Tot$veghf10[startsWith(Tot$veghf10, "CC") & endsWith(Tot$veghf16, "0")]

table(wasHF=Tot$washf, isHF=Tot$ishf)
table(wasHF=Tot$washf, isHF=Tot$ishf, changed=Tot$changed)

## apply back the cleanings
Tot$veghfch <- paste0(Tot$veghf10, "->", Tot$veghf16)
trVeg2 <- groupSums(trVeg, 2, Tot$veghfch)

Tot <- data.frame(veghfch=colnames(trVeg2), area=colSums(trVeg2)/10^6, prop=colSums(trVeg2)/sum(trVeg2))
tmp <- strsplit(colnames(trVeg2), "->")
Tot$veghf10 <- sapply(tmp, "[[", 1)
Tot$veghf16 <- sapply(tmp, "[[", 2)
Tot$changed <- Tot$veghf10 != Tot$veghf16
Tot$washf <- Tot$veghf10 %in% hf
Tot$ishf <- Tot$veghf16 %in% hf
Tot$sect10 <- tv$Sector61[match(Tot$veghf10, rownames(tv))]
Tot$sect16 <- tv$Sector61[match(Tot$veghf16, rownames(tv))]
Tot$age10 <- tv$AGE[match(Tot$veghf10, rownames(tv))]
Tot$age16 <- tv$AGE[match(Tot$veghf16, rownames(tv))]
Tot$ageUnk10 <- endsWith(Tot$veghf10, "0")
Tot$ageUnk16 <- endsWith(Tot$veghf16, "0")

Tot$an10 <- factor(as.character(Tot$age10), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an10) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an10 <- as.integer(Tot$an10) - 1
Tot$an16 <- factor(as.character(Tot$age16), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an16) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an16 <- as.integer(Tot$an16) - 1
Tot$diff <- Tot$an16 - Tot$an10
summary(Tot)
plot(table(Tot$diff))
sum(Tot[Tot$diff > 1, "area"])

## if 16 is not recent forestry: apply 16 cond to 10
Tot[Tot$diff > 1 & startsWith(Tot$veghf16, "CC"), "veghf10"] <-
    Tot[Tot$diff > 1 & startsWith(Tot$veghf16, "CC"), "veghf16"]
## if 16 is not forestry: apply 10 cond to 16
Tot[Tot$diff > 1 & !startsWith(Tot$veghf16, "CC"), "veghf16"] <-
    Tot[Tot$diff > 1 & !startsWith(Tot$veghf16, "CC"), "veghf10"]

## recalculate fields
Tot$changed <- Tot$veghf10 != Tot$veghf16
Tot$washf <- Tot$veghf10 %in% hf
Tot$ishf <- Tot$veghf16 %in% hf
Tot$sect10 <- tv$Sector61[match(Tot$veghf10, rownames(tv))]
Tot$sect16 <- tv$Sector61[match(Tot$veghf16, rownames(tv))]
Tot$age10 <- tv$AGE[match(Tot$veghf10, rownames(tv))]
Tot$age16 <- tv$AGE[match(Tot$veghf16, rownames(tv))]
Tot$ageUnk10 <- endsWith(Tot$veghf10, "0")
Tot$ageUnk16 <- endsWith(Tot$veghf16, "0")

Tot$an10 <- factor(as.character(Tot$age10), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an10) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an10 <- as.integer(Tot$an10) - 1
Tot$an16 <- factor(as.character(Tot$age16), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an16) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an16 <- as.integer(Tot$an16) - 1
Tot$diff <- Tot$an16 - Tot$an10

## fire that is not a recent age class 0
sum(Tot$diff < 0 & Tot$an16 != 0)
sum(Tot[Tot$diff < 0 & Tot$an16 != 0,"area"])
Tot[Tot$diff < 0 & Tot$an16 != 0,1:2]
Tot$veghf16[Tot$diff < 0 & Tot$an16 != 0] <- Tot$veghf10[Tot$diff < 0 & Tot$an16 != 0]

summary(Tot)
plot(table(Tot$diff))
sum(Tot[Tot$diff > 1, "area"])

table(wasHF=Tot$washf, isHF=Tot$ishf)
table(wasHF=Tot$washf, isHF=Tot$ishf, changed=Tot$changed)

Tot$veghfch <- paste0(Tot$veghf10, "->", Tot$veghf16)

trVeg3 <- groupSums(trVeg2, 2, Tot$veghfch)

Tot <- data.frame(veghfch=colnames(trVeg3), area=colSums(trVeg3)/10^6, prop=colSums(trVeg3)/sum(trVeg3))
tmp <- strsplit(colnames(trVeg3), "->")
Tot$veghf10 <- sapply(tmp, "[[", 1)
Tot$veghf16 <- sapply(tmp, "[[", 2)
Tot$changed <- Tot$veghf10 != Tot$veghf16
Tot$washf <- Tot$veghf10 %in% hf
Tot$ishf <- Tot$veghf16 %in% hf
Tot$sect10 <- tv$Sector61[match(Tot$veghf10, rownames(tv))]
Tot$sect16 <- tv$Sector61[match(Tot$veghf16, rownames(tv))]
Tot$age10 <- tv$AGE[match(Tot$veghf10, rownames(tv))]
Tot$age16 <- tv$AGE[match(Tot$veghf16, rownames(tv))]
Tot$ageUnk10 <- endsWith(Tot$veghf10, "0")
Tot$ageUnk16 <- endsWith(Tot$veghf16, "0")

Tot$an10 <- factor(as.character(Tot$age10), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an10) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an10 <- as.integer(Tot$an10) - 1
Tot$an16 <- factor(as.character(Tot$age16), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an16) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an16 <- as.integer(Tot$an16) - 1
Tot$diff <- Tot$an16 - Tot$an10

summary(Tot)

Tot$wascc <- startsWith(Tot$veghf10, "CC")
Tot$iscc <- startsWith(Tot$veghf16, "CC")
Tot$wasforest <- Tot$veghf10 %in% rownames(tv)[tv$is_forest]
Tot$isforest <- Tot$veghf16 %in% rownames(tv)[tv$is_forest]

Tot$type <- factor("Other", c("NewFor", "OldFor", "NewHF", "OldHF", "Fire", "Aging0", "Aging1", "Other"))
Tot$type[Tot$isforest & Tot$wasforest] <- "Aging0"
Tot$type[Tot$ishf & !Tot$washf] <- "NewHF"
Tot$type[Tot$ishf & Tot$washf] <- "OldHF"
Tot$type[Tot$iscc & !Tot$wascc] <- "NewFor"
Tot$type[Tot$iscc & Tot$wascc] <- "OldFor"
Tot$type[!Tot$ishf & !Tot$washf & Tot$diff < 0] <- "Fire"
Tot$type[!Tot$ishf & !Tot$washf & Tot$diff > 0] <- "Aging1"

table(Tot$type, useNA="a")
table(Tot$type, Tot$ishf, Tot$iscc)

aa <- sum_by(100*Tot$prop, Tot$type)
round(aa, 2)

round(as.matrix(Xtab(area ~ an10 + an16, Tot)))

save(trVeg3, Tot, file="d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-2010-2016-transitions.RData")

## reference

veghf_2010 <- groupSums(trVeg3, 2, Tot$veghf10)
veghf_2016 <- groupSums(trVeg3, 2, Tot$veghf16)
veghf_ref <- Xtab(Shape_Area ~ GRID_LABEL + veghf0, dd)
veghf_ref <- veghf_ref[rownames(veghf_2010),]
veghf_2016 <- veghf_2016[,colnames(veghf_2010)]

compare_sets(colnames(veghf_2010), colnames(veghf_2016))
compare_sets(colnames(veghf_ref), colnames(veghf_2016))
setdiff(colnames(veghf_ref), colnames(veghf_2016))
setdiff(colnames(veghf_2016), colnames(veghf_ref))

save(veghf_2010, veghf_2016, veghf_ref, file="d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-ref-2010-2016-veghf.RData")

## checking 0 and R

library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")

load_common_data()
ST <- get_species_table()
taxa <- rev(levels(ST$taxon))

load("d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-ref-2010-2016-veghf.RData")

PT <- get_id_table()

rt <- .read_raster_template()

for (AGE in c("R", 0:9)) {
#AGE <- "R"
cn <- colnames(veghf_2010)[endsWith(colnames(veghf_2010), AGE)]
v10 <- rowSums(veghf_2010[,cn]) / rowSums(veghf_2010)
v10 <- v10[match(rownames(PT), names(v10))]
v10[is.na(v10)] <- -1
summary(v10)
v16 <- rowSums(veghf_2016[,cn]) / rowSums(veghf_2016)
v16 <- v16[match(rownames(PT), names(v16))]
v16[is.na(v16)] <- -1
summary(v16)

r10 <- .make_raster(v10, PT, rt)
values(r10)[!is.na(values(r10)) & values(r10) < 0] <- NA
r10 <- trim(r10)
r16 <- .make_raster(v16, PT, rt)
values(r16)[!is.na(values(r16)) & values(r16) < 0] <- NA
r16 <- trim(r16)
rdf <- r16 - r10

col1 <- rev(viridis::viridis(100))[1:50]
col2 <- rev(viridis::cividis(100))
png(paste0("d:/tmp/alpac-ages/age-", AGE, ".png"), height=600, width=400*3, pointsize=20)
op <- par(mfrow=c(1,3), mar=c(6, 1, 4, 4))
plot(r10, main=paste("2010: age", AGE), box=FALSE, axes=FALSE, col=col1)
plot(r16, main=paste("2016: age", AGE), box=FALSE, axes=FALSE, col=col1)
plot(rdf, main=paste("2016-2010: age", AGE), box=FALSE, axes=FALSE, col=col2)
par(op)
dev.off()
}


## running transition predictions for species
## (use s:/ instead of d:/abmi/)

library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()

## load soil & veg
load("d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-ref-2010-2016-veghf.RData")
load("d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-2010-2016-transitions.RData")
rt <- .read_raster_template()
PT <- get_id_table()
XY <- get_id_locations()
ID <- rownames(trVeg3)
XY <- XY[ID,]
v <- ifelse(rownames(PT) %in% ID, 1, 0)
r <- .make_raster(v, PT, rt)
values(r)[!is.na(values(r)) & values(r) < 1] <- NA
#r <- trim(r)
trVeg3 <- trVeg3 / rowSums(trVeg3)
Tot0 <- Tot

#tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
tv <- read.csv("https://raw.githubusercontent.com/psolymos/abmianalytics/master/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
hf <- rownames(tv)[tv$is_HF]

## dealing with unknown ages (making age class 9 to be 8)
h <- Tot$veghf10[endsWith(Tot$veghf10, "9")]
h <- paste0(substr(h, 1, nchar(h)-1), "8")
Tot$veghf10[endsWith(Tot$veghf10, "9")] <- h
h <- Tot$veghf10[endsWith(Tot$veghf10, "0")]
h <- paste0(substr(h, 1, nchar(h)-1), "8")
Tot$veghf10[endsWith(Tot$veghf10, "0")] <- h

h <- Tot$veghf16[endsWith(Tot$veghf16, "9")]
h <- paste0(substr(h, 1, nchar(h)-1), "8")
Tot$veghf16[endsWith(Tot$veghf16, "9")] <- h
h <- Tot$veghf16[endsWith(Tot$veghf16, "0")]
h <- paste0(substr(h, 1, nchar(h)-1), "8")
Tot$veghf16[endsWith(Tot$veghf16, "0")] <- h

## matching the coefs used in the package
Tot$cn10 <- tv[Tot$veghf10, "CoefTabs"]
Tot$cn16 <- tv[Tot$veghf16, "CoefTabs"]

compare_sets(get_levels()$veg, colnames(groupSums(trVeg3, 2, Tot$cn10)))
compare_sets(get_levels()$veg, colnames(groupSums(trVeg3, 2, Tot$cn16)))

## reassign 'type' variable after unknown age fix
Tot$changed <- Tot$veghf10 != Tot$veghf16
Tot$washf <- Tot$veghf10 %in% hf
Tot$ishf <- Tot$veghf16 %in% hf
Tot$sect10 <- tv$Sector61[match(Tot$veghf10, rownames(tv))]
Tot$sect16 <- tv$Sector61[match(Tot$veghf16, rownames(tv))]
Tot$age10 <- tv$AGE[match(Tot$veghf10, rownames(tv))]
Tot$age16 <- tv$AGE[match(Tot$veghf16, rownames(tv))]
Tot$ageUnk10 <- endsWith(Tot$veghf10, "0")
Tot$ageUnk16 <- endsWith(Tot$veghf16, "0")

Tot$an10 <- factor(as.character(Tot$age10), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an10) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an10 <- as.integer(Tot$an10) - 1
Tot$an16 <- factor(as.character(Tot$age16), c("", "R", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
levels(Tot$an16) <- c("0", "0", "10", "20", "40", "60", "80", "100", "120", "140", "140")
Tot$an16 <- as.integer(Tot$an16) - 1
Tot$diff <- Tot$an16 - Tot$an10

Tot$wascc <- startsWith(Tot$veghf10, "CC")
Tot$iscc <- startsWith(Tot$veghf16, "CC")
Tot$wasforest <- Tot$veghf10 %in% rownames(tv)[tv$is_forest]
Tot$isforest <- Tot$veghf16 %in% rownames(tv)[tv$is_forest]

Tot$type <- factor("Other", c("NewFor", "OldFor", "NewHF", "OldHF", "Fire", "Aging0", "Aging1", "Other"))
Tot$type[Tot$isforest & Tot$wasforest] <- "Aging0"
Tot$type[Tot$ishf & !Tot$washf] <- "NewHF"
Tot$type[Tot$ishf & Tot$washf] <- "OldHF"
Tot$type[Tot$iscc & !Tot$wascc] <- "NewFor"
Tot$type[Tot$iscc & Tot$wascc] <- "OldFor"
Tot$type[!Tot$ishf & !Tot$washf & Tot$diff < 0] <- "Fire"
Tot$type[!Tot$ishf & !Tot$washf & Tot$diff > 0] <- "Aging1"
Tot$type[Tot$wascc & Tot$ishf & !Tot$iscc] <- "NewHF"

table(Re0=Tot$type, In=Tot0$type)
write.csv(Tot, row.names=FALSE, file="d:/abmi/AB_data_v2019/data/analysis/alpac/Tot-table-lookup.csv")

## making matrices for each type
mats10 <- list()
mats16 <- list()
for (i in levels(Tot$type)) {
    j <- Tot$type == i
    mats10[[i]] <- groupSums(trVeg3[,j], 2, Tot$cn10[j])
    mats10[[i]] <- mats10[[i]][,!(colnames(mats10[[i]]) %in% c("Bare", "Water"))]
    mats16[[i]] <- groupSums(trVeg3[,j], 2, Tot$cn16[j])
    mats16[[i]] <- mats16[[i]][,!(colnames(mats16[[i]]) %in% c("Bare", "Water"))]
}

rfun <- function(x) {
    x <- rowSums(pr10)
    x <- x[match(rownames(PT), ID)]
    x[is.na(x)] <- 0
    r10 <- .make_raster(x, PT, rt)
    trim(mask(r10, r))
}

## this is subset: change the subset vector here
tmp <- read.csv("d:/abmi/AB_data_v2019/data/analysis/alpac/AlpacFMA-LinkID_2019.csv")
ss <- list(
    AEI=rownames(trVeg3),
    FMA=as.character(tmp$LinkID)
)

## preallocating for storing results
tmp1 <- read.csv("d:/abmi/AB_data_v2019/data/analysis/alpac/ALPAC_AEI_2019-12-19.csv")
tmp2 <- read.csv("d:/abmi/AB_data_v2019/data/analysis/alpac/ALPAC_FMA_2019-12-19.csv")
## MountainChickadee should also be part of AEI not just FMA
SPP <- list(
    AEI=c("MountainChickadee", as.character(tmp1$SpeciesID[tmp1$Pred_Keep & tmp1$Model_north])),
    FMA=as.character(tmp2$SpeciesID[tmp2$Pred_Keep & tmp2$Model_north])
)
str(SPP)

Sums2010_AEI <- matrix(0, length(SPP$AEI), nlevels(Tot$type))
dimnames(Sums2010_AEI) <- list(SPP$AEI, levels(Tot$type))
Sums2016_AEI <- Sums2010_AEI

Sums2010_FMA <- matrix(0, length(SPP$FMA), nlevels(Tot$type))
dimnames(Sums2010_FMA) <- list(SPP$FMA, levels(Tot$type))
Sums2016_FMA <- Sums2010_FMA

## this is a template for storing predictions before summing
pr <- matrix(0, nrow(trVeg3), nlevels(Tot$type))
dimnames(pr) <- list(rownames(trVeg3), levels(Tot$type))

## run this in a separate process for speedup
for (spp in SPP$AEI) {
    cat("2010:", spp, which(SPP$AEI==spp), "/", length(SPP$AEI), as.character(Sys.time()), "\n")
    flush.console()
    object <- load_spclim_data(spp)
    pr10 <- pr
    for (i in levels(Tot$type)) {
        pr10[,i] <- rowSums(predict_mat(object, XY, mats10[[i]], method="bilinear")$veg)
        pr10[is.na(pr10[,i]),i] <- 0
        q <- quantile(pr10[,i], 0.99)
        pr10[pr10[,i]>q,i] <- q
    }
    Sums2010_AEI[spp,] <- colSums(pr10[ss$AEI,,drop=FALSE])
    if (spp %in% SPP$FMA)
        Sums2010_FMA[spp,] <- colSums(pr10[ss$FMA,,drop=FALSE])
}
save(Sums2010_AEI, Sums2010_FMA,
    file="d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-tr-results-2010.RData")

## run this in a separate process for speedup
for (spp in SPP$AEI) {
    cat("2016:", spp, which(SPP$AEI==spp), "/", length(SPP$AEI), as.character(Sys.time()), "\n")
    flush.console()
    object <- load_spclim_data(spp)
    pr16 <- pr
    for (i in levels(Tot$type)) {
        pr16[,i] <- rowSums(predict_mat(object, XY, mats16[[i]], method="bilinear")$veg)
        pr16[is.na(pr16[,i]),i] <- 0
        q <- quantile(pr16[,i], 0.99)
        pr16[pr16[,i]>q,i] <- q
    }
    Sums2016_AEI[spp,] <- colSums(pr16[ss$AEI,,drop=FALSE])
    if (spp %in% SPP$FMA)
        Sums2016_FMA[spp,] <- colSums(pr16[ss$FMA,,drop=FALSE])
}
save(Sums2016_AEI, Sums2016_FMA,
    file="d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-tr-results-2016.RData")


predict_mat.c4ispclim <-
function(object, xy, veg, soil, method="simple", ...)
{
    xy <- .tr_xy(xy)
    ## coefs in object are on log/logit scale, need linkinv
    fi <- if (object$taxon == "birds")
        poisson("log")$linkinv else binomial("logit")$linkinv
    if (missing(veg) && missing(soil))
        stop("veg or soil must be provided")
    xy <- spTransform(xy, proj4string(.read_raster_template()))
    if (!missing(veg) && !is.null(veg)) {
        if (is.null(object$cveg)) {
            warning(sprintf("veg based estimates are unavailable for %s", object$species))
            Nveg <- NULL
        } else {
            if (nrow(veg) != nrow(coordinates(xy)))
                stop("nrow(veg) must equal number of points in xy")
            .check(as.factor(colnames(veg)), names(object$cveg))
            iveg <- extract(object$rveg, xy, method)
            imatv <- t(array(iveg, dim(veg), dimnames(veg)))
            mveg <- object$cveg[match(colnames(veg), names(object$cveg))]
            Nveg <- fi(t(mveg + imatv)) * veg
            if (any(colnames(veg) == "SoftLin") && object$taxon == "birds" && object$version == "2017") {
                warning("veg contained SoftLin: check your assumptions")
                Nveg[,colnames(veg) == "SoftLin"] <- NA
            }
        }
    } else {
        Nveg <- NULL
    }
    if (!missing(soil) && !is.null(soil)) {
        if (is.null(object$csoil)) {
            warning(sprintf("soil based estimates are unavailable for %s", object$species))
            Nsoil <- NULL
        } else {
            if (nrow(soil) != nrow(coordinates(xy)))
                stop("nrow(veg) must equal number of points in xy")
            .check(as.factor(colnames(soil)), names(object$csoil))
            isoil <- extract(object$rsoil, xy, method)
            ## pAspen here used as habitat covariate NOT as North weight
            rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))
            ipa <- extract(rpa, xy, method)
            imats <- t(array(object$caspen * ipa + isoil, dim(soil), dimnames(soil)))
            msoil <- object$csoil[match(colnames(soil), names(object$csoil))]
            Nsoil <- fi(t(msoil + imats)) * soil
            if (any(colnames(soil) == "SoftLin") && object$taxon == "birds" && object$version == "2017") {
                warning("soil contained SoftLin: check your assumptions")
                Nsoil[,colnames(soil) == "SoftLin"] <- NA
            }
        }
    } else {
        Nsoil <- NULL
    }
    OUT <- list(veg=Nveg, soil=Nsoil)
    class(OUT) <- c("c4ippredmat")
    OUT
}


## ignore this part inside the loop
#spp <- "CommonYellowthroat"
#spp = "Ovenbird"
for (spp in SPP) {

    cat(spp, which(SPP==spp), "/", length(SPP), as.character(Sys.time()), "\n")
    flush.console()

    object <- load_spclim_data(spp)

    pr10 <- matrix(0, nrow(trVeg3), nlevels(Tot$type))
    dimnames(pr10) <- list(rownames(trVeg3), levels(Tot$type))
    pr16 <- pr10

    for (i in levels(Tot$type)) {
        pr10[,i] <- rowSums(predict_mat(object, XY, mats10[[i]])$veg)
        pr16[,i] <- rowSums(predict_mat(object, XY, mats16[[i]])$veg)
    }
    pr10[is.na(pr10)] <- 0
    pr16[is.na(pr16)] <- 0

    Sums2010[spp,] <- colSums(pr10)
    Sums2016[spp,] <- colSums(pr16)

    if (FALSE) {
        prdf <- pr16 - pr10

        z <- rbind(pr2010=colSums(pr10),
            pr2016=colSums(pr16))
        z <- rbind(z, perc_change=100*(z[2,]-z[1,])/z[1,])
        round(z, 1)

        summary(pr10)
        summary(pr16)
        summary(prdf)

        r10 <- rfun(rowSums(pr10))
        r16 <- rfun(rowSums(pr16))
        rdf <- rfun(rowSums(prdf))
        col1 <- rev(viridis::viridis(100))[1:50]
        col2 <- rev(viridis::cividis(100))
        png(paste0("d:/tmp/alpac-ages/", spp, ".png"), height=600, width=400*3, pointsize=20)
        op <- par(mfrow=c(1,3), mar=c(6, 1, 4, 4))
        plot(r10, main=paste("2010:", spp), box=FALSE, axes=FALSE, col=col1)
        plot(r16, main=paste("2016:", spp), box=FALSE, axes=FALSE, col=col1)
        plot(rdf, main=paste("2016-2010:", spp), box=FALSE, axes=FALSE, col=col2)
        par(op)
        dev.off()
    }
}


## save into Excel file: spp table, 2010, 2016

load("d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-tr-results-MAX.RData")
load("d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-tr-results-2010.RData")
load("d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-tr-results-2016.RData")
ST <- get_species_table(mregion="north")
all(rownames(ST)==SPP)
A <- 100*colMeans(groupSums(trVeg3, 2, Tot$type))
A <- A[colnames(Sums2010)]

ST$MAX <- MAX
ST$MEAN <- 0.5 * (rowSums(Sums2010)/nrow(trVeg3) + rowSums(Sums2016)/nrow(trVeg3))
ST$IN <- ST$MEAN > ST$MAX * 0.01

library(openxlsx)
l <- list(
    Perc_Area=data.frame(TransitionType=names(A), Perc_Area=A),
    Species=ST[ST$IN,],
    Sums_2010=data.frame(ST[,1:2], round(Sums2010, 6))[ST$IN,],
    Sums_2016=data.frame(ST[,1:2], round(Sums2016, 6))[ST$IN,],
    Perc_Ch_by_type=data.frame(ST[,1:2],
        round(100*(Sums2016-Sums2010)/Sums2010, 6))[ST$IN,],
    Perc_Ch_total=data.frame(ST[,1:2],
        round(100*(Sums2016-Sums2010)/rowSums(Sums2010), 6))[ST$IN,])
write.xlsx(l, "d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-tr-results.xlsx")

## some plots

z1 <- (100*(Sums2016-Sums2010)/rowSums(Sums2010))[ST$IN,]
z2 <- (100*(Sums2016-Sums2010)/Sums2010)[ST$IN,]
z3 <- 100 * t(t(z1) / A[colnames(z1)])

op <- par(mfrow=c(3,1))
boxplot(z1, ylim=c(-30, 30), col="#0000ff88", main="% of 2010 total")
abline(h=0, col=2)
boxplot(z2, ylim=c(-100, 300), col="#0000ff88", main="% of 2010 by type (~under HF)")
abline(h=0, col=2)
boxplot(z3, ylim=c(-300, 300), col="#0000ff88", main="unit effect")
abline(h=0, col=2)
abline(h=c(-100, 100), col=2, lty=2)
par(op)


x <- rowSums(Sums2016)/rowSums(Sums2010)
v <- z1[,"Fire"]
head(data.frame(v=v,ST[ST$IN,c("taxon","native")])[order(v, decreasing=TRUE),])

