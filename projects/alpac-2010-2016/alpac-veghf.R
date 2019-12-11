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
dd$veghf10 <- as.character(dd2010$VEGHFAGEclass)
dd$veghf16 <- as.character(dd2016$VEGHFAGEclass)
#dd$veghfch <- dd$veghf10
#i <- dd$veghf10 != dd$veghf16
#table(i) / 10^6
#dd$veghfch[i] <- paste0(dd$veghf10[i], "->", dd$veghf16[i])
dd$veghfch <- paste0(dd$veghf10, "->", dd$veghf16)

## using AlPac FMA only
trVeg <- Xtab(Shape_Area ~ GRID_LABEL + veghfch, dd[d$FMA_NAME=="Alberta-Pacific Forest Industries Inc.",])
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

Tot$type <- factor("Other", c("NewFor", "OldFor", "NewHF", "OldHF", "Fire", "Aging", "Other"))
Tot$type[Tot$isforest & Tot$wasforest] <- "Aging"
Tot$type[Tot$ishf & !Tot$washf & !Tot$iscc] <- "NewHF"
Tot$type[Tot$ishf & Tot$washf] <- "OldHF"
Tot$type[Tot$iscc & !Tot$wascc & !Tot$washf] <- "NewFor"
Tot$type[Tot$iscc & Tot$wascc] <- "OldFor"
Tot$type[!Tot$ishf & !Tot$washf & Tot$diff < 0] <- "Fire"
Tot$type[!Tot$ishf & !Tot$washf & Tot$diff > 0] <- "Aging"

table(Tot$type, useNA="a")

aa <- sum_by(100*Tot$prop, Tot$type)
round(aa, 2)

save(trVeg3, Tot, file="d:/abmi/AB_data_v2019/data/analysis/alpac/alpac-2010-2016-transitions.RData")
