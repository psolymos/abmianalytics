od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)

ROOT <- "s:/AB_data_v2020/data/raw/veghf/reports/norbord"
f <- file.path(ROOT,
  "IC_Reporting_NorbordReportingAreaFinal_3regions_10TM_Reporting_Vegetation_HF2018SC_rawdata.csv")

d0 <- read.csv(f)
d0$Shape_Area <- d0$Area_ha * 100^2
d0$Soil_Type_1 <- "UNK"
d0$FEATURE_TY[d0$FEATURE_TY == "WELL-UNKNOWN"] <- "WELL-OTHER"
d0$Origin_Year <- 2020-80
str(d0)

d <- make_vegHF_wide_v6(d0,
    col.label="Region",
    col.year=2018,
    col.HFyear="YEAR",
    col.HABIT="Vegetation",
    col.SOIL="Soil_Type_1",
    HF_fine=TRUE, wide=TRUE)$veg_current
data.frame(A=colSums(d))

d <- d[,!(colnames(d) %in%
  c("Decid0", "DecidR", "Decid1", "Decid2", "Decid3", "Decid5",
  "Decid6", "Decid7", "Decid8", "Decid9", "Mixedwood0", "MixedwoodR",
  "Mixedwood1", "Mixedwood2", "Mixedwood3", "Mixedwood5", "Mixedwood6",
  "Mixedwood7", "Mixedwood8", "Mixedwood9", "Pine0", "PineR", "Pine1",
  "Pine2", "Pine3", "Pine5", "Pine6", "Pine7", "Pine8", "Pine9",
  "Spruce0", "SpruceR", "Spruce1", "Spruce2", "Spruce3", "Spruce5",
  "Spruce6", "Spruce7", "Spruce8", "Spruce9", "TreedBog0", "TreedBogR",
  "TreedBog1", "TreedBog2", "TreedBog3", "TreedBog5", "TreedBog6",
  "TreedBog7", "TreedBog8", "TreedBog9", "TreedFen0", "TreedFenR",
  "TreedFen1", "TreedFen2", "TreedFen3", "TreedFen5", "TreedFen6",
  "TreedFen7", "TreedFen8", "TreedFen9", "TreedSwamp0", "TreedSwampR",
  "TreedSwamp1", "TreedSwamp2", "TreedSwamp3", "TreedSwamp5", "TreedSwamp6",
  "TreedSwamp7", "TreedSwamp8", "TreedSwamp9", "CutBlocks"))]

for (i in c("Decid", "Mixedwood", "Pine", "Spruce", "TreedBog", "TreedFen",
"TreedSwamp")) {
  colnames(d)[colnames(d) == paste0(i, "4")] <- i
}

write.csv(d, file=file.path(ROOT, "norbord-veghf2018.csv"))

