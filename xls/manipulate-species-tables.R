library(xlsx)
#fn <- "e:/peter/sppweb2016/round01/tables/ABMI-species-v3.2_mammals.xlsx"
fn <- "e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_lichens.xlsx"
#fn <- "e:/peter/sppweb2016/round01/tables/ABMI-species-v3.2_vplants-nonnative.xlsx"
#fn <- "e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_mites.xlsx"
#fn <- "e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_vplants.xlsx"
#fn <- "e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_mosses.xlsx"

x1 <- read.xlsx(fn, sheetName="VegetationNorth")
x2 <- read.xlsx(fn, sheetName="LinearNorth")
keep1 <- which(rowSums(is.na(x1)) != ncol(x1))
keep2 <- which(rowSums(is.na(x2)) != ncol(x2))
nam_avg <- c("WhiteSpruce0", "WhiteSpruce10", "WhiteSpruce20",
    "WhiteSpruce40", "WhiteSpruce60", "WhiteSpruce80", "WhiteSpruce100",
    "WhiteSpruce120", "WhiteSpruce140", "Pine0", "Pine10", "Pine20",
    "Pine40", "Pine60", "Pine80", "Pine100", "Pine120", "Pine140",
    "Deciduous0", "Deciduous10", "Deciduous20", "Deciduous40", "Deciduous60",
    "Deciduous80", "Deciduous100", "Deciduous120", "Deciduous140",
    "Mixedwood0", "Mixedwood10", "Mixedwood20", "Mixedwood40", "Mixedwood60",
    "Mixedwood80", "Mixedwood100", "Mixedwood120", "Mixedwood140",
    "BlackSpruce0", "BlackSpruce10", "BlackSpruce20", "BlackSpruce40",
    "BlackSpruce60", "BlackSpruce80", "BlackSpruce100", "BlackSpruce120",
    "BlackSpruce140", "Larch0", "Larch10", "Larch20", "Larch40",
    "Larch60", "Larch80", "Larch100", "Larch120", "Larch140",
    "Swamp", "WetGrass", "WetShrub", "Shrub", "GrassHerb")
x3 <- x1[keep1,]
x4 <- x2[keep2,]
x4 <- x4[,c(colnames(x4)[1], "AverageCoef", "SoftLin", "SoftLin.LCI", "SoftLin.UCI",
    "HardLin", "HardLin.LCI", "HardLin.UCI")]
x4$AverageCoef <- rowMeans(x3[,nam_avg], na.rm=TRUE)
for (i in c("SoftLin", "SoftLin.LCI", "SoftLin.UCI", "HardLin", "HardLin.LCI", "HardLin.UCI")) {
    x4[,i] <- 0.1*x4[,i]+0.9*x4$AverageCoef
}
write.xlsx(x4, fn, sheetName="LinearNorth2", append=TRUE, row.names=FALSE)

x1 <- read.xlsx(fn, sheetName="SoilNontreedSouth")
x2 <- read.xlsx(fn, sheetName="LinearSouth")
keep1 <- which(rowSums(is.na(x1)) != ncol(x1))
keep2 <- which(rowSums(is.na(x2)) != ncol(x2))
nam_avg <- c("Productive", "Clay", "Saline", "RapidDrain")
x3 <- x1[keep1,]
x4 <- x2[keep2,]
x4 <- x4[,c(colnames(x4)[1], "AverageCoef", "SoftLin", "SoftLin.LCI", "SoftLin.UCI",
    "HardLin", "HardLin.LCI", "HardLin.UCI")]
x4$AverageCoef <- rowMeans(x3[,nam_avg], na.rm=TRUE)
for (i in c("SoftLin", "SoftLin.LCI", "SoftLin.UCI", "HardLin", "HardLin.LCI", "HardLin.UCI")) {
    x4[,i] <- 0.1*x4[,i]+0.9*x4$AverageCoef
}
write.xlsx(x4, fn, sheetName="LinearSouth2", append=TRUE, row.names=FALSE)

