library(cure4insect)
load_common_data()
.c4if=cure4insect:::.c4if
.c4is=cure4insect:::.c4is

rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))
KT <- .c4if$KT
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
load("d:/abmi/AB_data_v2019/misc/overlap/OverlapReg.RData")
rownames(OverlapReg) <- OverlapReg$LinkID
all(rownames(KT) == rownames(kgrid))

OverlapReg$pAspen <- kgrid[rownames(OverlapReg), "pAspen"]
OverlapReg$wN <- OverlapReg$pAspen / (OverlapReg$pAspen + (1-OverlapReg$pForest))

kgrid$wN <- ifelse(KT$reg_nr == "Grassland", 0, 1)
kgrid[rownames(OverlapReg), "wN"] <- OverlapReg$wN

rol <- .make_raster(kgrid$wN, kgrid, rpa)
plot(rol)

writeRaster(rol, "~/repos/cure4insect/inst/extdata/wNorth.tif")

#kgrid$wN2 <- ifelse(KT$reg_nr == "Grassland", 0, 1)
#kgrid[rownames(OverlapReg), "wN2"] <- OverlapReg$pAspen
#rol2 <- .make_raster(kgrid$wN2, kgrid, rpa)
#plot(rol2)


## Drat
library(drat)
options("dratRepo"="~/repos/ABbiodiversity-drat")
# only on setup!
# pruneRepo(remove=TRUE)

insertPackage("d:/abmi/AB_data_v2019/misc/cure4insect_0.1-1.tar.gz")
