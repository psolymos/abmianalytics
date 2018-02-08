## processing common data (kgrid, areas, species lookups)

library(mefa4)
library(pbapply)
## HF 2012
load("e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata")

load("e:/peter/AB_data_v2017/data/analysis/kgrid_table_km.Rdata")
stopifnot(all(rownames(kgrid) == rownames(dd1km_pred[[1]])))

library(cure4insect)
load_common_data()
.c4if=cure4insect:::.c4if

## do vegs
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
veg_cf <- .c4if$CFbirds$marginal$veg
veg1 <- groupSums(dd1km_pred[[1]], 2, tv0[colnames(dd1km_pred[[1]]), "UseInAnalysisCoef"])
veg0 <- groupSums(dd1km_pred[[2]], 2, tv0[colnames(dd1km_pred[[2]]), "UseInAnalysisCoef"])
veg_cf[is.na(veg_cf)] <- 0
veg_cf[,"SoftLin"] <- log(rowMeans(exp(veg_cf[,c("Shrub", "GrassHerb")])))
veg_cf[,"HardLin"] <- -10
cn1 <- intersect(colnames(veg1), colnames(veg_cf))
cn0 <- intersect(colnames(veg0), colnames(veg_cf))
veg1 <- veg1[,cn1]
veg0 <- veg0[,cn0]
rs1 <- rowSums(veg1)
rs1[rs1==0] <- 1
rs0 <- rowSums(veg0)
rs0[rs0==0] <- 1
veg1 <- veg1/rs1
veg0 <- veg0/rs1
vcf1 <- exp(veg_cf[,colnames(veg1)])
vcf0 <- exp(veg_cf[,colnames(veg0)])


mat_veg1 <- pbapply(vcf1, 1, function(z) as.numeric(veg1 %*% z))
mat_veg0 <- pbapply(vcf0, 1, function(z) as.numeric(veg0 %*% z))

for (i in 1:ncol(mat_veg1)) {
    q <- quantile(mat_veg1[,i], 0.99)
    mat_veg1[mat_veg1[,i] > q,i] <- q
    q <- quantile(mat_veg0[,i], 0.99)
    mat_veg0[mat_veg0[,i] > q,i] <- q
}

## do soils
ts0 <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf.csv")
soil_cf <- .c4if$CFbirds$marginal$soil
pA_cf <- .c4if$CFbirds$marginal$paspen
pA <- kgrid$pAspen
soil1 <- groupSums(dd1km_pred[[3]], 2, ts0[colnames(dd1km_pred[[3]]), "UseInAnalysis"])
soil0 <- groupSums(dd1km_pred[[4]], 2, ts0[colnames(dd1km_pred[[4]]), "UseInAnalysis"])

soil_cf[,"SoftLin"] <- log(rowMeans(exp(soil_cf), na.rm=TRUE)) # SoftLin is NA
soil_cf[,"HardLin"] <- -10
soil_cf[is.na(soil_cf)] <- 0
cn1 <- intersect(colnames(soil1), colnames(soil_cf))
cn0 <- intersect(colnames(soil0), colnames(soil_cf))
soil1 <- soil1[,cn1]
soil0 <- soil0[,cn0]
rs1 <- rowSums(soil1)
rs1[rs1==0] <- 1
rs0 <- rowSums(soil0)
rs0[rs0==0] <- 1
soil1 <- soil1/rs1
soil0 <- soil0/rs1
scf1 <- exp(soil_cf[,colnames(soil1)])
scf0 <- exp(soil_cf[,colnames(soil0)])

is_soil <- kgrid$pSoil > 0
pAterm <- exp(pbsapply(pA_cf, function(z) pA[is_soil] * z))
mat_soil1 <- pbapply(scf1, 1, function(z) as.numeric(soil1[is_soil,] %*% z)) * pAterm
mat_soil0 <- pbapply(scf0, 1, function(z) as.numeric(soil0[is_soil,] %*% z)) * pAterm

for (i in 1:ncol(mat_soil1)) {
    q <- quantile(mat_soil1[,i], 0.99)
    mat_soil1[mat_soil1[,i] > q,i] <- q
    q <- quantile(mat_soil0[,i], 0.99)
    mat_soil0[mat_soil0[,i] > q,i] <- q
}

save(mat_veg1, mat_veg0, mat_soil1, mat_soil0,
    file="e:/peter/AB_data_v2017/data/analysis/birds-hab-only-km-preds.Rdata")
