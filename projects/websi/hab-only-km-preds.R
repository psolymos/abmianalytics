## processing common data (kgrid, areas, species lookups)

library(mefa4)
library(pbapply)
## HF 2012
load("e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata")
tv0 <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")

library(cure4insect)
load_common_data()
.c4if=cure4insect:::.c4if
veg_cf <- .c4if$CFbirds$marginal$veg
soil_cf <- .c4if$CFbirds$marginal$soil

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
