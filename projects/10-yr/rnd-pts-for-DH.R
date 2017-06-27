load("e:/peter/AB_data_v2017/data/analysis/Random samples info for North and South.RData")

lt <- read.csv("c:/Users/Peter/repos/abmispecies/_data/birds.csv")

SPP <- as.character(lt$AOU[lt$veghf.north])
rnd_table_curr_N <- matrix(NA, nrow(Random.N), length(SPP))
dimnames(rnd_table_curr_N) <- list(rownames(Random.N), SPP)
rnd_table_ref_N <- rnd_table_curr_N
for (spp in SPP) {

load(paste0("e:/peter/AB_data_v2016/out/birds/pred1combined-as-ES/", spp, ".Rdata"))
NN <- km2[rownames(Random.N),]
rnd_table_curr_N[,spp] <- NN[,"Curr"]
rnd_table_ref_N[,spp] <- NN[,"Ref"]

}
colnames(rnd_table_curr_N) <- lt[SPP, "sppid"]
colnames(rnd_table_ref_N) <- lt[SPP, "sppid"]

SPP <- as.character(lt$AOU[lt$soilhf.south])
rnd_table_curr_S <- matrix(NA, nrow(Random.S), length(SPP))
dimnames(rnd_table_curr_S) <- list(rownames(Random.S), SPP)
rnd_table_ref_S <- rnd_table_curr_S
for (spp in SPP) {

load(paste0("e:/peter/AB_data_v2016/out/birds/pred1combined-as-ES/", spp, ".Rdata"))
NN <- km2[rownames(Random.S),]
rnd_table_curr_S[,spp] <- NN[,"Curr"]
rnd_table_ref_S[,spp] <- NN[,"Ref"]

}
colnames(rnd_table_curr_S) <- lt[SPP, "sppid"]
colnames(rnd_table_ref_S) <- lt[SPP, "sppid"]

save(rnd_table_curr_N, rnd_table_ref_N, rnd_table_curr_S, rnd_table_ref_S,
    file="e:/peter/sppweb2017/birds-rnd-1000-pixel.Rdata")
