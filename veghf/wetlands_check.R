library(mefa4)

## 2014 lookup and data
e1 <- new.env()
load("c:/p/AB_data_v2014/R/veghf_abmi_allscales.Rdata",
    envir=e1)
lt1 <- read.csv("c:/p/AB_data_v2014/lookup/VEG_HF_interim2.csv")
names(e1)
names(lt1)

## 2015 lookup and data
e2 <- new.env()
load("c:/p/AB_data_v2015/out/abmi_onoff/veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0.Rdata",
    envir=e2)
lt2 <- read.csv("c:/Users/Peter/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
names(e2)
names(lt2)

## matching

dd1 <- as.matrix(e1$dd1ha$veg_current)
dd2 <- as.matrix(e2$dd1ha$veg_current)

n1 <- strsplit(rownames(dd1), "_")
n2 <- strsplit(rownames(dd2), "_")
n2 <- sapply(n2, function(z) paste(z[4], z[5], sep="_"))
rownames(dd2) <- n2

compare.sets(rownames(dd1), rownames(dd2))
nn <- intersect(rownames(dd1), rownames(dd2))

dd1 <- dd1[nn,]
dd2 <- dd2[nn,]

## wetland classes

colnames(dd1)
colnames(dd2)

compare.sets(colnames(dd1), lt1$VEGHFAGE)
compare.sets(colnames(dd2), lt2$VEGHFAGE)

ddd1 <- groupSums(dd1, 2, lt1$Levels8)
ddd2 <- groupSums(dd2, 2, lt2$LevsForComp)

colnames(ddd2)[colnames(ddd2) == "Conif"] <- "WhiteSpruce"
colnames(ddd2)[colnames(ddd2) == "Decid"] <- "Deciduous"
colnames(ddd2)[colnames(ddd2) == "Mixwood"] <- "Mixedwood"
colnames(ddd2)[colnames(ddd2) == "BSpr"] <- "BlackSpruce"
colnames(ddd2)[colnames(ddd2) == "Larch"] <- "LarchFen"
colnames(ddd2)[colnames(ddd2) == "Shrub"] <- "Shrubland"
colnames(ddd2)[colnames(ddd2) == "GrassHerb"] <- "Grassland"

ddd1x <- ddd1[,intersect(colnames(ddd1), colnames(ddd2))]
ddd2x <- ddd2[,intersect(colnames(ddd1), colnames(ddd2))]

m1 <- cbind(v14=100*colSums(ddd1x)/sum(ddd1),
    v15=100*colSums(ddd2x)/sum(ddd2))
barplot(m1, beside=TRUE)

m2 <- cbind(v14=100*colSums(ddd1[,setdiff(colnames(ddd1), colnames(ddd2))])/sum(ddd1),
    v15=NA)
m3 <- cbind(v14=NA,
    v15=100*colSums(ddd2[,setdiff(colnames(ddd2), colnames(ddd1))])/sum(ddd2))

m4 <- rbind(m1, m2, m3)

mm <- data.frame(
    up15=rowSums(ddd2[,c("WhiteSpruce", "Deciduous", "Mixedwood", "Pine")]),
    twet15=rowSums(ddd2[,c("BlackSpruce", "LarchFen")]),
    wet15=rowSums(ddd2[,c("WetBare", "WetGrassHerb", "WetShrub")]),
    up14=rowSums(ddd1[,c("WhiteSpruce", "Deciduous", "Mixedwood", "Pine")]),
    twet14=rowSums(ddd1[,c("BlackSpruce", "LarchFen")]),
    wet14=rowSums(ddd1[,c("Bog", "Fen", "Swamp", "Marsh")]))
summary(mm)
plot(mm)
plot(mm[mm$wet15==0,])




write.csv(m4, file="abmi-1ha-comparison-v2014-2015.csv")
write.csv(ddd1, file="v2014.csv")
write.csv(ddd2, file="v2015.csv")

