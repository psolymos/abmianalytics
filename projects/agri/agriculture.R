## verified info for 4 Cult classes:
## - 2013-2016 (2017?) uses RGB 1.5m res
## - 1999-2003 uses B/W 2.5m res (Crop & Tame not separated well)

library(mefa4)
load("e:/peter/AB_data_v2018/data/analysis/site/veg-hf_SiteCenter_v6verified.Rdata")
load("e:/peter/AB_data_v2018/data/inter/species/vplants_2018-05-23.Rdata")

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-V6.csv")

dd_1ha[[1]] <- dd_1ha[[1]] / rowSums(dd_1ha[[1]])
dd_qha[[1]] <- dd_qha[[1]] / rowSums(dd_qha[[1]])

compare_sets(colnames(dd_1ha[[1]]), rownames(tv))
setdiff(colnames(dd_1ha[[1]]), rownames(tv))
setdiff(rownames(tv), colnames(dd_1ha[[1]]))

tv <- tv[colnames(dd_1ha[[1]]),]

x <- groupSums(dd_1ha[[1]], 2, tv$VEGHF_FINE)
x2 <- groupSums(dd_qha[[1]], 2, tv$VEGHF_FINE)
rn <- sapply(sapply(rownames(x2), strsplit, "_"), function(z)
    paste0(z[1], "_", z[2]))
x2 <- groupSums(x2, 1, rn) / 4
ii <- intersect(rownames(x), rownames(x2))
x2 <- x2[ii,]
x <- x[ii,]

d <- data.frame(P=100*colMeans(x))
round(d, 2)
#write.csv(d, file="v6-fineHF_percent-in-1ha.csv")

pdf("v6-fineHF_percent-in-1ha.pdf", onefile=TRUE, width=12, height=8)
op <- par(mfrow=c(2,3))
for (i in colnames(x)) {
    hist(100*x[,i], xlab="Percent in 1 ha", main=i, col="tomato", xlim=c(0,100))
    hist(100*x2[,i], xlab="Percent in 1/4 ha", main="with 0s", col="gold", xlim=c(0,100))
    plot(100*x[,i], 100*x2[,i], xlim=c(0,100), ylim=c(0,100),
        xlab="Percent in 1 ha", ylab="Percent in 1/4 ha", col=2)
    abline(0,1,lty=2)

    if (sum(x[,i]>0) == 0)
        plot.new() else hist(100*x[x[,i]>0,i], xlab="Percent in 1 ha",
            main=i, col="tomato", xlim=c(0,100))
    if (sum(x2[,i]>0) == 0)
        plot.new() else hist(100*x2[x2[,i]>0,i], xlab="Percent in 1/4 ha",
            main="without 0s", col="gold", xlim=c(0,100))
    plot.new()
}
par(op)
dev.off()

## indval

meta <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
mm <- m_site_year_sub
samp(mm)$nr <- meta$NATURAL_REGIONS[match(samp(mm)$nearest, meta$SITE_ID)]
compare_sets(rownames(mm), rownames(dd_qha[[1]]))
jj <- intersect(rownames(mm), rownames(dd_qha[[1]]))

z <- dd_qha[[1]][jj,]
mm <- mm[jj,]
mx <- find_max(z)
samp(mm)$hab <- mx$index
samp(mm)$val <- mx$value

k <- c("GrassHerb", "CultivationRoughPasture", "CultivationTamePasture", "CultivationCrop")
kk <- c("Grass", "Rough", "Tame", "Cult")
mm <- mm[samp(mm)$hab %in% k & samp(mm)$val >= 0.75,]
mm <- mm[samp(mm)$nr %in% c("Grassland", "Parkland"),]
summary(samp(mm))
mm <- mm[,colSums(xtab(mm)) >= 20 & taxa(mm)$TaxonomicResolution == "Species"]
samp(mm)$hab <- factor(as.character(samp(mm)$hab), k)
levels(samp(mm)$hab) <- kk

library(opticut)

y <- as.matrix(mm)
oc <- opticut(y ~ 1, strata=samp(mm)$hab, dist="binomial")

plot(oc, sort=1)

df <- as.data.frame(summary(oc))
tx <- taxa(mm)[rownames(df),1:3]
write.csv(data.frame(tx, df), row.names=FALSE, file="cultivation-indicators.csv")
