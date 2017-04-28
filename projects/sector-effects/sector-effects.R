library(opticut)
library(XLConnect)
library(mefa4)
library(vegan)

fl <- c(#mammals="e:/peter/sppweb2016/round01/tables/ABMI-species-v3.2_mammals.xlsx",
    birds="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_birds.xlsx",
    lichens="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_lichens.xlsx",
    mites="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_mites.xlsx",
    mosses="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_mosses.xlsx",
    vplants="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_vplants.xlsx")
cn <- c("WhiteSpruce0", "WhiteSpruce10", "WhiteSpruce20",
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
    "Larch60", "Larch80", "Larch100", "Larch120", "Larch140", "WhiteSpruceCC0",
    "WhiteSpruceCC10", "WhiteSpruceCC20", "WhiteSpruceCC40", "WhiteSpruceCC60",
    "PineCC0", "PineCC10", "PineCC20", "PineCC40", "PineCC60", "DeciduousCC0",
    "DeciduousCC10", "DeciduousCC20", "DeciduousCC40", "DeciduousCC60",
    "MixedwoodCC0", "MixedwoodCC10", "MixedwoodCC20", "MixedwoodCC40",
    "MixedwoodCC60", "Swamp", "WetGrass", "WetShrub", "Shrub", "GrassHerb",
    "Cult", "UrbInd")
cn2 <- c("C", "C", "C",
    "C", "C", "C", "C",
    "C", "C", "C", "C", "C",
    "C", "C", "C", "C", "C", "C",
    "D", "D", "D", "D", "D",
    "D", "D", "D", "D",
    "D", "D", "D", "D", "D",
    "D", "D", "D", "D",
    "B", "B", "B", "B",
    "B", "B", "B", "B",
    "B", "B", "B", "B", "B",
    "B", "B", "B", "B", "B", "C",
    "C", "C", "C", "C",
    "C", "C", "C", "C", "C", "D",
    "D", "D", "D", "D",
    "D", "D", "D", "D",
    "D", "W", "W", "W", "O", "O",
    "H", "H")

load("e:/peter/AB_data_v2016/out/abmi_onoff/veg-hf-clim-reg_abmi-onoff_siteCentre_incl2015.Rdata")

vt <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
compare_sets(vt$UseInAnalysisCoef, cn)
all(colnames(dd1ha[[1]])==rownames(vt))

aveg <- groupSums(dd1ha[[1]], 2, as.character(vt$UseInAnalysisCoef))
pveg <- as.matrix(aveg)
pveg <- pveg[climSite$NRNAME %in% c("Boreal", "Canadian Shield", "Foothills", "Parkland"),cn]
pveg <- pveg / rowSums(pveg)
pveg <- pveg[rowSums(is.na(pveg))==0,]
mveg <- find_max(pveg)
vveg <- as.character(mveg$index)
nn <- colMeans(pveg)
pveg2 <- pveg
pveg2[] <- 0
for (ii in 1:nrow(mveg))
    pveg2[ii, vveg[ii]] <- 1

#Population effect of agriculture sector (% of study area = 50.75%)
#Population effect of forestry sector (% of study area = 0.71%)
#Population effect of energy sector (% of study area = 1.99%)
#Population effect of rural/urban sector (% of study area = 1.88%)
#Population effect of transportation sector (% of study area = 2.55%)

veg <- list()
sef <- list()
for (i in names(fl)) {
    cat(i, "\n");flush.console()
    wb <- loadWorkbook(fl[i])
    veg[[i]] <- readWorksheet(wb, sheet = "VegetationNorth", header = TRUE)
    rownames(veg[[i]]) <- veg[[i]]$Species
    veg[[i]] <- veg[[i]][,cn]
    sef[[i]] <- readWorksheet(wb, sheet = "SectorNorth", header = TRUE)
    rownames(sef[[i]]) <- sef[[i]]$Species
    sef[[i]]$Species <- NULL
    sef[[i]] <- sef[[i]][rownames(veg[[i]]),]
}

pred <- list()
pred2 <- list()
lc <- list()
lcmat <- list()
for (i in names(fl)) {
    pred[[i]] <- apply(veg[[i]], 1, function(z) colSums(t(pveg) * z, na.rm=TRUE))
    pred2[[i]] <- apply(veg[[i]], 1, function(z) colSums(t(pveg2) * z, na.rm=TRUE))
    lc[[i]] <- t(sapply(1:nrow(veg[[i]]), function(j)
        unclass(summary(lorenz(as.numeric(ifelse(is.na(veg[[i]][j,]), 0, veg[[i]][j,])), nn)))))
    lcmat[[i]] <- ifelse(ifelse(is.na(as.matrix(veg[[i]])), 0,
        as.matrix(veg[[i]])) > lc[[i]][,"x[t]"], 1, 0)
}

nn <- table(mveg$index)/length(mveg$index)
mu <- aggregate(x, list(g), mean)

## use lc on predicted stuff at 1ha or pc level to

COL <- rep("grey", length(cn))
names(COL) <- cn
COL[c("WhiteSpruce0", "WhiteSpruce10", "WhiteSpruce20",
    "WhiteSpruce40", "WhiteSpruce60", "WhiteSpruce80", "WhiteSpruce100",
    "WhiteSpruce120", "WhiteSpruce140", "Pine0", "Pine10", "Pine20",
    "Pine40", "Pine60", "Pine80", "Pine100", "Pine120", "Pine140")] <- "green1"
COL[c("WhiteSpruceCC0",
    "WhiteSpruceCC10", "WhiteSpruceCC20", "WhiteSpruceCC40", "WhiteSpruceCC60",
    "PineCC0", "PineCC10", "PineCC20", "PineCC40", "PineCC60")] <- "turquoise1"
COL[c("Deciduous0", "Deciduous10", "Deciduous20", "Deciduous40", "Deciduous60",
    "Deciduous80", "Deciduous100", "Deciduous120", "Deciduous140",
    "Mixedwood0", "Mixedwood10", "Mixedwood20", "Mixedwood40", "Mixedwood60",
    "Mixedwood80", "Mixedwood100", "Mixedwood120", "Mixedwood140")] <- "green4"
COL[c("DeciduousCC0",
    "DeciduousCC10", "DeciduousCC20", "DeciduousCC40", "DeciduousCC60",
    "MixedwoodCC0", "MixedwoodCC10", "MixedwoodCC20", "MixedwoodCC40",
    "MixedwoodCC60")] <- "turquoise4"
COL[c("BlackSpruce0", "BlackSpruce10", "BlackSpruce20", "BlackSpruce40",
    "BlackSpruce60", "BlackSpruce80", "BlackSpruce100", "BlackSpruce120",
    "BlackSpruce140", "Larch0", "Larch10", "Larch20", "Larch40",
    "Larch60", "Larch80", "Larch100", "Larch120", "Larch140", "WhiteSpruceCC0")] <- "darkblue"
COL[c("Swamp", "WetGrass", "WetShrub")] <- "cyan"
COL[c("Shrub", "GrassHerb")] <- "orange"
COL[c("Cult", "UrbInd")] <- "red"

COL2 <- c(C="green1", D="green4", B="darkblue", W="cyan", O="orange", H="red")

k <- 1
obj <- t(lcmat[[k]])
obj2 <- t(groupMeans(obj, 1, cn2))
wm <- find_max(obj2)
#obj <- t(veg[[k]])
se <- t(as.matrix(sef[[k]])[,6:10])
rownames(se) <- c("Agr", "For", "En", "Urb", "Tran")


rownames(obj) <- make.cepnames(rownames(obj))
colnames(obj) <- make.cepnames(colnames(obj))

o <- rda(obj)
o2 <- rda(se)
s <- scores(o)
s2 <- scores(o2)
o3 <- rda(t(obj), obj2)
o4 <- rda(t(obj), t(se))

s[[1]][,1] <- s[[1]][,1]/max(abs(s[[1]][,1]))
s[[1]][,2] <- s[[1]][,2]/max(abs(s[[1]][,2]))
s[[2]][,1] <- s[[2]][,1]/max(abs(s[[2]][,1]))
s[[2]][,2] <- s[[2]][,2]/max(abs(s[[2]][,2]))
s2[[1]][,1] <- s2[[1]][,1]/max(abs(s2[[1]][,1]))
s2[[1]][,2] <- s2[[1]][,2]/max(abs(s2[[1]][,2]))
s2[[2]][,1] <- s2[[2]][,1]/max(abs(s2[[2]][,1]))
s2[[2]][,2] <- s2[[2]][,2]/max(abs(s2[[2]][,2]))

par(mfrow=c(1,3))
plot(0, type="n", ylim=c(-1,1), xlim=c(-1,1), ann=FALSE)
title(main="Land cover classes")
abline(h=0,v=0,lty=3,col="lightgrey")
text(s$sites, labels=rownames(s$sites), col=COL, cex=0.6)

plot(0, type="n", ylim=c(-1,1), xlim=c(-1,1), ann=FALSE)
title(main="Species")
abline(h=0,v=0,lty=3,col="lightgrey")
text(s$species, labels=rownames(s$species), col=COL2[as.character(wm$index)], cex=0.6)

plot(0, type="n", ylim=c(-1,1), xlim=c(-1,1), ann=FALSE)
title(main="Setctors")
abline(h=0,v=0,lty=3,col="lightgrey")
text(s2$sites, labels=rownames(s2$sites), col=2, cex=1)
text(s2$species, labels=rownames(s2$species), col=COL2[as.character(wm$index)], cex=0.6)


plot(o4, scaling="sites")
text(scores(o4)$sites, labels=rownames(scores(o4)$sites), col=COL, cex=0.6)

plot(o4, scaling="none")
text(scores(o4)$species, labels=rownames(scores(o4)$species),
    col=COL2[as.character(wm$index)], cex=0.6)

lcm <- t(obj)
sem <- tanh(t(se)/100)
m2 <- rda(lcm ~ ., data.frame(sem))
plot(m2, scaling=3, type="none")
text(m2, "species", col=COL2[as.character(wm$index)], cex=0.4, scaling=3)
text(m2, "sites", col=COL, cex=0.6, scaling=3)
text(m2, display="bp", col=4, cex=0.8, scaling=3)

par(mfrow=c(1,2))
plot(m2, scaling=3, type="none")
text(m2, "species", col=COL, cex=0.4, scaling=3)
text(m2, display="bp", col=4, cex=0.8, scaling=3)

plot(m2, scaling=3, type="none")
text(m2, "sites", col=COL2[as.character(wm$index)], cex=0.6, scaling=3)
text(m2, display="bp", col=4, cex=0.8, scaling=3)
