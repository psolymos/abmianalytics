library(XLConnect)
library(mefa4)

fl <- c(#mammals="e:/peter/sppweb2016/round01/tables/ABMI-species-v3.2_mammals.xlsx",
    birds="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_birds.xlsx",
    lichens="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_lichens.xlsx",
    mites="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_mites.xlsx",
    mosses="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_mosses.xlsx",
    vplants="e:/peter/sppweb2016/round01/tables/ABMI-species-v4.0_vplants.xlsx")

upfor <- c("WhiteSpruce0", "WhiteSpruce10", "WhiteSpruce20",
    "WhiteSpruce40", "WhiteSpruce60", "WhiteSpruce80", "WhiteSpruce100",
    "WhiteSpruce120", "WhiteSpruce140", "Pine0", "Pine10", "Pine20",
    "Pine40", "Pine60", "Pine80", "Pine100", "Pine120", "Pine140",
    "Deciduous0", "Deciduous10", "Deciduous20", "Deciduous40", "Deciduous60",
    "Deciduous80", "Deciduous100", "Deciduous120", "Deciduous140",
    "Mixedwood0", "Mixedwood10", "Mixedwood20", "Mixedwood40", "Mixedwood60",
    "Mixedwood80", "Mixedwood100", "Mixedwood120", "Mixedwood140")
lowfor <- c("BlackSpruce0", "BlackSpruce10", "BlackSpruce20", "BlackSpruce40",
    "BlackSpruce60", "BlackSpruce80", "BlackSpruce100", "BlackSpruce120",
    "BlackSpruce140", "Larch0", "Larch10", "Larch20", "Larch40",
    "Larch60", "Larch80", "Larch100", "Larch120", "Larch140")
open <- c("Swamp", "WetGrass", "WetShrub", "Shrub", "GrassHerb")

all <- list()
for (i in names(fl)) {
    cat(i, "\n");flush.console()
    wb <- loadWorkbook(fl[i])
    d0 <- readWorksheet(wb, sheet = "VegetationNorth", header = TRUE)
    d <- d0[,!grepl("[[:punct:]]", colnames(d0))]
    rownames(d0) <- d0$Species
    rownames(d) <- d$Species
    dc <- d[,grepl("WhiteSpruce", colnames(d))]
    dd <- d[,grepl("Deciduous", colnames(d))]
    dcf <- dc[,!grepl("CC", colnames(dc))]
    dcc <- dc[,grepl("CC", colnames(dc))]
    ddf <- dd[,!grepl("CC", colnames(dd))]
    ddc <- dd[,grepl("CC", colnames(dd))]
    all[[i]] <- list(hab=d0, conif_fire=dcf, decid_fire=ddf,
        conif_cc=dcc, decid_cc=ddc, spp=data.frame(Species=d[,1]))
}

## establish age preference

for (i in names(fl)) {
    maxc <- find_max(all[[i]]$conif_fire)
    maxd <- find_max(all[[i]]$decid_fire)
    maxc$max <- as.integer(gsub("[[:alpha:]]", "", as.character(maxc$index)))
    all[[i]]$spp$esc <- maxc$max < 60
    maxd$max <- as.integer(gsub("[[:alpha:]]", "", as.character(maxd$index)))
    all[[i]]$spp$esd <- maxd$max < 60
    all[[i]]$spp$upfor <- rowMeans(all[[i]]$hab[,colnames(all[[i]]$hab) %in% upfor], na.rm=TRUE)
    all[[i]]$spp$lowfor <- rowMeans(all[[i]]$hab[,colnames(all[[i]]$hab) %in% lowfor], na.rm=TRUE)
    all[[i]]$spp$open <- rowMeans(all[[i]]$hab[,colnames(all[[i]]$hab) %in% open], na.rm=TRUE)
    tmp <- find_max(all[[i]]$spp[,c("upfor","lowfor","open")])
    all[[i]]$spp$class <- tmp$index
    all[[i]]$spp$keep <- all[[i]]$spp$esc & all[[i]]$spp$esd & all[[i]]$spp$class == "upfor"
}

for (i in names(fl)) {
    tmp <- cbind(
        cf=all[[i]]$conif_fire[,"WhiteSpruce0"],
        cc=all[[i]]$conif_cc[,"WhiteSpruceCC0"],
        df=all[[i]]$decid_fire[,"Deciduous0"],
        dc=all[[i]]$decid_cc[,"DeciduousCC0"])
    tmp <- tmp[all[[i]]$spp$keep,,drop=FALSE]
    all[[i]]$yr0 <- data.frame(taxa=i,
        spp=droplevels(all[[i]]$spp$Species[all[[i]]$spp$keep]), tmp)
}

dat <- do.call(rbind, lapply(all, "[[", "yr0"))
dat$rc <- tanh(0.5 * log(dat$cc / dat$cf))
dat$rd <- tanh(0.5 * log(dat$dc / dat$df))

par(mfrow=c(1,2))
boxplot(rc ~ taxa, dat, ylab="Association (+likes, -avoids CC)", main="Upland Conifer", col="pink")
abline(h=0, col=4)
boxplot(rd ~ taxa, dat, ylab="Association (+likes, -avoids CC)", main="Deciduous", col="pink")
abline(h=0, col=4)



