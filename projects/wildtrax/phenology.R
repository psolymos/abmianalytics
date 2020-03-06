library(mgcv)
library(mefa4)
library(pbapply)
ROOT <- "d:/abmi/AB_data_v2019"

load(file.path(ROOT, "data", "inter", "species", "birds-revisit.Rdata"))

str(xt_vis)
str(x_vis)
head(x_vis)


SPP <- c(
    "AlderFlycatcher", "AmericanBittern", "AmericanCoot", "AmericanCrow",
    "AmericanGoldfinch", "AmericanRedstart", "AmericanRobin", "AmericanThreetoedWoodpecker",
    "AmericanTreeSparrow", "AmericanWigeon", "BairdsSparrow", "BaltimoreOriole",
    "BarnSwallow", "BarredOwl", "BaybreastedWarbler",
    "BlackbilledMagpie", "BlackcappedChickadee", "BlackthroatedGreenWarbler",
    "BlackAndWhiteWarbler", "BlackTern",
    "BlackpollWarbler", "BlueheadedVireo", "BluewingedTeal", "BlueJay",
    "BonapartesGull", "BorealChickadee", "BorealChorusFrog", "BorealOwl",
    "BrewersBlackbird", "BrownheadedCowbird", "BrownCreeper", "CanadaGoose",
    "CanadaWarbler", "CanadianToad", "CapeMayWarbler", "CedarWaxwing",
    "ChestnutcollaredLongspur", "ChippingSparrow", "ClaycoloredSparrow",
    "CommonGrackle", "CommonLoon", "CommonNighthawk", "CommonRaven",
    "CommonRedpoll", "CommonYellowthroat", "ConnecticutWarbler", "DarkeyedJunco", "DownyWoodpecker",
    "EasternKingbird", "EuropeanStarling", "EveningGrosbeak",
    "FoxSparrow", "FranklinsGull", "Gadwall", "GoldencrownedKinglet",
    "GrasshopperSparrow", "GrayCatbird", "GrayJay", "GreatHornedOwl",
    "GreaterWhitefrontedGoose", "GreaterYellowlegs", "GreenwingedTeal",
    "HairyWoodpecker", "HermitThrush", "HornedLark",
    "HouseSparrow", "HouseWren", "Killdeer", "LeContesSparrow", "LeastFlycatcher",
    "LesserYellowlegs", "LincolnsSparrow", "LongbilledCurlew", "MagnoliaWarbler",
    "Mallard", "MarbledGodwit", "MarshWren", "MourningDove",
    "MourningWarbler", "NelsonsSparrow", "NorthernFlicker", "NorthernPintail",
    "NorthernSawwhetOwl", "NorthernShoveler", "NorthernWaterthrush",
    "OlivesidedFlycatcher", "OrangecrownedWarbler", "Ovenbird", "PalmWarbler",
    "PerchingBirdsSongbirds", "PhiladelphiaVireo", "PiedbilledGrebe",
    "PileatedWoodpecker", "PineSiskin", "PurpleFinch", "RedbreastedNuthatch",
    "RedeyedVireo", "RedneckedGrebe", "RedneckedPhalarope", "RedwingedBlackbird",
    "RedSquirrel", "RingbilledGull", "RingneckedPheasant", "RockPigeon",
    "RosebreastedGrosbeak", "RubycrownedKinglet", "RuffedGrouse",
    "RustyBlackbird", "SandhillCrane",
    "SavannahSparrow", "SedgeWren", "SharptailedGrouse",
    "SnowGoose", "SolitarySandpiper", "SongSparrow", "Sora",
    "SpottedSandpiper", "SpraguesPipit", "SwainsonsThrush", "SwampSparrow",
    "TennesseeWarbler", "TreeSwallow",
    "UplandSandpiper", "VariedThrush",
    "VesperSparrow", "Vireos", "WarblingVireo",
    "WesternMeadowlark", "WesternTanager", "WesternToad", "WesternWoodpewee",
    "WhitecrownedSparrow", "WhitethroatedSparrow", "WhitewingedCrossbill",
    "Willet", "WilsonsSnipe", "WilsonsWarbler", "WinterWren", "WoodFrog",
    "YellowbelliedFlycatcher", "YellowbelliedSapsucker",
    "YellowheadedBlackbird", "YellowrumpedWarbler", "YellowRail",
    "YellowWarbler")

xt_vis[xt_vis > 0] <- 1
xt_vis <- xt_vis[,SPP]

dim(xt_vis)

gamfun <- function(spp, sub=FALSE) {
    x <- data.frame(spp=as.numeric(xt_vis[,spp]),
        ToY=x_vis$ToY,
        site=x_vis$ABMISite)
    if (sub) {
        z <- aggregate(x$spp, list(site=x$site), sum)
        x <- x[x$site %in% rownames(z)[z$x > 0],]
    }

    m <- gam(spp ~ s(ToY), data=x, family="binomial")
    xn <- data.frame(ToY=seq(0, 365, 1))
    xn$p <- predict(m, newdata=xn, type="response")
    m[["_predicted"]] <- xn
    m
}

M <- pblapply(colnames(xt_vis), gamfun)
names(M) <- colnames(xt_vis)

Mx <- pblapply(colnames(xt_vis), gamfun, sub=TRUE)
names(Mx) <- colnames(xt_vis)

glmfun <- function(spp, sub=FALSE) {
    x <- data.frame(spp=as.numeric(xt_vis[,spp]),
        ToY=x_vis$ToY,
        site=x_vis$ABMISite)
    x$ToY2 <- x$ToY^2
    x$ToY3 <- x$ToY^3
    x$ToY4 <- x$ToY^4
    if (sub) {
        z <- aggregate(x$spp, list(site=x$site), sum)
        x <- x[x$site %in% rownames(z)[z$x > 0],]
    }
    m <- list(m0 = glm(spp ~ 1, data=x, family="binomial"),
        m1 = update(m0, . ~ . + ToY),
        m2 = update(m1, . ~ . + ToY2),
        m3 = update(m2, . ~ . + ToY3),
        m4 = update(m3, . ~ . + ToY4))
    m <- m[[which.min(sapply(list(m0, m1, m2, m3, m4), BIC))]]

    xn <- data.frame(ToY=seq(0, 365, 1))
    xn$ToY2 <- xn$ToY^2
    xn$ToY3 <- xn$ToY^3
    xn$ToY4 <- xn$ToY^4
    xn$p <- predict(m, newdata=xn, type="response")
    m[["_predicted"]] <- xn[,c("ToY", "p")]
    m
}

M2 <- pblapply(colnames(xt_vis), glmfun)
names(M2) <- colnames(xt_vis)
M2x <- pblapply(colnames(xt_vis), glmfun, sub=TRUE)
names(M2x) <- colnames(xt_vis)


ss <- M[[1]][["_predicted"]][,1] >= 70 & M[[1]][["_predicted"]][,1] <= 208


pdf("phenol.pdf",onefile=TRUE, width=8, height=10)
for (i in names(sort(colSums(xt_vis)))) {
    Max <- max(M[[i]][["_predicted"]][ss,2], M2[[i]][["_predicted"]][ss,2],
        Mx[[i]][["_predicted"]][ss,2], M2x[[i]][["_predicted"]][ss,2])
    op <- par(mfrow=c(2,2))
    plot(M[[i]][["_predicted"]], type="l", main=paste(i, "GAM"), sub=paste("n =", sum(xt_vis[,i])),
        col=2, ylim=c(0,Max))
    lines(M[[i]][["_predicted"]][ss,])
    rug(jitter(x_vis$ToY[xt_vis[,i] > 0]))
    plot(M2[[i]][["_predicted"]], type="l", main="GLM", sub=paste("n =", sum(xt_vis[,i])),
        col=2, ylim=c(0,Max))
    lines(M2[[i]][["_predicted"]][ss,])
    rug(jitter(x_vis$ToY[xt_vis[,i] > 0]))

    plot(Mx[[i]][["_predicted"]], type="l", main=paste(i, "GAM"), sub=paste("n =", sum(xt_vis[,i])),
        col=2, ylim=c(0,Max))
    lines(Mx[[i]][["_predicted"]][ss,])
    rug(jitter(x_vis$ToY[xt_vis[,i] > 0]))
    plot(M2x[[i]][["_predicted"]], type="l", main="GLM", sub=paste("n =", sum(xt_vis[,i])),
        col=2, ylim=c(0,Max))
    lines(M2x[[i]][["_predicted"]][ss,])
    rug(jitter(x_vis$ToY[xt_vis[,i] > 0]))
    par(op)
}
dev.off()

