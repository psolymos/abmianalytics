library(parallel)
library(mgcv)
library(mefa4)
library(pbapply)
ROOT <- "d:/abmi/AB_data_v2019"

load(file.path(ROOT, "data", "inter", "species", "birds-revisit.Rdata"))

str(xt_vis)
str(x_vis)
head(x_vis)

xy <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
z <- x_vis$ABMISite
z <- gsub("B", "", z)
z <- strsplit(z, "-")
z <- sapply(z, function(a) if (length(a) > 1) a[3] else a[1])
z <- as.integer(z)
#sort(unique(z))
#range(z)
x_vis <- data.frame(x_vis, xy[match(z, xy$SITE_ID),])

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
    "CommonYellowthroat", "ConnecticutWarbler", "DarkeyedJunco", "DownyWoodpecker",
    "EasternKingbird", "EuropeanStarling", "EveningGrosbeak",
    "FoxSparrow", "FranklinsGull", "Gadwall", "GoldencrownedKinglet",
    "GrasshopperSparrow", "GrayCatbird", "GrayJay", "GreatHornedOwl",
    "GreaterYellowlegs", "GreenwingedTeal",
    "HairyWoodpecker", "HermitThrush", "HornedLark",
    "HouseSparrow", "HouseWren", "Killdeer", "LeContesSparrow", "LeastFlycatcher",
    "LesserYellowlegs", "LincolnsSparrow", "LongbilledCurlew", "MagnoliaWarbler",
    "Mallard", "MarbledGodwit", "MarshWren", "MourningDove",
    "MourningWarbler", "NelsonsSparrow", "NorthernFlicker", "NorthernPintail",
    "NorthernSawwhetOwl", "NorthernShoveler", "NorthernWaterthrush",
    "OlivesidedFlycatcher", "OrangecrownedWarbler", "Ovenbird", "PalmWarbler",
    "PhiladelphiaVireo", "PiedbilledGrebe",
    "PileatedWoodpecker", "PineSiskin", "PurpleFinch", "RedbreastedNuthatch",
    "RedeyedVireo", "RedneckedGrebe", "RedneckedPhalarope", "RedwingedBlackbird",
    "RedSquirrel", "RingbilledGull", "RingneckedPheasant", "RockPigeon",
    "RosebreastedGrosbeak", "RubycrownedKinglet", "RuffedGrouse",
    "RustyBlackbird", "SandhillCrane",
    "SavannahSparrow", "SedgeWren", "SharptailedGrouse",
    "SolitarySandpiper", "SongSparrow", "Sora",
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

#sort(colSums(xt_vis[x_vis$ToY > 150,]))

xt_vis[xt_vis > 0] <- 1
xt_vis <- xt_vis[,SPP]

dim(xt_vis)

cl <- makeCluster(8)
clusterEvalQ(cl, library(mefa4))
clusterEvalQ(cl, library(mgcv))
clusterExport(cl, c("xt_vis", "x_vis"))

gamfun <- function(spp, sub=FALSE) {
    x <- data.frame(spp=as.numeric(xt_vis[,spp]),
        ToY=x_vis$ToY,
        site=x_vis$ABMISite)
    if (sub) {
        j <- x$ToY > 150
        z <- aggregate(x$spp[j], list(site=x$site[j]), sum)
        x <- x[x$site %in% rownames(z)[z$x > 0],]
    }

    m <- gam(spp ~ s(ToY), data=x, family="binomial")
    xn <- data.frame(ToY=seq(70, 208, 1))
    xn$p <- predict(m, newdata=xn, type="response")
    m[["_predicted"]] <- xn
    m
}

glmfun <- function(spp, sub=FALSE) {
    x <- data.frame(spp=as.numeric(xt_vis[,spp]),
        ToY=x_vis$ToY,
        site=x_vis$ABMISite)
    x$ToY2 <- x$ToY^2
    x$ToY3 <- x$ToY^3
    x$ToY4 <- x$ToY^4
    if (sub) {
        j <- x$ToY > 150
        z <- aggregate(x$spp[j], list(site=x$site[j]), sum)
        x <- x[x$site %in% rownames(z)[z$x > 0],]
    }
    m0 <- glm(spp ~ 1, data=x, family="binomial")
    m1 <- update(m0, . ~ . + ToY)
    m2 <- update(m1, . ~ . + ToY2)
    m3 <- update(m2, . ~ . + ToY3)
    m4 <- update(m3, . ~ . + ToY4)
    m <- list(m0=m0, m1=m1, m2=m2, m3=m3, m4=m4)
    b <- sapply(m, BIC)
    if (!m0$converged) b[1] <- Inf
    if (!m1$converged) b[2] <- Inf
    if (!m2$converged) b[3] <- Inf
    if (!m3$converged) b[4] <- Inf
    if (!m4$converged) b[5] <- Inf
    m <- m[[which.min(b)]]

    xn <- data.frame(ToY=seq(70, 208, 1))
    xn$ToY2 <- xn$ToY^2
    xn$ToY3 <- xn$ToY^3
    xn$ToY4 <- xn$ToY^4
    xn$p <- predict(m, newdata=xn, type="response")
    m[["_predicted"]] <- xn[,c("ToY", "p")]
    m
}

M <- pblapply(colnames(xt_vis), gamfun, cl=cl)
names(M) <- colnames(xt_vis)
#Mx <- pblapply(colnames(xt_vis), gamfun, sub=TRUE, cl=cl)
#names(Mx) <- colnames(xt_vis)

M2 <- pblapply(colnames(xt_vis), glmfun, cl=cl)
names(M2) <- colnames(xt_vis)
#M2x <- pblapply(colnames(xt_vis), glmfun, sub=TRUE, cl=cl)
#names(M2x) <- colnames(xt_vis)

dfun <- function(x, y, v) {
    md <- floor(max(diff(sort(x[y > 0])))*2)
    v <- max(v, md)
    sig <- v/4
    r <- range(x)
    k <- r[1]:(r[2]-v)
    xo <- 0.5 * (k + (k+v))
    yo <- xo
    yo[] <- NA
    for (i in seq_along(k)) {
        j <- which(x >= k[i] & x <= k[i]+v)
        if (length(j) & sum(y[j])>0) {
            w <- dnorm(x[j], xo[i], sig) / dnorm(xo[i], xo[i], sig)
            #yo[i] <- sum(y[j])/length(j)
            yo[i] <- sum(w*y[j])/sum(w)
        }
    }
    data.frame(x=xo, y=yo)
}

M3 <- pblapply(colnames(xt_vis), function(spp) {
    dfun(x = x_vis$ToY,
        y = as.numeric(xt_vis[,spp]),
        v = 14)
})
names(M3) <- colnames(xt_vis)

pdf("phenol.pdf",onefile=TRUE, width=6, height=4)
for (i in names(sort(colSums(xt_vis)))) {
    Max <- max(M[[i]][["_predicted"]][,2], M3[[i]][,2], na.rm=TRUE)
    op <- par(mfrow=c(1,1))
    plot(M[[i]][["_predicted"]], type="l", main=i, sub=paste("n =", sum(xt_vis[,i])),
        col=1, ylim=c(0, Max))
    lines(M2[[i]][["_predicted"]], col=2)
    lines(M3[[i]], col=4)
    rug(jitter(x_vis$ToY[xt_vis[,i] > 0]))
    legend("topleft", lty=1, col=c(1,2,4), bty="n", legend=c("GAM", "GLM", "KDE"))
    par(op)
}
dev.off()



gamfun2 <- function(spp, B=99) {
    x <- data.frame(spp=as.numeric(xt_vis[,spp]),
        ToY=x_vis$ToY)
    m <- gam(spp ~ s(ToY), data=x, family="binomial")
    xn <- data.frame(ToY=seq(70, 208, 1))
    p <- matrix(0, nrow(xn), B+1)
    j <- 1
    p[,1] <- predict(m, newdata=xn, type="response")
    while (j <= 99) {
        x1 <- x[sample(nrow(x), nrow(x), replace=TRUE),]
        if (sum(x1$spp) > 0) {
            m1 <- update(m, data=x1)
            p[,j+1] <- predict(m1, newdata=xn, type="response")
            j <- j + 1
            if (j > 200) break
        }
    }
    p <- p[,1:(j-1)]
    q <- t(apply(p, 1, quantile, c(0.5, 0.05, 0.95)))
    Mean <- rowMeans(p)
    z1 <- q[,3]-q[,2]
    z2 <- min(z1)/z1
    y <- Mean*z2+q[,2]*(1-z2)
    list(x=xn$ToY, y=y, m=Mean, p=p, q=q)
}

M4 <- pblapply(colnames(xt_vis), gamfun2, cl=cl)
names(M4) <- colnames(xt_vis)


pdf("phenol.pdf",onefile=TRUE, width=7, height=5)
for (i in names(sort(colSums(xt_vis)))) {
    z <- M4[[i]]
    matplot(z$x, z$p, type="l", lty=1, col="#00000022", main=i,
        ylim=c(0, max(z$q[z$x > 100 & z$x < 180,])))
    lines(z$x, z$m, col=4,lwd=2)
    z1 <- z$q[,3]-z$q[,2]
    z2 <- min(z1)/z1
    y <- z$q[,1]*z2+z$q[,2]*(1-z2)
    lines(spline(z$x, y, n=50), col=2,lwd=2)
    rug(jitter(x_vis$ToY[xt_vis[,i] > 0]), col="grey")
}
dev.off()

stopCluster(cl)



spp="CanadaWarbler"
x <- data.frame(spp=as.numeric(xt_vis[,spp]),
    ToY=x_vis$ToY,
    Yr=x_vis$Year,
    Yrf=as.factor(x_vis$Year),
    x=x_vis$PUBLIC_LONGITUDE,
    y=x_vis$PUBLIC_LATTITUDE)
m0 <- gam(spp ~ 1, data=x, family="binomial")
m1 <- gam(spp ~ Yr, data=x, family="binomial")
m1f <- gam(spp ~ Yrf, data=x, family="binomial")
m2 <- gam(spp ~ s(ToY), data=x, family="binomial")
m3 <- gam(spp ~ Yr + s(ToY), data=x, family="binomial")
m3f <- gam(spp ~ Yrf + s(ToY), data=x, family="binomial")

m4 <- gam(spp ~ Yrf + s(ToY) + s(x))
m5 <- gam(spp ~ Yrf + s(ToY) + s(y))
m6 <- gam(spp ~ Yrf + s(ToY) + s(x) + s(y))
m7 <- gam(spp ~ Yrf + s(ToY) + s(x, y))

AIC(m0, m1, m2, m3, m1f, m3f, m4)

AIC(m4, m5, m6, m7)

plot(m7)
points(y ~ x, x)

summary(m4)
summary(m3f)

x$ToY2 <- x$ToY^2
x$ToY3 <- x$ToY^3
x$ToY4 <- x$ToY^4
x$x2 <- x$x^2
x$x3 <- x$x^3
x$x4 <- x$x^4
x$y2 <- x$y^2
x$y3 <- x$y^3
x$y4 <- x$y^4


spp="CanadaWarbler"
pdf("phenol2.pdf",onefile=TRUE, width=7, height=5)
for (spp in SPP) {
x <- data.frame(spp=as.numeric(xt_vis[,spp]),
    ToY=x_vis$ToY,
    Yr=x_vis$Year,
    Yrf=as.factor(x_vis$Year),
    x=x_vis$PUBLIC_LONGITUDE,
    y=x_vis$PUBLIC_LATTITUDE)
x$ToY2 <- x$ToY^2
x$ToY3 <- x$ToY^3
x$ToY4 <- x$ToY^4
m <- glm(spp ~ y * (ToY + ToY2 + ToY3), x, family=binomial)
nd <- expand.grid(ToY=80:200, y=c(50, 54, 58))
nd$ToY2 <- nd$ToY^2
nd$ToY3 <- nd$ToY^3
nd$ToY4 <- nd$ToY^4
nd$p <- predict(m, newdata=nd, type="response")

plot(p ~ ToY, nd[nd$y==50,], type="l", ylim=c(0, max(nd$p)), main=spp)
lines(p ~ ToY, nd[nd$y==54,], col=2)
lines(p ~ ToY, nd[nd$y==58,], col=4)
legend("topleft", lty=1,col=c(1,2,4), legend=c("lat=50","lat=54","lat=58"))
}
dev.off()
