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
cn1 <- c("ConR", "ConR", "ConY",
    "ConY", "ConY", "ConO", "ConO",
    "ConO", "ConO", "ConR", "ConR", "ConY",
    "ConY", "ConY", "ConO", "ConO",
    "ConO", "ConO",
    "DecR", "DecR", "DecY", "DecY", "DecY",
    "DecO", "DecO", "DecO", "DecO",
    "DecR", "DecR", "DecY", "DecY", "DecY",
    "DecO", "DecO", "DecO", "DecO",
    "BSpR", "BSpR", "BSpY", "BSpY",
    "BSpY", "BSpO", "BSpO", "BSpO",
    "BSpO", "BSpR", "BSpR", "BSpY", "BSpY",
    "BSpY", "BSpO", "BSpO", "BSpO", "BSpO", "CCConR",
    "CCConR", "CCConY", "CCConY", "CCConY",
    "CCConR", "CCConR", "CCConY", "CCConY", "CCConY", "CCDecR",
    "CCDecR", "CCDecY", "CCDecY", "CCDecY",
    "CCDecR", "CCDecR", "CCDecY", "CCDecY",
    "CCDecY", "Wet", "Wet", "Wet", "Shrub", "Grass",
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
save(sef, veg, vt, aveg, cn, cn2, fl, climSite,
    file="~/Dropbox/abmi/10yr/sector/sector.Rdata")

load("~/Dropbox/abmi/10yr/sector/sector.Rdata")

pveg3 <- as.matrix(aveg)
pveg3 <- pveg3[climSite$NRNAME %in% c("Boreal", "Canadian Shield", "Foothills", "Parkland"),cn]
pveg3 <- pveg3 / rowSums(pveg3)
pveg3 <- pveg3[rowSums(is.na(pveg3))==0,]
pveg3 <- groupSums(pveg3, 2, cn1)

pveg <- as.matrix(aveg)
pveg <- pveg[climSite$NRNAME %in% c("Boreal", "Canadian Shield", "Foothills", "Parkland"),cn]
pveg <- pveg / rowSums(pveg)
pveg <- pveg[rowSums(is.na(pveg))==0,]
#pveg3 <- pveg3[rownames(pveg),]
all(rownames(pveg)==rownames(pveg3))
mveg <- find_max(pveg)
vveg <- as.character(mveg$index)
nn <- colMeans(pveg)
pveg2 <- pveg
pveg2[] <- 0
for (ii in 1:nrow(mveg))
    pveg2[ii, vveg[ii]] <- 1
mveg2 <- find_max(pveg3)
vveg2 <- as.character(mveg2$index)
nn2 <- colMeans(pveg3)
pveg4 <- pveg3
pveg4[] <- 0
for (ii in 1:nrow(mveg))
    pveg4[ii, vveg2[ii]] <- 1

#Population effect of agriculture sector (% of study area = 50.75%)
#Population effect of forestry sector (% of study area = 0.71%)
#Population effect of energy sector (% of study area = 1.99%)
#Population effect of rural/urban sector (% of study area = 1.88%)
#Population effect of transportation sector (% of study area = 2.55%)


pred <- list()
pred2 <- list()
lc <- list()
lcmat <- list()
lcmat2 <- list()
for (i in names(fl)) {
    pred[[i]] <- apply(veg[[i]], 1, function(z) colSums(t(pveg) * z, na.rm=TRUE))
    pred2[[i]] <- apply(veg[[i]], 1, function(z) colSums(t(pveg2) * z, na.rm=TRUE))
    lc[[i]] <- t(sapply(1:nrow(veg[[i]]), function(j)
        unclass(summary(lorenz(as.numeric(ifelse(is.na(veg[[i]][j,]), 0, veg[[i]][j,])), nn)))))
    lcmat[[i]] <- ifelse(ifelse(is.na(as.matrix(veg[[i]])), 0,
        as.matrix(veg[[i]])) > lc[[i]][,"x[t]"], 1, 0)
    lcmat2[[i]] <- t(sapply(1:nrow(veg[[i]]), function(j) {
        tmp <- as.numeric(ifelse(is.na(veg[[i]][j,]), 0, veg[[i]][j,]))
        tmp <- sum_by(tmp, cn1)
        lc <- unclass(summary(lorenz(tmp[,"x"], nn2)))
        ifelse(tmp[,"x"] > lc["x[t]"], 1, 0)
    }))
    rownames(lcmat2[[i]]) <- rownames(lcmat[[i]])
}

Pairs <- function (x, ...)
{
    #y <- data.frame(x)
    y <- x
    fun.lower <- function(x1, x2, ...) {
        COR <- cor(x1, x2)
        text(mean(range(x1, na.rm = TRUE)), mean(range(x2, na.rm = TRUE)),
            round(COR, 2), cex = 0.5 + 2*abs(COR))
        box(col="grey")
    }
    fun.upper <- function(x1, x2, ...) {
        if (is.factor(x1)) {
            x1 <- as.integer(x1)
        }
        if (is.factor(x2)) {
            x1 <- as.integer(x2)
        }
        abline(h=0, v=0, lty=1, col="grey")
        points(x1, x2, col=ColSp)
        #LM <- lm(x2 ~ x1)
        #abline(LM, col="#D7191C")
        box(col="#80B9D8")
    }
    panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$density
        Max <- max(y)
        y <- y/Max
        rect(breaks[-nB], 0, breaks[-1], y, col = "#FDC980", border = "#F07C4A",
            ...)
        abline(v=0, lty=1, col="grey")
        den <- density(x)
        den$y <- den$y/Max
        lines(den, col="#F07C4A")
        box(col="#F07C4A")
    }
    pairs.default(y, lower.panel = fun.lower, upper.panel = fun.upper,
        diag.panel = panel.hist)
    invisible(NULL)
}
fs <- function(x) 0.5*(1+x/max(abs(x)))

k <- 1
for (k in 1:5) {
NAM <- names(fl)[k]
obj <- t(lcmat2[[k]])
rownames(obj) <- make.cepnames(rownames(obj))
#colnames(obj) <- make.cepnames(colnames(obj))
se <- t(as.matrix(sef[[k]])[,6:10]) # unit
#se <- t(as.matrix(sef[[k]])[,1:5]) # total
rownames(se) <- c("Agr", "For", "En", "Urb", "Tran")

lcm <- t(obj)
sem <- tanh(t(se)/100)
m <- rda(lcm ~ ., data.frame(sem))
s <- scores(m, 1:3)
ColSi <- 2
ColSp <- rgb(red=fs(s$sites[,1]),
    green=fs(s$sites[,2]), blue=fs(s$sites[,3]))

pdf(paste0("~/Dropbox/abmi/10yr/sector/triplot-", NAM, ".pdf"))
plot(m, scaling=3, type="none", main=NAM)
text(m, "sites", col=ColSp, cex=0.5, scaling=3)
#points(m, "sites", col=ColSp, cex=0.5, scaling=3)
text(m, "species", col=1, cex=1, scaling=3)
text(m, display="bp", col=4, cex=0.8, scaling=3)
dev.off()

pdf(paste0("~/Dropbox/abmi/10yr/sector/seff-", NAM, ".pdf"))
Pairs(data.frame(sem))
dev.off()
}

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

sty <- apply(lcmat2[[1]], 1, function(z)
    paste(ifelse(z==1, colnames(lcmat2[[1]]), "-"), collapse=""))
data.frame(table(sty))
h <- grepl("H", sty)
COLS <- rep("grey", length(sty))
## D
COLS[!h & grepl("D", sty)] <- "green1"
## C & CD
COLS[!h & (grepl("C", sty) | grepl("CD", sty))] <- "green4"
## H+
COLS[h] <- "red"
## no H & W+
COLS[!h & grepl("W", sty)] <- "cyan"
## no H & O+
COLS[!h & grepl("O", sty)] <- "orange"
## no H & B+ but not CDB
COLS[!h & grepl("B", sty)] <- "darkblue"
## CDB
COLS[!h & grepl("CDB", sty)] <- "green"


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
text(m2, "species", col=COL, cex=0.4, scaling=3)
text(m2, "sites", col=COLS, cex=0.6, scaling=3)
#text(m2, "sites", col=ct, cex=0.6, scaling=3)
text(m2, display="bp", col=4, cex=0.8, scaling=3)

par(mfrow=c(1,2))
plot(m2, scaling=3, type="none")
text(m2, "species", col=COL, cex=0.4, scaling=3)
text(m2, display="bp", col=4, cex=0.8, scaling=3)

plot(m2, scaling=3, type="none")
text(m2, "sites", col=COLS, cex=0.6, scaling=3)
text(m2, display="bp", col=4, cex=0.8, scaling=3)
