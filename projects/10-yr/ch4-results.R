library(mefa4)
ROOT <- "v:/contents/2017/species"
TAXA <- c("mites", "mosses", "lichens", "vplants")
AllIn <- list()
vegcols <- read.csv("c:/Users/Peter/repos/abmianalytics/projects/10-yr/veg-cols.csv")

TAX <- "birds"

lt0 <- read.csv(file.path(ROOT, TAX, "lookup.csv"))
linn0 <- read.csv(file.path(ROOT, TAX, "lin10-north.csv"))
lins0 <- read.csv(file.path(ROOT, TAX, "lin10-south.csv"))
sectn0 <- read.csv(file.path(ROOT, TAX, "sector-north.csv"))
sects0 <- read.csv(file.path(ROOT, TAX, "sector-south.csv"))
soilnt0 <- read.csv(file.path(ROOT, TAX, "soilhf-nontreed-south.csv"))
soiltr0 <- read.csv(file.path(ROOT, TAX, "soilhf-treed-south.csv"))
usen0 <- read.csv(file.path(ROOT, TAX, "useavail-north.csv"))
uses0 <- read.csv(file.path(ROOT, TAX, "useavail-south.csv"))
veg0 <- read.csv(file.path(ROOT, TAX, "veghf-north.csv"))
sects0$PopEffect.Forestry <- sects0$UnitEffect.Forestry <- NULL
AllIn[["birds"]] <- list(lt=lt0, usen=usen0, uses=uses0,
    veg=veg0, linn=linn0, soilnt=soilnt0, soiltr=soiltr0,
    lins=lins0, sectn=sectn0, sects=sects0)

#TAX <- "mites"
for (TAX in TAXA) {

lt <- read.csv(file.path(ROOT, TAX, "lookup.csv"))
linn <- read.csv(file.path(ROOT, TAX, "lin10-north.csv"))
lins <- read.csv(file.path(ROOT, TAX, "lin10-south.csv"))
paspen <- read.csv(file.path(ROOT, TAX, "paspen.csv"))
sectn <- read.csv(file.path(ROOT, TAX, "sector-north.csv"))
sects <- read.csv(file.path(ROOT, TAX, "sector-south.csv"))
veg <- read.csv(file.path(ROOT, TAX, "veghf-north.csv"))
soil <- read.csv(file.path(ROOT, TAX, "soilhf-south.csv"))
usen <- read.csv(file.path(ROOT, TAX, "useavail-north.csv"))
uses <- read.csv(file.path(ROOT, TAX, "useavail-south.csv"))

## lookup
SPP <- as.character(lt$Analysis_Name)
N <- length(SPP)
fl <- gsub(".png", "", list.files(file.path(ROOT, TAX, "map-det")))
compare_sets(SPP, fl)
lt1 <- with(lt, data.frame(
    SpeciesID=Analysis_Name,
    Species=Scientific_Name_Analysis,
    CommonName=NA,
    ScientificName=Scientific_Name_Analysis,
    SppCode=SppCode,
    TSNID=TSNID,
    nSites=nSites,
    nQuadrant=Occurrences,
    map.det=fl %in% SPP,
    veghf.north=Analysis.North,
    soilhf.south=Analysis.South,
    map.pred=Maps,
    useavail.north=UseAvailability.North,
    useavail.south=UseAvailability.South,
    comments=NA
))
lt1$useavail.north[lt1$veghf.north] <- FALSE
lt1$useavail.south[lt1$soilhf.south] <- FALSE
table(mod=lt1$veghf.north, use=lt1$useavail.north)
table(mod=lt1$soilhf.south, use=lt1$useavail.south)
stopifnot(sum(diag(table(modAll=lt1$veghf.north | lt1$soilhf.south, map=lt1$map.pred))) == N)

## veghf
veg1 <- veg
veg1 <- veg1[,colnames(veg1) %in% vegcols$In]
#colnames(veg1) <- vegcols$Out[match(colnames(veg1), vegcols$In)]
veg1 <- veg1[,match(vegcols$In, colnames(veg1))]
colnames(veg1) <- vegcols$Out
compare_sets(lt1$Species, veg1$Species)

## lin N
linn1 <- data.frame(Species=veg[,"Species"], linn,
    veg[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
colnames(linn1) <- c("Species", "AverageCoef", "SoftLin10", "HardLin10", "SoftLin",
    "SoftLin.LCL", "SoftLin.UCL", "HardLin", "HardLin.LCL", "HardLin.UCL")
compare_sets(lt1$Species, linn1$Species)

## soil, nontreed and treed
soilnt1 <- soil
soilnt1$X <- NULL
colnames(soilnt1) <- gsub(".LCI", ".LCL", colnames(soilnt1))
colnames(soilnt1) <- gsub(".UCI", ".UCL", colnames(soilnt1))
soilnt1 <- soilnt1[,colnames(soilnt0)]
soiltr1 <- data.frame(Species=soilnt1$Species,
    plogis(qlogis(as.matrix(soilnt1[,-1])) + paspen$pAspen))
compare_sets(lt1$Species, soilnt1$Species)

## lin S
lins1 <- data.frame(Species=soil[,"Species"], lins,
    soil[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
colnames(lins1) <- c("Species", "AverageCoef", "SoftLin10", "HardLin10", "SoftLin",
    "SoftLin.LCL", "SoftLin.UCL", "HardLin", "HardLin.LCL", "HardLin.UCL")
compare_sets(lt1$Species, lins1$Species)

## use-avail
usen1 <- usen
usen1$Species <- usen1$SpLabel
usen1$SpLabel <- NULL
usen1 <- usen1[,!grepl("\\.WRSI", colnames(usen1))]
colnames(usen1) <- gsub("\\.rWRSI", "", colnames(usen1))
#compare_sets(lt1$Species, usen1$Species)

uses1 <- uses
uses1$Species <- uses1$SpLabel
uses1$SpLabel <- NULL
uses1 <- uses1[,!grepl("\\.WRSI", colnames(uses1))]
colnames(uses1) <- gsub("\\.rWRSI", "", colnames(uses1))

## sector
sectn1 <- sectn
sectn1 <- data.frame(Species=sectn1$Sp, sectn1)
sectn1$Sp <- sectn1$pcTotalProvAbundance <- sectn1$Reg_Model <- NULL
colnames(sectn1) <- gsub("pcTotalEffect", "PopEffect", colnames(sectn1))
colnames(sectn1) <- gsub("pcUnitEffect", "UnitEffect", colnames(sectn1))
sectn1 <- sectn1[,colnames(sectn0)]
compare_sets(lt1$Species, sectn1$Species)

## note: this does not have Forestry (Grassland!)
sects1 <- sects
sects1 <- data.frame(Species=sects1$Sp, sects1)
sects1$Sp <- sects1$pcTotalProvAbundance <- sects1$Reg_Model <- NULL
#sects1$pcTotalEffect.Forestry <- sects1$pcUnitEffect.Forestry <- 0
#sects1 <- sects1[,colnames(sectn1)]
colnames(sects1) <- gsub("pcTotalEffect", "PopEffect", colnames(sects1))
colnames(sects1) <- gsub("pcUnitEffect", "UnitEffect", colnames(sects1))
sects1 <- sects1[,colnames(sects0)]
compare_sets(lt1$Species, sects1$Species)

AllIn[[TAX]] <- list(lt=lt1, usen=usen1, uses=uses1,
    veg=veg1, linn=linn1, soilnt=soilnt1, soiltr=soiltr1,
    lins=lins1, sectn=sectn1, sects=sects1)
}
save(AllIn, file="~/Dropbox/abmi/10yr/R/AllTables.Rdata")


## making sense

library(mefa4)
load("~/Dropbox/abmi/10yr/R/AllTables.Rdata")

get_stuff0 <- function(TAX, TABLE, COL, north=TRUE) {
    SPP <- AllIn[[TAX]]$lt$Species
    i <- if (north)
        AllIn[[TAX]]$lt[,"veghf.north"] else AllIn[[TAX]]$lt[,"soilhf.south"]
    SPP <- as.character(SPP[i])
    z <- AllIn[[TAX]][[TABLE]]
    rownames(z) <- z$Species
    #compare_sets(rownames(z), SPP)
    z[SPP, COL]
}
rugged_mat <- function(x) {
    n <- sapply(x, length)
    out <- matrix(NA, max(n), length(x))
    colnames(out) <- names(x)
    for (i in 1:length(x))
        out[1:length(x[[i]]),i] <- x[[i]]
    out
}
get_stuff <- function(TABLE, COL, north=TRUE)
    rugged_mat(lapply(structure(names(AllIn), names=names(AllIn)),
        get_stuff0, TABLE=TABLE, COL=COL, north=north))
vp1 <- function(x, p=c(0, 1), ...) {
    x2 <- x[!is.na(x)]
    q <- quantile(x2, p)
    x2 <- x2[x2 > q[1] & x2 < q[2]]
    d <- density(x2, ...)
    data.frame(h=d$x, w=d$y/max(d$y))
}
vp <- function(x, p=c(0,1),
col="#FF664D80", border="#FF664D",
#col="#40E0D180", border="#40E0D1",
ylim, ylab="", xlab="", main="", ...) {
    xx <- apply(x, 2, vp1, p=p, ...)
    if (missing(ylim))
        ylim <- range(x, na.rm=TRUE)
    plot(0, type="n", axes=FALSE, ann=FALSE, xlim=c(0.5, ncol(x)+0.5),
        ylim=ylim)
    axis(1, 1:ncol(x), colnames(x), lwd=0)
    axis(2)
    title(ylab=ylab, xlab=xlab, main=main)
    for (i in 1:ncol(x)) {
        polygon(0.45*c(-xx[[i]]$w, rev(xx[[i]]$w))+i,
            c(xx[[i]]$h, rev(xx[[i]]$h)), col=col, border=border)
        #lines(c(i,i), range(xx[[i]]$h), col=col, lwd=2)
        lines(c(i-0.2, i+0.2), rep(median(x[,i], na.rm=TRUE), 2), lwd=3)
    }
    #points(1:ncol(x), colMeans(x, na.rm=TRUE), pch=21, cex=1.5)
    invisible(NULL)
}
fu <- function(x) sign(x) * plogis(log(abs(x/100)))

## --- linear

softn <- get_stuff("linn", "SoftLin10", TRUE) /
    get_stuff("linn", "AverageCoef", TRUE)
hardn <- get_stuff("linn", "HardLin10", TRUE) /
    get_stuff("linn", "AverageCoef", TRUE)
softs <- get_stuff("lins", "SoftLin10", FALSE) /
    get_stuff("lins", "AverageCoef", FALSE)
hards <- get_stuff("lins", "HardLin10", FALSE) /
    get_stuff("lins", "AverageCoef", FALSE)
softn[softn > 10] <- NA
softs[softs > 10] <- NA
hardn[hardn > 10] <- NA
hards[hards > 10] <- NA

## might need to calculate mean density for birds instead of current AvgCoefficient
par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(0, 3)
p <- c(0.025, 0.975)
vp(softn, main="Soft Linear, North", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)
vp(hardn, main="Hard Linear, North", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)
vp(softs, main="Soft Linear, South", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)
vp(hards, main="Hard Linear, South", ylim=ylim, p=p, ylab="Std. Effect")
abline(h=1, lty=2)

## --- sector

tnagr <- get_stuff("sectn", "PopEffect.Agriculture", TRUE)
tnfor <- get_stuff("sectn", "PopEffect.Forestry", TRUE)
tneng <- get_stuff("sectn", "PopEffect.Energy", TRUE)
tnurb <- get_stuff("sectn", "PopEffect.RuralUrban", TRUE)
tntra <- get_stuff("sectn", "PopEffect.Transportation", TRUE)

tsagr <- get_stuff("sects", "PopEffect.Agriculture", FALSE)
tsfor <- get_stuff("sects", "PopEffect.Forestry", FALSE)
tseng <- get_stuff("sects", "PopEffect.Energy", FALSE)
tsurb <- get_stuff("sects", "PopEffect.RuralUrban", FALSE)
tstra <- get_stuff("sects", "PopEffect.Transportation", FALSE)

unagr <- get_stuff("sectn", "UnitEffect.Agriculture", TRUE)
unfor <- get_stuff("sectn", "UnitEffect.Forestry", TRUE)
uneng <- get_stuff("sectn", "UnitEffect.Energy", TRUE)
unurb <- get_stuff("sectn", "UnitEffect.RuralUrban", TRUE)
untra <- get_stuff("sectn", "UnitEffect.Transportation", TRUE)

usagr <- get_stuff("sects", "UnitEffect.Agriculture", FALSE)
usfor <- get_stuff("sects", "UnitEffect.Forestry", FALSE)
useng <- get_stuff("sects", "UnitEffect.Energy", FALSE)
usurb <- get_stuff("sects", "UnitEffect.RuralUrban", FALSE)
ustra <- get_stuff("sects", "UnitEffect.Transportation", FALSE)

par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(-1, 1)
p <- c(0, 1)
vp(fu(uneng), main="Energy, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unagr), main="Algiculture, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unurb), main="Rural/Urban, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(untra), main="Transportation, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(unfor), main="Forestry, North", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
plot.new()

par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(-1, 1)
p <- c(0, 1)
vp(fu(useng), main="Energy, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(usagr), main="Algiculture, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(usurb), main="Rural/Urban, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")
vp(fu(ustra), main="Transportation, South", ylim=ylim, p=p, ylab="Unit Effect")
abline(h=0, lty=2)
abline(h=c(-0.5,0.5), lty=3, col="darkgrey")


par(mfrow=c(3,2), las=1, yaxs="i")
ylim <- c(-10, 10)
p <- c(0.025, 0.975)
vp(tneng, main="Energy, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tnagr, main="Algiculture, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tnurb, main="Rural/Urban, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tntra, main="Transportation, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tnfor, main="Forestry, North", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
plot.new()

par(mfrow=c(2,2), las=1, yaxs="i")
ylim <- c(-20, 20)
p <- c(0.025, 0.975)
vp(tseng, main="Energy, South", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tsagr, main="Algiculture, South", ylim=c(-100, 100), p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tsurb, main="Rural/Urban, South", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)
vp(tstra, main="Transportation, South", ylim=ylim, p=p, ylab="Total Effect")
abline(h=0, lty=2)

## --- ordination

ord_fun <- function(TAX, north=TRUE, scaling=2, ...) {
    Excl <- c("WhiteSpruceCC20", "WhiteSpruceCC40", "WhiteSpruceCC60",
        "PineCC20", "PineCC40", "PineCC60",
        "DeciduousCC20",   "DeciduousCC40", "DeciduousCC60",
        "MixedwoodCC20", "MixedwoodCC40",   "MixedwoodCC60")
    z <- if (north)
        AllIn[[TAX]]$veg else AllIn[[TAX]]$soilnt
    rownames(z) <- z[,1]
    z <- as.matrix(z[,-1])
    z <- z[,!grepl("\\.", colnames(z))]
    z <- z[,!(colnames(z) %in% Excl)]
    cn <- colnames(z)
    cn <- gsub("[[:digit:]]", "", cn)
    #cn <- gsub("CC", "", cn)
    zz <- groupMeans(z, 2, cn)
    rownames(zz) <- make.cepnames(rownames(zz))
    colnames(zz) <- make.cepnames(colnames(zz))

    library(vegan)
    fs <- function(x) 0.5*(1+x/max(abs(x)))
    yy <- (t(zz))
    m <- cca(yy)
    s <- scores(m, 1:3)
    ColSi <- 2
    ColSp <- rgb(red=fs(s$species[,1]),
        green=fs(s$species[,2]), blue=fs(s$species[,3]))
    plot(m, scaling=scaling, type="none", ...)
    #text(m, "species", col=ColSp, cex=0.5, scaling=3)
    points(m, "species", col=ColSp, cex=1, pch=19, scaling=3)
    text(m, "sites", col=1, cex=1, scaling=3)
    invisible(NULL)
}

par(mfrow=c(2,2), las=1)
ord_fun("birds", TRUE, main="Birds, North")
ord_fun("birds", FALSE, main="Birds, South")
ord_fun("mites", TRUE, main="Mites, North")
ord_fun("mites", FALSE, main="Mites, South")

par(mfrow=c(2,2), las=1)
ord_fun("lichens", TRUE, main="Lichens, North")
ord_fun("lichens", FALSE, main="Lichens, South")
ord_fun("mosses", TRUE, main="Bryophytes, North")
ord_fun("mosses", FALSE, main="Bryophytes, South")

par(mfrow=c(1,2), las=1)
ord_fun("vplants", TRUE, main="Vascular Plants, North")
ord_fun("vplants", FALSE, main="Vascular Plants, South")

## --- SI distributions: violplots, Prov/N/S  by taxon

## --- SI vs HF: gam/loess splines Prov/N/S (thf, alien, succ by taxon)


## OCCC: 10km Rdata file Coef.bs (R object Mite intactness North by 10x10km unit BS 2017 10Km2.Rdata)

