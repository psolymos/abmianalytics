library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"

level <- 0.9

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
#    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

## names etc
e <- new.env()
#load(file.path(ROOT, "data", "data-full-withrevisit.Rdata"), envir=e)
load(file.path(ROOT, "data", "data-wrsi.Rdata"), envir=e)
tax <- droplevels(e$TAX)
tax$Fn <- droplevels(tax$English_Name)
levels(tax$Fn) <- nameAlnum(levels(tax$Fn), capitalize="mixed", collapse="")
yy <- e$YYw
yy01 <- yy
yy01[yy01>0] <- 1

load(file.path(ROOT, "data", "data-cawa.Rdata"))
y_cawa <- YY[,"CAWA"]
off_cawa <- OFF[,"CAWA"]
rm(e, OFF, BB)

## terms and design matrices
Terms <- getTerms(mods, "list")
Xn <- model.matrix(getTerms(mods, "formula"), DAT)
colnames(Xn) <- fixNames(colnames(Xn))

names(mods)
# [1] "Hab"      "Age"      "CC"       "Contrast" "ARU"      "Space"   
# [7] "Wet"      "Dec"      "HF"       "Year" 
stage_hab <- 5
spp <- "CAWA"
NAM <- as.character(tax[spp, "English_Name"])

load(file.path(ROOT, "results", "cawa", "birds_abmi-cawa_CAWA.Rdata"))
#load(file.path(ROOT, "results", "north", "birds_abmi-north_CAWA.Rdata"))

getFancyMidTab(res, mods)
printCoefmat(getSummary(res))

est_hab <- getEst(res, stage=stage_hab, na.out=FALSE, Xn)
res_coef <- pred_veghf(est_hab, Xn, burn_included=FALSE)
MAX <- fig_veghf_ymax(res_coef)
NDAT <- sum(y_cawa > 0)
#fname <- file.path(ROOT, "figs", "veghf-north", 
#    paste0(as.character(tax[spp, "Spp"]), ".png"))
#png(file=fname,width=1500,height=700)
fig_veghf(res_coef, paste0(NAM, " (n = ", NDAT, " detections)"), ymax=MAX)
#dev.off()
## linear
#fname <- file.path(ROOT, "figs", "linear-north", 
#    paste0(as.character(tax[spp, "Spp"]), ".png"))
#png(file=fname,width=350,height=400)
fig_linear(attr(res_coef, "linear"), paste0(NAM, "\nNorth (n = ", NDAT, " det.)"))
#dev.off()

## climate & surrounding hf tables, climate surface maps

est_hf <- getEst(res, stage=length(mods)-1, na.out=FALSE, Xn)
fig_hf_noremn(est_hf, Xn, LAB=paste0(NAM, ", North"))



## -- old

cn <- c("xPET", "xMAT", "xAHM", "xFFP", 
    "xMAP", "xMWMT", "xMCMT", "xlat", "xlong", "xlat2", "xlong2", 
    "THF_KM", "Lin_KM", "Nonlin_KM", "Succ_KM", "Alien_KM", "Noncult_KM", 
    "Cult_KM", "THF2_KM", "Nonlin2_KM", "Succ2_KM", "Alien2_KM", 
    "Noncult2_KM")
transform_CLIM <- function(x, ID="PKEY") {
    z <- x[,ID,drop=FALSE]
    z$xlong <- (x$POINT_X - (-113.7)) / 2.15
    z$xlat <- (x$POINT_Y - 53.8) / 2.28
    z$xAHM <- (x$AHM - 0) / 50
    z$xPET <- (x$PET - 0) / 800
    z$xFFP <- (x$FFP - 0) / 130
    z$xMAP <- (x$MAP - 0) / 2200
    z$xMAT <- (x$MAT - 0) / 6
    z$xMCMT <- (x$MCMT - 0) / 25
    z$xMWMT <- (x$MWMT - 0) / 20
    z
}
xclim <- transform_CLIM(kgrid, "Row_Col")
xclim$xlat2 <- xclim$xlat^2
xclim$xlong2 <- xclim$xlong^2
ffTerms <- getTerms(modsn["Space"], "formula", intercept=FALSE)
Xclim <- model.matrix(ffTerms, xclim)
colnames(Xclim) <- fixNames(colnames(Xclim))
excln <- kgrid$NRNAME %in% c("Rocky Mountain", "Grassland")
excls <- rep(TRUE, nrow(kgrid))
excls[kgrid$NRNAME %in% c("Grassland", "Parkland")] <- FALSE
excls[kgrid$NSRNAME %in% c("Dry Mixedwood")] <- FALSE
clim_n <- list()
clim_s <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(tax[spp, "English_Name"])
if (tax[spp, "surroundinghf_north"]) {
    resn <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_sp <- getEst(resn, stage=stage_hab_n + 2, na.out=FALSE, Xnn)
    sp_n <- colMeans(estn_sp[,cn])
    clim_n[[spp]] <- sp_n

    fname <- file.path(ROOT, "figs", "climate-north", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    ## quick and dirty
    pr <- exp(drop(Xclim %*% colMeans(estn_sp[,colnames(Xclim)])))
    ## bootstrap based and correct
#    pr <- rowMeans(exp(apply(estn_sp[,colnames(Xclim)], 1, function(z) drop(Xclim %*% z))))
    q <- quantile(pr, 0.99)
    pr[pr > q] <- q
    pr <- pr/max(pr)
    pr[excln] <- NA
    qq <- quantile(pr, seq(0.1, 0.9, 0.1), na.rm=TRUE)
    z <- cut(pr, c(-1, unique(qq), 2))
    Col <- rev(terrain.colors(nlevels(z)))
	png(file=fname, width=600, height=1000)

    plot(kgrid$X, kgrid$Y, pch=15, cex=0.2, col=Col[z], axes=FALSE, ann=FALSE)
    points(kgrid$X[excln], kgrid$Y[excln], pch=15, cex=0.2, col="darkgrey")
    mtext(paste0(NAM, ", North"), line=2, side=3, adj=0.5, cex=1.4, col="grey40")
    points(xyw, pch=15, cex=0.2, col=rgb(0.3,0.45,0.9))
    points(city, pch=18, col="grey10")
    text(city, rownames(city), cex=0.8, adj=-0.1, col="grey10")
    legend("bottomleft", col=rev(Col), fill=rev(Col),
        legend=c("High", rep("", length(Col)-2), "Low"), bty="n")

	dev.off()
}
if (tax[spp, "surroundinghf_south"]) {
    ress <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-south_", spp, ".Rdata")))
    ests_sp <- getEst(ress, stage=stage_hab_s + 2, na.out=FALSE, Xns)
    sp_s <- colMeans(ests_sp[,cn])
    clim_s[[spp]] <- sp_s

    fname <- file.path(ROOT, "figs", "climate-south", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    ## quick and dirty
    pr <- exp(drop(Xclim %*% colMeans(ests_sp[,colnames(Xclim)])))
    ## bootstrap based and correct
#    pr <- rowMeans(exp(apply(ests_sp[,colnames(Xclim)], 1, function(z) drop(Xclim %*% z))))
    q <- quantile(pr, 0.99)
    pr[pr > q] <- q
    pr <- pr/max(pr)
    pr[excls] <- NA
    qq <- quantile(pr, seq(0.1, 0.9, 0.1), na.rm=TRUE)
    z <- cut(pr, c(-1, unique(qq), 2))
    Col <- rev(terrain.colors(nlevels(z)))
	png(file=fname, width=600, height=1000)

    plot(kgrid$X, kgrid$Y, pch=15, cex=0.2, col=Col[z], axes=FALSE, ann=FALSE)
    points(kgrid$X[excls], kgrid$Y[excls], pch=15, cex=0.2, col="darkgrey")
    mtext(paste0(NAM, ", South"), line=2, side=3, adj=0.5, cex=1.4, col="grey40")
    points(xyw, pch=15, cex=0.2, col=rgb(0.3,0.45,0.9))
    points(city, pch=18, col="grey10")
    text(city, rownames(city), cex=0.8, adj=-0.1, col="grey10")
    legend("bottomleft", col=rev(Col), fill=rev(Col),
        legend=c("High", rep("", length(Col)-2), "Low"), bty="n")

	dev.off()
}
}

clim_N <- data.frame(tax[names(clim_n), c("English_Name","Scientific_Name")], 
    do.call(rbind, clim_n))
clim_S <- data.frame(tax[names(clim_s), c("English_Name","Scientific_Name")], 
    do.call(rbind, clim_s))
clim_N <- clim_N[rownames(slt)[slt$modelN],]
clim_S <- clim_S[rownames(slt)[slt$modelS],]

write.csv(clim_N, file=file.path(ROOT, "figs", "climatehf-north.csv"))
write.csv(clim_S, file=file.path(ROOT, "figs", "climatehf-south.csv"))


## surroundinghf-north
## surroundinghf-south

for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(tax[spp, "English_Name"])
if (tax[spp, "surroundinghf_north"]) {
    resn <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_sp <- getEst(resn, stage=stage_hab_n + 2, na.out=FALSE, Xnn)
    fname <- file.path(ROOT, "figs", "surroundinghf-north", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    png(file=fname, width=7.5, height=5.7, units="in", res=300)
    op <- par(mai=c(0.9,1,0.2,0.3))
    fig_hf_noremn(estn_sp, Xnn, LAB=paste0(NAM, ", North"))
    par(op)
    dev.off()
}
if (tax[spp, "surroundinghf_south"]) {
    ress <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-south_", spp, ".Rdata")))
    ests_sp <- getEst(ress, stage=stage_hab_s + 2, na.out=FALSE, Xns)
    fname <- file.path(ROOT, "figs", "surroundinghf-south", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    png(file=fname, width=7.5, height=5.7, units="in", res=300)
    op <- par(mai=c(0.9,1,0.2,0.3))
    fig_hf_noremn(ests_sp, Xns, LAB=paste0(NAM, ", North"))
    par(op)
    dev.off()
}
}

## trend

res_trend <- matrix(NA, nrow(tax), 10)
colnames(res_trend) <- c("Mean_North","Median_North","LCL_North","UCL_North","n_North",
    "Mean_South","Median_South","LCL_South","UCL_South","n_South")
res_trend[,5] <- tax$ndet_n
res_trend[,10] <- tax$ndet_s
rownames(res_trend) <- rownames(tax)
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
if (tax[spp, "trend_north"]) {
    resn <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_yr <- getEst(resn, stage=stage_hab_n + 3, na.out=FALSE, Xnn)
    yr_n <- 100 * (exp(estn_yr[,"YR"]) - 1)
    res_trend[spp, 1:4] <- fstat(yr_n)
    NDATN <- sum(yyn[,spp] > 0)
    NN <- aggregate(yyn[,spp], list(year=xnn$YEAR), mean)
}
if (tax[spp, "trend_south"]) {
    ress <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-south_", spp, ".Rdata")))
    ests_yr <- getEst(ress, stage=stage_hab_s + 3, na.out=FALSE, Xns)
    yr_s <- 100 * (exp(ests_yr[,"YR"]) - 1)
    res_trend[spp, 6:9] <- fstat(yr_s)
    NDATS <- sum(yys[,spp] > 0)
    NS <- aggregate(yys[,spp], list(year=xns$YEAR), mean)
}
if (tax[spp, "trend_north"] | tax[spp, "trend_south"]) {
    NAM <- as.character(tax[spp, "English_Name"])
    fname <- file.path(ROOT, "figs", "trend", 
        paste0(as.character(tax[spp, "file"]), ".png"))
    png(file=fname, width=600, height=600)
    op <- par(mfrow=c(2,2), cex=0.8)
    if (tax[spp, "trend_north"]) {
        plot(NN, ylab="Annual Mean Abundance Index", xlab="Year", 
            type="b", col=1, pch=19,
            main=paste0(NAM, ", North (n = ", NDATN, " detections)"))
        abline(lm(x ~ year, NN), col="red4", lty=1, lwd=2)
        hist(yr_n, col="gold", xlab="Decadal Trend (%)", main="")
        abline(v=fstat(yr_n)[1], col="red4", lty=1, lwd=2)
        abline(v=fstat(yr_n)[3:4], col="red4", lty=2, lwd=1)
    } else {
        plot.new()
        plot.new()
    }
    if (tax[spp, "trend_south"]) {
        plot(NS, ylab="Annual Mean Abundance Index", xlab="Year", 
            type="b", col=1, pch=19,
            main=paste0(NAM, ", South (n = ", NDATS, " detections)"))
        abline(lm(x ~ year, NS), col="red4", lty=1, lwd=2)
        hist(yr_n, col="gold", xlab="Decadal Trend (%)", main="")
        abline(v=fstat(yr_s)[1], col="red4", lty=1, lwd=2)
        abline(v=fstat(yr_s)[3:4], col="red4", lty=2, lwd=1)
    } else {
        plot.new()
        plot.new()
    }
    par(op)
    dev.off()
}
}
res_trend2 <- data.frame(tax[,c("English_Name","Scientific_Name")], res_trend)
write.csv(res_trend2, file=file.path(ROOT, "figs", "trend.csv"))

rank_fun(res_trend$Mean_North, res_trend$LCL_North, res_trend$UCL_North,
    n=res_trend$n_North, col=1, lab = rownames(res_trend))

rank_fun(res_trend$Mean_South, res_trend$LCL_South, res_trend$UCL_South,
    n=res_trend$n_South, col=1, lab = rownames(res_trend))

## ARU effect
res_aru <- list()
for (spp in rownames(tax[tax$surroundinghf_north,])) {
    pres <- sum(yyn[substr(rownames(yyn), 1, 5) == "EMCLA",spp] > 0)
    if (pres > 0) {
        cat(spp, "\n");flush.console()
        resn <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-north_", spp, ".Rdata")))
        estn <- getEst(resn, stage=6, na.out=FALSE, Xnn)
        aru <- estn[,"ARU"]
        res_aru[[spp]] <- c(aru, pres)
    }
}
res_aru <- do.call(rbind, res_aru)
tmp <- res_aru[res_aru[,241]>19,-241]
tmp[tmp == 0] <- NA
rowMeans(exp(tmp), na.rm=TRUE)


## Linear features coefficients

spp <- "BTNW"
xlin <- nonDuplicated(xnn[,c("ROAD01","hab1","hab_lcc","hab_lcc2","hab_lcc3")],
    hab1, TRUE)
xlin <- xlin[c("Decid", "Mixwood", "Conif", "Pine", "BSpr", "Larch", 
    "Decid", "Mixwood", "Conif", "Pine", "BSpr", "Larch", 
    "GrassHerb", "Shrub", "Wetland", "Cult", "UrbInd"),]
rownames(xlin)[1:12] <- paste0(rep(rownames(xlin)[1:6], 2), 
    rep(c("0-40","40+"), each=6))
xlin$ROAD01 <- 1
xlin$SoftLin_PC <- 0
xlin$hab_lcc[] <- c(4,4, 3,3,3,3, 2,2, 1,1,1,1, 5,5,5,5,5)
xlin$hab_lcc3 <- xlin$hab_lcc
levels(xlin$hab_lcc3) <- c("1", "1", "2", "2", "3")
xlin$hab_lcc2 <- xlin$hab_lcc
levels(xlin$hab_lcc2) <- c("1", "1", "1", "1", "2")


Xlin <- model.matrix(getTerms(modsn["Contrast"], "formula", intercept=TRUE), xlin)
colnames(Xlin) <- fixNames(colnames(Xlin))
Xlin <- Xlin[,-1]


res_soft <- list()
res_hard <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(tax[spp, "English_Name"])
if (tax[spp, "veghf_north"]) {
    resn <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_lin <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    colnames(estn_lin) <- fixNames(colnames(estn_lin))
    estn_lin2 <- estn_lin[,colnames(Xlin)]
    pr <- apply(estn_lin2, 1, function(z) Xlin %*% z)
    rownames(pr) <- rownames(xlin)
    tab <- t(apply(exp(pr), 1, quantile, c(0.5, 0.05, 0.95)))
    res_hard[[spp]] <- data.frame(Species=spp, Habitat=rownames(tab), tab)
    res_soft[[spp]] <- quantile(estn_lin[,"SoftLin_PC"], c(0.5, 0.05, 0.95))
}
}
## note: roadside stuff is exponentiated, but soft lin is not,
## because it is exp(x * est)

softlin <- data.frame(Species=tax[names(res_soft), "English_Name"], do.call(rbind, res_soft))
hardlin <- do.call(rbind, res_hard)
hardlin$Species <- tax[as.character(hardlin$Species), "English_Name"]

softlin <- droplevels(softlin[rownames(slt)[slt$veghf.north],])
hardlin <- droplevels(hardlin[hardlin$Species %in% softlin$Species,])

write.csv(softlin, row.names=FALSE, 
    file=file.path(ROOT, "figs", "soft-linear-coefs-2015.csv"))
write.csv(hardlin, row.names=FALSE, 
    file=file.path(ROOT, "figs", "hard-linear-EXPcoefs-2015.csv"))
    
softlin2 <- softlin[c("BTNW","BBWA","OVEN","BRCR","CAWA"),]
hardlin2 <- do.call(rbind, res_hard[c("BTNW","BBWA","OVEN","BRCR","CAWA")])
hardlin2$Species <- tax[as.character(hardlin2$Species), "English_Name"]

write.csv(softlin2, row.names=FALSE, 
    file=file.path(ROOT, "figs", "soft-linear-coefs-2015-5spp.csv"))
write.csv(hardlin2, row.names=FALSE, 
    file=file.path(ROOT, "figs", "hard-linear-EXPcoefs-2015-5spp.csv"))

## upland/lowland classification of species

tax2 <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(tax2) <- tax2$AOU
tax3 <- read.csv("~/repos/abmianalytics/lookup/vertebrate-guilds.csv")
rownames(tax3) <- tax3$AOU.Code
setdiff(tax2$AOU[tax2$map.pred], tax3$AOU.Code)
setdiff(tax2$AOU[tax2$map.pred], tax$Species_ID)
SPP <- intersect(tax2$AOU[tax2$map.pred], tax3$AOU.Code)
tax2 <- droplevels(tax2[SPP,])
tax3 <- droplevels(tax3[SPP,])
native <- tax3[,grep("Native.to.", colnames(tax3))]
native[is.na(native)] <- 0
native[native > 0] <- 1
wet <- tax3[,c("General.Habitat.Category.Bog", "General.Habitat.Category.WetAq",
    "Wetland.Types.Wet_NestTerrOrWet", "Wetland.Types.Aq_NestTerrOrWet")]
wet[is.na(wet)] <- 0

tax2$native <- ifelse(rowSums(native)>0, 1, 0)
tax2 <- cbind(tax2, wet)

dat2 <- dat[dat$useOK & dat$keep,]
wetcl <- c("BSpr","Larch","Wetland")
dat2$strat <- as.factor(ifelse(dat2$hab1 %in% wetcl, "lowland", "upland"))
yy2 <- as.matrix(yy[rownames(dat2), SPP])
off2 <- e$OFFmean[rownames(dat2)]

table(dat2$strat, dat2$pWater >0.5)
dat2$strat[dat2$pWater >0.5] <- "lowland"

library(opticut)

XXX <- model.matrix(~ ROAD01 + SoftLin_PC, dat2)
oc1 <- opticut1(yy2[,1], XXX, dat2$strat, dist="poisson")

oc <- opticut(yy2 ~ ROAD01 + SoftLin_PC, dat2, strata=dat2$strat, 
    offset=off2, dist="poisson", comb="rank")
os <- summary(oc)$summary
os <- os[SPP,]

tax2v <- data.frame(tax2[SPP,], os[SPP,])
tax2v$w <- NULL
tax2v$ndet_n <- NULL
tax2v$ndet_s <- NULL
tax2v$ndet_ns <- NULL
tax2v$map.det <- NULL
tax2v$veghf.north <- NULL
tax2v$soilhf.south <- NULL
tax2v$map.pred <- NULL
tax2v$useavail.north <- NULL
tax2v$useavail.south <- NULL
tax2v$lablo <- NULL
tax2v$labhi <- NULL

#levels(tax2v$split) <- c("lowland", "upland", "nopref")
#tax2v$split[tax2v$logLR < 2] <- "nopref"
table(tax2v$split)

tax2v$order <- tax3[SPP, "Order"]
tax2v$split2 <- as.character(tax2v$split)
tax2v$split2[] <- ""
tax2v$split2[tax2v$General.Habitat.Category.Bog + 
    tax2v$General.Habitat.Category.WetAq > 0 & tax2v$split == "lowland"] <- "lowland"
tax2v$split2[tax2v$General.Habitat.Category.Bog + 
    tax2v$General.Habitat.Category.WetAq == 0 & tax2v$split == "upland"] <- "upland"
tax2v$split2[tax2v$order %in% c("ANSERIFORMES","CHARADRIIFORMES","CICONIIFORMES",
    "PODICIPEDIFORMES","PELECANIFORMES","GAVIIFORMES","GRUIFORMES")] <- "lowland"
tax2v$split2[tax2v$order %in% c("COLUMBIFORMES","FALCONIFORMES",
    "GALLIFORMES","PICIFORMES","STRIGIFORMES")] <- "upland"
tax2v$split2[tax2v$native == 0] <- "nonnative"

table(tax2v$order,tax2v$split)
table(tax2v$split2)
write.csv(tax2v, file="~/birds-upland-lowland-classification.csv", row.names=FALSE)

