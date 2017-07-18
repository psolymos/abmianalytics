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
TAX <- droplevels(e$TAX)
TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")
yy <- e$YYw
yy01 <- yy
yy01[yy01>0] <- 1

yrf <- read.csv("e:\\peter\\AB_data_v2016\\oracle\\birds-rf.csv")
ysm <- read.csv("e:\\peter\\AB_data_v2016\\oracle\\birds-aru.csv")
yrfsm <- data.frame(name=c(as.character(yrf$COMMON_NAME), as.character(ysm$COMMON_NAME)),
    tsn=c(as.character(yrf$TSNID), as.character(ysm$TSNID)))
yrfsm <- nonDuplicated(yrfsm, name)
levels(yrfsm$name)[levels(yrfsm$name) == "Black-crowned Night Heron"] <- "Black-crowned Night-Heron"
levels(yrfsm$name)[levels(yrfsm$name) == "Black and White Warbler"] <- "Black-and-white Warbler"
levels(yrfsm$name)[levels(yrfsm$name) == "Macgillivray Warbler"] <- "MacGillivray's Warbler"
compare_sets(TAX$English_Name, yrfsm$name)
setdiff(TAX$English_Name, yrfsm$name)
TAX$TSNID <- yrfsm$tsn[match(TAX$English_Name, yrfsm$name)]

en <- new.env()
load(file.path(ROOT, "data", "data-north.Rdata"), envir=en)
xnn <- en$DAT
modsn <- en$mods
yyn <- en$YY

es <- new.env()
load(file.path(ROOT, "data", "data-south.Rdata"), envir=es)
xns <- es$DAT
modss <- es$mods
yys <- es$YY

rm(e, en, es)

compare_sets(rownames(TAX), colnames(yyn))
compare_sets(rownames(TAX), colnames(yys))

#yyn <- yy[rownames(yyn0),]
#yys <- yy[rownames(yys0),]

## model for species
fln <- list.files(file.path(ROOT, "results", "north"))
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- list.files(file.path(ROOT, "results", "south"))
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)

SPP <- sort(unique(c(fln, fls, colnames(yy))))
tax <- droplevels(TAX[SPP, c("Spp","English_Name","Scientific_Name","TSNID","Family_Sci")])
tax$AOUcode <- rownames(tax)
tax$modelN <- rownames(tax) %in% fln
tax$modelS <- rownames(tax) %in% fls
tax$ndet_n <- colSums(yyn>0)[match(rownames(tax), colnames(yyn))]
tax$ndet_s <- colSums(yys>0)[match(rownames(tax), colnames(yys))]
tax$ndet_n[is.na(tax$ndet_n)] <- 0
tax$ndet_s[is.na(tax$ndet_s)] <- 0

## tax placeholders for all the output
tax$ndet_ns <- colSums(yy01[,rownames(tax)])
tax$ndet <- tax$ndet_ns
tax$map_det <- tax$ndet_ns > 0
tax$useavail_north <- tax$ndet_n > 3
tax$useavail_south <- tax$ndet_s > 3
tax$trend_north <- tax$modelN
tax$trend_south <- tax$modelS
tax$veghf_north <- tax$modelN
tax$soilhf_nontreed_south <- tax$modelS
tax$soilhf_treed_south <- tax$modelS
tax$linear_north <- tax$modelN
tax$linear_south <- tax$modelS
tax$surroundinghf_north <- tax$modelN
tax$surroundinghf_south <- tax$modelS

## this is hand edited version
spp_OK <- read.csv("e:/peter/AB_data_v2016/out/birds/tables/birds-lookup-final.csv")
rownames(spp_OK) <- spp_OK$AOU

## species lookup table for web
slt <- data.frame(sppid=tax$Spp,
    species=tax$English_Name,
    scinam=tax$Scientific_Name,
    tsnid=tax$TSNID,
    AOU=tax$AOUcode,
    family=tax$Family_Sci,
    tax[,c("ndet","modelN","modelS","ndet_n","ndet_s", "ndet_ns")])
slt$modelN <- spp_OK[rownames(slt),"showN"]
slt$modelS <- spp_OK[rownames(slt),"showS"]
slt$known_issues <- spp_OK[rownames(slt),"Known_Problems"]

slt$map.det <- tax$map_det
#slt$veghf.north <- tax$modelN & tax$ndet_n > 99
#slt$soilhf.south <- tax$modelS & tax$ndet_s > 49
slt$veghf.north <- slt$modelN
slt$soilhf.south <- slt$modelS
slt$map.pred <- slt$veghf.north | slt$soilhf.south
slt$useavail.north=tax$useavail_north & !slt$veghf.north
slt$useavail.south=tax$useavail_south & !slt$soilhf.south
slt <- slt[slt$map.det,]

gl <- read.csv("~/repos/abmianalytics/lookup/vertebrate-guilds.csv")
intersect(gl$AOU.Code, rownames(slt))
slt[setdiff(rownames(slt), gl$AOU.Code),1:3]

slt$oldforest <- gl$Forest.Types.Old[match(rownames(slt), gl$AOU.Code)]
slt$oldforest[is.na(slt$oldforest)] <- 0
slt$oldforest[slt$oldforest > 0] <- 1

if (FALSE) {
    sb <- read.csv("~/repos/bamanalytics/lookup/singing-species.csv")
    rownames(sb) <- sb$Species_ID
    compare_sets(rownames(tax), rownames(sb))
    sb2 <- tax[setdiff(rownames(tax), rownames(sb)),c("English_Name",
        "Family_Sci")]
    sb2 <- sb2[rownames(sb2) != "NONE",]
    sb2$Singing_birds <- NA
    sb3 <- rbind(sb[intersect(rownames(tax), rownames(sb)),colnames(sb2)], sb2)
    write.csv(sb3, file=file.path(ROOT, "tables", "birds-lookup-temp.csv"))
    #write.csv(sb3, row.names=F, file="Bird-spp.csv")
}

sb <- read.csv("~/repos/abmianalytics/lookup/singing-species-alberta.csv")
slt$singing <- sb$Singing_birds[match(rownames(slt), sb$Species_ID)]
#write.csv(slt, row.names=FALSE, file="~/repos/abmispecies/_data/birds.csv")
#write.csv(slt, row.names=FALSE,
#    file=file.path(ROOT, "tables", "birds-lookup.csv"))

slt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(slt) <- slt$AOU
slt$comments <- NULL

## terms and design matrices
nTerms <- getTerms(modsn, "list")
sTerms <- getTerms(modss, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), xnn)
colnames(Xnn) <- fixNames(colnames(Xnn))
Xns <- model.matrix(getTerms(modss, "formula"), xns)
colnames(Xns) <- fixNames(colnames(Xns))

stage_hab_n <- 5
stage_hab_s <- 3

## computing times

t_n <- numeric(length(fln))
names(t_n) <- fln
for (spp in fln) {
    cat(spp, "\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "results", "north", paste0("birds_abmi-north_", spp, ".Rdata")))
    t_n[spp] <- attr(resn, "ncl") * abs(attr(resn, "timing")["elapsed"])
}
t_s <- numeric(length(fls))
names(t_s) <- fls
for (spp in fls) {
    cat(spp, "\n");flush.console()
    ress <- loadSPP(file.path(ROOT, "results", "south", paste0("birds_abmi-south_", spp, ".Rdata")))
    t_s[spp] <- attr(ress, "ncl") * abs(attr(ress, "timing")["elapsed"])
}
summary(t_n)
summary(t_s)
((sum(c(t_n, t_s))/60)/60)/24 # 475.5 days

## EC2: 1.68 US / h, 19170.64
1.68 * (sum(c(t_n, t_s))/60)/60

## spp specific output

spp <- "BTNW"


## veghf-north
## linear-north
## table: veghf-north

for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(tax[spp, "English_Name"])
if (tax[spp, "veghf_north"]) {
    resn <- loadSPP(file.path(ROOT, "results", "north", paste0("birds_abmi-north_", spp, ".Rdata")))
    estn7 <- getEst(resn, stage=7, na.out=FALSE, Xnn)
    fname <- file.path(ROOT, "tables", "coefs",
        paste0(as.character(tax[spp, "Spp"]), "_Stage7_coefs.csv"))
    write.csv(estn7, row.names=FALSE, file=fname)
}
}

## all N & S data for plots
res_coef <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    res_coef[[spp]] <- list(veg=NULL, soil=NULL, max=c(NA, NA))
    if (tax[spp, "veghf_north"]) {
        resn <- loadSPP(file.path(ROOT, "results", "north",
            paste0("birds_abmi-north_", spp, ".Rdata")))
        estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
        tmp1 <- pred_veghf(estn_hab, Xnn, burn_included=FALSE)
        res_coef[[spp]]$veg <- tmp1
        res_coef[[spp]]$max[1] <- fig_veghf_ymax(tmp1)
    }
    if (tax[spp, "soilhf_treed_south"]) {
        ress <- loadSPP(file.path(ROOT, "results", "south",
            paste0("birds_abmi-south_", spp, ".Rdata")))
        ests_hab <- getEst(ress, stage=stage_hab_s, na.out=FALSE, Xns)
        tmp2 <- pred_soilhf(ests_hab, Xns)
        res_coef[[spp]]$soil <- tmp2
        res_coef[[spp]]$max[2] <- max(fig_soilhf_ymax(tmp2$treed), fig_soilhf_ymax(tmp2$nontreed))
    }
}
save(res_coef, file=file.path(ROOT, "tables", "res_coef.Rdata"))

## get all kinds of linear coefs

get_lin <- function(tp, what="veg") {
    if (what=="veg")
        pr <- unname(attr(tp$veg, "linear"))
    if (what=="soil")
        pr <- unname(tp$soil$linear)
    c(AverageCoef=pr[1],
        SoftLin10=pr[1] * exp(0.1*log(pr[2])),
        HardLin10=pr[1] * pr[5],
        SoftLin	=pr[2],
        SoftLin.LCI	=pr[3],
        SoftLin.UCI	=pr[4],
        HardLin	=pr[5],
        HardLin.LCI	=pr[6],
        HardLin.UCI=pr[7])
}

lin1_n <- lin1_s <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    tp <- res_coef[[spp]]
    if (!all(is.na(tp$max))) {
        if (tax[spp, "veghf_north"])
            lin1_n[[spp]] <- get_lin(tp, "veg")
        if (tax[spp, "soilhf_treed_south"])
            lin1_s[[spp]] <- get_lin(tp, "soil")
    }
}

lin1_n <- do.call(rbind, lin1_n)
lin1_s <- do.call(rbind, lin1_s)
lin2_n <- data.frame(tax[rownames(lin1_n), c("English_Name","Scientific_Name","TSNID")],
    lin1_n)
lin2_s <- data.frame(tax[rownames(lin1_s), c("English_Name","Scientific_Name","TSNID")],
    lin1_s)
write.csv(lin2_n,
    file=file.path(ROOT, "tables", "birds-linearfull-north.csv"), row.names=FALSE)
write.csv(lin2_s,
    file=file.path(ROOT, "tables", "birds-linearfull-south.csv"), row.names=FALSE)


## plots for N & S
lin_n <- lin_s <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    tp <- res_coef[[spp]]
    if (!all(is.na(tp$max))) {
        NAM <- as.character(tax[spp, "English_Name"])
        if (max(tp$max, na.rm=TRUE) > 3*min(tp$max, na.rm=TRUE)) {
            MAXn <- tp$max[1]
            MAXs <- tp$max[2]
        } else {
            MAXn <- max(tp$max, na.rm=TRUE)
            MAXs <- max(tp$max, na.rm=TRUE)
        }
        if (is.na(MAXn))
            MAXn <- max(tp$max, na.rm=TRUE)
        if (is.na(MAXs))
            MAXs <- max(tp$max, na.rm=TRUE)

        if (tax[spp, "veghf_north"]) {
            prn <- tp$veg
            NDAT <- sum(yyn[,spp] > 0)
            ## veghf
            fname <- file.path(ROOT, "figs", "veghf-north",
                paste0(as.character(tax[spp, "Spp"]), ".png"))
            png(file=fname,width=1500,height=700)
            fig_veghf(prn,
                paste0(NAM, " (n = ", NDAT, " detections)"), ymax=MAXn)
            dev.off()
            ## linear
            fname <- file.path(ROOT, "figs", "linear-north",
                paste0(as.character(tax[spp, "Spp"]), ".png"))
            png(file=fname,width=350,height=400)
            lin_n[[spp]] <- fig_linear(attr(prn, "linear"),
                paste0(NAM, "\nNorth (n = ", NDAT, " det.)"))
            dev.off()
        }
        if (tax[spp, "soilhf_treed_south"]) {
            prs <- tp$soil
            NDAT <- sum(yys[,spp] > 0)
            ## treed
            fname <- file.path(ROOT, "figs", "soilhf-treed-south",
                paste0(as.character(tax[spp, "Spp"]), ".png"))
            png(file=fname,width=500,height=450)
            fig_soilhf(prs$treed,
                paste0(NAM, ", South, Treed (n = ", NDAT, " detections)"),
                ymax=MAXs)
            dev.off()
            ## nontreed
            fname <- file.path(ROOT, "figs", "soilhf-nontreed-south",
                paste0(as.character(tax[spp, "Spp"]), ".png"))
            png(file=fname,width=500,height=450)
            fig_soilhf(prs$nontreed,
                paste0(NAM, ", South, Non-treed (n = ", NDAT, " detections)"),
                ymax=MAXs)
            dev.off()
            ## linear
            fname <- file.path(ROOT, "figs", "linear-south",
                paste0(as.character(tax[spp, "Spp"]), ".png"))
            png(file=fname,width=350,height=400)
            lin_s[[spp]] <- fig_linear(prs$linear,
                paste0(NAM, "\nSouth (n = ", NDAT, " det.)"))
            dev.off()
        }
    }
}

lin_n <- do.call(rbind, lin_n)
lin_s <- do.call(rbind, lin_s)

f1 <- function(spp) {
    x <- res_coef[[spp]]$veg
    pr <- attr(x, "linear")
    rr <- c(pr[1], # Average
        pr[1] * exp(0.1*log(pr[2:4])), # SoftLin
        pr[1] * pr[5:7]) # HardLin
    names(rr) <- c("AverageCoef", "SoftLinear", "SoftLinear.LCL", "SoftLinear.UCL",
        "HardLinear", "HardLinear.LCL", "HardLinear.UCL")
    #x <- x[rownames(x) != "Burn",c(2,3,4)]
    x <- x[rownames(x) != "Burn",c(2,3,4)]
    rownames(x) <- gsub(" ", "", rownames(x))
    xx <- t(x)
    dim(xx) <- NULL
    names(xx) <- paste0(rep(rownames(x), each=3), c("", ".LCL", ".UCL"))
    c(xx, rr)
}
vhf <- t(sapply(fln, f1))
vhf2 <- data.frame(tax[rownames(vhf), c("English_Name","Scientific_Name","TSNID")],
    vhf)
vhf2 <- vhf2[rownames(slt)[slt$veghf.north],]
write.csv(vhf2,
    file=file.path(ROOT, "tables", "birds-veghf-north.csv"), row.names=FALSE)

if (FALSE) {
SPP <- rownames(slt)[slt$veghf.north]

vhf2 <- read.csv(file.path(ROOT, "figs", "birds-veghf-north.csv"))
rownames(vhf2) <- vhf2$X
vhf2 <- vhf2[,-(1:3)]

excl <- c(grep(".LCL", colnames(vhf2)), grep(".UCL", colnames(vhf2)))
vhf2 <- as.matrix(vhf2[,-excl])
Max <- apply(vhf2[,1:(ncol(vhf2)-2)], 1, max)
vhf2 <- vhf2 / Max
vhf2[vhf2[,"HardLinear"] > 5,"HardLinear"] <- 5
}

f2 <- function(spp) {
    x <- res_coef[[spp]]$soil
    xTreed <- mean((x$treed / x$nontreed)[,2])
    pr <- x$linear
    x <- x$nontreed
    rownames(x) <- gsub(" ", "", rownames(x))

    rr <- c(pr[1], # Average
        pr[1] * exp(0.1*log(pr[2:4])), # SoftLin
        pr[1] * pr[5:7]) # HardLin
    names(rr) <- c("AverageCoef", "SoftLinear", "SoftLinear.LCL", "SoftLinear.UCL",
        "HardLinear", "HardLinear.LCL", "HardLinear.UCL")

    xx <- t(x[,2:4])
    dim(xx) <- NULL
    names(xx) <- paste0(rep(rownames(x), each=3), c("", ".LCL", ".UCL"))
    c(xx, rr, xTreed=xTreed)
}
soil <- t(sapply(fls, f2))
soil2 <- data.frame(tax[rownames(soil), c("English_Name","Scientific_Name","TSNID")],
    soil)
soil2 <- soil2[rownames(slt)[slt$soilhf.south],]
write.csv(soil2,
    file=file.path(ROOT, "tables", "birds-soilhf-south.csv"), row.names=FALSE)

## sector effects tables
load(file.path(ROOT, "tables", "sector-effects.Rdata"))

nres <- list()
sres <- list()
for (spp in names(seff_res)) {
    nres[[spp]] <- 100*c(PopEffect=seff_res[[spp]]$N[,2], UnitEffect=seff_res[[spp]]$N[,3])
    sres[[spp]] <- 100*c(PopEffect=seff_res[[spp]]$S[,2], UnitEffect=seff_res[[spp]]$S[,3])
}
nres <- do.call(rbind, nres)
sres <- do.call(rbind, sres)
nres <- data.frame(tax[names(seff_res), c("English_Name","Scientific_Name","TSNID")],
    round(nres, 6))[fln,]
sres <- data.frame(tax[names(seff_res), c("English_Name","Scientific_Name","TSNID")],
    round(sres, 6))[fls,]

round(sapply(seff_res[[1]], function(z) 100*z[,1]), 2)

nres <- nres[rownames(slt)[slt$veghf.north],]
sres <- sres[rownames(slt)[slt$soilhf.south],]

write.csv(nres, row.names=FALSE,
    file=file.path(ROOT, "tables", "Birds_SectorEffects_North.csv"))
write.csv(sres, row.names=FALSE,
    file=file.path(ROOT, "tables", "Birds_SectorEffects_South.csv"))



## climate & surrounding hf tables, climate surface maps

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
    resn <- loadSPP(file.path(ROOT, "results", "north",
        paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_sp <- getEst(resn, stage=stage_hab_n + 3, na.out=FALSE, Xnn)
    fname <- file.path(ROOT, "figs", "surroundinghf-north",
        paste0(as.character(tax[spp, "Spp"]), ".png"))
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

res_lin <- list()
for (spp in rownames(tax)) {
    cat(spp, "\n");flush.console()
    NAM <- as.character(tax[spp, "English_Name"])
if (tax[spp, "veghf_north"]) {
#    resn <- loadSPP(file.path(ROOT, "results", paste0("birds_abmi-north_", spp, ".Rdata")))
    resn <- loadSPP(file.path(ROOT, "results", "north",
            paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_lin <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    colnames(estn_lin) <- fixNames(colnames(estn_lin))
    ehl0 <- exp(estn_lin[,c("ROAD01")])
    ehl1 <- exp(estn_lin[,c("ROAD01")] + estn_lin[,c("habCl:ROAD01")])
    sl <- estn_lin[,c("SoftLin_PC")]
    res_lin[[spp]] <- c(expHardLin_Open=fstat(ehl0),
        expHardLin_Closed=fstat(ehl1),
        SoftLin=fstat(sl))
}
}
## note: roadside stuff is exponentiated, but soft lin is not,
## because it is exp(x * est)

## habCl is 0/1:
## 0 = 40yrs or older forest % in 150m radius buffer <= 50%
## 1 = 40yrs or older forest % in 150m radius buffer > 50%

lin <- data.frame(Species=tax[names(res_lin), "English_Name"], do.call(rbind, res_lin))

write.csv(lin, row.names=FALSE,
    file=file.path(ROOT, "tables", "linear-coefs-2016.csv"))

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

## desperate attempt to clean up labeling mess
library(mefa4)
v <- read.csv("e:/peter/sppweb2016/round01/tables/col-xwalk.csv")

## mosses
tax1 <- read.csv("e:/peter/sppweb2016/round01/tables/mosses-tax-eta.csv")
tax2 <- read.csv("e:/peter/sppweb2016/round01/tables/mosses-tax-tsn.csv")
compare_sets(tax1$Species, tax2$SCIENTIFIC_NAME)
tax1$TSNID <- tax2$TSN_ID[match(tax1$Species, tax2$SCIENTIFIC_NAME)]
write.csv(tax1, file="e:/peter/sppweb2016/round01/tables/mosses-tax-eta-tsn.csv")

y <- read.csv("e:/peter/sppweb2016/round01/tables/mosses-veghf.csv")
z <- data.frame(Species=y$Species)
for (i in 1:nrow(v)) {
    if (as.character(v$eta)[i] != "") {
        z[[as.character(v$ref[i])]] <- y[[as.character(v$eta)[i]]]
    } else {
        z[[as.character(v$ref[i])]] <- NA
    }
}
write.csv(z, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/mosses-veghf2.csv")

yy <- as.matrix(y[,-1])
yy <- yy[,!grepl(".LCI", colnames(yy))]
yy <- yy[,!grepl(".UCI", colnames(yy))]
yy <- yy[,!(colnames(yy) %in% c("SoftLin", "HardLin"))]
yyy <- data.frame(Species=y$Species, AverageCoef=rowMeans(yy, na.rm=TRUE), y[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
write.csv(yyy, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/mosses-lin2.csv")

## NN plants
y <- read.csv("e:/peter/sppweb2016/round01/tables/nnplants-veghf.csv")
y <- y[1,,drop=FALSE]
z <- data.frame(Species=y$Species)
for (i in 1:nrow(v)) {
    if (as.character(v$eta)[i] != "") {
        z[[as.character(v$ref[i])]] <- y[[as.character(v$eta)[i]]]
    } else {
        z[[as.character(v$ref[i])]] <- NA
    }
}
write.csv(z, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/nnplants-veghf2.csv")

## mites
y <- read.csv("e:/peter/sppweb2016/round01/tables/mites-veghf.csv")
z <- data.frame(Species=y$Species)
for (i in 1:nrow(v)) {
    if (as.character(v$eta)[i] != "") {
        z[[as.character(v$ref[i])]] <- y[[as.character(v$eta)[i]]]
    } else {
        z[[as.character(v$ref[i])]] <- NA
    }
}
write.csv(z, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/mites-veghf2.csv")

## mammals
y <- read.csv("e:/peter/sppweb2016/round01/tables/mammals-veghf.csv")
z <- data.frame(SpeciesID=y$SpeciesID)
for (i in 1:nrow(v)) {
    if (as.character(v$eta)[i] != "") {
        z[[as.character(v$ref[i])]] <- y[[as.character(v$eta)[i]]]
    } else {
        z[[as.character(v$ref[i])]] <- NA
    }
}
write.csv(z, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/mammals-veghf2.csv")

## lichens
tax1 <- read.csv("e:/peter/sppweb2016/round01/tables/lichens-tax-eta.csv")
#tax2 <- read.csv("e:/peter/AB_data_v2016/data/species/OUT_Lichens_Species_Taxa_2016-04-28.csv")
tax3 <- read.csv("e:/peter/AB_data_v2016/oracle/lichens-0308.csv")
tax4 <- read.csv("e:/peter/AB_data_v2016/oracle/lichens-09.csv")

tax2 <- rbind(nonDuplicated(tax3[,c("SCIENTIFIC_NAME", "TSNID")], tax3$SCIENTIFIC_NAME),
    nonDuplicated(tax4[,c("SCIENTIFIC_NAME", "TSNID")], tax4$SCIENTIFIC_NAME))
tax2 <- nonDuplicated(tax2, SCIENTIFIC_NAME)
tax2 <- tax2[order(tax2[,1]),]

compare_sets(tax1$scinam, tax2$SCIENTIFIC_NAME)
tax1$TSNID <- tax2$TSNID[match(tax1$scinam, tax2$SCIENTIFIC_NAME)]
write.csv(tax1, row.names=FALSE, na="", file="e:/peter/sppweb2016/round01/tables/lichens-tax-eta-tsn.csv")

write.csv(tax2, row.names=FALSE, na="", file="e:/peter/sppweb2016/round01/tables/lichens-tax-tsn.csv")

x <- read.csv("e:/peter/sppweb2016/round01/tables/lichens-use-s.csv")
x$Species <- tax1$scinam[match(x$SpLabel, tax1$sppid)]
write.csv(x, row.names=FALSE, na="", file="e:/peter/sppweb2016/round01/tables/lichens-use-s2.csv")

y <- read.csv("e:/peter/sppweb2016/round01/tables/lichens-veghf.csv")
z <- data.frame(Species=y$Species)
for (i in 1:nrow(v)) {
    if (as.character(v$eta)[i] != "") {
        z[[as.character(v$ref[i])]] <- y[[as.character(v$eta)[i]]]
    } else {
        z[[as.character(v$ref[i])]] <- NA
    }
}
write.csv(z, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/lichens-veghf2.csv")

yy <- as.matrix(y[,-1])
yy <- yy[,!grepl(".LCI", colnames(yy))]
yy <- yy[,!grepl(".UCI", colnames(yy))]
yy <- yy[,!(colnames(yy) %in% c("SoftLin", "HardLin"))]
yyy <- data.frame(Species=y$Species, AverageCoef=rowMeans(yy, na.rm=TRUE), y[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
write.csv(yyy, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/lichens-lin2.csv")

x <- read.csv("e:/peter/sppweb2016/round01/tables/lichens-soil.csv")
y <- read.csv("e:/peter/sppweb2016/round01/tables/lichens-paspen.csv")[,1:2]
y$Species <- tax1$scinam[match(y$Sp, tax1$sppid)]
all(as.character(y$Species) == as.character(x$Species))
xx <- x[,c("Productive", "Productive.LCI", "Productive.UCI",
"Clay", "Clay.LCI", "Clay.UCI", "Saline", "Saline.LCI", "Saline.UCI",
"RapidDrain", "RapidDrain.LCI", "RapidDrain.UCI", "Cult", "Cult.LCI",
"Cult.UCI", "UrbInd", "UrbInd.LCI", "UrbInd.UCI")]
xx <- data.frame(Species=x$Species, plogis(qlogis(as.matrix(xx))+y$pAspen))
write.csv(xx, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/lichens-s-paspen.csv")

yy <- as.matrix(x[,-1])
yy <- yy[,!grepl(".LCI", colnames(yy))]
yy <- yy[,!grepl(".UCI", colnames(yy))]
yy <- yy[,!(colnames(yy) %in% c("SoftLin", "HardLin"))]
yyy <- data.frame(Species=x$Species, AverageCoef=rowMeans(yy, na.rm=TRUE), x[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
write.csv(yyy, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/lichens-lin2south.csv")


## vplants
tax1 <- read.csv("e:/peter/sppweb2016/round01/tables/vplants-tax-eta.csv")
tax2 <- read.csv("e:/peter/AB_data_v2016/data/species/OUT_VPlants_Species_Taxa_2016-08-18.csv")

compare_sets(tax1$scinam, tax2$SPECIES_OLD)
tax1$TSNID <- tax2$TSN_ID[match(tax1$scinam, tax2$SPECIES_OLD)]
write.csv(tax1, row.names=FALSE, na="", file="e:/peter/sppweb2016/round01/tables/vplants-tax-eta-tsn.csv")

y <- read.csv("e:/peter/sppweb2016/round01/tables/vplants-veghf.csv")
z <- data.frame(Species=y$Species)
for (i in 1:nrow(v)) {
    if (as.character(v$eta)[i] != "") {
        z[[as.character(v$ref[i])]] <- y[[as.character(v$eta)[i]]]
    } else {
        z[[as.character(v$ref[i])]] <- NA
    }
}
write.csv(z, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/vplants-veghf2.csv")

yy <- as.matrix(y[,-1])
yy <- yy[,!grepl(".LCI", colnames(yy))]
yy <- yy[,!grepl(".UCI", colnames(yy))]
yy <- yy[,!(colnames(yy) %in% c("SoftLin", "HardLin"))]
yyy <- data.frame(Species=y$Species, AverageCoef=rowMeans(yy, na.rm=TRUE), y[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
write.csv(yyy, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/vplants-lin2.csv")

x <- read.csv("e:/peter/sppweb2016/round01/tables/vplants-soil.csv")
y <- read.csv("e:/peter/sppweb2016/round01/tables/vplants-paspen.csv")[,1:2]
y$Species <- tax1$scinam[match(y$Sp, tax1$sppid)]
y <- y[!is.na(y$Species),]
all(as.character(y$Species) == as.character(x$Species))
xx <- x[,c("Productive", "Productive.LCI", "Productive.UCI",
"Clay", "Clay.LCI", "Clay.UCI", "Saline", "Saline.LCI", "Saline.UCI",
"RapidDrain", "RapidDrain.LCI", "RapidDrain.UCI", "Cult", "Cult.LCI",
"Cult.UCI", "UrbInd", "UrbInd.LCI", "UrbInd.UCI")]
xx <- data.frame(Species=x$Species, plogis(qlogis(as.matrix(xx))+y$pAspen))
write.csv(xx, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/vplants-s-paspen.csv")

yy <- as.matrix(x[,-1])
yy <- yy[,!grepl(".LCI", colnames(yy))]
yy <- yy[,!grepl(".UCI", colnames(yy))]
yy <- yy[,!(colnames(yy) %in% c("SoftLin", "HardLin"))]
yyy <- data.frame(Species=x$Species, AverageCoef=rowMeans(yy, na.rm=TRUE), x[,c("SoftLin","SoftLin.LCI","SoftLin.UCI","HardLin","HardLin.LCI","HardLin.UCI")])
write.csv(yyy, row.names=FALSE, na = "", file="e:/peter/sppweb2016/round01/tables/vplants-lin2south.csv")

## GoF AUC ROC R2 and stuff

#library(pROC)
library(ResourceSelection)

pr_fun_for_gof <- function(est, X, off=0) {
    if (is.null(dim(est))) {
        mu0 <- drop(X %*% est)
        exp(mu0 + off)
    } else {
        mu0 <- apply(est, 1, function(z) X %*% z)
        rowMeans(exp(mu0 + off))
    }
}
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}
#system.time(roc12 <- roc(y10sp[ss], pr12[ss]))
#system.time(roczz <- simple_roc(y10sp[ss], pr12[ss]))
#as.numeric(auc(roc12))
#simple_auc(roczz)
#plot(roc12)
#lines(1-roczz[,"FPR"], roczz[,"TPR"], col=2)


en <- new.env()
load(file.path(ROOT, "data", "data-north.Rdata"), envir=en)
OFFn <- en$OFF
OFFmn <- en$OFFmean
BBn <- en$BB
DATn <- en$DAT
rm(en)
## out-of-sample set
bunique <- unique(BBn)
INTERNAL <- 1:nrow(DATn) %in% bunique
ss <- which(INTERNAL)
ss1 <- which(!INTERNAL)
bid <- DATn$bootid
levels(bid) <- sapply(strsplit(levels(bid), "\\."), "[[", 1)

spp <- "OVEN"
all_acc_N <- list()
for (spp in fln) {
cat(spp, "N\n");flush.console()

## Deviance based pseudo R^2 (only internal)
y1sp <- yyn[,spp]
y10sp <- ifelse(y1sp > 0, 1, 0)
off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
resn <- loadSPP(file.path(ROOT, "results", "north", paste0("birds_abmi-north_", spp, ".Rdata")))
est5 <- getEst(resn, stage=5, na.out=FALSE, Xnn)
est6 <- getEst(resn, stage=6, na.out=FALSE, Xnn)
## 1st run
pr11 <- pr_fun_for_gof(est5[1,], Xnn, off=off1sp)
pr12 <- pr_fun_for_gof(est6[1,], Xnn, off=off1sp)
## Null model: intercept and offset
ll0 <- sum(dpois(y1sp[ss], mean(y1sp[ss]), log=TRUE))
## Saturated: one parameter per observation
lls <- sum(dpois(y1sp[ss], y1sp[ss], log=TRUE))
## Full: our 1st prediction
ll11 <- sum(dpois(y1sp[ss], pr11[ss], log=TRUE))
ll12 <- sum(dpois(y1sp[ss], pr12[ss], log=TRUE))
R211 <- 1-(lls-ll11)/(lls-ll0)
R212 <- 1-(lls-ll12)/(lls-ll0)

## ROC
#roc11 <- roc(y10sp[ss], pr11[ss])
#roc12 <- roc(y10sp[ss], pr12[ss])
#auc11 <- as.numeric(roc11$auc)
#auc12 <- as.numeric(roc12$auc)
roc11 <- simple_roc(y10sp[ss], pr11[ss])
roc12 <- simple_roc(y10sp[ss], pr12[ss])
auc11 <- simple_auc(roc11)
auc12 <- simple_auc(roc12)
spp1res <- c(R211=R211, R212=R212,
    AUC11=auc11, AUC12=auc12, prevalence=mean(y10sp))
all_acc_N[[spp]] <- spp1res
}
all_acc_N <- do.call(rbind, all_acc_N)

es <- new.env()
load(file.path(ROOT, "data", "data-south.Rdata"), envir=es)
OFFs <- es$OFF
OFFms <- es$OFFmean
BBs <- es$BB
DATs <- es$DAT
rm(es)
## out-of-sample set
bunique <- unique(as.numeric(BBs))
INTERNAL <- 1:nrow(DATs) %in% bunique
ss <- which(INTERNAL)
ss1 <- which(!INTERNAL)
bid <- DATs$bootid
levels(bid) <- sapply(strsplit(levels(bid), "\\."), "[[", 1)

all_acc_S <- list()
for (spp in fls) {
cat(spp, "S\n");flush.console()

## Deviance based pseudo R^2 (only internal)
y1sp <- yys[,spp]
y10sp <- ifelse(y1sp > 0, 1, 0)
off1sp <- if (spp %in% colnames(OFFs)) OFFs[,spp] else OFFms
ress <- loadSPP(file.path(ROOT, "results", "south", paste0("birds_abmi-south_", spp, ".Rdata")))

est5 <- getEst(ress, stage=3, na.out=FALSE, Xns)
est6 <- getEst(ress, stage=4, na.out=FALSE, Xns)
## 1st run
pr11 <- pr_fun_for_gof(est5[1,], Xns, off=off1sp)
pr12 <- pr_fun_for_gof(est6[1,], Xns, off=off1sp)
## Null model: intercept and offset
ll0 <- sum(dpois(y1sp[ss], mean(y1sp[ss]), log=TRUE))
## Saturated: one parameter per observation
lls <- sum(dpois(y1sp[ss], y1sp[ss], log=TRUE))
## Full: our 1st prediction
ll11 <- sum(dpois(y1sp[ss], pr11[ss], log=TRUE))
ll12 <- sum(dpois(y1sp[ss], pr12[ss], log=TRUE))
R211 <- 1-(lls-ll11)/(lls-ll0)
R212 <- 1-(lls-ll12)/(lls-ll0)

## ROC
#roc11 <- roc(y10sp[ss], pr11[ss])
#roc12 <- roc(y10sp[ss], pr12[ss])
#auc11 <- as.numeric(roc11$auc)
#auc12 <- as.numeric(roc12$auc)
roc11 <- simple_roc(y10sp[ss1], pr11[ss1])
roc12 <- simple_roc(y10sp[ss1], pr12[ss1])
auc11 <- simple_auc(roc11)
auc12 <- simple_auc(roc12)
spp1res <- c(R211=R211, R212=R212,
    AUC11=auc11, AUC12=auc12, prevalence=mean(y10sp))
all_acc_S[[spp]] <- spp1res

}
all_acc_S <- do.call(rbind, all_acc_S)

## OCCC metrics
library(epiR)
occc_res_N <- list()
boot_res_N <- list()
for (spp in fln) {
    cat(spp, "N\n");flush.console()
    resn <- loadSPP(file.path(ROOT, "results", "north",
        paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    pr <- exp(pred_veghf(estn_hab, Xnn, burn_included=FALSE, raw=TRUE))
    occc_res_N[[spp]] <- epi.occc(pr)
    boot_res_N[[spp]] <- pr
}
occc_res_S <- list()
boot_res_S <- list()
for (spp in fls) {
    cat(spp, "S\n");flush.console()
    ress <- loadSPP(file.path(ROOT, "results", "south",
        paste0("birds_abmi-south_", spp, ".Rdata")))
    ests_hab <- getEst(ress, stage=stage_hab_s, na.out=FALSE, Xns)
    pr <- exp(pred_soilhf(ests_hab, Xns, raw=TRUE))
    occc_res_S[[spp]] <- epi.occc(pr)
    boot_res_S[[spp]] <- pr
}

save(all_acc_S, all_acc_N, occc_res_N, boot_res_N, occc_res_S, boot_res_S, slt,
    file="e:/peter/sppweb2017/birds-r2-auc-occc-boot.Rdata")
if (FALSE) {

ss <- as.character(slt[slt$veghf.north,"AOU"])
ren <- as.character(slt[slt$veghf.north,"sppid"])
all_acc_N <- all_acc_N[ss,]
rownames(all_acc_N) <- ren
occc_res_N <- occc_res_N[ss]
names(occc_res_N) <- ren
boot_res_N <- boot_res_N[ss]
names(boot_res_N) <- ren

ss <- as.character(slt[slt$soilhf.south,"AOU"])
ren <- as.character(slt[slt$soilhf.south,"sppid"])
all_acc_S <- all_acc_S[ss,]
rownames(all_acc_S) <- ren
occc_res_S <- occc_res_S[ss]
names(occc_res_S) <- ren
boot_res_S <- boot_res_S[ss]
names(boot_res_S) <- ren

save(all_acc_S, all_acc_N, occc_res_N, boot_res_N, occc_res_S, boot_res_S, slt,
    file="e:/peter/sppweb2017/birds-r2-auc-occc-boot-subset.Rdata")

}

SPPn <- as.character(slt[slt$veghf.north,"AOU"])
SPPs <- as.character(slt[slt$soilhf.south,"AOU"])
load("e:/peter/sppweb2017/birds-r2-auc-occc-boot.Rdata")

z <- all_acc_N
SPP <- SPPn
with(data.frame(z), plot(R211, R212, type="n", ylim=c(0,1), xlim=c(0,1)))
abline(0,1, col="grey")
with(data.frame(z[rownames(z) %in% SPP,]), points(R211, R212))
with(data.frame(z[!(rownames(z) %in% SPP),]), points(R211, R212, col=2))

load("e:/peter/sppweb2017/birds-r2-auc-occc-boot-subset.Rdata")
occn <- t(sapply(occc_res_N, "[", 1:3))
occs <- t(sapply(occc_res_S, "[", 1:3))

par(mfrow=c(3,2))
with(data.frame(all_acc_N), plot(R211, R212, ylim=c(0,1), xlim=c(0,1),
    main="R^2, North", xlab="Habitat", ylab="Habitat + Climate"))
abline(0,1, col="grey")
with(data.frame(all_acc_S), plot(R211, R212, ylim=c(0,1), xlim=c(0,1),
    main="R^2, South", xlab="Habitat", ylab="Habitat + Climate"))
abline(0,1, col="grey")

with(data.frame(all_acc_N), plot(AUC11, AUC12, ylim=c(0,1), xlim=c(0,1),
    main="AUC, North", xlab="Habitat", ylab="Habitat + Climate"))
abline(0,1, col="grey")
abline(h=0.5, v=0.5, col="grey", lty=2)
with(data.frame(all_acc_S), plot(AUC11, AUC12, ylim=c(0,1), xlim=c(0,1),
    main="AUC, South", xlab="Habitat", ylab="Habitat + Climate"))
abline(0,1, col="grey")
abline(h=0.5, v=0.5, col="grey", lty=2)

with(data.frame(occn), plot(oprec, oaccu , ylim=c(0,1), xlim=c(0,1),
    main="AUC, North", xlab="Location", ylab="Scale"))
abline(0,1, col="grey")
with(data.frame(occs), plot(oprec, oaccu , ylim=c(0,1), xlim=c(0,1),
    main="AUC, North", xlab="Location", ylab="Scale"))
abline(0,1, col="grey")





## means
flam <- function(ss1) {
    c(lam_obs=mean(y1sp[ss1]),
    lam_est=mean(prf2[ss1]),
    p_obs=mean(y10sp[ss1]),
    p_est=mean(1-exp(-prf2[ss1])))
}
lam_all <- list()
for (spp in fln) {
cat(spp, "\n");flush.console()
y1sp <- yyn[,spp]
y10sp <- ifelse(y1sp > 0, 1, 0)
off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
resn <- loadSPP(file.path(ROOT, "results", "north", paste0("birds_abmi-north_", spp, ".Rdata")))
#est5 <- getEst(resn, stage=5, na.out=FALSE, Xnn)
est6 <- getEst(resn, stage=6, na.out=FALSE, Xnn)
#prf1 <- pr_fun_for_gof(est5, Xnn, off=off1sp)
prf2 <- pr_fun_for_gof(est6, Xnn, off=off1sp)

## regional ROC analysis
regres <- list()
#REG <- "Foothills"
for (REG in levels(bid)) {
    regres[[REG]] <- flam(bid == REG & !INTERNAL)
}
regres <- do.call(rbind, regres)
regres <- regres[c("NW", "NE", "SW", "SE", "Foothills", "Parkland", "Rocky Mountain"),]

lam_all[[spp]] <- rbind(All=flam(ss1), regres)
}

save(all_acc, occc_res, lam_all, file=file.path(ROOT, "tables", "res_acc.Rdata"))

load(file=file.path(ROOT, "tables", "res_acc.Rdata"))

pdet <- sapply(fln, function(z) sum(yyn[ss,z]>0)/length(ss))
overall <- t(sapply(all_acc, function(z) z$overall))
kreg <- t(sapply(all_acc, function(z) z$regions[,"kf2"]))
occc <- t(sapply(occc_res, function(z) unlist(z[1:3])))
lam <- t(sapply(lam_all, function(z) z["All", 1:2]))

blt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(blt) <- blt$AOU
nmok <- blt$map.pred[match(fln, blt$AOU)]
sing <- blt$singing[match(fln, blt$AOU)]
names(sing) <- names(nmok) <- fln
table(sing, nmok)
typ <- factor(rep("OK", length(SPP)), c("OK", "Prblm", "NSng"))
names(typ) <- SPP
typ[!sing] <- "NSng"
typ[sing & !nmok] <- "Prblm"

plot(overall[,c("R2f2", "R212")], xlim=c(0,1), ylim=c(0,1))
abline(0,1)

overall[overall[,"R2f2"]<0,]

plot(overall[,c("R211", "R212")], xlim=c(0,1), ylim=c(0,1),
    cex=1+2*sqrt(pdet), col=as.integer(typ)+1,
    xlab="Pseudo R^2, Habitat", ylab="Pseudo R^2, Habitat+Space")
abline(0,1, lty=2, col="grey")

plot(lam, xlim=c(0,0.5), ylim=c(0,0.5),
    cex=1+2*sqrt(pdet), col=as.integer(typ)+1,
    xlab="Mean Observed Count", ylab="Mean Predicted # Det.")
abline(0,1, lty=2, col="grey")

plot(pdet*length(ss), overall[,"R212"], xlim=c(0,200))

plot(overall[,c("kf1", "kf2")], xlim=c(-1,1), ylim=c(-1,1),
    cex=1+2*sqrt(pdet), col=as.integer(typ)+1,
    xlab="CAUC, Habitat", ylab="CAUC, Habitat+Space")
abline(0,1, lty=2, col="grey")
abline(h=0,v=0, lty=3, col="grey")

plot(overall[,c("AUCf1", "AUCf2")], xlim=c(0,1), ylim=c(0,1), cex=1+2*sqrt(pdet))
abline(0,1)

boxplot(kreg[rowSums(is.na(kreg))==0,], ylim=c(-1,1), ylab="CAUC, Habitat+Space")
abline(h=0,lty=2, col="grey")

library(plotrix)
ladderplot(kreg[rowSums(is.na(kreg))==0,], ylim=c(-1,1), pch=NA)
abline(h=0)

layout(matrix(c(1,1,1,2,1,1,1,3,1,1,1,4), 3, 4, byrow=TRUE))
plot(occc[,2], occc[,3], xlim=c(0,1), ylim=c(0,1),
    cex=1+2*sqrt(pdet), col=as.integer(typ)+1,
    xlab="Overall Precision", ylab="Overall Accuracy")
abline(0,1, lty=2, col="grey")
hist(occc[,1], xlim=c(0,1), main="Overall Concordance")
hist(occc[,2], xlim=c(0,1), main="Overall Precision")
hist(occc[,3], xlim=c(0,1), main="Overall Accuracy")

occcx <- occc[rowSums(is.na(occc))==0,]
occcx[occcx[,2] < 0.2 & occcx[,3] < 0.2,,drop=FALSE]
occcx[occcx[,2] < 0.4 & occcx[,3] > 0.6,,drop=FALSE]
occcx[occcx[,2] > 0.8 & occcx[,3] < 0.2,,drop=FALSE]

foccc <- function(spp, nn=10) {
    resn <- loadSPP(file.path(ROOT, "results", "north",
        paste0("birds_abmi-north_", spp, ".Rdata")))
    estn_hab <- getEst(resn, stage=stage_hab_n, na.out=FALSE, Xnn)
    pr <- exp(pred_veghf(estn_hab, Xnn, burn_included=FALSE, raw=TRUE))
    ii <- sample(ncol(pr), min(nn, ncol(pr)))

    prr <- pr[order(rowMeans(pr)),ii]
    prr <- t(t(prr)/apply(prr,2,max))
    matplot(prr, type="l", col=sample(rainbow(length(ii))), lty=1,
        axes=FALSE, xlab="Land cover types", ylab="Relative abundance",
        main=blt[spp,"species"], ylim=c(0,1.2))
    text(1, 1.1, paste("Overall Concordance =", round(occc_res[[spp]]$occc, 3),
        "\nOverall Precision =", round(occc_res[[spp]]$oprec, 3),
        "\nOverall Accuracy =", round(occc_res[[spp]]$oaccu, 3)), pos=4)
    box()
}
par(mfrow=c(2,2))
foccc("AMKE")
foccc("BOCH")
foccc("CITE")
foccc("PUMA")





par(mfrow=c(2,3))
ymax <- max(pr0,prf1, prf2)
ResourceSelection:::.mep(y1sp[ss], pr0[ss], link="log", type="unique", level=0.9,
    main="Internal Null", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss], prf1[ss], link="log", type="unique", level=0.9,
    main="Internall Local", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss], prf2[ss], link="log", type="unique", level=0.9,
    main="Internal Space", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss1], pr0[ss1], link="log", type="unique", level=0.9,
    main="External Null", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss1], prf1[ss1], link="log", type="unique", level=0.9,
    main="Externall Local", ylim=c(0,ymax))
ResourceSelection:::.mep(y1sp[ss1], prf2[ss1], link="log", type="unique", level=0.9,
    main="External Space", ylim=c(0,ymax))

#ss <- lapply(1:240, function(i) which(!(1:nrow(DAT) %in% BB[,i])))
#ssd <- data.frame(id=unlist(ss), iter=unlist(lapply(1:240, function(i) rep(i, length(ss[[i]])))))
#ssx <- Xtab(~iter + id, ssd)
#ssi <- which(colSums(ssx) == 240)
#ssc <- as.integer(colnames(ssx)[ssi])
#compare_sets(ssc, ss1)

Col <- c("blue","darkgreen","red")

pdf(file.path(OUTDIR, "cawa-fig-roc.pdf"))
op <- par(las=1)
plot(Show[[1]], col=Col[1], lty=1)
for (i in 2:length(Show))
    lines(Show[[i]], col=Col[i], lty=1)
aucx <- sapply(Show, function(z) as.numeric(z$auc))
#txt <- paste0(names(aucx), " (AUC = ", round(aucx, 3), ")")
txt <- paste0(c("Local (stages 1-6)",
    "Stand level (stages 1-9)","Year (stages 1-10)"), " (AUC = ", round(aucx, 3), ")")
legend("bottomright", bty="n", col=rev(Col),
    lty=1, lwd=2, legend=rev(txt))
dev.off()
png(file.path(OUTDIR, "cawa-fig-roc.png"))
op <- par(las=1)
plot(Show[[1]], col=Col[1], lty=1)
for (i in 2:length(Show))
    lines(Show[[i]], col=Col[i], lty=1)
aucx <- sapply(Show, function(z) as.numeric(z$auc))
#txt <- paste0(names(aucx), " (AUC = ", round(aucx, 3), ")")
txt <- paste0(c("Local","Spatial","Year"), " (AUC = ", round(aucx, 3), ")")
legend("bottomright", bty="n", col=rev(Col),
    lty=1, lwd=2, legend=rev(txt))
dev.off()

## mean plots

REGNAMS <- substr(rownames(lam_all[[1]])[-1], 1, pmin(nchar(rownames(lam_all[[1]])[-1]), 5))
pdf(file.path("gof-figures.pdf"), width=12, height=6, onefile=TRUE)
for (spp in fln) {
cat(spp, "\n");flush.console()
op <- par(las=1, mfrow=c(1,2))
MAX <- max(lam_all[[spp]][,1:2])*1.2
plot(lam_all[[spp]][-1,1:2], pch=21, ylim=c(0, MAX), xlim=c(0, MAX),
    col=4, cex=2, main=blt[spp,"species"],
    xlab="Mean Observed Count", ylab="Mean Expected Number of Detections")
abline(0,1)
abline(v=lam_all[[spp]][1,1], h=lam_all[[spp]][1,2], col=4)
text(lam_all[[spp]][-1,1], lam_all[[spp]][-1,2],
    round(all_acc[[spp]]$regions[,"kf2"], 2), pos=4, cex=0.6)
text(lam_all[[spp]][-1,1], lam_all[[spp]][-1,2], REGNAMS, pos=3, cex=0.4)
STAT <- all_acc[[spp]]$overall[c("R2f2", "kf2")]
text(0, MAX*0.9, paste0("Deviance R^2 = ", max(0, round(STAT[1],3)),
    "\nCorrected AUC = ", round(STAT[2],3)), pos=4)
foccc(spp, 25)
par(op)
}
dev.off()

## estimating overall ARU effect

table(xnn$ARU3)

aru <- list()
for (spp in colnames(OFFn)) {
    cat(spp, "\n");flush.console()
    y1sp <- yyn[,spp]
    off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
    m0 <- glm(y1sp ~ ARU3-1, xnn, family=poisson, offset=off1sp)
    aru[[spp]] <- exp(coef(m0))
}
aru <- do.call(rbind, aru)

nn <- list()
for (spp in colnames(OFFn)) {
    cat(spp, "\n");flush.console()
    nn[[spp]] <- table(xnn$ARU3, factor(ifelse(yyn[,spp]>0, 1, 0), 0:1))[,"1"]
}
nn <- do.call(rbind, nn)
smin <- apply(nn, 1, min)

boxplot(aru[smin >= 5, ])
plotrix::ladderplot(aru[smin >= 5, ], pch=NA)

MIN <- 10
plot(density(aru[smin >= MIN,"ARU3SM"]/aru[smin >= MIN,"ARU3RF"]))
summary(aru[smin >= MIN,"ARU3SM"]/aru[smin >= MIN,"ARU3RF"])

## this is too small -- use SVW and YIP

## Geo XV geographic cross valudation

en <- new.env()
load(file.path(ROOT, "data", "data-north-geoxv.Rdata"), envir=en)
OFFn <- en$OFF
OFFmn <- en$OFFmean
BBn <- en$BB
DATn <- en$DAT
yyn <- en$YY
SPP <- colnames(yyn)
modsn <- en$mods

nTerms <- getTerms(modsn, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), DATn)
colnames(Xnn) <- fixNames(colnames(Xnn))

rm(en)

pr_fun_for_geoxv <-
function(est, X, off=0) {
    mu0 <- apply(est, 1, function(z) X %*% z)
    exp(mu0 + off)
}

spp <- "OVEN"
geoxv_N <- list()
for (spp in SPP) {
cat(spp, "N\n");flush.console()

    y1sp <- yyn[,spp]
    y10sp <- ifelse(y1sp > 0, 1, 0)
    off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
    resn <- loadSPP(file.path(ROOT, "results", "north-geoxv",
        paste0("birds_abmi-north-geoxv_", spp, ".Rdata")))
    resn[[1]] <- NULL # iteration 0 (full)
    est6 <- getEst(resn, stage=6, na.out=FALSE, Xnn)
    ## 1st run
    pr <- pr_fun_for_geoxv(est6, Xnn, off=off1sp)
    auc <- matrix(NA, ncol(pr), 4)
    dimnames(auc) <- list(1:ncol(pr), c("n_in", "n_out", "AUC_in", "AUC_out"))

    for (jj in 1:ncol(pr)) {
        INTERNAL <- BBn[,1] %in% BBn[BBn[,2]!=jj,1]
        ss <- which(INTERNAL)
        ss1 <- which(!INTERNAL)
        roc_in <- simple_roc(y10sp[ss], pr[ss,jj])
        roc_out <- simple_roc(y10sp[ss1], pr[ss1,jj])
        auc[jj,"n_in"] <- sum(y10sp[ss])
        auc[jj,"n_out"] <- sum(y10sp[ss1])
        auc[jj,"AUC_in"] <- simple_auc(roc_in)
        auc[jj,"AUC_out"] <- simple_auc(roc_out)
    }
    geoxv_N[[spp]] <- auc
}

en <- new.env()
load(file.path(ROOT, "data", "data-south-geoxv.Rdata"), envir=en)
OFFn <- en$OFF
OFFmn <- en$OFFmean
BBn <- en$BB
DATn <- en$DAT
yyn <- en$YY
SPP <- colnames(yyn)
modsn <- en$mods

nTerms <- getTerms(modsn, "list")
Xnn <- model.matrix(getTerms(modsn, "formula"), DATn)
colnames(Xnn) <- fixNames(colnames(Xnn))

rm(en)

#spp <- "OVEN"
geoxv_S <- list()
for (spp in SPP) {
cat(spp, "S\n");flush.console()

    y1sp <- yyn[,spp]
    y10sp <- ifelse(y1sp > 0, 1, 0)
    off1sp <- if (spp %in% colnames(OFFn)) OFFn[,spp] else OFFmn
    resn <- loadSPP(file.path(ROOT, "results", "south-geoxv",
        paste0("birds_abmi-north-geoxv_", spp, ".Rdata")))
    resn[[1]] <- NULL # iteration 0 (full)
    est6 <- getEst(resn, stage=4, na.out=FALSE, Xnn)
    ## 1st run
    pr <- pr_fun_for_geoxv(est6, Xnn, off=off1sp)
    auc <- matrix(NA, ncol(pr), 4)
    dimnames(auc) <- list(1:ncol(pr), c("n_in", "n_out", "AUC_in", "AUC_out"))

    for (jj in 1:ncol(pr)) {
        INTERNAL <- BBn[,1] %in% BBn[BBn[,2]!=jj,1]
        ss <- which(INTERNAL)
        ss1 <- which(!INTERNAL)
        roc_in <- simple_roc(y10sp[ss], pr[ss,jj])
        roc_out <- simple_roc(y10sp[ss1], pr[ss1,jj])
        auc[jj,"n_in"] <- sum(y10sp[ss])
        auc[jj,"n_out"] <- sum(y10sp[ss1])
        auc[jj,"AUC_in"] <- simple_auc(roc_in)
        auc[jj,"AUC_out"] <- simple_auc(roc_out)
    }
    geoxv_S[[spp]] <- auc
}


gxn <- do.call(rbind, lapply(1:length(geoxv_N), function(i) {
    data.frame(Species=names(geoxv_N)[i], Cluster=rownames(geoxv_N[[i]]), geoxv_N[[i]])
}))
gxnin <- data.frame(What=factor("In", c("In", "Out")), gxn[,c(1,2,3,5)])
gxnout <- data.frame(What=factor("Out", c("In", "Out")), gxn[,c(1,2,4,6)])
colnames(gxnin) <- colnames(gxnout) <- c("What", "Species", "Cluster", "n", "AUC")
gxn <- rbind(gxnin, gxnout)
gxn$Cross <- interaction(gxn$What, gxn$Cluster)
gxn$AUC[gxn$n < 3] <- -1
mln <- as.matrix(Xtab(AUC ~ Species + Cross, gxn))
mln[mln < 0] <- NA
gxn$AUC[gxn$n < 3] <- NA

gxs <- do.call(rbind, lapply(1:length(geoxv_S), function(i) {
    data.frame(Species=names(geoxv_S)[i], Cluster=rownames(geoxv_S[[i]]), geoxv_S[[i]])
}))
gxsin <- data.frame(What=factor("In", c("In", "Out")), gxs[,c(1,2,3,5)])
gxsout <- data.frame(What=factor("Out", c("In", "Out")), gxs[,c(1,2,4,6)])
colnames(gxsin) <- colnames(gxsout) <- c("What", "Species", "Cluster", "n", "AUC")
gxs <- rbind(gxsin, gxsout)
gxs$Cross <- interaction(gxs$What, gxs$Cluster)
gxs$AUC[gxs$n < 3] <- -1
mls <- as.matrix(Xtab(AUC ~ Species + Cross, gxs))
mls[mls < 0] <- NA
gxs$AUC[gxs$n < 3] <- NA

par(mfrow=c(2,1),las=1)
boxplot(AUC ~ Cross, gxn, range=0, ylab="AUC", ylim=c(0,1),
    main="Birds, Geographical Cross-validation, North")
matlines(1:nlevels(gxn$Cross), t(mln), col=rgb(0,0,0,0.25), lty=1)
boxplot(AUC ~ Cross, gxn, range=0, add=TRUE, col=rgb(1,0.5,0.5,0.4))
abline(h=0.5, lty=2, col="grey")
boxplot(AUC ~ Cross, gxs, range=0, ylab="AUC", ylim=c(0,1),
    main="Birds, Geographical Cross-validation, South")
matlines(1:nlevels(gxs$Cross), t(mls), col=rgb(0,0,0,0.25), lty=1)
boxplot(AUC ~ Cross, gxs, range=0, add=TRUE, col=rgb(1,0.5,0.5,0.4))
abline(h=0.5, lty=2, col="grey")

save(geoxv_N, geoxv_S, gxn, mln, gxs, mls,
    file="e:/peter/sppweb2017/birds-auc-geoxv.Rdata")


