## Sector effects

## check partial backill classes to match sectors used here!!!


library(mefa4)
library(intrval)
library(raster)
library(rgdal)
library(sp)

INDIR <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/sector")

sectors_all <- c("Agr", "Energy", "For", "Urban","Transp")

load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$X <- kgrid$POINT_X
kgrid$Y <- kgrid$POINT_Y


## ch2soil ch2veg trSoil trVeg
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v6hf2016v3noDistVeg.Rdata")
stopifnot(all(rownames(kgrid) == rownames(trVeg)))
stopifnot(all(rownames(kgrid) == rownames(trSoil)))

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]
tv <- droplevels(tv[!endsWith(rownames(tv), "0"),])
tv$ao <- as.factor(paste0(as.character(tv[, "UseInAnalysisFine"]), ifelse(tv[, "MatureOld"], "O", "")))

compare_sets(ch2veg$cr, rownames(tv))
setdiff(ch2veg$cr, rownames(tv))
setdiff(rownames(tv), ch2veg$cr)

ch2veg$rf2 <- tv$UseInAnalysisFineAge[match(ch2veg$rf, rownames(tv))]
ch2veg$cr2 <- tv$UseInAnalysisFineAge[match(ch2veg$cr, rownames(tv))]
ch2veg$sector <- tv$Sector61[match(ch2veg$cr, rownames(tv))]
ch2veg$rf3 <- tv$ao[match(ch2veg$rf, rownames(tv))]
ch2veg$cr3 <- tv$ao[match(ch2veg$cr, rownames(tv))]

with(ch2veg, table(cr3, sector))

EXCL <- c("HWater", "SoilUnknown", "SoilWater", "Water")
MODIF <- c("SoftLin", "Well", "EnSoftLin", "TrSoftLin", "Seismic")

keepn <- rownames(ch2veg)[!(ch2veg$cr2 %in% EXCL) & !(ch2veg$rf2 %in% EXCL)]
trVeg <- trVeg[,keepn]
ch2veg <- ch2veg[keepn,]
ch2veg$modif <- ch2veg$cr2 %in% MODIF
rsn <- rowSums(trVeg)
rsn[rsn==0] <- 1
trVeg <- trVeg / rsn


## define OSR
ogrListLayers("d:/spatial/Oilsands-Boundaries.gdb")
pl <- readOGR("d:/spatial/Oilsands-Boundaries.gdb", "OilsandRegionDissolve10TM")
xy <- kgrid[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- ~ POINT_X + POINT_Y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(pl))
o <- over(xy, pl)
#plot(xy, pch=".", col=ifelse(is.na(o$FIELDCODE), 1, 4))

ss <- !is.na(o$FIELDCODE)
trVegSS <- trVeg[ss,]
AVegSS <- colSums(trVegSS)

levs <- c(
    "MineSite",
    "Pipeline",
    "TransmissionLine",
    "RailHardSurface",
    "RailVegetatedVerge",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge",
    "SeismicLineNarrow",
    "SeismicLineWide",
    "IndustrialSiteRural",
    "UrbanIndustrial",
    "WellSite")

sect <- as.character(ch2veg$sector)
sect[ch2veg$cr %in% levs] <- as.character(ch2veg$cr)[ch2veg$cr %in% levs]
sect[sect %in% c("SeismicLineNarrow", "SeismicLineWide")] <- "SeismicLine"
sect[sect %in% c("Pipeline", "TransmissionLine")] <- "PipeTransLine"
sect[sect %in% c("RailHardSurface",
    "RailVegetatedVerge",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge")] <- "RoadRailVerge"
sect[sect %in% c("IndustrialSiteRural", "UrbanIndustrial")] <- "Industrial"
sect[sect == "Energy"] <- "OtherEnergy"
data.frame(x=table(sect))


SPP <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATTW", "BAOR", "BARS",
    "BAWW", "BBMA", "BBWA", "BBWO", "BCCH", "BHCO", "BHVI", "BLJA",
    "BLPW", "BOCH", "BRBL", "BRCR", "BTNW", "CAWA", "CCSP", "CEDW",
    "CHSP", "CMWA", "CONW", "CORA", "COYE", "DEJU", "DOWO", "EAKI",
    "EVGR", "FOSP", "GCKI", "GRAJ", "GRCA", "HAWO", "HETH", "HOLA",
    "HOWR", "LCSP", "LEFL", "LISP", "MAWA", "MODO", "MOWA", "NOFL",
    "NOWA", "OCWA", "OSFL", "OVEN", "PAWA", "PHVI", "PISI", "PIWO",
    "PUFI", "RBGR", "RBNU", "RCKI", "REVI", "RUBL", "RWBL", "SAVS",
    "SOSP", "SWSP", "SWTH", "TEWA", "TRES", "VATH", "VEER", "VESP",
    "WAVI", "WBNU", "WCSP", "WETA", "WEWP", "WIWA", "WIWR", "WTSP",
    "WWCR", "YBFL", "YBSA", "YEWA", "YRWA")

save(SPP, AVegSS, ch2veg,
    file=paste0(INDIR, "/sector-tools.RData"))

## summarize species

library(mefa4)
library(intrval)
library(ggplot2)

INDIR <- paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/sector")
#INDIR <- paste0("~/GoogleWork/tmp/sector")
load(paste0(INDIR, "/sector-tools.RData"))

sectors_list <- list(
    Agr=c("Crop", "RoughP", "TameP"),
    Transp=c("HardLin", "TrSoftLin"),
    EnSoft=c("EnSoftLin", "Seismic"),
    EnHard=c("Mine", "Well", "Industrial"),
    Urb=c("Rural", "Urban"),
    For=c("ForHarv"))
ch2veg$sector <- as.character(ch2veg$sector)
ch2veg$sector2 <- as.character(ch2veg$sector)
ch2veg$sector2[ch2veg$cr3 %in% sectors_list$EnSoft] <- "EnSoft"
ch2veg$sector2[ch2veg$cr3 %in% sectors_list$EnHard] <- "EnHard"
## note: predictions treated Rural/Industrial as Energy and not RuralUrban, peat mine also Energy
ch2veg$sector[ch2veg$sector=="Misc" & ch2veg$sector2=="EnHard"] <- "Energy"
ch2veg$sector[ch2veg$sector=="RuralUrban" & ch2veg$sector2=="EnHard"] <- "Energy"
table(ch2veg$sector2, ch2veg$sector)


A <- groupSums(matrix(AVegSS, ncol=1), 1, ch2veg$sector)
A <- A[,1]/sum(A)
A2 <- groupSums(matrix(AVegSS, ncol=1), 1, ch2veg$sector2)
A2 <- A2[,1]/sum(A2)

load_spp <- function(spp, EN=FALSE) {
    if (EN) {
        SECS <- c("All", "For", "Transp", "EnS", "EnH", "Urban", "Agr")
        sv <- ch2veg$sector2
    } else {
        SECS <- c("All", "For", "Transp", "Energy", "Urban", "Agr")
        sv <- ch2veg$sector
    }
    out <- list()
    for (SEC in SECS) {
        e <- new.env()
        load(paste0(INDIR, "/", spp, "/", spp, "-HF-", toupper(SEC), ".RData"), envir=e)
        out[[SEC]] <- list(
            HCR = groupSums(e$HABCR, 1, sv),
            HRF = groupSums(e$HABRF, 1, sv))
    }
    out
}

se_fun <- function(i, H) {
    SECS <- names(H)
    Nsect <- list()
    for (SEC in SECS) {
        Nsect[[SEC]] <- cbind(
            cr=H[[SEC]]$HCR[,i],
            rf=H[[SEC]]$HRF[,i])
    }
    Nref <- Nsect[[1]][,2]
    Diffs <- sapply(Nsect, function(z) z[,1]-Nref)
    Diffs <- Diffs[,SECS]
    Df2 <- cbind(Native=0,Diffs[,-1])
    Df2 <- Df2[,c("Native", SECS[-1])]
    Tots <- data.frame(#All=Diffs[,1],
        Joint=Diffs[rownames(Diffs) != "Misc","All"], # joint
        Direct=diag(Df2), # pbf direct (pbf=partial backfilled)
        Spillover=colSums(Df2)-diag(Df2)) # indirect effects: pbf all - pbf direct
    Tots <- 100*Tots/sum(Nref)
    Tots$Indirect <- Tots[,"Joint"]-Tots[,"Direct"]
    as.matrix(Tots)
}

sector <- function(spp, level=0.95, EN=FALSE) {
    a <- c((1-level)/2, 1-((1-level)/2))
    H <- load_spp(spp, EN=EN)
    B <- ncol(H[[1]][[1]])
    tmp <- se_fun(1, H)
    S <- array(0, c(dim(tmp), B))
    dimnames(S) <- list(rownames(tmp), colnames(tmp), NULL)
    for (i in seq_len(B))
        S[,,i] <- se_fun(i, H)
    OUT <- list(estimate=tmp, lower=tmp, upper=tmp)
    for (j in seq_len(nrow(tmp))) {
        for (k in seq_len(ncol(tmp))) {
            z <- quantile(S[j,k,], c(0.5, a))
            OUT$estimate[j,k] <- z[1]
            OUT$lower[j,k] <- z[2]
            OUT$upper[j,k] <- z[3]
        }
    }
    OUT
}

sector_res <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    sector_res[[spp]] <- list(coarse=sector(spp), fine=sector(spp, EN=TRUE))
}
save(sector_res, file=paste0(INDIR, "/sector-res.RData"))


load(paste0(INDIR, "/sector-res.RData"))

v1 <- lapply(1:6, function(i)
    t(sapply(sector_res, function(z) z$coarse$estimate[i,1:2])))
names(v1) <- rownames(sector_res[[1]]$coarse$estimate)
v2 <- lapply(1:7, function(i)
    t(sapply(sector_res, function(z) z$fine$estimate[i,1:2])))
names(v2) <- rownames(sector_res[[1]]$fine$estimate)

pdf(paste0(INDIR, "/sector-all.pdf"), width=12, height=8)
op <- par(mfrow=c(2,3))
v <- v1
hist(v[[1]][,1], main="Native", sub="", breaks=20, xlim=c(-50,250),
    xlab="Indirect", col="#00000066")
abline(v=0, col=3, lty=1, lwd=2)
for (i in 2:6) {
    vv <- v[[i]]
    vv <- vv[order(abs(vv[,1] - vv[,2])),]
    lim <- c(-10, 25)
    if (i == 6)
        lim <- c(-50, 300)
    col <- ifelse(vv[,1] > vv[,2], "#ff000066", "#0000ff66")
    col[abs(vv[,1] - vv[,2]) %[]% c(0, 1)] <- "#00000066"
    plot(vv, col=col, main=names(v)[i],
        xlim=lim, ylim=lim, pch=19, cex=1.5)
    abline(0,1, col="grey")
    abline(h=0,v=0, col=3, lty=1)
}
par(op)
dev.off()

pdf(paste0(INDIR, "/sector-all-en.pdf"), width=16, height=8)
op <- par(mfrow=c(2,4))
v <- v2
hist(v[[1]][,1], main="Native", sub="", breaks=20, xlim=c(-50,250),
    xlab="Indirect", col="#00000066")
abline(v=0, col=3, lty=1, lwd=2)
for (i in 2:7) {
    vv <- v[[i]]
    vv <- vv[order(abs(vv[,1] - vv[,2])),]
    lim <- c(-10, 25)
    if (i == 6)
        lim <- c(-50, 300)
    col <- ifelse(vv[,1] > vv[,2], "#ff000066", "#0000ff66")
    col[abs(vv[,1] - vv[,2]) %[]% c(0, 1)] <- "#00000066"
    plot(vv, col=col, main=names(v)[i],
        xlim=lim, ylim=lim, pch=19, cex=1.5)
    abline(0,1, col="grey")
    abline(h=0,v=0, col=3, lty=1)
}
par(op)
dev.off()

v <- NULL
for (i in names(v1)) {
    AA <- if (i == "Native")
        sum(A) else A[i]
    v <- rbind(v, data.frame(Sector=i,
        Species=rownames(v1[[i]]),
        Type=rep(colnames(v1[[i]]), each=nrow(v1[[i]])),
        Total=c(v1[[i]][,1], v1[[i]][,2]),
        Unit=c(v1[[i]][,1], v1[[i]][,2])/AA))
}
v1x <- v
ggplot(v[v$Total <= 25,], aes(x=Sector, y=Total, fill=Type)) +
    geom_boxplot() + geom_abline(intercept = 0, slope=0)
ggsave(paste0(INDIR, "/sector-box-total.pdf"))
ggplot(v[v$Unit %[]% c(-200, 200),], aes(x=Sector, y=Unit, fill=Type)) +
    geom_boxplot() + geom_abline(intercept = 0, slope=0) +
    geom_abline(intercept = c(-100, 100), slope=0, lty=2)
ggsave(paste0(INDIR, "/sector-box-unit.pdf"))

v <- NULL
for (i in names(v2)) {
    AA <- if (i == "Native")
        sum(A2) else A2[i]
    v <- rbind(v, data.frame(Sector=i,
        Species=rownames(v2[[i]]),
        Type=rep(colnames(v2[[i]]), each=nrow(v2[[i]])),
        Total=c(v2[[i]][,1], v2[[i]][,2]),
        Unit=c(v2[[i]][,1], v2[[i]][,2])/AA))
}
v2x <- v
ggplot(v[v$Total <= 25,], aes(x=Sector, y=Total, fill=Type)) +
    geom_boxplot() + geom_abline(intercept = 0, slope=0)
ggsave(paste0(INDIR, "/sector-box-total-en.pdf"))
ggplot(v[v$Unit %[]% c(-200, 200),], aes(x=Sector, y=Unit, fill=Type)) +
    geom_boxplot() + geom_abline(intercept = 0, slope=0) +
    geom_abline(intercept = c(-100, 100), slope=0, lty=2)
ggsave(paste0(INDIR, "/sector-box-unit-en.pdf"))

write.csv(v1x, row.names = FALSE, file=paste0(INDIR, "/sector-results.csv"))
write.csv(v2x, row.names = FALSE, file=paste0(INDIR, "/sector-results-en.csv"))



plot_sect <- function(spp, unit=FALSE, EN=FALSE) {

    z <- if (EN) sector_res[[spp]]$fine else sector_res[[spp]]$coarse
    AA <- if (EN) A2 else A
    z[[1]][1,"Direct"] <- sum(z[[1]][,"Indirect"])
    z[[2]][1,"Direct"] <- sum(z[[2]][,"Indirect"])
    z[[3]][1,"Direct"] <- sum(z[[3]][,"Indirect"])
    if (unit) {
        z[[1]] <- z[[1]]/AA
        z[[2]] <- z[[2]]/AA
        z[[3]] <- z[[3]]/AA
    }
    rr1 <- range(z[[1]][,1:2])
    rr2 <- range(z[[1]][,1:2], z[[2]][,1:2], z[[3]][,1:2])
    rr <- rr2
    rr[1] <- if (abs(rr2[1] - rr1[1])/diff(rr1) > 0.5)
        rr1[1]-0.5*diff(rr1) else rr2[1]
    rr[2] <- if (abs(rr2[2] - rr1[2])/diff(rr1) > 0.5)
        rr1[2]+0.5*diff(rr1) else rr2[2]
    #rr <- rr - (rr %% 10)
    #rr <- rr + c(-10, 10)

    YLAB <- if (unit)
        "Unit effect" else "Regional effect"
    plot(0, type="n", xlim=c(0.5, length(AA)+0.5), ylim=rr, ann=FALSE, axes=FALSE)
    w <- 0.4
    Cols <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')
    if (EN)
        Cols <- Cols[c(1,2,3,4,4,5,6)]
    for (i in 1:length(AA)) {
        COL <- paste0(Cols[i], '66')
        polygon(
            i + c(w, -w, -w, w),
            c(0, 0, z$estimate[i, "Direct"], z$estimate[i, "Direct"]),
            col=paste0(Cols[i], '44'), border=paste0(Cols[i], 'aa'))
        NONSIG <- 0 %[]% c(z$lower[i, "Direct"], z$upper[i, "Direct"])
        lines(
            c(i,i),
            c(z$lower[i, "Direct"], z$upper[i, "Direct"]),
            lwd=if (NONSIG) 1.5 else 3,lend=1,
            col=if (NONSIG) paste0(Cols[i], 'aa') else paste0(Cols[i], 'ff'))
        lines(
            i + c(-w, w),
            c(z$estimate[i, "Direct"], z$estimate[i, "Direct"]),
            lwd=3, lend=1,
            col=if (NONSIG) paste0(Cols[i], 'aa') else paste0(Cols[i], 'ff'))
    }
    title(ylab=YLAB)
    axis(2, col="grey")
    abline(h=0, col="grey")
    if (unit)
        abline(h=c(-100, 100), lty=2)
    LAB <- c("Indirect", "Forestry", "Transp.", "Energy", "Urban", "Agricult.")
    if (EN)
        LAB <- c("Indirect", "Forestry", "Transp.", "En. Lin.", "En. M&W", "Urban", "Agricult.")
    axis(1, 1:length(AA), LAB, tick=FALSE)
    invisible(z)
}

pdf(paste0(INDIR, "/sector-results.pdf"), onefile=TRUE, height=12*0.8, width=16*0.8)
for (spp in SPP) {
    op <- par(mfrow=c(2,2), las=1, mar=c(4,4,4,1))
    plot_sect(spp)
    title(main=spp)
    plot_sect(spp, EN=TRUE)
    plot_sect(spp, TRUE)
    plot_sect(spp, TRUE, EN=TRUE)
    par(op)
}
dev.off()

vvv <- data.frame(sapply(v, function(z) z[,1] - z[,2]))
plot(vvv)

op <- par(mfrow=c(2,3))
for (i in 1:6) {
hist(vvv[,i], main=colnames(vvv)[i], sub="", col="#00000066",
    xlab="Indirect", breaks=20)
abline(v=0, col=3, lty=1, lwd=2)
}
par(op)

## looking at pixel level effects

As <- groupSums(trVegSS, 2, ch2veg$sector)
A2s <- groupSums(trVegSS, 2, ch2veg$sector2)

load_spp_pixel <- function(spp) {
    SECS <- c("All", "For", "Transp", "EnS", "EnH", "Urban", "Agr")
    sv <- ch2veg$sector2
    CN <- c(
        "Agr"="Agriculture",
        "Transp"="Transportation",
        "Energy"="Energy",
        "EnS"="EnSoft",
        "EnH"="EnHard",
        "Urban"="RuralUrban",
        "For"="Forestry")
    CR <- list()
    for (SEC in SECS) {
        e <- new.env()
        load(paste0(INDIR, "/", spp, "/", spp, "-HF-", toupper(SEC), "-pixel.RData"), envir=e)
        CR[[SEC]] <- as.matrix(e$SCR)
        if (SEC == "All") {
            RF <- as.matrix(e$SRF)
            XD <- RF
            XD[] <- 0
        } else {
            XD[,CN[SEC]] <- CR[[SEC]][,CN[SEC]]
        }
    }
    XI <- CR[["All"]] - XD
    X <- XD
    colnames(X)[1L] <- "Indirect"
    X[,"Indirect"] <- rowSums(XI)
    X
}


library(mgcv)
k <- 5

pdf(paste0(INDIR, "/sector-indirect.pdf"), onefile=TRUE, height=8, width=15)
spp <- "OVEN"
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    X <- load_spp_pixel(spp)

    d <- data.frame(Indirect=as.numeric(X[,1]), as.matrix(A2s))
    ss <- sample(which(d$Native < 0.95), 5000)
    sss <- sample(ss, 1000)
    #hist(d[ss,1])
    #plot(d[ss,])

    #library(mgcv)
    q <- quantile(d[ss,1], c(0.01, 0.99))

    m1 <- gam(Indirect ~ s(Native, k=k), data=d[ss,])
    m2 <- gam(Indirect ~ s(Forestry, k=k), data=d[ss,])
    m3 <- gam(Indirect ~ s(Transportation, k=k), data=d[ss,])
    m4 <- gam(Indirect ~ s(EnSoft, k=k), data=d[ss,])
    m5 <- gam(Indirect ~ s(EnHard, k=k), data=d[ss,])
    m6 <- gam(Indirect ~ s(RuralUrban, k=k), data=d[ss,])
    m7 <- gam(Indirect ~ s(Agriculture, k=k), data=d[ss,])
    mm <- list(m1, m2, m3, m4, m5, m6, m7)
    names(mm) <- colnames(d)[-1]
    a <- sapply(mm, AIC)
    a <- a - min(a)
    best <- names(mm)[which.min(a)]

    op <- par(mfrow=c(2,4))
    for (i in colnames(d)[-1]) {
        nd <- data.frame(x=seq(0, quantile(d[,i], 0.99), length.out = 200))
        colnames(nd) <- i
        p <- predict(mm[[i]], newdata=nd)
        plot(d[sss,i], d[sss,1], xlab=i, ylab="Indirect", ylim=q,
            col=if (i == best) "#0000ff22" else "#00000022",
            sub=paste("dAIC =", round(a[i], 2)))
        if (i == "Native")
            title(main=spp)
        abline(h=0,col=1)
        #lines(lowess(d[ss,i], d[ss,1]), col=2, lwd=2, lty=2)
        lines(nd[,1], p, col=2, lwd=2, lty=1)
    }
    par(op)
}
dev.off()


aic <- list()
best <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    X <- load_spp_pixel(spp)
    d <- data.frame(Indirect=as.numeric(X[,1]), as.matrix(A2s))
    ss <- sample(which(d$Native < 0.95), 5000)
    m1 <- gam(Indirect ~ s(Native, k=k), data=d)
    m2 <- gam(Indirect ~ s(Forestry, k=k), data=d)
    m3 <- gam(Indirect ~ s(Transportation, k=k), data=d)
    m4 <- gam(Indirect ~ s(EnSoft, k=k), data=d)
    m5 <- gam(Indirect ~ s(EnHard, k=k), data=d)
    m6 <- gam(Indirect ~ s(RuralUrban, k=k), data=d)
    m7 <- gam(Indirect ~ s(Agriculture, k=k), data=d)
    mm <- list(m1, m2, m3, m4, m5, m6, m7)
    names(mm) <- colnames(d)[-1]
    a <- sapply(mm, AIC)
    a <- a - min(a)
    aic[[spp]] <- a
    best[[spp]] <- names(mm)[which.min(a)]
}

aic <- do.call(rbind, aic)
best <- unlist(best)
data.frame(best=table(best))

#    best.best best.Freq
#1 Agriculture        22
#2      EnSoft         2
#3    Forestry         1
#4      Native        59
#5  RuralUrban         1
#> best[best=="Forestry"]
#      WIWR
#"Forestry"
#> best[best=="EnSoft"]
#    BHCO     BTNW
#"EnSoft" "EnSoft"
#> best[best=="RuralUrban"]
#        AMCR
#"RuralUrban"

