library(mefa4)
library(data.table)
getOption("datatable.fread.datatable")
options(datatable.fread.datatable=FALSE)
getOption("datatable.fread.datatable")

## HF 2014

load("e:/peter/AB_data_v2016/out/kgrid/kgrid_table.Rdata")

fv5 <- "e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata"
e <- new.env()
load(fv5, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
v5 <- e$dd1km_pred[[1]]
rm(e)

x <- fread("e:/peter/AB_data_v2016/raw_new/hf2014/hf2014_union_1by1.csv")
str(x)

lt <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type-v2014.csv")
lt$Description <- NULL
lt$ref <- ifelse(lt$Refined=="", as.character(lt$HF_GROUP), as.character(lt$Refined))

yrft <- Xtab(SHAPE_Area ~ x$YEAR + x$FEATURE_TY, x)
round(100 * rowSums(yrft) / sum(yrft),2)

round(100*yrft[,"CULTIVATION_ABANDONED"][yrft[,"CULTIVATION_ABANDONED"] > 0] / sum(yrft[,"CULTIVATION_ABANDONED"]), 2)

rn <- rownames(yrft)
rn[-1] <- "1"
yrft2 <- groupSums(yrft, 1, rn)
yrft2 <- as.matrix(t(yrft2) / colSums(yrft))
yr0p <- 100 * yrft["0",] / colSums(yrft)

yrft2 <- round(100 * yrft2[-1,], 2)
colnames(yrft2) <- c("PercentNoAge", "PercentWithAge")
write.csv(yrft2, file="e:/peter/AB_data_v2016/raw_new/hf2014/feature-types-ages.csv")

## reclass
## compare to 2012
## map

table(x$YEAR)
compare_sets(x$FEATURE_TY, lt$FEATURE_TY)
setdiff(x$FEATURE_TY, lt$FEATURE_TY)

map12 <- cbind(
    c("", as.character(lt$FEATURE_TY)),
    c("", as.character(lt$HF_GROUP)))
map14 <- cbind(
    c("", as.character(lt$FEATURE_TY)),
    c("", as.character(lt$ref)))
x$FT12 <- reclass(x$FEATURE_TY, map12)
x$FT14 <- reclass(x$FEATURE_TY, map14)

y <- Xtab(SHAPE_Area ~ Row_Col + FEATURE_TY, x)
y <- y[rownames(kgrid),]

compare_sets(levels(x$FEATURE_TY), lt$FEATURE_TY)
setdiff(levels(x$FEATURE_TY), lt$FEATURE_TY)

y12 <- Xtab(SHAPE_Area ~ Row_Col + FT12, x)
y14 <- Xtab(SHAPE_Area ~ Row_Col + FT14, x)
y12 <- y12[rownames(kgrid),]
y14 <- y14[rownames(kgrid),]

cn <- colnames(v5)
cn[substr(cn, 1, 2) == "CC"] <- "CutBlocks"
cn[cn %notin% colnames(y12)] <- ""
vv5 <- groupSums(v5, 2, cn)
compare_sets(colnames(vv5), colnames(y12))
vv5 <- vv5[,colnames(y12)]

## subregional comparisons

vv5ns <- groupSums(vv5, 1, kgrid$NSRNAME)
y12ns <- groupSums(y12, 1, kgrid$NSRNAME)

col <- "grey"
cols <- colorRampPalette(c("blue","red"))(5)

pdf("e:/peter/AB_data_v2016/data/kgrid-V6/nsr-level-comparison-HF.pdf",
    width=5, height=5.5, onefile=TRUE)
for (i in 1:ncol(vv5)) {
    j <- colnames(vv5)[i]
    cr5 <- vv5ns[,j] / rowSums(vv5ns) # 2012 HF
    cr6 <- y12ns[,j] / rowSums(y12ns) # 2014 HF reclasses for 12
    if (j == "") {
        j <- "Total HF"
        cr5 <- 1 - cr5
        cr6 <- 1 - cr6
    }
    lim <- c(0, max(0.1, cr5, cr6))
    d <- ceiling(100*abs(cr5 - cr6))
    d0 <- d
    d[d < 1] <- 1
    d[d >= 5] <- 5
    plot(cr5, cr6, main=j,
        xlim=lim, ylim=lim,
        xlab="HF 2012", ylab="HF 2014", col=cols[d], pch=19)
    text(cr5, cr6, ifelse(d0 > 0.5, names(cr5), ""),
        col=cols[d], pos=3, cex=0.4, xpd=NA)
    abline(0,1,col=col,lty=1)
    abline(-0.01,1,col=col,lty=2)
    abline(0.01,1,col=col,lty=2)
    abline(-0.05,1,col=col,lty=3)
    abline(0.05,1,col=col,lty=3)
}
dev.off()

## pixel level distribution

cex <- 0.25
legcex <- 1.5
C1 <- colorRampPalette(c("#edf8e9","#bae4b3","#74c476","#31a354","#006d2c"))(101)
C2 <- colorRampPalette(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6"))(101)
C3 <- C1[ceiling(100*sqrt(seq(0,1,len=101)))+1] # sqrt transformed
CW <- rgb(0.4,0.3,0.8) # water
CE <- "lightcyan4" # exclude
library(raster)
library(sp)
library(rgdal)
city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
city <- spTransform(city, CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
city <- as.data.frame(city)

cr5a <- vv5 / rowSums(vv5)
cr6a <- y12 / rowSums(y12)

y14re <- y14 / rowSums(y14)
y14re <- y14re[,colnames(y14re) %notin% colnames(y12)]

for (i in 1:ncol(cr5a)) {

    j <- colnames(cr5a)[i]
    cr5 <- cr5a[,j]
    cr6 <- cr6a[,j]
    if (j == "") {
        j <- "TotalHF"
        cr5 <- 1 - cr5
        cr6 <- 1 - cr6
    }
    cat(j, "\n");flush.console()

    png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/hfmaps/HF12-HF14_", j, ".png"),
        width=3*600, height=1000)
    op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,3))

    MAX <- max(cr5, cr6)

    iii <- as.integer(pmin(101, round(100*cr5/MAX)+1))
    plot(kgrid$X, kgrid$Y, col=C3[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste("HF2012:", j),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 5700000, paste0(round(MAX*100),"%"))

    iii <- as.integer(pmin(101, round(100*cr6/MAX)+1))
    plot(kgrid$X, kgrid$Y, col=C3[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste("HF2014:", j),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 5700000, paste0(round(MAX*100),"%"))

    iii <- as.integer(pmin(101, round(50+50*(cr6-cr5))+1))
    plot(kgrid$X, kgrid$Y, col=C2[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste("HF2014-HF2012:", j),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C2[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "+100%")
    text(240000, 0.5*(5450000 + 5700000), "0%")
    text(240000, 5700000, "-100%")

    par(op)
    dev.off()
}

## map the refined classes only

for (i in 1:ncol(y14re)) {

    j <- colnames(y14re)[i]
    cr5 <- y14re[,j]
    if (j == "") {
        j <- "TotalHF"
        cr5 <- 1 - cr5
    }
    cat(j, "\n");flush.console()

    png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/hfmaps/HF14refined_", j, ".png"),
        width=1*600, height=1000)
    op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,1))

    MAX <- max(cr5)

    iii <- as.integer(pmin(101, round(100*cr5/MAX)+1))
    plot(kgrid$X, kgrid$Y, col=C3[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste("HF2012:", j),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 5700000, paste0(round(MAX*100),"%"))
    par(op)
    dev.off()
}

## all feature types

rs <- rowSums(y)
for (i in 1:ncol(y)) {

    j <- colnames(y)[i]
    cr5 <- y[,j] / rs
    if (j != "") {
    cat(j, "\n");flush.console()

    png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/hfmaps/HF14ft_", j, ".png"),
        width=1*600, height=1000)
    op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,1))

    MAX <- max(cr5)

    iii <- as.integer(pmin(101, round(100*cr5/MAX)+1))
    plot(kgrid$X, kgrid$Y, col=C3[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste("HF2012:", j),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 5700000, paste0(round(MAX*100),"%"))
    par(op)
    dev.off()
    }
}
