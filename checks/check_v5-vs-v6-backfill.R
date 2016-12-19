library(mefa4)
load("e:/peter/AB_data_v2016/out/kgrid/kgrid_table.Rdata")
fv5 <- "e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata"
fv6 <- "e:/peter/AB_data_v2016/data/kgrid-V6/veg-hf_1kmgrid_v6.Rdata"

e <- new.env()
load(fv5, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
v5 <- e$dd1km_pred[1:2]

e <- new.env()
load(fv6, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
v6 <- e$dd1km_pred[1:2]

rm(e)

##

compare_sets(colnames(v5[[1]]), colnames(v6[[1]]))
intersect(colnames(v5[[1]]), colnames(v6[[1]]))
setdiff(colnames(v5[[1]]), colnames(v6[[1]]))
setdiff(colnames(v6[[1]]), colnames(v5[[1]]))

#cn <- data.frame(v6=colnames(v6[[1]]))
#cn$v5 <- colnames(v5[[1]])[match(colnames(v6[[1]]), colnames(v5[[1]]))]
#write.csv(cn, row.names=FALSE,
#    file="e:/peter/AB_data_v2016/data/kgrid-V6/xwalk-v5-v6.csv")

compare_sets(colnames(v5[[2]]), colnames(v6[[2]]))
intersect(colnames(v5[[2]]), colnames(v6[[2]]))
setdiff(colnames(v5[[2]]), colnames(v6[[2]]))
setdiff(colnames(v6[[2]]), colnames(v5[[2]]))

xt <- read.csv("e:/peter/AB_data_v2016/data/kgrid-V6/xwalk-v5-v6.csv")
setdiff(colnames(v5[[1]]), xt$v5)
setdiff(colnames(v5[[2]]), xt$v5)
setdiff(colnames(v6[[1]]), xt$v6)
setdiff(colnames(v6[[2]]), xt$v6)

xt$v5 <- as.character(xt$v5)
xt$v6 <- as.character(xt$v6)
xt$v5[is.na(xt$v5)] <- "XXX"
xt$v6[is.na(xt$v6)] <- "XXX"

v5[[1]] <- cBind(v5[[1]], "XXX"=0)
v5[[2]] <- cBind(v5[[2]], "XXX"=0)
v6[[1]] <- cBind(v6[[1]], "XXX"=0)
v6[[2]] <- cBind(v6[[2]], "XXX"=0)

v5[[1]] <- v5[[1]] / ifelse(rowSums(v5[[1]])==0, 1, rowSums(v5[[1]]))
v5[[2]] <- v5[[2]] / ifelse(rowSums(v5[[2]])==0, 1, rowSums(v5[[2]]))
v6[[1]] <- v6[[1]] / ifelse(rowSums(v6[[1]])==0, 1, rowSums(v6[[1]]))
v6[[2]] <- v6[[2]] / ifelse(rowSums(v6[[2]])==0, 1, rowSums(v6[[2]]))

## NSR level comparison

cr5 <- groupSums(v5[[1]], 2, xt$v5[match(colnames(v5[[1]]), xt$v5)])
rf5 <- groupSums(v5[[2]], 2, xt$v5[match(colnames(v5[[2]]), xt$v5)])

cr6 <- groupSums(v6[[1]], 2, xt$v5[match(colnames(v6[[1]]), xt$v6)])
rf6 <- groupSums(v6[[2]], 2, xt$v5[match(colnames(v6[[2]]), xt$v6)])

compare_sets(colnames(cr5), colnames(cr6))
compare_sets(colnames(rf5), colnames(rf6))

intersect(colnames(cr5), colnames(cr6))
setdiff(colnames(cr5), colnames(cr6))
setdiff(colnames(cr6), colnames(cr5))

colSums(cr5[,setdiff(colnames(cr5), colnames(cr6))])
cr5 <- cr5[,colnames(cr6)]
colSums(rf5[,setdiff(colnames(rf5), colnames(rf6))])
rf5 <- rf5[,colnames(rf6)]

#cr5 <- cr5 / ifelse(rowSums(cr5)==0, 1, rowSums(cr5))
#cr6 <- cr6 / ifelse(rowSums(cr6)==0, 1, rowSums(cr6))
#rf5 <- rf5 / ifelse(rowSums(rf5)==0, 1, rowSums(rf5))
#rf6 <- rf6 / ifelse(rowSums(rf6)==0, 1, rowSums(rf6))

if (TRUE) {
z <- colnames(cr5)
z1 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
z <- colnames(rf5)
z2 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)

cr5 <- groupSums(cr5, 2, z1)
cr6 <- groupSums(cr6, 2, z1)
rf5 <- groupSums(rf5, 2, z2)
rf6 <- groupSums(rf6, 2, z2)
}

cr5 <- groupMeans(cr5, 1, kgrid$NSRNAME)
cr6 <- groupMeans(cr6, 1, kgrid$NSRNAME)
rf5 <- groupMeans(rf5, 1, kgrid$NSRNAME)
rf6 <- groupMeans(rf6, 1, kgrid$NSRNAME)

col <- "grey"
cols <- colorRampPalette(c("blue","red"))(5)
pdf("e:/peter/AB_data_v2016/data/kgrid-V6/nsr-level-comparison-noage.pdf",
    width=10, height=5.5, onefile=TRUE)
for (i in 1:ncol(cr6)) {
    j <- colnames(cr6)[i]
    lim <- c(0, max(0.1, cr5[,j], cr6[,j]))
    if (j %in% colnames(rf6))
        lim <- c(0, max(lim, rf5[,j], rf6[,j]))
    op <- par(mfrow=c(1,2))
    d <- ceiling(100*abs(cr5[,j] - cr6[,j]))
    d[d < 1] <- 1
    d[d >= 5] <- 5
    plot(cr5[,j], cr6[,j], main=j,
        xlim=lim, ylim=lim,
        xlab="current 1km mean v5", ylab="current 1km mean v6", col=cols[d], pch=19)
    text(cr5[,j], cr6[,j], ifelse(d > 1, rownames(cr5), ""),
        col=cols[d], pos=3, cex=0.4, xpd=NA)
        abline(0,1,col=col,lty=1)
        abline(-0.01,1,col=col,lty=2)
        abline(0.01,1,col=col,lty=2)
        abline(-0.05,1,col=col,lty=3)
        abline(0.05,1,col=col,lty=3)

    if (j %in% colnames(rf6)) {
        d <- ceiling(100*abs(rf5[,j] - rf6[,j]))
        d[d < 1] <- 1
        d[d >= 5] <- 5
        plot(rf5[,j], rf6[,j], main="",
            xlim=lim, ylim=lim,
            xlab="reference 1km mean v5", ylab="reference 1km mean v6", col=cols[d], pch=19)
        text(rf5[,j], rf6[,j], ifelse(d > 1, rownames(cr5), ""),
            col=cols[d], pos=3, cex=0.4, xpd=NA)
        abline(0,1,col=col,lty=1)
        abline(-0.01,1,col=col,lty=2)
        abline(0.01,1,col=col,lty=2)
        abline(-0.05,1,col=col,lty=3)
        abline(0.05,1,col=col,lty=3)
    } else {
        plot.new()
    }
    par(op)
}
dev.off()

## only v6
crS <- rowSums(v6[[1]][,setdiff(colnames(v6[[1]]), colnames(v5[[1]]))]) / ifelse(rowSums(v6[[1]])==0, 1, rowSums(v6[[1]]))

## maps

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

cr5 <- groupSums(v5[[1]], 2, xt$v5[match(colnames(v5[[1]]), xt$v5)])
rf5 <- groupSums(v5[[2]], 2, xt$v5[match(colnames(v5[[2]]), xt$v5)])

cr6 <- groupSums(v6[[1]], 2, xt$v5[match(colnames(v6[[1]]), xt$v6)])
rf6 <- groupSums(v6[[2]], 2, xt$v5[match(colnames(v6[[2]]), xt$v6)])

cr5 <- cr5[,colnames(cr6)]
rf5 <- rf5[,colnames(rf6)]

z <- colnames(cr5)
z1 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
z <- colnames(rf5)
z2 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)

cr5 <- groupSums(cr5, 2, z1)
cr6 <- groupSums(cr6, 2, z1)
rf5 <- groupSums(rf5, 2, z2)
rf6 <- groupSums(rf6, 2, z2)

for (i in 1:ncol(rf6)) {

j <- colnames(rf6)[i]
cat(j, "\n");flush.console()
png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/maps/v5v6_", j, ".png"),
    width=3*600, height=1000)
op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,3))

MAX <- max(rf5[,j], rf6[,j])

iii <- as.integer(pmin(101, round(100*rf5[,j]/MAX)+1))
plot(kgrid$X, kgrid$Y, col=C1[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
mtext(side=3,paste(j, "v5 reference"),col="grey30", cex=legcex)
points(city, pch=18, cex=cex*2)
text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
for (ii in 1:101) {
    jj <- ii * abs(diff(c(5450000, 5700000)))/100
    segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
}
text(240000, 5450000, "0%")
text(240000, 5700000, paste0(round(MAX*100),"%"))

iii <- as.integer(pmin(101, round(100*rf6[,j]/MAX)+1))
plot(kgrid$X, kgrid$Y, col=C1[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
mtext(side=3,paste(j, "v6 reference"),col="grey30", cex=legcex)
points(city, pch=18, cex=cex*2)
text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
for (ii in 1:101) {
    jj <- ii * abs(diff(c(5450000, 5700000)))/100
    segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
}
text(240000, 5450000, "0%")
text(240000, 5700000, paste0(round(MAX*100),"%"))

iii <- as.integer(pmin(101, round(50+50*(rf6[,j]-rf5[,j]))+1))
plot(kgrid$X, kgrid$Y, col=C2[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
mtext(side=3,paste(j, "v6-v5 reference"),col="grey30", cex=legcex)
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

## v5 only

cr5 <- v5[[1]]
rf5 <- v5[[2]]
z <- colnames(cr5)
z1 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
z <- colnames(rf5)
z2 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
cr5 <- groupSums(cr5, 2, z1)
rf5 <- groupSums(rf5, 2, z2)


for (i in 1:ncol(cr5)) {

j <- colnames(cr5)[i]
cat(j, "\n");flush.console()
png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/maps/v5_", j, ".png"),
    width=3*600, height=1000)
op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,3))

MAX <- if (j %in% colnames(rf5))
    max(cr5[,j], rf5[,j]) else max(cr5[,j])

iii <- as.integer(pmin(101, round(100*cr5[,j]/MAX)+1))
plot(kgrid$X, kgrid$Y, col=C1[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
mtext(side=3,paste(j, "v5 current"),col="grey30", cex=legcex)
points(city, pch=18, cex=cex*2)
text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
for (ii in 1:101) {
    jj <- ii * abs(diff(c(5450000, 5700000)))/100
    segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
}
text(240000, 5450000, "0%")
text(240000, 5700000, paste0(round(MAX*100),"%"))

if (j %in% colnames(rf5)) {
    iii <- as.integer(pmin(101, round(100*rf5[,j]/MAX)+1))
    plot(kgrid$X, kgrid$Y, col=C1[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste(j, "v5 reference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 5700000, paste0(round(MAX*100),"%"))

    iii <- as.integer(pmin(101, round(50+50*(cr5[,j]-rf5[,j]))+1))
    plot(kgrid$X, kgrid$Y, col=C2[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste(j, "current-reference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C2[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "+100%")
    text(240000, 0.5*(5450000 + 5700000), "0%")
    text(240000, 5700000, "-100%")

} else {
    plot.new()
    plot.new()
}
par(op)
dev.off()
}

## v6 only

cr6 <- v6[[1]]
rf6 <- v6[[2]]
z <- colnames(cr6)
z1 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
z <- colnames(rf6)
z2 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
cr6 <- groupSums(cr6, 2, z1)
rf6 <- groupSums(rf6, 2, z2)


for (i in 1:ncol(cr6)) {

j <- colnames(cr6)[i]
cat(j, "\n");flush.console()
png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/maps/v6col_", j, ".png"),
    width=3*600, height=1000)
op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,3))

MAX <- if (j %in% colnames(rf6))
    max(cr6[,j], rf6[,j]) else max(cr6[,j])

iii <- as.integer(pmin(101, round(100*cr6[,j]/MAX)+1))
plot(kgrid$X, kgrid$Y, col=C3[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
mtext(side=3,paste(j, "v6 current"),col="grey30", cex=legcex)
points(city, pch=18, cex=cex*2)
text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
for (ii in 1:101) {
    jj <- ii * abs(diff(c(5450000, 5700000)))/100
    segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
}
text(240000, 5450000, "0%")
text(240000, 5700000, paste0(round(MAX*100),"%"))

if (j %in% colnames(rf6)) {
    iii <- as.integer(pmin(101, round(100*rf6[,j]/MAX)+1))
    plot(kgrid$X, kgrid$Y, col=C3[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste(j, "v6 reference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C1[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "0%")
    text(240000, 5700000, paste0(round(MAX*100),"%"))

    iii <- as.integer(pmin(101, round(50+50*(cr6[,j]-rf6[,j]))+1))
    plot(kgrid$X, kgrid$Y, col=C2[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,paste(j, "current-reference"),col="grey30", cex=legcex)
    points(city, pch=18, cex=cex*2)
    text(city[,1], city[,2], rownames(city), cex=0.8, adj=-0.1, col="grey10")
    for (ii in 1:101) {
        jj <- ii * abs(diff(c(5450000, 5700000)))/100
        segments(190000, 5450000+jj, 220000, 5450000+jj, col=C2[ii], lwd=2, lend=2)
    }
    text(240000, 5450000, "+100%")
    text(240000, 0.5*(5450000 + 5700000), "0%")
    text(240000, 5700000, "-100%")

} else {
    plot.new()
    plot.new()
}
par(op)
dev.off()
}

## total area check

e <- new.env()
load(fv5, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
ta5 <- sapply(e$dd1km_pred[1:4], rowSums)

e <- new.env()
load(fv6, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
ta6 <- sapply(e$dd1km_pred[1:4], rowSums)

rm(e)

all(rownames(kgrid) == rownames(ta5))
all(rownames(kgrid) == rownames(ta6))

colnames(ta5) <- paste0(colnames(ta5), "_v5")
colnames(ta6) <- paste0(colnames(ta6), "_v6")

ta <- cbind(ta5, ta6)

for (i in 1:8) {

cat(i, "\n");flush.console()
png(paste0("e:/peter/AB_data_v2016/data/kgrid-V6/maps/Atot_", colnames(ta)[i], ".png"),
    width=600, height=1000)
op <- par(mar=c(0, 0, 4, 0) + 0.1, mfrow=c(1,1))

    iii <- as.integer(pmin(101, round(50+50*(ta[,i]/10^6-kgrid$Area_km2))+1))
    plot(kgrid$X, kgrid$Y, col=C2[iii], pch=15, cex=cex, ann=FALSE, axes=FALSE)
    with(kgrid[kgrid$pWater > 0.99,], points(X, Y, col=CW, pch=15, cex=cex))
    mtext(side=3,colnames(ta)[i],col="grey30", cex=legcex)
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

