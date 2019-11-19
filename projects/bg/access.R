## big grid test

library(DBI)
library(RPostgres)
con <- dbConnect(RPostgres::Postgres(),
    dbname = 'bis',
    host = 'prod.wildtrax.ca',
    port = 5432,
    user = 'science_center',
    password = 'b1rd2Science')
(tl <- dbListTables(con))

q <- "SELECT DISTINCT
data_set,site,station,latitude,longitude,
SUBSTRING(recording_date::VARCHAR FROM 1 FOR 4) AS year_deploy,
is_buffered_location
FROM recording
where data_set like 'BG%'"

x <- dbGetQuery(con, q)

dbDisconnect(con)



library(sp)
library(raster)
x <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/BGxy.csv")
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

coordinates(x) <- ~ longitude + latitude
proj4string(x) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
x <- spTransform(x, proj4string(rt))

x@data$gid <- with(x@data, interaction(data_set, site, sep="_"))
x@data$stid <- with(x@data, interaction(data_set, site, station, sep="_"))
levels(x@data$gid)

dim(x)

out <- NULL
for (i in c(2,4:18)) {
    z <- x[x@data$gid == levels(x@data$gid)[i],]
    if (i==2 || i==12) {
        s <- c(100, 1)
    } else {
        s <- c(91, 10)
    }
    if (i == 14)
        s[2] <- 19
    xy1 <- coordinates(z)[z@data$station == s[1],]
    xy2 <- coordinates(z)[z@data$station == s[2],]
    d <- sqrt(sum((xy1-xy2)^2))
    de <- sqrt((9*600)^2 + (9*600)^2)
    v <- data.frame(
        gid=paste0("BG_", c(i, i)),
        stid=paste0("BG_", c(i, i), "_", s),
        x=c(xy1[1], xy2[1]),
        y=c(xy1[2], xy2[2]))
    out <- rbind(out, v)
}
write.csv(out, row.names=FALSE, file="BG-corner-points-2019-10-29.csv")

op <- par(mfrow=c(3,3), mar=rep(1,4)+0.1)
for (i in c(2,4:9)) {
    z <- x[x@data$gid == levels(x@data$gid)[i],]
    if (i==2) {
        xy1 <- coordinates(z)[z@data$station == 100,]
        xy2 <- coordinates(z)[z@data$station == 1,]
    } else {
        xy1 <- coordinates(z)[z@data$station == 91,]
        xy2 <- coordinates(z)[z@data$station == 10,]
    }
    d <- sqrt(sum((xy1-xy2)^2))
    de <- sqrt((9*600)^2 + (9*600)^2)
    m <- paste(levels(x@data$gid)[i], round(d-de))
    plot(z, main=m)
    if (i==2) {
        points(z[z@data$station == 100,], pch=19, col=2)
        points(z[z@data$station == 1,], pch=19, col=4)
    } else {
        points(z[z@data$station == 91,], pch=19, col=2)
        points(z[z@data$station == 10,], pch=19, col=4)
    }
}
par(op)

op <- par(mfrow=c(3,3), mar=rep(1,4)+0.1)
for (i in c(10:18)) {
    z <- x[x@data$gid == levels(x@data$gid)[i],]
    if (i==12) {
        xy1 <- coordinates(z)[z@data$station == 100,]
        xy2 <- coordinates(z)[z@data$station == 1,]
    } else {
        xy1 <- coordinates(z)[z@data$station == 91,]
        xy2 <- coordinates(z)[z@data$station == 10,]
    }
    d <- sqrt(sum((xy1-xy2)^2))
    de <- sqrt((9*600)^2 + (9*600)^2)
    m <- paste(levels(x@data$gid)[i], round(d-de))
    plot(z, main=m)
    if (i==12) {
        points(z[z@data$station == 100,], pch=19, col=2)
        points(z[z@data$station == 1,], pch=19, col=4)
    } else {
        points(z[z@data$station == 91,], pch=19, col=2)
        if (i == 14) {
            points(z[z@data$station == 19,], pch=19, col=4)
        } else {
            points(z[z@data$station == 10,], pch=19, col=4)
        }
    }
}
par(op)

f <- function(i) {
    z <- x[x@data$gid == levels(x@data$gid)[i],]
    h <- coordinates(z)
    plot(h, type="n", axes=FALSE, ann=FALSE, main=i)
    text(h[,1], h[,2], z@data$station, cex=0.6)
}


xy <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/BGxy.csv")
str(xy)
head(xy)
table(xy$site)


## -- summarizing big grids

#e:/peter/AB_data_v2018/data/raw/veghf/abmi

HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))

## 600x600
d <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/20191029_grids_600sqm_veg61HF2016.csv")
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "600m x 600m BC cells [V61 backfilled + 2016 HFI]"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd16 <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

d <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/20191031_grids_600sqm_veg61HF2017.csv")
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "600m x 600m BC cells [V61 backfilled + 2016 HFI]"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dd17 <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

all(rownames(dd17[[1]]) == rownames(dd16[[1]]))
dif <- dd17[[1]] - dd16[[1]]
sum(rowSums(abs(dif)) > 0)

## 150m buffers
d <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/20191118_grids_points_buf150m_veg61hf2016.csv")
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "150m radius circle buffer in BG [V61 backfilled + 2016 HFI]"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dx$NSRNAME[dx$NSRNAME == ""] <- "Central Mixedwood"
ddp16 <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

d <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/20191118_grids_points_buf150m_veg61hf2017.csv")
dd <- make_vegHF_wide_v6(d,
    col.label="UID",
    col.year=2016,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE) # use refined classes
dd$scale <- "150m radius circle buffer in BG [V61 backfilled + 2016 HFI]"
dx <- nonDuplicated(d, UID, TRUE)[rownames(dd[[1]]),]
dx$NSRNAME[dx$NSRNAME == ""] <- "Central Mixedwood"
ddp17 <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

all(rownames(dd16[[1]]) == rownames(ddp16[[1]]))
all(rownames(ddp17[[1]]) == rownames(ddp16[[1]]))

save(dd16, dd17, ddp16, ddp17, file="d:/abmi/AB_data_v2019/data/misc/bg/veghf-summaries.RData")


## this is the GIS based index
dx1 <- data.frame(id_gis=dx$UID, x_gis=dx$POINT_X, y_gis=dx$POINT_Y)
dx1$site_gis <- as.integer(sapply(strsplit(as.character(dx1$id_gis), "_"), "[[", 2))
dx1$station_gis <- as.integer(sapply(strsplit(as.character(dx1$id_gis), "_"), "[[", 3))
dx1 <- dx1[order(dx1$site_gis, dx1$station_gis),]
m1 <- matrix(1:100, 10, 10)
m2 <- t(m1[,10:1])
t(m2[10:1,])
as.integer(t(m1[,10:1]))
dx1$station_ideal <- dx1$station_gis
for (i in unique(dx1$site_gis)) {
    if (!(i %in% c(2, 12))) {
        dx1$station_ideal[dx1$site_gis==i] <- as.integer(t(m1[,10:1]))
    }
}
f <- function(i, x) {
    x <- x[x[,4] == i,]
    h <- x[,2:3]
    plot(h, type="l", axes=FALSE, main=i, col="grey")
    text(h[,1], h[,2], x[,5], cex=0.6)
}
sort(unique(dx1$site_gis))
f(2, dx1)# bottomleft
f(4, dx1)
f(5, dx1)
f(6, dx1)
f(7, dx1)
f(9, dx1)
f(10, dx1)
f(11, dx1)
f(12, dx1)# bottomleft
f(13, dx1)
f(14, dx1)
f(15, dx1)
f(16, dx1)
f(17, dx1)
f(18, dx1)

## this is the original sample id
dx2 <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/BGxy.csv")
dx2$gid <- with(dx2, interaction(data_set, site, sep="_"))
dx2$stid <- with(dx2, interaction(data_set, site, station, sep="_"))
rownames(dx1) <- dx1[,1]

dx2 <- dx2[,c("stid", "longitude", "latitude", "site", "station")]
colnames(dx2) <- c("id_brd", "x_brd", "y_brd", "site_brd", "station_brd")
dx2 <- dx2[!duplicated(dx2[,1]),]
rownames(dx2) <- dx2[,1]

library(sp)
library(raster)
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

xy <- dx2
coordinates(xy) <- ~ x_brd + y_brd
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, proj4string(rt))
head(coordinates(xy))
dx2$x_brd <- coordinates(xy)[,1]
dx2$y_brd <- coordinates(xy)[,2]

dx2 <- dx2[order(dx2$site_brd, dx2$station_brd),]
dx2 <- droplevels(dx2[dx2$site_brd %in% unique(dx1$site_gis),])

sort(unique(dx2$site_brd))
f(2, dx2)
f(4, dx2)
f(5, dx2)
f(6, dx2)
f(7, dx2)
f(9, dx2)
f(10, dx2)
f(11, dx2)
f(12, dx2)
f(13, dx2)
f(14, dx2) # weird - we assume points are where they should be
f(15, dx2)
f(16, dx2)
f(17, dx2)
f(18, dx2)

compare_sets(dx1$id_gis, dx2$id_brd)

rownames(dx1) <- paste0("BG_", dx1$site_gis, "_", dx1$station_ideal)
dx0 <- data.frame(dx1, dx2[match(rownames(dx1), dx2$id_brd),])



pc <- read.csv("d:/abmi/AB_data_v2019/data/misc/bg/allgrids-bird-counts.csv")
pcx <- pc[,1:26]
pcy <- pc[,27:271]
pcx$site_pc <- as.integer(sapply(strsplit(as.character(pcx$Site), "-"), "[[", 2))
pcx$station_pc <- as.integer(sapply(strsplit(as.character(pcx$StationKey), "-"), "[[", 3))
table(pcx$site_pc)
table(pcx$station_pc)
pcx$id_pc <- paste0("BG_", pcx$site_pc, "_", pcx$station_pc)
dx3 <- nonDuplicated(pcx, id_pc, TRUE)
dx3 <- droplevels(dx3[dx3$site_pc %in% unique(dx1$site_gis),])


dx3 <- dx3[,c("id_pc", "Longitude", "Latitude", "site_pc", "station_pc", "Year")]
colnames(dx3) <- c("id_pc", "x_pc", "y_pc", "site_pc", "station_pc", "year")
rownames(dx3) <- dx3[,1]

sort(unique(dx3$site_pc))
f(2, dx3)
f(4, dx3)
f(5, dx3)
f(6, dx3)
f(7, dx3)
f(9, dx3)
f(10, dx3)
#f(11, dx3)
#f(12, dx3)
f(13, dx3)
f(14, dx3)

f(15, dx3)
f(16, dx3)
f(17, dx3)
f(18, dx3)

compare_sets(dx1$id_gis, dx3$id_pc)

dx0$had_brd_data <- dx0$id_gis %in% dx3$id_pc
dx0$year_sampled <- dx3$year[match(dx0$id_brd, dx3$id_pc)]

dd16cr <- as.matrix(dd16[[1]])
dd16cr <- dd16cr[match(dx0$id_gis, rownames(dd16cr)),]
dd17cr <- as.matrix(dd17[[1]])
dd17cr <- dd17cr[match(dx0$id_gis, rownames(dd17cr)),]
ddcr <- dd17cr
ii <- !is.na(dx0$year_sampled) & dx0$year_sampled < 2017
ddcr[ii,] <- dd16cr[ii,]

rownames(ddcr) <- rownames(dx0)


#write.csv(ddcr, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-from-peter-VEGHF-2019-11-14.csv")
#write.csv(dx0, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-from-peter-IDS-2019-11-14.csv")

grids <- data.frame(station=1:100)

m2 <- matrix(0, 10, 10)
z <- 0
for (i in 1:10) {
    if (i %% 2) {
        m2[i,] <- rep((z+1):(z+10/2), each=2)
    } else {
        m2[i,] <- m2[i-1,]
        z <- z + 10/2
    }
}
grids$g2x2 <- as.integer(m2)

m3 <- matrix(0, 10, 10)
z <- 0
for (i in 1:9) {
    if (i %% 3) {
        m3[i,] <- c(rep((z+1):(z+3), each=3), 0)
    } else {
        m3[i,] <- m3[i-1,]
        z <- z + 3
    }
}
grids$g3x3 <- as.integer(m3)

m4 <- matrix(0, 10, 10)
z <- 0
for (i in 1:8) {
    if (i %% 4) {
        m4[i,] <- c(rep((z+1):(z+2), each=4), 0, 0)
    } else {
        m4[i,] <- m4[i-1,]
        z <- z + 2
    }
}
grids$g4x4 <- as.integer(m4)

m5 <- matrix(0, 10, 10)
z <- 0
for (i in 1:10) {
    if (i %% 5) {
        m5[i,] <- rep((z+1):(z+2), each=5)
    } else {
        m5[i,] <- m5[i-1,]
        z <- z + 2
    }
}
grids$g5x5 <- as.integer(m5)

write.csv(grids, file="d:/abmi/AB_data_v2019/data/misc/bg/bg-from-peter-GRIDS-2019-11-14.csv")
