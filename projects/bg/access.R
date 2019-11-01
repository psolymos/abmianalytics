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
x <- read.csv("d:/abmi/AB_data_v2019/data/raw/species/bg/BGxy.csv")
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


xy <- read.csv("d:/abmi/AB_data_v2019/data/raw/species/bg/BGxy.csv")
str(xy)
head(xy)
table(xy$site)


## -- summarizing big grids

#e:/peter/AB_data_v2018/data/raw/veghf/abmi

HF_VERSION <- "2016_fine"
source("~/repos/abmianalytics/veghf/veghf-setup.R")
load(file.path(ROOT, VER, "data", "analysis", "ages-by-nsr.Rdata"))

d <- read.csv("d:/abmi/AB_data_v2019/data/analysis/bg/20191029_grids_600sqm_veg61HF2016.csv")
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

d <- read.csv("d:/abmi/AB_data_v2019/data/analysis/bg/20191031_grids_600sqm_veg61HF2017.csv")
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

