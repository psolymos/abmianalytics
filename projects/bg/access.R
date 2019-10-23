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

x <- read.csv("d:/abmi/AB_data_v2019/data/raw/species/bg/BGxy.csv")


library(sp)
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))

coordinates(x) <- ~ longitude + latitude
proj4string(x) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
x <- spTransform(x, proj4string(rt))

x@data$gid <- with(x@data, interaction(data_set, site, sep="_"))
x@data$stid <- with(x@data, interaction(data_set, site, station, sep="_"))
levels(x@data$gid)

dim(x)

op <- par(mfrow=c(3,3), mar=rep(1,4)+0.1)
for (i in 1:9)
    plot(x[x@data$gid == levels(x@data$gid)[i],], main=levels(x@data$gid)[i])
par(op)

op <- par(mfrow=c(3,3), mar=rep(1,4)+0.1)
for (i in 10:18)
    plot(x[x@data$gid == levels(x@data$gid)[i],], main=levels(x@data$gid)[i])
par(op)


xy <- read.csv("d:/abmi/AB_data_v2019/data/raw/species/bg/BGxy.csv")
str(xy)
head(xy)
table(xy$site)

