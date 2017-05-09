## processing XY for intersections etc

x <- read.csv("e:/peter/AB_data_v2017/data/raw/xy/ForPeterS-BU_May9-2017.csv")
x <- x[x$MaxOfYear <= 2017,]

x1 <- x[x$CountOfYear <= 2,]
x1$YEAR <- x1$MinOfYear
x2 <- x[x$CountOfYear == 2,]
x2$YEAR <- x2$MaxOfYear

xx <- droplevels(rbind(x1, x2))

xx$SS <- as.factor(with(xx, paste(ProjectID, Cluster, SITE, STATION, sep="::")))
xx <- droplevels(xx[,c("SS","YEAR","UTM_Zone","EASTING","NORTHING","Latitude","Longitude")])
write.csv(xx, row.names=FALSE,
    file="e:/peter/AB_data_v2017/data/raw/xy/xy-bu-bird-points-2017-05-09.csv")
