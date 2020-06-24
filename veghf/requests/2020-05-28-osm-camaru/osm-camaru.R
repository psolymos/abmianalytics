library(mefa4)

fl <- c("veg-hf_CAMARU-2019_Veg61-vHF.Rdata",
    "veg-hf_CameraARU_v6verif_2017-2018-sites.Rdata",
    "veg-hf_CameraARU_v6verif_Additional.Rdata",
    "veg-hf_CameraARU_v6verif_Additional-20190228.Rdata",
    "veg-hf_CameraARU_v6verified.Rdata")

root <- "d:/abmi/AB_data_v2020/misc/osm-camaru/"

e <- new.env()
load(paste0(root, fl[1]), envir=e)
d1 <- e$d_wide$veg_current

e <- new.env()
load(paste0(root, fl[5]), envir=e)
names(e)
d2 <- e$dd_150m$veg_current

e <- new.env()
load("s:/AB_data_v2019/data/misc/bg/veghf-summaries.RData", envir=e)
d3 <- e$ddp16$veg_current

d <- rbind(d1, d2[,colnames(d1)])
rn <- rownames(d)
tmp1 <- strsplit(rn, "_")
YEAR <- as.integer(sapply(tmp1, function(z) z[length(z)]))
table(YEAR)
tmp2 <- sapply(strsplit(rn, "-"), function(z) {
    a <- strsplit(z, "_")[[1]][1]
    if (a == "OG")
        z[3] else a
})

z <- data.frame(Site=as.integer(tmp2), Year=YEAR)
rownames(z) <- rn
d <- d[!is.na(z$Site),]
z <- z[!is.na(z$Site),]


s <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
z$X <- s$PUBLIC_LONGITUDE[match(z$Site, s$SITE_ID)]
z$Y <- s$PUBLIC_LATTITUDE[match(z$Site, s$SITE_ID)]

library(sp)
library(rgdal)
ply <- readOGR(dsn=system.file("extdata/OSA_bound.geojson", package="cure4insect"))

xy <- z
coordinates(xy) <- ~ X + Y
proj4string(xy) <- proj4string(ply)

o <- over(xy, as(ply, "SpatialPolygons"))
i <- !is.na(o)

z$insideOSA <- i

save(z, d, file= "d:/abmi/AB_data_v2020/misc/osm-camaru/osm-camaru-results-2020-06-26.RData")

