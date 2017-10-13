setwd("~/repos/abmianalytics/apps/hfchange")
source("globals.R")

## preprocess data

load(file.path("e:/peter/AB_data_v2017", "data", "analysis",
    "veg-hf_3x7_1999-2015-summaries.Rdata"))

if (FALSE) {
crs <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

od <- setwd("~/Dropbox/courses/st-johns-2017/data/NatRegAB")
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
AB <- spTransform(AB, crs)
ABnr <- gUnaryUnion(AB, AB@data$NRNAME) # natural regions
setwd(od)
ABnr <- gSimplify(ABnr, tol=500, topologyPreserve=TRUE)
pts <- gis
coordinates(pts) <- ~ PUBLIC_LONGITUDE + PUBLIC_LATTITUDE
proj4string(pts) <-
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
pts <- spTransform(pts, crs)
ab <- spTransform(ABnr, '+proj=longlat +datum=WGS84')
save(ab, ABnr, pts, file=file.path(DIR, "AB_NatReg_Simplified.Rdata"))
}
load(file.path(DIR, "AB_NatReg_Simplified.Rdata"))


f <- function(z) {
    z <- z[,colnames(z) != "NATIVE"]
#    z <- cbind(Total=rowSums(z), z)
    z <- z[,c("Agriculture", "Forestry", "Energy", "RuralUrban", "Transportation", "Misc")]
    colnames(z) <- c("Agriculture", "Forestry", "Energy",
        "Urban", "Transportation", "Other")
#    AB <- colMeans(z)
    NR <- groupMeans(z, 1, gis$NATURAL_REGIONS)
    NR <- NR[c("Grassland", "Parkland", "Foothills", "Boreal", "Canadian Shield", "Rocky Mountain"),]
#    rbind(Alberta=AB, NR)
    NR
}

hf <- lapply(veg3x7_sector, f)
HF <- array(NA, c(dim(hf[[1]]), length(veg3x7_sector)))
dimnames(HF) <- list(rownames(hf[[1]]), colnames(hf[[1]]), names(veg3x7_sector))
for (i in 1:length(hf))
    HF[,,i] <- hf[[i]]

## plot

r <- 2:3
c <- 2:4
br <- FALSE
d <- get_data0(r,c,br)

plot(gvisLineChart(data.frame(Year=d$x, d$y), options=list(gvis.editor="Edit me!")))

get_plot(r, c, br)
get_hplot(r, c, br)
get_rplot(r, c, br)
get_map(r)

save(ABnr, pts, HF, file="hfchange.rda")


# make the coordinates a numeric matrix
qk_mx <- data.matrix(gis[,c("PUBLIC_LONGITUDE", "PUBLIC_LATTITUDE")])
# convert the coordinates to a multipoint feature
qk_mp <- st_multipoint(qk_mx)
# convert the multipoint feature to sf
qk_sf <- st_sf(st_cast(st_sfc(qk_mp), "POINT"), gis, crs=4326)

