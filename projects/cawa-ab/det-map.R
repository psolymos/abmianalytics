library(mefa4)
library(RColorBrewer)

ROOT <- "e:/peter/AB_data_v2016/out/birds"

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
#up()
source("~/repos/abmianalytics/R/wrsi_functions.R")

load(file.path(ROOT, "data", "data-wrsi-all.Rdata"))

TAX$Fn <- droplevels(TAX$English_Name)
levels(TAX$Fn) <- nameAlnum(levels(TAX$Fn), capitalize="mixed", collapse="")
DATw$NSRxLUF <- interaction(DATw$NSRNAME, DATw$LUF_NAME, sep="_", drop=TRUE)

if (FALSE) {
xn <- DATw[DATw$useNorth,]
yn <- YYw[DATw$useNorth,]
yn <- as.matrix(yn)
xs <- DATw[DATw$useSouth,]
ys <- YYw[DATw$useSouth,]
ys <- as.matrix(ys)

pveghf <- pVeg[rownames(yn), c("Decid", "Mixwood", "Conif", "Pine",
    "BSpr", "Larch", "Open", "Shrub", "Wet",
    "HFor", "Cult", "UrbInd", "HardLin", "SoftLin")]
pveghf[,"Open"] <- pveghf[,"Open"] + pveghf[,"Shrub"]
pveghf <- pveghf[,colnames(pveghf) != "Shrub"]
colnames(pveghf) <- c("Deciduous", "Mixedwood", "White Spruce", "Pine",
    "Black Spruce", "Larch", "Open", "Wet",
    "Forestry", "Cultivated", "Urban/Industrial", "Hard Linear", "Soft Linear")

psoilhf <- pSoil[rownames(ys),]
psoilhf <- psoilhf[,c("Productive", "Clay", "Saline", "RapidDrain", "Cult", "UrbInd")]
colnames(psoilhf) <- c("Productive", "Clay", "Saline", "Rapid Drain",
    "Cultivated", "Urban/Industrial")

tab_srv_n <- as.matrix(groupSums(YYw[DATw$useNorth,], 1, DATw[DATw$useNorth,"NSRxLUF"]))
tab_srv_s <- as.matrix(groupSums(YYw[DATw$useSouth,], 1, DATw[DATw$useSouth,"NSRxLUF"]))
tab_srv_a <- as.matrix(groupSums(YYw, 1, DATw[,"NSRxLUF"]))
colnames(tab_srv_n) <- colnames(tab_srv_s) <- colnames(tab_srv_a) <- TAX$English_Name

tab_loc_n <- as.matrix(groupSums(ifelse(yn > 0, 1, 0), 1, DATw[DATw$useNorth,"NSRxLUF"]))
tab_loc_s <- as.matrix(groupSums(ifelse(ys > 0, 1, 0), 1, DATw[DATw$useSouth,"NSRxLUF"]))
tab_loc_a <- as.matrix(groupSums(ifelse(as.matrix(YYw) > 0, 1, 0), 1, DATw[,"NSRxLUF"]))
colnames(tab_loc_n) <- colnames(tab_loc_s) <- colnames(tab_loc_a) <- TAX$English_Name

write.csv(tab_srv_n, file=file.path(ROOT, "tables", "birds-number-of-surveys-north.csv"))
write.csv(tab_srv_s, file=file.path(ROOT, "tables", "birds-number-of-surveys-south.csv"))
write.csv(tab_srv_a, file=file.path(ROOT, "tables", "birds-number-of-surveys-all.csv"))
write.csv(tab_loc_n, file=file.path(ROOT, "tables", "birds-number-of-locations-north.csv"))
write.csv(tab_loc_s, file=file.path(ROOT, "tables", "birds-number-of-locations-south.csv"))
write.csv(tab_loc_a, file=file.path(ROOT, "tables", "birds-number-of-locations-all.csv"))

yn <- yn[,colSums(yn) > 0]
ys <- ys[,colSums(ys) > 0]
}

## spp specific output

spp <- "CAWA"

## map-det

load(file.path("e:/peter/AB_data_v2016/out", "kgrid", "kgrid_table.Rdata"))
col1 <- c("#C8FBC8","#C8E6FA","#F5E6F5","#FFDCEC","#FFE6CD","#FFF1D2")[match(kgrid$NRNAME,
    c("Boreal","Foothills","Rocky Mountain","Canadian Shield","Parkland","Grassland"))]

library(raster)
library(sp)
library(rgdal)
city <-data.frame(x = -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60,
    y = c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60)
rownames(city) <- c("Calgary","Edmonton","Lethbridge","Fort McMurray",
    "High Level","Grande Prairie")
coordinates(city) <- ~ x + y
proj4string(city) <- CRS(paste0("+proj=longlat +datum=WGS84 ",
    "+ellps=WGS84 +towgs84=0,0,0"))
city <- as.data.frame(spTransform(city, CRS(paste0("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 ",
    "+x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))))
xyw <- as.matrix(kgrid[kgrid$pWater >= 0.99,c("X","Y")])
blank <- matrix(0, 0, 2)



    ## mask (e.g. for CAWA project)
    #iii <- rep(TRUE, nrow(DATw))
    iii <- DATw$POINT_Y > 50 & DATw$NRNAME != "Grassland"

#    cat("map-det", spp, i, "/", ncol(YYw), "\n");flush.console()
    xy0 <- as.matrix(DATw[YYw[,spp] == 0 & iii, c("X","Y")])
    xy1 <- as.matrix(DATw[YYw[,spp] > 0 & iii, c("X","Y")])

    NAM <- as.character(TAX[spp, "English_Name"])
    (nsrv <- sum(YYw[iii,spp]))
    (nloc <- sum(YYw[iii,spp] > 0))
    #fname <- file.path(ROOT, "figs", "map-det",
    #    paste0(as.character(TAX[spp, "Fn"]), ".png"))
    fname <- "c:/Users/Peter/Dropbox/josm/cawa-jeff/revision/cawa-map-det.png"
	png(file=fname, width=600, height=1000)
#    postscript(paste0("e:/peter/", as.character(tax[spp, "file"]), "-detections.eps"),
#        horizontal = FALSE, onefile = FALSE, paper = "special",
#        width=6, height=10)
    plot(kgrid$X, kgrid$Y, pch=15, cex=0.2, col=col1, axes=FALSE, ann=FALSE)
    iiii <- !(kgrid$POINT_Y > 50 & kgrid$NRNAME != "Grassland")
    points(kgrid$X[iiii], kgrid$Y[iiii], pch=15, cex=0.2, col="grey")
    points(xyw, pch=15, cex=0.2, col=rgb(0.3,0.45,0.9))
    #points(xy0, pch=19, cex=0.5, col="red3")
    points(xy0, pch=3, cex=0.6, col="red2")
    points(xy1, pch=19, cex=1.6, col="red4")
    #mtext(paste0(NAM, "\n(detected at ", nloc, " locations, ", nsrv, " surveys)"), line=1,
    #    side=3, adj=0.5, cex=1.4, col="grey40")
    points(city, pch=18, col="grey10")
    text(city, rownames(city), cex=1, adj=-0.1, col="grey10")
    legend("bottomleft", pch=c(15,15,15,15,15,NA,15, NA, 3,19),
        col=c("#C8FBC8","#C8E6FA","#F5E6F5","#FFDCEC","#FFE6CD",NA, "grey", NA, "red2", "red4"),
        legend=c(
        "Boreal","Foothills","Rocky Mountain","Canadian Shield","Parkland",
        NA,"Outside of study area",
        NA, "Survey location","Detection"), title="Natural Regions",
        cex=1.2, pt.cex=c(4,4,4,4,4,NA,4, NA, 1.5,2), bty="n")
	dev.off()

