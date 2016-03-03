library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"

int <- read.csv(file.path(ROOT, "data", "aru-coni", 
    "2015_CONI_DetectionsByRecording_1minInterval.csv"))
tms <- read.csv(file.path(ROOT, "data", "aru-coni",
    "2015_CONIPeent3.4_30_70_truepositives_details.csv"))
tms$file.path <- NULL # same levels
sit <- read.csv(file.path(ROOT, "data", "aru-coni",
    "CONImodel_ARU_sitesv2.csv"))
loc <- read.csv(file.path(ROOT, "data", "aru-coni",
    "CONImodel_ARU_locations-proper-ID.csv"))

## dd150m, dd1km, points
load(file.path(ROOT, "out", "aru", 
    "veg-hf-clim-reg_aru-coni_fix-fire_fix-age0.Rdata"))
e <- new.env()
load(file=file.path(ROOT, "out/abmi_onoff", 
        "veg-hf-clim-reg_abmi-onoff_fix-fire_fix-age0_with2015.Rdata"), envir=e)
dd150mPT_2015 <- e$dd150mPT_2015
dd1kmPT_2015 <- e$dd1kmPT_2015
climPT_2015 <- e$climPT_2015
rm(e)

points$oldID <- rownames(dd150m[[1]])
climPT_2015$oldID <- rownames(dd150mPT_2015[[1]])

compare_sets(points$POINT_ID, loc$POINT_ID)
points$ID <- loc$ID[match(points$POINT_ID, loc$POINT_ID)]
rownames(points) <- points$ID

## setdiff must be due to faulty units
compare_sets(sit$ID, int$ID)
compare_sets(sit$ID, tms$ID)
compare_sets(sit$ID, points$ID)
setdiff(sit$ID, points$ID) # ABMI + some stuff??
setdiff(points$ID, sit$ID) # faulty units

climPT_2015$SITE000 <- sapply(climPT_2015$Nearest, function(z) {
    ch <- as.character(z)
    paste0(paste(rep(0, 4-nchar(ch)), collapse=""), ch)
})
climPT_2015$ID <- with(climPT_2015, paste("ABMI", SITE000, Cam_ARU_Location, sep="-"))
compare_sets(sit$ID, climPT_2015$ID)
i <- intersect(sit$ID, climPT_2015$ID)
climPT <- climPT_2015[climPT_2015$ID %in% i,]
climPT <- nonDuplicated(climPT, ID, TRUE)

compare_sets(colnames(points), colnames(climPT))
cn <- intersect(colnames(points), colnames(climPT_2015))
pt <- rbind(data.frame(points[,cn], part=1), data.frame(climPT[,cn], part=2))
compare_sets(sit$ID, pt$ID)

veghfpc <- rBind(dd150m[[1]], dd150mPT_2015[[1]])
vegpc <- rBind(dd150m[[2]], dd150mPT_2015[[2]])
veghfkm <- rBind(dd1km[[1]], dd1kmPT_2015[[1]])
vegkm <- rBind(dd1km[[2]], dd1kmPT_2015[[2]])

veghfpc <- veghfpc[match(pt$oldID, rownames(veghfpc)),]
rownames(veghfpc) <- pt$ID
vegpc <- vegpc[match(pt$oldID, rownames(vegpc)),]
rownames(vegpc) <- pt$ID
veghfkm <- veghfkm[match(pt$oldID, rownames(veghfkm)),]
rownames(veghfkm) <- pt$ID
vegkm <- vegkm[match(pt$oldID, rownames(vegkm)),]
rownames(vegkm) <- pt$ID

summary(rowSums(veghfpc))
summary(rowSums(veghfkm))
table(rowSums(veghfpc) > 151^2*pi)
table(rowSums(veghfkm) > 10^6)

sit2 <- nonDuplicated(sit, ID, TRUE)
sit2 <- sit2[rownames(pt),]
sit2$ID <- NULL
pt <- data.frame(pt, sit2)

## wrapping up covariates
save(pt, veghfpc, vegpc, veghfkm, vegkm,
    file=file.path(ROOT, "data", "aru-coni", 
    "coni-compiled-covariates.Rdata"))

## wrapping up detections

compare_sets(pt$ID, int$ID)
ii <- intersect(pt$ID, int$ID)
pt <- droplevels(pt[ii,])
veghfpc <- veghfpc[ii,]
vegpc <- vegpc[ii,]
veghfkm <- veghfkm[ii,]
vegkm <- vegkm[ii,]
int <- droplevels(int[int$ID %in% ii,])
## times are only positives
compare_sets(pt$ID, tms$ID)
tms <- droplevels(tms[tms$ID %in% ii,])

save(pt, veghfpc, vegpc, veghfkm, vegkm, tms, int,
    file=file.path(ROOT, "data", "aru-coni", 
    "coni-compiled-all.Rdata"))



