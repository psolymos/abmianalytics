library(sf)
library(mefa4)

## read in data (from Google drive)
ROOT <- "~/Desktop/x/Human Footprint Data Extraction"

a <- st_read(file.path(ROOT, "sta_latlong_10TMf_buff5km.gdb"))
b <- st_read(file.path(ROOT, "MAPS_birdVitals_HFI2018clp.gdb"))
b <- st_cast(b,  "MULTIPOLYGON")

## site locations
xy <- cbind(X=a$utm_e, Y=a$utm_n)
rownames(xy) <- a$name
z <- st_multipoint(xy)
z <- st_sfc(z)
st_crs(z) <- st_crs(a)

## total area in 5km radius buffers
Atot <- 5000^2*pi

## footprint lookup:
## https://github.com/psolymos/abmianalytics/blob/master/lookup/lookup-hf-type-v2014.csv
hf <- read.csv("~/repos/abmianalytics/lookup/lookup-hf-type-v2014.csv")
hf <- hf[!duplicated(hf[,1]),]
rownames(hf) <- hf[,1]

## loop over sites and intersect with buffer
## gdb did not have site ID, that's why we need this
i <- 1
out <- NULL
for (i in 1:nrow(xy)) {
    cat(i, "\n")
    flush.console()
    zi <- st_point(xy[i,])
    v <- st_buffer(zi, 5000)
    v <- st_sfc(v)
    st_crs(v) <- st_crs(a)

    ii <- st_intersection(b, v) # this is why we need as_cast above, otherwise this fails

    tmp <- data.frame(
        hf=ii$FEATURE_TY,
        p=st_area(ii)/Atot,
        site=factor(rownames(xy)[i], rownames(xy)))
    out <- rbind(out, tmp)
}

## one feature type category is undefined
str(out)
levels(out$hf)[levels(out$hf) == "CONVENTIONAL-SEISMIC"] <- "CUTLINE_TRAIL"

## pivot table
y <- Xtab(p ~ site + hf, out)

## check feature types: we aren't missing any
compare_sets(colnames(y), rownames(hf))
setdiff(colnames(y), rownames(hf))

## aggregate (add up columns in the same class)
y <- groupSums(y, 2, hf[colnames(y), "HF_GROUP"])
y <- as.matrix(y) # sparese to dense matrix coercion

## checks
rowSums(y)
summary(y)

## save
yy <- data.frame(Site=rownames(y), y)
write.csv(yy, row.names=FALSE, file=file.path(ROOT, "HF_class_proportions.csv"))
