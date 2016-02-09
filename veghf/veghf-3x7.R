source("~/repos/abmianalytics/veghf/veghf-setup.R")


## 3 x 7 yearly veg+HF

fl <- c("veg_hf_3x7_1999.csv","veg_hf_3x7_2001.csv",
    "veg_hf_3x7_2002.csv","veg_hf_3x7_2003.csv","veg_hf_3x7_2004.csv",
    "veg_hf_3x7_2005.csv","veg_hf_3x7_2006.csv","veg_hf_3x7_2007.csv",
    "veg_hf_3x7_2008.csv","veg_hf_3x7_2009.csv","veg_hf_3x7_2010.csv",
    "veg_hf_3x7_2011.csv","veg_hf_3x7_2012.csv","veg_hf_3x7_2013.csv",
    "veg_hf_3x7_2014.csv")
yr <- c(1999, 2001, 2002:2014)

yearly_vhf0 <- list()

#i <- 1
for (i in 1:length(fl)) {
    cat(i, "\n");flush.console()
    f <- file.path(ROOT, VER, "data/veghf/3x7", fl[i])
    d <- read.csv(f)

    d$inventory_year <- yr[i]
    #head(d)
    if (!("year" %in% colnames(d)))
        colnames(d)[colnames(d) == "YEAR"] <- "year"
    dd <- make_vegHF_wide(d, col.label = "ABMI", 
        col.year="inventory_year", col.HFyear="year")
    dd$scale <- "3x7 km rectangle around NSF site"
    dd$sample_year <- yr[i]
    #str(dd)

    yearly_vhf0[[as.character(yr[i])]] <- dd
}

## age0

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
rownames(gis) <- gis$SITE_ID

Sites <- rownames(gis)
NSR <- gis$NATURAL_SUBREGIONS

yearly_vhf <- yearly_vhf0

for (i in 1:length(fl)) {
    yearly_vhf0[[i]][[1]] <- yearly_vhf0[[i]][[1]][Sites,]
    yearly_vhf0[[i]][[2]] <- yearly_vhf0[[i]][[2]][Sites,]
    yearly_vhf0[[i]][[3]] <- yearly_vhf0[[i]][[3]][Sites,]
    yearly_vhf0[[i]][[4]] <- yearly_vhf0[[i]][[4]][Sites,]
    yearly_vhf[[i]] <- fill_in_0ages(yearly_vhf0[[i]], NSR)
}

lapply(yearly_vhf, function(zz) 
    sapply(zz[1:4], function(z) range(rowSums(z) / (3*7*10^6))))
lapply(yearly_vhf0, function(zz) sum(zz[[1]][,Target0]))
lapply(yearly_vhf, function(zz) sum(zz[[1]][,Target0]))

save(yearly_vhf,
    file=file.path(ROOT, VER, "out/3x7", 
    "veg-hf_3x7_fix-fire_fix-age0.Rdata"))

for (i in yr) {
cat(i, "\n");flush.console()
a1 <- data.frame(as.matrix(yearly_vhf[[as.character(i)]][[1]]))
a2 <- data.frame(as.matrix(yearly_vhf[[as.character(i)]][[2]]))
a3 <- data.frame(as.matrix(yearly_vhf[[as.character(i)]][[3]]))
a4 <- data.frame(as.matrix(yearly_vhf[[as.character(i)]][[4]]))
write.csv(a1, file=file.path(ROOT, VER, "out/3x7", paste0("veg-hf_3x7_", i, "_veg_current.csv")))
write.csv(a2, file=file.path(ROOT, VER, "out/3x7", paste0("veg-hf_3x7_", i, "_veg_reference.csv")))
write.csv(a3, file=file.path(ROOT, VER, "out/3x7", paste0("veg-hf_3x7_", i, "_soil_current.csv")))
write.csv(a4, file=file.path(ROOT, VER, "out/3x7", paste0("veg-hf_3x7_", i, "_soil_reference.csv")))
}
