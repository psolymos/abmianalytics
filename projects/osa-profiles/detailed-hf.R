library(mefa4)

e <- new.env()
load("e:/peter/sppweb2018/OSMR-kgrid_2017-11-14.Rdata", envir=e)
osa3 <- e[["km2.final.os.area"]]
e <- new.env()
load("e:/peter/sppweb2018/OSMR-minable-kgrid_2017-11-14.Rdata", envir=e)
mineable <- e[["km2.final.os.area"]]
e <- new.env()
load("e:/peter/sppweb2018/OilsandRegionDissolve.RData", envir=e)
ids <- list(merge=e$id)
ids$mineable <- as.character(mineable$LinkID)
ids$athabasca <- as.character(osa3$LinkID[osa3$Point.Poly == "Athabasca Oilsand Area"])
ids$coldlake <- as.character(osa3$LinkID[osa3$Point.Poly == "Cold Lake Oilsand Area"])
ids$peaceriver <- as.character(osa3$LinkID[osa3$Point.Poly == "Peace River Oilsand Area"])
rm(e)

if (FALSE) {
load(file.path("e:/peter/AB_data_v2018", "data", "analysis", "kgrid_table_km.Rdata")) # kgrid
load(file.path("e:/peter/AB_data_v2018", "data", "analysis", "grid",
    "veg-hf_grid_v6hf2016v3.Rdata")) # dd_kgrid
veg_cr <- t(sapply(ids, function(z) colSums(dd_kgrid[[1]][z,])))
veg_rf <- t(sapply(ids, function(z) colSums(dd_kgrid[[2]][z,])))
write.csv(veg_cr, file="Wall-to-wall_vegHF_for_OSR-current.csv")
write.csv(veg_rf, file="Wall-to-wall_vegHF_for_OSR-reference.csv")
}

load(file.path("e:/peter/AB_data_v2016", "out", "kgrid", "kgrid_table.Rdata")) # kgrid
load(file.path("e:/peter/AB_data_v2016",
    "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv$Sector2 <- tv$SectorRefined
vhf <- groupSums(dd1km_pred[[1]], 2, tv[colnames(dd1km_pred[[1]]), "Sector2"])
colnames(vhf)[colnames(vhf) == "NATIVE"] <- "Native"
AA <- t(sapply(ids, function(z) colSums(vhf[z,])/sum(colSums(vhf[z,]))))

CN <- c("CultivationCropPastureBareground",
    "Urban",
    "RuralResidentialIndustrial",
    "IndustrialSiteRural",
    "MineSite",
    "WellSite",
    "Pipeline",
    "SeismicLine",
    "TransmissionLine",
    "RoadHardSurface",
    "RoadTrailVegetated",
    "RoadVegetatedVerge",
    "RailHardSurface",
    "RailVegetatedVerge",
    "CCDecid",
    "CCMixwood",
    "CCConif",
    "CCPine")
CN2 <- c("Arg",
    "Urban",
    "RuralResidentialIndustrial",
    "IndustrialSiteRural",
    "MineSite",
    "WellSite",
    "Pipeline",
    "SeismicLine",
    "TransmissionLine",
    "RoadRail",
    "RoadRail",
    "RoadRail",
    "RoadRail",
    "RoadRail",
    "CCDecid",
    "CCMixwood",
    "CCConif",
    "CCPine")

fl <- list.files("e:/peter/sppweb2018/detailed-hf")

res <- list()

#i <- 1
for (i in 1:length(fl)) {
    cat(fl[i], "\n");flush.console()
    e <- new.env()
    load(file.path("e:/peter/sppweb2018/detailed-hf", fl[i]), envir=e)

    MAX <- max(quantile(rowSums(e$SA.Curr), 0.99), quantile(rowSums(e$SA.Ref), 0.99))
    cr0 <- rowSums(e$SA.Curr[intersect(ids$merge, rownames(e$SA.Curr)),])
    rf0 <- rowSums(e$SA.Ref[intersect(ids$merge, rownames(e$SA.Ref)),])
    MEAN_cr <- mean(cr0[cr0 <= quantile(cr0, 0.99)])
    MEAN_rf <- mean(rf0[rf0 <= quantile(rf0, 0.99)])
    MEAN <- max(MEAN_cr, MEAN_rf)

    cr <- t(sapply(ids, function(z) colSums(e$SA.Curr[intersect(z, rownames(e$SA.Curr)),])))
    rf <- t(sapply(ids, function(z) colSums(e$SA.Ref[intersect(z, rownames(e$SA.Curr)),])))
    rftot <- rowSums(rf)
    df <- cr - rf
    SEtot <- 100 * groupSums(df[,CN], 2, CN2) / rftot
    SEuhf <- 100 * groupSums(df[,CN], 2, CN2) / groupSums(rf[,CN], 2, CN2)
    #SEuhf[is.na(SEuhf)] <- 0
    SEuni <- SEtot / groupSums(AA[,CN], 2, CN2)
    SEuni[is.na(SEuni)] <- 0
    res[[substr(fl[i], 1, nchar(fl[i])-6)]] <- list(cr=cr, rf=rf,
        tot=SEtot, uhf=SEuhf, uni=SEuni, MAX=MAX, MEAN=MEAN)
}
save(res, file="e:/peter/sppweb2018/detailed-hf-birds.RData")


library(openxlsx)

keep <- (sapply(res, "[[", "MEAN") / sapply(res, "[[", "MAX")) >= 0.01
names(res)[!keep]

i <- "merge"
j <- "uni"
f <- function(i, j) {
    zz <- t(sapply(res[keep], function(z) z[[j]][i,]))
    zz <- data.frame(Species=rownames(zz), zz)
    zz[is.na(zz)] <- 0
    zz
}

xlslist <- list(
    OSA_Total=f("merge", "tot"),
    Mineable_Total=f("mineable", "tot"),
    Athabasca_Total=f("athabasca", "tot"),
    ColdLake_Total=f("coldlake", "tot"),
    PeaceRiver_Total=f("peaceriver", "tot"),
    OSA_UnderHF=f("merge", "uhf"),
    Mineable_UnderHF=f("mineable", "uhf"),
    Athabasca_UnderHF=f("athabasca", "uhf"),
    ColdLake_UnderHF=f("coldlake", "uhf"),
    PeaceRiver_UnderHF=f("peaceriver", "uhf"))

write.xlsx(xlslist, file="e:/peter/sppweb2018/detailed-hf-birds.xlsx", overwrite=TRUE)
