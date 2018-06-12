library(mefa4)
#library(cure4insect)

ROOT <- "e:/peter/AB_data_v2016"
#INDIR <- paste0("e:/peter/josm/2018/earlyseralTRUE")
#INDIR <- paste0("e:/peter/josm/2018/earlyseralFALSE")
INDIR <- paste0("e:/peter/josm/2018/hshfix")

sectors_all <- c("Agriculture", "EnergyLin", "EnergyMW", "Forestry", "Misc",
    "RuralUrban","Transportation")

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata")) # kgrid
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv$SectorForSeMs <- factor(ifelse(is.na(tv$SectorForSeMs), "NATIVE", as.character(tv$SectorForSeMs)),
    c("NATIVE", "Agriculture", "EnergyLin", "EnergyMW",
    "Forestry", "Misc", "RuralUrban","Transportation"))
stopifnot(all(rownames(kgrid) == rownames(dd1km_pred$veg_current)))
stopifnot(all(colnames(dd1km_pred$veg_current) == rownames(tv)))
HFarea <- groupSums(dd1km_pred$veg_current, 2, tv$SectorForSeMs)
rm(dd1km_pred)

## subset definition
kgrid$subset <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"
#kgrid$subset <- kgrid$POINT_Y > 50 & kgrid$NRNAME != "Grassland"

#rt <- .read_raster_template()
#rsub <- .make_raster(as.integer(kgrid$subset), kgrid, rt)
#plor(rsub)

## model for species
fln <- list.files(file.path(ROOT, "out", "birds", "results", "josmshf"))
fln <- sub("birds_abmi-josmshf_", "", fln)
fln <- sub(".Rdata", "", fln)

SPP <- fln
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL","RWBL")

spp <- "CAWA"

Ahf <- colSums(HFarea[kgrid$subset,])
Nsect <- list()
for (sect in c("All", sectors_all)) {
    e <- new.env()
    load(file.path(INDIR, sect, paste0(spp, ".RData")), envir=e)
    #stopifnot(all(rownames(kgrid) == rownames(e$SA.Curr)))
    Nsect[[sect]] <- cbind(cr=colSums(e$SA.Curr[kgrid$subset,]),
        rf=colSums(e$SA.Ref[kgrid$subset,]))
}

## reference abundance should not differ (much?)
round(100*sapply(Nsect, function(z) (z[,2]-Nsect[[1]][,2])/sum(Nsect[[1]][,2])), 3)
max(abs(round(100*sapply(Nsect, function(z) (z[,2]-Nsect[[1]][,2])/sum(Nsect[[1]][,2])), 3)))

Nref <- Nsect[[1]][,2]
Diffs <- sapply(Nsect, function(z) z[,1]-Nref)

## treating all indirects (over native and not)
#Tots <- cbind(All=Diffs[,1],
#    Conversion=c(0, diag(Diffs[-1,-1])), # direct
#    Disturbance=colSums(Diffs)-c(0, diag(Diffs[-1,-1]))) # indirect
## treating only native indirects
Tots <- cbind(All=Diffs[,1],
    Conversion=c(0, diag(Diffs[-1,-1])), # direct
    Disturbance=c(0, Diffs[1,-1])) # indirect
Tots["Native","Disturbance"] <- sum(Tots[,"All"])-sum(Tots[-1,-1]) # synergy
Tots <- rbind(Tots, Total=colSums(Tots))
round(100*Tots/sum(Nref), 2)

