source("~/repos/abmianalytics/veghf/veghf-setup.R")

### BAM+BBS bird points, 150 m radius buffer --------------------------------

## processing csv files
fl <- list.files(file.path(ROOT, VER, "data", "veghf", "bammbbs150m"))
tmplist <- list()
for (fn in fl) {
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "veghf", "bammbbs150m", fn)
    d <- read.csv(f)
    hfc <- "year"
    if (!(hfc %in% colnames(d)))
        hfc <- "YEAR"
    dd <- make_vegHF_wide(d, col.label="PKEY", 
        col.year="YEAR_", col.HFyear=hfc, sparse=TRUE)
    tmplist[[fn]] <- dd
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
}

## assembling return object
dd150m_bambbs <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year=NA,
    scale = "150 m radius circle around bird points")

### BAM+BBS bird points, 1 km^2 buffer

## processing csv files
fl <- list.files(file.path(ROOT, VER, "data", "veghf", "bammbbs564m"))
tmplist <- list()
for (fn in fl) {
    cat(which(fl == fn), "/", length(fl), "--", fn, "\n");flush.console()
    f <- file.path(ROOT, VER, "data", "veghf", "bammbbs564m", fn)
    d <- read.csv(f)
    hfc <- "year"
    if (!(hfc %in% colnames(d)))
        hfc <- "YEAR"
    dd <- make_vegHF_wide(d, col.label="PKEY", 
        col.year="YEAR_", col.HFyear=hfc, sparse=TRUE)
    tmplist[[fn]] <- dd
}

## binding together the pieces
veg_current <- tmplist[[1]]$veg_current
veg_reference <- tmplist[[1]]$veg_reference
soil_current <- tmplist[[1]]$soil_current
soil_reference <- tmplist[[1]]$soil_reference
for (j in 2:length(tmplist)) {
    cat("binding", j-1, "&", j, "/", length(tmplist), "\n");flush.console()
    veg_current <- bind_fun2(veg_current, tmplist[[j]]$veg_current)
    veg_reference <- bind_fun2(veg_reference, tmplist[[j]]$veg_reference)
    soil_current <- bind_fun2(soil_current, tmplist[[j]]$soil_current)
    soil_reference <- bind_fun2(soil_reference, tmplist[[j]]$soil_reference)
}

## assembling return object
dd1km_bambbs <- list(
    veg_current = veg_current,
    veg_reference = veg_reference,
    soil_current = soil_current,
    soil_reference = soil_reference,
    sample_year=NA,
    scale = "564 m radius circle around bird points")

climPoint_bambbs <- read.csv(
    file.path(ROOT, VER, "data/climate", "AllBird_fromPeter_april2015_climates.csv"))
colnames(climPoint_bambbs)[colnames(climPoint_bambbs) == "Eref"] <- "PET"
colnames(climPoint_bambbs)[colnames(climPoint_bambbs) == 
    "Populus_tremuloides_brtpred_nofp"] <- "pAspen"
climPoint_bambbs$OBJECTID <- NULL
colnames(climPoint_bambbs)[colnames(climPoint_bambbs) == "YEAR_"] <- "YEAR"
rownames(climPoint_bambbs) <- climPoint_bambbs$PKEY

compare.sets(rownames(dd150m_bambbs[[1]]), rownames(dd1km_bambbs[[1]]))
compare.sets(rownames(climPoint_bambbs), rownames(dd150m_bambbs[[1]]))

all(rownames(dd150m_bambbs[[1]]) == rownames(dd1km_bambbs[[1]]))
climPoint_bambbs <- climPoint_bambbs[rownames(dd150m_bambbs[[1]]),]
all(rownames(dd150m_bambbs[[1]]) == rownames(climPoint_bambbs))

## fix age 0 in saved files -----------------------------
load(file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

## dd150m_bambbs, dd1km_bambbs -- need NSR from previous climate table

load(file.path(ROOT, VER, "out/bambbs", "veg-hf_bambbs_fix-fire.Rdata"))

sum(dd150m_bambbs[[1]][,Target0])
dd150m_bambbs <- fill_in_0ages(dd150m_bambbs, climPoint_bambbs$NSRNAME)
sum(dd150m_bambbs[[1]][,Target0])

sum(dd1km_bambbs[[1]][,Target0])
dd1km_bambbs <- fill_in_0ages(dd1km_bambbs, climPoint_bambbs$NSRNAME)
sum(dd1km_bambbs[[1]][,Target0])

if (SAVE)
    save(dd150m_bambbs, dd1km_bambbs, climPoint_bambbs,
        file=file.path(ROOT, VER, "out/bambbs", "veg-hf_bambbs_fix-fire_fix-age0.Rdata"))

