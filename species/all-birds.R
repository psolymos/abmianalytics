
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!! this file focuses on revisits only for now !!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(mefa4)
source("~/repos/abmianalytics/species/abmi-r-api.R")
#data.frame(table=get_table_names())

## settings
taxon <- "birds"
ROOT <- "d:/abmi/AB_data_v2018"

## common stuff
DATE <- as.Date(Sys.time(), tz=Sys.timezone(location = TRUE))
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
add_labels <- function(res, sub_col) {
    res$offgrid <- startsWith(as.character(res$ABMISite), "OG")
    res$subunit <- res[[sub_col]]
    res[[sub_col]] <- NULL
    res$site_year <- interaction(res$ABMISite, res$Year, drop=TRUE, sep="_")
    res$site_year_sub <- interaction(res$ABMISite, res$Year, res$subunit, drop=TRUE, sep="_")

    tmp <- strsplit(as.character(res$ABMISite), "-")
    res$nearest <- sapply(tmp, function(z) {
        zz <- if (length(z) > 1) z[3] else z[1]
        as.integer(gsub("\\D+", "", zz))
    })
    res
}
normalize_species <- function(res) {
#    res$ScientificName0 <- res$ScientificName
    levels(res$ScientificName) <- gsub("X ", "", levels(res$ScientificName))
    levels(res$ScientificName) <- gsub(" x ", " ", levels(res$ScientificName))
    levels(res$ScientificName) <- sapply(strsplit(levels(res$ScientificName), " "), function(z) {
        paste(z[1:min(2, length(z))], collapse=" ")
    })

    levels(res$TaxonomicResolution)[levels(res$TaxonomicResolution) %in%
        c("Subspecies", "Variety")] <- "Species"

    res$SpeciesID <- res$ScientificName
    levels(res$SpeciesID) <- nameAlnum(levels(res$SpeciesID), capitalize="mixed", collapse="")
    res$SpeciesID <- droplevels(res$SpeciesID)

    res
}
find_max <- function(x) {
    if (is.null(dim(x)))
        stop("x must be matrix like object with dim attribute")
    if (is.null(colnames(x)))
        colnames(x) <- paste0("X", seq_len(ncol(x)))
    m <- ncol(x)
    idx <- integer(nrow(x))+1L
    val <- x[,1L]
    for (j in seq_len(m)[-1L]) {
        s <- x[,j] > val
        val[s] <- x[s,j]
        idx[s] <- j
    }
    i <- factor(idx, levels = seq_len(m))
    levels(i) <- colnames(x)
    out <- data.frame(index = i, value = val)
    rownames(out) <- rownames(x)
    out
}
find_min <- function (x)
{
    out <- find_max(-1 * x)
    out$value <- -1 * out$value
    out
}
cn1 <- c("ABMISite", "Year", "subunit", "site_year", "site_year_sub", "offgrid", "nearest")
cn2 <- c("SpeciesID", "CommonName", "ScientificName", "TaxonomicResolution",
    "UniqueTaxonomicIdentificationNumber")

resRF0 <- get_table("T26A")
resSM0 <- get_table("T26B")
resRF <- resRF0
resSM <- resSM0
table(resRF$Year)
table(resSM$Year)

save(resRF0, resSM0, file=file.path(ROOT, "data", "raw", "species",
    paste0(taxon, "_", DATE, ".Rdata")))

## River Forks, 3 intervals ------------------------------

colnames(resRF) <- gsub(" ", "", colnames(resRF))
resRF <- add_labels(resRF, sub_col="PointCountStation")
#resRF <- normalize_species(resRF)
resRF$SpeciesID <- resRF$CommonName
levels(resRF$SpeciesID) <- nameAlnum(levels(resRF$SpeciesID), capitalize="mixed", collapse="")
resRF$SpeciesID <- droplevels(resRF$SpeciesID)

resRF$TBB_TIME_1ST_DETECTED <- resRF[["TimeFirstDetected(seconds)"]]
resRF$TBB_INTERVAL_1 <- resRF[["Interval1(0-200seconds)"]]
resRF$TBB_INTERVAL_2 <- resRF[["Interval2(201-400seconds)"]]
resRF$TBB_INTERVAL_3 <- resRF[["Interval3(401-600seconds)"]]
resRF$TBB_TIME_1ST_DETECTED <- as.character(resRF$TBB_TIME_1ST_DETECTED)
table(resRF$TBB_TIME_1ST_DETECTED,resRF$Year)
resRF$TBB_TIME_1ST_DETECTED[resRF$TBB_TIME_1ST_DETECTED %in%
    c("DNC", "NONE", "VNA")] <- NA
resRF$TBB_TIME_1ST_DETECTED <- as.numeric(as.character(resRF$TBB_TIME_1ST_DETECTED))
resRF$period1st <- as.numeric(cut(resRF$TBB_TIME_1ST_DETECTED, c(-1, 200, 400, 600)))

resRF$keep <-
    resRF$TBB_INTERVAL_1 %in% c("0","1") &
    resRF$TBB_INTERVAL_2 %in% c("0","1") &
    resRF$TBB_INTERVAL_3 %in% c("0","1")
#resRF <- resRF[resRF$TBB_INTERVAL_1 %in% c("0","1"),]
#resRF <- resRF[resRF$TBB_INTERVAL_2 %in% c("0","1"),]
#resRF <- resRF[resRF$TBB_INTERVAL_3 %in% c("0","1"),]
resRF$TBB_INTERVAL_1 <- as.integer(resRF$TBB_INTERVAL_1) - 1L
resRF$TBB_INTERVAL_2 <- as.integer(resRF$TBB_INTERVAL_2) - 1L
resRF$TBB_INTERVAL_3 <- as.integer(resRF$TBB_INTERVAL_3) - 1L

tmp <- col(matrix(0,sum(resRF$keep),3)) *
    resRF[resRF$keep,c("TBB_INTERVAL_1","TBB_INTERVAL_2","TBB_INTERVAL_3")]
tmp[tmp==0] <- NA
tmp <- cbind(999,tmp)
resRF$period123 <- 999
resRF$period123[resRF$keep] <- apply(tmp, 1, min, na.rm=TRUE)
resRF$period123[resRF$period123 > 3] <- 999
with(resRF, table(period1st, period123))
resRF$period1 <- pmin(resRF$period1st, resRF$period123)
with(resRF, table(period1st, period1))
with(resRF, table(period123, period1))


## SM units --------------------------------------------

colnames(resSM) <- gsub(" ", "", colnames(resSM))

#cn <- c("ABMISite", "Year", "Quadrant", "Method",
#    "Interval1(1minute)", "Interval2(1minute)", "Interval3(1minute)",
#    "CommonName", "RecordingDate", "RecordingTime")
#resSM <- resSM[,cn]
tmp <- paste(resSM$RecordingDate, resSM$RecordingTime)
resSM$Start <- strptime(tmp, "%d-%b-%y %H:%M:%S")

resSM$Duration <- NA
resSM$Duration[resSM$Method %in% c("11", "14")] <- 3
resSM$Duration[resSM$Method %in% c("12", "13")] <- 1

resSM <- add_labels(resSM, sub_col="Quadrant")
#resSM <- normalize_species(resSM)
resSM$SpeciesID <- resSM$CommonName
levels(resSM$SpeciesID) <- nameAlnum(levels(resSM$SpeciesID), capitalize="mixed", collapse="")
resSM$SpeciesID <- droplevels(resSM$SpeciesID)

## first detection interval
resSM$int1 <- ifelse(resSM$`Interval1(1minute)` == "VNA", NA, as.integer(resSM$`Interval1(1minute)`))
resSM$int2 <- ifelse(resSM$`Interval2(1minute)` == "VNA", NA, as.integer(resSM$`Interval2(1minute)`))
resSM$int3 <- ifelse(resSM$`Interval3(1minute)` == "VNA", NA, as.integer(resSM$`Interval3(1minute)`))
tmp <- col(resSM[,c("int1", "int2", "int3")])
tmp[is.na(resSM[,c("int1", "int2", "int3")])] <- Inf
tmp2 <- find_min(tmp)
tmp2$value[is.infinite(tmp2$value)] <- NA
resSM$res1 <- tmp2$value

resSM$ToY <- resSM$Start$yday
resSM$ToYc <- as.integer(cut(resSM$ToY, c(0, 105, 120, 140, 150, 160, 170, 180, 365)))

resSM$visit <- interaction(resSM$site_year_sub, resSM$ToYc, drop=TRUE)

resSM$ToD <- resSM$Start$hour + resSM$Start$min / 60
resSM$ToDx <- round(resSM$ToD, 0)
resSM$ToDc <- as.factor(ifelse(resSM$ToDx == 0, "Midnight", "Morning"))

## crosstabs

cd <- c("NONE","SNI", "VNA", "DNC", "PNA")

xt_10 <- Xtab(~ site_year_sub + SpeciesID,
    resRF[resRF$period123 <= 3,], cdrop=cd)

xt_pts <- Xtab(~ site_year_sub + SpeciesID + period123, resRF, cdrop=cd)[1:3]
xt_pts[[1]] <- as.matrix(xt_pts[[1]])
xt_pts[[2]] <- as.matrix(xt_pts[[2]])
xt_pts[[3]] <- as.matrix(xt_pts[[3]])
x_pts <- nonDuplicated(resRF, site_year_sub, TRUE)
x_pts <- x_pts[rownames(xt_pts[[1]]),]

xt_stn <- as.matrix(Xtab(~ site_year_sub + SpeciesID, resSM, cdrop=cd))
x_stn <- nonDuplicated(resSM, site_year_sub, TRUE)
x_stn <- x_stn[rownames(xt_stn),]

keep <- !is.na(resSM$visit)
xt_vis <- as.matrix(Xtab(~ visit + SpeciesID, resSM[keep,], cdrop=cd))
x_vis <- nonDuplicated(resSM[keep,], visit, TRUE)
x_vis <- x_vis[rownames(xt_vis),]

save(xt_pts, xt_stn, xt_vis, x_pts, x_stn, x_vis, xt_10, resSM, resRF,
    file=file.path(ROOT, "data", "inter", "species", "birds-revisit.Rdata"))

## trend piece --------------------------------------


compare_sets(x_pts$ABMISite, x_stn$ABMISite)
setdiff(x_stn$ABMISite, x_pts$ABMISite)
setdiff(x_pts$ABMISite, x_stn$ABMISite)

SITES <- intersect(droplevels(x_vis$ABMISite), x_pts$ABMISite)
SPP <- intersect(colnames(xt_vis), colnames(xt_pts[[1]]))
BB <- cbind(1:length(SITES), replicate(199, sample(length(SITES), replace=TRUE)))


#site <- "1"
#spp <- "Ovenbird"
one_site <- function(site, spp) {
    i1 <- !is.na(x_pts$ABMISite) & x_pts$ABMISite == site
    x1 <- x_pts[i1,]
    y1 <- as.matrix(xt_pts[[1]])[i1,spp]
    t1 <- min(unique(x1$Year))
    y1 <- y1[x1$Year == t1]

    i2 <- x_vis$ABMISite == site & x_vis$Method %in% c("11", "14")
    x2 <- x_vis[i2,]
    y2 <- xt_vis[i2,spp]
    t2 <- max(unique(x2$Year))
    y2 <- y2[x2$Year == t2]

    data.frame(site=site, species=spp,
        yr1=t1, yr2=t2, dt=t2-t1, y1=mean(y1), y2=mean(y2))
}
get_pac <- function(zz) {
    TC <- with(zz, sum(y2) / sum(y1))
    WM <- with(zz, sum(y1 * dt) / sum(y1))
    lam <- TC^(1/WM)
    c(TC=TC, WM=WM, lambda=lam, pac=100*(lam-1))
}
one_species <- function(spp) {
    zz <- do.call(rbind, lapply(SITES, one_site, spp=spp))
    zzz <- t(apply(BB, 2, function(z) get_pac(zz[z,])))
    cbind(Estimate=zzz[1,], t(apply(zzz, 2, quantile, c(0.05, 0.95), na.rm=TRUE)))
}

#for (i in SPP) one_species(i)
all <- pbapply::pblapply(SPP, one_species)
names(all) <- SPP

save(all, file="~/GoogleWork/josm/2018/trend/josm-trend-results-2018july.RData")

pacs <- t(sapply(all, function(z) z["pac",]))
pacs <- pacs[!is.na(pacs[,3]) & pacs[,1] > -100 & pacs[,1] < 100,]
pacs <- pacs[order(pacs[,1]),]

write.csv(pacs, file="~/GoogleWork/josm/2018/trend/josm-trend-results-2018july.csv")

library(ggplot2)
d <- data.frame(Species=as.factor(rownames(pacs)), pacs)
colnames(d)[3:4] <- c("lwr", "upr")
p <- ggplot(aes(x=Estimate, y=Species), data=d) +
    geom_point() +
    NULL


