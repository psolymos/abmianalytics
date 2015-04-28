ROOT <- "y:/Oracle_access_2015"
getwd()
source("R/00globalvars.R")

T <- "Birds"
if (do.prof) {
    proffile <- paste(D, "/OUT_", tolower(T), d, ".Rprof",sep="")
    Rprof(proffile)
}

if (FALSE) {## extablish connection
con <- odbcConnectAccess2007(DBVERSION)
## queries
res <- sqlFetch(con, "CSVDOWNLOAD_A_RT_BIRD_COUNT_V")
lookup <- sqlFetch(con, "RAWDATA_G_OFFGRID_SITE_LABEL")
taxo <- sqlFetch(con, "PUBLIC_ACCESS_PUBLIC_DETAIL_TAXONOMYS")
gis <- sqlFetch(con, "METADATA_GIS_SITE_SUMMARY")
## close connection
close(con)
}

res <- read.csv(file.path(ROOT, "data/birds.csv"))
gis <- read.csv(file.path(ROOT, "data/sitemetadata.csv"))
taxo <- read.csv(file.path(ROOT, "data/taxonomy.csv"))


## Labels etc
rr <- LabelFun(res)
rr <- nonDuplicated(rr, Label, TRUE)

res <- data.frame(res, rr[match(res$SITE_LABEL, rr$Label),])

## add in XY
res$PUBLIC_X <- gis$PUBLIC_LONGITUDE[match(res$ClosestABMISite, gis$SITE_ID)]
res$PUBLIC_Y <- gis$PUBLIC_LATTITUDE[match(res$ClosestABMISite, gis$SITE_ID)]

## add XY if required
if (FALSE) {
con <- odbcConnectAccess2007(DBVERSION)
gis <- sqlFetch(con, "METADATA_GIS_SITE_SUMMARY")
close(con)
res$PUBLIC_X <- gis$PUBLIC_LONGITUDE[match(res$ClosestABMISite, gis$SITE_ID)]
res$PUBLIC_Y <- gis$PUBLIC_LATTITUDE[match(res$ClosestABMISite, gis$SITE_ID)]
resxx <- res[!is.na(res$PUBLIC_X),c("SITE", "YEAR", "COMMON_NAME", 
    "SCIENTIFIC_NAME", "RANK_NAME", "TSNID", "Label", "PUBLIC_X", "PUBLIC_Y")]
summary(resxx)
write.csv(resxx, file="c:/p/GoA_ABMI_data_toTrish_20140113.csv")
}



table(res$BEHAVIOUR)

## exclude parts here if necessary
keep <- res$BEHAVIOUR %in% levels(res$BEHAVIOUR)
res <- res[keep,]

## this should all apply at point level! this is not 1 ha plot
## unique labels where species is DNC
tmp001 <- unique(res[res$SCIENTIFIC_NAME %in% c("DNC"),c("Label")])
## unique labels where species is NOT DNC
tmp002 <- unique(res[!(res$SCIENTIFIC_NAME %in% c("DNC")),c("Label")])
## use the collection status (does not work for 03-08 data)
#tmp003 <- unique(res[res$TMLIR_COLLECTION_STATUS %in% c("DNC","PNA"),c("Label")])
## unique labels where only DNC species occur (no data collected)
pcs.to.exclude <- setdiff(tmp001, tmp002)
## unique labels where collection status indicates so
#pcs.to.exclude <- union(pcs.to.exclude, tmp003)
#qs.to.exclude <- unique(as.character(res$Label[res$Label2 %in% pcs.to.exclude]))

## crosstab
res$SPECIES_OLD <- res$COMMON_NAME
levels(res$COMMON_NAME) <- nameAlnum(levels(res$COMMON_NAME), capitalize="mixed", collapse="")
res$COMMON_NAME <- droplevels(res$COMMON_NAME)

xt <- Xtab(~ Label + COMMON_NAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"), 
    rdrop=pcs.to.exclude, drop.unused.levels = FALSE)
#xt <- Xtab(~ Label + SCIENTIFIC_NAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"), 
#    subset = !(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")), drop.unused.levels = TRUE)
## get taxonomy
z <- nonDuplicated(res, res$COMMON_NAME, TRUE)
## add here higher taxa too
z <- z[colnames(xt),c("TSNID","COMMON_NAME","SCIENTIFIC_NAME","RANK_NAME","SPECIES_OLD")]
#z[] <- lapply(z, function(z) z[drop=TRUE])
#summary(z)
## trim trailing white space -- just to make sure
#for (i in 1:ncol(z))
#    if (is.factor(z[[i]]))
#        levels(z[[i]]) <- sub(' +$', '', levels(z[[i]]))
z <- droplevels(z)

## crosstab
#xt <- Xtab(~ Label + SCIENTIFIC_NAME, res)
## get taxonomy
#z <- nonDuplicated(res, res$SCIENTIFIC_NAME, TRUE)
## add here higher taxa too
#z <- z[,c("TSNID","COMMON_NAME","SCIENTIFIC_NAME","RANK_NAME")]

x <- nonDuplicated(res[,c("Label", "Label2", "ROTATION", "SITE", "YEAR", "ADATE", "CREW", 
    "TBB_POINT_COUNT", "WIND_CONDITION", "PRECIPTATION", 
    "TBB_START_TIME", "TBB_END_TIME", "OnOffGrid", 
    "SiteLabel", "DataProvider", "Visit", "ClosestABMISite", "PUBLIC_X", "PUBLIC_Y")], res$Label, TRUE)

## crosstab on PC level
m <- Mefa(xt, x, z)
## exclude not species level taxa
table(m@taxa$RANK_NAME) # no sub-specific levels for birds
#m <- m[,m@taxa$RANK_NAME == "Species"] # -- Do not, do it by Order
## exclude squirrels and assoc. -- Don't
#m <- m[,m@taxa$CLASS_NAME == "Aves"]

## excl. waterfowls and pirds of prey
#m <- m[,!(toupper(m@taxa$ORDER_NAME) %in% c("ANSERIFORMES","CHARADRIIFORMES",
#    "FALCONIFORMES","PELECANIFORMES","PODICIPEDIFORMES","STRIGIFORMES",
#    "GAVIIFORMES","CICONIIFORMES"))]

#aaa <- lapply(structure(levels(m@taxa$FAMILY),
#    names=toupper(levels(m@taxa$FAMILY))), function(i) {
#    tmp <- m@taxa$COMMON_NAME[m@taxa$FAMILY == i]
#    tmp[!is.na(tmp)]
#})

if (FALSE) {
m <- m[,!(toupper(m@taxa$FAMILY) %in% c(
    "ACCIPITRIDAE",
    "ANATIDAE",
    "ARDEIDAE",
    "CANIDAE", # Coyote
    "CAPRIMULGIDAE",
    "FALCONIDAE",
    "GAVIIDAE",
    "GRUIDAE",
    "LARIDAE",
    "PELECANIDAE",
    "PHALACROCORACIDAE",
    "PODICIPEDIDAE",
    "STRIGIDAE",
    "TYTONIDAE")) | m@taxa$COMMON_NAME == "Sandhill Crane"]
m@taxa[] <- lapply(m@taxa, function(z) z[drop=TRUE])
}

## site level info
m2 <- groupSums(m, 1, m@samp$Label2)
#tmp <- Xtab(~Label2 + TBB_POINT_COUNT, m@samp)

## this is 0/1 at PC level
m01 <- m
xtab(m01)[xtab(m01)>0] <- 1
m012 <- groupSums(m01, 1, m@samp$Label2)

## crosstabs
if (!combine.tables) {
    res1 <- as.matrix(m)
    res2 <- as.matrix(m2)
    res3 <- as.matrix(m012)
    rntw <- TRUE
} else {
    res1 <- data.frame(samp(m), as.matrix(m))
    mmm <- Mefa(xtab(m2), nonDuplicated(samp(m)[,-which(colnames(samp(m))=="Label")], Label2, TRUE))
    res2 <- data.frame(samp(mmm), as.matrix(mmm))
    mmm01 <- Mefa(xtab(m012), nonDuplicated(samp(m)[,-which(colnames(samp(m))=="Label")], Label2, TRUE))
    tmp <- table(m@samp$Label2)
    samp(mmm01)$TotalNoOfPC <- tmp[match(samp(mmm01)$Label2, names(tmp))]
    res3 <- data.frame(samp(mmm01), as.matrix(mmm01))
    rntw <- FALSE
}
tax <- taxa(m)
#sam1 <- samp(m)
#sam2 <- samp(m)
tax[] <- lapply(tax, function(z) z[drop=TRUE])
str(res1)
str(res2)
str(res3)
m
m2
str(tax)
#str(sam1)
#str(sam2)
range(xtab(m))
range(xtab(m2))
range(xtab(m012))

if (check.completeness) {
    mln <- paste(D, "/OUT_", "MasterList", "_SpringAndSummerSite", d, ".csv",sep="")
    if (file.exists(mln)) {
        ml <- read.csv(mln)
        ml$SiteFound <- ml$Label %in% rownames(res2)
        cat("\nSites per year:\n\n")
        print(with(ml, table(Year=Year, SiteFound=SiteFound)))
        cat("\nSites per year and per data provider:\n\n")
        print(with(ml, table(Year=Year, SiteFound=SiteFound, DataProvider=DataProvider)))
        cat("\nSites not found in ", T, "data set:\n\n")
        print(data.frame(NotFound=sort(as.character(ml$Label)[!ml$SiteFound])))
        if (sum(!ml$SiteFound)) {
            con <- odbcConnectAccess2007(DBVERSION)
            ## queries
            id0 <- sqlQuery(con, paste("SELECT * FROM SITE_MONITOR_SITE_SUMMARY_PUBLIC_V"))
            ## close connection
            close(con)
            PROT <- c("T26 Breeding Birds")
            id0 <- id0[id0$PROTOCOL_NAME %in% c("Terrestrial Spring","Terrestrial Summer"),]
            id0$SiteYear <- as.factor(with(id0, paste(SITE, SYML_YEAR, sep="_")))
            id <- id0[id0$SiteYear %in% as.character(ml$SiteYear)[!ml$SiteFound],]
            id <- id[id$DPL_COLLECTION_NAME %in% PROT,]
            id <- data.frame(id, ml[match(id$SiteYear, ml$SiteYear),])
            write.csv(id, file=paste(D, "PROBLEM_SITES_", T, d, ".csv",sep=""), row.names = rntw)
        }
    } else {
        cat("no master list found \n\n")
    }
}

## write static files into csv
write.csv(res1, file=paste(D, "/OUT_", T, "_Species_PC-Counts", d, ".csv",sep=""), row.names = rntw)
write.csv(res2, file=paste(D, "/OUT_", T, "_Species_Site-Counts", d, ".csv",sep=""), row.names = rntw)
write.csv(res3, file=paste(D, "/OUT_", T, "_Species_Site-Binomial", d, ".csv",sep=""), row.names = rntw)
write.csv(tax, file=paste(D, "/OUT_", T, "_Species_Taxa", d, ".csv",sep=""), row.names = TRUE)
write.csv(data.frame(Excluded_Points=sort(pcs.to.exclude)),
    file=paste(D, "/OUT_", T, "_Excluded_Sites", d, ".csv",sep=""), row.names = FALSE)

if (do.prof)
    summaryRprof(proffile)
if (do.image)
    save.image(paste(D, "/OUT_", tolower(T), d, ".Rdata",sep=""))
## quit without saving workspace
quit(save="no")

## compare abundances across years
tmp <- groupMeans(xtab(m), 1, samp(m)$YEAR)
tmp <- t(as.matrix(tmp))
matplot(log(t(tmp)+1), type="l", col=1, lty=1)
tmp2 <- ifelse(tmp>0, 1, 0)
tmp <- tmp[order(rowMeans(tmp), decreasing=TRUE),order(as.integer(colnames(tmp)))]

round(tmp[tmp[,"2013"]==0,],2)
aa <- droplevels(res[res$YEAR==2013,])
bb <- droplevels(res[res$YEAR!=2013,])
setdiff(levels(aa$SPECIES_OLD), levels(bb$SPECIES_OLD))
setdiff(levels(bb$SPECIES_OLD), levels(aa$SPECIES_OLD))
intersect(levels(bb$SPECIES_OLD), levels(aa$SPECIES_OLD))

## add BAM stuff here ----------------------------------

#sppc <- read.csv("y:/Oracle_access/src2/species_codes.csv")

library(RODBC)
con <- odbcConnectAccess2007("c:/bam/BAM_NATIONAL_V2_2011_TEAM_MARCH11.accdb")
taxo2 <- sqlQuery(con, "SELECT * FROM [DD_Combined Life History];")
close(con)

mm <- m[, !is.na(m@taxa$CLASS_NAME)]
mm <- mm[, mm@taxa$CLASS_NAME=="Aves"]
setdiff(colnames(mm), taxo2[,"SCIENTIFIC NAME"])



## detectability
