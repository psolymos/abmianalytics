ROOT <- "e:/peter/AB_data_v2016/oracle"
OUTDIR <- "e:/peter/AB_data_v2016/data/species"
getwd()
if (interactive())
    source("~/repos/abmianalytics/species/00globalvars.R") else source("00globalvars.R")

T <- "Mites"
if (do.prof) {
    proffile <- paste(D, "/OUT_", tolower(T), d, ".Rprof",sep="")
    Rprof(proffile)
}

if (FALSE) {## extablish connection
con <- odbcConnectAccess2007(DBVERSION)
## queries
res <- sqlQuery(con, paste("SELECT * FROM CSVDOWNLOAD_A_RT_SOIL_MITES_V"))
#lookup <- sqlQuery(con, paste("SELECT * FROM RAWDATA_G_OFFGRID_SITE_LABEL"))
taxo <- sqlQuery(con, paste("SELECT * FROM PUBLIC_ACCESS_PUBLIC_DETAIL_TAXONOMYS"))
gis <- sqlFetch(con, "METADATA_GIS_SITE_SUMMARY")
#check <- sqlQuery(con, paste("SELECT * FROM CSVDOWNLOAD_RT_SOIL_ARTHROPOD_C_DOWNLOAD"))
## close connection
close(con)
}

res <- read.csv(file.path(ROOT, "mites.csv"))
#gis <- read.csv(file.path(ROOT, "data/sitemetadata.csv"))
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
taxo <- read.csv(file.path(ROOT, "taxonomy.csv"))

if (FALSE) {
## mite check, Dave Apr 25, 2014

##- Two sites were in the mite file but not in the master list: 
##  OG-ABMI-15-1 had a row in the mite file but no mites, so I just deleted it.  
#new site where no species was found (all NONE for all 4 quadrants)
##  OG-ABMI-1280-1B has a row, with 21 mite occurrences.  
##  The results are very similar but not identical to OG-ABMI-1280-1.  
##  So, not sure what the "B" site is.
#Ask IC

##- All sampled sites that were missing mite data seemed to be in lakes, so that's fine.
#OK
##- There is one species Anachipteria.sp..1.DEW and one Anachipteria.sp.1.DEW 
##  (two dots after "sp" versus one dot).  The second one is the "non-standard" name, 
##  because it is missing the "." after "sp".  
#keep "Anachipteria sp. 1 DEW", merge "Anachipteria sp.1 DEW" 

##- There is a Cepheus.sp..2.DEW with lots of records, and a Cepheus.sp..2B.DEW 
##  with one record.  Not sure if that is supposed to be a separate new species?
#Ask PC/IC how to treat  "Cepheus sp. 2 DEW" and "Cepheus sp. 2B DEW"

##- There is a species "Eueremaeus.marshalli...Eueremaeus.quadrilamellatus" 
##  don't know if that's a typo, or just that the two options couldn't be distinguished.
#"Eueremaeus marshalli / Eueremaeus quadrilamellatus" -- ask IC/PC is this should be 
#treated as Eueremaeus marshalli

##- There is a species "Eupterotegaeus.rostratus." as well as the normal 
##  name without the dot at the end.  Must be a trailing blank space.
#keep "Eupterotegaeus rostratus", merge "Eupterotegaeus rostratus "

##- There is a Tectoribates.sp..1.DEW with 1 record, and a 
## Tectoribates.sp..B.DEW with 6 records.  Not sure if that is intentional?
#Ask IC/PC how to treat "Tectoribates sp. 1 DEW" and "Tectoribates sp. B DEW"

}


## Labels etc
rr <- LabelFun(res)
rr <- nonDuplicated(rr, Label, TRUE)

res <- data.frame(res, rr[match(res$SITE_LABEL, rr$Label),])
## add in XY
res$PUBLIC_X <- gis$PUBLIC_LONGITUDE[match(res$ClosestABMISite, gis$SITE_ID)]
res$PUBLIC_Y <- gis$PUBLIC_LATTITUDE[match(res$ClosestABMISite, gis$SITE_ID)]

## exclude parts here if necessary
keep <- 1:nrow(res)
res <- res[keep,]

## this should all apply at site level!
## unique labels where species is DNC
tmp001 <- unique(res[res$SCIENTIFIC_NAME %in% c("DNC"),c("Label2")])
## unique labels where species is NOT DNC
tmp002 <- unique(res[!(res$SCIENTIFIC_NAME %in% c("DNC")),c("Label2")])
## use the collection status (does not work for 03-08 data)
#tmp003 <- unique(res[res$TMLIR_COLLECTION_STATUS %in% c("DNC","PNA"),c("Label2")])
## unique labels where only DNC species occur (no data collected)
pcs.to.exclude <- setdiff(tmp001, tmp002)
## unique labels where collection status indicates so
#pcs.to.exclude <- union(pcs.to.exclude, tmp003)
qs.to.exclude <- unique(as.character(res$Label[res$Label2 %in% pcs.to.exclude]))

## crosstab
res$count <- as.character(res$COUNTNUM)
res$count[res$count %in% c("NONE","SNI", "VNA", "DNC", "PNA")] <- "0"
res$count <- as.integer(res$count)

res$SPECIES_OLD <- res$SCIENTIFIC_NAME
levels(res$SCIENTIFIC_NAME) <- nameAlnum(levels(res$SCIENTIFIC_NAME), capitalize="mixed", collapse="")
res$SCIENTIFIC_NAME <- droplevels(res$SCIENTIFIC_NAME)

xt <- Xtab(count ~ Label + SCIENTIFIC_NAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"), 
    rdrop=qs.to.exclude, drop.unused.levels = FALSE)
#xt <- Xtab(~ Label + SCIENTIFIC_NAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"), 
#    subset = !(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")), drop.unused.levels = TRUE)

## get taxonomy
z <- nonDuplicated(res[!(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")),],
    res$SCIENTIFIC_NAME[!(res$SCIENTIFIC_NAME %in% c("VNA", "DNC"))], TRUE)
## add here higher taxa too
z <- z[,c("TSNID","SCIENTIFIC_NAME","SPECIES_OLD","RANK_NAME")]
z2 <- taxo[taxo$SCIENTIFIC_NAME %in% z$SCIENTIFIC_NAME,]
z <- data.frame(z, z2[match(z$SCIENTIFIC_NAME, z2$SCIENTIFIC_NAME),setdiff(colnames(taxo), colnames(z))])
#z[] <- lapply(z, function(z) z[drop=TRUE])
#summary(z)
## trim trailing white space -- just to make sure
for (i in 1:ncol(z))
    if (is.factor(z[[i]]))
        levels(z[[i]]) <- sub(' +$', '', levels(z[[i]]))

x <- nonDuplicated(res[,c("Label", "Label2", "ROTATION", "SITE", "YEAR", "ADATE", "CREWS", 
    "TSM_QUADRANT", "OnOffGrid", 
#    "OGLabel", 
    "SiteLabel", "DataProvider", "Visit", "ClosestABMISite", "PUBLIC_X", "PUBLIC_Y")], res$Label, TRUE)

## crosstab on PC level
m <- Mefa(xt, x, z)
## exclude not species level taxa
table(m@taxa$RANK_NAME) # here are sub-specific levels
## here seems that all subspecies are different species, but just to be sure...
#m <- m[,m@taxa$RANK_NAME %in% c("Species", "Subspecies")]
m <- m[,m@taxa$RANK_NAME %in% c("Genus", "Species", "Subspecies")]

#tmp <- Xtab(~Label2 + TSM_QUADRANT, m@samp)

## here I want to collapse sub-specific taxa, but TAX WB is not that clean to map these
#zz <- m@taxa[!is.na(m@taxa$SPECIES),]
#zz <- nonDuplicated(zz, zz$SPECIES, TRUE)
#m <- groupSums(m[,!is.na(m@taxa$SPECIES)], 2, m@taxa$SPECIES[!is.na(m@taxa$SPECIES)],
#    replace=zz)

## site level info
m2 <- groupSums(m, 1, m@samp$Label2)

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
    samp(mmm01)$TotalNoOfQU <- tmp[match(samp(mmm01)$Label2, names(tmp))]
    res3 <- data.frame(samp(mmm01), as.matrix(mmm01))
    rntw <- FALSE
}
tax <- taxa(m)
tax[] <- lapply(tax, function(z) z[drop=TRUE])
str(res1)
str(res2)
str(res3)
m
m2
str(tax)
range(xtab(m))
range(xtab(m2))
range(xtab(m012))

if (check.completeness) {
    mln <- paste(D, "OUT_", "MasterList", "_SpringAndSummerSite", d, ".csv",sep="")
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
            PROT <- c("T23 Soil Arthropods (Springtails and Mites) Collection",
                "T24A Soil Arthropods (Mites) Identification")
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
write.csv(res1, file=paste(D, "/OUT_", T, "_Species_QU-Counts", d, ".csv",sep=""), row.names = rntw)
write.csv(res2, file=paste(D, "/OUT_", T, "_Species_Site-Counts", d, ".csv",sep=""), row.names = rntw)
write.csv(res3, file=paste(D, "/OUT_", T, "_Species_Site-Binomial", d, ".csv",sep=""), row.names = rntw)
write.csv(tax, file=paste(D, "/OUT_", T, "_Species_Taxa", d, ".csv",sep=""), row.names = TRUE)
write.csv(data.frame(Excluded_Sites=sort(pcs.to.exclude)),
    file=paste(D, "/OUT_", T, "_Excluded_Sites", d, ".csv",sep=""), row.names = FALSE)

if (do.prof)
    summaryRprof(proffile)
if (do.image)
    save.image(paste(D, "/OUT_", tolower(T), d, ".Rdata",sep=""))
## quit without saving workspace
quit(save="no")

