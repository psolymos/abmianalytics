source("y:/Oracle_access/src2/00globalvars.R")
T <- "AqPlants"
if (do.prof) {
    proffile <- paste(D, "OUT_", tolower(T), d, ".Rprof",sep="")
    Rprof(proffile)
}
## load mefa4 package
stopifnot(require(mefa4))
## load package RODBC
stopifnot(require(RODBC))
## extablish connection
con <- odbcConnectAccess2007(DBVERSION)
## queries
res <- sqlQuery(con, paste("SELECT * FROM CSVDOWNLOAD_A_RW_VASCULAR_PLANT"))
lookup <- sqlQuery(con, paste("SELECT * FROM RAWDATA_G_OFFGRID_SITE_LABEL"))
taxo <- sqlQuery(con, paste("SELECT * FROM PUBLIC_ACCESS_PUBLIC_DETAIL_TAXONOMYS"))
## close connection
close(con)

## Labels etc
res$OnOffGrid <- ifelse(res$SITE > 1656, "OGW", "RT")
res$OGLabel <- lookup$SITE_LABEL[match(res$SITE, lookup$SITE_ID)]
res$SiteLabel <- toupper(ifelse(res$SITE > 1656, as.character(res$OGLabel), as.character(res$SITE)))
tmp <- lapply(res$SiteLabel, function(z) strsplit(z, split="-")[[1]])
res$DataProvider <- sapply(tmp, function(z) ifelse(length(z) == 1, "ABMI", z[2]))
#res$Visit <- sapply(tmp, function(z) ifelse(length(z) == 1, 1L, as.integer(z[4])))
res$Visit <- 1
res$ClosestABMISite <- sapply(1:length(tmp), function(i) ifelse(length(tmp[[i]]) == 1, res$SITE[i], as.integer(tmp[[i]][3])))
res$Label <- with(res, paste(OnOffGrid, DataProvider, SiteLabel, YEAR, Visit, "WZ", ZONE1, sep="_"))
res$Label2 <- with(res, paste(OnOffGrid, DataProvider, SiteLabel, YEAR, Visit, sep="_"))
rm(tmp)

## exclude parts here if necessary
## remove bad wetlands
source("y:/Oracle_access/src2/00globalvars_wetland.R")
keep <- !(res$SITE %in% REJECT)
res <- res[keep,]

## crosstab
xt <- Xtab(~ Label + SCIENTIFICNAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"), 
    subset = !(res$SCIENTIFICNAME %in% c("VNA", "DNC")), drop.unused.levels = TRUE)
## get taxonomy
z <- nonDuplicated(res[!(res$SCIENTIFICNAME %in% c("VNA", "DNC")),],
    res$SCIENTIFICNAME[!(res$SCIENTIFICNAME %in% c("VNA", "DNC"))], TRUE)
## add here higher taxa too
z <- z[,c("TSNID","COMMONNAME","SCIENTIFICNAME","TAXONOMICRESOLUTION")]
z2 <- taxo[taxo$SCIENTIFIC_NAME %in% z$SCIENTIFICNAME,]
z <- data.frame(z, z2[match(z$SCIENTIFICNAME, z2$SCIENTIFIC_NAME),setdiff(colnames(taxo), colnames(z))])
#z[] <- lapply(z, function(z) z[drop=TRUE])
#summary(z)

x <- nonDuplicated(res[,c("Label", "Label2", "ROTATION", "SITE", "YEAR", "FIELDDATE", "CREWMEMBER", 
    "ZONE1", "OnOffGrid", 
    "OGLabel", "SiteLabel", "DataProvider", "Visit", "ClosestABMISite")], res$Label, TRUE)

## crosstab on PC level
m <- Mefa(xt, x, z)
## exclude not species level taxa
table(m@taxa$RANK_NAME) # here are sub-specific levels
m <- m[,m@taxa$RANK_NAME %in% c("Genus", "Species", "Subspecies", "Variety")]
## exclude not species level taxa
#m <- m[,m@taxa$TAXONOMICRESOLUTION == "Species"]
xtab(m) <- as(xtab(m) > 0, "dgCMatrix")

## site level info
m2 <- groupSums(m, 1, m@samp$Label2)
xtab(m2) <- as(xtab(m2) > 0, "dgCMatrix")

## crosstabs
res1 <- as.matrix(m)
res2 <- as.matrix(m2)
tax <- taxa(m)
tax[] <- lapply(tax, function(z) z[drop=TRUE])
str(res1)
str(res2)
m
m2
str(tax)

## write static files into csv
write.csv(res1, file=paste(D, "OUT_", T, "_Species_WZ", d, ".csv",sep=""))
write.csv(res2, file=paste(D, "OUT_", T, "_Species_Site", d, ".csv",sep=""))
write.csv(tax, file=paste(D, "OUT_", T, "_Species_Taxa", d, ".csv",sep=""))

if (do.prof)
    summaryRprof(proffile)
if (do.image)
    save.image(paste(D, "OUT_", tolower(T), d, ".Rdata",sep=""))
## quit without saving workspace
quit(save="no")

sss <- read.csv("y:/Oracle_access/src2/sitesummary_download.csv")[,c(1,2,3,5)]
setdiff(samp(m)$SiteLabel[samp(m)$YEAR==2011],sss[sss$Wetland=="Completed","SiteLabel"])
setdiff(sss[sss$Wetland=="Completed","SiteLabel"],samp(m)$SiteLabel[samp(m)$YEAR==2011])

setdiff(setdiff(sss[sss$Wetland=="Completed","SiteLabel"],samp(m)$SiteLabel[samp(m)$YEAR==2011]), REJECT)

"OGW-ABMI-1058-2" is missing from summaries but expected based on SSWB

