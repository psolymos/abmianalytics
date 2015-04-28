ROOT <- "y:/Oracle_access_2015"
getwd()
if (interactive())
    source("~/repos/abmianalytics/species/R/00globalvars.R") else source("R/00globalvars.R")

T <- "Lichens"
if (do.prof) {
    proffile <- paste(D, "/OUT_", tolower(T), d, ".Rprof",sep="")
    Rprof(proffile)
}

if (FALSE) {## extablish connection
con <- odbcConnectAccess2007(DBVERSION)
## queries
resa0 <- sqlFetch(con, "CSVDOWNLOAD_A_RT_LAB_LICHEN_03_08_V")
resb0 <- sqlFetch(con, "CSVDOWNLOAD_A_RT_LAB_LICHEN_CURRENT_V")
#resb0 <- sqlFetch(con, "PUBLIC_ACCESS_MOSS_LICHEN_DATA_V")
#lookup <- sqlQuery(con, paste("SELECT * FROM RAWDATA_G_OFFGRID_SITE_LABEL"))
taxo <- sqlQuery(con, paste("SELECT * FROM PUBLIC_ACCESS_PUBLIC_DETAIL_TAXONOMYS"))
## close connection
close(con)
}

resa0 <- read.csv(file.path(ROOT, "data/lichens0308.csv"))
resb0 <- read.csv(file.path(ROOT, "data/lichens09.csv"))
taxo <- read.csv(file.path(ROOT, "data/taxonomy.csv"))
cs <- read.csv(file.path(ROOT, "data/moss-lichen-collstatus.csv"))

## Labels etc
sum(grepl("ALPAC-SK", resa0$SITE_LABEL))
resa0 <- resa0[!grepl("ALPAC-SK", resa0$SITE_LABEL),]
rra <- LabelFun(resa0)
rra <- nonDuplicated(rra, Label, TRUE)
resa <- data.frame(resa0, rra[match(resa0$SITE_LABEL, rra$Label),])
#resa <- data.frame(resa0, rra[match(resa0$SITE_LABEL, rownames(rra)),])

sum(grepl("ALPAC-SK", resb0$SITE_LABEL))
resb0 <- resb0[!grepl("ALPAC-SK", resb0$SITE_LABEL),]
rrb <- LabelFun(resb0)
rrb <- nonDuplicated(rrb, Label, TRUE)
resb <- data.frame(resb0, rrb[match(resb0$SITE_LABEL, rrb$Label),])

## identification issues to exclude
resb <- resb[resb$QUALIFIER != "cf.",]

if (FALSE) {
tmp <- t(sapply(as.character(resb0$MS_SITE), function(z) {
    zz <- strsplit(z, "-")[[1]]
    if (length(zz) == 4)
        zz else c("IG", "ABMI", zz, NA)
}))
resb0$SITE_LABEL <- as.factor(with(resb0, paste0("T_", tmp[,1], "_", tmp[,2],
    "_", MS_SITE, "_", MS_YEAR, "_1_PL-MH_", MS_QUADRANT, "-XX")))
sum(grepl("ALPAC-SK", resb0$SITE_LABEL))
resb0 <- resb0[!grepl("ALPAC-SK", resb0$SITE_LABEL),]
rrb <- LabelFun(resb0)
rrb <- nonDuplicated(rrb, Label, TRUE)
resb <- data.frame(resb0, rrb[match(resb0$SITE_LABEL, rrb$Label),])
#resb <- data.frame(resb0, rrb[match(resb0$SITE_LABEL, rownames(rrb)),])
}

sum(grepl("ALPAC-SK", resa0$SITE_LABEL))
sum(grepl("ALPAC-SK", resb0$SITE_LABEL))
sum(grepl("ALPAC-SK", resa$Label2))
sum(grepl("ALPAC-SK", resb$Label2))
table(resa0$YEAR)
table(resb0$YEAR)

#resa$TMLIR_PLOT <- resa$TLML_PLOT
#resa$TMLIR_STRATUM <- resa$TLML_MICROHABITAT_TYPE
#resb$TSN_ID <- resb$TSNID
#resb$SCIENTIFIC_NAME <- resb$MLS_SCIENTIFIC_NAME
#resb$YEAR <- resb$MS_YEAR

resa$TMLIR_COLLECTION_STATUS <- factor("C", levels(resb$TMLIR_COLLECTION_STATUS))
#ccol <- c("SITE_LABEL","ROTATION", "SITE", "YEAR", "ADATE", "CREWS", "TMLIR_PLOT", 
#    "TMLIR_STRATUM", "SCIENTIFIC_NAME", "COMMON_NAME", "RANK_NAME", "TSN_ID",
#    "OnOffGrid", "SiteLabel", "DataProvider", 
#    "SubType", "SubTypeID", "OGSeqenceID",
#    "Visit", "ClosestABMISite", "Label", "Label2")
ccol <- c("SITE_LABEL","YEAR", "SCIENTIFIC_NAME", 
    "OnOffGrid", "SiteLabel", "DataProvider", 
    "SubType", "SubTypeID", "OGSeqenceID",
    "Visit", "ClosestABMISite", "Label", "Label2","TMLIR_COLLECTION_STATUS")
res <- rbind(resa[,ccol], resb[,ccol])
#res <- rbind(resa, resb)

sum(grepl("ALPAC-SK", res$SITE_LABEL))

## exclude parts here if necessary
#keep <- 1:nrow(res)
#res <- res[keep,]

tmp <- as.factor(substr(res$SubTypeID, 1, 2))
levels(tmp)[levels(tmp)=="DN"] <- "DNC"
levels(tmp)[levels(tmp)=="VN"] <- "VNA"
levels(tmp)[!(levels(tmp) %in% c("DNC", "VNA", "NW","NE","SW","SE"))] <- "Plot1ha"
tmp[tmp != "Plot1ha" & res$YEAR < 2009] <- "Plot1ha"
table(tmp, res$YEAR)
res$Quadrant <- tmp
res$Label <- paste(res$Label2, "PL", res$Quadrant, sep="_")

## this should all apply at site level!
## unique labels where species is DNC
tmp001 <- unique(res[res$SCIENTIFIC_NAME %in% c("DNC"),c("Label2")])
## unique labels where species is NOT DNC
tmp002 <- unique(res[!(res$SCIENTIFIC_NAME %in% c("DNC")),c("Label2")])
## use the collection status (does not work for 03-08 data)
tmp003 <- unique(res[res$TMLIR_COLLECTION_STATUS %in% c("DNC","PNA"),c("Label2")])
## unique labels where only DNC species occur (no data collected)
pcs.to.exclude <- setdiff(tmp001, tmp002)
## unique labels where collection status indicates so
pcs.to.exclude <- union(pcs.to.exclude, tmp003)
qs.to.exclude <- unique(as.character(res$Label[res$Label2 %in% pcs.to.exclude]))

## crosstab
res$sppnam <- res$SCIENTIFIC_NAME
levels(res$sppnam) <- nameAlnum(levels(res$sppnam),"first",collapse="")
res$sppnam <- droplevels(res$sppnam)

xt <- Xtab(~ Label + sppnam, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"), 
    rdrop=qs.to.exclude, drop.unused.levels = FALSE)
## get rid of them here, because drop=FALSE
sum(grepl("ALPAC-SK", rownames(xt)))
xt <- xt[!grepl("ALPAC-SK", rownames(xt)),]
sum(grepl("ALPAC-SK", rownames(xt)))

## get taxonomy
z <- nonDuplicated(res[!(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")) & !is.na(res$SCIENTIFIC_NAME),],
    res$sppnam[!(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")) & !is.na(res$SCIENTIFIC_NAME)], TRUE)
if (FALSE) {
## add here higher taxa too
z <- z[,c("TSN_ID","COMMON_NAME","SCIENTIFIC_NAME","RANK_NAME")]
z2 <- taxo[taxo$SCIENTIFIC_NAME %in% z$SCIENTIFIC_NAME,]
z <- data.frame(z, z2[match(z$SCIENTIFIC_NAME, z2$SCIENTIFIC_NAME),setdiff(colnames(taxo), colnames(z))])
#z[] <- lapply(z, function(z) z[drop=TRUE])
#summary(z)
## trim trailing white space -- just to make sure
for (i in 1:ncol(z))
    if (is.factor(z[[i]]))
        levels(z[[i]]) <- sub(' +$', '', levels(z[[i]]))
}

#x <- nonDuplicated(res[,c("Label", "Label2", "ROTATION", "SITE", "YEAR", "ADATE", "CREWS", 
#    "OnOffGrid", 
#    "SiteLabel", "DataProvider", "Visit", "ClosestABMISite",
#    "SubType", "SubTypeID", "OGSeqenceID", "Quadrant")], res$Label, TRUE)
x <- nonDuplicated(res[,c("Label", "Label2", "YEAR", 
    "OnOffGrid", 
    "SiteLabel", "DataProvider", "Visit", "ClosestABMISite",
    "SubType", "SubTypeID", "OGSeqenceID", "Quadrant")], res$Label, TRUE)

## crosstab on PC level
m <- Mefa(xt, x, z)
## exclude not species level taxa
#table(m@taxa$RANK_NAME) # here are sub-specific levels
#m <- m[,taxa(m)$RANK_NAME %in% c("Genus", "Species", "Subspecies", "Variety")]

## here I want to collapse sub-specific taxa, but TAX WB is not that clean to map these
if (FALSE) {
tmp <- m@taxa$SCIENTIFIC_NAME
tmp <- sapply(strsplit(as.character(tmp)," "), function(z) {
    if (z[1]=="X") {
        dd <- if (length(z) == 2)
            paste(z[2], "spp.") else paste(z[2], z[3])
        return(dd)
    }
    if (length(z) == 1)
        paste(z[1], "spp.") else paste(z[1], z[2])
    })
taxa(m)$SpeciesFinal <- as.factor(tmp)
zz <- nonDuplicated(m@taxa, m@taxa$SpeciesFinal, TRUE)
m <- groupSums(m, 2, m@taxa$SpeciesFinal, replace=zz)
}
xtab(m) <- as(xtab(m) > 0, "dgCMatrix") # can be >1 due to COMM_DOM
#table(m@taxa$RANK_NAME) # here are sub-specific levels

## site level info
mmm <- m[samp(m)$Quadrant != "DNC",]
m2 <- groupSums(mmm, 1, mmm@samp$Label2)
samp(m2) <- nonDuplicated(samp(mmm), Label2, TRUE)
#xtab(m2) <- as(xtab(m2) > 0, "dgCMatrix") # it might be 5 due to QU=DNC

## this is 0/1 at PC level
m01 <- m
xtab(m01)[xtab(m01)>0] <- 1
m012 <- groupSums(m01, 1, m@samp$Label2)
samp(m012) <- nonDuplicated(samp(m01), Label2, TRUE)

aa=apply(xtab(m012),1,max)
table(aa,samp(m012)$YEAR)

samp(m012)[samp(m012)$YEAR < 2009 & aa > 1,]
xtab(m012)[samp(m012)$YEAR < 2009 & aa > 1,]

xtab(m2)[xtab(m2)>0] <- 1

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
            PROT <- c("T17 Moss and Lichen Microhabitat", "T18 Moss Lichen Search", "T19 Moss Identification")
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
write.csv(res1, file=paste(D, "/OUT_", T, "_Species_QU-PA", d, ".csv",sep=""), row.names = rntw)
write.csv(res2, file=paste(D, "/OUT_", T, "_Species_Site-PA", d, ".csv",sep=""), row.names = rntw)
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
