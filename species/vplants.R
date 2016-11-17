ROOT <- "e:/peter/AB_data_v2016/oracle"
OUTDIR <- "e:/peter/AB_data_v2016/data/species"
getwd()
if (interactive())
    source("~/repos/abmianalytics/species/00globalvars.R") else source("00globalvars.R")

T <- "VPlants"

com_dom <- FALSE

if (do.prof) {
    proffile <- paste(D, "OUT_", tolower(T), d, ".Rprof",sep="")
    Rprof(proffile)
}

if (FALSE) {## extablish connection
con <- odbcConnectAccess2007(DBVERSION)
## queries
res <- sqlQuery(con, paste("SELECT * FROM CSVDOWNLOAD_A_RT_VASCULAR_PLANT_V"))
#lookup <- sqlQuery(con, paste("SELECT * FROM RAWDATA_G_OFFGRID_SITE_LABEL"))
taxo <- sqlQuery(con, paste("SELECT * FROM PUBLIC_ACCESS_PUBLIC_DETAIL_TAXONOMYS"))
## close connection
close(con)
}

#res <- read.csv(file.path(ROOT, "data", "plants2015",
#    "A_T15_Vascular_Plants_AllYear_AB_689250437025811679.csv"), row.names=NULL)
res <- read.csv(file.path(ROOT, "vplants.csv"))
#gis <- read.csv(file.path(ROOT, "data/sitemetadata.csv"))
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
taxo <- read.csv(file.path(ROOT, "taxonomy.csv"))

## Labels etc
rr <- LabelFun(res)
rr <- nonDuplicated(rr, Label, TRUE)

res <- data.frame(res, rr[match(res$SITE_LABEL, rr$Label),])

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

## using species only
levels(res$SCIENTIFIC_NAME) <- gsub("X ", "", levels(res$SCIENTIFIC_NAME))
levels(res$SCIENTIFIC_NAME) <- gsub(" x ", " ", levels(res$SCIENTIFIC_NAME))
levels(res$SCIENTIFIC_NAME) <- sapply(strsplit(levels(res$SCIENTIFIC_NAME), " "), function(z) {
    paste(z[1:min(2, length(z))], collapse=" ")
})

res$SPECIES_OLD <- res$SCIENTIFIC_NAME
levels(res$SCIENTIFIC_NAME) <- nameAlnum(levels(res$SCIENTIFIC_NAME), capitalize="mixed", collapse="")
res$SCIENTIFIC_NAME <- droplevels(res$SCIENTIFIC_NAME)

xt <- Xtab(~ Label + SCIENTIFIC_NAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"),
    rdrop=qs.to.exclude, drop.unused.levels = FALSE)
xt[xt>0] <- 1

if (com_dom) {
    res$comdom <- res$COMM_DOM
    levels(res$comdom)[!(levels(res$comdom) %in% c("Dominant", "Common"))] <- "Uncommon"
    table(res$comdom)
    xt0 <- Xtab(~ Label + SCIENTIFIC_NAME + comdom, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"),
        rdrop=pcs.to.exclude, drop.unused.levels = FALSE)
    xt1 <- xt0$Uncommon
    xt1[xt1>0] <- 1
    xt10 <- xt0$Common * 10
    xt10[xt10>0] <- 10
    xt100 <- xt0$Dominant * 100
    xt100[xt100>0] <- 100

    xt <- xt1
    xt[xt10>0] <- xt10[xt10>0]
    xt[xt100>0] <- xt100[xt100>0]
}
#xt <- Xtab(~ Label + SCIENTIFIC_NAME, res, cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"),
#    subset = !(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")), drop.unused.levels = TRUE)
## get taxonomy
z <- nonDuplicated(res[!(res$SCIENTIFIC_NAME %in% c("VNA", "DNC")),],
    res$SCIENTIFIC_NAME[!(res$SCIENTIFIC_NAME %in% c("VNA", "DNC"))], TRUE)
## add here higher taxa too
z <- z[,c("TSN_ID","COMMON_NAME","SCIENTIFIC_NAME","RANK_NAME","SPECIES_OLD")]
z2 <- taxo[taxo$SCIENTIFIC_NAME %in% z$SCIENTIFIC_NAME,]
z <- data.frame(z, z2[match(z$SCIENTIFIC_NAME, z2$SCIENTIFIC_NAME),setdiff(colnames(taxo), colnames(z))])
levels(z$RANK_NAME)[levels(z$RANK_NAME) %in% c("Subspecies","Variety")] <- "Species"
#z[] <- lapply(z, function(z) z[drop=TRUE])
#summary(z)
## trim trailing white space -- just to make sure
#for (i in 1:ncol(z))
#    if (is.factor(z[[i]]))
#        levels(z[[i]]) <- sub(' +$', '', levels(z[[i]]))

x <- nonDuplicated(res[,c("Label", "Label2", "ROTATION", "SITE", "YEAR", "GSA_DATE", "GCL_NAMES",
    "SubType", "SubTypeID", "OnOffGrid",
    "SiteLabel", "DataProvider", "Visit", "ClosestABMISite")], res$Label, TRUE)

## crosstab on PC level
m <- Mefa(xt, x, z)
## exclude not species level taxa
#table(m@taxa$RANK_NAME) # here are sub-specific levels
m <- m[,m@taxa$RANK_NAME %in% c("Genus", "Species")]

## here I want to collapse sub-specific taxa, but TAX WB is not that clean to map these
## this is already taken care of by truncating species names
if (FALSE) {
tmp <- m@taxa$SCIENTIFIC_NAME_OLD
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

#xtab(m) <- as(xtab(m) > 0, "dgCMatrix") # can be >1 due to COMM_DOM -- handled above
table(m@taxa$RANK_NAME) # here are sub-specific levels

## site level info
mmm <- m[samp(m)$SubTypeID != "DNC",]
m2 <- groupSums(mmm, 1, mmm@samp$Label2)
no_q <- table(m[samp(m)$SubTypeID != "DNC",]@samp$Label2)
no_q <- no_q[rownames(m2)]
samp(m2) <- data.frame(NoQuadrantsPerSite = no_q)
#xtab(m2) <- as(xtab(m2) > 0, "dgCMatrix") # it might be 5 due to QU=DNC

range(xtab(m))
range(xtab(m2))

## crosstabs
if (!combine.tables) {
    res1 <- as.matrix(m)
    res2 <- as.matrix(m2)
    rntw <- TRUE
} else {
    res1 <- data.frame(samp(m), as.matrix(m))
    mmm <- Mefa(xtab(m2), data.frame(nonDuplicated(samp(m)[,-which(colnames(samp(m))=="Label")],
        Label2, TRUE)[rownames(m2),], samp(m2)))
    res2 <- data.frame(samp(mmm), as.matrix(mmm))
    rntw <- FALSE
}
tax <- taxa(m)
tax[] <- lapply(tax, function(z) z[drop=TRUE])
str(res1)
str(res2)
m
m2
str(tax)
range(xtab(m))
range(xtab(m2))

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
            PROT <- c("T15 Vascular Plants")
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
write.csv(res1, file=paste(D, "/OUT_", T, "_Species_QU-PA", d,
    ifelse(com_dom, "_ComDom", ""), ".csv",sep=""), row.names = rntw)
write.csv(res2, file=paste(D, "/OUT_", T, "_Species_Site-Binomial", d,
    ifelse(com_dom, "_ComDom", ""), ".csv",sep=""), row.names = rntw)
write.csv(tax, file=paste(D, "/OUT_", T, "_Species_Taxa", d,
    ifelse(com_dom, "_ComDom", ""), ".csv",sep=""), row.names = TRUE)
#write.csv(rr, file=paste(D, "/OUT_", T, "_SampleInfo", d,
#    ifelse(com_dom, "_ComDom", ""), ".csv",sep=""))
write.csv(data.frame(Excluded_Sites=sort(pcs.to.exclude)),
    file=paste(D, "/OUT_", T, "_Excluded_Sites", d,
    ifelse(com_dom, "_ComDom", ""), ".csv",sep=""), row.names = FALSE)


if (do.prof)
    summaryRprof(proffile)
if (do.image)
    save.image(paste(D, "/OUT_", tolower(T), d,
        ifelse(com_dom, "_ComDom", ""), ".Rdata",sep=""))
## quit without saving workspace
quit(save="no")

