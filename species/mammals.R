#from:	 QiongYan Fang <qfang@ualberta.ca>
#to:	 Shawn Morrison <shawn.morrison@dryasresearch.com>
#cc:	 Peter Solymos <solymos@ualberta.ca>
#date:	 Wed, Jun 12, 2013 at 10:54 AM
#subject:	 RE: updated 2002-09 Mammal Transect Segment numbering
#Hi Shawn,
#I have two tables updated:
#EXT_ILM_MAMMAL (species table)
#EXT_ILM_MAMMAL_HU (habitat table)
#These two following tables are unchanged:
#EXT_ILM_MAMMAL_COORDINATES (has segment information, but it is already 1-40.  I don’t have coordinates 
#information to update it. I can renumber it, as the segment already within ranges.)
#EXT_ILM_MAMMAL_SNOW (doesn’t contain any segment information)
#These two following tables are to be removed, but I leave them in, in case you want to compare them.
#EXT_ILM_MAMMAL_HU_OLD
#EXT_ILM_MAMMAL_OLD

ROOT <- "y:/Oracle_access_2015"
getwd()
if (interactive())
    source("~/repos/abmianalytics/species/00globalvars.R") else source("00globalvars.R")

T <- "Mammals"
if (do.prof) {
    proffile <- paste(D, "/OUT_", tolower(T), d, ".Rprof",sep="")
    Rprof(proffile)
}

if (FALSE) {
## load package RODBC
stopifnot(require(RODBC))
## extablish connection
con <- odbcConnectAccess2007(DBVERSION)

## ABMI data
abmi1 <- sqlFetch(con, "CSVDOWNLOAD_A_RT_SNOWTRACK_TRANSECT") # species table
abmi2 <- sqlFetch(con, "CSVDOWNLOAD_A_RT_SNOWTRACKING_HU") # habitat table
abmi3 <- sqlFetch(con, "CSVDOWNLOAD_A_RT_SNOWTRACKING_SITE_V") # snow etc
#abmi4 <- read.csv("y:/Oracle_access/mamTransect.csv")

## ILM data
ilm1 <- sqlFetch(con, "PUBLIC_ACCESS_EXT_ILM_MAMMAL") # species table
ilm2 <- sqlFetch(con, "PUBLIC_ACCESS_EXT_ILM_MAMMAL_HU") # habitat table
ilm3 <- sqlFetch(con, "PUBLIC_ACCESS_EXT_ILM_MAMMAL_SNOW") # snow etc
ilm4 <- sqlFetch(con, "PUBLIC_ACCESS_EXT_ILM_MAMMAL_COORDINATES") # xy

close(con)
}

## ABMI data
abmi1 <- read.csv(file.path(ROOT, "data/snowtracking-transect.csv"))
abmi2 <- read.csv(file.path(ROOT, "data/snowtracking-hu.csv"))
abmi3 <- read.csv(file.path(ROOT, "data/snowtracking-site-v.csv"))

## ILM data
ilm1 <- read.csv(file.path(ROOT, "data/snowtracking-ilm.csv"))
ilm2 <- read.csv(file.path(ROOT, "data/snowtracking-ilm-hu.csv"))
ilm3 <- read.csv(file.path(ROOT, "data/snowtracking-ilm-snow.csv"))
ilm4 <- read.csv(file.path(ROOT, "data/snowtracking-ilm-coord.csv"))


abmi1$Date <- strptime(abmi1$ADATE, "%d-%b-%y") # 27-Jan-05
abmi1$YearSampling <- abmi1$Date$year + 1900
table(ABMI=abmi1$YEAR, Date=abmi1$YearSampling)
stopifnot(all(abmi1$YEAR <= abmi1$YearSampling))

table(ABMI=ilm1$ABMIYEAR, Date=ilm1$SAMPLEYEAR)
stopifnot(all(ilm1$ABMIYEAR <= ilm1$SAMPLEYEAR))

## labels

## transects with 2 visits
tmp <- abmi3[abmi3$SITE_LABEL %in% as.character(abmi3$SITE_LABEL)[abmi3$TMTC_VISIT > 1],]
with(tmp, table(as.character(SITE_LABEL), YEAR))
vis2 <- as.character(tmp$SITE_LABEL)

abmi1$INTER_SEG <- ifelse(abmi1$YEAR < 2005, abmi1$TMTT_SEGMENT, 
    ((abmi1$TMTT_SEGMENT - 1) %/% 4) + 1)
table(abmi1$TMTT_SEGMENT, abmi1$INTER_SEG)
abmi1$label_tr <- with(abmi1, paste("T", "IG", "ABMI", SITE, YEAR, 1, sep="_"))
abmi1$label_int <- with(abmi1, paste("T", "IG", "ABMI", SITE, YEAR, 1, "SI", INTER_SEG, sep="_"))
#tmp <- abmi1[abmi1$label_tr %in% vis2,]
#with(tmp, table(as.character(label_tr), INTER_SEG))
#with(abmi1, table(as.character(label_tr), INTER_SEG))
levels(abmi1$TSNID)[levels(abmi1$TSNID) == "VNA"] <- NA
abmi1$TSNID <- as.integer(as.character(abmi1$TSNID))

#ilm1$method <- as.factor(ifelse(ilm1$ABMIYEAR < 2005, "triangle", "linear"))
table(ilm1$SEGMENT, ilm1$ABMIYEAR)
ilm1$INTER_SEG <- ifelse(ilm1$ABMIYEAR < 2005, ilm1$SEGMENT, 
    ((ilm1$SEGMENT - 1) %/% 4) + 1)
table(ilm1$INTER_SEG, ilm1$ABMIYEAR)
ilm1$label_tr <- with(ilm1, paste("T", "OG", "ILM", ABMI_SITE, 
    ABMIYEAR, 1, sep="_"))
ilm1$label_int <- with(ilm1, paste("T", "OG", "ILM", ABMI_SITE, 
    ABMIYEAR, 1, "SI", INTER_SEG, sep="_"))
ilm1$length <- ifelse(ilm1$ABMIYEAR < 2005, 1000, 250)
ilm1 <- ilm1[ilm1$SPP_COUNT != "DNC",]

ilm3 <- ilm3[ilm3$SOURCES=="ILM",]
ilm3[] <- lapply(ilm3, function(z) z[drop=TRUE])
ilm3$label_tr <- with(ilm3, paste("T", "OG", "ILM", SITE, ABMIYEAR, 1, sep="_"))

res1 <- abmi1[,c("label_tr", "label_int", "COMMON_NAME", "SCIENTIFIC_NAME", "TSNID")]
res2 <- ilm1[,c("label_tr", "label_int", "COMMON_NAME","SCIENTIFIC_NAME", "TSN")]
colnames(res2) <- colnames(res1)
res1$TSNuse <- res1$TSNID
res2$TSNuse <- res2$TSNID

cid <- intersect(res1$TSNID, res2$TSNID)
cid <- cid[!is.na(cid)]
tmp1 <- res1[res1$TSNID %in% cid,3:5]
tmp1 <- nonDuplicated(tmp1, TSNID, TRUE)[,-3]
tmp2 <- res2[res2$TSNID %in% cid,3:5]
tmp2 <- nonDuplicated(tmp2, TSNID, TRUE)[,-3]
cbind(tmp1[as.character(cid),], tmp2[as.character(cid),])

setdiff(res1$TSNID, res2$TSNID)
tmp3 <- res1[res1$TSNID %in% setdiff(res1$TSNID, res2$TSNID),3:5]
(tmp3 <- nonDuplicated(tmp3, TSNID, TRUE))

setdiff(res2$TSNID, res1$TSNID)
tmp4 <- res2[res2$TSNID %in% setdiff(res2$TSNID, res1$TSNID),3:5]
(tmp4 <- nonDuplicated(tmp4, TSNID, TRUE))

res1$TSNuse[res1$TSNID %in% c(180604, 180607)] <- 180603 # foxes
res1$TSNuse[res1$TSNID %in% c(677542)] <- 175803  # lagopus

res <- rbind(res1, res2)
res$TSNuse[is.na(res$TSNuse)] <- -1
res$TSNuse[res$SCIENTIFIC_NAME == "DNC"] <- NA

tax0 <- nonDuplicated(res[,3:6], TSNID, FALSE)
tax02 <- nonDuplicated(res[,3:6], COMMON_NAME, FALSE)
tax02 <- tax02[is.na(tax02$TSNuse),]
tax0 <- rbind(tax0, tax02) 
tax0$COMMON_NAME <- as.character(tax0$COMMON_NAME)
tax0$COMMON_NAME[tax0$SCIENTIFIC_NAME == "Erethizon dorsatum"] <- "Porcupine"
tax0$COMMON_NAME[tax0$SCIENTIFIC_NAME == "Vulpes velox"] <- "Foxes"
tax0$COMMON_NAME[tax0$SCIENTIFIC_NAME == "Vulpes vulpes"] <- "Foxes"
tax0$COMMON_NAME[tax0$SCIENTIFIC_NAME == "Rangifer tarandus caribou"] <- "Caribou"

## dds
sam1 <- abmi3[,c("SITE_LABEL", "TMTC_DSS", "YEAR")]#, "TMTC_SNOW_CONDITION", "TMTC_TEMPERATURE")]
sam2 <- ilm3[,c("label_tr", "DSS", "ABMIYEAR")]#, "SNOW_CONDITION", "TEMP")]
colnames(sam1) <- colnames(sam2)
levels(sam2$DSS)[levels(sam2$DSS) == "3 to 6"] <- "4"
sam2$DSS <- as.integer(as.character(sam2$DSS))
sam <- rbind(sam1, sam2)
sam <- nonDuplicated(sam, label_tr, TRUE)
#sam$Length_assumed <- ifelse(sam$ABMIYEAR < 2005, 9000, 10000)

len1 <- abmi1[,c("label_tr", "TMTC_TOTAL_TRANSECT_LEN", "YEAR")]
len1 <- nonDuplicated(len1, label_tr, TRUE)
levels(len1$TMTC_TOTAL_TRANSECT_LEN)[levels(len1$TMTC_TOTAL_TRANSECT_LEN) == "VNA"] <- NA
len1$TMTC_TOTAL_TRANSECT_LEN <- as.integer(as.character(len1$TMTC_TOTAL_TRANSECT_LEN))
sam$Length_measured <- len1$TMTC_TOTAL_TRANSECT_LEN[match(rownames(sam), rownames(len1))]
#plot(sam$Length_measured, sam$Length_assumed);abline(0,1)

## no HUMAN and domestic here (inconsistent across data sets it seems)
xt <- Xtab(~ label_int + TSNuse, res, drop.unused.levels=FALSE,
    subset = !is.na(res$TSNuse) & res$TSNuse != (-96000), cdrop="-1")
xt[xt>0] <- 1

x <- nonDuplicated(res[,1:2], label_int, TRUE)
tmp <- data.frame(t(sapply(rownames(x), function(z) strsplit(z, "_")[[1]])))
colnames(tmp) <- c("Protocol", "OnOffGrid", "DataProvider", "Site", "Year", "Visit", "SubType", "InterSegID")
x <- data.frame(x, tmp[rownames(x),])
z <- nonDuplicated(res[!is.na(res$TSNuse),c(6,3,4)], TSNuse, TRUE)
m <- Mefa(xt, x, z)
m
colnames(m) <- tax0$COMMON_NAME[match(colnames(m), as.character(tax0$TSNuse))]
colnames(m) <- nameAlnum(colnames(m),"first",collapse="")


m2 <- groupSums(m, 1, samp(m)$label_tr, ra.rm=TRUE)
x2 <- nonDuplicated(samp(m), label_tr, TRUE)
x2 <- data.frame(x2, sam[rownames(x2),-(1)])
x2$label_int <- NULL
tmp <- table(samp(m)$label_tr)
x2$n_inter <- tmp[match(rownames(x2), names(tmp))]
x2$Length_assumed <- x2$n_inter * 1000
#plot(x2$Length_assumed, x2$Length_measured)

samp(m2) <- x2

resm1 <- as.matrix(m)
resm2 <- as.matrix(m2)
tax <- NULL
samp1 <- samp(m)
samp2 <- samp(m2)
str(resm1)
str(samp1)
str(samp2)
m
m2
str(tax)

stopifnot(all(apply(resm2, 1, max) <= x2$n_inter))

if (FALSE) {
aaa <- read.csv("y:/Oracle_access/mammals/Segments with no mammal data.csv")
aaa$onoff <- ifelse(substr(aaa$Site, 1, 3) == "ILM", "OG", "IG")
aaa$dprov <- ifelse(substr(aaa$Site, 1, 3) == "ILM", "ILM", "ABMI")
aaa$Lab <- with(aaa, paste("T", onoff, dprov, Site, Year, "1", "SI", Seg, sep="_"))
bbb <- rownames(resm1)[rowSums(resm1) == 0]
intersect(aaa$Lab, bbb)
zzz2 <- data.frame(InterSegment=setdiff(aaa$Lab, bbb), Type="InGIS_NotInSpp")
zzz3 <- data.frame(InterSegment=setdiff(bbb, aaa$Lab), Type="NotInGIS_InSpp")
zzz <- rbind(zzz2, zzz3)
write.csv(zzz,file="y:/Oracle_access/mammals/ProblemSegments_Mammals_Oct18-2013.csv")

}

## write static files into csv
if (combine.tables) {
    out1 <- data.frame(samp1, resm1)
    out2 <- data.frame(samp2, resm2)
    write.csv(out1, file=paste(D, "/OUT_", T, "_Species_InterSegment", d, 
        ".csv",sep=""), row.names = FALSE)
    write.csv(out2, file=paste(D, "/OUT_", T, 
        "_Species_Transect-Binomial-Length-DSS", d, ".csv",sep=""), 
        row.names = FALSE)
} else {
    write.csv(resm1, file=paste(D, "/OUT_", T, "_Species_InterSegment", 
        d, ".csv",sep=""), row.names = TRUE)
    write.csv(resm2, file=paste(D, "/OUT_", T, "_Species_Transect-Binomial", 
        d, ".csv",sep=""), row.names = TRUE)
    write.csv(samp1, file=paste(D, "/OUT_", T, "_Species_SegmentInfo", 
        d, ".csv",sep=""), row.names = TRUE)
    write.csv(samp2, file=paste(D, "/OUT_", T, "_Species_TransectInfo-Length-DSS", 
        d, ".csv",sep=""), row.names = TRUE)
}
tax <- data.frame(Species=colnames(m), taxa(m))
write.csv(tax, file=paste(D, "/OUT_", T, "_Species_Taxa", d, ".csv",sep=""), 
    row.names = FALSE)

if (do.prof)
    summaryRprof(proffile)
if (do.image)
    save.image(paste(D, "/OUT_", tolower(T), d, ".Rdata",sep=""))
## quit without saving workspace
quit(save="no")
