library(mefa4)

ROOT <- "e:/peter/AB_data_v2016/data/aru-raw"

issu <- read.csv(file.path(ROOT, "2015-ABMI-MC-SC-AudioRecordingTranscription-issues.csv"))
depl <- read.csv(file.path(ROOT, "2015-ABMI-MC-SC-AudioRecordingTranscription-depl.csv"))
spp <- read.csv(file.path(ROOT, "2015-ABMI-MC-SC-AudioRecordingTranscription-spp.csv"))
det <- read.csv(file.path(ROOT, "2015-ABMI-MC-SC-AudioRecordingTranscription-det.csv"))

## resolve duration
det$Duration <- NA
det$Duration[det$Method %in% c(11, 14)] <- 3
det$Duration[det$Method %in% c(12, 13)] <- 1

## resolve species, none, unkn
levels(det$AOU_Code)
det$Spp <- det$AOU_Code
det$Spp[det$AOU_Code == ""] <- "NONE"
det$Spp[grepl("Unknown", as.character(det$ENGLISH.NAME))] <- "NONE"
det$Spp[spp$ORDER[match(det$AOU_Code, spp$CODE)] == "ABIOTIC"] <- "NONE"
det$Spp <- droplevels(det$Spp)
levels(det$Spp)

## format date/time
tmp <- paste(det$RECORDING_DATE, det$RECORD_TIME)
det$Start <- strptime(tmp, "%d-%b-%y %H:%M:%S")

det <- droplevels(det[det$Replicate == 1, ]) # ~40 rows
#det <- det[det$Spp != "NONE", ]

## first detection interval
tmp <- col(det[,c("X0min", "X1min", "X2min")])
tmp[is.na(det[,c("X0min", "X1min", "X2min")])] <- Inf
tmp2 <- find_min(tmp)
tmp2$value[is.infinite(tmp2$value)] <- NA
det$Det1 <- tmp2$value

## make sure not double counted: indiv_id # ~60 rows
tmp <- paste(det$RecordingKey, det$Spp, det$INDIV_ID)
tmp2 <- paste(det$RecordingKey, det$Spp)
dc <- names(table(tmp))[table(tmp) > 1]
zz <- det[tmp %in% dc,]
zz <- zz[zz$Spp != "NONE",]
zz[,c("RecordingKey", "Spp","INDIV_ID", "X0min", "X1min", "X2min")]

#rec <- nonDuplicated(det[det$Replicate == 1, ], RecordingKey, TRUE)
rec <- nonDuplicated(det, RecordingKey, TRUE)
rec <- rec[,c("RecordingKey", "ProjectID", "Cluster", "SITE", "STATION", 
    "Year", "Round", "FileName", "RECORDING_DATE", "RECORD_TIME", 
    "Replicate", "Observer", "Rain", "Wind", "Industry", "Noise", 
    "Microphone", "ProsTime", "Method", "Comment.Recording", "Duration", "Start")]
rec$ToY <- rec$Start$yday
rec$ToD <- rec$Start$hour + rec$Start$min / 60
rec$ToDx <- round(rec$ToD, 0)

xt <- as.matrix(Xtab(~ ToY + ToDx, rec))
xt1 <- Xtab(~ ToY + ToDx + Duration, rec)
xt2 <- as.matrix(Xtab(~ Duration + ToDx, rec))

hist(rec$ToY)
hist(rec$ToD)
hist(rec$ToD[rec$Duration == 1])
hist(rec$ToD[rec$Duration == 3])

barplot(xt2)

aa=table(round(rec$ProsTime),rec$Duration)
aa <- aa[-(1:2),]
barplot(t(aa),beside=T)

table(det$Replicate)
table(det$Replicate, det$Method)

table(det$SITE, det$Duration)


table(rowSums(!is.na(det[det$Duration == 3,c("X0min", "X1min", "X2min")])))

