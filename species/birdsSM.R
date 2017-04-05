ROOT   <- "e:/peter/AB_data_v2017/data/raw/species"
OUTDIR <- "e:/peter/AB_data_v2017/data/analysis/species"
getwd()
if (interactive())
    source("~/repos/abmianalytics/species/00globalvars.R") else source("00globalvars.R")

T <- "BirdsSM"

det <- read.csv(file.path(ROOT, "birds-sm-20170404.csv"))

## resolve duration
det$Duration <- NA
det$Duration[det$METHOD %in% c("11", "14")] <- 3
det$Duration[det$METHOD %in% c("12", "13")] <- 1

## format date/time
tmp <- paste(det$RECORDING_DATE, det$RECORDING_TIME)
det$Start <- strptime(tmp, "%d-%b-%y %H:%M:%S")

#det <- det[det$Spp != "NONE", ]

## first detection interval
det$int1 <- ifelse(det$MIN_1 == "VNA", NA, as.integer(det$MIN_1))
det$int2 <- ifelse(det$MIN_2 == "VNA", NA, as.integer(det$MIN_2))
det$int3 <- ifelse(det$MIN_3 == "VNA", NA, as.integer(det$MIN_3))
tmp <- col(det[,c("int1", "int2", "int3")])
tmp[is.na(det[,c("int1", "int2", "int3")])] <- Inf
tmp2 <- find_min(tmp)
tmp2$value[is.infinite(tmp2$value)] <- NA
det$Det1 <- tmp2$value

f <- function(x) {
    ifelse(x == "VNA", NA, as.numeric(as.character(x)))
}
det$RAIN <- f(det$RAIN)
det$WIND <- f(det$WIND)
det$INDUSTRY <- f(det$INDUSTRY)
det$NOISE <- f(det$NOISE)
det$MICROPHONE <- f(det$MICROPHONE)

spp_keep <- unique(as.character(det$COMMON_NAME)[det$RANK_NAME == "Species"])
spp_keep <- spp_keep[spp_keep != "VNA"]
det$Spp <- det$COMMON_NAME
levels(det$Spp)[!(levels(det$Spp) %in% spp_keep)] <- "NONE"

levels(det$Spp) <- nameAlnum(levels(det$Spp), capitalize="mixed", collapse="")

## make sure not double counted: indiv_id # ~60 rows
tmp <- paste(det$RECORDING_KEY, det$Spp, det$INDIVIDUAL_ID)
tmp2 <- paste(det$RECORDING_KEY, det$Spp)
dc <- names(table(tmp))[table(tmp) > 1]
zz <- det[tmp %in% dc,]
zz <- zz[zz$Spp != "NONE",]
#zz[,c("RECORDING_KEY", "Spp","INDIVIDUAL_ID", "int1", "int2", "int3")]
## leave it for now (until it is resolved at BU end)

det$site_stn <- interaction(det$SITE, det$STATION, drop=TRUE)

det$ToY <- det$Start$yday
det$ToYc <- as.integer(cut(det$ToY, c(0, 105, 120, 140, 150, 160, 170, 180, 365)))
det$visit <- interaction(det$site_stn, det$ToYc, drop=TRUE)

det$ToD <- det$Start$hour + det$Start$min / 60
det$ToDx <- round(det$ToD, 0)
det$ToDc <- as.factor(ifelse(det$ToDx == 0, "Midnight", "Morning"))


xt_stn <- as.matrix(Xtab(~ site_stn + Spp, det, cdrop="NONE"))
xt_vis <- as.matrix(Xtab(~ visit + Spp, det, cdrop="NONE"))

xt_tod <- data.frame(as.matrix(Xtab(~ Spp + ToDc, det, rdrop="NONE")))
xt_tod$MidP <- round(xt_tod$Midnight / (xt_tod$Midnight + xt_tod$Morning), 4)
xt_tod[order(xt_tod$MidP),]

xt_toy <- as.matrix(Xtab(~ Spp + ToYc, det, rdrop="NONE"))

Class <- nonDuplicated(det[!is.na(det$visit),], visit, TRUE)
Class <- Class[rownames(xt_vis),]
Class$STR2 <- factor(NA, c("A_Early", "B_Mid", "C_Late"))
Class$STR2[Class$ToYc %in% 1:3] <- "A_Early"
Class$STR2[Class$ToYc %in% 4:7] <- "B_Mid"
Class$STR2[Class$ToYc %in% 8] <- "C_Late"
table(Class$STR2, Class$ToYc)

table(det$ToYc, det$Duration)

## crosstab for all-in-one models

keep <- det$Duration == 3 & det$ToDc == "Morning"
keep[is.na(keep)] <- FALSE
det2 <- det[keep,]
det2$PKEY <- interaction(det2$SITE_LABEL, "_", det2$ToY, ":",
    det2$Start$hour, ":", det2$Start$min, sep="", drop=TRUE)
xt <- as.matrix(Xtab(~ PKEY + Spp, det2, cdrop="NONE"))
x <- nonDuplicated(det2, PKEY, TRUE)
x <- x[rownames(xt), c("PKEY", "SITE_LABEL", "ROTATION",
    "SITE", "YEAR", "STATION", "RAIN", "WIND", "INDUSTRY", "NOISE", "MICROPHONE",
    "Start", "ToY", "ToD")]

d <- paste("_", Sys.Date(), sep="")
save(det, xt, x, file=paste(OUTDIR, "/OUT_", tolower(T), d, ".Rdata",sep=""))

## Use the 1st 1-minute segment only

keep <- det$MIN_1 != "VNA"
det3 <- det[keep,]
det3$PKEY <- interaction(det3$SITE_LABEL, "_", det3$ToY, ":",
    det3$Start$hour, ":", det3$Start$min, sep="", drop=TRUE)
xt <- as.matrix(Xtab(~ PKEY + Spp, det3, cdrop="NONE"))
x <- nonDuplicated(det3, PKEY, TRUE)

x$MONTH <- x$Start$mon+1
x$MDAY <- x$Start$mday
x$HOUR <- x$Start$hour
x$MINUTE <- x$Start$min
x$YDAY <- x$Start$yday

x <- x[rownames(xt), c("PKEY", "SITE_LABEL", "ROTATION",
    "SITE", "YEAR", "STATION", "RAIN", "WIND", "INDUSTRY", "NOISE", "MICROPHONE",
    "Start", "ToY", "ToD", "ToYc", "ToDc",
    "MONTH", "MDAY", "HOUR", "MINUTE", "YDAY")]
xt <- xt[,colSums(xt)>0]

save(xt, x, file=paste(OUTDIR, "/OUT_", tolower(T), d, "_1stMinOnly.Rdata",sep=""))
#save(xt, x, file="~/Dropbox/collaborations/opticut/R/abmi-aru-1min.Rdata")


