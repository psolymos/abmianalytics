library(mefa4)
library(pbapply)
library(maptools)

ROOT <- "e:/peter/AB_data_v2016"

tms <- read.csv(file.path(ROOT, "data", "aru-coni",
    "2015_CONIPeent3.4_30_70_truepositives_details.csv"))

tmp <- with(tms, paste0("2015-", month, "-", day, " ", hour, ":", minute, ":", second))
tms$FileStart <- strptime(tmp, format="%Y-%m-%e %H:%M:%S")

tmp <- strsplit(as.character(tms$offset), ":")
OffMin <- sapply(tmp, function(z) as.numeric(z[length(z)-1]))
OffSec <- sapply(tmp, function(z) as.numeric(z[length(z)]))
summary(OffMin)
summary(OffSec)
## start of event relative to FileStart, in minutes
tms$EventStart <- OffMin + OffSec / 60

fls <- nonDuplicated(tms, file.path, FALSE)

## events list for files
evt <- list()
for (j in seq_len(nrow(fls))) {
    i <- as.character(fls$file.path[j])
    #cat(i, "of", nrow(fls));flush.console()
    ss <- tms[tms$file.path == i,]
    ss <- nonDuplicated(ss, offset, FALSE)
    ss <- ss[order(ss$EventStart),]
    if (nrow(ss) < 2) {
        evt[[j]] <- cbind(ID=j, DIFF=NA)
    } else {
        evt[[j]] <- cbind(ID=j, DIFF=diff(ss$EventStart))
    }
}
evt <- do.call(rbind, evt)
evt2 <- evt[!is.na(evt[,2]),]
summary(evt[,2] * 60)
plot(density(evt2[,2]))
polygon(density(evt2[,2]),col="grey")


## Julian day
fls$JULIAN <- fls$FileStart$yday # this is kept as original
fls$JDAY <- fls$JULIAN / 365
summary(fls$JULIAN)
summary(fls$JDAY)

## need XY for tssr -- do not have it yet
## prevent too far extrapolation
PKEY$JDAY[PKEY$JDAY < 0.35 | PKEY$JDAY > 0.55] <- NA
## TSSR = time since sunrise
Coor <- as.matrix(cbind(as.numeric(SS$X),as.numeric(SS$Y)))[match(PKEY$SS, rownames(SS)),]
JL <- as.POSIXct(DD)
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
PKEY$srise <- NA
PKEY$srise[subset] <- sr
PKEY$start_time <- PKEY$HOUR + PKEY$MIN/60
TZ <- SS$TZONE[match(PKEY$SS, rownames(SS))]
lttz <- read.csv("~/repos/bamanalytics//lookup/tzone.csv")
lttz <- nonDuplicated(lttz, Timezone, TRUE)
PKEY$MDT_offset <- lttz$MDT_offset[match(TZ, rownames(lttz))]
table(TZ, PKEY$MDT_offset)
PKEY$TSSR <- (PKEY$start_time - PKEY$srise + PKEY$MDT_offset) / 24
PKEY$TSSR_orig <- PKEY$TSSR # keep a full copy
PKEY$TSSR[PKEY$start_time > 12] <- NA ## after noon
summary(PKEY$TSSR)
summary(PKEY$start_time)

## sin/cos for TSSR
pkDur$TSSRsin <- sin(pkDur$TSSR_orig * 2 * pi)
pkDur$TSSRcos <- cos(pkDur$TSSR_orig * 2 * pi)
pkDur$TSSRsin2 <- sin(pkDur$TSSR_orig * 2 * pi)^2
pkDur$TSSRcos2 <- cos(pkDur$TSSR_orig * 2 * pi)^2

ff <- list(
    "1"= ~ JDAY + JDAY2 + TSSRsin + TSSRcos,
    "2"= ~ JDAY + JDAY2 + TSSRsin + TSSRcos + TSSRsin2 + TSSRcos2,
    "3"= ~ JDAY + JDAY2 + JDAY3 + TSSRsin + TSSRcos,
    "4"= ~ JDAY + JDAY2 + JDAY3 + TSSRsin + TSSRcos + TSSRsin2 + TSSRcos2,
    "5"= ~ JDAY + JDAY2 + JDAY3,
    "6"= ~ TSSRsin + TSSRcos + TSSRsin2 + TSSRcos2)

