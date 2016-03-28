library(mefa4)
library(maptools)

ROOT <- "e:/peter/AB_data_v2016"

load(file.path(ROOT, "data", "aru-coni", 
    "coni-compiled-all.Rdata"))
## int: all intervals from all files
## pt: point level info
## tms: events from positive files
## veghfkm: veg & HF composition in 1 km^2 buffer
## veghfpc: veg & HF composition in 150 m radius circle
## vegkm: backfilled veg composition in 1 km^2 buffer
## vegpc: backfilled veg composition in 150 m radius circle




tmp <- with(tms, paste0("2015-", month, "-", day, " ", hour, ":", minute, ":", second))
tms$FileStart <- strptime(tmp, format="%Y-%m-%e %H:%M:%S")

tmp <- strsplit(as.character(tms$offset), ":")
OffMin <- sapply(tmp, function(z) as.numeric(z[length(z)-1]))
OffSec <- sapply(tmp, function(z) as.numeric(z[length(z)]))
summary(OffMin)
summary(OffSec)
## start of event relative to FileStart, in minutes
tms$EventStart <- OffMin + OffSec / 60

tmp <- strsplit(as.character(int$file.interval), " - ")
int$file.name <- sapply(tmp, "[[", 1)
tmp <- nonDuplicated(int, file.name, FALSE)
tmp$Duration <- ifelse(tmp$Size < 30*10^6, 3, 10)
tms$Duration <- tmp$Duration[match(tms$file.name, tmp$file.name)]
tms <- tms[!is.na(tms$Duration),]

fls <- nonDuplicated(tms, file.name, FALSE)
compare_sets(fls$file.name, tms$file.name)

int2 <- nonDuplicated(int, file.name, TRUE)
tmp <- with(int2, paste0("2015-", month, "-", day, " ", hour, ":", minute, ":", second))
int2$FileStart <- strptime(tmp, format="%Y-%m-%e %H:%M:%S")
int2$JULIAN <- int2$FileStart$yday # this is kept as original
int2$JDAY <- int2$JULIAN / 365
int2$TOD <- (int2$hour + (int2$minute / 60)) / 24

xt <- Xtab(CONI.hit ~ file.name + interval, int)
xt <- xt[rownames(int2),]
plot(TOD*24 ~ JULIAN, int2, pch=".")
points(TOD*24 ~ JULIAN, int2[rowSums(xt)>0,], pch=19, col=2)

if (FALSE) {
## events list for files
evt <- list()
for (j in seq_len(nrow(fls))) {
    i <- as.character(fls$file.name[j])
    #cat(i, "of", nrow(fls));flush.console()
    ss <- tms[tms$file.name == i,]
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
}

## survival model compatible data
evt <- list()
for (j in seq_len(nrow(fls))) {
    i <- as.character(fls$file.name[j])
    #cat(i, "of", nrow(fls));flush.console()
    ss <- tms[tms$file.name == i,]
    ss <- nonDuplicated(ss, offset, FALSE)
    ss <- ss[order(ss$EventStart),]
    ## last event is taken right censored
    ## it is kind of equivalent of left censoring the 1st, but easier to code
    if (nrow(ss) == 1) {
        evt[[j]] <- cbind(ID=j, 
            Start=ss$Duration - ss$EventStart, 
            Event=0)
    } else {
        evt[[j]] <- cbind(ID=j, 
            Start=c(diff(ss$EventStart), ss$Duration[1] - ss$EventStart[nrow(ss)]),
            Event=c(rep(1, nrow(ss)-1), 0))
    }
}
evt <- do.call(rbind, evt)


## Julian day
fls$JULIAN <- fls$FileStart$yday # this is kept as original
fls$JDAY <- fls$JULIAN / 365
summary(fls$JULIAN)
summary(fls$JDAY)

## TOD time-of-day
fls$TOD <- (fls$hour + (fls$minute / 60)) / 24
summary(fls$TOD)

plot(TOD*24 ~ JULIAN, fls)

fls$JDAY2 <- fls$JDAY^2
fls$sin <- sin(fls$TOD * 2 * pi)
fls$cos <- cos(fls$TOD * 2 * pi)
fls$sin2 <- sin(fls$TOD * 2 * pi)^2
fls$cos2 <- cos(fls$TOD * 2 * pi)^2
par(mfrow=c(2,2))
plot(sin ~ TOD, fls)
plot(cos ~ TOD, fls)
plot(sin2 ~ TOD, fls)
plot(cos2 ~ TOD, fls)

fls2 <- fls[evt[,"ID"],]
fls2$Time <- evt[,"Start"] # badly named but it is inter-event time
fls2$Event <- evt[,"Event"] # 1=dead, 0=alive (censored)
str(fls2)

library(survival)
sv <- with(fls2, Surv(Time, Event))
head(sv)
head(evt)


m0 <- survreg(sv ~ 1, fls2, dist="exponential")
m1 <- survreg(sv ~ JDAY, fls2, dist="exponential")
m2 <- survreg(sv ~ JDAY + JDAY2, fls2, dist="exponential")
m3 <- survreg(sv ~ sin + cos + sin2 + cos2, fls2, dist="exponential")
m4 <- survreg(sv ~ JDAY + sin + cos + sin2 + cos2, fls2, dist="exponential")
m5 <- survreg(sv ~ JDAY + JDAY2 + sin + cos + sin2 + cos2, fls2, dist="exponential")
m6 <- survreg(sv ~ JDAY + sin + cos + sin2 + cos2 +
    JDAY:sin + JDAY:sin2 + JDAY:cos + JDAY:cos2, fls2, dist="exponential")
m7 <- survreg(sv ~ JDAY + JDAY2 + sin + cos + sin2 + cos2 +
    JDAY:sin + JDAY:sin2 + JDAY:cos + JDAY:cos2 +
    JDAY:sin + JDAY:sin2 + JDAY:cos + JDAY:cos2, fls2, dist="exponential")

m0w <- survreg(sv ~ 1, fls2, dist="weibull")
m1w <- survreg(sv ~ JDAY, fls2, dist="weibull")
m2w <- survreg(sv ~ JDAY + JDAY2, fls2, dist="weibull")
m3w <- survreg(sv ~ sin + cos + sin2 + cos2, fls2, dist="weibull")
m4w <- survreg(sv ~ JDAY + sin + cos + sin2 + cos2, fls2, dist="weibull")
m5w <- survreg(sv ~ JDAY + JDAY2 + sin + cos + sin2 + cos2, fls2, dist="weibull")
m6w <- survreg(sv ~ JDAY + sin + cos + sin2 + cos2 +
    JDAY:sin + JDAY:sin2 + JDAY:cos + JDAY:cos2, fls2, dist="weibull")
m7w <- survreg(sv ~ JDAY + JDAY2 + sin + cos + sin2 + cos2 +
    JDAY:sin + JDAY:sin2 + JDAY:cos + JDAY:cos2 +
    JDAY:sin + JDAY:sin2 + JDAY:cos + JDAY:cos2, fls2, dist="weibull")

aic <- AIC(m0, m1, m2, m3, m4, m5, m6, m7, m0w, m1w, m2w, m3w, m4w, m5w, m6w, m7w)
aic$dAIC <- aic$AIC - min(aic$AIC)

mb <- m7
mbw <- m7w

vjd <- seq(min(fls2$JDAY), max(fls2$JDAY), len=51*10)
vtd <- seq(0, 23/24, len=24*10)
pr <- expand.grid(JDAY=vjd, TOD=vtd)
pr$JDAY2 <- pr$JDAY^2
pr$sin <- sin(pr$TOD * 2 * pi)
pr$cos <- cos(pr$TOD * 2 * pi)
pr$sin2 <- sin(pr$TOD * 2 * pi)^2
pr$cos2 <- cos(pr$TOD * 2 * pi)^2

fit <- predict(mb, newdata=pr)
z <- matrix(fit, length(vjd), length(vtd))
fitw <- predict(mbw, newdata=pr)
zw <- matrix(fitw, length(vjd), length(vtd))

par(mfrow=c(2,2))
plot(TOD*24 ~ JULIAN, int2, pch=21, cex=0.6, col="grey", main="Exp, rate",
    xlab="Julian day", ylab="Hour")
points(TOD*24 ~ JULIAN, int2[rowSums(xt)>0,], pch=19, cex=0.6, col=1)
contour(vjd*365, vtd*24, 1/z, add=TRUE, col=2)
plot(TOD*24 ~ JULIAN, int2, pch=21, cex=0.6, col="grey", main="Weibull, rate",
    xlab="Julian day", ylab="Hour")
points(TOD*24 ~ JULIAN, int2[rowSums(xt)>0,], pch=19, cex=0.6, col=1)
contour(vjd*365, vtd*24, 1/zw, add=TRUE, col=2)

plot(TOD*24 ~ JULIAN, int2, pch=21, cex=0.6, col="grey", main="Exp, p10",
    xlab="Julian day", ylab="Hour")
points(TOD*24 ~ JULIAN, int2[rowSums(xt)>0,], pch=19, cex=0.6, col=1)
contour(vjd*365, vtd*24, 1-exp(-10*1/z), add=TRUE, col=2)
plot(TOD*24 ~ JULIAN, int2, pch=21, cex=0.6, col="grey", main="Weibull, p10",
    xlab="Julian day", ylab="Hour")
points(TOD*24 ~ JULIAN, int2[rowSums(xt)>0,], pch=19, cex=0.6, col=1)
contour(vjd*365, vtd*24, 1-exp(-(10*1/zw)^(1/mbw$scale)), add=TRUE, col=2)

## offsets for files (int2)
int2$JDAY2 <- int2$JDAY^2
int2$sin <- sin(int2$TOD * 2 * pi)
int2$cos <- cos(int2$TOD * 2 * pi)
int2$sin2 <- sin(int2$TOD * 2 * pi)^2
int2$cos2 <- cos(int2$TOD * 2 * pi)^2
int2$Duration <- ifelse(int2$Size < 30*10^6, 3, 10)

int2$Rate <- 1 / predict(mb, newdata=int2)
int2$p <- 1-exp(-int2$Duration * int2$Rate)

#int2$Ratew <- 1 / predict(mbw, newdata=int2)
#int2$pw <- 1-exp(-(int2$Duration * int2$Ratew)^(1/mbw$scale))

int3 <- int2[,c("file.name","ID","JULIAN","TOD","Rate","Duration","p")]
int3$Nhits <- rowSums(xt)
int3$hit01 <- ifelse(int3$Nhits > 0, 1, 0)
boxplot(p ~ hit01, int3)

plot(p ~ pw, int2, ylim=c(0,1), xlim=c(0,1))
abline(0,1)

write.csv(int3, row.names=FALSE, 
    file=file.path(ROOT, "data", "aru-coni", "offsets.csv"))

coni_data <- data.frame(int3,
    pt[match(int3$ID, pt$ID), c("NRNAME", "NSRNAME", "AHM", "FFP", "MAP", "MAT", "MCMT", "MWMT", 
    "ProjectID", "Cluster", "SITE", "STATION", 
    "ARU", "Year", "DEPLOY_DATE", "RETRIEVE_DATE", "UTM_Zone", "EASTING", 
    "NORTHING")])

lt <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
lt$ec <- paste0(ifelse(lt$TypeCC == "CC", "CC", ""), lt$Type, lt$EC_AGE)
lt$ec[lt$WATER == 1] <- "Water"
lt$ec[lt$ec == "XXX"] <- "Bare"


vpt <- as.matrix(veghfpc)
vpt <- groupSums(vpt, 2, lt[colnames(vpt), "ec"])
vpt <- vpt[match(coni_data$ID, rownames(vpt)),]
colnames(vpt) <- paste0(colnames(vpt), "_pt")

vkm <- as.matrix(veghfkm)
vkm <- groupSums(vkm, 2, lt[colnames(vkm), "ec"])
vkm <- vpt[match(coni_data$ID, rownames(vkm)),]
colnames(vkm) <- paste0(colnames(vkm), "_km")

veghf <- data.frame(vpt, vkm)
rownames(veghf) <- rownames(coni_data)

coni_data <- data.frame(coni_data, veghf)

save(coni_data, 
    file=file.path(ROOT, "data", "aru-coni", "coni-analysis.Rdata"))


## modeling starts here ---------------------------------------------------

coni_use <- coni_data
coni_use$hit01 <- NULL

coni_use$file.name <- NULL
coni_use$ID <- NULL
coni_use$JULIAN <- NULL
coni_use$TOD <- NULL
coni_use$RATE <- NULL
coni_use$Duration <- NULL
coni_use$NRNAME <- NULL
coni_use$NSRNAME <- NULL
coni_use$ProjectID <- NULL
coni_use$Cluster <- NULL
coni_use$SITE <- NULL
coni_use$STATION <- NULL
coni_use$ARU <- NULL
coni_use$YEAR <- NULL
coni_use$DEPLOY_DATE <- NULL
coni_use$RETRIEVE_DATE <- NULL
coni_use$UTM_Zone <- NULL
coni_use$EASTING <- NULL
coni_use$NORTHING <- NULL

library(dismo)

## GLM, Poisson model (use)
mod1 <- glm(Nhits ~ MAT, data=coni_use, offset=log(coni_use$p), family=poisson)

## BRT

mod2 <- gbm.step(data=coni_use, 
    gbm.x=colnames(coni_use)[c(colnames(coni_use) %in% c("p","Nhits"))],
    gbm.y="Nhits",
    offset=coni_data$p,
    family="poisson")
## follow these steps:
## https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf

## other things to try:
## * restrict data set to early morning times (use TOD variable)
##   to mimic that data set has been collected using regular field obs times
## * maybe take a random sample matching that sample size to eliminate 
##   sample size related differences in model.
## * this can give you something for pitching for ARUs


