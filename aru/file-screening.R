library(pbapply)
library(tuneR)
library(seewave)

## problem list
x <- read.csv(file.path("e:/peter/AB_data_v2016/data/aru-raw",
    "2015-ABMI-MC-SC-AudioRecordingTranscription-issues.csv"))

DIRIN <- "u:\\2015_AllData\\ABMI-Core\\ABMI-0032\\ABMI-0032-NW"
#DIRIN <- "u:\\2015_AllData\\ABMI-Core\\ABMI-0235\\ABMI-0235-SW"
DIROUT <- "e:\\Peter\\tmp"
CPROG <- "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe"

aru_convert <-
function(DIRIN, DIROUT, FNAME, CPROG)
{
    ## fin & fout must be full path
    fin <- file.path(DIRIN, FNAME, fsep="\\")
    fout <- file.path(DIROUT, gsub(".wac", ".wav", FNAME), fsep="\\")
    CALL <- paste(CPROG, fin, fout)
    system(CALL)
    invisible(fout)
}
aru_parse_location <-
function(f)
{
    f <- gsub(".wac", "", f)
    f <- sapply(strsplit(f, "\\+"), "[[", 1L)
    substr(f, 1, nchar(f)-2)
}
aru_parse_date <-
function(f, type=c("full", "date", "time"))
{
    type <- match.arg(type)
    f <- gsub(".wac", "", f)
    if (type == "full") {
        f <- substr(sapply(strsplit(f, "\\+"), "[[", 2L), 3, nchar(f))
        out <- as.POSIXlt(strptime(f, "%Y%m%d_%H%M%S"))
    } else {
        out <- sapply(strsplit(f, "_"), "[[", if (type == "date") 3L else 4L)
    }
    out
}
aru_calculateH <-
function (wave, f, wl = 512, envt = "hil", msmooth = NULL, ksmooth = NULL)
{
    input <- inputw(wave = wave, f = f)
    wave <- input$w
    f <- input$f
    rm(input)
    spec <- meanspec(wave = wave, f = f, wl = wl, plot = FALSE)
    SH <- sh(spec)
    enve <- env(wave = wave, f = f, envt = envt, msmooth = msmooth,
        ksmooth = ksmooth, plot = FALSE)
    TH <- th(enve)
    H <- SH * TH
    return(c(H=H, SH=SH, TH=TH))
}
aru_leftright <- function(object) {
    d <- object@left - object@right
    sdl <- sd(object@left)
    sdr <- sd(object@right)
    c(sd_left=sdl, sd_right=sdr,
        mean_left=mean(object@left), mean_right=mean(object@right),
        min_left=min(object@left), min_right=min(object@right),
        max_left=max(object@left), max_right=max(object@right),
        std_diff = mean(d / (0.5*(sdl+sdr))),
        log_sd_ratio = log(sd(object@left) / sd(object@right)))
}
aru_summarize <-
function(DIRIN, DIROUT, FNAME, CPROG, keep_wav=FALSE, quick=FALSE, return_object=TRUE)
{
    t0 <- proc.time()
    if (!file.exists(file.path(DIRIN, FNAME, fsep="\\")))
        stop("file does not exist")
    FNWAV <- try(aru_convert(DIRIN, DIROUT, FNAME, CPROG))
    if (inherits(FNWAV, "try-error"))
        stop("wac2wav error")
    object <- readWave(FNWAV)
    #object <- if (sample_p < 1)
    #    downsample(object0, samp.rate=sample_p*object@samp.rate) else object0
    #SPECTRO <- spectro(object, plot=FALSE)
    #CSH <- csh(object, plot=FALSE) # range, mean etc ?
    H <- if (quick)
        NULL else aru_calculateH(object)
    out <- list(file=FNAME,
        timer = unname(proc.time() - t0)[3L],
        location = aru_parse_location(FNAME),
        date = aru_parse_date(FNAME),
        duration = length(object@left)/object@samp.rate,
        size_wac = file.size(file.path(DIRIN, FNAME, fsep="\\")),
        size_wav = file.size(FNWAV),
        entropy=H,
        channels=aru_leftright(object),
        #sample_p=sample_p,
        object=if (return_object) object else NULL)
    if (!keep_wav)
        file.remove(FNWAV)
    out
}

if (FALSE) {
FLWAC <- list.files(DIRIN)
FLWAC <- FLWAC[grep(".wac", FLWAC)]

LOC <- aru_parse_location(FLWAC[1L])
DATES <- aru_parse_date(FLWAC)
SIZES <- pbsapply(FLWAC, function(z) file.size(file.path(DIRIN, z, fsep="\\")))

#FNAME <- FLWAC[which.min(SIZES)]
FNAME <- "ABMI-0235-SW_0+1_20150522_055400.wac"
#FNAME <- "ABMI-0032-NW_0+1_20150312_142828.wac"
#FNAME <- "ABMI-0032-NW_0+1_20150312_150000.wac"
#FNAME <- "ABMI-0032-NW_A_Summary.txt"


res <- aru_summarize(
    FNAME = "ABMI-0032-NW_0+1_20150723_020255.wac",
    DIRIN = "u:\\2015_AllData\\ABMI-Core\\ABMI-0032\\ABMI-0032-NW",
    DIROUT = "e:\\Peter\\tmp",
    CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
    keep_wav = TRUE)

system.time(res <- aru_summarize(
    FNAME = "ABMI-0388-SW_0+1_20150626_071600.wac",
    DIRIN = "u:\\2015_AllData\\ABMI-Core\\ABMI-0388\\ABMI-0388-SW",
    DIROUT = "e:\\Peter\\tmp",
    CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
    keep_wav = TRUE))

## --



FWAV <- aru_convert(DIRIN, DIROUT, FNAME, CPROG)
fi <- file.info(FWAV)
fs <- file.size(FWAV)

#FNAME <- FLWAC[which.min(SIZES)]
FNAME <- "ABMI-0235-SW_0+1_20150522_055400.wac"
#FNAME <- "ABMI-0032-NW_0+1_20150312_142828.wac"
#FNAME <- "ABMI-0032-NW_0+1_20150312_150000.wac"
#FNAME <- "ABMI-0032-NW_A_Summary.txt"

## too small
"ABMI-0032-NW_0+1_20150723_020255.wac"
## very busy
"ABMI-0034-NE_0+1_20150520_055700.wac"
## good practice rec
"ABMI-0235-SW_0+1_20150522_055400.wac"
## full of birds
"ABMI-0388-SW_0+1_20150626_071600.wac"
## Right microphone already dead on deployment recording
"ABMI-0331-NW_0+1_20150301_090000.wac"

## short file
res0 <- aru_summarize(
    FNAME = "ABMI-0032-NW_0+1_20150723_020255.wac",
    DIRIN = "u:\\2015_AllData\\ABMI-Core\\ABMI-0032\\ABMI-0032-NW",
    DIROUT = "e:\\Peter\\tmp",
    CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
    keep_wav = TRUE, quick=TRUE)
## Right microphone already dead on deployment recording
system.time(res1 <- aru_summarize(
    FNAME = "ABMI-0331-NW_0+1_20150301_090000.wac",
    DIRIN = "u:\\2015_AllData\\ABMI-Core\\ABMI-0331\\ABMI-0331-NW",
    DIROUT = "e:\\Peter\\tmp",
    CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
    keep_wav = TRUE, quick=TRUE))
## OK file
system.time(res2 <- aru_summarize(
    FNAME = "ABMI-0388-SW_0+1_20150626_071600.wac",
    DIRIN = "u:\\2015_AllData\\ABMI-Core\\ABMI-0388\\ABMI-0388-SW",
    DIROUT = "e:\\Peter\\tmp",
    CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
    keep_wav = TRUE, quick=TRUE))

object <- res1$object
## downsample
#o2 <- downsample(object, samp.rate=0.2*object@samp.rate)

aru_calculateH(object[,1])
aru_calculateH(object[,2])
}

## (1) define duration wac_size relationship based on random sample

## randomly picking 100 files
dir <- c("u:\\2015_AllData\\ABMI-Core\\ABMI-0001",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0002",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0003",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0013",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0014",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0015",
    #"u:\\2015_AllData\\ABMI-Core\\ABMI-0032",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0033",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0034",
    #"u:\\2015_AllData\\ABMI-Core\\ABMI-0234",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0235",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0236",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0264",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0265",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0266",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0294",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0295",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0296",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0329",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0330",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0331",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0358",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0359",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0360",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0388",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0389",
    #"u:\\2015_AllData\\ABMI-Core\\ABMI-0390", # wav in the mix
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0412",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0413",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0414",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0442",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0443",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0444",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0472",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0473",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0474",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0512",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0514",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0543",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0544",
    #"u:\\2015_AllData\\ABMI-Core\\ABMI-0545",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0574",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0575",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0576",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0599",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0600",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0601",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0631",
    #"u:\\2015_AllData\\ABMI-Core\\ABMI-0632",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0633",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0663",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0664",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0665",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0791",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0792",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0793",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0794",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0823",
    #"u:\\2015_AllData\\ABMI-Core\\ABMI-0824",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0825",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0855",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0856",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0857",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0989",
#    "u:\\2015_AllData\\ABMI-Core\\ABMI-0990", # unexpected dir name
    "u:\\2015_AllData\\ABMI-Core\\ABMI-0991",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-1022",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-1023",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-1024",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-1056",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-1057",
    "u:\\2015_AllData\\ABMI-Core\\ABMI-1058")

set.seed(1)
#rnd <- c("NE","NW","SE","SW")[sample(1:4, length(dir), replace=TRUE)]
#site <- sapply(strsplit(dir, "\\", fixed=TRUE), "[[", 4)
#dirs <- paste0(dir, "\\", site, "-", rnd)
dirs <- rep(dir, each=4)
site <- sapply(strsplit(dirs, "\\", fixed=TRUE), "[[", 4)
dirs <- paste0(dirs, "\\", site, "-", c("NE","NW","SE","SW"))

files <- character(length(dirs))
for (i in 1:length(dirs)) {
    cat(i, "/", length(dirs), "\n");flush.console()
    FLWAC <- list.files(dirs[i])
    FLWAC <- FLWAC[grep(".wac", FLWAC)]
    files[i] <- sample(FLWAC, 1)
}

save(dirs, files, file="e:\\peter\\AB_data_v2016\\out\\aru-filtering\\dirs-files.Rdata")

RES <- list()
for (i in 1:length(dirs)) {
    cat(i, "/", length(dirs), "\n");flush.console()
    res <- aru_summarize(
        FNAME = files[i],
        DIRIN = dirs[i],
        DIROUT = "e:\\Peter\\tmp",
        CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
        keep_wav = FALSE, quick=TRUE, return_object=FALSE)
    RES[[i]] <- res
}

save(RES, file="e:\\peter\\AB_data_v2016\\out\\aru-filtering\\sample.Rdata")

z <- data.frame(t(sapply(1:length(RES), function(i) unlist(RES[[i]][c(2,5,6,7)]))))
plot(z)

## wac is 34-75% compression of original wav
summary(z$size_wac / z$size_wav)

## assume 3 min as minimum
## some <180, but mostly 180 and 600 sec
table(z$duration)
(z0 <- z[z$duration < 180,])
summary(z[z$duration >= 180,])
summary(z[z$duration == 180,])
summary(z[z$duration == 600,])

#(WAC_SIZE_MIN <- max(min(z$size_wac[z$duration >= 180]), max(z$size_wac[z$duration < 180])))
(WAC_SIZE_MIN0 <- min(z$size_wac[z$duration >= 180]))
WAC_SIZE_MIN <- WAC_SIZE_MIN0

## (2) assess 1 directory at a time and find small files

#WAC_SIZE_MIN <- 11849974 # 3 minutes minimum
WAC_SIZE_MIN <- 17825792 # some 3 mins in (17Mb)
RESbyDIR <- list()
ii <- 1:length(dirs)
#ii <- matrix(1:length(dirs), nrow=4)[1,]
for (i in ii) {
    DIR <- dirs[i]
    cat(DIR, "---", which(ii==i), "/", length(ii), "\n");flush.console()
    FLWAC <- list.files(DIR)
    FLWAC <- FLWAC[grep(".wac", FLWAC)]

    #LOC <- aru_parse_location(FLWAC[1L])
    dd <- data.frame(file=FLWAC,
        full_date=aru_parse_date(FLWAC, type="full"),
        date=aru_parse_date(FLWAC, type="date"),
        time=aru_parse_date(FLWAC, type="time"),
        size_wac=pbsapply(FLWAC, function(z) file.size(file.path(DIR, z, fsep="\\"))))
    rownames(dd) <- NULL

    si <- dd$size_wac / 1024^2
    names(si) <- dd$file
    si <- sort(si[si < 40])
    sid <- cbind(Theoretical=qnorm(seq(0.0001,0.9999,len=length(si))),
        Empirical=(si-mean(si))/sd(si),
        pval=pnorm(-abs((si-mean(si))/sd(si)))) # 1-sided
    rownames(sid) <- names(si)
    sid[sid[,3] < 0.05,]
    #dd$pval <- NA
    dd$pval <- sid[match(dd$file, rownames(sid)),3]
    dd$flag <- FALSE
    dd$flag[!is.na(dd$pval) & dd$pval < 0.05] <- TRUE
    dd <- dd[order(dd$size_wac),]
    dd$flag[which.max(dd$size_wac[dd$flag])+1] <- TRUE
    dd <- dd[order(dd$full_date),]

#    dd$small <- dd$size_wac < WAC_SIZE_MIN
    dd$small <- dd$flag
    dd$duration <- NA


    tmp <- list()
    k <- 1
    sss <- which(dd$small)
    if (length(sss) > 20) {
        sss <- sample(sss, 20)
        #sss <- c(sss[1:10], sss[(length(sss)-9):length(sss)])
    }
    for (j in sss) {
        cat(k, "/", length(sss), "\n");flush.console()
        res <- try(aru_summarize(
            FNAME = as.character(dd$file[j]),
            #FNAME = "ABMI-0264-NW_0+1_20150421_195400.wac",
            DIRIN = DIR,
            DIROUT = "e:\\Peter\\tmp",
            CPROG = "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe",
            keep_wav = FALSE, quick=TRUE, return_object=FALSE))
        if (!inherits(res, "try-error")) {
#            tmp[[j]] <- res
            dd$duration[j] <- res$duration
        }
        k <- k + 1
        rm(res)
    }
    print(table(dd$duration))

    ## update min
#    ddd <- dd[!is.na(dd$duration),]
#    ddd <- ddd[ddd$duration >= 180,]
#    WAC_SIZE_MIN_NEW <- min(ddd$size_wac)
#    WAC_SIZE_MIN <- min(WAC_SIZE_MIN, WAC_SIZE_MIN_NEW)

    RESbyDIR[[DIR]] <- dd
}

save(RESbyDIR, file="e:\\peter\\AB_data_v2016\\out\\aru-filtering\\filter-dirs.Rdata")

zz1 <- do.call(rbind, lapply(RESbyDIR, function(z) {
    z <- z[!is.na(z$duration),,drop=FALSE]
    z[z$duration < 180,,drop=FALSE]
}))
zz2 <- do.call(rbind, lapply(RESbyDIR, function(z) {
    z <- z[!is.na(z$duration),,drop=FALSE]
    z[z$duration >= 180,,drop=FALSE]
}))

range_fun <- function(d1) {
    d1$durc <- 10
    d1$durc[d1$size_wac/1024^2 < 40] <- 3
    d1$durc[!is.na(d1$duration) & d1$duration < 180] <- 0
    table(d1$durc)
    t(sapply(c(0, 3, 10), function(z) range(d1$size_wac[d1$durc==z]/1024^2)))
}
ranges <- lapply(RESbyDIR, range_fun)
r0 <- t(sapply(ranges, function(z) z[1,]))
r0[!is.finite(r0)] <- NA
r3 <- t(sapply(ranges, function(z) z[2,]))
r3[!is.finite(r3)] <- NA
r10 <- t(sapply(ranges, function(z) z[3,]))
r10[!is.finite(r10)] <- NA

with(rbind(zz1, zz2), plot(size_wac/1024^2, duration,
    col=c(rep(2, nrow(zz1)), rep(1, nrow(zz2)))))

min(zz2$size_wac/1024^2) # 11.2
max(zz1$size_wac/1024^2) # 17.7
table(zz2$size_wac > max(zz1$size_wac))
table(zz1$size_wac < min(zz2$size_wac))

## confusion tables
zz <- rbind(zz1, zz2)
True <- zz$duration < 180
Min3min <- zz$size_wac < min(zz2$size_wac)
MaxSmall <- zz$size_wac < max(zz1$size_wac)
n <- nrow(zz)

##  True pos           | False neg (Type II)
## --------------------+---------------------
##  False pos (Type I) | True neg
(cm1 <- round((table(True=True, Pred=Min3min)/n)[c(2,1),c(2,1)],2))
##        Pred
## True    TRUE FALSE
##   TRUE  0.56  0.04
##   FALSE 0.00  0.40
(cm2 <- round((table(True=True, Pred=MaxSmall)/n)[c(2,1),c(2,1)],2))
##        Pred
## True    TRUE FALSE
##   TRUE  0.60  0.00
##   FALSE 0.05  0.36

## Accuracy = (True pos + True neg) / n
(Acc1 <- sum(diag(cm1)))
(Acc2 <- sum(diag(cm2)))

hist(zz1$size_wac/1024^2, col="grey")
rug(zz1$size_wac/1024^2)
abline(v=min(zz2$size_wac/1024^2), col=2)

hist(zz2$size_wac/1024^2, col="grey")
rug(zz2$size_wac/1024^2)
abline(v=max(zz1$size_wac/1024^2), col=2)

plot(0, type="n", xlim=c(0, max(r10, na.rm=TRUE)), ylim=c(0, length(RESbyDIR)+1),
    xlim="WAC size (Mb)", ylim="ARU unit")
for (i in 1:length(RESbyDIR)) {
    lines(r0[i,], c(i,i), col=2)
    lines(r3[i,], c(i,i), col=1)
    lines(r10[i,], c(i,i), col=4)
}




## have a look at site level variation in file sizes that are OK (3 & 10 min)
## also check min 3 minutes limit

## (3) plot

#source("~/repos/abmianalytics/aru/qcc_functions.R")

i <- 1
    DIR <- dirs[i]
    cat(DIR, "---", i, "/", length(dirs), "\n");flush.console()
    FLWAC <- list.files(DIR)
    FLWAC <- FLWAC[grep(".wac", FLWAC)]

    LOC <- aru_parse_location(FLWAC[1L])
    dd <- data.frame(file=FLWAC,
        full_date=aru_parse_date(FLWAC, type="full"),
        date=aru_parse_date(FLWAC, type="date"),
        time=aru_parse_date(FLWAC, type="time"),
        size_wac=pbsapply(FLWAC, function(z) file.size(file.path(DIR, z, fsep="\\"))))
    rownames(dd) <- NULL
    dd$small <- dd$size_wac < WAC_SIZE_MIN0
    dd$yday <- as.POSIXlt(dd$full_date)$yday

    dd$date2 <- as.POSIXlt(strptime(dd$date, "%Y%m%d"))
    dd$time2 <- as.POSIXlt(strptime(dd$time, "%H%M%S"))


si <- dd$size_wac / 1024^2
names(si) <- dd$file
si <- sort(si[si < 40])
sid <- cbind(Theoretical=qnorm(seq(0.0001,0.9999,len=length(si))),
    Empirical=si,
    pval=pnorm(-abs((si-mean(si))/sd(si)))) # 1-sided
rownames(sid) <- names(si)
sid[sid[,3] < 0.05,]
#dd$pval <- NA
dd$pval <- sid[match(dd$file, rownames(sid)),3]
dd$flag <- FALSE
dd$flag[!is.na(dd$pval) & dd$pval < 0.05] <- TRUE
dd <- dd[order(dd$size_wac),]
dd$flag[which.max(dd$size_wac[dd$flag])+1] <- TRUE
dd <- dd[order(dd$full_date),]

plot(sid[,-3], pch=19, col=ifelse(sid[,3] < 0.05, 2, 1))
qqline(sid[,2], col=4)
#qqplot(qnorm(seq(0.0001,0.9999,len=length(si))), si)
#qqnorm(si)


hist(dd$size_wac/1024^2, xlab="WAC file size (Mb)", main=LOC, col="grey")
rug(dd$size_wac/1024^2)
abline(v=WAC_SIZE_MIN0/1024^2, col=2)

with(dd, plot(full_date, size_wac/1024^2, type="l"))
abline(h=WAC_SIZE_MIN0/1024^2, col=2)

#mav <- function(x, n = 8) filter(x, rep(1 / n, n), sides = 2)
#dd$mav <- mav(dd$size_wac/1024^2)
tmp <- as.matrix(aggregate(dd$size_wac, list(dd$yday), range))
tmp <- tmp[match(as.character(dd$yday), as.character(tmp[,1])),]
dd$dmin <- tmp[,2]
dd$dmax <- tmp[,3]
tmp <- as.matrix(aggregate(dd$size_wac, list(dd$yday), mean))
tmp <- tmp[match(as.character(dd$yday), as.character(tmp[,1])),]
dd$dmean <- tmp[,2]

op <- par(mfrow=c(2,2))
plot(dd$date2, dd$time2, pch=19, col=1, main=LOC,
    xlab="Date", ylab="Time")
with(dd[dd$pval < 0.05,], points(date2, time2, pch=19, col=2))

plot(sid[,-3], type="l")
points(sid[,-3], pch=19, col=ifelse(sid[,3] < 0.05, 2, 1))
qqline(sid[,2], col=4)

with(dd[dd$yday <= min(dd$yday)+7,], plot(full_date, size_wac/1024^2, type="n",
    main="", ylab="WAC file size (Mb)", xlab="1 week after deployment",
    ylim=c(0, max(dd$size_wac/1024^2))))
polygon(c(dd$full_date, rev(dd$full_date)),
    c(dd$dmin/1024^2, rev(dd$dmax/1024^2)), col="lightblue", border=NA)
with(dd, lines(full_date, size_wac/1024^2, col=1))
with(dd, points(full_date, size_wac/1024^2, pch=19,
    col=ifelse(pval < 0.05, 2, 1)))
abline(h=WAC_SIZE_MIN0/1024^2, col=2)
#lines(dd$full_date, dd$dmean/1024^2, col=4)

with(dd[dd$yday >= max(dd$yday)-7,], plot(full_date, size_wac/1024^2, type="n",
    main="", ylab="WAC file size (Mb)", xlab="1 week before retrieval",
    ylim=c(0, max(dd$size_wac/1024^2))))
polygon(c(dd$full_date, rev(dd$full_date)),
    c(dd$dmin/1024^2, rev(dd$dmax/1024^2)), col="lightblue", border=NA)
with(dd, lines(full_date, size_wac/1024^2, col=1))
with(dd, points(full_date, size_wac/1024^2, pch=19,
    col=ifelse(pval < 0.05, 2, 1)))
abline(h=WAC_SIZE_MIN0/1024^2, col=2)
#lines(dd$full_date, dd$dmean/1024^2, col=4)

par(op)

plot_control_chart1(dd$full_date, dd$size_wac/1024^2)
type=c("qcc","cusum","ewma"), n=8, nsigmas=4.5,
main, use.date=TRUE, offset=0.2, lambda=0.2, sort=TRUE)

with(dd[dd$yday >= max(dd$yday)-7,], plot_control_chart1(full_date, size_wac/1024^2,
    "cusum", n=floor(0.9*nrow(dd)), use.date=FALSE))




