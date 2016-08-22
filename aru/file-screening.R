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
function(f)
{
    f <- gsub(".wac", "", f)
    f <- substr(sapply(strsplit(f, "\\+"), "[[", 2L), 3, nchar(f))
    as.POSIXlt(strptime(f, "%Y%m%d_%H%M%S"))
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
