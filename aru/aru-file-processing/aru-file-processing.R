## install R and the following packages: pbapply, tuneR, seewave
## set variables in the preamble

## --- preamble starts here ---
## set path to the directory you want to screen
## program used value defined here if called interactively
## otherwise this should be 1st argument passed from command line
DIRIN <- "u:\\2015_AllData\\ABMI-Core\\ABMI-0032\\ABMI-0032-NW"
## set path for temporary/intermediate/output files
TMPDIR <- "e:\\Peter\\tmp"
## set path to the exacutable wac2wav converter
CPROG <- "C:\\Users\\Peter\\Dropbox\\abmi\\aru\\wac2wav\\wac2wav.exe"
## set if converted wav files should be kept
KEEP_WAV <- FALSE
## set if conversion needs to be made, else the unverified summaries are returned
## conversion is needed for L/R channel testing
CONVERT <- TRUE
## number of files to check for duration (usually the last NMAX in time sequence)
NMAX <- 20
## --- preamble ends here ---


## --- function and final variable definitions start here ---
#library(pbapply)
library(tuneR)
library(seewave)
ARGS <- if (interactive())
    NULL else commandArgs(trailingOnly = TRUE)
DIRIN <- if (interactive())
    DIRIN else ARGS[1]

aru_convert <-
function(DIRIN, DIROUT, FNAME, CPROG)
{
    ## fin & fout must be full path
    fin <- file.path(DIRIN, FNAME, fsep="\\")
    ## check extension
    EXT <- tolower(substr(FNAME, nchar(FNAME)-3, nchar(FNAME)))
    if (!(EXT %in% c(".wac",".wav")))
        stop("input file must be WAC or WAV")
    if (EXT == ".wac") { # convert
        fout <- file.path(DIROUT, gsub(".wac", ".wav", FNAME), fsep="\\")
        CALL <- paste(CPROG, fin, fout)
        system(CALL)
    } else { # or copy
        fout <- file.path(DIROUT, FNAME, fsep="\\")
        file.copy(fin, fout, overwrite=TRUE)
    }
    invisible(fout)
}
aru_parse_location <-
function(f)
{
    #f <- gsub(".wac", "", f)
    f <- substr(f, 1, nchar(f)-4)
    f <- sapply(strsplit(f, "\\+"), "[[", 1L)
    substr(f, 1, nchar(f)-2)
}
aru_parse_date <-
function(f, type=c("full", "date", "time"))
{
    type <- match.arg(type)
    #f <- gsub(".wac", "", f)
    f <- substr(f, 1, nchar(f)-4)
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
aru_parse_directory <-
function(DIRIN)
{
    FLIST <- list.files(DIRIN)
    EXT <- tolower(substr(FLIST, nchar(FLIST)-3, nchar(FLIST)))
    KEEP <- grepl(".wac", tolower(FLIST)) | grepl(".wav", tolower(FLIST))
    FLIST <- FLIST[KEEP]
    dd <- data.frame(file=FLIST,
        full_date=aru_parse_date(FLIST, type="full"),
        date=aru_parse_date(FLIST, type="date"),
        time=aru_parse_date(FLIST, type="time"),
        size_wac=sapply(FLIST, function(z) file.size(file.path(DIRIN, z, fsep="\\"))))
    dd <- dd[order(dd$full_date),]
    dd$yday <- as.POSIXlt(dd$full_date)$yday
    dd$wday <- as.POSIXlt(dd$full_date)$wday + 1
    dd$mon <- as.POSIXlt(dd$full_date)$mon + 1
    dd$mday <- as.POSIXlt(dd$full_date)$mday
    dd$chrtime <- factor(NA, c("Midnight","2am","SunRise1","SunRise2",
        "Noon","3pm","SunSet1","SunSet2","OUT_OF_SEQUENCE"))
    dd$chrtime[dd$time=="000000"] <- "Midnight"
    dd$chrtime[dd$time=="020000"] <- "2am"
    dd$chrtime[dd$time=="120000"] <- "Noon"
    dd$chrtime[dd$time=="150000"] <- "3pm"
    IDS <- which(is.na(dd$chrtime))
    IDS <- IDS[IDS <= nrow(dd)-2]
    IDS <- IDS[IDS != 1]
    for (i in IDS) {
        if (!is.na(dd$chrtime[i+2])) {
            if (dd$chrtime[i-1] %in% c("First", "2am") && dd$chrtime[i+2] == "Noon")
                dd$chrtime[i] <- "SunRise1"
            if (dd$chrtime[i-1] %in% c("First", "SunRise1") && dd$chrtime[i+1] == "Noon")
                dd$chrtime[i] <- "SunRise2"
            if (dd$chrtime[i-1] %in% c("First", "3pm") && dd$chrtime[i+2] == "Midnight")
                dd$chrtime[i] <- "SunSet1"
            if (dd$chrtime[i-1] %in% c("First", "SunSet1") && dd$chrtime[i+1] == "Midnight")
                dd$chrtime[i] <- "SunSet2"
        }
    }
    dd$chrtime[is.na(dd$chrtime)] <- "OUT_OF_SEQUENCE"
    dd$chrtime[c(1,nrow(dd)-1, nrow(dd))] <- NA
    #table(dd$chrtime,useNA="a")
    dd$duration_class <- factor(NA, c("0-180sec","180-600sec",
        "600sec"))
    Mb <- dd$size_wac / 1024^2
    dd$duration_class[Mb > 50] <- "600sec"
    dd$duration_class[Mb > 30 & Mb <= 50] <- "180-600sec"
    dd$duration_class[Mb <= 30] <- "0-180sec"
    #table(dd$chrtime,dd$duration_class,useNA="a")
    dd
}
aru_flag_files <-
function(dat, nmax=20)
{
    Mb <- dat$size_wac / 1024^2
    names(Mb) <- dat$file
    Mb <- sort(Mb[Mb < 40])
    sid <- cbind(Theoretical=qnorm(seq(0.0001,0.9999,len=length(Mb))),
        Empirical=(Mb-mean(Mb))/sd(Mb),
        pval=pnorm(-abs((Mb-mean(Mb))/sd(Mb)))) # 1-sided
    rownames(sid) <- names(Mb)
    dat$pval <- sid[match(dat$file, rownames(sid)),3]
    dat$to_check <- 0L
    dat$to_check[!is.na(dat$pval) & dat$pval < 0.05] <- -1L
    dat <- dat[order(dat$size_wac),]
#    dat$to_check[which.max(dat$size_wac[dat$to_check])+1] <- TRUE
    dat <- dat[order(dat$full_date),]
    OK <- which(dat$to_check >= 0 & dat$chrtime == "SunRise1" &
        dat$duration_class == "600sec")
    OKJune <- which(dat$to_check >= 0 & dat$chrtime == "SunRise1" &
        dat$duration_class == "600sec" & dat$mon == 6)
    ID <- unique(sort(c(1, OK[1], OKJune[1], OK[length(OK)])))
    dat$to_check[ID] <- 1L
    if (sum(dat$to_check < 0) > nmax) {
        tmp <- which(dat$to_check < 0)
#        tmp1 <- tmp[1:(length(tmp)-nmax)]
#        tmp2 <- tmp[(length(tmp)-nmax+1):length(tmp)]
        dat$to_check[tmp[1:(length(tmp)-nmax)]] <- 0L
    }
    dat
}
aru_plot <-
function(dat) {
    dat$date2 <- as.POSIXlt(strptime(dat$date, "%Y%m%d"))
    dat$time2 <- as.POSIXlt(strptime(dat$time, "%H%M%S"))
    LOC <- aru_parse_location(rownames(dat)[1])
    PROBLM <- dat$to_check < 0 | dat$chrtime == "OUT_OF_SEQUENCE"
    xx <- dat[!is.na(dat$channel_left_mean),,drop=FALSE]
    xxx <- xx[xx$to_check > 0,,drop=FALSE]
    #xxx <- xx
    xval <- 1:nrow(xxx)
    #xval <- xxx$yday

    left <- cbind(xxx$channel_left_mean,
        xxx$channel_left_mean-xxx$channel_left_sd,
        xxx$channel_left_mean+xxx$channel_left_sd)
    right <- cbind(xxx$channel_right_mean,
        xxx$channel_right_mean-xxx$channel_right_sd,
        xxx$channel_right_mean+xxx$channel_right_sd)

    op <- par(mfrow=c(2,2))

    plot(dat$date2, dat$time2, pch=19, col=1, main=LOC,
        xlab="Date", ylab="Time", cex=1)
    with(dat[PROBLM,], points(date2, time2, pch=19, col=2))
    with(dat, points(date2[1], time2[1], pch=19, col=4))

    plot(xval-0.1, left[,1], col=2, pch=19, ylim=range(left,right),
        xlim=c(min(xval)-0.5, max(xval)+0.5), type="n",
        ylab="Amplitude (dB) [mean +/- SD]", xlab="Julian day", axes=FALSE)
    axis(2)
    axis(1, at=xval, xxx$yday)
    box()
    segments(xval-0.1, y0=left[,2], y1=left[,3], col="tomato", lwd=6)
    segments(xval+0.1, y0=right[,2], y1=right[,3], col="lightblue", lwd=6)
    points(xval-0.1, left[,1], col=2, pch=19, ylim=range(left,right), cex=1.5)
    points(xval+0.1, right[,1], col=4, pch=19, ylim=range(left,right), cex=1.5)
    lines(xval-0.1, left[,1], col=2, lty=1, lwd=2)
    lines(xval+0.1, right[,1], col=4, lty=1, lwd=2)
    legend("bottomright", pch=19, lty=1, lwd=2, col=c(2,4),
        legend=c("Left","Right"), bty="n")

    dd1 <- dat[!is.na(dat$yday) & dat$yday <= min(dat$yday,na.rm=TRUE)+7,]
    if (nrow(dd1) < 10)
        dd1 <- dat[!is.na(dat$yday) & dat$yday <= min(dat$yday,na.rm=TRUE)+14,]

    with(dd1, plot(full_date, size_wac/1024^2, type="n",
        main="", ylab="WAC file size (Mb)", xlab="1 week after deployment",
        ylim=c(0, max(dat$size_wac/1024^2))))
    with(dat, lines(full_date, size_wac/1024^2, col=1))
    with(dat[!PROBLM,], points(full_date, size_wac/1024^2, pch=19, col=1))
    with(dat[PROBLM,], points(full_date, size_wac/1024^2, pch=19, col=2))
    with(dat[1,,drop=FALSE], points(full_date, size_wac/1024^2, pch=19, col=4))

    dd2 <- dat[!is.na(dat$yday) & dat$yday >= max(dat$yday,na.rm=TRUE)-7,]
    if (nrow(dd2) < 10)
        dd2 <- dat[!is.na(dat$yday) & dat$yday >= max(dat$yday,na.rm=TRUE)-14,]

    with(dd2, plot(full_date, size_wac/1024^2, type="n",
        main="", ylab="WAC file size (Mb)", xlab="1 week before retrieval",
        ylim=c(0, max(dat$size_wac/1024^2))))
    with(dat, lines(full_date, size_wac/1024^2, col=1))
    with(dat[!PROBLM,], points(full_date, size_wac/1024^2, pch=19, col=1))
    with(dat[PROBLM,], points(full_date, size_wac/1024^2, pch=19, col=2))

    par(op)
    invisible(NULL)
}
## --- function and final variable definitions end here ---

dat1 <- aru_parse_directory(DIRIN)
dat2 <- aru_flag_files(dat1, nmax=NMAX)
dat3 <- dat2
dat3$duration <- NA
dat3$size_wav <- NA
dat3$channel_left_mean <- NA
dat3$channel_right_mean <- NA
dat3$channel_left_sd <- NA
dat3$channel_right_sd <- NA

if (CONVERT) {
    files <- rownames(dat2)[dat2$to_check != 0]
    RES <- list()
    cat("Directory:", DIRIN, "\nStarted at", date(), "\n")
    for (i in seq_len(length(files))) {
        cat(i, "/", length(files), "\n")
        flush.console()
        #OBJECT <- dat2[files[i], "to_check"] > 0
        OBJECT <- FALSE
        RES[[files[i]]] <- aru_summarize(
            FNAME = files[i],
            DIRIN = DIRIN,
            DIROUT = TMPDIR,
            CPROG = CPROG,
            keep_wav = KEEP_WAV,
            quick=TRUE,
            return_object=OBJECT)
        dat3[files[i], "duration"] <- RES[[files[i]]]$duration
        dat3[files[i], "size_wav"] <- RES[[files[i]]]$size_wav
        dat3[files[i], "channel_left_mean"] <- RES[[files[i]]]$channels["mean_left"]
        dat3[files[i], "channel_right_mean"] <- RES[[files[i]]]$channels["mean_right"]
        dat3[files[i], "channel_left_sd"] <- RES[[files[i]]]$channels["sd_left"]
        dat3[files[i], "channel_right_sd"] <- RES[[files[i]]]$channels["sd_right"]
    }
    cat("Finished at", date(), "\n\n")
}
DAT <- dat3
LOC <- aru_parse_location(rownames(dat3)[1])
#csv_out <- gsub("\\", "+", DIRIN, fixed=TRUE)
#csv_out <- gsub(":", "", csv_out, fixed=TRUE)
write.csv(DAT, row.names=FALSE, file=file.path(TMPDIR, paste0(LOC, ".csv")))
save(DIRIN, TMPDIR, CPROG, KEEP_WAV, CONVERT, NMAX, RES, DAT, LOC,
    file=file.path(TMPDIR, paste0(LOC, ".RData")))
pdf(file.path(TMPDIR, paste0(LOC, ".pdf")), width=10, height=10)
aru_plot(DAT)
dev.off()


do_all <- function(DIRIN) {
    dat1 <- aru_parse_directory(DIRIN)
    dat2 <- aru_flag_files(dat1, nmax=NMAX)
    dat3 <- dat2
    dat3$duration <- NA
    dat3$size_wav <- NA
    dat3$channel_left_mean <- NA
    dat3$channel_right_mean <- NA
    dat3$channel_left_sd <- NA
    dat3$channel_right_sd <- NA

    if (CONVERT) {
        files <- rownames(dat2)[dat2$to_check != 0]
        RES <- list()
        cat("Directory:", DIRIN, "\nStarted at", date(), "\n")
        for (i in seq_len(length(files))) {
            cat(i, "/", length(files), "\n")
            flush.console()
            #OBJECT <- dat2[files[i], "to_check"] > 0
            OBJECT <- FALSE
            RES[[files[i]]] <- aru_summarize(
                FNAME = files[i],
                DIRIN = DIRIN,
                DIROUT = TMPDIR,
                CPROG = CPROG,
                keep_wav = KEEP_WAV,
                quick=TRUE,
                return_object=OBJECT)
            dat3[files[i], "duration"] <- RES[[files[i]]]$duration
            dat3[files[i], "size_wav"] <- RES[[files[i]]]$size_wav
            dat3[files[i], "channel_left_mean"] <- RES[[files[i]]]$channels["mean_left"]
            dat3[files[i], "channel_right_mean"] <- RES[[files[i]]]$channels["mean_right"]
            dat3[files[i], "channel_left_sd"] <- RES[[files[i]]]$channels["sd_left"]
            dat3[files[i], "channel_right_sd"] <- RES[[files[i]]]$channels["sd_right"]
        }
        cat("Finished at", date(), "\n\n")
    }
    DAT <- dat3
    LOC <- aru_parse_location(rownames(dat3)[1])
    #csv_out <- gsub("\\", "+", DIRIN, fixed=TRUE)
    #csv_out <- gsub(":", "", csv_out, fixed=TRUE)
    write.csv(DAT, row.names=FALSE, file=file.path(TMPDIR, paste0(LOC, ".csv")))
    save(DIRIN, TMPDIR, CPROG, KEEP_WAV, CONVERT, NMAX, RES, DAT, LOC,
        file=file.path(TMPDIR, paste0(LOC, ".RData")))
    pdf(file.path(TMPDIR, paste0(LOC, ".pdf")), width=10, height=10)
    on.exit(dev.off())
    aru_plot(DAT)
    invisible(NULL)
}

load("e:\\peter\\AB_data_v2016\\out\\aru-filtering\\dirs-files.Rdata")
for (i in 1:length(dirs)) {
    cat("\n\n--->>>", i, "/", length(dirs), "<<<---\n")
    try(do_all(dirs[i]))
}
