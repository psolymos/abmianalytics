ROOT   <- "e:/peter/AB_data_v2017/data/raw/species"
OUTDIR <- "e:/peter/AB_data_v2017/data/analysis/species"
library(mefa4)

x <- read.csv(file.path(ROOT, "vplant-loweg-20170404.csv"))
#levels(x$YEAR) <- gsub(",", "", levels(x$YEAR))
#x$YEAR <- as.integer(as.character(x$YEAR))

table(x$TLVC_COVER_PCT)
x$TLVC_COVER_PCT[x$TLVC_COVER_PCT %in% c("DNC", "PNA", "VNA", "SNI")] <- NA
x$TLVC_COVER_PCT <- droplevels(x$TLVC_COVER_PCT)
levels(x$TLVC_COVER_PCT)[levels(x$TLVC_COVER_PCT) == "<1"] <- "1"
x$TLVC_COVER_PCT <- as.integer(as.character(x$TLVC_COVER_PCT))

slugify <- function(x, ...) {
    if (!is.factor(x))
        x <- as.character(x)
    l0 <- if (is.factor(x))
        levels(x) else unique(x)
    l <- l0

    l <- gsub("[[:punct:]]", "", l)
    l <- gsub("[[:space:]]", "-", l)
    l <- tolower(l)

    if (is.factor(x)) {
        levels(x) <- l
    } else {
        x <- l[match(x, l0)]
    }
    x
}

#x$Species <- slugify(x$SCIENTIFIC_NAME)

x$SPECIES_OLD <- x$SCIENTIFIC_NAME
levels(x$SCIENTIFIC_NAME) <- nameAlnum(levels(x$SCIENTIFIC_NAME), capitalize="mixed", collapse="")
x$SiteYear <- interaction(x$SITE, x$YEAR, x$TLVC_PLOT, drop=TRUE, sep="_")

y <- Xtab(TLVC_COVER_PCT ~ SiteYear + SCIENTIFIC_NAME, x,
    cdrop=c("NONE", "DNC", "PNA", "VNA", "SNI"))

cn1 <- c("SiteYear", "ROTATION", "SITE", "YEAR", "FIELDDATE", "FIELDCREW", "TLVC_ID_DATE",
    "TLVC_ID_CREW", "TLVC_PLOT")
cn2 <- c("SCIENTIFIC_NAME", "COMMON_NAME", "SPECIES_OLD", "RANK_NAME", "TSNID")
s <- nonDuplicated(x[,cn1], SiteYear, TRUE)
t <- nonDuplicated(x[,cn2], SCIENTIFIC_NAME, TRUE)
m <- Mefa(y, s, t, "inner")

df1 <- data.frame(samp(m), as.matrix(xtab(m)))
df2 <- taxa(m)
write.csv(df1, row.names=FALSE, file=file.path(OUTDIR,
    paste0("LowVegCover_out_", as.Date(Sys.time()), ".csv")))
write.csv(df2, row.names=FALSE, file=file.path(OUTDIR,
    paste0("LowVegCover_spp_", as.Date(Sys.time()), ".csv")))
