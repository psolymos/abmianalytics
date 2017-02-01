ROOT <- "e:/peter/AB_data_v2017"
library(mefa4)

x <- read.csv(file.path(ROOT, "data", "raw", "species", "LowVegCover-PSolymos.csv"))
levels(x$YEAR) <- gsub(",", "", levels(x$YEAR))
x$YEAR <- as.integer(as.character(x$YEAR))

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

x$Species <- slugify(x$SCIENTIFIC_NAME)

y <- Xtab(TLVC_COVER_PCT ~ SITE_LABEL + Species, x,
    cdrop=tolower(c("NONE", "DNC", "PNA", "VNA", "SNI")))
s <- nonDuplicated(x[,1:9], SITE_LABEL, TRUE)
t <- nonDuplicated(x[,c(15, 10:13)], Species, TRUE)
m <- Mefa(y, s, t, "inner")

df1 <- data.frame(samp(m), as.matrix(xtab(m)))
df2 <- taxa(m)
write.csv(df1, row.names=FALSE, file=file.path(ROOT, "data", "inter", "species",
    paste0("LowVegCover_out_", as.Date(Sys.time()), ".csv")))
write.csv(df2, row.names=FALSE, file=file.path(ROOT, "data", "inter", "species",
    paste0("LowVegCover_spp_", as.Date(Sys.time()), ".csv")))
