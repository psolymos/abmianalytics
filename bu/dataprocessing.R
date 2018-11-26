library(odbc)
library(DBI)
source("~/.ssh/boreal")

con <- dbConnect(odbc::odbc(),
    dsn = "BOREAL",
    database = "EMCLA_SQL_Database",
    driver = "SQL Server",
    uid = ..boreal_db_access$uid,
    pwd = ..boreal_db_access$pwd)

dbl <- dbListTables(con)

d0 <- dbReadTable(con, "ViewSpecies_LongFormCount")

dbDisconnect(con)

for (i in 1:ncol(d0))
    if (is.character(d0[,i]))
        d0[,i] <- as.factor(d0[,i])

d <- droplevels(d0[d0$Method %in% as.character(c(0, 8, 11:14)) & d0$Replicate == 1,])
d$MaxDur <- 3
d$MaxDur[d$Method %in% c("12", "13")] <- 1
d$MaxDur[d$Method == "0"] <- 10

## what is UID

d$SS <- with(d, interaction(
    ProjectID,
    Cluster,
    SITE,
    STATION,
    sep="::", drop=TRUE))
tmp1 <- as.character(d$RECORDING_DATE)
tmp2 <- sapply(strsplit(as.character(d$RECORD_TIME), " "), "[[", 2)
d$yr <- as.numeric(substr(tmp1, 1, 4))
d$dati <- as.POSIXlt(paste(tmp1, tmp2))
d$PKEY <- as.factor(paste0(
    as.character(d$SS),
    "_",
    tmp1,
    "-",
    tmp2))

load(file.path("d:/abmi/AB_data_v2018", "data", "analysis", "site",
    "veg-hf_BAM-BBS-BU_v6verified.Rdata"))

mefa4::compare_sets(d$SS, dd_point$SS)

