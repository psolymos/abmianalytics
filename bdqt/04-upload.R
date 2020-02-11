library(cure4insect)
set_options(path="d:/abmi/reports")
load_common_data()

library(DBI)
source("~/.ssh/postgres")
data.frame(..postgres_science)

con <- dbConnect(
    odbc::odbc(),
    driver   = "PostgreSQL Unicode(x64)",
    server   = ..postgres_science$host,
    database = ..postgres_science$database,
    uid      = ..postgres_science$username,
    pwd      = ..postgres_science$password,
    port     = ..postgres_science$port)
(dbl <- dbListTables(con))

#tmp <- data.frame(x=1:3, y=5:7)
#dbWriteTable(con, "test", tmp)
#dbSendQuery(con, "drop table test")

#d0 <- dbReadTable(con, "map_YellowheadedBlackbird")
#dbSendQuery(con, "DROP TABLE map_YellowheadedBlackbird;")

#dbl <- dbl[grepl("map_", dbl)]
#for (i in dbl) {
#    dbSendQuery(con, paste0("DROP TABLE \"", i, "\";"))
#}

dbDisconnect(db)

SP <- get_species_table()
XY <- get_id_table()


dbWriteTable(con, "species_map_id", XY)
dbWriteTable(con, "species_map_lookup", SP)

for (spp in rownames(tax)[48:173]) {

    cat(spp, "\n");flush.console()
#    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-04-01/", spp, ".RData"))
    load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-05-14/", spp, ".RData"))

    TYPE <- "C" # combo
    if (tax[spp, "ModelSouth"] && !tax[spp, "ModelNorth"])
        TYPE <- "S"
    if (!tax[spp, "ModelSouth"] && tax[spp, "ModelNorth"])
        TYPE <- "N"

    for (i in 1:ncol(CURRB)) {
        q <- quantile(CURRB[,i], 0.99, na.rm=TRUE)
        CURRB[CURRB[,i] > q,i] <- q
        q <- quantile(REFB[,i], 0.99, na.rm=TRUE)
        REFB[REFB[,i] > q,i] <- q
    }
    Dcr <- rowMeans(CURRB, na.rm=TRUE)
    Drf <- rowMeans(REFB, na.rm=TRUE)
    SD <- apply(CURRB, 1, sd, na.rm=TRUE)

#    Dcr <- rowSums(Curr)
#    q <- quantile(Dcr, 0.99)
#    Dcr[Dcr > q] <- q
#    summary(Dcr)
#    Drf <- rowSums(Ref)
#    q <- quantile(Drf, 0.99)
#    Drf[Drf > q] <- q
#    summary(Drf)
    MAX <- max(Dcr, Drf)

    df <- (Dcr-Drf) / MAX
    df <- sign(df) * abs(df)^0.5
    df <- pmin(200, ceiling(99 * df)+100)
    df[df==0] <- 1
    cr <- pmin(100, ceiling(99 * sqrt(Dcr / MAX))+1)
    rf <- pmin(100, ceiling(99 * sqrt(Drf / MAX))+1)
    crsd <- 100 * SD / mean(Dcr)
    si <- 100 * pmin(Dcr, Drf)/pmax(Dcr, Drf)
    si[is.na(si)] <- 100
    si[si==0] <- 1

    if (SAVE) {
        tmp <- data.frame(
            ID=rownames(kgrid),
            Current=Dcr,
            Reference=Drf)
            #Color_Current=col1[cr],
            #Color_Reference=col1[rf],
            #Color_Difference=col3[df]

        if (TYPE == "S")
            tmp <- tmp[kgrid$useS,]
        if (TYPE == "N")
            tmp <- tmp[kgrid$useN,]
        dbWriteTable(con, "test_num", tmp,
            overwrite=TRUE, row.names=FALSE)
    }

