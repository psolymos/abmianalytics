## read in species
## O - N+S averaged
## N & S: no averaging just veg or soil
## birds: p-transformation
## store results in matrix that is saved

## preamble
#library(mefa4)
taxa <- c("mammals", "birds", "lichens", "mites", "mosses", "vplants")
CH <- c("S1", "O1", "O2", "O3", paste0("N", 1:11))

#ch="O1"
#tax="mammals"
for (ch in CH) {
    for (tax in taxa) {
        cat(ch, tax, "\n")
        flush.console()
        gc()
        load(paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, ".RData")) # OUT
        if (tax == "mammals") {
            sr <- OUT[,1:length(taxa)]
            colnames(sr) <- taxa
            sr[] <- 0
        }
        sr[,tax] <- rowSums(OUT)
        rm(OUT)
        gc()
    }
    save(sr, file=paste0("s:/AB_data_v2019/bdqt/final/richness/sr-", ch, ".RData"))
}

## NN and OF

library(cure4insect)
set_options(path = "d:/abmi/reports")
load_common_data()
ST <- get_species_table()
ST2 <- read.csv("~/repos/abmispecies/_data/birds.csv")
nnp <- rownames(ST)[!ST$native & ST$taxon == "vplants"]
ofb <- as.character(ST2$sppid[ST2$oldforest == 1])

for (ch in CH[1:14]) {
    cat(ch, "\n")
    flush.console()
    load(paste0("s:/AB_data_v2019/bdqt/inter/", "vplants", "/", ch, ".RData")) # OUT
    nn <- rowSums(OUT[,colnames(OUT) %in% nnp])
    rm(OUT)
    gc()
    load(paste0("s:/AB_data_v2019/bdqt/inter/", "birds", "/", ch, ".RData")) # OUT
    of <- rowSums(OUT[,colnames(OUT) %in% ofb])
    rm(OUT)
    gc()
    save(nn, of, file=paste0("s:/AB_data_v2019/bdqt/final/richness/nn-of-", ch, ".RData"))
}

## stitching the pieces

taxa <- c("mammals", "birds", "lichens", "mites", "mosses", "vplants")
CH <- c("S1", "O1", "O2", "O3", paste0("N", 1:11))

load("s:/AB_data_v2019/bdqt/bdqt-poly-hab_2019-12-10.RData")
SR <- matrix(0, nrow(x), length(taxa)+3)
rownames(SR) <- x$UID
rm(x)
colnames(SR) <- c(taxa, "total", "nnplants", "ofbirds")

for (ch in CH) {
    gc()
    cat(ch, "\n")
    flush.console()
    load(paste0("s:/AB_data_v2019/bdqt/final/richness/sr-", ch, ".RData"))
    SR[rownames(sr),taxa] <- sr
}
summary(SR[,taxa])

for (ch in CH) {
    gc()
    cat(ch, "\n")
    flush.console()
    load(paste0("s:/AB_data_v2019/bdqt/final/richness/nn-of-", ch, ".RData"))
    SR[names(nn),"nnplants"] <- nn
    SR[names(of),"ofbirds"] <- of
}

## standardizing richness to 0-1
for (tax in taxa) {
    cat(tax, "\n")
    flush.console()
    M <- max(SR[,tax])
    SR[,tax] <- SR[,tax] / M
}
SR[,"total"] <- rowMeans(SR[,taxa])
M <- max(SR[,"total"])
SR[,"total"] <- SR[,"total"] / M
summary(SR)

M <- max(SR[,"nnplants"])
SR[,"nnplants"] <- SR[,"nnplants"] / M
M <- max(SR[,"ofbirds"])
SR[,"ofbirds"] <- SR[,"ofbirds"] / M
summary(SR)

save(SR, file="s:/AB_data_v2019/bdqt/bdqt-richness_2019-12-18.RData")

## uploading layers to DB

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

load("s:/AB_data_v2019/bdqt/bdqt-richness_2019-12-18.RData")



dbDisconnect(con)

