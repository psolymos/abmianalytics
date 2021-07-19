## REVI tagging

## Find 125 presence + 125 absence
## filters: 3min, Maj 25 - June 25, ABMI or BU data

library(mefa4)
library(intrval)
load("d:/abmi/AB_data_v2020/data/analysis/species/birds/ab-birds-all-2020-09-23.RData")

LIMIT <- as.POSIXlt(c("2020-05-25", "2020-06-25"))$yday

#SPP <- c("REVI", "WTSP", "YEWA", "CCSP")

x <- dd[dd$PCODE == "ABMISM" & dd$MAXDUR == 3 & dd$DATI$yday %[]% LIMIT,]
dim(x)

y <- data.frame(as.matrix(yy[rownames(x), ]))
cs <- colSums(yy>0)
SPP <- colnames(y)[cs >= 125]
SPP <- SPP[!(SPP %in% "MOTR")]
y <- y[,SPP]

cn <- c("PCODE", "SS", "SSYR", "PKEY", "YEAR", "DATE", "DATI",
    "X", "Y", "NRNAME", "NSRNAME", "LUF_NAME")
x <- droplevels(x[,cn])
rownames(x) <- NULL

## p is presence (TRUE) or absence (FALSE)
rnd_fun <- function(p) {
    n0 <- sum(!p)
    n1 <- sum(p)
    v <- rep(NA_integer_, length(p))
    v[!p] <- -sample(n0, n0)
    v[p] <- sample(n1, n1)
    v
}

set.seed(1)
for (spp in SPP) {
    x[[paste0("count_", spp)]] <- as.numeric(y[,spp])
    x[[paste0("select_", spp)]] <- rnd_fun(y[,spp] > 0)
}

#x$REVI_select <- rnd_fun(x$REVI > 0)
#x$WTSP_select <- rnd_fun(x$WTSP > 0)
#x$YEWA_select <- rnd_fun(x$YEWA > 0)
#x$CCSP_select <- rnd_fun(x$CCSP > 0)

write.csv(x, row.names=FALSE, file="~/Desktop/bird-presence-absence-for-tagging-allspp.csv")

table(random_subset=abs(x$select_REVI) <= 125, presence_absence=x$count_REVI > 0)
