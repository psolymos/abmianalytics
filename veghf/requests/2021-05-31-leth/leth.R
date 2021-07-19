sites <- c("OG-ABMI-1600-71-1", "OG-ABMI-1600-71-4", "OG-ABMI-1635-71-1",
    "OG-ABMI-1600-71-11", "OG-ABMI-1600-71-12", "OG-ABMI-1600-71-7",
    "OG-ABMI-1633-71-7", "OG-ABMI-1633-71-6", "OG-ABMI-1633-71-5",
    "OG-ABMI-1633-71-4", "OG-ABMI-1633-71-1", "OG-ABMI-1617-72-42",
    "OG-ABMI-1617-72-40", "OG-ABMI-1617-72-31", "OG-ABMI-1617-72-20",
    "OG-ABMI-1617-72-10", "OG-ABMI-1600-72-68", "OG-ABMI-1600-72-52",
    "OG-ABMI-1600-72-4", "OG-ABMI-1600-72-34", "OG-ABMI-1600-72-33",
    "OG-ABMI-1600-72-12", "OG-ABMI-1600-71-10", "OG-ABMI-1617-72-57",
    "OG-ABMI-1600-72-27", "OG-ABMI-1600-72-26", "OG-ABMI-1600-72-19")


library(mefa4)
load("d:/abmi/AB_data_v2020/data/analysis/species/birds/ab-birds-all-2020-09-23.RData")


s <- dd[dd$PCODE %in% c("ABMIRF", "ABMISM"),]

z <- sites[1]
sum(grepl("OG", rownames(s)))

sapply(sites, function(z) sum(grepl(z, rownames(s))))


h <- rownames(dd)[grepl("OG-ABMI", rownames(dd))]
table(substr(h, 1, 8))
h1 <- h[substr(h, 1, 8) == "OG-ABMI-"]
h2 <- h[substr(h, 1, 8) == "OG-ABMI:"]

hh1 <- sapply(strsplit(h1, "_"), "[[", 1)
hh2 <- strsplit(h2, "_")
hh2 <- strsplit(sapply(hh2, "[[", 1), "::")
hh2 <- sapply(hh2, function(z) paste0(z[1], "-", z[3]))
hh <- h
names(hh) <- h
hh[substr(h, 1, 8) == "OG-ABMI-"] <- hh1
hh[substr(h, 1, 8) == "OG-ABMI:"] <- hh2

sites2 <- unique(sapply(strsplit(sites, "-"), function(z) paste0(z[1:4], collapse="-")))


sum(grepl(sites2[1], hh))
sum(grepl(sites2[2], hh))
sum(grepl(sites2[3], hh))
sum(grepl(sites2[4], hh))
sum(grepl(sites2[5], hh))

