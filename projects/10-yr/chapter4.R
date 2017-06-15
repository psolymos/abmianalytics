library(mefa4)
fl <- c(
aqplants="e:/peter/AB_data_v2017/data/analysis/species/OUT_aqplants_2017-04-07.Rdata",
birdsrf="e:/peter/AB_data_v2017/data/analysis/species/OUT_birdsrf_2017-04-26.Rdata",
lichens="e:/peter/AB_data_v2017/data/analysis/species/OUT_lichens_2017-05-29.Rdata",
mites="e:/peter/AB_data_v2017/data/analysis/species/OUT_mites_2017-04-27.Rdata",
vplants="e:/peter/AB_data_v2017/data/analysis/species/OUT_vplants_2017-04-13.Rdata",
mosses="e:/peter/AB_data_v2017/data/analysis/species/OUT_mosses_2017-04-05.Rdata")

u <- list()
for (i in 1:length(fl)) {
e <- new.env()
load(fl[i], envir=e)
u[[names(fl)[i]]] <- list(nsite=nrow(e$res2), nspp=ncol(e$m), siteyr=table(e$res2$YEAR))
}

e <- new.env()
load("e:/peter/AB_data_v2017/data/analysis/species/OUT_birdsrf_2017-04-26.Rdata", envir=e)
load("e:/peter/AB_data_v2017/data/analysis/species/OUT_birdssm_2017-04-05.Rdata")

length(union(colnames(e$m), colnames(xt)))
length(setdiff(colnames(xt), colnames(e$m)))
setdiff(colnames(xt), colnames(e$m))
ncol(xt)

xx <- nonDuplicated(x, SITE, TRUE)
table(xx$YEAR)

m <- read.csv("e:/peter/AB_data_v2015/out/species/OUT_Mammals_Species_Transect-Binomial-Length-DSS_2015-06-01.csv")
table(m$OnOffGrid, m$Year)