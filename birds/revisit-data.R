library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")

load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData")
## raw RF and SM data
load(file.path("d:/abmi/AB_data_v2018", "data", "inter", "species", "birds-revisit.Rdata"))
ls()


ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds"
ee <- new.env()
load(file.path(ROOT, "ab-birds-all-2018-11-29.RData"), envir=ee)
TAX <- ee$tax
rm(ee)

en <- new.env()
#load(file.path(ROOT, "data", "ab-birds-north-2018-12-07.RData"), envir=en)
load(file.path(ROOT, "data", "ab-birds-north-2019-07-29-reference.RData"), envir=en)
es <- new.env()
#load(file.path(ROOT, "data", "ab-birds-south-2018-12-07.RData"), envir=es)
load(file.path(ROOT, "data", "ab-birds-south-2019-07-29-reference.RData"), envir=es)

Xn <- get_model_matrix(en$DAT, en$mods)
Xs <- get_model_matrix(es$DAT, es$mods)

Xn <- cbind(Xn, "vegcCrop"=0, "vegcIndustrial"=0, "vegcMine"=0, "vegcRoughP"=0, "vegcRural"=0,
    "vegcTameP"=0, "vegcUrban"=0)

Xage <- as.matrix(read.csv("~/repos/abmianalytics/lookup/Xn-veg-v61.csv"))
colnames(Xage) <- colnames(Xn)[match(colnames(Xage), make.names(colnames(Xn)))]

#keep <- en$DAT$PCODE %in% c("ABMIRF", "ABMISM")
keep <- en$DAT$PCODE %in% c("ABMIRF")
table(keep)

d <- droplevels(en$DAT[keep,])
tmp <- strsplit(as.character(d$SS), "_")
d$ABMIsite <- sapply(tmp, "[[", 1)
d$ABMIbirdpt <- sapply(tmp, "[[", 2)
d$ABMIsiteYear <- paste0(d$ABMIsite, "_", d$YEAR)
rev <- table(d$ABMIsite, d$YEAR)
rev[rev > 0] <- 1
str(d)
table(rowSums(rev), d$NRNAME[match(names(rowSums(rev)), d$ABMIsite)])

Dat <- nonDuplicated(d, d$ABMIsiteYear, TRUE)[,c("PCODE",
    "SS", "SSYR", "PKEY", "YEAR", "DATE", "DATI",
    "X", "Y", "NRNAME", "NSRNAME", "LUF_NAME",
    "ABMIsite", "ABMIbirdpt", "ABMIsiteYear")]

y <- en$YY[keep,]
y[y > 0] <- 1
y <- y[,colSums(y) > 0]

Y <- groupSums(cbind(nBirdPoints=1, y), 1, d$ABMIsiteYear)
range(Y)
table(Y[,"nBirdPoints"])
Dat$nBirdPoints <- Y[,1]
Y <- as.matrix(Y[,-1])[rownames(Dat),]
D <- Y
D[] <- 0

#spp <- "ALFL"
for (spp in colnames(D)) {
    cat(spp, "\n");flush.console()
    resn <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    if (!is.null(resn)) {
        lam <- rowMeans(exp(predict_with_SSH(resn, Xn[keep,], stage="Space")))
        tmp <- sum_by(lam, d$ABMIsiteYear)[rownames(Dat),]
        D[,spp] <- tmp[,"x"] / tmp[,"by"]
    }
}

colnames(Y) <- colnames(D) <- as.character(TAX[colnames(Y), "sppid"])

z <- nonDuplicated(resRF, resRF$site_year, TRUE)[rownames(Dat),c(
    "IdentificationAnalyst", "WindConditions", "Precipitation")]
Dat <- data.frame(Dat, z)

v <- read.csv("s:/reports/2018/data/species-info.csv")
rownames(v) <- v$SpeciesID
sp <- intersect(v$SpeciesID[v$model_north], colnames(Y))
str(sp)

Y <- Y[,sp]
D <- D[,sp]

range(D)
for (i in colnames(D)) {
    q <- quantile(D[,i], 0.999)
    D[D[,i] > q,i] <- q
}
range(D)
apply(D, 2, max)

save(Dat, Y, D, file="RFdata-toDave-20190724.RData", version=2)

write.csv(D, file="reference-bird-densities-20190730.csv")

