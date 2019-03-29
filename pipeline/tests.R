library(mefa4)

rm(list=ls())

test_names <- function(x, ref, MARGIN, details=FALSE) {
    if (details) {
        compare_sets(dimnames(ref)[[MARGIN]], dimnames(x)[[MARGIN]])
    } else {
        all(dimnames(ref)[[MARGIN]] %in% dimnames(x)[[MARGIN]])
    }
}


f0 <- "d:/abmi/sppweb2018/c4i/tables/StandardizedOutput.RData"
tax <- read.csv("d:/abmi/AB_data_v2018/data/raw/species/taxonomy.csv")

tx <- "vplants"
#tx <- "mites"
#tx <- "mosses"
#tx <- "lichens"

if (tx == "vplants") {
    fn <- paste0("s:/Result from Ermias_2018/vplants/north/ABMI-vplant-north-results-2019-03-20.RData")
    fs <- paste0("s:/Result from Ermias_2018/vplants/south/ABMI-vplant-south-results-2019-03-20.RData")
    fl <- paste0("s:/Result from Ermias_2018/vplants/Species lookup for Plants 2018.RData")
}
if (tx == "mites") {
    fn <- paste0("s:/Result from Ermias_2018/mites/north/ABMI-mites-north-results-2019-03-20.RData")
    fs <- paste0("s:/Result from Ermias_2018/mites/south/ABMI-mites-south-results-2019-03-20.RData")
    fl <- paste0("s:/Result from Ermias_2018/mites/Species lookup for Mites 2018.RData")
}
if (tx == "mosses") {
    fn <- paste0("s:/Result from Ermias_2018/mosses/north/ABMI-moss-north-results-2019-03-20.RData")
    fs <- paste0("s:/Result from Ermias_2018/mosses/south/ABMI-moss-south-results-2019-03-20.RData")
    fl <- paste0("s:/Result from Ermias_2018/mosses/Species lookup for Moss 2018.RData")
}
if (tx == "lichens") {
    fn <- paste0("s:/Result from Ermias_2018/lichens/north/ABMI-lichen-north-results-2019-03-20.RData")
    fs <- paste0("s:/Result from Ermias_2018/lichens/south/ABMI-lichen-south-results-2019-03-20.RData")
    fl <- paste0("s:/Result from Ermias_2018/lichens/Species lookup for Lichens 2018.RData")
}

stopifnot(file.exists(fn))
stopifnot(file.exists(fs))
stopifnot(file.exists(fl))

e0 <- new.env()
en <- new.env()
es <- new.env()
el <- new.env()
load(f0, envir=e0)
load(fn, envir=en)
load(fs, envir=es)
load(fl, envir=el)

sort(names(e0))
sort(names(en))
sort(names(es))
sort(names(el))

## Lookup

test_names(el$Lookup, e0$Lookup, 2, details=TRUE)


el$Lookup$LinkHabitat <- "binomial_logit"
el$Lookup$LinkSpclim <- "binomial_logit"
el$Lookup$AUCNorth <- en$ModEval.North[rownames(el$Lookup), "AUC.All"]
el$Lookup$AUCSouth <- es$ModEval.South[rownames(el$Lookup), "AUC.All"]
if (tax == "vplants") {
    el$Lookup$Nonnative[el$Lookup$Nonnative == "Native/Exotic"] <- "TRUE"
    el$Lookup$Nonnative <- as.logical(el$Lookup$Nonnative)
} else {
    el$Lookup$Nonnative <- FALSE
}

for (i in colnames(e0$Lookup))
    if (!(i %in% colnames(el$Lookup)))
        el$Lookup[[i]] <- NA
Lookup <- el$Lookup[,colnames(e0$Lookup)]
rownames(Lookup) <- Lookup$SpeciesID
Lookup$CommonName <- tax$COMMON_NAME[match(Lookup$ScientificName, tax$SCIENTIFIC_NAME)]
Lookup$CommonName[Lookup$CommonName == ""] <- NA
Lookup$CommonName <- droplevels(Lookup$CommonName)
Lookup$UseavailNorth <- rownames(Lookup) %in% rownames(en$UseavailNorth)# & !Lookup$ModelNorth
Lookup$UseavailSouth <- rownames(Lookup) %in% rownames(es$UseavailSouth)# & !Lookup$ModelSouth
str(Lookup)
summary(Lookup)

## UseavailSouth

UseavailSouth <- array(NA,
    c(nrow(es$UseavailSouth), ncol(e0$UseavailSouth)),
    dimnames=list(rownames(es$UseavailSouth), colnames(e0$UseavailSouth)))
for (i in colnames(UseavailSouth))
    if (i %in% colnames(es$UseavailSouth))
        UseavailSouth[,i] <- es$UseavailSouth[,i]

## UseavailNorth

UseavailNorth <- array(NA,
    c(nrow(en$UseavailNorth), ncol(e0$UseavailNorth)),
    dimnames=list(rownames(en$UseavailNorth), colnames(e0$UseavailNorth)))
for (i in colnames(UseavailNorth))
    if (i %in% colnames(en$UseavailNorth))
        UseavailNorth[,i] <- en$UseavailNorth[,i]

## CoefSouth + CI

CoefSouth <- array(NA,
    c(nrow(es$CoefSouth), ncol(e0$CoefSouth), dim(es$CoefSouth.bs)[3]+1),
    dimnames=list(rownames(es$CoefSouth), colnames(e0$CoefSouth), NULL))
test_names(es$CoefSouth, e0$CoefSouth, 2, details=TRUE)

cn <- intersect(colnames(es$CoefSouth), colnames(e0$CoefSouth))
CoefSouth[,cn,1] <- es$CoefSouth[dimnames(CoefSouth)[[1]],cn]
CoefSouth[,cn,-1] <- es$CoefSouth.bs[dimnames(CoefSouth)[[1]],cn,]
Lookup$ModelSouth <- rownames(Lookup) %in% rownames(CoefSouth)
LowerSouth <- CoefSouth[,,1]
LowerSouth[] <- NA
LowerSouth[,cn] <- es$LowerSouth[dimnames(CoefSouth)[[1]],cn]
UpperSouth <- CoefSouth[,,1]
UpperSouth[] <- NA
UpperSouth[,cn] <- es$UpperSouth[dimnames(CoefSouth)[[1]],cn]

## CoefNorth + CI

CoefNorth <- array(NA,
    c(nrow(en$CoefNorth), ncol(e0$CoefNorth), dim(en$CoefNorth.bs)[3]+1),
    dimnames=list(rownames(en$CoefNorth), colnames(e0$CoefNorth), NULL))
test_names(en$CoefNorth, e0$CoefNorth, 2, details=TRUE)

cn <- intersect(colnames(en$CoefNorth), colnames(e0$CoefNorth))
CoefNorth[,cn,1] <- en$CoefNorth[dimnames(CoefNorth)[[1]],cn]
CoefNorth[,cn,-1] <- en$CoefNorth.bs[dimnames(CoefNorth)[[1]],cn,]
Lookup$ModelNorth <- rownames(Lookup) %in% rownames(CoefNorth)
LowerNorth <- CoefNorth[,,1]
LowerNorth[] <- NA
LowerNorth[,cn] <- en$LowerNorth[dimnames(CoefNorth)[[1]],cn]
UpperNorth <- CoefNorth[,,1]
UpperNorth[] <- NA
UpperNorth[,cn] <- en$UpperNorth[dimnames(CoefNorth)[[1]],cn]

## LinearSouth

LinearSouth <- as.matrix(es$LinearHF.10[rownames(CoefSouth),
    c("AverageCoef", "SoftLin10", "HardLin10")])

## LinearNorth

LinearNorth <- as.matrix(en$LinearHF.10[rownames(CoefNorth),
    c("AverageCoef", "SoftLin10", "HardLin10")])

## SpclimSouth

SpclimSouth <- array(NA,
    c(nrow(es$SpclimSouth), ncol(e0$SpclimSouth), dim(es$SpclimSouth.bs)[3]+1),
    dimnames=list(rownames(es$SpclimSouth), colnames(e0$SpclimSouth), NULL))
test_names(es$SpclimSouth, e0$SpclimSouth, 2, details=TRUE)

cn <- intersect(colnames(es$SpclimSouth), colnames(e0$SpclimSouth))
SpclimSouth[,cn,1] <- as.matrix(es$SpclimSouth[dimnames(SpclimSouth)[[1]],cn])
SpclimSouth[,cn,-1] <- es$SpclimSouth.bs[dimnames(SpclimSouth)[[1]],cn,]

## SpclimNorth

SpclimNorth <- array(NA,
    c(nrow(en$SpclimNorth), ncol(e0$SpclimNorth), dim(en$SpclimNorth.bs)[3]+1),
    dimnames=list(rownames(en$SpclimNorth), colnames(e0$SpclimNorth), NULL))
test_names(en$SpclimNorth, e0$SpclimNorth, 2, details=TRUE)

cn <- intersect(colnames(en$SpclimNorth), colnames(e0$SpclimNorth))
SpclimNorth[,cn,1] <- as.matrix(en$SpclimNorth[dimnames(SpclimNorth)[[1]],cn])
SpclimNorth[,cn,-1] <- en$SpclimNorth.bs[dimnames(SpclimNorth)[[1]],cn,]


toSave <- c("Lookup",
    "CoefNorth", "CoefSouth",
    "SpclimNorth", "SpclimSouth",
    "UseavailNorth", "UseavailSouth",
    "LinearNorth", "LinearSouth",
    "UpperNorth", "UpperSouth",
    "LowerNorth", "LowerSouth")
for (i in toSave) {
    cat("\n\n-------------\n", i, "\n\n")
    print(str(get(i)))
}

table(rowSums(Lookup[,c("ModelNorth", "ModelSouth", "UseavailNorth", "UseavailSouth")]))
colSums(Lookup[,c("ModelNorth", "ModelSouth", "UseavailNorth", "UseavailSouth")])
table(N=Lookup$ModelNorth, S=Lookup$ModelSouth)

save(list=toSave,
    file=paste0("d:/abmi/sppweb2018/c4i/tables/StandardizedOutput-", tx, ".RData"))

