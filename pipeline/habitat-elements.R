## processing habitat element data to be displayed on bio browser
library(readxl)

ROOT <- "~/GoogleWork/abmi/habelem2020"

e <- list()
for (i in c("Species",                    # OK
            "VeghfNorth",                 #
            "LinearNorth",                #
            "SoilhfSouthNontreed",        # OK
            "SoilhfSouthTreed",           # OK
            "LinearSouth")                # OK
     ) {
    e[[i]] <- as.data.frame(
        read_xlsx(file.path(ROOT, "DataPortalUpdate_2019-04-10.xlsx"), i))
}
lapply(e, colnames)

## lookup

s <- read.csv(file.path(ROOT, "Species lookup for habitat June 2017.csv"))
s$Sp <- gsub("BA and volume summary ", "", as.character(s$Sp))

t1 <- data.frame(
    SpeciesID = s$Sp,
    ScientificName = NA,
    TSNID = NA,
    CommonName = s$SpName,
    ModelNorth = s$Analysis.North,
    ModelSouth = s$Analysis.South,
    UseavailNorth = FALSE,
    UseavailSouth = FALSE, SizeNorth=NA, SizeSouth=NA,Nonnative=FALSE,
    LinkHabitat=s$Link,
    LinkSpclim=s$Link,
    AUCNorth=NA, AUCSouth=NA,
    Comments="",
    Group="habitatelements",
    Unit=s$Units
)
rownames(t1) <- t1$SpeciesID

## south
es <- new.env()
load(file.path(ROOT,
    "Habitat coefficients South May 2020 OFFICIAL coefficients.Rdata"), envir=es)
names(es)
l10s <- read.csv(file.path(ROOT,
    "Linear 10pc figure values for habitat South.csv"))
pA <- read.csv(file.path(ROOT,
    "Habitat coefficients South May 2020 pApen coefficients.csv"))
pA <- structure(pA$pAspen, names=as.character(pA$Sp))
xs <- es$Coef.official
rownames(xs)[rownames(xs) == "CanopyCover"] <- "CanopyClosure"
xsl <- es$Coef.official.lci
rownames(xsl)[rownames(xsl) == "CanopyCover"] <- "CanopyClosure"
xsu <- es$Coef.official.uci
rownames(xsu)[rownames(xsu) == "CanopyCover"] <- "CanopyClosure"
names(pA)[names(pA) == "CanopyCover"] <- "CanopyClosure"
pA <- pA[rownames(xs)]
all(names(pA)==rownames(xs))

mefa4::compare_sets(rownames(t1), rownames(xs))
mefa4::compare_sets(colnames(e$SoilhfSouthNontreed), colnames(xs))
setdiff(colnames(xs), colnames(e$SoilhfSouthNontreed))
cn <- intersect(colnames(xs), colnames(e$SoilhfSouthNontreed))
tmp0 <- xs[,cn]
tmp1 <- xsl[,cn]
colnames(tmp1) <- paste0("Lower_", cn)
tmp2 <- xsu[,cn]
colnames(tmp2) <- paste0("Upper_", cn)
t4 <- data.frame(
    t1[rownames(xs),c("SpeciesID", "ScientificName",
                      "CommonName", "TSNID", "Group")],
    tmp0,
    tmp1,
    tmp2
)
LINK <- t1$LinkHabitat
levels(LINK) <- c("logit", "log", "log", "identity")
LINK <- as.character(LINK)
names(LINK) <- rownames(t1)

tmp0t <- tmp0
tmp1t <- tmp1
tmp2t <- tmp2
for (i in rownames(xs)) {
    f <- binomial(LINK[i])$linkfun
    fi <- binomial(LINK[i])$linkinv
    tmp0t[i,] <- fi(f(tmp0[i,]) + pA[i])
    tmp1t[i,] <- fi(f(tmp1[i,]) + pA[i])
    tmp2t[i,] <- fi(f(tmp2[i,]) + pA[i])
}

t5 <- data.frame(
    t1[rownames(xs),c("SpeciesID", "ScientificName",
                      "CommonName", "TSNID", "Group")],
    tmp0t,
    tmp1t,
    tmp2t
)

rownames(l10s) <- l10s$Element
rownames(l10s)[rownames(l10s) == "CanopyCover"] <- "CanopyClosure"
l10s2 <- l10s[rownames(xs),c("lin10.mean", "lin10.soft", "lin10.hard")]
colnames(l10s) <- c("AverageCoef", "SoftLin10", "HardLin10")
t6 <- data.frame(
    t1[rownames(xs),c("SpeciesID", "ScientificName",
                      "CommonName", "TSNID", "Group")],
    l10s2
)


## north

en <- new.env()
load(file.path(ROOT,
    "Habitat coefficients North May 2020 OFFICIAL coefficients ALL.Rdata"), envir=en)
names(en)
mefa4::compare_sets(colnames(en$Coef.official.ALL), colnames(e$VeghfNorth))
setdiff(colnames(en$Coef.official.ALL), colnames(e$VeghfNorth))

xn <- en$Coef.official.ALL
xnl <- en$Coef.official.ALL.lci
xnu <- en$Coef.official.ALL.uci

rownames(xn)[rownames(xn) == "CanopyCover"] <- "CanopyClosure"
rownames(xnl)[rownames(xnl) == "CanopyCover"] <- "CanopyClosure"
rownames(xnu)[rownames(xnu) == "CanopyCover"] <- "CanopyClosure"

rownames(xn)[rownames(xn) == "Live Decid BA"] <- "Live Deciduous BA"
rownames(xnl)[rownames(xnl) == "Live Decid BA"] <- "Live Deciduous BA"
rownames(xnu)[rownames(xnu) == "Live Decid BA"] <- "Live Deciduous BA"

mefa4::compare_sets(rownames(xn), s$Sp)

cn <- intersect(colnames(xn), colnames(e$VeghfNorth))
tmp0 <- xn[,cn]
tmp1 <- xnl[,cn]
colnames(tmp1) <- paste0("Lower_", cn)
tmp2 <- xnu[,cn]
colnames(tmp2) <- paste0("Upper_", cn)
t2 <- data.frame(
    t1[rownames(xn),c("SpeciesID", "ScientificName",
                      "CommonName", "TSNID", "Group")],
    tmp0,
    tmp1,
    tmp2
)

l10n <- read.csv(file.path(ROOT,
    "Linear 10pc figure values.csv"))
l10n2 <- read.csv(file.path(ROOT,
    "Linear 10pc figure values BA.csv"))
colnames(l10n2) <- colnames(l10n)
l10n <- rbind(l10n, l10n2)

rownames(l10n) <- l10n$Element
rownames(l10n)[rownames(l10n) == "CanopyCover"] <- "CanopyClosure"
rownames(l10n)[rownames(l10n) == "Live Decid BA"] <- "Live Deciduous BA"

stopifnot(all(rownames(l10n) == rownames(t2)))

l10n <- l10n[rownames(xn),c("lin10.mean", "lin10.soft", "lin10.hard")]
colnames(l10n) <- c("AverageCoef", "SoftLin10", "HardLin10")

t3 <- data.frame(
    t1[rownames(xn),c("SpeciesID", "ScientificName",
                      "CommonName", "TSNID", "Group")],
    l10n
)

library(openxlsx)

zn <- rownames(t1) %in% rownames(xn)
zs <- rownames(t1) %in% rownames(xs)
table(zn, t1$ModelNorth)
table(zs, t1$ModelSouth)

OUT <- list(t1, t2, t3, t4, t5, t6)
names(OUT) <- names(e)

write.xlsx(OUT, file=file.path(ROOT, "DataPortalUpdate_2020-09-14_habitatelements.xlsx"))



## writing csv files with current and reference abundances

fl <- list.files("d:/abmi/AB_data_v2020/misc/Habitat elements May 2020/Km2 summaries/")
a <- gsub(".csv", "",
    gsub("Km2 North reference and current May 2020 ", "", fl, fixed=TRUE), fixed=TRUE)
mefa4::compare_sets(a, rownames(t1))
setdiff(a, rownames(t1))
setdiff(rownames(t1), a)

for (spp in a) {

    cat(spp, "\n");flush.console()


    fin <- paste0("d:/abmi/AB_data_v2020/misc/Habitat elements May 2020/Km2 summaries/",
        "Km2 North reference and current May 2020 ", spp, ".csv")

    fout <- paste0("d:/abmi/AB_data_v2020/misc/Habitat elements May 2020/_dataportal/normalized_maps/", spp, ".csv")
    y <- read.csv(fin)
    print(round(range(y$Curr), 4))
    if (spp == "pH") {
        MIN <- min(y$Ref, y$Curr)
        y$Ref <- y$Ref+abs(MIN)
        y$Curr <- y$Curr+abs(MIN)
    }

    Curr <- y$Curr
    Ref <- y$Ref
    Dcr <- Curr
    q <- quantile(Dcr, 0.99)
    Dcr[Dcr > q] <- q
    Drf <- Ref
    q <- quantile(Drf, 0.99)
    Drf[Drf > q] <- q

    d <- data.frame(ID=y$LinkID,
        Current=round(Dcr, 6), Reference=round(Drf, 6))
    d$Current[d$Current < 10^-6] <- 0
    d$Reference[d$Reference < 10^-6] <- 0
    write.csv(d, row.names=FALSE, file=fout)

}



library(mefa4)
library(sf)

ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds" # change this bit
load(file.path(ROOT, "data", "ab-birds-north-2019-01-30.RData"))

bound <- st_read("s:/GC_eric/FromEric/Reporting_Boundaries/20191029_Boundaries_2019.gdb",
    'FMA_Current_ALPAC')

d <- DAT[,c("PKEY","SS","X", "Y", "vegc", "vegw")]
a <- c(0, 10, 20, 40, 60, 80, 100, 120, Inf)
d$Age <- DAT$wtAge*200
d$agec <- cut(d$Age, a, include.lowest=TRUE)
table(d$agec, useNA="a")
levels(d$agec) <- paste0(a[-length(a)], "-", a[-1])
d$agec <- as.character(d$agec)
d$agec[!(d$vegc %in% c("Decid", "BSpr", "Mixedwood", "Pine", "Spruce"))] <- ""
table(d$agec, useNA="a")

d$stand_age <- paste0(ifelse(DAT$fCC2 != 0, "CC", ""), as.character(d$vegc), "_", d$agec)
table(d$stand_age)


d <- st_as_sf(d, coords=c("X", "Y"), crs=4326)
d <- st_transform(d, st_crs(bound))

plot(d$geometry)
plot(bound$SHAPE, add=TRUE, col=2)

dd <- st_intersection(d, bound)
dd$Age <- round(dd$Age)

plot(dd$geometry)
plot(bound$SHAPE, add=TRUE, col=2)

table(dd$vegc, dd$agec)


dd <- as.data.frame(dd[,1:8])
dd$geometry <- NULL
save(dd, file="alpac-bird-points.RData")
