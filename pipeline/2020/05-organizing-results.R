## organizing files for sc-dev website and BB
library(mefa4)
library(jsonlite)

Taxa <- c("lichens", "mites", "mosses", "vplants", "birds",
    "mammals", "habitats", "nnplants")
names(Taxa) <- Taxa

ROOT <- "s:/AB_data_v2020/Results"
labs <- c("det", "useavail-north", "useavail-south",
    "veghf", "lin-north",
    "soilhf", "lin-south",
    "map", "sector-north", "sector-south")

RES <- list()
for (taxon in Taxa) {
    SPP <- list.dirs(file.path(ROOT, "web", taxon),
        full.names = FALSE, recursive = FALSE)
    TAB <- NULL
    for (spp in SPP) {
        X <- gsub("\\.png", "", list.files(file.path(ROOT, "web", taxon, spp)))
        D <- data.frame(
            taxon=factor(taxon, Taxa),
            species=factor(spp, SPP),
            plot=factor(X[X %in% labs], labs))
        TAB <- rbind(TAB, D)
    }
    RES[[taxon]] <- TAB
}


XT <- lapply(Taxa, function(i) {
    y <- as.matrix(Xtab(~ species + plot, RES[[i]]))
    y[y[,"veghf"]>0,"useavail-north"] <- 0
    y[y[,"soilhf"]>0,"useavail-south"] <- 0
    q <- y > 0
    as.data.frame(q)
})
XT$habitats$`sector-north` <- FALSE
XT$habitats$`sector-south` <- FALSE
XT$nnplants$`sector-north` <- FALSE
XT$nnplants$`sector-south` <- FALSE
sapply(XT, function(z) range(rowSums(z)))

addmargins(sapply(XT, function(z) colSums(z)))

toJSON(head(XT$lichens), pretty=TRUE)

# append linear image to coef
if (FALSE) {
library(magick)
RT <- "s:"
#RT <- "d:/abmi"
for (taxon in names(XT)) {
    for (spp in rownames(XT[[taxon]])) {
        gc(full=TRUE)
        cat(taxon, spp, "\n")
        flush.console()
        xt <- XT[[taxon]][spp,]
        if (xt$veghf) {
            i1 <- image_append(image_read(
                paste0(RT, "/AB_data_v2020/Results/web/", taxon, "/", spp,
                    c("/veghf.png", "/lin-north.png"))))
            image_write(i1,
                paste0(RT, "/AB_data_v2020/Results/web/", taxon, "/", spp,
                    "/coef-north.png"))
        }
        if (xt$soilhf) {
            i1 <- image_append(image_read(
                paste0(RT, "/AB_data_v2020/Results/web/", taxon, "/", spp,
                    c("/soilhf.png", "/lin-south.png"))))
            image_write(i1,
                paste0(RT, "/AB_data_v2020/Results/web/", taxon, "/", spp,
                    "/coef-south.png"))
        }
    }
}
}

## create JSON objects

xx <- list(
    creator="ABMI Science Centre",
    url="https://science.abmi.ca",
    version="2020",
    taxa=list(
        list(
            id="lichens",
            name="Lichens"
        ),
        list(
            id="mosses",
            name="Bryophytes"
        ),
        list(
            id="vplants",
            name="Vascular plants"
        ),
        list(
            id="mites",
            name="Soil mites"
        ),
        list(
            id="birds",
            name="Birds"
        ),
        list(
            id="mammals",
            name="Mammals (camera)"
        ),
        list(
            id="habitats",
            name="Habitat elements"
        ),
        list(
            id="nnplants",
            name="Non-native plant richness"
        )
    )
)

toJSON(xx, auto_unbox=TRUE, pretty=TRUE)
#writeLines(
#    toJSON(xx, auto_unbox=TRUE),
#    "s:/AB_data_v2020/Results/web/index.json")


load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
COEFS$mammals <- COEFS2$mammals
COEFS$habitats <- COEFS2$habitats
COEFS$nnplants <- COEFS2$nnplants

COEFS$mammals$species["Mountaingoat", "CommonName"] <- "Mountaingoat"
COEFS$mammals$species["RiverOtter", "CommonName"] <- "River Otter"
COEFS$vplants$species["Petunia.x.atkinsiana","ScientificName"] <- "Petunia x atkinsiana"

#taxon <- "habitats"
RES <- NULL
ckeep <- c("id", "scientific", "common", "display",
        "det", "useavailnorth",
        "useavailsouth", "coefnorth", "coefsouth", "map", "sectornorth",
        "sectorsouth", "idprev", "idnext")
for (taxon in names(COEFS)) {
    cat("\n\n\n", taxon, "---------------\n\n")

    s <- COEFS[[taxon]]$species
    s$SpeciesID[endsWith(s$SpeciesID, ".")] <- substr(
        s$SpeciesID[endsWith(s$SpeciesID, ".")], 1,
            nchar(s$SpeciesID[endsWith(s$SpeciesID, ".")])-1)
    s$id <- s$SpeciesID
    s$scientific <- s$ScientificName
    s$common <- s$CommonName

    rownames(s) <- s$id
    if (all(is.na(s$common))) {
        s$display <- s$scientific
    } else {
        s$display <- ifelse(is.na(s$scientific), s$common,
            paste0(s$common, " (", s$scientific, ")"))
    }
    if (taxon=="vplants") {
        s$display <- ifelse(is.na(s$common), s$scientific,
            paste0(s$scientific, " - ", s$common))
    }

    #"Do.not.analyze" to exclude
    x <- XT[[taxon]]
    #x <- x[rownames(x) != "Do.not.analyze", ]
    print(compare_sets(rownames(s), rownames(x)))
    print(setdiff(rownames(s), rownames(x)))
    print(setdiff(rownames(x), rownames(s)))

    s <- s[rownames(x),]

    s$det <- x$det
    s$useavailnorth <- x$`useavail-north`
    s$useavailsouth <- x$`useavail-south`
    s$coefnorth <- x$veghf
    s$coefsouth <- x$soilhf
    s$map <- x$map
    s$sectornorth <- x$`sector-north`
    s$sectorsouth <- x$`sector-south`

    s <- s[order(as.character(s$display)),]
    s$idprev <- c(s$id[length(s$id)], s$id[-length(s$id)])
    s$idnext <- c(s$id[-1], s$id[1])

    o <- xx$taxa[sapply(xx$taxa, "[[", "id") == taxon][[1]]
    o$species <- s[,ckeep]

if (FALSE) {
    #toJSON(o, auto_unbox=TRUE, pretty=TRUE)
    writeLines(
        toJSON(o, auto_unbox=TRUE),
        paste0("s:/AB_data_v2020/Results/web/", taxon, "/index.json"))

    for (spp in rownames(s)) {

        oo <- as.list(s[spp,ckeep])
        oo$taxonid <- xx$taxa[sapply(xx$taxa, "[[", "id") == taxon][[1]]$id
        oo$taxonname <- xx$taxa[sapply(xx$taxa, "[[", "id") == taxon][[1]]$name

        toJSON(oo, auto_unbox=TRUE, pretty=TRUE)
        writeLines(
            toJSON(oo, auto_unbox=TRUE),
            paste0("s:/AB_data_v2020/Results/web/", taxon, "/", spp, "/index.json"))

    }
}
    sss <- s
    sss$taxonid <- o$id
    sss$taxonname <- o$name
    RES <- rbind(RES, sss)


}

ok <- RES$ModelNorth | RES$ModelSouth
data.frame(table(RES$Group[ok]))
any(is.na(RES$display))

VER <- read.csv("pipeline/2020/VER.csv")

cn <- c("id", "scientific", "common", "display",  "taxonid", "taxonname")
xxx <- xx
xxx$info <- VER
xxx$species <- RES[,cn]

writeLines(
    toJSON(xxx, auto_unbox=TRUE),
    "s:/AB_data_v2020/Results/web/index.json")



## BioBrowser updates

library(openxlsx)

RT <- "s:/AB_data_v2020/Results"
load(file.path(RT, "BB-ESTIMATES.RData")) # "RESULTS" "SPECIES"
load(file.path(RT, "UseAvail-all.RData")) # UA
SEff <- list()
load(file.path(RT, "SEffect-birds.RData")) # SE
#tmp <- rownames(RES)[RES$Group=="birds"]
#table(names(SE) %in% tmp)
#SE <- SE[names(SE) %in% tmp]
#save(SE, file=file.path(RT, "SEffect-birds.RData"))
SEff$birds <- SE
load(file.path(RT, "SEffect-habitats.RData"))
SEff$habitats <- SE
load(file.path(RT, "SEffect-lichens.RData"))
SEff$lichens <- SE
load(file.path(RT, "SEffect-mammals.RData"))
SEff$mammals <- SE
load(file.path(RT, "SEffect-mites.RData"))
SEff$mites <- SE
load(file.path(RT, "SEffect-mosses.RData"))
SEff$mosses <- SE
load(file.path(RT, "SEffect-nnplants.RData"))
SEff$nnplants <- SE
load(file.path(RT, "SEffect-vplants.RData"))
SEff$vplants <- SE
SE <- SEff
rm(SEff)
cn0 <- c("SpeciesID", "ScientificName",  "CommonName", "TSNID", "Group")

# useavail north
v <- do.call(rbind, lapply(UA, function(z) z$north))
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
p <- RES[RES$useavailnorth, cn0]
compare_sets(rownames(p), rownames(v))
UAn <- data.frame(p, v[rownames(p),])

# useavail south
v <- do.call(rbind, lapply(UA, function(z) z$south))
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
p <- RES[RES$useavailsouth, cn0]
compare_sets(rownames(p), rownames(v))
UAs <- data.frame(p, v[rownames(p),])

fse <- function(z, north=TRUE) {
    zz <- if (north)
        z$north else z$south
    if (is.null(zz))
        return(NULL)
    c(area=zz$unit["area",],
        regional=zz$regional, underhf=zz$underhf, unit=zz$unit["unit",])
}
# sector north
v <- do.call(rbind, lapply(SE, function(tx) do.call(rbind, lapply(tx, fse, north=TRUE))))
p <- RES[RES$sectornorth, cn0]
compare_sets(rownames(p), rownames(v))
SEn <- data.frame(p, v[rownames(p),])

# sector south
v <- do.call(rbind, lapply(SE, function(tx) do.call(rbind, lapply(tx, fse, north=FALSE))))
p <- RES[RES$sectorsouth, cn0]
compare_sets(rownames(p), rownames(v))
SEs <- data.frame(p, v[rownames(p),])

# veghf north
v <- do.call(rbind, lapply(RESULTS$veghf, function(z) {
    c0 <- z$Estimate
    c1 <- z$LCL
    c2 <- z$UCL
    names(c0) <- rownames(z)
    names(c1) <- paste0(rownames(z), ".LCL")
    names(c2) <- paste0(rownames(z), ".UCL")
    c(c0, c1, c2)
}))
p <- RES[RES$coefnorth, cn0]
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
compare_sets(rownames(p), rownames(v))
CFn <- data.frame(p, v[rownames(p),])

# veghf north
v <- do.call(rbind, lapply(RESULTS$linn, function(z) {
    o <- c(z$First, z$Lower, z$Upper)
    names(o) <- c(rownames(z), paste0(rownames(z), ".LCL"), paste0(rownames(z), ".UCL"))
    o
}))
p <- RES[RES$coefnorth, cn0]
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
compare_sets(rownames(p), rownames(v))
CFLn <- data.frame(p, v[rownames(p),])

# veghf south
v <- do.call(rbind, lapply(RESULTS$lins, function(z) {
    o <- c(z$First, z$Lower, z$Upper)
    names(o) <- c(rownames(z), paste0(rownames(z), ".LCL"), paste0(rownames(z), ".UCL"))
    o
}))
p <- RES[RES$coefsouth, cn0]
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
compare_sets(rownames(p), rownames(v))
CFLs <- data.frame(p, v[rownames(p),])

# soilhfnontreed
v <- do.call(rbind, lapply(RESULTS$soilhfnontreed, function(z) {
    c0 <- z$Estimate
    c1 <- z$LCL
    c2 <- z$UCL
    names(c0) <- rownames(z)
    names(c1) <- paste0(rownames(z), ".LCL")
    names(c2) <- paste0(rownames(z), ".UCL")
    c(c0, c1, c2)
}))
p <- RES[RES$coefsouth, cn0]
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
compare_sets(rownames(p), rownames(v))
CFNTs <- data.frame(p, v[rownames(p),])

# soilhftreed
v <- do.call(rbind, lapply(RESULTS$soilhftreed, function(z) {
    c0 <- z$Estimate
    c1 <- z$LCL
    c2 <- z$UCL
    names(c0) <- rownames(z)
    names(c1) <- paste0(rownames(z), ".LCL")
    names(c2) <- paste0(rownames(z), ".UCL")
    c(c0, c1, c2)
}))
p <- RES[RES$coefsouth, cn0]
rownames(v)[endsWith(rownames(v), ".")] <-
    substr(rownames(v)[endsWith(rownames(v), ".")], 1, nchar(rownames(v)[endsWith(rownames(v), ".")])-1)
compare_sets(rownames(p), rownames(v))
CFTs <- data.frame(p, v[rownames(p),])

z <- RESULTS$veghf[[1]]
tmp <- c(cn0, rownames(z), paste0(rownames(z), ".LCL"), paste0(rownames(z), ".UCL"))
colnames(CFn) <- tmp

z <- RESULTS$soilhfnontreed[[1]]
tmp <- c(cn0, rownames(z), paste0(rownames(z), ".LCL"), paste0(rownames(z), ".UCL"))
colnames(CFNTs) <- colnames(CFTs) <- tmp

meta <- read.csv("pipeline/2020/meta.csv")


fcheck <- function(zz) {
    zzz <- zz[,6:ncol(zz)]
    print(zzz[rowSums(is.na(zzz))>0,])
    cat("\n")
    print(range(zzz, na.rm=TRUE))
    invisible()
}
fcheck(UAn)
fcheck(UAs)
fcheck(CFn)
fcheck(CFLn)
fcheck(CFNTs)
fcheck(CFTs)
fcheck(CFLs)
fcheck(SEn)
fcheck(SEs)

ffix <- function(zz) {
    zzz <- zz[,!(colnames(zz) %in% cn0)]
    zzz[is.na(zzz)] <- 0
    zzz[zzz > 10^4] <- 10^4
    data.frame(zz[,cn0], zzz)
}
fcheck(ffix(UAn))
fcheck(ffix(UAs))
fcheck(ffix(CFn))
fcheck(ffix(CFLn))
fcheck(ffix(CFNTs))
fcheck(ffix(CFTs))
fcheck(ffix(CFLs))
fcheck(ffix(SEn))
fcheck(ffix(SEs))

LIST <- list(
    Info=VER,
    Species=RES[,c("SpeciesID", "ScientificName",  "CommonName", "TSNID", "Group",
        "det", "useavailnorth", "useavailsouth",
        "coefnorth", "coefsouth", "map", "sectornorth", "sectorsouth",
        "Occurrences", "nSites", "Nonnative", "LinkHabitat",
        "LinkSpclim", "AUCNorth", "AUCSouth", "Comments")],
    UseavailNorth=ffix(UAn),
    UseavailSouth=ffix(UAs),
    VeghfNorth=ffix(CFn),
    LinearNorth=ffix(CFLn),
    SoilhfSouthNontreed=ffix(CFNTs),
    SoilhfSouthTreed=ffix(CFTs),
    LinearSouth=ffix(CFLs),
    SectorNorth=ffix(SEn),
    SectorSouth=ffix(SEs),
    Metadata=meta)


write.xlsx(LIST, file.path(RT, "DataPortalUpdate_2021-01-18.xlsx"))

## csv version
lf <- list.files(file.path(RT, "normalized-maps"), recursive = TRUE, full.names = TRUE)
lf2 <- lf
lf2 <- gsub("normalized-maps", "normalized-maps-csv", lf2)
lf2 <- gsub("RData", "csv", lf2)
if (FALSE) {
for (i in 1:length(lf2)) {
    cat(i, "\n")
    flush.console()
    local({
        load(lf[i])
        out <- data.frame(LinkID=rownames(out), out)
        write.csv(out, row.names = FALSE, file=lf2[i])
    })
}
}

