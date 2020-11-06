## organizing files for sc-dev website and BB
library(mefa4)
library(jsonlite)

Taxa <- c("lichens", "mites", "mosses", "vplants", "birds",
    "mammals", "habitats")
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
sapply(XT, function(z) range(rowSums(z)))

addmargins(sapply(XT, function(z) colSums(z)))

toJSON(head(XT$lichens), pretty=TRUE)

# append linear image to coef
if (FALSE) {
library(magick)
for (taxon in names(XT)) {
    for (spp in rownames(XT[[taxon]])) {
        cat(taxon, spp, "\n")
        flush.console()
        xt <- XT[[taxon]][spp,]
        if (xt$veghf) {
            i1 <- image_append(image_read(
                paste0("s:/AB_data_v2020/Results/web/", taxon, "/", spp,
                    c("/veghf.png", "/lin-north.png"))))
            image_write(i1,
                paste0("s:/AB_data_v2020/Results/web/", taxon, "/", spp,
                    "/coef-north.png"))
        }
        if (xt$soilhf) {
            i1 <- image_append(image_read(
                paste0("s:/AB_data_v2020/Results/web/", taxon, "/", spp,
                    c("/soilhf.png", "/lin-south.png"))))
            image_write(i1,
                paste0("s:/AB_data_v2020/Results/web/", taxon, "/", spp,
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
            name="Vasculer plants"
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
        )
    )
)

toJSON(xx, auto_unbox=TRUE, pretty=TRUE)
writeLines(
    toJSON(x, auto_unbox=TRUE),
    "s:/AB_data_v2020/Results/web/index.json")


load("s:/AB_data_v2020/Results/COEFS-ALL.RData")
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
COEFS$mammals <- COEFS2$mammals
COEFS$habitats <- COEFS2$habitats


#taxon <- "habitats"
for (taxon in names(COEFS)) {
    cat(taxon, "\n")


    s0 <- COEFS[[taxon]]$species

    if (!(taxon %in% c("habitats", "birds", "mammals"))) {
        s <- s0
        s <- data.frame(
                id=as.character(s$SpeciesID),
                scientific=as.character(s[,grep("Scient", colnames(s))]),
                common=NA_character_,
                stringsAsFactors = FALSE)
        #s$id[grep("\\.\\.", s$id)] <- gsub("\\.\\.", "\\.", s$id[grep("\\.\\.", s$id)])
        s$id[endsWith(s$id, ".")] <- substr(s$id[endsWith(s$id, ".")], 1,
            nchar(s$id[endsWith(s$id, ".")])-1)
    }
    if (taxon=="habitats") {
        s <- s0
        s <- data.frame(
            id=as.character(s$sppid),
            scientific=as.character(s$ScientificName),
            common=as.character(s$CommonName),
            stringsAsFactors = FALSE)
    }
    if (taxon=="mammals") {
        s <- s0
        s <- data.frame(
            id=as.character(s$SpeciesID),
            scientific=as.character(s$ScientificName),
            common=as.character(s$CommonName),
            stringsAsFactors = FALSE)
    }
    if (taxon=="birds") {
        ROOT <- "d:/abmi/AB_data_v2020/data/analysis/species/birds"
        ee <- new.env()
        load(file.path(ROOT, "ab-birds-all-2020-09-23.RData"), envir=ee)
        TAX <- ee$tax
        rm(ee)
        s <- TAX
        s <- data.frame(
            id=as.character(s$sppid),
            scientific=as.character(s$scinam),
            common=as.character(s$species),
            stringsAsFactors = FALSE)
    }

    rownames(s) <- s$id
    if (all(is.na(s$common))) {
        s$display <- s$scientific
    } else {
        s$display <- ifelse(is.na(s$scientific), s$common,
            paste0(s$common, " (", s$scientific, ")"))
    }

    #"Do.not.analyze" to exclude
    x <- XT[[taxon]]
    x <- x[rownames(x) != "Do.not.analyze", ]
    compare_sets(rownames(s), rownames(x))
    setdiff(rownames(s), rownames(x))
    setdiff(rownames(x), rownames(s))

    s <- s[rownames(x),]

    s$det <- x$det
    s$useavailnorth <- x$`useavail-north`
    s$useavailsouth <- x$`useavail-south`
    s$coefnorth <- x$veghf
    s$coefsouth <- x$soilhf
    s$map <- x$map
    s$sectornorth <- x$`sector-north`
    s$sectorsouth <- x$`sector-south`

    s <- s[order(s$display),]
    s$idprev <- c(s$id[length(s$id)], s$id[-length(s$id)])
    s$idnext <- c(s$id[-1], s$id[1])

    o <- xx$taxa[sapply(xx$taxa, "[[", "id") == taxon][[1]]
    o$species <- s

    #toJSON(o, auto_unbox=TRUE, pretty=TRUE)
    writeLines(
        toJSON(o, auto_unbox=TRUE),
        paste0("s:/AB_data_v2020/Results/web/", taxon, "/index.json"))

    for (spp in rownames(s)) {

        oo <- as.list(s[spp,])
        oo$taxonid <- xx$taxa[sapply(xx$taxa, "[[", "id") == taxon][[1]]$id
        oo$taxonname <- xx$taxa[sapply(xx$taxa, "[[", "id") == taxon][[1]]$name

        toJSON(oo, auto_unbox=TRUE, pretty=TRUE)
        writeLines(
            toJSON(oo, auto_unbox=TRUE),
            paste0("s:/AB_data_v2020/Results/web/", taxon, "/", spp, "/index.json"))

    }

}
