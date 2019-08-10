library(mefa4)
source("~/repos/abmianalytics/species/abmi-r-api.R")
#data.frame(table=get_table_names())

## settings
TAXA <- c("vplants", "mites")
ROOT <- "d:/abmi/AB_data_v2019"

## common stuff
DATE <- as.Date(Sys.time(), tz=Sys.timezone(location = TRUE))
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
add_labels <- function(res, sub_col) {
    res$offgrid <- startsWith(as.character(res$ABMISite), "OG")
    res$subunit <- res[[sub_col]]
    res[[sub_col]] <- NULL
    res$site_year <- interaction(res$ABMISite, res$Year, drop=TRUE, sep="_")
    res$site_year_sub <- interaction(res$ABMISite, res$Year, res$subunit, drop=TRUE, sep="_")

    tmp <- strsplit(as.character(res$ABMISite), "-")
    res$nearest <- sapply(tmp, function(z) {
        zz <- if (length(z) > 1) z[3] else z[1]
        as.integer(gsub("\\D+", "", zz))
    })
    res
}
normalize_species <- function(res, spgen_only=TRUE) {
    if (spgen_only) {
#        res$ScientificName0 <- res$ScientificName
        levels(res$ScientificName) <- gsub("X ", "", levels(res$ScientificName))
        levels(res$ScientificName) <- gsub(" x ", " ", levels(res$ScientificName))
        levels(res$ScientificName) <- sapply(strsplit(levels(res$ScientificName), " "), function(z) {
            paste(z[1:min(2, length(z))], collapse=" ")
        })

        levels(res$TaxonomicResolution)[levels(res$TaxonomicResolution) %in%
            c("Subspecies", "Variety")] <- "Species"
    }
    res$SpeciesID <- res$ScientificName
    levels(res$SpeciesID) <- nameAlnum(levels(res$SpeciesID), capitalize="mixed", collapse="")
    res$SpeciesID <- droplevels(res$SpeciesID)

    res
}
cn1 <- c("ABMISite", "Year", "subunit", "site_year", "site_year_sub", "offgrid", "nearest")
cn2 <- c("SpeciesID", "CommonName", "ScientificName", "TaxonomicResolution",
    "UniqueTaxonomicIdentificationNumber")

for (taxon in TAXA) {
    cat("taxon:", taxon, "\n  - pulling and normalizing data")
    flush.console()

    ## vascular plants -------------------------
    if (taxon == "vplants") {
        tab <- "T15"
        sub_col <- "Quadrant"
        allowed_subunits <- c("NE", "NW", "SE", "SW")
        allowed_resolution <- c("Genus", "Species")
        sub_max <- 4

        res <- get_table(tab)
        res0 <- res
        save_list <- "res0"

        colnames(res) <- gsub(" ", "", colnames(res))
        res <- add_labels(res, sub_col=sub_col)
        res <- normalize_species(res)
    }

    ## mosses -------------------------
    if (taxon == "mosses") {
        tab1 <- "T19A" # T19A Moss Identification (2003-2008)
        tab2 <- "T19B" # T19B Moss Identification (since 2009)
        sub_col <- "Quadrant"
        allowed_subunits <- c("NE", "NW", "SE", "SW", "1ha")
        allowed_resolution <- c("Genus", "Species")
        sub_max <- 4

        res1 <- get_table(tab1)
        res2 <- get_table(tab2)
        res01 <- res1
        res02 <- res2
        save_list <- c("res01", "res02")

        colnames(res1) <- gsub(" ", "", colnames(res1))
        res1[[sub_col]] <- as.factor("1ha")
        res1 <- add_labels(res1, sub_col=sub_col)
        res1 <- normalize_species(res1)

        colnames(res2) <- gsub(" ", "", colnames(res2))
        res2 <- add_labels(res2, sub_col=sub_col)
        res2 <- normalize_species(res2)

        tmp <- intersect(colnames(res1), colnames(res2))
        res <- rbind(res1[,tmp], res2[,tmp])
    }

    ## mites -------------------------
    if (taxon == "mites") {
        tab <- "T24A"
        sub_col <- "Quadrant"
        allowed_subunits <- c("NE", "NW", "SE", "SW")
        allowed_resolution <- c("Genus", "Species", "Subspecies")
        sub_max <- 4

        res <- get_table(tab)
        res0 <- res
        save_list <- "res0"

        colnames(res) <- gsub(" ", "", colnames(res))
        res <- add_labels(res, sub_col=sub_col)
        res <- normalize_species(res, spgen_only=FALSE) # keep spp names as is
        res$CommonName <- NA
    }

    cat(" --- OK\n  - processing attributes and x-tabs")
    flush.console()
    ## sample attributes
    x_site_year <- nonDuplicated(res, res$site_year, TRUE)[,cn1]
    x_site_year$subunit <- x_site_year$site_year_sub <- NULL
    rownames(x_site_year) <- x_site_year$site_year
    x_site_year_sub <- nonDuplicated(res, res$site_year_sub, TRUE)[,cn1]
    rownames(x_site_year_sub) <- x_site_year_sub$site_year_sub

    ## species attributes
    z <- nonDuplicated(res, res$SpeciesID, TRUE)[,c(cn2)]
    rownames(z) <- z$SpeciesID

    ## sample-species cross tabs
    y_site_year_sub <- Xtab(~ site_year_sub + SpeciesID, res,
        cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"))
    y_site_year_sub01 <- y_site_year_sub
    y_site_year_sub01[y_site_year_sub01 > 0] <- 1
    if (taxon %in% c("vplants", "mosses", "lichens"))
        y_site_year_sub <- y_site_year_sub01

    ## mefa bundles for sample/subunits
    m_site_year_sub <- Mefa(y_site_year_sub, x_site_year_sub, z)
    m_site_year_sub <- m_site_year_sub[,taxa(m_site_year_sub)$TaxonomicResolution %in%
        allowed_resolution]
    m_site_year_sub01 <- Mefa(y_site_year_sub01, x_site_year_sub, z)
    m_site_year_sub01 <- m_site_year_sub01[,taxa(m_site_year_sub)$TaxonomicResolution %in%
        allowed_resolution]

    ## aggregated cross tabs for binomial tables
    tmp <- m_site_year_sub01[samp(m_site_year_sub01)$subunit %in% allowed_subunits]
    nn <- sum_by(rep(1, nrow(tmp)), droplevels(samp(tmp)$site_year))
    y_site_year <- groupSums(xtab(tmp), 1, droplevels(samp(tmp)$site_year))
    stopifnot(max(y_site_year) == sub_max)

    ## aggregated mefa bundles for samples
    m_site_year <- Mefa(y_site_year, x_site_year, z)
    samp(m_site_year)$nQuadrant <- nn[match(rownames(m_site_year), rownames(nn)),"by"]

    ## csv view
    out_site_year <- data.frame(samp(m_site_year), as.matrix(xtab(m_site_year)))
    out_site_year_sub <- data.frame(samp(m_site_year_sub), as.matrix(xtab(m_site_year_sub)))
    out_species <- taxa(m_site_year)

    cat(" --- OK\n  - saving files")
    flush.console()

    ## save raw input
    save(list=save_list, file=file.path(ROOT, "data", "raw", "species",
        paste0(taxon, "_", DATE, ".Rdata")))

    ## save bundles
    save(m_site_year, m_site_year_sub, file=file.path(ROOT, "data", "inter", "species",
        paste0(taxon, "_", DATE, ".Rdata")))

    ## write csv & binary
    write.csv(out_site_year, row.names=FALSE,
        file=file.path(ROOT, "data", "analysis", "species",
        paste0(taxon, "_SiteBinom_", DATE, ".csv")))
    write.csv(out_site_year_sub, row.names=FALSE,
        file=file.path(ROOT, "data", "analysis", "species",
        paste0(taxon, "_Quadrant_", DATE, ".csv")))
    write.csv(out_species, row.names=FALSE,
        file=file.path(ROOT, "data", "analysis", "species",
        paste0(taxon, "_Species_", DATE, ".csv")))
    save(out_site_year, out_site_year_sub, out_species,
        file=file.path(ROOT, "data", "analysis", "species",
        paste0(taxon, "_out_", DATE, ".Rdata")))

    cat(" --- OK\n")
    flush.console()
}

if (FALSE) {

    tn <- c("T19A Moss Identification (2003-2008)",
        "T19B Moss Identification (since 2009)",
        "T20A Lichen Identification (2003-2008)",
        "T20B Lichen Identification (since 2009)",
        "T24A Soil Arthropods (Mites) Identification",
        "T15 Vascular Plants",

        "T26A Breeding Birds",
        "T26B Breeding Birds (ARUs)",
        "T26C ARU Deployment and Retrieval",
        "T26D Breeding Birds (ARU) Abiotic")
    names(tn) <- tn
    tn <- sapply(strsplit(tn, " "), "[[", 1)

    x <- list()
    for (tt in tn) {
        cat(tt, "\n");flush.console()
        try(x[[tt]] <- get_table(table = tt))
    }

    lapply(x, colnames)
    lapply(x, function(z) max(z$Year))

}
