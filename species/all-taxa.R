library(mefa4)
source("~/repos/rabmi/R/rabmi.R")

ROOT <- "e:/peter/AB_data_v2018"
DATE <- as.Date(Sys.time(), tz=Sys.timezone(location = TRUE))
gis <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")

## vascular plants -------------------------

taxon <- "vplants"


tab <- "T15"
sub_col <- "Quadrant"

res <- get_table(tab)
save(res, file=file.path(ROOT, "data", "raw", "species",
    paste0(taxon, "_", DATE, ".Rdata")))

colnames(res) <- gsub(" ", "", colnames(res))
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

levels(res$ScientificName) <- gsub("X ", "", levels(res$ScientificName))
levels(res$ScientificName) <- gsub(" x ", " ", levels(res$ScientificName))
levels(res$ScientificName) <- sapply(strsplit(levels(res$ScientificName), " "), function(z) {
    paste(z[1:min(2, length(z))], collapse=" ")
})

res$SpeciesID <- res$ScientificName
levels(res$SpeciesID) <- nameAlnum(levels(res$SpeciesID), capitalize="mixed", collapse="")
res$SpeciesID <- droplevels(res$SpeciesID)


cn1 <- c("ABMISite", "Year", "subunit", "site_year", "site_year_sub", "offgrid", "nearest")
cn2 <- c("SpeciesID", "CommonName", "ScientificName", "TaxonomicResolution",
    "UniqueTaxonomicIdentificationNumber")
x_site_year <- nonDuplicated(res, res$site_year, TRUE)[,cn1]
x_site_year$subunit <- x_site_year$site_year_sub <- NULL
rownames(x_site_year) <- x_site_year$site_year
x_site_year_sub <- nonDuplicated(res, res$site_year_sub, TRUE)[,cn1]
rownames(x_site_year_sub) <- x_site_year_sub$site_year_sub

z <- nonDuplicated(res, res$SpeciesID, TRUE)[,c(cn2)]
rownames(z) <- z$SpeciesID

y_site_year_sub <- Xtab(~ site_year_sub + SpeciesID, res,
    cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"))
y_site_year_sub[y_site_year_sub > 0] <- 1

m_site_year_sub <- Mefa(y_site_year_sub, x_site_year_sub, z)
m_site_year_sub <- m_site_year_sub[,taxa(m_site_year_sub)$TaxonomicResolution %in%
    c("Genus", "Species")]

tmp <- m_site_year_sub[samp(m_site_year_sub)$subunit %in%
    c("NE", "NW", "SE", "SW")]
nn <- sum_by(rep(1, nrow(tmp)), droplevels(samp(tmp)$site_year))
y_site_year <- groupSums(xtab(tmp), 1, droplevels(samp(tmp)$site_year))
stopifnot(max(y_site_year) == 4)

m_site_year <- Mefa(y_site_year, x_site_year, z)
samp(m_site_year)$nQuadrant <- nn[match(rownames(m_site_year), rownames(nn)),"by"]

save(m_site_year, m_site_year_sub, file=file.path(ROOT, "data", "inter", "species",
    paste0(taxon, "_", DATE, ".Rdata")))


out_site_year <- data.frame(samp(m_site_year), as.matrix(xtab(m_site_year)))
out_site_year_sub <- data.frame(samp(m_site_year_sub), as.matrix(xtab(m_site_year_sub)))
out_species <- taxa(m_site_year)

write.csv(out_site_year, row.names=FALSE,
    file=file.path(ROOT, "data", "analysis", "species",
    paste0(taxon, "_SiteBinom_", DATE, ".csv")))

write.csv(out_site_year_sub, row.names=FALSE,
    file=file.path(ROOT, "data", "analysis", "species",
    paste0(taxon, "_Quadrant_", DATE, ".csv")))

write.csvout_species(, row.names=FALSE,
    file=file.path(ROOT, "data", "analysis", "species",
    paste0(taxon, "_Species_", DATE, ".csv")))

save(out_site_year, out_site_year_sub, out_species,
    file=file.path(ROOT, "data", "analysis", "species",
    paste0(taxon, "_out_", DATE, ".Rdata")))

