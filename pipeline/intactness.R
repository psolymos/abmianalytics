library(mefa4)

load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
b <- read.csv("d:/abmi/sppweb2018/c4i/tables/StandardizedOutput-birds-final-lookup.csv")


#tx <- "birds"
allSI <- list()
for (tx in c("lichens", "mites", "mosses", "vplants", "birds")) {
    rt <- paste0("s:/reports/2018/results/", tx, "/boot/")
    fl <- list.files(rt)
    if (tx == "birds")
        fl <- paste0(as.character(b[b$ModelNorth | b$ModelSouth,"SpeciesID"]), ".RData")

    res <- list()
    for (i in fl) {
        cat(tx, i, "\n")
        flush.console()

        e <- new.env()
        load(paste0(rt, i), envir=e)

        NC <- e$Curr.Boot[match(kgrid$Row10_Col10, rownames(e$Curr.Boot)),]
        NC[is.na(NC)] <- 0
        NR <- e$Ref.Boot[match(kgrid$Row10_Col10, rownames(e$Ref.Boot)),]
        NR[is.na(NR)] <- 0

        NC <- groupSums(NC, 1, kgrid$LUF_NAME)
        NR <- groupSums(NR, 1, kgrid$LUF_NAME)

        SI <- NC / NR
        SI[is.na(SI)] <- 1
        SI[SI > 1] <- (NR / NC)[SI > 1]
        SI <- 100*SI

        res[[gsub(".RData", "", i)]] <- list(NC=NC, NR=NR, SI=SI)
    }
    allSI[[tx]] <- res
}

## evaluate percent abundance in LUF region to screen out spp not there
## use the subset and average SI over those
## save allSI and taxa results
## calculate mean + CI and present some histograms/vases?

save(allSI, file="d:/abmi/AB_data_v2018/data/analysis/intactness-by-LUF.Rdata")

load("d:/abmi/AB_data_v2018/data/analysis/intactness-by-LUF.Rdata")

TX <- c("lichens", "mites", "mosses", "vplants", "birds")


#tx <- "birds"
SI <- NULL
for (tx in TX) {

    a <- array(NA, c(7, 100, length(allSI[[tx]])))
    dimnames(a) <- list(rownames(allSI[[1]][[1]]$SI), NULL, names(allSI[[tx]]))

    for (i in names(allSI[[tx]])) {
        tmp <- allSI[[tx]][[i]]
        nc <- tmp$NC[,1]
        nc <- nc*100/sum(nc)
        si <- tmp$SI
        si[nc < 1,] <- NA
        a[,,i] <- si
    }

    z <- sapply(dimnames(a)[[1]], function(i) rowMeans(a[i,,], na.rm=TRUE))
    rownames(z) <- 1:100
    z <- Melt(z)
    colnames(z) <- c("iter", "LUF_Region", "Intactness")
    SI <- rbind(SI, data.frame(Taxonomic_Group=tx, z))
}
str(SI)

library(ggplot2)

p <- ggplot(SI, aes(x=Taxonomic_Group, y=Intactness)) +
    geom_boxplot(aes(fill=Taxonomic_Group)) +
    facet_grid(. ~ LUF_Region) + ylim(0, 100) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


aa <- aggregate(SI$Intactness, list(LUF=SI$LUF_Region, Taxon=SI$Taxonomic_Group), quantile, c(0.5, 0.05, 0.95))
write.csv(aa, row.names=FALSE, file="intactness-by-LUF.csv")


