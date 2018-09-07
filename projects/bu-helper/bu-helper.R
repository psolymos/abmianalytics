library(mefa4)
library(jsonlite)
library(cure4insect)
set_options(verbose=0, path = "w:/reports")
load_common_data()

Species <- get_all_species("birds")
TAX <- droplevels(cure4insect:::.c4if$SP[Species,])
df <- cure4insect:::.c4if$KT
Area <- 1 # keep it as density

lh <- read.csv("~/GoogleWork/bam/lhreg2/bna-life-history-info.csv")
lh$id <- lh$species
levels(lh$id) <- nameAlnum(levels(lh$id), capitalize="mixed", collapse="")
TAX$Hab <- lh$habitat[match(rownames(TAX), lh$id)]
tmp <- c("BairdsSparrow"="Grassland",
    "BlackAndWhiteWarbler"="Forest",
    "Bobolink"="Grassland",
    "BrewersSparrow"="Scrub",
    "ChestnutcollaredLongspur"="Grassland",
    "DuskyFlycatcher"="Open Woodland",
    "GrasshopperSparrow"="Grassland",
    "LarkBunting"="Grassland",
    "LarkSparrow"="Grassland",
    "LongbilledCurlew"="Grassland",
    "MacgillivrayWarbler"="Open Woodland",
    "MarbledGodwit"="Marsh",
    "McCownsLongspur"="Grassland",
    "RockPigeon"="Town",
    "SaysPhoebe"="Grassland",
    "SharptailedGrouse"="Grassland",
    "SpraguesPipit"="Grassland",
    "WesternKingbird"="Grassland",
    "WesternMeadowlark"="Grassland",
    "Willet"="Shore-line",
    "WilsonsPhalarope"="Lake/Pond")
TAX[names(tmp), "Hab"] <- tmp

#species <- Species[1]
p <- list()
for (species in Species) {
    cat(species, "\n")
    flush.console()
    y <- load_species_data(species)
    df$D <- 0
    df[rownames(y$SA.Curr), "D"] <- rowSums(y$SA.Curr)
    x <- aggregate(list(p=df$D), list(NSR=df$reg_nsr), mean)
    p[[species]] <- structure(x$p, names=levels(x$NSR))
}

save(TAX, p, file="~/GoogleWork/abmi/bu-helper/bu-helper.RData")

all <- list()
for (species in Species) {
    all[[species]] <- list(
        common_name=as.character(TAX[species, "CommonName"]),
        scientific_name=as.character(TAX[species, "ScientificName"]),
        tsn_id=as.character(TAX[species, "TSNID"]),
        habitat=as.character(TAX[species, "Hab"]),
        nsr=as.list(round(p[[species]], 4)))
}

out <- data.frame(t(sapply(all, unlist)))
outj <- toJSON(all, pretty=TRUE)

write.csv(out, file="~/GoogleWork/abmi/bu-helper/bu-helper.csv")
writeLines(outj, con="~/GoogleWork/abmi/bu-helper/bu-helper.json")
