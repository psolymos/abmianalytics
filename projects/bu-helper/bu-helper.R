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


## -- time of year processing from ARU data

library(mefa4)
library(cure4insect)
set_options(verbose=0, path = "d:/abmi/reports")
load_common_data()
spt <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(spt) <- spt$AOU
library(lhreg)
data("lhreg_data")
rt <- .read_raster_template()
KT <- cure4insect:::.c4if$KT

rasterize_results_cr <- function (y, current=TRUE)
{
    NC <- if (current)
        rowSums(y$SA.Curr) else rowSums(y$SA.Ref)
    KT$NC <- NC[match(rownames(KT), names(NC))]
    KT$NC[is.na(KT$NC)] <- 0
    r <- .make_raster(KT$NC, rc = KT, rt = rt)
    r
}


SPP <- intersect(lhreg_data$spp, spt$AOU[spt$map.pred])
WR <- as.character(lhreg_data$spp)[lhreg_data$Mig == "WR"]
SD <- as.character(lhreg_data$spp)[lhreg_data$Mig == "SD"]
zz <- matrix(0, 12, length(SPP))
dimnames(zz) <- list(1:12, SPP)
zz[5:9,] <- 1
zz[,colnames(zz) %in% WR] <- 1
zz[4,colnames(zz) %in% SD] <- 1
zz <- t(zz)
zz

Species <- as.character(spt[SPP, "sppid"])
plist <- list()
for (spp in Species) {
    cat(spp, "\n")
    y <- load_species_data(spp)
    r <- rasterize_results_cr(y)
    p <- p_bird(r, "ha", 7)
    plist[[spp]] <- p
}
t(sapply(plist, function(z) range(values(z),na.rm=T)))

for (spp in Species) {
    p <- plist[[spp]]
    MAX <- max(values(p), na.rm=TRUE)
    p <- p / MAX
    plist[[spp]] <- p
}

names(plist) <- SPP
rr <- stack(plist)

writeRaster(rr, "d:/abmi/AB_data_v2018/data/analysis/birds/BU-helper.tif", overwrite=TRUE)
write.csv(zz,file="d:/abmi/AB_data_v2018/data/analysis/birds/BU-helper.csv")
