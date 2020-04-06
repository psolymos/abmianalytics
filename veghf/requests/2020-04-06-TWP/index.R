## Need c4i based suitabilities for spp in:
## Twp 102,103, Rg 9,10

library(mefa4)
library(cure4insect)
set_options(path = "d:/abmi/reports")
load_common_data()

## QS IDs are built up as MER-RGE-TWP-SEC-QS
## all should be 36*4*4=576
vv <- expand.grid(MER=6, RGE=9:10, TWP=102:103, SEC=1:36, QS=c("NE", "NW", "SE", "SW"))
vv$QSID <- paste(vv$MER, vv$RGE, vv$TWP, vv$SEC, vv$QS, sep="-")
rownames(vv) <- vv$QSID
vv$RT <- paste(vv$MER, vv$RGE, vv$TWP, sep="-")
vv$KM <- as.character(cure4insect:::.c4if$QT2KT[rownames(vv)])
## check
#ID <- sort(qs2km(rownames(vv)))
#z <- sort(unique(vv$KM))
#all(ID==z)

## save results (spp x id) + metadata

ST <- get_species_table()
SPP <- rownames(ST)

CR <- matrix(0, length(SPP), length(unique(vv$RT)))
rownames(CR) <- SPP
colnames(CR) <- sort(unique(vv$RT))

#spp <- SPP[1]
for (spp in SPP) {
    cat(spp, which(SPP==spp), "/", length(SPP), "\n")
    flush.console()
    y <- load_species_data(spp, boot=FALSE)
    names(y)
    cr <- rowMeans(groupSums(y$SA.Curr[vv$KM,], 1, vv$RT))
    CR[spp, names(cr)] <- cr
}

out <- data.frame(ST, CR)
write.csv(out, row.names=FALSE,
    file="GRE-9-10-TWP-102-103_ABMI-species_2020-04-06.csv")

