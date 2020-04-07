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

## IDs in the North
IDN <- get_all_id("north")

## save results (spp x id) + metadata

ST <- get_species_table(mregion="north")
SPP <- rownames(ST)

## matrices to store quartiles
CR1 <- matrix(0, length(SPP), length(unique(vv$RT)))
rownames(CR1) <- SPP
colnames(CR1) <- sort(unique(vv$RT))
CR2 <- CR3 <- CR4 <- CR1

#spp <- SPP[1]
for (spp in SPP) {
    cat(spp, which(SPP==spp), "/", length(SPP), "\n")
    flush.console()
    y <- load_species_data(spp, boot=FALSE)
    Nkm <- rowSums(y$SA.Curr)
    q <- quantile(Nkm, seq(0, 1, 0.25))
    Qkm <- cut(Nkm, q, include.lowest=TRUE)
    levels(Qkm) <- paste0("Q", 1:4)
    names(Qkm) <- names(Nkm)
    Qkm <- Qkm[vv$KM]
    tab <- Xtab(~ RT + Qkm, vv)
    tab <- 100 * tab/sum(tab)
    CR1[spp, rownames(tab)] <- as.numeric(tab[,"Q1"])
    CR2[spp, rownames(tab)] <- as.numeric(tab[,"Q2"])
    CR3[spp, rownames(tab)] <- as.numeric(tab[,"Q3"])
    CR4[spp, rownames(tab)] <- as.numeric(tab[,"Q4"])
}

range(CR1)
range(CR2)
range(CR3)
range(CR4)

out <- data.frame(ST[,c(1:5, 11)],
    Q1=round(CR1, 2),Q2=round(CR2, 2),Q3=round(CR3, 2),Q4=round(CR4, 2))

write.csv(out, row.names=FALSE,
    file="GRE-9-10-TWP-102-103_ABMI-species_2020-04-07.csv")

