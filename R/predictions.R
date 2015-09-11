library(mefa4)

ROOT <- "c:/p/AB_data_v2015"
OUTDIR1 <- "e:/peter/sppweb2015/birds-pred-1/"
OUTDIRB <- "e:/peter/sppweb2015/birds-pred-B/"

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
#source("~/repos/bragging/R/glm_skeleton.R")
#source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
regs <- levels(kgrid$LUFxNSR)
useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
useS <- kgrid$NRNAME == "Grassland"

e <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-full-withrevisit.Rdata"), envir=e)
tax <- e$TAX
rm(e)
tax$file <- nameAlnum(as.character(tax$English_Name), "mixed", "")

## model for species
fl <- list.files(file.path(ROOT, "out", "birds", "results"))
fln <- fl[grep("-north_", fl)]
fln <- sub("birds_abmi-north_", "", fln)
fln <- sub(".Rdata", "", fln)
fls <- fl[grep("-south_", fl)]
fls <- sub("birds_abmi-south_", "", fls)
fls <- sub(".Rdata", "", fls)
SPP <- union(fln, fls)

load(file.path(ROOT, "out", "transitions", paste0(regs[1], ".Rdata")))
Aveg <- rbind(colSums(trVeg))
rownames(Aveg) <- regs[1]
colnames(Aveg) <- colnames(trVeg)
Asoil <- rbind(colSums(trSoil))
rownames(Asoil) <- regs[1]
colnames(Asoil) <- colnames(trSoil)

for (i in 2:length(regs)) {
    cat(regs[i], "\n");flush.console()
    load(file.path(ROOT, "out", "transitions", paste0(regs[i], ".Rdata")))
    Aveg <- rbind(Aveg, colSums(trVeg))
    rownames(Aveg) <- regs[1:i]
    Asoil <- rbind(Asoil, colSums(trSoil))
    rownames(Asoil) <- regs[1:i]
}

spp <- "CAWA"

for (spp in SPP) {

cat(spp, "\n");flush.console()

load(file.path(OUTDIR1, spp, paste0(regs[1], ".Rdata")))
rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
rownames(pxScr1) <- rownames(pxSrf1) <- names(Cells)
pxNcr <- pxNcr1
pxNrf <- pxNrf1
pxScr <- pxScr1
pxSrf <- pxSrf1
pSoil <- pSoil1
for (i in 2:length(regs)) {
    cat(spp, regs[i], "\n");flush.console()
    load(file.path(OUTDIR1, spp, paste0(regs[i], ".Rdata")))
    rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
    rownames(pxScr1) <- rownames(pxSrf1) <- names(Cells)
    pxNcr <- rbind(pxNcr, pxNcr1)
    pxNrf <- rbind(pxNrf, pxNrf1)
    pxScr <- rbind(pxScr, pxScr1)
    pxSrf <- rbind(pxSrf, pxSrf1)
    pSoil <- c(pSoil, pSoil1)
}

pxNcr <- pxNcr[rownames(kgrid),]
pxNrf <- pxNrf[rownames(kgrid),]
pxScr <- pxScr[rownames(kgrid),]
pxSrf <- pxSrf[rownames(kgrid),]
pSoil <- pSoil[rownames(kgrid)]

wS <- 1-kgrid$pAspen
if (!NSest["north"])
    wS[] <- 1
if (!NSest["south"])
    wS[] <- 0

pxcr <- wS * pxScr + (1-wS) * pxNcr
pxcr[useN] <- pxNcr[useN]
pxcr[useS] <- pxScr[useS]
pxrf <- wS * pxSrf + (1-wS) * pxNrf
pxrf[useN] <- pxNrf[useN]
pxrf[useS] <- pxSrf[useS]

if (!NSest["north"]) {
    pxcr[useN] <- -1#NA
    pxrf[useN] <- -1#NA
}
if (!NSest["south"]) {
    pxcr[useS] <- -1#NA
    pxrf[useS] <- -1#NA
}

km <- data.frame(LinkID=kgrid$Row_Col,
    Ref=pxrf,
    Curr=pxcr)
NAM <- as.character(tax[spp, "English_Name"])

write.csv(km, row.names=FALSE,
    paste0("e:/peter/sppweb2015/birds-pred/", paste0(as.character(tax[spp, "file"]), ".csv")))
}


dim(pxNcr)
dim(kgrid)

with(kgrid, plot(X, Y, col=ifelse(pSoil == 1 & kgrid$pWater < 0.99, 2, 1), pch="."))

* save csv
* ??? fix PUMA
* plot from csv
* sector effects
* sector effects table
* regional abund
* habitat based table -- with upland/lowland
* create guild classification
* check how guild figures are made (reuse Dave code)
* produce old forest guild figures
* CoV map from Boot (lastly)

