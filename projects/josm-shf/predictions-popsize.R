library(mefa4)

ROOT <- "e:/peter/AB_data_v2016"
STAGE <- list(veg = 7) # hab=5, hab+clim=6, hab+clim+shf=7

OUTDIR1 <- paste0("e:/peter/josm/2017/stage", STAGE$veg, "/pred1")
OUTDIRB <- paste0("e:/peter/josm/2017/stage", STAGE$veg, "/predB")


load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata"))
#source("~/repos/bragging/R/glm_skeleton.R")
#source("~/repos/abmianalytics/R/results_functions.R")
#source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/abmianalytics/R/maps_functions.R")
regs <- levels(kgrid$LUFxNSR)
kgrid$useN <- !(kgrid$NRNAME %in% c("Grassland", "Parkland") | kgrid$NSRNAME == "Dry Mixedwood")
kgrid$useN[kgrid$NSRNAME == "Dry Mixedwood" & kgrid$POINT_Y > 56.7] <- TRUE
kgrid$useS <- kgrid$NRNAME == "Grassland"
kgrid$useBCR6 <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"

e <- new.env()
load(file.path(ROOT, "out", "birds", "data", "data-josmshf.Rdata"), envir=e)
xn <- e$DAT
mods <- e$mods
yy <- e$YY
#BB <- e$BB
#tax <- droplevels(TAX[colnames(yyn),])
rm(e)

fln <- list.files(file.path(ROOT, "out", "birds", "results", "josmshf"))
fln <- sub("birds_abmi-josmshf_", "", fln)
fln <- sub(".Rdata", "", fln)

SPP <- fln
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL")

PREDS <- matrix(0, sum(kgrid$useBCR6), length(SPP))
rownames(PREDS) <- rownames(kgrid)[kgrid$useBCR6]
colnames(PREDS) <- SPP
PREDS0 <- PREDS

AREA_ha <- (1-kgrid$pWater) * kgrid$Area_km2 * 100
AREA_ha <- AREA_ha[kgrid$useBCR6]

for (spp in SPP) {
    cat(spp, "--------------------------------------\n");flush.console()
    fl <- list.files(file.path(OUTDIR1, spp))
    ssRegs <- gsub("\\.Rdata", "", fl)
    pxNcr <- pxNrf <- NULL
    for (i in ssRegs) {
        cat(spp, i, "\n");flush.console()
        load(file.path(OUTDIR1, spp, paste0(i, ".Rdata")))
        rownames(pxNcr1) <- rownames(pxNrf1) <- names(Cells)
        pxNcr <- rbind(pxNcr, pxNcr1)
        pxNrf <- rbind(pxNrf, pxNrf1)
    }
    PREDS[,spp] <- pxNcr[rownames(PREDS),]
    PREDS0[,spp] <- pxNrf[rownames(PREDS0),]
}
N <- colSums(PREDS*AREA_ha) / 10^6
save(AREA_ha, N, PREDS, PREDS0, file=file.path(OUTDIR1, "predictions.Rdata"))

## looking at results

e <- new.env()
load("e:/peter/josm/2017/stage6/pred1/predictions.Rdata", envir=e)
N6 <- e$N
e <- new.env()
load("e:/peter/josm/2017/stage7/pred1/predictions.Rdata", envir=e)
N7 <- e$N
rm(e)

## ---

N <- colSums(PREDS*AREA_ha) / 10^6
#N <- N[N < max(N)]
summary(N)

## PIF table
pif <- read.csv("~/Dropbox/bam/PIF-AB/popBCR-6AB_v2_22-May-2013.csv")
mefa4::compare_sets(tax$English_Name, pif$Common_Name)
setdiff(tax$English_Name, pif$Common_Name)
pif <- pif[match(tax$English_Name, pif$Common_Name),]


## roadside_bias
load(file.path(ROOT, "out", "birds", "josmshf", "roadside_bias.Rdata"))

load(file.path(ROOT, "out", "birds", "data", "mean-qpad-estimates.Rdata"))
qpad_vals <- qpad_vals[rownames(tax),]

## roadside avoidance
library(mefa4)
load(file.path(ROOT, "out", "birds", "josmshf", "roadside_avoidance.Rdata"))
tmp <- cbind(ROAD=rai_data$ROAD, rai_pred)
rai <- groupSums(tmp[BBn[,1],], 1, rai_data$HAB[BBn[,1]], TRUE)
rai <- t(t(rai) / colSums(rai))
RAI <- 1 - colSums(rai[,1] * rai)
summary(RAI)
RAIc <- RAI-RAI["ROAD"]

#yy <- cbind(ALL=1, ROAD=xnn[BBn[,1],"ROAD01"],
#    ifelse(as.matrix(yyn[BBn[,1],]) > 0, 1, 0))
#rai <- groupSums(yy, 1, xnn$hab1[BBn[,1]], TRUE)
#n <- rai[,"ALL"]
#rai <- rai[,-1]
#rai <- t(t(rai) / colSums(rai))
#sai <- groupSums(yy, 1, xnn$hab1[BBn[,1]], TRUE)
#RAI <- 1 - colSums(rai[,1] * rai)

pop <- tax[,c("Species_ID", "English_Name", "Scientific_Name", "Spp")]
pop$RAI <- RAI[match(rownames(pop), names(RAI))]
pop$RAIc <- RAIc[match(rownames(pop), names(RAIc))]
pop$RAIroad <- RAI["ROAD"]
pop$Don <- roadside_bias[rownames(pop), "on"]
pop$Doff <- roadside_bias[rownames(pop), "off"]
pop$DeltaRoad <- roadside_bias[rownames(pop), "onoff"]
pop$Nqpad <- colSums(PREDS*AREA_ha) / 10^6 # M males
pop$Nqpad[pop$Nqpad > 1000] <- NA
pop$Npif <- (pif$Pop_Est / pif$Pair_Adjust) / 10^6 # M males
pop$DeltaObs <- pop$Nqpad / pop$Npif
pop$TimeAdj <- pif$Time_Adjust
pop$MDD <- pif$Detection_Distance_m
pop$p3 <- 1-exp(-3 * qpad_vals$phi0)
pop$EDR <- qpad_vals$phi0 * 100
pop$DeltaTime <- (1/pop$p3)/pop$TimeAdj
pop$DeltaDist <- pop$MDD^2 / pop$EDR^2
pop$DeltaExp <- pop$DeltaRoad * pop$DeltaTime * pop$DeltaDist
pop$DeltaRes <- pop$DeltaObs / pop$DeltaExp
pop <- pop[rowSums(is.na(pop))==0,]

#write.csv(pop, row.names=FALSE, file="~/Dropbox/bam/PIF-AB/qpad-pif-results.csv")

boxplot(log(pop[,c("DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")]))
abline(h=0, col=2)

boxplot(log(pop[,c("DeltaObs", "DeltaExp")]))
abline(h=0, col=2)

mat <- log(pop[,c("DeltaObs", "DeltaExp", "DeltaRoad", "DeltaTime", "DeltaDist", "DeltaRes")])
rnd <- runif(nrow(pop), -0.1, 0.1)
boxplot(mat, range=0)
for (i in 2:ncol(mat))
    segments(x0=i+rnd-1, x1=i+rnd, y0=mat[,i-1], y1=mat[,i], col="lightgrey")
for (i in 1:ncol(mat))
    points(i+rnd, mat[,i], col="darkgrey", pch=19)
abline(h=0, col=2, lwd=2)
boxplot(mat, range=0, add=TRUE)

with(pop, plot(RAI, log(DeltaRes), type="n"))
abline(h=0, v=RAI["ROAD"], col=2, lwd=2)
with(pop, text(RAI, log(DeltaRes), rownames(pop), cex=0.75))

boxplot(pop[,c("Npif", "Nqpad")], ylim=c(0,10))


