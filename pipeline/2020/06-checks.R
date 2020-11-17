# chacking mammal predictions: exact and approximate

library(mefa4)
library(ggplot2)
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
ROOT <- "s:/AB_data_v2020/Results/pred"

taxon <- "mammals"
SPPn <- dimnames(COEFS2[[taxon]]$north)[[1]]
SPPs <- dimnames(COEFS2[[taxon]]$south)[[1]]
SPP <- sort(unique(c(SPPn, SPPs)))

spp <- "Coyote"

res <- NULL
for (spp in SPP) {
    cat(spp, "\n")

load(file.path(ROOT, taxon, paste0(spp, ".RData"))) # Ncr, Nrf
z <- data.frame(
    LinkID=rownames(Ncr),
    Ref=rowSums(Nrf),
    Curr=rowSums(Ncr))

x <- read.csv(paste0(
    "s:/AB_data_v2020/Results/",
    "Camera mammal models revised June 2020/",
    "Km2 summaries/", spp, ".csv"))
rownames(x) <- x$LinkID

s <- intersect(rownames(x), rownames(z))

u <- data.frame(
    rf_exact=x[s,"Ref"],
    rf_approx=z[s,"Ref"],
    cr_exact=x[s,"Curr"],
    cr_approx=z[s,"Curr"])

res <- rbind(res, c(rf=cor(u[,1:2])[2,1], cr=cor(u[,3:4])[2,1]))
}
rownames(res) <- SPP

k <- lm(rf_approx ~ rf_exact, u)
p <- ggplot(u, aes(x=rf_exact, y=rf_approx)) +
    geom_bin2d() +
    geom_abline(slope=1, col=2) +
    geom_abline(intercept=coef(k)[1], slope=coef(k)[2], col=2, lty=2) +
    labs(title=paste(spp, "Ref")) +
    theme_minimal()

k <- lm(cr_approx ~ cr_exact, u)
p <- ggplot(u, aes(x=cr_exact, y=cr_approx)) +
    geom_bin2d() +
    geom_abline(slope=1, col=2) +
    geom_abline(intercept=coef(k)[1], slope=coef(k)[2], col=2, lty=2) +
    labs(title=paste(spp, "Curr")) +
    theme_minimal()

head(x)
head(z)
compare_sets(x$LinkID, z$LinkID)


## check predictive maps

load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
#load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_2018_wide.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]


xvr <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "rf"])
xvc <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "cr"])
setdiff(colnames(xvr), colnames(xvc))
setdiff(colnames(xvc), colnames(xvr))
i <- intersect(colnames(xvc), colnames(xvr))
d <- (xvr[,i] - xvc[,i]) / 10^6
range(d)
summary(rowSums(d))
table(rowSums(d < 0))
j <- which(rowSums(d < 0) > 13)
u <- d[j,]
u <- u[,colSums(u<0) > 1]

xvc <- dd_2018$veg_current
xvr <- dd_2018$veg_reference
d <- (xvr[,i] - xvc[,i]) / 10^6
range(d)

load("s:/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v61hf2016v3WildFireUpTo2016.Rdata")
