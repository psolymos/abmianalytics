#devtools::install_github("ABbiodiversity/cure4insect")
library(cure4insect)
library(knitr)
library(rgdal)
#set_options(verbose=0)
load_common_data()

load(system.file("extdata/raw_all.rda", package="cure4insect"))
z <- do.call(rbind, lapply(res, flatten))
class(z) <- c("c4idf", class(z))

sectors <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation")

A <- colMeans(z[,paste0("Area_", sectors)], na.rm=TRUE)

load("~/repos/abmianalytics/apps/hfchange/hfchange.rda")
str(HF)
dimnames(HF)[[2]] <- c("Agriculture","Forestry","Energy","RuralUrban","Transportation","Misc")
hf <- matrix(NA, 5, 2015-1999+1)
rownames(hf) <- sectors
yrs <- 1999:2015
colnames(hf) <- yrs
for (i in dimnames(HF)[[3]])
    hf[,i] <- colMeans(HF[,sectors,i])

mods <- lapply(sectors, function(s) lm(hf[s,] ~ I((yrs-1999)/10)))
names(mods) <- sectors
sapply(sectors, function(s) unname(coef(mods[[s]])[2]))
Apr1 <- sapply(sectors, function(s) unname(hf[s,"2015"]+coef(mods[[s]])[2]))

library(forecast)
mods <- lapply(sectors, function(s) ets(ts(hf[s,],start=1999)))
names(mods) <- sectors
preds <- lapply(mods, predict, 10)
Apr2 <- sapply(preds, function(p) p$mean[10])

U <- z[,paste0("Unit_", sectors)]
#Tot0 <- t(t(U)*A/100)
#head(Tot0)
#head(z[,paste0("Total_", sectors)])
Tot0 <- z[,paste0("Total_", sectors)]
Tot1 <- t(t(U)*Apr1/100)
Tot2 <- t(t(U)*Apr2/100)

Tab <- round(data.frame(MeanTotal=colMeans(Tot0, na.rm=TRUE),
    yr10lm=colMeans(Tot1, na.rm=TRUE),
    yr10ets=colMeans(Tot2, na.rm=TRUE),
    A=A, Alm10=Apr1, Aets10=Apr2), 3)


h <- function(tmp) {
    v <- tmp[species, sectors]
    pmax(-100, pmin(100, v))
}
op <- par(mfrow=c(3,1))
tmp <- plot_sector(z, type="regional", main="Regional")
points(1:5, h(tmp), pch=4, col=1, cex=3)

tmp <- plot_sector(z, type="underhf", main="Under HF")
points(1:5, h(tmp), pch=4, col=1, cex=3)

tmp <- plot_sector(z, type="unit", main="Unit", ylim=c(-200,200))
points(1:5, h(tmp), pch=4, col=1, cex=3)
par(op)


