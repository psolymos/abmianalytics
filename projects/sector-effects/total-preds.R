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

TT <- 50

mods <- lapply(sectors, function(s) lm(hf[s,] ~ I((yrs-1999)/TT)))
names(mods) <- sectors
sapply(sectors, function(s) unname(coef(mods[[s]])[2]))
Apr1 <- sapply(sectors, function(s) unname(hf[s,"2015"]+coef(mods[[s]])[2]))

library(forecast)
mods <- lapply(sectors, function(s) ets(ts(hf[s,],start=1999)))
names(mods) <- sectors
preds <- lapply(mods, predict, TT)
Apr2 <- sapply(preds, function(p) p$mean[TT])

U <- z[,paste0("Unit_", sectors)]
#Tot0 <- t(t(U)*A/100)
#head(Tot0)
#head(z[,paste0("Total_", sectors)])
Tot0 <- z[,paste0("Total_", sectors)]
Tot1 <- t(t(U)*Apr1/100)
Tot2 <- t(t(U)*Apr2/100)
colnames(Tot1) <- colnames(Tot0)
colnames(Tot2) <- colnames(Tot0)

Tab <- round(data.frame(MeanTotal=colMeans(Tot0, na.rm=TRUE),
    yr10lm=colMeans(Tot1, na.rm=TRUE),
    yr10ets=colMeans(Tot2, na.rm=TRUE),
    A=A, Alm10=Apr1, Aets10=Apr2), 3)


Tot01 <- rbind(Tot0, Tot1)
Tot01$status <- rep(0:1, each=nrow(Tot0))
mb <- glm(status ~ ., Tot01, family=binomial)
summary(mb)
round(coef(mb),4)

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


z2 <- read.csv("~/GoogleWork/abmi/Custom_Report_2018-01-19.csv")
z2$X <- NULL
rownames(z2) <- z2$SpeciesID
class(z2) <- class(z)
all(rownames(z)==rownames(z2))

par(mfrow=c(3,2))
plot_sector(z, "underhf")
plot_sector(z2, "underhf")
plot_sector(z, "regional")
plot_sector(z2, "regional")
plot_sector(z, "unit", ylim=c(-200,200))
plot_sector(z2, "unit", ylim=c(-200,200))


library(plotrix)

ss <- !is.na(z$Total_Energy) & !is.na(z2$Total_Energy) & z$Taxon=="birds"
WHAT <- "Total_Energy"
ladderplot(cbind(AB=z[ss,WHAT],OSA=z2[ss,WHAT]))
a <- z2[ss,WHAT]-z[ss,WHAT]
names(a) <- rownames(z)[ss]
head(a[order(abs(a), decreasing=TRUE)], 20)

plot(density(z$Total_Energy[ss]), xlim=c(-100,100))
lines(density(z2$Total_Energy[ss]), col=2)
summary(z2$Total_Energy[ss]-z$Total_Energy[ss])




opar <- set_options(path = "w:/reports")
getOption("cure4insect")
load_common_data()
species <- c("Carex.trisperma", "Hymenoxys.richardsonii", "Lithospermum.ruderale",
"Zigadenus.elegans")
subset_common_data(id=get_all_id(),
    species=species)
## see how these compare
system.time(res <- report_all())

z <- do.call(rbind, lapply(res, flatten))
class(z) <- c("c4idf", class(z))

y <- load_species_data(species[1])
x <- calculate_results(y)
x
flatten(x)

