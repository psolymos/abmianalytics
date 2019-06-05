#' # Comparing HFI mine poly feature types based on NDVI
#'
#' Sentinel - NDVI (Normalized Difference Vegetation Index)
#'
#' This most known and used vegetation index is a simple, but effective VI for
#' quantifying green vegetation. It normalizes green leaf scattering in the
#' Near Infra-red wavelength and chlorophyll absorption in the red wavelength.
#'
#' Values description: The value range of an NDVI is -1 to 1.
#' Negative values of NDVI (values approaching -1) correspond to water.
#' Values close to zero (-0.1 to 0.1) generally correspond to barren areas of rock, sand, or snow.
#' Low, positive values represent shrub and grassland (approximately 0.2 to 0.4),
#' while high values indicate temperate and tropical rainforests (values approaching 1).
#'
x <- read.csv("d:/abmi/AB_data_v2019/data/raw/HFI_Mines_NDVI.csv")
summary(x)
hist(x$NDVImean)
table(x$Mines_FEATURE_TY, is.na(x$NDVImean), useNA="a")
x <- x[!is.na(x$NDVImean),]

m <- lm(NDVImean ~ Mines_FEATURE_TY-1, x)
summary(m)

boxplot(NDVImean ~ Mines_FEATURE_TY-1, x)
abline(h=c(-0.1,0, 0.1), col=2)

a <- aggregate(x$NDVImean, list(FTY=x$Mines_FEATURE_TY), quantile, c(0.5, 0, 1))
a <- a[order(a$x[,"50%"]),]
a

x$FEATURE_TY <- as.character(x$Mines_FEATURE_TY)
x$FEATURE_TY <- factor(x$FEATURE_TY, as.character(a$FTY))
#levels(x$FEATURE_TY) <- paste(LETTERS[1:nlevels(x$FEATURE_TY)], x$FEATURE_TY)

library(ggplot2)

p <- ggplot(x, aes(y=NDVImean, x=FEATURE_TY, fill=FEATURE_TY, color=FEATURE_TY)) +
  geom_violin() + coord_flip() + theme(legend.position = "none") +
  geom_hline(yintercept = 0) + geom_hline(yintercept = -0.1, lty=2) + geom_hline(yintercept = 0.1, lty=2)
p
