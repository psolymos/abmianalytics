#remotes::install_github("mabecker89/abmi.themes")

library(ggplot2)
library(mefa4)
library(sp)
library(raster)
library(cure4insect)
library(rgdal)
library(ggplot2)
library(abmi.themes)
library(intrval)
#set_options(path = "d:/abmi/reports")
load_common_data()
source("~/repos/abmianalytics/birds/00-functions.R")

ROOT <- "~/GoogleWork/abmi/validation"

load(file.path(ROOT, "bg-data-package.RData"))
load(file.path(ROOT, "bg-predictions-ALL.RData"))
load(file.path(ROOT, "bg-predictions-ALLWB.RData"))

rt <- .read_raster_template()

b <- readOGR(dsn=system.file("extdata/AB_bound.geojson", package="cure4insect"))
b <- spTransform(b, proj4string(xy))

load(file.path(ROOT, "veghf-summaries.RData"))
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]

vt <- ddp17$veg_current
vt <- row_std(groupSums(vt, 2, tv[colnames(vt), "Sector61"]))
vt <- as.matrix(vt)
vt1 <- groupMeans(vt, 1, vv$site)
vt1 <- vt1[order(vt1[,1], decreasing=TRUE),]

bc <- c("#576a26", "#919e39", "#b1ab3e", "#d6c350", "#eece5a",
        "#d1a123", "#966521", "#5c361f", "#261d19")


## map with BGs
rownames(xy@coords) <- vv$id_final
tmp <- xy[!duplicated(vv$site),]

png(file.path(ROOT, "fig1-map.png"), width=1000, height=1500)
op <- par(mar=c(0,0,0,0)+0.2)
plot(b, col=bc[5], border=bc[7])
plot(xy, col="white", add=TRUE, pch=".")
set.seed(2)
text(jitter(coordinates(tmp), 10^2),
    label=vv$site[!duplicated(vv$site)], cex=1.5)
par(op)
dev.off()



## BGs with points
png(file.path(ROOT, "fig2-bg-pts.png"), width=1600, height=1600, pointsize=24)
col <- colorRampPalette(bc)(100)
op <- par(mfrow=c(4,4), mar=c(3,3,3,3))
for (i in as.integer(rownames(vt1)))
    plot(xy[vv$site==i,], col=col[ceiling(100*(1-vt[vv$site==i,1]))],
        cex=2, pch=19, main=paste("Site", i))
plot(0, type="n", xlim=c(0, 100), ylim=c(0, 100), axes=FALSE, ann=FALSE,
    xaxs = "i", yaxs = "i")
for (i in 1:100)
    polygon(c(i-1, i, i, i-1), c(0, 0, 20, 20), col=col[i], border=col[i])
axis(1, tick=FALSE, line=-1)
text(50, 30, "Total human footprint (%)", cex=1.5)
par(op)
dev.off()

## THF in BGs table
df1 <- data.frame(Site=rownames(vt1),
    round(vt1[,-1]*100,2),
    TotalHF=round(100*(1-vt1[,1]), 2))
knitr::kable(df1, digits=2, row.names = FALSE)
write.csv(df1, row.names=FALSE, file=file.path(ROOT, "tab1-bg-hf.csv"))





## styled ggplot
if (FALSE) {
p <- ggplot(mammals, aes(x = common_name, y = images, fill = common_name)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Total number of images captured",
       subtitle = "For eight of the most common mammal species",
       caption = "Based on sampling done between 2013 and 2018.",
       y = "Images",
       x = "")

# Change to ABMI theme, remove legend, and use the "main" palette
p1 <- p +
    theme_abmi(font = "montserrat") +
    theme(legend.position = "none") +
    scale_fill_manual(values=bc)

# Finally add ABMI logo; defaults to acronym version
add_logo(p1)


library(colorspace)

cc <- bc
cc <- rev(bc[c(5:1)])
cc <- rev(bc[c(5:9)])
demoplot(cc, "bar")
demoplot(cc, "heatmap")
hclplot(cc)
specplot(cc, type = "o")


}



theme_new <- theme_abmi(font = "montserrat")
theme_new$panel.grid.major.x <- ggplot2::element_line(color = "grey80")



SPP <- names(ALL)

tab <- NULL
for (spp in names(ALL)) {
    for (stage in c("ARU", "Space", "HF")) {
        for (i in names(ALL[[spp]][[stage]])) {
            tmp <- data.matrix(ALL[[spp]][[stage]][[i]]$stats)
            tab <- rbind(tab, data.frame(Species=spp, Stage=stage, Scale=i, t(tmp)))
        }
    }
}
summary(tab)

po <- colMeans(YY[,SPP]>0)
tab$pocc <- po[match(tab$Species, names(po))]

tmp1 <- tab[,c("Species", "Stage", "Scale", "COR", "pocc")]
colnames(tmp1)[4:5] <- c("Value", "Pocc")
tmp1$Statistic <- factor("r", c("r", "R2"))

tmp2 <- tab[,c("Species", "Stage", "Scale", "R2", "pocc")]
colnames(tmp2)[4:5] <- c("Value", "Pocc")
tmp2$Statistic <- factor("R2", c("r", "R2"))

df2 <- data.frame(rbind(tmp1, tmp2))

p1 <- ggplot(df2[df2$Scale != "10",],
             aes(x = Scale, y = Value, fill=Scale)) +
    geom_boxplot() +
    facet_wrap(vars(Statistic, Stage)) +
    labs(title = "Validation statistics",
        subtitle = "Between predicted abundance and observed counts",
        caption = "Statistics: r = rank correlation, R2 = coefficient of determination.",
        y = "Statistic value",
        x = "Scale (x600 m)") +
    theme_abmi(font = "montserrat") +
    theme(legend.position = "none") +
    scale_fill_manual(values=bc[5:1])
#add_logo(p1)
ggsave(file.path(ROOT, "fig3-cor-by.png"))

p2 <- ggplot(tab[tab$Scale!="10" & b01$stage=="HF",],
             aes(x = COR, y = R2, color=Scale, group=Scale)) +
    geom_point() +
    #geom_smooth(se=FALSE) +
    labs(title = "Validation statistics",
        subtitle = "Between predicted abundance and observed counts",
        caption = "Model stages: local + space + landscape.",
        x = "Rank correlation",
        y = "Coefficient of determination") +
    theme_new +
    #theme(legend.position = "none") +
    scale_color_manual(values=bc[5:1])
#add_logo(p2)
ggsave(file.path(ROOT, "fig4-cor.png"))

## all effects move in the positive directions
## interactions are not important
mr <- lm(Value ~ (Stage + Scale + Pocc)^3, tmp1)
mr <- step(mr)
summary(mr)
knitr::kable(coefficients(summary(mr)))

mr2 <- lm(Value ~ (Stage + Scale + Pocc)^3, tmp2)
mr2 <- step(mr2)
summary(mr2)
knitr::kable(coefficients(summary(mr2)))

df2 <- data.frame(r=coefficients(summary(mr)), R2=coefficients(summary(mr)))
df2 <- data.frame(Terms=rownames(df2), df2)
write.csv(df2, row.names=FALSE, file=file.path(ROOT, "tab2-lm.csv"))



## regression is based on mean (per PC value)
## otherwise intercepts are
c_oef <- function(spp, stage="HF") {
    c_oef1 <- function(k, x) {
        xx <- x[[as.character(k)]]$results
        xx$yobs <- xx$yobs / xx$n
        xx$yhat <- xx$yhat / xx$n
        object <- lm(yobs ~ yhat, xx)
        cf <- unname(coef(object))
        ses <- unname(sqrt(diag(vcov(object))))
        level <- 0.95
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        fac <- qt(a, object$df.residual)
        ci <- cf + ses %o% fac
        out <- c(b0=cf[1], b1=cf[2], is00=0 %)(% ci[1,], is11=1 %)(% ci[2,])
        if (0 %)(% ci[1,] & cf[1] < 0)
          out[3] <- -1 * out[3]
        if (0 %)(% ci[2,] & cf[1] < 0)
          out[4] <- -1 * out[4]
        out
    }
    x <- ALL[[spp]][[stage]]
    out <- data.frame(t(sapply(names(x), c_oef1, x=x)))
    out$spp <- factor(spp, SPP)
    out$stage <- factor(stage, c("ARU", "Space", "HF"))
    out$scale <- factor(rownames(out), rownames(out))
    out
}

v <- expand.grid(spp=SPP, stage=c("ARU", "Space", "HF"))
b01 <- do.call(rbind, mapply(c_oef, spp=v$spp, stage=v$stage,
    SIMPLIFY=FALSE))
summary(b01[b01$scale != "10",])
with(b01, table(int=is00, slp=is11))
with(b01, table(int=is00, slp=is11, scale))

with(b01[b01$scale != "10",], table(int=is00, slp=is11, scale))

df3 <- data.frame(with(b01, table(slp=is11, int=is00)))
rownames(df3) <- paste(df3$int, df3$slp)
df3 <- df3[,2:1]
df3$scale1 <- 0
df3$scale2 <- 0
df3$scale3 <- 0
df3$scale4 <- 0
df3$scale5 <- 0

for ( i in c("1", "2", "3", "4", "5")) {
  tmp <- data.frame(with(b01[b01$scale == i & b01$stage=="HF",],
                         table(slp=is11, int=is00)))
  rownames(tmp) <- paste(tmp$int, tmp$slp)
  df3[rownames(tmp), paste0("scale", i)] <- tmp$Freq
}
write.csv(df3, row.names=FALSE, file=file.path(ROOT, "tab3-3x3.csv"))


str(tab)
rownames(tab) <- with(tab, paste(Species, Stage, Scale))
str(b01)
rownames(b01) <- with(b01, paste(spp, stage, scale))

tax <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(tax) <- tax$AOU
tax <- tax[SPP,]


tab2 <- data.frame(tax[match(tab$Species, tax$AOU),c(2,3,6)],
                  tab,
                  b01[rownames(tab), 1:4])
tab2 <- tab2[tab2$Scale != "10",]
tab2$CORU <- tab2$PRC <- tab2$ACC <- NULL
write.csv(tab2, row.names=FALSE, file=file.path(ROOT, "Appendix.csv"))

boxplot(b0 ~ scale, b01,ylim=c(-2,3))
abline(h=c(0,1),col=2:3)
boxplot(b1 ~ scale, b01,ylim=c(-2,3))
abline(h=c(0,1),col=2:3)

plot(b1 ~ b0, b01, pch=19, col="#00000022")
abline(h=0, v=0, col=2)
abline(h=1,col=2, lty=2)

p3 <- ggplot(b01[b01$scale != "10" & b01$stage=="HF",],
            aes(x = b1, y = b0, color=scale, group=scale)) +
    geom_point() +
    geom_smooth(se=FALSE) +
    labs(title = "Validation statistics",
        subtitle = "Between predicted abundance and observed counts",
        caption = "Model stages: local + space + landscape.",
        x = "Slope",
        y = "Intercept") +
    theme_new +
    #theme(legend.position = "none") +
    scale_color_manual(values=bc[5:1]) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=1)
#add_logo(p3)
ggsave(file.path(ROOT, "fig5-b01.png"))




p_lot <- function(x, k, ...) {
    xx <- x[[as.character(k)]]$results
    Max <- max(xx$yobs, xx$yhat)
    plot(yhat ~ jitter(yobs), xx, ylim=c(0, Max), xlim=c(0, Max),
        pch=19, col="#65d7fc", xlab="Observed", ylab="Predicted", ...)
    abline(0,1,lty=2)
    abline(lm(yhat ~ yobs, xx))
    s <- x[[as.character(k)]]$stats
    g <- switch(as.character(k),
        "1"="1x1",
        "2"="2x2",
        "3"="3x3",
        "4"="4x4",
        "5"="5x5",
        "10"="10x10")
    txt <- paste0(c(g, paste0(names(s), " = ", round(s, 3))), collapse="\n")
    mtext(txt, at=0, adj=0, padj=1, line=-1, cex=0.75)
    invisible()
}


