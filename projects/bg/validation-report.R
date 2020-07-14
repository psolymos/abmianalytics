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
load(file.path(ROOT, "bg-predictions-ALL-2020-06-23.RData"))
load(file.path(ROOT, "bg-predictions-ALLWB-2020-06-23.RData"))
load(file.path(ROOT, "North-xy.RData"))

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
theme_new <- theme_abmi(font = "montserrat")
theme_new$panel.grid.major.x <- ggplot2::element_line(color = "grey80")


## map with BGs
rownames(xy@coords) <- vv$id_final
tmp <- xy[!duplicated(vv$site),]
xy2 <- DATss[,c("X","Y")]
coordinates(xy2) <- ~ X + Y
proj4string(xy2) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy2 <- spTransform(xy2, proj4string(xy))

png(file.path(ROOT, "fig1-map.png"), width=1000, height=1500)
op <- par(mar=c(0,0,0,0)+0.2)
plot(b, col=bc[5], border=bc[7])
plot(xy, col="white", add=TRUE, pch=".")
set.seed(2)
text(jitter(coordinates(tmp), 10^2),
    label=vv$site[!duplicated(vv$site)], cex=1.5)
par(op)
dev.off()

png(file.path(ROOT, "fig0-map-train.png"), width=1000, height=1500)
op <- par(mar=c(0,0,0,0)+0.2)
plot(b, col=bc[5], border=bc[7])
plot(xy2, col=bc[9], add=TRUE, pch=19, cex=0.4)
set.seed(2)
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




## evaluate statistics

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
levels(df2$Scale) <- paste0(levels(df2$Scale), "x", levels(df2$Scale))

p1 <- ggplot(df2[df2$Scale != "10x10",],
             aes(x = Scale, y = Value, fill=Scale)) +
    geom_boxplot() +
    facet_grid(rows=vars(Statistic), cols=vars(Stage)) +
    labs(title = "Validation statistics",
        subtitle = "Between predicted abundance and observed counts",
        caption = "Statistics: r = rank correlation, R2 = coefficient of determination.",
        y = "Statistic value",
        x = "Scale") +
    theme_abmi(font = "montserrat") +
    theme(legend.position = "none") +
    scale_fill_manual(values=bc[5:1])
#add_logo(p1)
p1
ggsave(file.path(ROOT, "fig3-cor-by.png"), p1)

levels(tab$Scale) <- paste0(levels(tab$Scale), "x", levels(tab$Scale))

p2 <- ggplot(tab[tab$Scale!="10x10" & tab$Stage=="HF",],
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
p2
ggsave(file.path(ROOT, "fig4-cor.png"), p2, width=7, height=5)

## all effects move in the positive directions
## interactions are not important
#mr <- lm(Value ~ (Stage + Scale + Pocc)^3, tmp1[tmp1$Scale != "10",])
#mr <- step(mr)
mr <- lm(Value ~ Stage + Scale + Pocc, tmp1[tmp1$Scale != "10",])
summary(mr)
knitr::kable(coefficients(summary(mr)))

#mr2 <- lm(Value ~ (Stage + Scale + Pocc)^3, tmp2[tmp2$Scale != "10",])
#mr2 <- step(mr2)
mr2 <- lm(Value ~ Stage + Scale + Pocc, tmp2[tmp2$Scale != "10",])
summary(mr2)
knitr::kable(coefficients(summary(mr2)))

df2 <- data.frame(r=coefficients(summary(mr)), R2=coefficients(summary(mr)))
df2 <- data.frame(Terms=rownames(df2), df2)
write.csv(df2, row.names=FALSE, file=file.path(ROOT, "tab2-lm.csv"))

## try poisson
## Over: OBS < PRED, Under: OBS > PRED
c_oefp <- function(spp, stage="HF") {
    c_oefp1 <- function(k, x, n=10^3) {
        xx <- x[[as.character(k)]]$results
        level <- 0.95
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        PI <- matrix(0, nrow(xx), 2)
        for (i in 1:nrow(xx))
            PI[i,] <- quantile(rpois(n, xx$yhat[i]), a)
        out <- structure(
          as.numeric(table(factor(xx$yobs %[c]% PI, c(-1L, 0L, 1L)))),
          names=c("Over", "Equal", "Under" ))
        out
    }
    x <- ALL[[spp]][[stage]]
    out <- data.frame(t(sapply(names(x), c_oefp1, x=x)))
    out$spp <- factor(spp, SPP)
    out$stage <- factor(stage, c("ARU", "Space", "HF"))
    out$scale <- factor(rownames(out), rownames(out))
    out
}

n <- 10^3
level <- 0.95
a <- (1 - level)/2
a <- c(a, 1 - a)
for (spp in names(ALL)) {
    for (stage in c("ARU", "Space", "HF")) {
        for (i in names(ALL[[spp]][[stage]])) {
            xx <- ALL[[spp]][[stage]][[i]]$results
            PI <- matrix(0, nrow(xx), 2)
            for (j in 1:nrow(xx))
                PI[j,] <- quantile(rpois(n, xx$yhat[j]), a)
            ALL[[spp]][[stage]][[i]]$results$inPI <- xx$yobs %[c]% PI
        }
    }
}



## regression is based on mean (per PC value)
## otherwise intercepts are
c_oef <- function(spp, stage="HF") {
    c_oef1 <- function(k, x) {
        xx <- x[[as.character(k)]]$results
        xx$yobs <- xx$yobs / xx$n
        xx$yhat <- xx$yhat / xx$n
        object <- lm(yobs ~ yhat, xx)
        sig <- sigma(object)
        cf <- unname(coef(object))
        ses <- unname(sqrt(diag(vcov(object))))
        level <- 0.95
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        fac <- qt(a, object$df.residual)
        ci <- cf + ses %o% fac
        out <- c(b0=cf[1], b1=cf[2], is00=0 %)(% ci[1,], is11=1 %)(% ci[2,], sig=sig)
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
b01 <- do.call(rbind, pbapply::pbmapply(c_oef, spp=v$spp, stage=v$stage,
    SIMPLIFY=FALSE))
b01$Variance <- b01$sig^2
summary(b01[b01$scale != "10",])
with(b01, table(int=is00, slp=is11))
with(b01, table(int=is00, slp=is11, scale))

p01 <- do.call(rbind, pbapply::pbmapply(c_oefp, spp=v$spp, stage=v$stage,
    SIMPLIFY=FALSE))
p01$pEqual <- p01$Equal/(p01$Over+p01$Equal+p01$Under)
p01$pOver <- p01$Over/(p01$Over+p01$Equal+p01$Under)
p01$pUnder <- p01$Under/(p01$Over+p01$Equal+p01$Under)

p01$Scale <- paste0(p01$scale, "x", p01$scale)
with(b01[b01$scale != "10",], table(int=is00, slp=is11, scale))

df3 <- rbind(df2[,colnames(df2) != "Pocc",],
                  data.frame(Species=p01$spp,
                             Stage=p01$stage,
                             Scale=paste0(p01$scale, "x", p01$scale),
                             Value=p01$pEqual,
                             Statistic="Coverage"))

tmp3 <- df3[df3$Scale != "10x10" & df3$Statistic=="Coverage",]
tmp3$Pocc <- po[match(tmp3$Species, names(po))]

mr3 <- lm(Value ~ Stage + Scale + Pocc, tmp3)
#mr3 <- step(mr3)
summary(mr3)

p1x <- ggplot(df3[df3$Scale != "10x10",],
             aes(x = Scale, y = Value, fill=Scale)) +
    geom_boxplot() +
    facet_grid(rows=vars(Statistic), cols=vars(Stage)) +
    labs(title = "Validation statistics",
        subtitle = "Between predicted abundance and observed counts",
        caption = "Statistics: r = rank correlation, R2 = coefficient of determination.",
        y = "Statistic value",
        x = "Scale") +
    theme_abmi(font = "montserrat") +
    theme(legend.position = "none") +
    scale_fill_manual(values=bc[5:1])
p1x
ggsave(file.path(ROOT, "fig3x-cor-by.png"), p1x)


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

## explanatory figure

png(file.path(ROOT, "fig5-3x3.png"), width=1200, height=1200, pointsize=24)
op <- par(mfrow=c(3,3), mar=c(4,4,4,4))
for (i in 1:3) {
  Int <- c(0.25, 0, -0.25)[i]
  for (j in 1:3) {
    Slp <- c(0.5, 1, 2)[j]
    bg <- paste0(bc[7], "88")
    if (i==1 && j==1)
      bg <- "lightgrey"
    if (i==2 && j==2)
      bg <- bc[5]
    if (i==3 && j %in% 2:3)
      bg <- "lightgrey"
    if (i %in% 2:3 && j==1)
      bg <- paste0(bc[3], "88")
    plot(0, type="n", axes=FALSE, ann=FALSE, xlim=c(0,1), ylim=c(0,1),
         xaxs = "i", yaxs = "i", asp=1)
    polygon(c(0,1,1,0), c(0,0,1,1), col=bg, border=NA)
    abline(0,1,lty=2,lwd=7, col="white")
    abline(Int, Slp, lty=1,lwd=9, col=bc[1])

    if (i==1)
      switch(j,
             "1"=mtext(expression(beta[1]<1), 3),
             "2"=mtext(expression(beta[1]==1), 3),
             "3"=mtext(expression(beta[1]>1), 3))
    if (j==3)
      switch(i,
             "1"=mtext(expression(beta[0]>0), 4, line=1),
             "2"=mtext(expression(beta[0]==0), 4, line=1),
             "3"=mtext(expression(beta[0]<0), 4, line=1))
    if (i==2 & j==1)
      mtext("Overestimate\nObserved < Predicted", 1, col=bc[3], line=4)
    if (i==1 & j==3)
      mtext("Underestimate\nObserved > Predicted", 1, col=bc[7], line=4)
  }
}
par(op)
title(xlab="Predicted abundance")
title(ylab="Observed abundance")
dev.off()


## writing appendix
str(tab)
rownames(tab) <- with(tab, paste(Species, Stage, Scale))
str(b01)
rownames(b01) <- with(b01, paste(spp, stage, scale))
str(p01)
rownames(p01) <- with(p01, paste(spp, stage, scale))

tax <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(tax) <- tax$AOU
tax <- tax[SPP,]

tab2 <- data.frame(tax[match(tab$Species, tax$AOU),c(2,3,6)],
                  tab,
                  b01[rownames(tab), 1:4],
                  p01[rownames(tab), 7:9])
tab2 <- tab2[tab2$Scale != "10",]
#tab2$CORU <- tab2$PRC <- tab2$ACC <- NULL
write.csv(tab2, row.names=FALSE, file=file.path(ROOT, "Appendix.csv"))

ggplot(tab2, aes(x=Scale, y=pEqual)) +
  geom_boxplot() +
  facet_wrap(vars(Stage))

ggplot(tab2[tab2$Stage=="HF",],
       aes(x=R2, y=pEqual)) +
  geom_point() +
  facet_wrap(vars(Scale)) +
  geom_smooth(method="lm")

ggplot(tab2[tab2$Stage=="HF",],
       aes(x=R2, y=pUnder, group=Scale, color=Scale)) +
  geom_point() +
  scale_color_manual(values=bc[1:5])

round(as.dist(cor(tab2[,sapply(tab2, is.numeric)])),2)

ps <- ggplot(b01[b01$scale != "10" & b01$stage=="HF",],
            aes(x = scale, y=sig, color=scale, group=scale, fill=scale)) +
  geom_violin() +
  scale_fill_manual(values=bc[5:1]) +
  scale_color_manual(values=bc[5:1]) +
    labs(title = "Residual variance",
        subtitle = "Between predicted abundance and observed counts",
        caption = "Model stages: local + space + landscape.",
        x = "Scale",
        y = "Variance") +
    theme_new
ps
#ggsave(file.path(ROOT, "fig5b-sigma.png"))

tmp <- droplevels(p01[p01$scale != "10" & p01$stage=="HF",])
ddd <- data.frame(scale=rep(tmp$scale, 3),
                  Bias=factor(rep(c("OBS > PRED", "OBS = PRED", "OBS < PRED"),
                    each=nrow(tmp)), c("OBS > PRED", "OBS = PRED", "OBS < PRED")),
                  value=c(tmp$pUnder, tmp$pEqual, tmp$pOver))
#tmp <- droplevels(p01[p01$scale != "10",])
#ddd <- data.frame(scale=tmp$scale,
#                  stage=tmp$stage,
#                  Bias=factor(rep(c("OBS > PRED", "OBS = PRED", "OBS < PRED"),
#                    each=nrow(tmp)), c("OBS > PRED", "OBS = PRED", "OBS < PRED")),
#                  value=c(tmp$pUnder, tmp$pEqual, tmp$pOver))

levels(ddd$scale) <- paste0(levels(ddd$scale), "x", levels(ddd$scale))
ps <- ggplot(ddd,
            aes(x = scale, y=value, fill=Bias)) +
    geom_bar(position="fill", stat="identity") +
    labs(title = "Prediction interval",
        subtitle = "Interval should contain observed counts",
        caption = "Model stages: local + space + landscape.",
        x = "Scale",
        y = "Coverage (%)") +
    theme_new +
    scale_fill_manual(values=bc[c(8,5,2)]) +
    scale_y_continuous(labels = scales::percent_format())
ps
ggsave(file.path(ROOT, "fig5b-PI.png"), width=6, height=5)

levels(b01$scale) <- paste0(levels(b01$scale), "x", levels(b01$scale))
p3 <- ggplot(b01[b01$scale != "10x10" & b01$stage=="HF",],
            aes(x = b1, y = b0, color=scale, group=scale)) +
    #geom_point() +
    geom_density2d() +
    facet_wrap(vars(scale)) +
    labs(title = "Regression summaries",
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
p3
ggsave(file.path(ROOT, "fig5-b01.png"))


## species specific results

#spp <- "ALFL"
stdfun <- function(x, w, s) {
      x <- data.frame(t(x))
      colnames(x) <- c("Value", "Med", "Lwr", "Upr")
      x$Scale <- as.factor(sc2)
      x$Stat <- factor(rep(s, length(sc1)), c("r", "R2"))
      x$Type <- factor(rep(w, length(sc1)), c("Between", "Within"))
      x$Species <- factor(rep(spp, 5), SPP)
      x$Signif <- 0
      x$Signif[x$Value < x$Lwr] <- -1
      x$Signif[x$Value > x$Upr] <- 1
      x
}
sc1 <- c("1", "2", "3", "4", "5")
sc2 <- c("1x1", "2x2", "3x3", "4x4", "5x5")
#sc1 <- c("1", "2", "3", "4", "5", "10")
#sc2 <- c("1x1", "2x2", "3x3", "4x4", "5x5", "10x10")

for (spp in SPP) {
  cat(spp, "\n")
  vals <- NULL
  for (stage in names(ALL[[spp]])) {
    for (scale in sc1) {
      tmp <- data.frame(ALL[[spp]][[stage]][[scale]]$results)
      tmp$Scale <- factor(paste0(scale, "x", scale), sc2)
      tmp$Stage <- factor(stage, names(ALL[[spp]]))
      vals <- rbind(vals, tmp)
    }
  }
  p40 <- ggplot(vals, aes(x = yhat, y = yobs)) +
    geom_point(color=ifelse(vals$inPI==0,paste0(bc[1],"80"), paste0(bc[5],"80"))) +
    geom_smooth(method="lm", color=bc[8], se=FALSE) +
    #facet_grid(rows=vars(Scale), cols=vars(Stage), scales="free") +
    facet_wrap(. ~ Scale + Stage, scales="free", nrow=5, ncol=3) +
    labs(title = as.character(tax[spp,"species"]),
            subtitle = "Regression summaries by scale and stage",
            caption = "Model stages: local + space + landscape.",
            x = "Predicted",
            y = "Observed") +
    theme_new +
    scale_color_manual(values=bc[length(sc1):1]) +
    geom_abline(intercept=0, slope = 1, col=1, lty=2)
  p40
  ggsave(file.path(ROOT, "species", paste0("all-", spp, ".png")), p40,
         height=10, width=7)
}

for (spp in SPP) {
    cat(spp, "\n")
    vals <- NULL
    for (scale in sc1) {
      tmp <- data.frame(ALL[[spp]][["HF"]][[scale]]$results)
      tmp$scale <- factor(paste0(scale, "x", scale),
                          sc2)
      vals <- rbind(vals, tmp)
    }

    s <- t(sapply(ALL[[spp]][["HF"]], "[[", "stats"))
    COR <- s[,"COR"][1:length(sc1)]
    R2 <- s[,"R2"][1:length(sc1)]
    STATS <- data.frame(Value=c(COR, COR, R2, R2),
        Scale=sc2,
        Stat=factor(rep(c("r", "r", "R2", "R2"), each=length(sc1)), c("r", "R2")),
        Type=factor(rep(c("Within","Between","Within","Between"), each=length(sc1)),
            c("Between", "Within")))
    STATS$What <- "Actual"

    tmp <- t(sapply(ALLW[[spp]], function(z) do.call(rbind, z)[,"COR"]))
    tmp <- tmp[,sc1]
    QCORW <- rbind(Value=COR, apply(tmp, 2, quantile, c(0.5, 0.05, 0.95)))
    tmp <- data.frame(Value=as.numeric(tmp),
        Scale=rep(sc2, each=nrow(tmp)))
    CORW <- tmp
    CORW$Stat <- factor("r", c("r", "R2"))
    CORW$Type <- factor("Within", c("Between", "Within"))

    tmp <- t(sapply(ALLB[[spp]], function(z) do.call(rbind, z)[,"COR"]))
    tmp <- tmp[,sc1]
    QCORB <- rbind(Value=COR, apply(tmp, 2, quantile, c(0.5, 0.05, 0.95)))
    tmp <- data.frame(Value=as.numeric(tmp),
        Scale=rep(sc2, each=nrow(tmp)))
    CORB <- tmp
    CORB$Stat <- factor("r", c("r", "R2"))
    CORB$Type <- factor("Between", c("Between", "Within"))

    tmp <- t(sapply(ALLW[[spp]], function(z) do.call(rbind, z)[,"R2"]))
    tmp <- tmp[,sc1]
    QR2W <- rbind(Value=R2, apply(tmp, 2, quantile, c(0.5, 0.05, 0.95)))
    tmp <- data.frame(Value=as.numeric(tmp),
        Scale=rep(sc2, each=nrow(tmp)))
    R2W <- tmp
    R2W$Stat <- factor("R2", c("r", "R2"))
    R2W$Type <- factor("Within", c("Between", "Within"))

    tmp <- t(sapply(ALLB[[spp]], function(z) do.call(rbind, z)[,"R2"]))
    tmp <- tmp[,sc1]
    QR2B <- rbind(Value=R2, apply(tmp, 2, quantile, c(0.5, 0.05, 0.95)))
    tmp <- data.frame(Value=as.numeric(tmp),
        Scale=rep(sc2, each=nrow(tmp)))
    R2B <- tmp
    R2B$Stat <- factor("R2", c("r", "R2"))
    R2B$Type <- factor("Between", c("Between", "Within"))

    V <- rbind(CORW, CORB, R2W, R2B)
    V$What <- "Randomized"

if (FALSE) {
    p4 <- ggplot(vals,
                aes(x = yhat, y = yobs, group=scale)) +
        geom_point(color=ifelse(vals$inPI==0,
                                paste0(bc[1],"80"), paste0(bc[5],"80"))) +
        geom_smooth(method="lm", color=bc[8]) +
        facet_wrap(vars(scale), scales="free") +
        labs(title = as.character(tax[spp,"species"]),
            subtitle = "Regression summaries by scale",
            caption = "Model stages: local + space + landscape.",
            x = "Predicted",
            y = "Observed") +
        theme_new +
        scale_color_manual(values=bc[length(sc1):1]) +
        geom_abline(intercept=0, slope = 1, col=1, lty=2)
    p4
    ggsave(file.path(ROOT, "species", paste0("regr-", spp, ".png")), p4)

    p5 <- ggplot(V,
                aes(x = Scale, y = Value, fill=Scale, color=Scale)) +
        geom_violin() +
        facet_grid(rows=vars(Stat), cols=vars(Type), scales="free") +
        labs(title = as.character(tax[spp,"species"]),
            subtitle = "Randomization summaries by scale",
            caption = "Model stages: local + space + landscape.",
            x = "Scale",
            y = "Value") +
        theme_new +
        scale_color_manual(values=bc[length(sc1):1], guide=FALSE) +
        scale_fill_manual(values=bc[length(sc1):1], guide=FALSE) +
        geom_point(data=STATS, aes(y=Value), col=1)
    p5
    ggsave(file.path(ROOT, "species", paste0("rand-", spp, ".png")), p5)
}

    if (spp == SPP[1]) {
        ALLRES <- rbind(STATS, V)
        ALLRES[[spp]] <- ALLRES$Value
        ALLRES <- ALLRES[,c("What", "Scale", "Stat", "Type")]
        ALLRES$What <- as.factor(ALLRES$What)
        QQ <- rbind(stdfun(QCORW, "Within", "r"), stdfun(QCORB, "Between", "r"),
                    stdfun(QR2W, "Within", "R2"), stdfun(QR2B, "Between", "R2"))
    }else {
        ALLRES[[spp]] <- rbind(STATS, V)$Value
        QQ <- rbind(QQ,
                    stdfun(QCORW, "Within", "r"), stdfun(QCORB, "Between", "r"),
                    stdfun(QR2W, "Within", "R2"), stdfun(QR2B, "Between", "R2"))
    }

}

write.csv(ALLRES, row.names = FALSE, file=file.path(ROOT, "Allres.csv"))
write.csv(QQ, row.names = FALSE, file=file.path(ROOT, "RndSignif.csv"))

library(magick)
for (spp in SPP) {
    cat(spp, "\n")
    f1 <- file.path(ROOT, "species", paste0("regr-", spp, ".png"))
    f2 <- file.path(ROOT, "species", paste0("rand-", spp, ".png"))
    i1 <- image_read(f1)
    i2 <- image_read(f2)
#    i1 <- image_resize(i1, "1000x")
#    i2 <- image_resize(i2, "1000x")
    i12 <- image_append(c(i1, i2), stack=TRUE)
    i12 <- image_resize(i12, "x2000")
    image_write(i12, file.path(ROOT, "species", paste0("comb-", spp, ".png")))
}

for (spp in SPP) {
    cat(spp, "\n")
    f1 <- file.path(ROOT, "species", paste0("all-", spp, ".png"))
    f2 <- file.path(ROOT, "species", paste0("rand-", spp, ".png"))
    i1 <- image_read(f1)
    i2 <- image_read(f2)
    i1 <- image_resize(i1, "1000x")
    i2 <- image_resize(i2, "1000x")
    image_write(i1, file.path(ROOT, "species", paste0("all-", spp, ".png")))
    image_write(i2, file.path(ROOT, "species", paste0("rand-", spp, ".png")))
}

for (spp in SPP) {
  cat(spp, "\n")

  scale <- "10"
  vals <- data.frame(ALL[[spp]][["HF"]][[scale]]$results)
  vals$scale <- factor(paste0(scale, "x", scale),
                          sc2)
  vals$Footprint <- df1[rownames(vals), "TotalHF"]

  p4 <- ggplot(vals,
                aes(x = yhat, y = yobs)) +
        geom_point(aes(color=Footprint, size=Footprint)) +
        geom_smooth(method="lm") +
        labs(title = as.character(tax[spp,"species"]),
            subtitle = "Regression summaries by scale",
            caption = "Model stages: local + space + landscape.",
            x = "Predicted",
            y = "Observed") +
        theme_new +
        #scale_color_manual(values=bc[5:1]) +
        geom_abline(intercept=0, slope = 1, col=1, lty=2)
    p4
    ggsave(file.path(ROOT, "species", paste0("thf-", spp, ".png")), p4)
}

px <- ggplot(df1,
    aes(x = Energy, y = Forestry, color=RuralUrban, size=Transportation)) +
    geom_point() +
    labs(title = "Footprint %",
            subtitle = "By industrial sectors",
            caption = "Representing 15 big grids.") +
    theme_new +
    scale_color_gradientn(colours=bc[5:1])
px
ggsave(file.path(ROOT, "figX-hf.png"), px, width=6, height=5)


table(QQ$Scale, QQ$Signif)

QQ$y <- QQ$Signif==1
mm <- glm(y ~ (Scale + Stat + Type)^2, data=QQ, family=binomial)
mm <- step(mm)
summary(mm)

nd <- QQ[QQ$Species=="OVEN",c("Scale", "Stat", "Type")]
nd$p <- predict(mm, newdata=nd, type="response")

p6 <- ggplot(nd, aes(x=Scale, y=p, color=Type, group=Type)) +
  geom_line(lwd=1.5) +
  facet_wrap(vars(Stat)) +
        labs(title = "Significance",
            subtitle = "Statistics higher than randomized",
            caption = "Model stages: local + space + landscape.",
            x = "Scale",
            y = "P(Stat > Random)") +
        theme_new +
        scale_color_manual(values=bc[c(5,1)]) +
        scale_fill_manual(values=bc[c(5,1)])
p6
ggsave(file.path(ROOT, "fig6-prob-rnd.png"), p6, width=7, height=5)


## fancy graphics for illustration



elmat %>%
  sphere_shade(texture = "desert") %>%
  add_water(detect_water(elmat), color = "desert") %>%
  add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
  add_shadow(ambient_shade(elmat), 0) %>%
  plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
Sys.sleep(0.2)
render_snapshot()


library(leaflet)

df <- as(spTransform(xy,
    CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    "SpatialPointsDataFrame")
df@data <- vv[rownames(coordinates(df)),]
leaflet(df) %>% addTiles() %>% addCircleMarkers()

m <- leaflet(df) %>%
  addTiles() %>%
  #fitBounds(0, 40, 10, 50) %>%
  #setView(-110.4768, 54.56633, zoom = 12) %>%
  addTiles(urlTemplate = "https://mts1.google.com/vt/lyrs=s&hl=en&src=app&x={x}&y={y}&z={z}&s=G", attribution = 'Google') %>%
  addCircleMarkers(
    radius = 2,
    color="red",
    stroke = FALSE, fillOpacity = 0.6)
m

