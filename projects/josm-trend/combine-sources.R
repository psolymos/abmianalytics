library(mefa4)
ROOT <- "~/Dropbox/josm/2017"

tax <- read.csv(file.path(ROOT, "Table-tax.csv"))
rownames(tax) <- tax$AOUcode
yr <- read.csv(file.path(ROOT, "Table-yr.csv"))
rownames(yr) <- yr$AOUcode
bbs1 <- read.csv(file.path(ROOT, "Table-bbs.csv"))
rownames(bbs1) <- bbs1$AOUcode
res1 <- read.csv(file.path(ROOT, "Table-res1.csv"))
rownames(res1) <- res1$AOUcode
res2 <- read.csv(file.path(ROOT, "Table-res2.csv"))
rownames(res2) <- res2$AOUcode
dh <- read.csv(file.path(ROOT, "Dave-numbers.csv"))
dh <- droplevels(dh[dh$AOUCode != "",])
rownames(dh) <- dh$AOUCode

bbs2 <- read.csv(file.path(ROOT, "ron_bbs_tc_e20170828-ABbcr6.csv"))
compare_sets(bbs2$commonName, tax$English_Name)
bbs2$AOUcode <- rownames(tax)[match(bbs2$commonName, tax$English_Name)]
bbs2 <- droplevels(bbs2[!is.na(bbs2$AOUcode),])
bbs2lt <- bbs2[bbs2$timeFrame=="long-term",]
rownames(bbs2lt) <- bbs2lt$AOUcode
bbs2st <- bbs2[bbs2$timeFrame=="short-term",]
rownames(bbs2st) <- bbs2st$AOUcode

SPP <- as.character(res1$AOUcode)
tax <- droplevels(tax[SPP,])
yr <- droplevels(yr[SPP,])
bbs1 <- droplevels(bbs1[SPP,])
bbs2lt <- droplevels(bbs2lt[SPP,])
bbs2st <- droplevels(bbs2st[SPP,])
dh <- droplevels(dh[SPP,])
res1 <- droplevels(res1[SPP,])
res2 <- droplevels(res2[SPP,])


x <- data.frame(
    Residual=round(res1$offroad_Dhf_Median, 3),
    Revisit=dh$Median,
    BBSst12=bbs1$ShortTerm.annualTrend,
    BBSlt12=bbs1$LongTerm.annualTrend,
    BBSst15=bbs2st$annualTrend,
    BBSlt15=bbs2lt$annualTrend)
rownames(x) <- SPP
x <- x[rowSums(is.na(x))==0,]
x <- x[apply(abs(x), 1, max) < 20,]
#write.csv(x, file.path(ROOT, "all.csv"))

pif <- read.csv("c:/Users/Peter/Dropbox/bam/PIF-AB/qpad-pif-results.csv")
rownames(pif) <- pif$Species_ID
pif <- droplevels(pif[rownames(x),])

f <- function(B) 1+B/100
xx <- data.frame(R=(f(x$Revisit) / f(x$BBSst15)),
    Roadside=pif$DeltaR, Habitat=1/pif$H)

Pairs <- function (x, ...)
{
    #y <- data.frame(x)
    y <- x
    fun.lower <- function(x1, x2, ...) {
        COR <- cor(x1, x2)
        text(mean(range(x1, na.rm = TRUE)), mean(range(x2, na.rm = TRUE)),
            round(COR, 2), cex = 0.5 + 2*abs(COR))
        box(col="grey")
    }
    fun.upper <- function(x1, x2, ...) {
        if (is.factor(x1)) {
            x1 <- as.integer(x1)
        }
        if (is.factor(x2)) {
            x1 <- as.integer(x2)
        }
        abline(h=0, v=0, lty=1, col="grey")
        points(x1, x2, col = "#2C7BB680", pch=19)
        LM <- lm(x2 ~ x1)
        abline(LM, col="#D7191C")
        box(col="#80B9D8")
    }
    panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$density
        Max <- max(y)
        y <- y/Max
#        rect(breaks[-nB], 0, breaks[-1], y, col = "#FDC980", border = "#F07C4A",
#            ...)
        den <- density(x)
        den$y <- den$y/Max
        polygon(den, col="#F07C4A80", border="#F07C4A")
        abline(v=0, lty=1, col="grey")
        box(col="#F07C4A")
    }
    pairs.default(y, lower.panel = fun.lower, upper.panel = fun.upper,
        diag.panel = panel.hist)
    invisible(NULL)
}

Pairs(x)
Pairs(xx[rowSums(is.na(xx))==0,])

