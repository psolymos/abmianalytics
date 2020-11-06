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
