library(opticut)
library(cure4insect)
opar <- set_options(path = "w:/reports")
load_common_data()

SPP <- get_all_species("birds")
subset_common_data(id=NULL, species=SPP)

level <- 0.8

res <- list()
for (spp in SPP) {
    cat(spp, "\n")
    y <- load_species_data(spp, boot=FALSE)
    lc0 <- lorenz(rowSums(y$SA.Ref))
    lc1 <- lorenz(rowSums(y$SA.Curr))
    xt0 <- quantile(lc0, 1-level, type="L")
    pt0 <- iquantile(lc0, xt0, type="p")
    xt1 <- quantile(lc1, 1-level, type="L")
    pt1 <- iquantile(lc1, xt1, type="p")
    res[[spp]] <- rbind(
        ref= c(unclass(summary(lc0)), pt0=unname(pt0)),
        curr=c(unclass(summary(lc1)), pt1=unname(pt1)))
}

pts <- t(sapply(res, function(z) z[,"pt0"]))
plot(pts, xlim=c(0.4,1), ylim=c(0.4,1))

ind <- (1-pts)/(level)
summary(t(sapply(res,function(z) z[,"G"])))


res0 <- list()
res1 <- list()
q <- c(0.05, 0.25, 0.5, 0.75, 0.95)
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    y <- load_species_data(spp, boot=FALSE)
    r <- rasterize_results(y)
    D <- values(r[["NC"]])
    D <- D[!is.na(D)]
    lc <- lorenz(D)
    Dmax <- max(D)
    xt <- quantile(lc, 1-level, type="L")
    pt <- iquantile(lc, xt, type="p")
    Lt_half <- iquantile(lc, Dmax*q, type="L")
    pt_half <- iquantile(lc, Dmax*q, type="p")
    tmp <- rbind(Lt=Lt_half, pt=pt_half)
    colnames(tmp) <- q
    res0[[spp]] <- c(unclass(summary(lc)), pt=unname(pt))
    res1[[spp]] <- tmp
}

plot(lc)
abline(h=Lt_half)
abline(v=pt_half)

pc <- t(sapply(res1, function(z) z[,"0.5"]))
