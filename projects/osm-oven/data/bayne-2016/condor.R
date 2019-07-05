library(intrval)
x <- read.csv("data/bayne-2016/condor-suppl.csv", stringsAsFactors=FALSE)
x1 <- x[,1:4]
x1[] <- lapply(x1, as.factor)
rownames(x1) <- x1$AOUcode
x2 <- x[,5:13]
x3 <- x[,14:22]

f1 <- function(i) {
    z <- t(sapply(strsplit(x2[,i], "/"), as.integer))
    colnames(z) <- paste0(colnames(x2)[i], c("_Impact", "_Control"))
    z
}
x2 <- do.call(cbind, lapply(1:ncol(x2), f1))
rownames(x2) <- x1$AOUcode

f21 <- function(z) {
    if (z == "")
        return(c(NA, NA, NA))
    zz <- strsplit(z, " ")
    zzz <- strsplit(zz[[1]][2], "-")
    as.numeric(gsub("[^0-9\\.]", "", c(zz[[1]][1], zzz[[1]])))
}
f2 <- function(i) {
    z <- t(sapply(x3[,i], f21))
    colnames(z) <- paste0(colnames(x3)[i], c("", "_Lo", "_Hi"))
    z
}
x3 <- do.call(cbind, lapply(1:ncol(x3), f2))
rownames(x3) <- x1$AOUcode
xx <- data.frame(x1, x2, x3)

spp <- "OVEN"

m <- matrix(x3[spp,], nrow=3)
m0 <- matrix(m[1,], 3, 3)
mL <- matrix(m[2,], 3, 3)
mU <- matrix(m[3,], 3, 3)
main <- paste(x1[spp, "CommonName"], x1[spp, "Habitat"])
colnames(m0) <- c("Seismic", "Pipeline", "Wellpad")
rownames(m0) <- c("50m", "100m", "Unlimited")
dimnames(mL) <- dimnames(mU) <- dimnames(m0)
m0[is.na(m0)] <- 0
mL[is.na(mL)] <- 0
mU[is.na(mU)] <- 0
M <- list(est=m0, lwr=mL, upr=mU)
dput(M)


v <- barplot(m0,beside=TRUE, main=main, col=col, ylim=c(0, min(10, max(mU))))
abline(h=1, col='#8da0cb')
for (k in 1:9) {
    int <- c(mL[k], mU[k])
    lwd <- if (1 %[]% int) 1 else 3
    lines(c(v[k], v[k]), int, lwd=lwd)
}


