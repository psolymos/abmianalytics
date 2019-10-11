library(epiR)

n <- 20
set.seed(1)
x <- runif(n)
e <- rnorm(n)
s <- runif(n)

f <- function(v=0.05, a=0, b=1, w=0) {
    exp(log((1-w)*x+w*s) + e*v) * b + a
}

vals <- expand.grid(
    v=c(0.05, 0.25, 0.5), # error
    a=c(0, 0.2), # intercept
    b=c(1, 0.5), # slope
    w=c(0, 0.5, 1)) # association

g <- function(v=0.05, a=0, b=1, w=0) {
    plot(x, f(v=v, a=a, b=b, w=w))
    abline(0, 1)
    invisible()
}

h <- function(x, y) {
    oc <- epiR::epi.occc(cbind(x, y))[1:3]
    names(oc) <- c("oc", "pr", "ac")
    oc$rs <- cor(x, y, method="spearman")
    oc$ru <- sum(x*y)/(sqrt(sum(x^2)*sum(y^2)))
    oc$r2 <- summary(lm(y ~ x))$r.squared
    unlist(oc)
}

res <- t(apply(vals, 1, function(z) h(x, f(v=z[1], a=z[2], b=z[3], w=z[4]))))


g2 <- function(i) {
    v <- vals$v[i]
    a <- vals$a[i]
    b <- vals$b[i]
    w <- vals$w[i]
    y <- f(v=v, a=a, b=b, w=w)
    plot(x, y, xlim=c(0, max(x, y)), ylim=c(0, max(x, y)),
        main=paste0(paste0(c("v", "a", "b", "w"), "=", c(v, a, b, w)), collapse=", "),
        sub=paste0(colnames(res), "=", round(res[i,], 2), collapse=", "))
    abline(0, 1, lty=2)
    abline(lm(y ~ x), col=2)
    abline(lm(y ~ x-1), col=2, lty=2)
    invisible()
}

pdf("assess.pdf", onefile=TRUE)
for (i in order(res[,1]))
    g2(i)
dev.off()

plot(data.frame(res))
round(cor(res), 2)
