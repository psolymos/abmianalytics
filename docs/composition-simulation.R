library(mefa4)
library(parallel)
set.seed(1)

m <- 1000

lamA <- 1
lamB <- 2

#p <- c(runif(m, 0, 0.2), runif(m, 0.2, 0.8), runif(m3, 0.8, 1))
#p <- c(runif(m, 0, 0.2), runif(m3, 0.8, 1))

o <- 0.1

f <- function(o=0, m=500, lamA=0.2, lamB=0.4) {
    p <- c(runif(m, 0, o), runif(m, 1-o, 1))
    z <- find_max(cbind(A=(1-p),B=p))
    x <- z$index
    w1 <- z$value
    w2 <- pmax(0, pmin(1, 2*z$value-0.5))

    lam <- (1-p)*lamA + p*lamB
    y <- rpois(length(p), lam)

    m <- list(
        m1 = glm(y ~ p, family=poisson),
        m2 = glm(y ~ x, family=poisson),
        m3 = glm(y ~ x, family=poisson, weights=w1),
        m4 = glm(y ~ x, family=poisson, weights=w2))

    cf <- sapply(m, coef)
    th <- rbind(lam1=exp(cf[1,]), lam2=exp(colSums(cf)))
    df <- th - c(lamA, lamB)
    list(cf=cf, th=th, df=df)
}


B <- 100
cl <- makeCluster(8)

res000 <- pbapply::pbreplicate(B, f(0), simplify = FALSE)
res005 <- pbapply::pbreplicate(B, f(0.05), simplify = FALSE)
res010 <- pbapply::pbreplicate(B, f(0.1), simplify = FALSE)
res020 <- pbapply::pbreplicate(B, f(0.2), simplify = FALSE)
res050 <- pbapply::pbreplicate(B, f(0.5), simplify = FALSE)

stopCluster(cl)

op <- par(mfrow=c(2,3))
boxplot(t(sapply(res000, function(z) as.numeric(z$df))), ylim=c(-0.5,0.5));abline(h=0)
boxplot(t(sapply(res005, function(z) as.numeric(z$df))), ylim=c(-0.5,0.5));abline(h=0)
boxplot(t(sapply(res010, function(z) as.numeric(z$df))), ylim=c(-0.5,0.5));abline(h=0)
boxplot(t(sapply(res020, function(z) as.numeric(z$df))), ylim=c(-0.5,0.5));abline(h=0)
boxplot(t(sapply(res050, function(z) as.numeric(z$df))), ylim=c(-0.5,0.5));abline(h=0)
par(op)
