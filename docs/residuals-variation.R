#' ---
#' title: "Demo of residual variation"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "`r as.Date(Sys.time())`"
#' output: pdf_document
#' ---
#'
n <- c(10, 100, 500, 1000, 2000, 5000, 10000)
mu <- 2
sig <- sqrt(0.5)

y1 <- lapply(n, rnorm, mean=mu, sd=sig)

mu_hat <- sapply(y1, mean)
sig_hat <- sapply(y1, sd)

plot(log(n), mu_hat)
abline(h=mu)

plot(log(n), sig_hat)
abline(h=sig)

# repeat this many times and compare sqrt(n)






# binomial weight

set.seed(1)
p <- 0.32
n <- 100
N <- 1
Y <- rbinom(n, N, p)
#W <- rep(1, n)
W <- rep(c(1,0.5), each=n/2)
pval <- seq(0,1,0.001)

nll <- function(p, w=1) {
    -sum(w * dbinom(Y, N, p, log = TRUE))
}

## no weights, response as 2 cols
o <- optimize(nll, c(0,1), w=W)
m <- glm(cbind(Y, N-Y) ~ 1, family=binomial, weights=W)
m <- glm(Y ~ 1, family=binomial, weights=W)

plot(pval, exp(-sapply(pval, nll)), type="l", xlab="p", ylab="L(Y;p)")
abline(v=o$minimum)

plogis(coef(m))
o$minimum

logLik(m)
-o$objective

## response as prop, weights=W
W <- rep(4, n)
o <- optimize(nll, c(0,1), w=W)
m <- glm((Y/N) ~ 1, family=binomial, weights=W)

plogis(coef(m))
o$minimum

logLik(m)
-o$objective

##
W1 <- rep(4, n)
W2 <- rep(c(1,0.5), each=n/2)
o <- optimize(nll, c(0,1), w=W2)
m <- glm((Y/N) ~ 1, family=binomial, weights=W1*W2)

plogis(coef(m))
o$minimum

logLik(m)
-o$objective

