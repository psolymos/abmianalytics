
# non-stratified sampling
library(mefa4)
set.seed(1)
K <- 5
N <- 1000
#Max <- runif(K, 0.6, 0.9)
#prop <- sapply(1:K, function(i) runif(N, 0.6, Max[i]))
prop <- t(sapply(1:N, function(i) rbeta(K, 0.2, 4)))
prop <- prop / rowSums(prop)
colnames(prop) <- paste0("H", 1:K)
summary(prop)
iv <- find_max(prop)
summary(iv)
iv$w <- 1- iv$value

Lam <- sort(rexp(K, 2))
names(Lam) <- paste0("indexH", 1:K)
Abu <- colSums(t(prop) * Lam)
Y <- rpois(N, Abu)
table(Y)

m0 <- glm(Y ~ index - 1, data=iv, family=poisson)
summary(m0)
(est0 <- exp(coef(m0))[match(names(Lam), names(coef(m0)))])
est0[is.na(est0)] <- 0
names(est0) <- names(Lam)

m1 <- glm(Y ~ index - 1, data=iv, family=poisson, weights=iv$w)
summary(m1)
(est1 <- exp(coef(m1))[match(names(Lam), names(coef(m1)))])
est1[is.na(est1)] <- 0
names(est1) <- names(Lam)

cbind(Lam=Lam, Hab=est0, Cont=est1)
