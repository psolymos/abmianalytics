library(opticut)
library(mefa4)

# simulate some data
n <- 1000 # sample size
K <- 10 # number of categories in x

x <- sample(LETTERS[1:K], n, replace=TRUE)
table(x)

lambda <- rexp(K)
names(lambda) <- LETTERS[1:K]
y <- rpois(n, lambda[x])

boxplot(y ~ x)
aggregate(y, list(x), mean)

# GLM
d <- data.frame(x=as.factor(x), y=y)
m <- glm(y ~ x, data=d, family=poisson)
summary(m)

## estimate lambda
nd <- data.frame(x=factor(levels(d$x), levels(d$x)))
nd$lambda <- lambda[levels(d$x)]
nd$lambda_hat <- predict(m, nd, type="response")
nd

## Lorenz-tangent approach to binarize a multi-level problem

## predict for each sample
fit <- fitted(m)
g <- sum_by(fit, d$x)
g # column 'x' is the sum, colum 'by' is the number of cases

## make a Lorenz curve out of this: mean vs size
## note: using the mean instead of the sum because
## not trying to identify the samples where abundance is concentrated
## but rather to find the largest differential in the ordered estimates
l <- lorenz(g[,"x"]/g[,"by"], g[,"by"])
l # see ?lorenz for explanation of terms

s <- summary(l)
## labels of 'good' classes
lab <- rownames(l[l[,"x"] >= s["x[t]"],])
lab

plot(l)
abline(0, 1, lty=2)
lines(rep(s["p[t]"], 2), c(s["p[t]"], s["L[t]"]), col=2)

