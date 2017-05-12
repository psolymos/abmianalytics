n <- 100

set.seed(1234)
Y <- rnorm(n, 2, 0.5)

fun <- function(p) {
    if (p<2)
        stop("p<2")
    X <- cbind(1, replicate(p-1, rnorm(n)))
    b <- lm.fit(X, Y)$coefficients
    yhat <- drop(X %*% b)
    SST <- sum((Y - mean(Y))^2)
    SSE <- sum((Y - yhat)^2)
    R2 <- 1-(SSE/SST)
    r2 <- cor(cbind(Y,yhat))[1,2]
    c(R2=R2, r2=r2)
}
fun2 <- function(p, R=100) {
    res <- t(replicate(R, fun(p)))
    x <- mean(res[,"R2"])
    c(p=p, 
        R=R, 
        R2mean = x,
        R2var = sd(res[,"R2"])^2,
        R2ex = p/(n-1),
        R2adj = x - (1-x) * (p/(n-p-1)))

}
res <- t(sapply(2^(1:6), fun2, R=200))

summary(res)
par(mfrow=c(1,4))
plot(res[,c(1,3)], type="b")
plot(res[,c(1,4)], type="b")
plot(res[,c(1,5)], type="b")
plot(res[,c(1,6)], type="b")
