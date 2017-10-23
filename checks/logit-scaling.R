x1 <- 0:100/100
x2 <- 1-x1

beta0 <- 1
beta1 <- +1
beta2 <- -1
betaC <- -0.5

P1 <- plogis(beta0 + beta1*x1 + beta2*x2 + betaC)
P2 <- plogis(beta0 + beta1*1 + betaC)*x1 + plogis(beta0 + beta2*1 + betaC)*x2

plot(P1, P2,type="l",col=2,main=paste("Mean P =", round(mean(P1),3)))
abline(0,1,lty=2)

f <- function(beta0=0, ...) {
    P1 <- plogis(beta0 + beta1*x1 + beta2*x2 + betaC)
    P2 <- plogis(beta0 + beta1*1 + betaC)*x1 + plogis(beta0 + beta2*1 + betaC)*x2
    lines(P1, P2, ...)
    invisible(cbind(P1, P2))
}

beta1 <- +1
beta2 <- -1
betaC <- -0.5
plot(0, type="n", ylim=c(0, 1), xlim=c(0, 1), xlab="Method 1", ylab="Method 2")
abline(0,1,lty=2)
f(-1, col=2)
f(0, col=3)
f(1, col=4)
