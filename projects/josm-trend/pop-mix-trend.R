## simulating equal numbers of sites for 2 populations (same species) with
## - different trend
## - different temporal variability
## then we sample different proportions of the 2 populations
## and see what kind of trend do we get

set.seed(1)

N <- 100
T <- 25
TT <- rep(seq_len(T), each=N)
NN <- rep(seq_len(N), T)

alpha1 <- 1
sigma1 <- 1
alpha2 <- 1
sigma2 <- 1

beta1 <- log(-5/100 + 1)
beta2 <- log(5/100 + 1)

tau1 <- 0.2
tau2 <- 0.2

eps1 <- rnorm(N, 0, sigma1)
eps2 <- rnorm(N, 0, sigma2)

nu1 <- rnorm(T, 0, tau1)
nu2 <- rnorm(T, 0, tau2)

lam1 <- exp(alpha1 + beta1*(TT-1) + eps1[NN] + nu1[TT])
lam2 <- exp(alpha2 + beta2*(TT-1) + eps2[NN] + nu2[TT])

lam1mat <- matrix(lam1, N, T)
lam2mat <- matrix(lam2, N, T)

p <- 0.5
pp <- round(N*p)
lam3mat <- rbind(lam1mat[1:pp,,drop=FALSE], lam2mat[(N-pp+1):N,,drop=FALSE])

plot(seq_len(T), colMeans(lam1mat), type="l", col=2,
    ylim=c(0, max(colMeans(lam1mat), colMeans(lam2mat))), lty=2,
    xlab="Years", ylab="Mean abundance")
abline(lm(colMeans(lam1mat) ~ seq_len(T)), col=2)
lines(seq_len(T), colMeans(lam2mat), col=4, lty=2)
abline(lm(colMeans(lam2mat) ~ seq_len(T)), col=4)
lines(seq_len(T), colMeans(lam3mat), col=3, lty=2)
abline(lm(colMeans(lam3mat) ~ seq_len(T)), col=3)
legend("topleft", bty="n", lty=1, col=c(2,4,3),
    legend=c("Pop #1", "Pop #2", "Mix"))

