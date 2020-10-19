library(mefa4)
library(mgcv)
load("d:/abmi/AB_data_v2020/data/analysis/species/birds/WildTrax_2015-2019_2020-09-23.RData")


y <- as.matrix(xt1min)
x <- xx1min
rm(list=setdiff(ls(), c("x", "y")))


str(x)

plot(x$ToY, x$ToD)
ss <- x$ToD < 12 & x$ToD > 2
table(ss)
with(x[ss,], plot(ToY, ToD))

x <- x[ss,]
y0 <- y[ss,]
y0[y0 > 0] <- 1
y <- y0[,colSums(y0) >= 100]

sort(colSums(y))

px <- data.frame(ToY=seq(min(x$ToY), max(x$ToY), 1))

f <- function(i, k=-1) {
    y01 <- y[,i]
    m <- mgcv::gam(y01 ~ s(ToY, k=k), data=x, family="binomial")
    predict(m, newdata=px, type="response")
}

py <- pbapply::pbsapply(colnames(y), f)
pym <- apply(py, 2, function(z) z/max(z))

matplot(px$ToY, py, type="l", lty=1, col="#00000044")


plot(px$ToY, py[,"BOCH"], ylim=c(0,1), type="l")
plot(px$ToY, py[,"TEWA"], ylim=c(0,1), type="l")
