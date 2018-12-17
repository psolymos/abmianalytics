library(mefa4)
library(intrval)
source("~/repos/abmianalytics/birds/00-functions.R")


## validation
ROOT <- "d:/abmi/AB_data_v2018/data/analysis/birds"
PROJ <- "validation"

load(file.path(ROOT, "data", "ab-birds-validation-2018-12-07.RData"))
X <- get_model_matrix(DAT, mods)

if (FALSE) {
    rn <- rownames(DAT)[DAT$ABMIsite!=""]
    DATv <- DAT[rn,]
    YYv <- YY[rn,]
    OFFv <- OFF[rn,]
    SSHv <- SSH[rn,]
}

Xv <- get_model_matrix(DATv, mods)

spp <- "WTSP"
spp <- "OVEN"

names(mods)
stage <- "Water"
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
est <- get_coef(res, X, stage=stage, na.out=FALSE)
stopifnot(all(colnames(X)==colnames(est)))

mui <- X %*% t(est)
pri <- apply(exp(mui), 1, median)
offi <- OFF[,spp]
evi <- apply(exp(mui+offi), 1, median)
yi <- as.numeric(YY[,spp])
table(yi)

muo <- Xv %*% t(est)
pro <- apply(exp(muo), 1, median)
offo <- OFFv[,spp]
evo <- apply(exp(muo+offo), 1, median)
yo <- as.numeric(YYv[,spp])
table(yo)

(AUCi <- simple_auc(simple_roc(ifelse(yi>0, 1, 0), evi)))
(AUCo <- simple_auc(simple_roc(ifelse(yo>0, 1, 0), evo)))

MAX <- max(quantile(evi, 0.99), quantile(evo, 0.99))
boxplot(evi ~ yi, main=spp, sub="internal", xlab="Observed count", ylab="Expected value", ylim=c(0,MAX))
boxplot(evo ~ yo, main=spp, sub="internal", xlab="Observed count", ylab="Expected value", ylim=c(0,MAX))


boxplot(pri ~ DAT$vegc, cex.axis=0.6)

Predm <- groupSums(exp(muo+offo), 1, DATv$ABMIsite)
Y <- sum_by(yo, DATv$ABMIsite)
all(rownames(Predm)==rownames(Y))

q <- 1
y9obs <- Y[rownames(Y) != "","x"]
y9pred <- apply(Predm[rownames(Predm) != "",], 1, median)
MAX <- max(quantile(y9obs, q), quantile(y9pred, q))

plot(y9obs, y9pred, xlim=c(0,MAX), ylim=c(0,MAX))
abline(0,1)
abline(lm(y9pred ~ y9obs), col=2)
abline(lm(y9pred ~ y9obs-1), col=2)

cor(y9obs, y9pred, method="spearman")




validate <- function(res, Xv, stage=NULL, offset=TRUE) {
    est <- get_coef(res, Xv, stage=stage, na.out=FALSE)
    c1 <- colSums(abs(est)) > 0
    if (any(c1[c("SSH_KM", "SSH05_KM")])) {
        essh <- est[,c("SSH_KM", "SSH05_KM")]
        c1[c("SSH_KM", "SSH05_KM")] <- FALSE # drop SSH
        muo <- Xv[,c1,drop=FALSE] %*% t(est[,c1,drop=FALSE])
        mussh <- muo
        mussh[] <- 0 # put HHS effects here
        for (i in seq_len(nrow(est))) {
            ssh <- res[[i]]$ssh
            v <- rowSums(SSHv[,ssh$labels])
            mussh[,i] <- essh[1,"SSH_KM"]*v + essh[1,"SSH05_KM"]*sqrt(v)
        }
        muo <- muo + mussh # add them up
    } else {
        muo <- Xv[,c1,drop=FALSE] %*% t(est[,c1,drop=FALSE])
    }
    pro <- apply(exp(muo), 1, median)
    offo <- if (offset)
        OFFv[,spp] else 0
    evo <- apply(exp(muo+offo), 1, median)
    yo <- as.numeric(YYv[,spp])
    Predm <- groupSums(exp(muo+offo), 1, DATv$ABMIsite)
    Y <- sum_by(yo, DATv$ABMIsite)
    stopifnot(all(rownames(Predm)==rownames(Y)))
    y9obs <- Y[rownames(Y) != "","x"]
    y9pred <- apply(Predm[rownames(Predm) != "",], 1, median)
    AUCo <- simple_auc(simple_roc(ifelse(yo>0, 1, 0), evo))
    CORo <- cor(y9obs, y9pred, method="spearman")
    list(AUC=AUCo, COR9=CORo,
        y1obs=yo, lam1=evo, y9obs=y9obs, lam9=y9pred)
}
plotOne <- function(spp, q=1, All) {
    One <- All[[spp]]
    lam1 <- One$HF$lam1
    y1 <- One$HF$y1obs
    lam9 <- One$HF$lam9
    y9 <- One$HF$y9obs

    op <- par(mfrow=c(2,2))
    plot(0:9, sapply(One, "[[", "COR9"), type="l", col="#0000FF80", ylim=c(-1,1),
        xlab="Model stages", ylab="Correlation", main=spp, lwd=2)
    abline(h=0, lty=2)
    text(0:9, rep(-0.1, 10), round(sapply(One, "[[", "COR9"), 3), cex=0.6)
    plot(0:9, sapply(One, "[[", "AUC"), type="l", col="#0000FF80", ylim=c(0,1),
        xlab="Model stages", ylab="AUC", lwd=2)
    abline(h=0.5, lty=2)
    text(0:9, rep(0.45, 10), round(sapply(One, "[[", "AUC"), 3), cex=0.6)

    boxplot(lam1 ~ y1, xlab="Observed count at point", ylab="Expected value",
        ylim=c(0,quantile(lam1, q)), col="#0000FF80")

    MAX <- max(quantile(lam9, q), quantile(y9, q))
    plot(jitter(y9), lam9, xlim=c(0,MAX), ylim=c(0,MAX), xlab="Observed count at ABMI site",
        ylab="Expected value", col="#0000FF80", pch=19)
    abline(0,1, lty=2)
    abline(lm(lam9 ~ y9), col=2)
    abline(lm(lam9 ~ y9-1), col=2, lty=2)
    par(op)
    invisible(NULL)
}


spp <- "WTSP"
res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
V1 <- V <- try(c(list(validate(res, Xv, 0)),
    lapply(names(mods)[1:9], function(z) validate(res, Xv, z, TRUE))))
V2 <- V <- try(c(list(validate(res, Xv, 0)),
    lapply(names(mods)[1:9], function(z) validate(res, Xv, z, FALSE))))
res2 <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
V3 <- V <- try(c(list(validate(res2, Xv, 0)),
    lapply(names(mods)[1:9], function(z) validate(res2, Xv, z, TRUE))))
names(V1) <- names(V2) <- names(V3) <- c("Null", names(mods)[1:9])
plotOne("WTSP", All=list(WTSP=V1))
plotOne("WTSP", All=list(WTSP=V2))
plotOne("WTSP", All=list(WTSP=V3))

SPP <- colnames(YYv[DATv$ABMIsite != "",colSums(YYv>0)>20])

All <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    res <- load_species(file.path(ROOT, "out", PROJ, paste0(spp, ".RData")))
    V <- try(c(list(validate(res, Xv, 0)),
        lapply(names(mods)[1:9], function(z) validate(res, Xv, z))))
    if (!inherits(V, "try-error")) {
        names(V) <- c("Null", names(mods)[1:9])
        All[[spp]] <- V
    }
}




save(All, file="d:/abmi/AB_data_v2018/data/analysis/birds/validation-quick-results-2018-12-14.RData")

pdf("d:/abmi/AB_data_v2018/data/analysis/birds/validation-quick-results-2018-12-14.pdf",
    onefile=TRUE, height=10, width=10)
for (spp in colnames(YYv[DATv$ABMIsite != "",colSums(YYv>0)>20]))
    plotOne(spp, q=1, All=All)
dev.off()



load("d:/abmi/AB_data_v2018/data/analysis/birds/validation-quick-results-2018-12-12.RData")

coef9 <- function(spp, int=TRUE) {
    One <- All[[spp]]
    lam9 <- One$HF$lam9
    y9 <- One$HF$y9obs
    if (int)
        coef(lm(lam9 ~ y9)) else coef(lm(lam9 ~ y9-1))
}

cf9 <- t(sapply(names(All), coef9))
cf9 <- cf9[cf9[,2] < 2,]
cf9x <- sapply(rownames(cf9), coef9, int=FALSE)

hist(cf9[,2])
hist(cf9x)

library(lhreg)
data("lhreg_data")

Mass <- lhreg_data$mass[match(rownames(cf9), lhreg_data$spp)]
EDR <- exp(lhreg_data$logtau[match(rownames(cf9), lhreg_data$spp)])

plot(cf9x ~ log(Mass))
abline(lm(cf9x ~ log(Mass)))

plot(cf9[,2] ~ log(Mass))
abline(lm(cf9[,2] ~ log(Mass)))


cor(cf9x, Mass)
cor(cf9[,2], Mass)
cor(cf9x, EDR)
cor(cf9[,2], EDR)

## compare SM/RF estimates between validation and north
SPP <- colnames(YY)
stage <- "Water"
#spp <- "ALFL"
res <- matrix(0, length(SPP), 4)
dimnames(res) <- list(SPP, c("N_SM", "N_RF", "V_SM", "V_RF"))
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    res1 <- load_species(file.path(ROOT, "out", "north", paste0(spp, ".RData")))
    est1 <- get_coef(res1, X, stage=stage, na.out=FALSE)
    res2 <- load_species(file.path(ROOT, "out", "validation", paste0(spp, ".RData")))
    est2 <- get_coef(res2, X, stage=stage, na.out=FALSE)
    res[spp, ] <- c(est1[1,c("CMETHODSM","CMETHODRF")], est2[1,c("CMETHODSM","CMETHODRF")])
}
res <- res[rownames(cf9),]

plot(res[,c(1,3)], xlim=c(-5,5), ylim=c(-5,5));abline(0,1)
plot(res[,c(2,4)], xlim=c(-5,5), ylim=c(-5,5));abline(0,1)

plot(cf9[,2] ~ I(res[,1]-res[,3]))
abline(lm(cf9[,2] ~ I(res[,1]-res[,3])))

plot(cf9[,2] ~ I(res[,4]-res[,2]), xlab="RF est diff: Validation-North",ylab="Validation slope (ABMI 9 pt)")
abline(lm(cf9[,2] ~ I(res[,4]-res[,2])))

plot(cf9[,2] ~ res[,2])
abline(lm(cf9[,2] ~ res[,2]))


