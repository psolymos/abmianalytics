# chacking mammal predictions: exact and approximate

library(mefa4)
library(ggplot2)
load("s:/AB_data_v2020/Results/COEFS-ALL2.RData")
ROOT <- "s:/AB_data_v2020/Results/pred"

taxon <- "mammals"
SPPn <- dimnames(COEFS2[[taxon]]$north)[[1]]
SPPs <- dimnames(COEFS2[[taxon]]$south)[[1]]
SPP <- sort(unique(c(SPPn, SPPs)))

spp <- "Coyote"

res <- NULL
for (spp in SPP) {
    cat(spp, "\n")

load(file.path(ROOT, taxon, paste0(spp, ".RData"))) # Ncr, Nrf
z <- data.frame(
    LinkID=rownames(Ncr),
    Ref=rowSums(Nrf),
    Curr=rowSums(Ncr))

x <- read.csv(paste0(
    "s:/AB_data_v2020/Results/",
    "Camera mammal models revised June 2020/",
    "Km2 summaries/", spp, ".csv"))
rownames(x) <- x$LinkID

s <- intersect(rownames(x), rownames(z))

u <- data.frame(
    rf_exact=x[s,"Ref"],
    rf_approx=z[s,"Ref"],
    cr_exact=x[s,"Curr"],
    cr_approx=z[s,"Curr"])

res <- rbind(res, c(rf=cor(u[,1:2])[2,1], cr=cor(u[,3:4])[2,1]))
}
rownames(res) <- SPP

k <- lm(rf_approx ~ rf_exact, u)
p <- ggplot(u, aes(x=rf_exact, y=rf_approx)) +
    geom_bin2d() +
    geom_abline(slope=1, col=2) +
    geom_abline(intercept=coef(k)[1], slope=coef(k)[2], col=2, lty=2) +
    labs(title=paste(spp, "Ref")) +
    theme_minimal()

k <- lm(cr_approx ~ cr_exact, u)
p <- ggplot(u, aes(x=cr_exact, y=cr_approx)) +
    geom_bin2d() +
    geom_abline(slope=1, col=2) +
    geom_abline(intercept=coef(k)[1], slope=coef(k)[2], col=2, lty=2) +
    labs(title=paste(spp, "Curr")) +
    theme_minimal()

head(x)
head(z)
compare_sets(x$LinkID, z$LinkID)


## check predictive maps

load("d:/abmi/AB_data_v2020/data/analysis/kgrid_table_km.RData") # kgrid
## chSoil/chVeg/trSoil/trVeg
load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_ref_2018_transitions_wide.RData")
#load("d:/abmi/AB_data_v2020/data/analysis/veghf/veghf_w2w_2018_wide.RData")
trVeg <- trVeg[rownames(kgrid),rownames(chVeg)]
trSoil <- trSoil[rownames(kgrid),rownames(chSoil)]


xvr <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "rf"])
xvc <- groupSums(trVeg, 2, chVeg[colnames(trVeg), "cr"])
setdiff(colnames(xvr), colnames(xvc))
setdiff(colnames(xvc), colnames(xvr))
i <- intersect(colnames(xvc), colnames(xvr))
d <- (xvr[,i] - xvc[,i]) / 10^6
range(d)
summary(rowSums(d))
table(rowSums(d < 0))
j <- which(rowSums(d < 0) > 13)
u <- d[j,]
u <- u[,colSums(u<0) > 1]

xvc <- dd_2018$veg_current
xvr <- dd_2018$veg_reference
d <- (xvr[,i] - xvc[,i]) / 10^6
range(d)

load("s:/AB_data_v2018/data/analysis/grid/veg-hf_transitions_v61hf2016v3WildFireUpTo2016.Rdata")


## checking overlap region issues with birds

# run parts from 03-maps-and-sector.R
taxon <- "birds"
spp <- "AmericanCrow"


SPPn <- if (taxon!="birds")
    dimnames(COEFS[[taxon]]$north)[[1]] else dimnames(COEFS[[taxon]]$north$joint)[[1]]
SPPs <- if (taxon!="birds")
    dimnames(COEFS[[taxon]]$south)[[1]] else dimnames(COEFS[[taxon]]$south$joint)[[1]]
SPP <- sort(union(SPPn, SPPs))

type <- "C" # combo species (N+S)
M <- list(N=spp %in% SPPn, S=spp %in% SPPs)
if (M$N & !M$S)
    type <- "N"
if (!M$N & M$S)
    type <- "S"
cfn <- if (type == "S")
    NULL else COEFS[[taxon]]$north$joint[spp,,]
cfs <- if (type == "N")
    NULL else COEFS[[taxon]]$south$joint[spp,,]
XclimS <- Xclim_bird_S
XclimN <- Xclim_bird_N
FUN <- function (eta)
    pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)


i <- 1
if (type != "N") {
                gc()
                ## south calculations for the i'th run
                #compare_sets(rownames(cfs), chSoil$cr2)
                bscr <- cfs[chScr$cr2, i] # current land cover
                bsrf <- cfs[chSrf$rf2, i] # reference land cover
                ## space-climate coefs
                bscl <- if (taxon == "birds")
                    cfs[colnames(Xclim_bird), i] else cfs[colnames(Xclim_nonb), i]
#                if (taxon == "birds") {
#                    bscl[bscl < -10] <- -10
#                    bscl[bscl > 10] <- 10
#                }
                bscl[is.na(bscl)] <- 0 # this happens for habitat elements
                bspa <- cfs["pAspen", i]
                ## additive components for south
                if (taxon=="mammals") {
                    ## space/clim includes presence-absence piece
                    #Total.Abundance.approx<- exp(log(TA.vegHF) + SC.pOcc + pAspen.pa + pAspen.agp)
                    pApa <- COEFS[[taxon]]$pAspenPA[spp]
                    muscl <- drop(cbind(XclimS, pA) %*% c(bscl, pApa))
                    if (spp == "Pronghorn")
                        muscl <- (muscl - mean(muscl))/1000
                } else {
                    muscl <- drop(XclimS %*% bscl)
                }
                muspa <- pA * bspa

                muscr <- matrix(muscl + muspa, nrow=nrow(Pscr), ncol=ncol(Pscr))
                muscr <- t(t(muscr) + bscr)
                NScr <- as.matrix(groupSums(Pscr * FUN(muscr), 2, chScr$sector_use))
                NScr <- cbind(NScr, Forestry=0)

                musrf <- matrix(muscl + muspa, nrow=nrow(Psrf), ncol=ncol(Psrf))
                musrf <- t(t(musrf) + bsrf)
                NSrf <- as.matrix(groupSums(Psrf * FUN(musrf), 2, chSrf$sector_use))
                NSrf <- cbind(NSrf, Forestry=0)
            } else {
                NScr <- NULL
                NSrf <- NULL
            }
            if (type != "S") {
                gc()
                ## north calculations for the i'th run
                #compare_sets(rownames(cfn), chVeg$cr2)
                tmpn <- c(cfn[,i], Bare=-10^4, SnowIce= -10^4)
                bncr <- tmpn[chVcr$cr2] # current land cover
                bnrf <- tmpn[chVrf$rf2] # reference land cover
                ## space-climate coefs
                bncl <- if (taxon == "birds")
                    cfn[colnames(Xclim_bird), i] else cfn[colnames(Xclim_nonb), i]
                bncl[is.na(bncl)] <- 0 # this happens for habitat elements
                ## additive components for north
                muncl <- drop(XclimN %*% bncl)

                muncr <- matrix(muncl, nrow=nrow(Pncr), ncol=ncol(Pncr))
                muncr <- t(t(muncr) + bncr)
                NNcr <- as.matrix(groupSums(Pncr * FUN(muncr), 2, chVcr$sector_use))

                munrf <- matrix(muncl, nrow=nrow(Pnrf), ncol=ncol(Pnrf))
                munrf <- t(t(munrf) + bnrf)
                NNrf <- as.matrix(groupSums(Pnrf * FUN(munrf), 2, chVrf$sector_use))
            } else {
                NNcr <- NULL
                NNrf <- NULL
}
if (type == "C") {
                # averaging comes here
                NNcr <- NNcr[match(rownames(kgrid), rownames(NNcr)),]
                NNrf <- NNrf[match(rownames(kgrid), rownames(NNrf)),]
                NNcr[is.na(NNcr)] <- 0
                NNrf[is.na(NNrf)] <- 0
                rownames(NNcr) <- rownames(NNrf) <- rownames(kgrid)

                NScr <- NScr[match(rownames(kgrid), rownames(NScr)),]
                NSrf <- NSrf[match(rownames(kgrid), rownames(NSrf)),]
                NScr[is.na(NScr)] <- 0
                NSrf[is.na(NSrf)] <- 0
                rownames(NScr) <- rownames(NSrf) <- rownames(kgrid)

                Ncr <- kgrid$wN * NNcr + (1-kgrid$wN) * NScr
                Nrf <- kgrid$wN * NNrf + (1-kgrid$wN) * NSrf
}

library(raster)
rt <- raster(system.file("extdata/AB_1km_mask.tif", package="cure4insect"))
make_raster <- function(value, rc, rt) {
    value <- as.numeric(value)
    r <- as.matrix(Xtab(value ~ Row + Col, rc))
    r[is.na(as.matrix(rt))] <- NA
    raster(x=r, template=rt)
}
f <- function(x) {
    v <- rowSums(x)
    q <- quantile(v, 0.99)
    v[v>q] <- q
    hist(v)
    make_raster(v, kgrid, rt)
}
l <- stack(list(south=f(NScr), north=f(NNcr), combo=f(Ncr)))

par(mar=c(2,2,3,10))
plot(l, col=hcl.colors(100)[30:100])

i <- sample.int(nrow(kgrid), 10^4)
plot(rowSums(NScr)[i], rowSums(Ncr)[i], pch=19,
    col=ifelse(kgrid$NRNAME=="Grassland", "#ff000022", "#22222222")[i])
abline(0,1,col=2)
v <- rowSums(NNcr)
q <- quantile(v, 0.99)
v[v>q] <- q
plot(make_raster(v, kgrid, rt), col=hcl.colors(25))


## bootstrap based efficiencies
## can we use 10% of the pixels to predict the rest?
## compare agains the 100% based timings

library(qs)
ROOT2 <- "s:/AB_data_v2020/Results/pred-boot"

taxon <- "birds"
SPP <- list.dirs(file.path(ROOT2, taxon), full.names=FALSE)
SPP <- SPP[-1]
B <- 100
ii <- 1:B
names(ii) <- 1:B
ii <- ii[order(names(ii))]

spp <- SPP[1]
fl <- list.files(file.path(ROOT2, taxon, spp), full.names = TRUE)
fl <- fl[order(ii)]
ft <- file.info(fl)$ctime
dt <- as.numeric(diff(ft))

i <- 1
qreadm(file.path(ROOT2, taxon, spp, paste0(spp, "-", i, ".qrda")))
M <- matrix(0, nrow(Ncr), B)
rownames(M) <- rownames(Ncr)
M[,1] <- rowSums(Ncr)

for (i in 2:B) {
    cat(i, "\n")
    flush.console()
    qreadm(file.path(ROOT2, taxon, spp, paste0(spp, "-", i, ".qrda")))
    M[,i] <- rowSums(Ncr)
}

p <- 0.1
Rnd <- sample.int(100, nrow(M), replace=TRUE)
Mt <- data.frame(M[Rnd <= round(p*100),])
Mv <- data.frame(M[Rnd > round(p*100),])

library(e1071)

system.time({
    m <- svm(X1 ~ ., Mt)
    pr <- predict(m, Mv)
})

lm(Mv$X1 ~ pr)
cor(cbind(Mv$X1, pr))[1,2]
# taking 2-3x longer than calculation

## see if SD's are unaffected

res <- NULL
for (j in 1:10) {
    MM <- M[,sample(100, 100)]
    SD <- matrix(0, nrow(Ncr), 10)
    for (i in 1:10) {
        cat(j, i, "\n")
        flush.console()
        SD[,i] <- apply(MM[,1:(i*10)], 1, sd)
    }

    L <- list()
    for (i in 1:9) {
        L[[i]] <- lm(SD[,i] ~ SD[,10])
    }

    U <- t(sapply(L, function(z) {
        c(coef(z), sigma(z))
    }))
    colnames(U) <- c("b0", "b1", "s")
    U

    plot(1:9*10, U[,1], type="b")
    plot(1:9*10, U[,2], type="b")
    plot(1:9*10, U[,3], type="b")


    B <- 1:9*10
    u <- data.frame(coef=rep(colnames(U), each=9), estimate=as.numeric(U), B=B)
    u <- rbind(u, data.frame(
        coef=colnames(U),
        estimate=c(0,1,0),
        B=c(100, 100, 100)
    ))
    u$run <- j
    res <- rbind(res, u)
}


library(ggplot2)
ggplot(res, aes(x=B, y=estimate, color=as.factor(run))) +
    geom_line() +
#    geom_smooth() +
    facet_wrap(vars(coef), scales="free_y")
