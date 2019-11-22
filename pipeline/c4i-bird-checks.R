library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()

## load soil & veg

load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")
load("d:/abmi/AB_data_v2018/data/analysis/grid/veg-hf_grid_v61hf2016v3WildFireUpTo2016.Rdata")

ts <- read.csv("~/repos/abmianalytics/lookup/lookup-soil-hf-v61.csv")
rownames(ts) <- ts[,1]

sval <- as.character(ts[colnames(dd_kgrid$soil_current), "UseInAnalysisCoarse"])
sval[sval %in% c("SoilWater", "HWater")] <- "Water"
sval[sval %in% c("SoilUnknown", "HFor")] <- "XXX"
ps <- as.matrix(groupSums(dd_kgrid$soil_current, 2, sval))
ps <- ps / ifelse(rowSums(ps)==0, 1, rowSums(ps))
summary(rowSums(ps))
ps <- ps[,colnames(ps) != "XXX"]
summary(rowSums(ps))
compare_sets(get_levels()$soil, colnames(ps))
setdiff(get_levels()$soil, colnames(ps))
setdiff(colnames(ps), get_levels()$soil)

tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age-v61.csv")
rownames(tv) <- tv[,1]

pv <- as.matrix(groupSums(dd_kgrid$veg_current, 2, tv[colnames(dd_kgrid$veg_current), "CoefTabs"]))
pv <- pv / ifelse(rowSums(pv)==0, 1, rowSums(pv))
summary(rowSums(pv))
compare_sets(get_levels()$veg, colnames(pv))

SPP <- get_all_species("birds")
XY <- get_id_locations()
all(rownames(get_id_table()) == rownames(pv))
ss <- sample(seq_len(nrow(kgrid)), 10^4)

## load species stuff

res <- list()

for (spp in SPP) {

    cat(spp, "\n")
    flush.console()

    y <- load_species_data(spp)
    r <- rasterize_results(y)


    v1 <- extract(r[["NC"]], XY)

    ## predict species stuff

    object <- load_spclim_data(spp)
    pr <- predict_mat(object, XY, pv, ps)
    #if (spp %in% c("PineSiskin", "SpottedSandpiper"))
    #    pr$soil <- NULL

    if (!is.null(pr$veg))
        pr$veg[is.na(pr$veg)] <- 0
    if (!is.null(pr$soil))
        pr$soil[is.na(pr$soil)] <- 0

    if (is.null(pr$veg)) {
        prv <- rep(0, nrow(pr$soil))
    } else {
        prv <- rowSums(pr$veg)
    }
    if (is.null(pr$soil)) {
        prs <- rep(0, nrow(pr$veg))
    } else {
        prs <- rowSums(pr$soil)
    }
    v2 <- combine_veg_soil(XY, prv, prs)
    #v2 <- ifelse(is.na(v2), 0, v2)
    q <- quantile(v2, 0.99, na.rm=TRUE)
    v3 <- v2
    v3[!is.na(v3) & v3 > q] <- q

    r2 <- .make_raster(v2, kgrid, .read_raster_template())
    r3 <- .make_raster(v3, kgrid, .read_raster_template())

    #summary(v1)
    #summary(v2)
    #summary(v3)
    ct2 <- cor.test(v1, v2, method = "spearman")
    ct3 <- cor.test(v1, v3, method = "spearman")
    mx1 <- max(v1, na.rm=TRUE)
    mx2 <- max(v2, na.rm=TRUE)
    mx3 <- max(v3, na.rm=TRUE)

    #plot(v1[ss], v2[ss]);abline(0,1,col=2)


    crr <- cut(ct3$est, c(-Inf, 0.25, 0.5, 0.75, 0.9, Inf))
    flag <- c("0000check-", "000check-", "00check-", "0check-",
        "")[crr]
    png(paste0("d:/tmp/birds-check/", flag, spp, ".png"), height=600, width=400*3, pointsize=20)
    op <- par(mfrow=c(1,3), mar=c(6, 1, 4, 4))
    plot(r, "NC", sub="stored", box=FALSE, axes=FALSE, main=spp)
    plot(r3, sub="pred + trunc", box=FALSE, axes=FALSE,
        main=paste0("rho=",round(ct3$est,3), "\nmax=",
            round(mx3/mx1, 3)))
    plot(r2, sub="predicted", box=FALSE, axes=FALSE,
        main=paste0("rho=",round(ct2$est,3), "\nmax=",
            round(mx2/mx1, 3)))
    par(op)
    dev.off()


    res[[spp]] <- c(
        max_stored=max(v1, na.rm=TRUE), max_prtr=max(v3, na.rm=TRUE),
        mean_stored=mean(v1, na.rm=TRUE), mean_prtr=mean(v3, na.rm=TRUE),
        rho_prtr=ct3$estimate)
}

## compare bootstrap mean and individual coefs for 100 runs

library(mefa4)
library(cure4insect)
load("d:/abmi/AB_data_v2018/data/analysis/kgrid_table_km.Rdata")

e <- new.env()
load("d:/abmi/AB_data_v2018/data/analysis/birds/ab-birds-all-2018-11-29.RData", envir=e)
tax <- e$tax
rm(e)

tf <- function(x, p=0.99) {
    q <- quantile(x, p, na.rm=TRUE)
    x[!is.na(x) & x > q] <- q
    x
}

#spp <- "Ovenbird"
spp <- "PineSiskin" # south model is wrong
#spp <- "BlackburnianWarbler" # something is funny ~ scaling + roads?
#spp <- "SpottedSandpiper" # south model is wrong
#spp <- "GrayCatbird"
#spp <- "CommonYellowthroat"


aou <- rownames(tax)[tax$sppid == spp]

eb <- new.env()
load(paste0("d:/abmi/AB_data_v2018/data/analysis/birds/pred/2019-09-20/", aou, ".RData"), envir=eb)

b <- eb$CURRB
b <- b[,colSums(is.na(b))==0]
for (i in 1:ncol(b))
    b[,i] <- tf(b[,i])
#mn <- tf(rowMeans(b))
mn <- apply(b, 1, median)
mnb <- tf(rowSums(eb$CURR))
mn1 <- tf(b[,1])

rn <- intersect(names(mn), names(mnb))
rn <- intersect(rn, names(mn1))


cor(cbind(mn[rn], mnb[rn], mn1[rn]))
corr <- cor(cbind(mn[rn], b[rn,]))[-1,1]
corr


mr <- .make_raster(mn, kgrid[names(mn),], .read_raster_template())
mrb <- .make_raster(mnb, kgrid[names(mnb),], .read_raster_template())
mr1 <- .make_raster(mn1, kgrid[names(mn1),], .read_raster_template())

op <- par(mfrow=c(1,3), mar=c(6, 1, 4, 4))
plot(mr, main=spp, box=FALSE, axes=FALSE, sub="mean now")
plot(mrb, box=FALSE, axes=FALSE, sub="bavg")
plot(mr1, box=FALSE, axes=FALSE, sub="1st")
par(op)




predict_mat.c4ispclim <-
function(object, xy, veg, soil, method="simple", ...)
{
    xy <- .tr_xy(xy)
    ## coefs in object are on log/logit scale, need linkinv
    fi <- if (object$taxon == "birds")
        poisson("log")$linkinv else binomial("logit")$linkinv
    if (missing(veg) && missing(soil))
        stop("veg or soil must be provided")
    xy <- spTransform(xy, proj4string(.read_raster_template()))
    if (!missing(veg) && !is.null(veg)) {
        if (is.null(object$cveg)) {
            warning(sprintf("veg based estimates are unavailable for %s", object$species))
            Nveg <- NULL
        } else {
            if (nrow(veg) != nrow(coordinates(xy)))
                stop("nrow(veg) must equal number of points in xy")
            .check(as.factor(colnames(veg)), names(object$cveg))
            iveg <- extract(object$rveg, xy, method)
            imatv <- t(array(iveg, dim(veg), dimnames(veg)))
            mveg <- object$cveg[match(colnames(veg), names(object$cveg))]
            Nveg <- fi(t(mveg + imatv)) * veg
            if (any(colnames(veg) == "SoftLin") && object$taxon == "birds" && object$version == "2017") {
                warning("veg contained SoftLin: check your assumptions")
                Nveg[,colnames(veg) == "SoftLin"] <- NA
            }
        }
    } else {
        Nveg <- NULL
    }
    if (!missing(soil) && !is.null(soil)) {
        if (is.null(object$csoil)) {
            warning(sprintf("soil based estimates are unavailable for %s", object$species))
            Nsoil <- NULL
        } else {
            if (nrow(soil) != nrow(coordinates(xy)))
                stop("nrow(veg) must equal number of points in xy")
            .check(as.factor(colnames(soil)), names(object$csoil))
            isoil <- extract(object$rsoil, xy, method)
            ## pAspen here used as habitat covariate NOT as North weight
            rpa <- raster(system.file("extdata/pAspen.tif", package="cure4insect"))
            ipa <- extract(rpa, xy, method)
            imats <- t(array(object$caspen * ipa + isoil, dim(soil), dimnames(soil)))
            msoil <- object$csoil[match(colnames(soil), names(object$csoil))]
            Nsoil <- fi(t(msoil + imats)) * soil
            if (any(colnames(soil) == "SoftLin") && object$taxon == "birds" && object$version == "2017") {
                warning("soil contained SoftLin: check your assumptions")
                Nsoil[,colnames(soil) == "SoftLin"] <- NA
            }
        }
    } else {
        Nsoil <- NULL
    }
    OUT <- list(veg=Nveg, soil=Nsoil)
    class(OUT) <- c("c4ippredmat")
    OUT
}
