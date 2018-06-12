library(mefa4)
library(opticut)

ROOT <- "e:/peter/AB_data_v2016"

up <- function() {
    source("~/repos/bragging/R/glm_skeleton.R")
    source("~/repos/abmianalytics/R/results_functions.R")
    source("~/repos/bamanalytics/R/makingsense_functions.R")
#    source("~/repos/abmianalytics/R/wrsi_functions.R")
#    source("~/repos/abmianalytics/R/results_functions1.R")
#    source("~/repos/abmianalytics/R/results_functions2.R")
    invisible(NULL)
}
up()

load(file.path(ROOT, "out", "birds", "data", "data-josmshf.Rdata"))
#Terms <- getTerms(mods, "list")
Xn <- model.matrix(getTerms(mods, "formula"), DAT)
colnames(Xn) <- fixNames(colnames(Xn))
fl <- list.files(file.path(ROOT, "out", "birds", "results", "josmshf"))
SPP <- gsub(".Rdata", "", gsub("birds_abmi-josmshf_", "", fl))

#spp <- "CAWA"
j <- 1
for (spp in SPP) {
    cat(spp);flush.console()

    fn <- file.path(ROOT, "out", "birds", "results", "josmshf",
        paste0("birds_abmi-josmshf_", spp, ".Rdata"))
    res <- loadSPP(fn)

    mid <- res[[j]]$mid
    names(mid) <- names(mods)
    mid <- mid[names(mid) != "Year"]
    f1 <- y ~ 1
    for (k in names(mid)) {
        if (mid[k] > 0)
            f1 <- update(f1, mods[[k]][[mid[k]]])
    }
    f2 <- update(f1, . ~ . + HSH05_KM)

    estHf <- suppressWarnings(getEst(res, stage=which(names(mods)=="HF"), na.out=FALSE, Xn))[j,]
    estSp <- suppressWarnings(getEst(res, stage=which(names(mods)=="Space"), na.out=FALSE, Xn))
    Xn2 <- Xn[BB[,j],]
    Xn2[,grep("ROAD01",colnames(Xn2))] <- 0

    prSp <- exp(drop(Xn2 %*% estSp[j,]))

    g <- sum_by(prSp, DAT[BB[,j],"hab1ec"])
    l <- lorenz(g[,"x"]/g[,"by"], g[,"by"])
    s <- summary(l)
    lab <- rownames(l[l[,"x"] >= s["x[t]"],])
    DAT$HSH_KM <- rowSums(HSH[,lab,drop=FALSE])
    #DAT$HSH2_KM <- DAT$HSH_KM^2
    DAT$HSH05_KM <- sqrt(DAT$HSH_KM)

    y <- as.numeric(YY[BB[,j], spp])
    M <- glm(f2,
        DAT[BB[,j],],
        family=poisson(),
        offset=OFF[BB[,j], spp],
        x=FALSE, y=FALSE, model=FALSE)
    mod <- glm_skeleton(M, CAICalpha=attr(res,"CAICalpha"))
    mod$hsh_coef <- mod$coef
    names(mod$hsh_coef) <- fixNames(names(mod$hsh_coef))
    mod$no_hsh_coef <- estHf[estHf != 0]
    mod$lc_input <- g
    mod$lc <- l
    mod$hsh_labels <- lab
    AIC0 <- if (mid[which(names(mods)=="HF")] == 0) {
        attr(res[[j]]$caic[[which(names(mods)=="HF")]], "StartCAIC")
    } else {
        res[[j]]$caic[[which(names(mods)=="HF")]][mid[which(names(mods)=="HF")]]
    }
    mod$delta <- mod$caic - AIC0
    cat("\t", mod$delta, "\n")
    cat(lab, "\n")

    save(mod, file=paste0("e:/peter/josm/2018/hsh-estimates/", spp, ".Rdata"))
}

mod$no_hsh_coef[grep("_KM", names(mod$no_hsh_coef))]
mod$coef[grep("_KM", names(mod$coef))]

