do_1spec1run_noW <- function(j, i, mods, 
silent=FALSE, hsh_name=NA, CAICalpha=1, method=c("oc","lt"))
{
    method <- match.arg(method)
    select_hsh <- !is.na(hsh_name)
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    if (select_hsh)
        hsh <- HSH[BB[,j],]
    ## empty objects for storing results
    nmods <- length(mods)
    nnmods <- sapply(mods, length)
    mid <- numeric(nmods)
    bestList <- vector("list", nmods)
    caicList <- vector("list", nmods)
    ## Null
    null <- glm_skeleton(glm(y ~ 1, 
        x, 
        family=poisson(), 
        offset=off, 
        #weights=w,
        x=FALSE, y=FALSE, model=FALSE), CAICalpha=CAICalpha)
    best <- null
    if (select_hsh) {
        HABV <- x[,hsh_name]
        ## opticut based approach for core habitat delineation
        if (method=="oc") {
            require(opticut)
            oc <- opticut(y ~ ROAD01, data=x, strata=HABV, dist="poisson",
                offset=off, comb="rank")
            part <- drop(bestpart(oc))
            habmod <- glm_skeleton(bestmodel(oc)[[1]], CAICalpha=CAICalpha)
            cv <- NULL
            lam <- NULL
            ocres <- drop(as.matrix(summary(oc)$summary[,c("I","beta0","beta1","logLR","w")]))
            Prob <- table(HABV, part)[,"1"]
            ## missing/dropped levels are NaN=0/0
            Prob[is.na(Prob)] <- 0
            Hi <- names(Prob)[Prob > 0]
        }
        ## Lorenz-tangent approach for core habitat delineation
        if (method=="lt") {
            habmod <- glm_skeleton(try(glm(y ~ HABV + ROAD01,
                x,
                family=poisson(), 
                offset=off, 
                #weights=w,
                x=FALSE, y=FALSE, model=FALSE), silent=silent), CAICalpha=CAICalpha)
            ## need to correct for linear effects
            ## so that we estimate potential pop in habitats (and not realized)
            XHSH <- model.matrix(~ HABV + ROAD, x)
            XHSH[,"ROAD"] <- 0 # not predicting edge effects
            ## some levels might be dropped (e.g. Marsh)
            XHSH <- XHSH[,names(habmod$coef)]
            lam <- exp(drop(XHSH %*% habmod$coef))
            cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
            Freq <- table(hab=HABV, lc=ifelse(lam >= cv$lam, 1, 0))
            Prob <- Freq[,"1"] / rowSums(Freq)
            ## missing/dropped levels are NaN=0/0
            Prob[is.na(Prob)] <- 0
            Hi <- names(Prob)[Prob > 0.5]
            ocres <- NULL
        }
        x$HSH <- unname(rowSums(hsh[, Hi, drop=FALSE]))
        x$HSH2 <- x$HSH^2
    } else {
        Hi <- NULL
        lam <- NULL
        cv <- NULL
        habmod <- NULL
        ocres <- NULL
    }
    ## looping through models list
    for (l1 in 1:nmods) {
        if (nnmods[l1] > 0) {
            mlist <- vector("list", nnmods[l1])
            glist <- vector("list", nnmods[l1])
            for (l2 in 1:nnmods[l1]) {
                mlist[[l2]] <- glm_skeleton(try(update(object=best, 
                    formula=mods[[l1]][[l2]]), silent=silent), CAICalpha=CAICalpha)
            }
            mcaic <- sapply(mlist, "[[", "caic")
            attr(mcaic, "StartCAIC") <- best$caic
            for (l2 in 1:length(mlist)) { # check convergence
                if (mlist[[l2]]$class != "try-error" && !mlist[[l2]]$converge)
                    mcaic[l2] <- 2*.Machine$double.xmax^(1/3)
            }
            dcaic <- mcaic - best$caic
            mmid <- which.min(dcaic)
            if (dcaic[mmid] < 0) {
                best <- mlist[[mmid]]
                mid[l1] <- mmid
                gofbest <- glist[[mmid]]
            }
            caicList[[l1]] <- mcaic
        }
        bestList[[l1]] <- best
    }
    ## final assembly
    out <- list(species=i, iteration=j,
        null=null$coef,
        null_caic=null$caic,
        caic=caicList,
        coef=lapply(bestList, "[[", "coef"),
        mid=mid,
        hi=Hi,
        lc=cv,
        alpha=CAICalpha,
        method=method,
        ocres=ocres,
        #nmax=nmax,
        #w_id=w_id,
        habmod=habmod$coef,
        hsh_name=hsh_name)
    out
}

## z: colname for strata
cut_1spec1run_noW <- function(j, i, z)
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    HABV <- x[,z]
    ## opticut based approach for core habitat delineation
    require(opticut)
    oc <- opticut(y ~ ROAD01, data=x, strata=HABV, dist="poisson",
        offset=off, comb="rank")
    part <- drop(bestpart(oc))
    habmod_oc <- glm_skeleton(bestmodel(oc)[[1]])
    ocres <- drop(as.matrix(summary(oc)$summary[,c("I","beta0","beta1","logLR","w")]))
    Prob <- table(HABV, part)[,"1"]
    ## missing/dropped levels are NaN=0/0
    Prob[is.na(Prob)] <- 0
    Hi_oc <- names(Prob)[Prob > 0]
    ## Lorenz-tangent approach for core habitat delineation
    habmod <- glm_skeleton(try(glm(y ~ HABV + ROAD01,
        x,
        family=poisson(), 
        offset=off, 
        #weights=w,
        x=FALSE, y=FALSE, model=FALSE)), CAICalpha=CAICalpha)
    ## need to correct for linear effects
    ## so that we estimate potential pop in habitats (and not realized)
    XHSH <- model.matrix(~ HABV + ROAD01, x)
    XHSH[,"ROAD01"] <- 0 # not predicting edge effects
    ## some levels might be dropped (e.g. Marsh)
    XHSH <- XHSH[,names(habmod$coef)]
    lam <- exp(drop(XHSH %*% habmod$coef))
    cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
    Freq <- table(hab=HABV, lc=ifelse(lam >= cv$lam, 1, 0))
    Prob <- Freq[,"1"] / rowSums(Freq)
    ## missing/dropped levels are NaN=0/0
    Prob[is.na(Prob)] <- 0
    Hi_lc <- names(Prob)[Prob > 0.5]
    ## final assembly
    out <- list(species=i, iteration=j,
        hi_oc=Hi_oc,
        hi_lc=Hi_lc,
        lc=cv,
        ocres=ocres,
        #nmax=nmax,
        #w_id=w_id,
        habmod_oc=habmod_oc$coef,
        habmod_lc=habmod$coef)
    out
}
#system.time(x <- cut_1spec1run_noW(j=1, i="OVEN", z="hab1ec"))

level_1spec1run_noW <- function(j, i, z)
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    require(opticut)
    Z <- model.matrix(~ ROAD01 - 1, x)
    Z[,1] <- Z[,1] - mean(Z[,1], na.rm=TRUE)
    optilevel(y=y, x=x[,z], z=Z, dist="poisson",
        offset=off)
}
#system.time(x <- level_1spec1run_noW(j=1, i="OVEN", z="hab1ec"))

