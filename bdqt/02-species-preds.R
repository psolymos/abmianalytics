## read in species
## O - N+S averaged
## N & S: no averaging just veg or soil
## birds: p-transformation
## store results in matrix that is saved

## preamble
library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()
ST <- get_species_table()
taxa <- levels(ST$taxon)
CH <- c("S1", "O1", "O2", "O3", paste0("N", 1:11))

#tax <- "birds"
for (tax in taxa) {
    cat("Starting", tax, "============================\n")

    ## which chunk and taxon
    #ch <- "S1"
    for (ch in CH) {

        cat("Loading", ch, "--------------------------\n")

        ## load chunk
        type <- substr(ch, 1, 1)
        load(paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-", ch, "_2019-12-10.RData"))
        n <- nrow(xi)

        ## define species list
        SPP <- switch(type,
            "S"=rownames(ST)[ST$taxon == tax & ST$model_south],
            "O"=rownames(ST)[ST$taxon == tax],
            "N"=rownames(ST)[ST$taxon == tax & ST$model_north])

        ## output object
        OUT <- matrix(0, n, length(SPP))
        rownames(OUT) <- xi$UID
        colnames(OUT) <- SPP

        #spp <- "Ovenbird"
        for (spp in SPP) {
            cat(tax, "-", spp, which(SPP == spp), "/", length(SPP),  "in", ch)
            flush.console()
            gc()
            t0 <- proc.time()

            ## species model: North, South, Combo
            model <- "C"
            if (ST[spp, "model_north"] && !ST[spp, "model_south"])
                model <- "N"
            if (!ST[spp, "model_north"] && ST[spp, "model_south"])
                model <- "S"
            ## load species results
            object <- load_spclim_data(spp)
            ## predict
            if (type == "S") {
                pr <- if (ST[spp, "model_south"])
                    predict(object, xyi, soil=xi$SOILHFclass)$soil else rep(0, n)
            }
            if (type == "N") {
                pr <- if (ST[spp, "model_north"])
                    predict(object, xyi, veg=xi$VEGHFAGEclass)$veg else rep(0, n)
            }
            if (type == "O") {
                #pr <- predict(object, xyi, veg=xi$VEGHFAGEclass, soil=xi$SOILHFclass)$comb # 31sec
                prs <- if (ST[spp, "model_south"])
                    predict(object, xyi, soil=xi$SOILHFclass)$soil else rep(0, n)
                prv <- if (ST[spp, "model_north"])
                    predict(object, xyi, veg=xi$VEGHFAGEclass)$veg else rep(0, n)
                pr <- switch(model,
                    "C" = xi$wNorth * prv + (1 - xi$wNorth) * prs, # 28sec
                    "N" = prv,
                    "S" = prs)
            }
            pr[is.na(pr)] <- 0
            summary(pr)
            ## transform lambda to p for birds
            if (tax == "birds")
                pr <- p_bird(pr, area="ha", pair_adj=2)

            ## store results
            OUT[,spp] <- pr
            cat("\t", proc.time()[3] - t0[3], "sec\n")
        }
        save(OUT, file=paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, ".RData"))
    }
}




## fixing overlap zone

## preamble
library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()
ST <- get_species_table()
taxa <- levels(ST$taxon)
CH <- c("O1", "O2", "O3")

for (ch in CH) {
    cat("Loading", ch, "--------------------------\n")
    type <- "O"
    load(paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-", ch, "_2020-09-16.RData"))
    ss <- is.na(xi$SOILHFclass)
    xi <- xi[ss,]
    xyi <- xyi[ss]
    n <- nrow(xi)

    for (tax in taxa) {
        cat("Starting", tax, "============================\n")


        ## define species list
        SPP <- rownames(ST)[ST$taxon == tax & ST$model_south & ST$model_north]


        OUT <- matrix(0, n, length(SPP))
        rownames(OUT) <- xi$UID
        colnames(OUT) <- SPP

        #spp <- "Ovenbird"
        for (spp in SPP) {
            cat(tax, "-", spp, which(SPP == spp), "/", length(SPP),  "in", ch)
            flush.console()
            gc()
            t0 <- proc.time()

            object <- load_spclim_data(spp)
            pr <- if (ST[spp, "model_north"])
                predict(object, xyi, veg=xi$VEGHFAGEclass)$veg else rep(0, n)

            pr[is.na(pr)] <- 0
            summary(pr)
            ## transform lambda to p for birds
            if (tax == "birds")
                pr <- p_bird(pr, area="ha", pair_adj=2)

            ## store results
            OUT[,spp] <- pr
            cat("\t", proc.time()[3] - t0[3], "sec\n")
        }
        save(OUT, file=paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, "-weight-fix.RData"))
    }
}


library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()
ST <- get_species_table()
taxa <- levels(ST$taxon)
CH <- c("O1", "O2", "O3")

for (ch in CH) {
    for (tax in taxa) {
        cat(tax, ch, "\n")
        flush.console()
        load(paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, "-weight-fix.RData"))
        OUT1 <- OUT
        load(paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, ".RData"))
        OUT[rownames(OUT1),colnames(OUT1)] <- OUT1
        #save(OUT, file=paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, "-weight-fix-back.RData"))
        save(OUT, file=paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/", ch, ".RData"))
    }
}



## fixing S1 chunk (missing soil info):
## need to remove these rows and create N12 instead

library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()
ST <- get_species_table()
taxa <- levels(ST$taxon)

load(paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-S1_2019-12-10.RData"))
ss <- is.na(xi$SOILHFclass)
xi <- xi[ss,]
xyi <- xyi[ss]
n <- nrow(xi)

for (tax in taxa) {
    cat("Starting", tax, "============================\n")
    ## define species list
    SPP <- rownames(ST)[ST$taxon == tax & ST$model_north]

    OUT <- matrix(0, n, length(SPP))
    rownames(OUT) <- xi$UID
    colnames(OUT) <- SPP

    #spp <- "AlderFlycatcher"
    for (spp in SPP) {
        cat(tax, "-", spp, which(SPP == spp), "/", length(SPP))
        flush.console()
        gc()
        t0 <- proc.time()
        object <- load_spclim_data(spp)
        pr <- predict(object, xyi, veg=xi$VEGHFAGEclass)$veg

        pr[is.na(pr)] <- 0
        ## transform lambda to p for birds
        if (tax == "birds")
            pr <- p_bird(pr, area="ha", pair_adj=2)

        ## store results
        OUT[,spp] <- pr
        cat("\t", proc.time()[3] - t0[3], "sec\n")
    }
    save(OUT, file=paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/N12-weight-fix.RData"))
}


## remove soil NA polys from S1 chunks
library(cure4insect)
library(mefa4)
set_options(path = "d:/abmi/reports")
load_common_data()
ST <- get_species_table()
taxa <- levels(ST$taxon)

load(paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-S1_2019-12-10.RData"))
ss <- !is.na(xi$SOILHFclass)

for (tax in taxa[-1]) {
    cat(tax, "\n")
    flush.console()
    f <- paste0("s:/AB_data_v2019/bdqt/inter/", tax, "/S1.RData")
    #table(ss)
    load(f)
    #dim(OUT)
    OUT <- OUT[ss,]
    #dim(OUT)
    save(OUT, file=f)
}

load(paste0("s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-S1_2019-12-10.RData"))
ss <- !is.na(xi$SOILHFclass)

xi0 <- xi
xyi0 <- xyi

xi <- xi0[ss,]
xyi <- xyi0[ss]
save(xi, xyi, file="s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-S1_2019-12-10.RData")

xi <- xi0[!ss,]
xyi <- xyi0[!ss]
save(xi, xyi, file="s:/AB_data_v2019/bdqt/chunks/bdqt-poly-hab-N12_2019-12-10.RData")

