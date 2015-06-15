## wrsi: weighted relative suitability index

wrsi <-
    function(Y, X)
    {
        Y <- as.integer(ifelse(Y > 0, 1L, 0L))
        X <- data.matrix(X)
        n <- length(Y)
        if (nrow(X) != n)
            stop("dimension mismatch: X and Y")
        K <- ncol(X)
        if (is.null(colnames(X)))
            colnames(X) <- paste0("V", seq_len(K))
        ## relative suitability index
        ## # of available units of type k / total # of available units (any type)
        Pavail <- colSums(X) / sum(X)
        ## # of used units of type k / total # of used units (any type)
        Xu <- X * Y
        ## sum(Xu) = sum(Y) except when rowsum != 1
        Pused <- colSums(Xu) / sum(Xu)
        ## crude weighted p-occ
        Pw <- colSums(Xu) / colSums(X)
        ## Weighted Relative Suitability Index
        WRSI <- Pused / Pavail
        #Var <- (1/colSums(Xu)) - (1/sum(Xu)) + (1/colSums(X)) - (1/sum(X))
        res <- data.frame(
            WRSI=WRSI,
            zWRSI=log(WRSI),
            rWRSI=(exp(2 * log(WRSI)) - 1)/(1 + exp(2 * log(WRSI))),
            Pused=Pused,
            Pavail=Pavail,
            Pw=Pw,
            u=colSums(Xu),
            a=colSums(X))
        rownames(res) <- colnames(X)
        class(res) <- c("wrsi", "data.frame")
        res
    }

if (FALSE) {
library(mefa4)
library(RColorBrewer)

source("~/repos/abmianalytics/R/wrsi_functions.R")

load("~/Dropbox/Public/OUT_mites_2015-05-22.Rdata")
load("~/Dropbox/Public/veg-hf-clim-reg_abmi-onoff.Rdata")
tveg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")


yyy <- as.matrix(xtab(m2))
hhh <- as.matrix(dd1ha$veg_current)
iii <- intersect(rownames(yyy), rownames(hhh))
yyy <- yyy[iii,]
hhh <- hhh[iii,]
vvv <- as.character(tveg$Broad)
vvv[is.na(vvv)] <- as.character(tveg$UseInAnalysis[is.na(vvv)])
vvv[!is.na(tveg$Broad) & tveg$Broad=="Forest"] <-
    as.character(tveg$Type)[!is.na(tveg$Broad) & tveg$Broad=="Forest"]

vvv[colnames(dd1ha$veg_current) %in% c("Wetland-GrassHerb","Wetland-Shrub")] <- "WetOpen"
vvv[vvv=="Mine"] <- "UrbInd"
vvv[vvv=="Wetland"] <- "WetTreed"
vvv[vvv %in% c("Water","HWater","NonVeg","Wetland-Bare")] <- "EXCLUDE"

hhh <- groupSums(hhh, 2, vvv)
hhh <- hhh[,colnames(hhh) != "EXCLUDE"]
hhh <- hhh / ifelse(rowSums(hhh) == 0, 1, rowSums(hhh))

hhh <- hhh[, c("Decid", "Mixwood", "Conif", "Pine", "Shrub", "GrassHerb", 
    "WetTreed", "WetOpen", "Swamp", "Cult", "UrbInd", "HardLin", "SoftLin")]

keep <- rowSums(hhh) > 0
yyy <- yyy[keep,]
hhh <- hhh[keep,]
yyy <- yyy[,colSums(yyy>0) > 1]

plot_wrsi <- 
function(Y, X, ...)
{
    w <- wrsi(Y, X)
    rw <- w$rWRSI
    names(rw) <- rownames(w)

    col <- brewer.pal(8, "Accent")[c(1,1,1,1, 3,3, 2,2,2, 5,5, 8,8)]
    op <- par(mar=c(6,4,2,2)+0.1, las=2)
    tmp <- barplot(rw, horiz=FALSE, ylab="Affinity", 
        space=NULL, col=col, border=col, #width=sqrt(w$Pavail),
        ylim=c(-1,1), axes=FALSE, ...)
    axis(2)
    abline(h=0, col="red4", lwd=2)
    par(op)
    invisible(w)
}

spp <- 1
plot_wrsi(yyy[,spp], hhh)

pdf("~/Dropbox/Public/mites-pw.pdf", onefile=TRUE)
for (i in 1:ncol(yyy)) {
    main <- paste0(as.character(taxa(m2)[colnames(yyy)[i], "SPECIES_OLD"]),
        " (", sum(yyy[,i] > 0), " detection", ifelse(sum(yyy[,i] > 0) > 2,
        "s)", ")"))
    plot_wrsi(yyy[,i], hhh, main=main)
}
dev.off()

}
