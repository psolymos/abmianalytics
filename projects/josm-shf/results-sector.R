library(mefa4)
#library(cure4insect)

ROOT <- "e:/peter/AB_data_v2016"
#INDIR <- paste0("e:/peter/josm/2018/earlyseralTRUE")
#INDIR <- paste0("e:/peter/josm/2018/earlyseralFALSE")
INDIR <- paste0("e:/peter/josm/2018/hshfix")

sectors_all <- c("Agriculture", "EnergyLin", "EnergyMW", "Forestry", "Misc",
    "RuralUrban","Transportation")

load(file.path(ROOT, "out", "kgrid", "kgrid_table.Rdata")) # kgrid
load(file.path(ROOT, "out", "kgrid", "veg-hf_1kmgrid_fix-fire_fix-age0.Rdata")) # dd1km_pred
tv <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
tv$SectorForSeMs <- factor(ifelse(is.na(tv$SectorForSeMs), "NATIVE", as.character(tv$SectorForSeMs)),
    c("NATIVE", "Agriculture", "EnergyLin", "EnergyMW",
    "Forestry", "Misc", "RuralUrban","Transportation"))
stopifnot(all(rownames(kgrid) == rownames(dd1km_pred$veg_current)))
stopifnot(all(colnames(dd1km_pred$veg_current) == rownames(tv)))
HFarea <- groupSums(dd1km_pred$veg_current, 2, tv$SectorForSeMs)
rm(dd1km_pred)

## subset definition
kgrid$subset <- kgrid$BCRCODE == "  6-BOREAL_TAIGA_PLAINS"
#kgrid$subset <- kgrid$POINT_Y > 50 & kgrid$NRNAME != "Grassland"

#rt <- .read_raster_template()
#rsub <- .make_raster(as.integer(kgrid$subset), kgrid, rt)
#plor(rsub)

Ahf <- colSums(HFarea[kgrid$subset,])
names(Ahf)[1] <- "Native"


fln <- list.files(file.path(ROOT, "out", "birds", "results", "josmshf"))
fln <- sub("birds_abmi-josmshf_", "", fln)
fln <- sub(".Rdata", "", fln)

SPP <- fln
#SPP <- c("BOCH","ALFL","BTNW","CAWA","OVEN","OSFL","RWBL")

spp <- "CAWA"
sector_res <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    Nsect <- list()
    for (sect in c("All", sectors_all)) {
        e <- new.env()
        load(file.path(INDIR, sect, paste0(spp, ".RData")), envir=e)
        #stopifnot(all(rownames(kgrid) == rownames(e$SA.Curr)))
        Nsect[[sect]] <- cbind(cr=colSums(e$SA.Curr[kgrid$subset,]),
            rf=colSums(e$SA.Ref[kgrid$subset,]))
    }

    ## reference abundance should not differ
    #round(100*sapply(Nsect, function(z) (z[,2]-Nsect[[1]][,2])/sum(Nsect[[1]][,2])), 3)
    #max(abs(round(100*sapply(Nsect, function(z) (z[,2]-Nsect[[1]][,2])/sum(Nsect[[1]][,2])), 3)))

    Nref <- Nsect[[1]][,2]
    Diffs <- sapply(Nsect, function(z) z[,1]-Nref)
    Ahf <- Ahf[rownames(Diffs)]

    Df2 <- cbind(Native=0,Diffs[,-1])
    Tots <- cbind(#All=Diffs[,1],
        Conversion=diag(Df2), # direct
        Disturbance=colSums(Df2)-diag(Df2)) # indirect
    Tots <- 100*Tots/sum(Nref)
    U <- (Tots/100) / (Ahf/sum(Ahf))
    attr(Tots, "Nref") <- sum(Nref)
    attr(Tots, "synergy") <- 100*(sum(Diffs[,1])-sum(Df2))/sum(Nref)

    sector_res[[spp]] <- list(diff=Diffs, total=Tots, unit=U, ref=Nref, ahf=Ahf)

    if (FALSE) {
    ## treating all indirects (over native and not)
    Tots <- cbind(All=Diffs[,1],
        Conversion=c(0, diag(Diffs[-1,-1])), # direct
        Disturbance=colSums(Diffs)-c(0, diag(Diffs[-1,-1]))) # indirect
    Tots["Native","Disturbance"] <- sum(Tots[,"All"])-sum(Tots[-1,-1]) # synergy
    Tots <- rbind(Tots, Total=colSums(Tots))
    Tots1 <- Tots
    ## treating only native indirects
    Tots <- cbind(All=Diffs[,1],
        Conversion=c(0, diag(Diffs[-1,-1])), # direct
        Disturbance=c(0, Diffs[1,-1])) # indirect
    Tots["Native","Disturbance"] <- sum(Tots[,"All"])-sum(Tots[-1,-1]) # synergy
    Tots <- rbind(Tots, Total=colSums(Tots))
    Tots2 <- Tots
    Tots1 <- 100*Tots1/sum(Nref)
    attr(Tots1, "Nref") <- sum(Nref)

    sector_res[[spp]] <- Tots1
    }
}

save(sector_res, Ahf, file=file.path(INDIR, "sector_res.RData"))

## looking at results

load(file.path(INDIR, "sector_res.RData"))

## total
op <- par(mfrow=c(3,1))
z1 <- t(sapply(sector_res, function(z) z$total[-1,1]))
boxplot(z1, main="Direct", ylim=c(-20,20))
abline(h=0, col=4)
z2 <- t(sapply(sector_res, function(z) z$total[-1,2]))
boxplot(z2, main="Indirect", ylim=c(-20,20))
abline(h=0, col=4)
z3 <- t(sapply(sector_res, function(z) rowSums(z$total[-1,])))
boxplot(z3, main="Combined", ylim=c(-20,20))
abline(h=0, col=4)
par(op)

zzz <- matrix(NA, nrow(z1), 3*7+6)
for (i in 1:3) {
    for (j in 1:7) {
        if (i==1)
            z <- z1
        if (i==2)
            z <- z2
        if (i==3)
            z <- z3
        if (j==1)
            k <- 1:3
        if (j==2)
            k <- 5:7
        if (j==3)
            k <- 9:11
        if (j==4)
            k <- 13:15
        if (j==5)
            k <- 17:19
        if (j==6)
            k <- 21:23
        if (j==7)
            k <- 25:27
        zzz[,k[i]] <- z[,j]
    }
}
summary(zzz)

plot(0, type="n", xlim=c(0.5,27.5), ylim=c(-20,20), axes=FALSE, ann=FALSE)
axis(2)
title(ylab="Regional total effect")
abline(h=0,col='#6a3d9a')#;abline(h=c(-0.5,0.5),col='#6a3d9a',lty=2)
boxplot(zzz, add=TRUE,axes=FALSE, #range=0,
    col=c('#1f78b4', '#e31a1c', '#fb9a99', '#a6cee3'))
axis(1, c(2,6,10,14,18,22,26), c("Agr", "EnLin", "EnMW", "For", "Misc",
    "RurUrb", "Transp"))
box()
legend("bottomleft",bty="n",pch=19,col=c('#1f78b4', '#e31a1c', '#fb9a99'),
    legend=c("direct","indirect","combined"))

## unit
fu <- function(x) sign(x) * plogis(log(abs(x)))
ppp <- (sector_res[[1]]$ahf/sum(sector_res[[1]]$ahf))

op <- par(mfrow=c(3,1))
z1 <- fu(t(sapply(sector_res, function(z) z$unit[-1,1])))
boxplot(z1, main="Direct", ylim=c(-1,1), col="gold")
abline(h=0,col=2);abline(h=c(-0.5,0.5),col=2,lty=2)
z2 <- fu(t(sapply(sector_res, function(z) z$unit[-1,2])))
boxplot(z2, main="Indirect", ylim=c(-1,1), col="gold")
abline(h=0,col=2);abline(h=c(-0.5,0.5),col=2,lty=2)
z3 <- fu(t(sapply(sector_res, function(z) (rowSums(z$total[-1,])/100)/ppp[-1])))
boxplot(z3, main="Combined", ylim=c(-1,1), col="gold")
abline(h=0,col=2);abline(h=c(-0.5,0.5),col=2,lty=2)
par(op)

zzz <- matrix(NA, nrow(z1), 3*7+6)
for (i in 1:3) {
    for (j in 1:7) {
        if (i==1)
            z <- z1
        if (i==2)
            z <- z2
        if (i==3)
            z <- z3
        if (j==1)
            k <- 1:3
        if (j==2)
            k <- 5:7
        if (j==3)
            k <- 9:11
        if (j==4)
            k <- 13:15
        if (j==5)
            k <- 17:19
        if (j==6)
            k <- 21:23
        if (j==7)
            k <- 25:27
        zzz[,k[i]] <- z[,j]
    }
}
summary(zzz)

plot(0, type="n", xlim=c(0.5,27.5), ylim=c(-1,1), axes=FALSE, ann=FALSE)
axis(2)
title(ylab="Scaled unit effect")
abline(h=0,col='#6a3d9a');abline(h=c(-0.5,0.5),col='#6a3d9a',lty=2)
boxplot(zzz, add=TRUE,axes=FALSE, col=c('#1f78b4', '#e31a1c', '#fb9a99', '#a6cee3'))
axis(1, c(2,6,10,14,18,22,26), c("Agr", "EnLin", "EnMW", "For", "Misc",
    "RurUrb", "Transp"))
box()
legend("bottomleft",bty="n",pch=19,col=c('#1f78b4', '#e31a1c', '#fb9a99'),
    legend=c("direct","indirect","combined"))

## interaction effects
effs <- data.frame(
    interact = sapply(sector_res, function(z) attr(z$total, "synergy")),
    direct = sapply(sector_res, function(z) sum(z$total[,1])),
    indirect = sapply(sector_res, function(z) sum(z$total[,2])))
effs$combined <- effs$direct + effs$indirect
effs$total <- effs$combined + effs$interact

plot(effs[abs(effs$total)<=100,])
summary(effs)


## species

Plot <- function(spp) {

    Col <- c('#d95f02','#1b9e77','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')
    f <- function(v, ...) {
        v[v > 3] <- 3.5
        v[v < -3] <- -3.5
        barplot(v, width=ppp[-1], space=0, border=NA, col=Col, axes=FALSE, ann=FALSE,
            names.arg=NA, ylim=c(-3,3), ...)
        abline(h=0, col=4)
        abline(h=c(-1,1), col=4, lty=2)
        box(col="lightgrey")
        axis(2)
        invisible()
    }
    op <- par(mfrow=c(2,2))
    f(sector_res[[spp]]$unit[-1,1], main="Direct") # direct
    f(sector_res[[spp]]$unit[-1,2], main="Indirect") # indirect
    f((rowSums(sector_res[[spp]]$total[-1,])/100)/ppp[-1], main="Combined") # combined
    plot.new()
    title(sub=spp)
    legend("topleft", bty="n", fill=Col, legend=rownames(z), border=NA,
        text.col=Col, cex=1.5)
    par(op)

    invisible()
}

pdf(file.path(INDIR, "species.pdf"), onefile=TRUE, width=8, height=8)
for (spp in names(sector_res))
    Plot(spp)
dev.off()
