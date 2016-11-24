library(mefa4)
load("e:/peter/AB_data_v2016/out/kgrid/kgrid_table.Rdata")
fv5 <- "e:/peter/AB_data_v2016/out/kgrid/veg-hf_1kmgrid_fix-fire.Rdata"
fv6 <- "e:/peter/AB_data_v2016/data/kgrid-V6/veg-hf_1kmgrid_v6.Rdata"

e <- new.env()
load(fv5, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
v5 <- e$dd1km_pred[1:2]

e <- new.env()
load(fv6, envir=e)
e$dd1km_pred$sample_year
e$dd1km_pred$scale
v6 <- e$dd1km_pred[1:2]

rm(e)

##

compare_sets(colnames(v5[[1]]), colnames(v6[[1]]))
intersect(colnames(v5[[1]]), colnames(v6[[1]]))
setdiff(colnames(v5[[1]]), colnames(v6[[1]]))
setdiff(colnames(v6[[1]]), colnames(v5[[1]]))

#cn <- data.frame(v6=colnames(v6[[1]]))
#cn$v5 <- colnames(v5[[1]])[match(colnames(v6[[1]]), colnames(v5[[1]]))]
#write.csv(cn, row.names=FALSE,
#    file="e:/peter/AB_data_v2016/data/kgrid-V6/xwalk-v5-v6.csv")

compare_sets(colnames(v5[[2]]), colnames(v6[[2]]))
intersect(colnames(v5[[2]]), colnames(v6[[2]]))
setdiff(colnames(v5[[2]]), colnames(v6[[2]]))
setdiff(colnames(v6[[2]]), colnames(v5[[2]]))

xt <- read.csv("e:/peter/AB_data_v2016/data/kgrid-V6/xwalk-v5-v6.csv")
setdiff(colnames(v5[[1]]), xt$v5)
setdiff(colnames(v5[[2]]), xt$v5)
setdiff(colnames(v6[[1]]), xt$v6)
setdiff(colnames(v6[[2]]), xt$v6)

xt$v5 <- as.character(xt$v5)
xt$v6 <- as.character(xt$v6)
xt$v5[is.na(xt$v5)] <- "XXX"
xt$v6[is.na(xt$v6)] <- "XXX"

v5[[1]] <- cBind(v5[[1]], "XXX"=0)
v5[[2]] <- cBind(v5[[2]], "XXX"=0)
v6[[1]] <- cBind(v6[[1]], "XXX"=0)
v6[[2]] <- cBind(v6[[2]], "XXX"=0)

v5[[1]] <- v5[[1]] / ifelse(rowSums(v5[[1]])==0, 1, rowSums(v5[[1]]))
v5[[2]] <- v5[[2]] / ifelse(rowSums(v5[[2]])==0, 1, rowSums(v5[[2]]))
v6[[1]] <- v6[[1]] / ifelse(rowSums(v6[[1]])==0, 1, rowSums(v6[[1]]))
v6[[2]] <- v6[[2]] / ifelse(rowSums(v6[[2]])==0, 1, rowSums(v6[[2]]))

cr5 <- groupSums(v5[[1]], 2, xt$v5[match(colnames(v5[[1]]), xt$v5)])
rf5 <- groupSums(v5[[2]], 2, xt$v5[match(colnames(v5[[2]]), xt$v5)])

cr6 <- groupSums(v6[[1]], 2, xt$v5[match(colnames(v6[[1]]), xt$v6)])
rf6 <- groupSums(v6[[2]], 2, xt$v5[match(colnames(v6[[2]]), xt$v6)])

compare_sets(colnames(cr5), colnames(cr6))
compare_sets(colnames(rf5), colnames(rf6))

intersect(colnames(cr5), colnames(cr6))
setdiff(colnames(cr5), colnames(cr6))
setdiff(colnames(cr6), colnames(cr5))

colSums(cr5[,setdiff(colnames(cr5), colnames(cr6))])
cr5 <- cr5[,colnames(cr6)]
colSums(rf5[,setdiff(colnames(rf5), colnames(rf6))])
rf5 <- rf5[,colnames(rf6)]

#cr5 <- cr5 / ifelse(rowSums(cr5)==0, 1, rowSums(cr5))
#cr6 <- cr6 / ifelse(rowSums(cr6)==0, 1, rowSums(cr6))
#rf5 <- rf5 / ifelse(rowSums(rf5)==0, 1, rowSums(rf5))
#rf6 <- rf6 / ifelse(rowSums(rf6)==0, 1, rowSums(rf6))

if (TRUE) {
z <- colnames(cr5)
z1 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)
z <- colnames(rf5)
z2 <- ifelse(substr(z, nchar(z),nchar(z)) %in% c("R",0:9), substr(z, 1,nchar(z)-1), z)

cr5 <- groupSums(cr5, 2, z1)
cr6 <- groupSums(cr6, 2, z1)
rf5 <- groupSums(rf5, 2, z2)
rf6 <- groupSums(rf6, 2, z2)
}

cr5 <- groupMeans(cr5, 1, kgrid$NSRNAME)
cr6 <- groupMeans(cr6, 1, kgrid$NSRNAME)
rf5 <- groupMeans(rf5, 1, kgrid$NSRNAME)
rf6 <- groupMeans(rf6, 1, kgrid$NSRNAME)

col <- "grey"
cols <- colorRampPalette(c("blue","red"))(5)
pdf("e:/peter/AB_data_v2016/data/kgrid-V6/nsr-level-comparison-noage.pdf",
    width=10, height=5.5, onefile=TRUE)
for (i in 1:ncol(cr6)) {
    j <- colnames(cr6)[i]
    lim <- c(0, max(0.1, cr5[,j], cr6[,j]))
    if (j %in% colnames(rf6))
        lim <- c(0, max(lim, rf5[,j], rf6[,j]))
    op <- par(mfrow=c(1,2))
    d <- ceiling(100*abs(cr5[,j] - cr6[,j]))
    d[d < 1] <- 1
    d[d >= 5] <- 5
    plot(cr5[,j], cr6[,j], main=j,
        xlim=lim, ylim=lim,
        xlab="current 1km mean v5", ylab="current 1km mean v6", col=cols[d], pch=19)
    text(cr5[,j], cr6[,j], ifelse(d > 1, rownames(cr5), ""),
        col=cols[d], pos=3, cex=0.4, xpd=NA)
        abline(0,1,col=col,lty=1)
        abline(-0.01,1,col=col,lty=2)
        abline(0.01,1,col=col,lty=2)
        abline(-0.05,1,col=col,lty=3)
        abline(0.05,1,col=col,lty=3)

    if (j %in% colnames(rf6)) {
        d <- ceiling(100*abs(rf5[,j] - rf6[,j]))
        d[d < 1] <- 1
        d[d >= 5] <- 5
        plot(rf5[,j], rf6[,j], main="",
            xlim=lim, ylim=lim,
            xlab="reference 1km mean v5", ylab="reference 1km mean v6", col=cols[d], pch=19)
        text(rf5[,j], rf6[,j], ifelse(d > 1, rownames(cr5), ""),
            col=cols[d], pos=3, cex=0.4, xpd=NA)
        abline(0,1,col=col,lty=1)
        abline(-0.01,1,col=col,lty=2)
        abline(0.01,1,col=col,lty=2)
        abline(-0.05,1,col=col,lty=3)
        abline(0.05,1,col=col,lty=3)
    } else {
        plot.new()
    }
    par(op)
}
dev.off()

## only v6
crS <- rowSums(v6[[1]][,setdiff(colnames(v6[[1]]), colnames(v5[[1]]))]) / ifelse(rowSums(v6[[1]])==0, 1, rowSums(v6[[1]]))

