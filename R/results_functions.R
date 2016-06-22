loadSPP <- function(path) {
    if (!file.exists(path)) {
        cat("species results file does not exist\n")
        return(NULL)
    }
    e <- new.env()
    load(path, envir=e)
    e$res
}

fstat <- function(x, level=0.95) {
    c(Mean=mean(x), Median=median(x), quantile(x, c((1-level)/2, 1 - (1-level)/2)))
}
fstatv <- function(x, level=0.95) {
    out <- data.frame(rowMeans(x),
        apply(x, 1, sd),
        x[,1],
        t(pbapply(x, 1, quantile, c(0.5, (1-level)/2, 1 - (1-level)/2))))
    colnames(out) <- c("Mean","SD","Run1","Median","LCL","UCL")
    out
}

## predict coefficients
predX_age_cc <-
function(ht=c("Decid", "Mixwood", "Conif", "Pine", "BSpr", "Larch"),
ages=c(0,10,20,40,60,80,100,120,140), CC=FALSE, Xn)
{
    ht <- match.arg(ht)
    if (CC)
        ages <- ages[ages <= 60]
    ages <- sort(ages)
    wtage <- ages / 200
    X <- Xn[1:length(ages),]
    X[,-1] <- 0
    X[,"wtAge05"] <- sqrt(wtage)
    X[,"wtAge"] <- wtage
    X[,"wtAge2"] <- wtage^2

    MAXFOR <- 50/200
    fCC1 <- if (CC) pmax(0, 1 - (wtage/MAXFOR)) else 0
    ## Dave`s recovery trajectories
    age <- c(0, 1:20*4) / 200
    conif <- 1-c(0, 1.3, 4.7, 10, 17.3, 26, 35.5, 45.3, 54.6, 63.1, 70.7, 77.3,
        82.7, 87, 90.1, 92.3, 94, 95.3, 96.7, 98.2, 100)/100
    decid <- 1-c(0, 6.5, 15.1, 25.2, 36.1, 47.2, 57.6, 66.7, 74.3, 80.4, 85,
        88.3, 90.5, 92, 93, 94, 95.1, 96.4, 97.6, 98.8, 100)/100

    fCC2 <- if (ht %in% c("Decid", "Mixwood"))
        approxfun(age, decid)(wtage) else approxfun(age, conif)(wtage)
    fCC2[is.na(fCC2)] <- 0
    if (!CC)
        fCC2 <- 0

#    X[,"fCC1"] <- fCC1
    X[,"fCC2"] <- fCC2

#    if (ht=="Decid") {
#    }
    if (ht=="Mixwood") {
        X[,"hab1Mixwood"] <- 1
        X[,"isMix:wtAge05"] <- sqrt(wtage)
        X[,"isMix:wtAge"] <- wtage
        X[,"isMix:wtAge2"] <- wtage^2
    }
    if (ht=="Conif") {
        X[,"hab1Conif"] <- 1
        X[,"isWSpruce:wtAge05"] <- sqrt(wtage)
        X[,"isWSpruce:wtAge"] <- wtage
        X[,"isWSpruce:wtAge2"] <- wtage^2
        X[,"isCon:wtAge05"] <- sqrt(wtage)
        X[,"isCon:wtAge"] <- wtage
        X[,"isCon:wtAge2"] <- wtage^2
        X[,"isUpCon:wtAge05"] <- sqrt(wtage)
        X[,"isUpCon:wtAge"] <- wtage
        X[,"isUpCon:wtAge2"] <- wtage^2
    }
    if (ht=="Pine") {
        X[,"hab1Pine"] <- 1
        X[,"isPine:wtAge05"] <- sqrt(wtage)
        X[,"isPine:wtAge"] <- wtage
        X[,"isPine:wtAge2"] <- wtage^2
        X[,"isCon:wtAge05"] <- sqrt(wtage)
        X[,"isCon:wtAge"] <- wtage
        X[,"isCon:wtAge2"] <- wtage^2
        X[,"isUpCon:wtAge05"] <- sqrt(wtage)
        X[,"isUpCon:wtAge"] <- wtage
        X[,"isUpCon:wtAge2"] <- wtage^2
    }
    if (ht %in% c("Larch", "BSpr")) {
        if (ht=="Larch")
            X[,"hab1Larch"] <- 1
        if (ht=="BSpr")
            X[,"hab1BSpr"] <- 1
        X[,"isBSLarch:wtAge05"] <- sqrt(wtage)
        X[,"isBSLarch:wtAge"] <- wtage
        X[,"isBSLarch:wtAge2"] <- wtage^2
        X[,"isCon:wtAge05"] <- sqrt(wtage)
        X[,"isCon:wtAge"] <- wtage
        X[,"isCon:wtAge2"] <- wtage^2
    }
    rownames(X) <- paste(ht, ifelse(CC, "CC", ""), ages)
    X
}

pred_veghf <-
function(est, Xn, burn_included=TRUE)
{
    X <- Xn[1:13,colnames(est)]
    X[,-1] <- 0 # this should take care of all the modifiers (soft/hard/ARU)
    diag(X) <- 1
    rownames(X) <- c("Decid","Mixwood", "Conif", "Pine", "BSpr", "Larch",
    "Swamp", "WetGrass", "WetShrub", "Shrub", "GrassHerb", "Cult", "UrbInd")

    if (burn_included) {
        X3 <- X[1,,drop=FALSE]
        rownames(X3) <- "Burn"
        X3[1,"hab1bBurn"] <- 1
    }

    X <- rbind(
        predX_age_cc("Conif", CC=FALSE, Xn=Xn),
        predX_age_cc("Pine", CC=FALSE, Xn=Xn),
        predX_age_cc("Decid", CC=FALSE, Xn=Xn),
        predX_age_cc("Mixwood", CC=FALSE, Xn=Xn),
        predX_age_cc("BSpr", CC=FALSE, Xn=Xn),
        predX_age_cc("Larch", CC=FALSE, Xn=Xn),
        predX_age_cc("Conif", CC=TRUE, Xn=Xn),
        predX_age_cc("Pine", CC=TRUE, Xn=Xn),
        predX_age_cc("Decid", CC=TRUE, Xn=Xn),
        predX_age_cc("Mixwood", CC=TRUE, Xn=Xn),
        X[7:13,])
    if (burn_included) {
        X[,"hab1bMixwood"] <- X[,"hab1Mixwood"]
        X[,"hab1bConif"] <- X[,"hab1Conif"]
        X[,"hab1bPine"] <- X[,"hab1Pine"]
        X[,"hab1bBSpr"] <- X[,"hab1BSpr"]
        X[,"hab1bLarch"] <- X[,"hab1Larch"]
        X[,"hab1bWetland"] <- X[,"hab1Wetland"]
        X[,"hab1bShrub"] <- X[,"hab1Shrub"]
        X[,"hab1bGrassHerb"] <- X[,"hab1GrassHerb"]
        X[,"hab1bCult"] <- X[,"hab1Cult"]
        X[,"hab1bUrbInd"] <- X[,"hab1UrbInd"]
        X <- rbind(X, X3)
    }

    pr <- predStat(X, est, level, n=0, ci=TRUE, raw=FALSE)
    out <- pr[,c(1,2,5,6)]
    ## Linear features
    MEAN <- mean(out[,"Median"])
#    Soft <- quantile(MEAN * exp(0.1*est[,"SoftLin_PC"]), c(0.5, (1-level)/2, 1-(1-level)/2))
#    Hard <- quantile(MEAN * exp(est[,"ROAD01"]), c(0.5, (1-level)/2, 1-(1-level)/2))
    Soft <- quantile(exp(est[,"SoftLin_PC"]), c(0.5, (1-level)/2, 1-(1-level)/2))
    Hard <- quantile(exp(est[,"ROAD01"]), c(0.5, (1-level)/2, 1-(1-level)/2))
    attr(out, "linear") <- c(Baseline=MEAN, Soft=Soft, Hard=Hard)
    ## burn should not be shown when it is not selected (i.e. when sum == 0)
    ## REALLY: burn should be just part of young age class, and not being on its own
    if (burn_included)
        attr(out, "burn") <- sum(abs(est[,"hab1bBurn"]))
    rownames(out) <- gsub("Conif", "WhiteSpruce", rownames(out))
    rownames(out) <- gsub("Decid", "Deciduous", rownames(out))
    rownames(out) <- gsub("BSpr", "BlackSpruce", rownames(out))
    rownames(out) <- gsub("Mixwood", "Mixedwood", rownames(out))
    out
}


pred_soilhf <-
function(est, Xn)
{
    X <- Xn[1:5,colnames(est)]
    X[,-1] <- 0
    diag(X) <- 1
    X <- X[c(1,3,3,2,4,5),]
    rownames(X) <- c("Productive", "Clay", "Saline", "RapidDrain", "Cult", "UrbInd")

    X0 <- X1 <- X
    X1[,"pAspen"] <- 1

    pr0 <- predStat(X0, est, level, n=0, ci=TRUE, raw=FALSE)
    pr1 <- predStat(X1, est, level, n=0, ci=TRUE, raw=FALSE)
    ## Linear features
    MEAN <- mean(pr0[,"Median"])
    Soft <- quantile(MEAN * exp(0.1*est[,"SoftLin_PC"]), c(0.5, (1-level)/2, 1-(1-level)/2))
    Hard <- quantile(MEAN * exp(est[,"ROAD01"]), c(0.5, (1-level)/2, 1-(1-level)/2))
    list(treed=pr1[,c(1,2,5,6)], nontreed=pr0[,c(1,2,5,6)],
        linear=c(Baseline=MEAN, Soft=Soft, Hard=Hard))
}

## this is required for getting the axes right across figures
fig_veghf_ymax <-
function(pr2)
{
        lci <- pr2[,3]
        uci <- pr2[,4]
        y1 <- pr2[,2]
        ymax <- max(min(max(uci),2*max(y1)),y1*1.02)
        ymax
}
fig_soilhf_ymax <-
function(pr)
{
        labs <- c("Productive", "Clay", "Saline", "RapidDrain",
            "Cult", "UrbInd")
        pr2 <- pr[labs,]
        lci <- pr2[,3]
        uci <- pr2[,4]
        y1 <- pr2[,2]
        x <- 1:6
        ymax <- max(min(max(uci[x]),2*max(y1)),y1*1.02)
        ymax
}
fig_soilhf <-
function(pr, LAB="", ymax=NULL)
{
        labs <- c("Productive", "Clay", "Saline", "RapidDrain",
            "Cult", "UrbInd")
        op <- par(mai=c(1.5,1,0.2,0.3))
        pr2 <- pr[labs,]
        lci <- pr2[,3]
        uci <- pr2[,4]
        y1 <- pr2[,2]
        x <- 1:6
        if (is.null(ymax))
            ymax <- max(min(max(uci[x]),2*max(y1)),y1*1.02)
        space <- c(1,x[-1]-x[-length(x)])-0.9
        col.r <- c(0,0.3,0.5,1,rep(0.2,2))  # The red part of the rgb
        col.g<-c(0.8,0.5,0,0.2,rep(0.2,2))  # The green part
        col.b<-c(0,0.5,0.5,0.2,rep(0.2,2))  # The blue part
        x1 <- barplot(y1[x], space=space, border="white", col=rgb(col.r,col.g,col.b),
            ylim=c(0,ymax), xlim=c(-0.5,7.2), xaxs="i", yaxt="n", ylab="Relative abundance",
            col.lab="grey50", cex.lab=1.2,axisnames=FALSE)[,1]
        ax <- axis(side=2,cex.axis=0.9,col.axis="grey50",col.ticks="grey50",las=2)
        abline(h=ax, col="grey80")
        x1 <- barplot(y1[x], space=space, border="white", col=rgb(col.r,col.g,col.b),
            ylim=c(0,ymax), xlim=c(-0.5,7.2), xaxs="i", yaxt="n", ylab="Relative abundance",
            col.lab="grey50", cex.lab=1.2,axisnames=FALSE, add=TRUE)[,1]
        box(bty="l",col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2),c(lci[i],y1[i]),col="grey90")
            lines(rep(x1[i],2),c(uci[i],y1[i]),col=rgb(col.r[i],col.g[i],col.b[i]))
        }
        mtext(side=1,at=x1[1:4],line=1.4,c("Productive","Clay","Saline","Rapid Drain"),
            col=rgb(col.r,col.g,col.b)[1:4],las=1)
        mtext(side=1,at=x1[5:6],line=0.7,c("Cultivated HF","Urban/Industry HF"),
            col=rgb(col.r,col.g,col.b)[5:6],las=2)
        mtext(side=3,at=0,adj=0,LAB,col="grey30")
        par(op)
    invisible(NULL)
}

fig_veghf <-
function(pr, LAB="", ymax, ylab="Relative abundance")
{

        op <- par(mai=c(1.5,1,0.2,0.3))
        labs <- c(
            "WhiteSpruce  0", "WhiteSpruce  10",
            "WhiteSpruce  20", "WhiteSpruce  40", "WhiteSpruce  60", "WhiteSpruce  80",
            "WhiteSpruce  100", "WhiteSpruce  120", "WhiteSpruce  140",

            "Pine  0", "Pine  10", "Pine  20", "Pine  40", "Pine  60", "Pine  80", "Pine  100",
            "Pine  120", "Pine  140",

            "Deciduous  0", "Deciduous  10", "Deciduous  20", "Deciduous  40",
            "Deciduous  60", "Deciduous  80", "Deciduous  100", "Deciduous  120",
            "Deciduous  140",

            "Mixedwood  0", "Mixedwood  10", "Mixedwood  20",
            "Mixedwood  40", "Mixedwood  60", "Mixedwood  80", "Mixedwood  100",
            "Mixedwood  120", "Mixedwood  140",

            "BlackSpruce  0", "BlackSpruce  10",
            "BlackSpruce  20", "BlackSpruce  40", "BlackSpruce  60", "BlackSpruce  80",
            "BlackSpruce  100", "BlackSpruce  120", "BlackSpruce  140",

            "Larch  0", "Larch  10", "Larch  20", "Larch  40", "Larch  60",
            "Larch  80", "Larch  100", "Larch  120", "Larch  140",

            "GrassHerb", "Shrub", "Swamp", "WetGrass", "WetShrub",
            "Cult", "UrbInd",

            "WhiteSpruce CC 0", "WhiteSpruce CC 10", "WhiteSpruce CC 20", "WhiteSpruce CC 40", "WhiteSpruce CC 60",
            "Pine CC 0", "Pine CC 10", "Pine CC 20", "Pine CC 40", "Pine CC 60",
            "Deciduous CC 0", "Deciduous CC 10", "Deciduous CC 20", "Deciduous CC 40", "Deciduous CC 60",
            "Mixedwood CC 0", "Mixedwood CC 10", "Mixedwood CC 20", "Mixedwood CC 40", "Mixedwood CC 60")

        pr2 <- pr[labs,]
        lci <- pr2[,3]
        uci <- pr2[,4]
        y1 <- pr2[,2]
        if (missing(ymax))
            ymax <- min(max(uci),2*max(y1))
        #x <- c(rep(1:9,6)+rep(seq(0,50,10),each=9), 61,63,65, 68,70)
        x <- c(rep(1:9,6)+rep(seq(0,50,10),each=9), 61,63,65,67,69, 72,74)
        space <- c(1,x[-1]-x[-length(x)])-0.99  # The spacing between bars
        col.r <- c(rep(0,9),seq(0.3,0.6,length.out=9),seq(0.5,1,length.out=9),
            seq(0.8,0.9,length.out=9),rep(0,9),rep(0,9),
            0.8,0.2,0,0,0, rep(0.2,2))  # The red part
        col.g <- c(seq(0.5,1,length.out=9),seq(0.4,0.8,length.out=9),seq(0.1,0.2,length.out=9),
            seq(0.4,0.8,length.out=9),seq(0.4,0.7,length.out=9),seq(0.15,0.5,length.out=9),
            0.8,0.8,0,0,0, rep(0.2,2))  # The green part
        col.b <- c(rep(0,9),rep(0,9),rep(0,9),seq(0.2,0.4,length.out=9),
            seq(0.2,0.6,length.out=9),seq(0.4,0.7,length.out=9),
            0,0,1,1,1, rep(0.2,2))  # The blue part
        idx <- 1:length(x)
        x1 <- barplot(y1[idx],
            space=space,
            border="white",
            col=rgb(col.r,col.g,col.b),
            ylim=c(0,ymax),
            #xlim=c(-0.5,81.5),
            xlim=c(-0.5,75.5),
            xaxs="i", yaxt="n",
            ylab=ylab,
            col.lab="grey50",
            cex.lab=1.2,axisnames=FALSE)[,1]
        ax <- axis(side=2,cex.axis=0.9,col.axis="grey50",col.ticks="grey50",las=2)
        abline(h=ax, col="grey80")
        x1 <- barplot(y1[idx],
            space=space,border="white",col=rgb(col.r,col.g,col.b),ylim=c(0,ymax),
            xaxs="i",yaxt="n",
            #ylab="Relative abundance",
            col.lab="grey50",
            cex.lab=1.2,axisnames=FALSE, add=TRUE)[,1]
        box(bty="l",col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[idx][i], y1[idx][i]),col="grey90")
            lines(rep(x1[i],2), c(uci[idx][i], y1[idx][i]),col=rgb(col.r[i],col.g[i],col.b[i]))
        }
        mtext(side=1,at=x1[c(5,14,23,32,41,50)],line=1.4,
            c("Upland Spruce","Pine","Deciduous","Mixedwood","Black Spruce","Larch"),
            col=rgb(col.r[c(5,14,23,32,41,50)],col.g[c(5,14,23,32,41,50)],
            col.b[c(5,14,23,32,41,50)]),las=1)
        at1<-rep(seq(1,9,2),6)+rep(c(0,9,18,27,36,45),each=5)
        mtext(side=1,at=x1[at1]-0.3,rep(c("0","20","60","100","140"),6),
            line=0.2,adj=0.5,cex=0.8,col=rgb(col.r[at1],col.g[at1],col.b[at1]))
        mtext(side=1,at=-0.25,adj=1,line=0.2,"Age:",col="grey40",cex=0.8)
        mtext(side=3,at=0,adj=0,LAB,col="grey30")
        mtext(side=1,at=x1[c(55,56,57,58,59)],
            #c("Grass","Shrub","Wetland"),
            c("GrassHerb", "Shrub", "Swamp", "WetGrass", "WetShrub"),
            col=rgb(col.r[c(55,56,57,58,59)],col.g[c(55,56,57,58,59)],
            col.b[c(55,56,57,58,59)]),
            las=2,adj=1.1)
        mtext(side=1,at=x1[c(60,61)],c("Cultivated HF","Urban/Industry HF"),
            col=rgb(col.r[c(60,61)],col.g[c(60,61)],col.b[c(60,61)]),las=2,adj=1.1)

        ## Add cutblock trajectories - upland conifer
        i1<-which(names(y1)=="WhiteSpruce CC 0"):which(names(y1)=="WhiteSpruce CC 60")
        x2<-x1[1:5]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="WhiteSpruce  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Pine
        i1<-which(names(y1)=="Pine CC 0"):which(names(y1)=="Pine CC 60")
        x2<-x1[10:15]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Pine  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Deciduous
        i1<-which(names(y1)=="Deciduous CC 0"):which(names(y1)=="Deciduous CC 60")
        x2<-x1[19:24]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Deciduous  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Mixed
        i1<-which(names(y1)=="Mixedwood CC 0"):which(names(y1)=="Mixedwood CC 60")
        x2<-x1[28:33]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey60")
        x3 <- which(names(y1)=="Mixedwood  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=2)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")

        par(op)

    invisible(pr2)
}
fig_veghf_black <-
function(pr, LAB="", ymax, ylab="Relative abundance")
{

        op <- par(mai=c(1.5,1,0.2,0.3))
        labs <- c(
            "WhiteSpruce  0", "WhiteSpruce  10",
            "WhiteSpruce  20", "WhiteSpruce  40", "WhiteSpruce  60", "WhiteSpruce  80",
            "WhiteSpruce  100", "WhiteSpruce  120", "WhiteSpruce  140",

            "Pine  0", "Pine  10", "Pine  20", "Pine  40", "Pine  60", "Pine  80", "Pine  100",
            "Pine  120", "Pine  140",

            "Deciduous  0", "Deciduous  10", "Deciduous  20", "Deciduous  40",
            "Deciduous  60", "Deciduous  80", "Deciduous  100", "Deciduous  120",
            "Deciduous  140",

            "Mixedwood  0", "Mixedwood  10", "Mixedwood  20",
            "Mixedwood  40", "Mixedwood  60", "Mixedwood  80", "Mixedwood  100",
            "Mixedwood  120", "Mixedwood  140",

            "BlackSpruce  0", "BlackSpruce  10",
            "BlackSpruce  20", "BlackSpruce  40", "BlackSpruce  60", "BlackSpruce  80",
            "BlackSpruce  100", "BlackSpruce  120", "BlackSpruce  140",

            "Larch  0", "Larch  10", "Larch  20", "Larch  40", "Larch  60",
            "Larch  80", "Larch  100", "Larch  120", "Larch  140",

            "GrassHerb", "Shrub", "Swamp", "WetGrass", "WetShrub",
            "Cult", "UrbInd",

            "WhiteSpruce CC 0", "WhiteSpruce CC 10", "WhiteSpruce CC 20", "WhiteSpruce CC 40", "WhiteSpruce CC 60",
            "Pine CC 0", "Pine CC 10", "Pine CC 20", "Pine CC 40", "Pine CC 60",
            "Deciduous CC 0", "Deciduous CC 10", "Deciduous CC 20", "Deciduous CC 40", "Deciduous CC 60",
            "Mixedwood CC 0", "Mixedwood CC 10", "Mixedwood CC 20", "Mixedwood CC 40", "Mixedwood CC 60")

        lwd <- 1
        pr2 <- pr[labs,]
        lci <- pr2[,3]
        uci <- pr2[,4]
        y1 <- pr2[,2]
        if (missing(ymax))
            ymax <- min(max(uci),2*max(y1))
        #x <- c(rep(1:9,6)+rep(seq(0,50,10),each=9), 61,63,65, 68,70)
        x <- c(rep(1:9,6)+rep(seq(0,50,10),each=9), 61,63,65,67,69, 72,74)
        space <- c(1,x[-1]-x[-length(x)])-0.99  # The spacing between bars
        col.r <- c(rep(0,9),seq(0.3,0.6,length.out=9),seq(0.5,1,length.out=9),
            seq(0.8,0.9,length.out=9),rep(0,9),rep(0,9),
            0.8,0.2,0,0,0, rep(0.2,2))  # The red part
        col.g <- c(seq(0.5,1,length.out=9),seq(0.4,0.8,length.out=9),seq(0.1,0.2,length.out=9),
            seq(0.4,0.8,length.out=9),seq(0.4,0.7,length.out=9),seq(0.15,0.5,length.out=9),
            0.8,0.8,0,0,0, rep(0.2,2))  # The green part
        col.b <- c(rep(0,9),rep(0,9),rep(0,9),seq(0.2,0.4,length.out=9),
            seq(0.2,0.6,length.out=9),seq(0.4,0.7,length.out=9),
            0,0,1,1,1, rep(0.2,2))  # The blue part
        idx <- 1:length(x)
        x1 <- barplot(y1[idx],
            space=space,
            border="white",
            col=rgb(col.r,col.g,col.b),
            ylim=c(0,ymax),
            #xlim=c(-0.5,81.5),
            xlim=c(-0.5,75.5),
            xaxs="i", yaxt="n",
            ylab=ylab,
            col.lab=1,
            cex.lab=1.2,axisnames=FALSE)[,1]
        ax <- axis(side=2,cex.axis=0.9,col.axis=1,col.ticks=1,las=2)
        #abline(h=ax, col=1)
        x1 <- barplot(y1[idx],
            space=space,border="white",col=rgb(col.r,col.g,col.b),ylim=c(0,ymax),
            xaxs="i",yaxt="n",
            #ylab="Relative abundance",
            col.lab=1,
            cex.lab=1.2,axisnames=FALSE, add=TRUE)[,1]
        box(bty="o",col=1)
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[idx][i], y1[idx][i]),col="white", lwd=lwd)
            lines(rep(x1[i],2), c(uci[idx][i], y1[idx][i]),col=rgb(col.r[i],col.g[i],col.b[i]), lwd=lwd)
        }
        mtext(side=1,at=x1[c(5,14,23,32,41,50)],line=1.4,
            c("Upland Spruce","Pine","Deciduous","Mixedwood","Black Spruce","Larch"),
            col=rgb(col.r[c(5,14,23,32,41,50)],col.g[c(5,14,23,32,41,50)],
            col.b[c(5,14,23,32,41,50)]),las=1)
        at1<-rep(seq(1,9,2),6)+rep(c(0,9,18,27,36,45),each=5)
        mtext(side=1,at=x1[at1]-0.3,rep(c("0","20","60","100","140"),6),
            line=0.2,adj=0.5,cex=0.8,col=rgb(col.r[at1],col.g[at1],col.b[at1]))
        mtext(side=1,at=-0.25,adj=1,line=0.2,"Age:",col="grey10",cex=0.8)
        mtext(side=3,at=0,adj=0,LAB,col=1)
        mtext(side=1,at=x1[c(55,56,57,58,59)],
            #c("Grass","Shrub","Wetland"),
            c("GrassHerb", "Shrub", "Swamp", "WetGrass", "WetShrub"),
            col=rgb(col.r[c(55,56,57,58,59)],col.g[c(55,56,57,58,59)],
            col.b[c(55,56,57,58,59)]),
            las=2,adj=1.1)
        mtext(side=1,at=x1[c(60,61)],c("Cultivated HF","Urban/Industry HF"),
            col=rgb(col.r[c(60,61)],col.g[c(60,61)],col.b[c(60,61)]),las=2,adj=1.1)

        ## Add cutblock trajectories - upland conifer
        i1<-which(names(y1)=="WhiteSpruce CC 0"):which(names(y1)=="WhiteSpruce CC 60")
        x2<-x1[1:5]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey10", lwd=lwd)
        x3 <- which(names(y1)=="WhiteSpruce  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=1, lwd=lwd)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Pine
        i1<-which(names(y1)=="Pine CC 0"):which(names(y1)=="Pine CC 60")
        x2<-x1[10:15]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey10", lwd=lwd)
        x3 <- which(names(y1)=="Pine  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=1, lwd=lwd)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Deciduous
        i1<-which(names(y1)=="Deciduous CC 0"):which(names(y1)=="Deciduous CC 60")
        x2<-x1[19:24]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey10", lwd=lwd)
        x3 <- which(names(y1)=="Deciduous  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=1, lwd=lwd)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")
        ## Mixed
        i1<-which(names(y1)=="Mixedwood CC 0"):which(names(y1)=="Mixedwood CC 60")
        x2<-x1[28:33]+0.15*(x1[2]-x1[1])
        for (j in 1:5)
            lines(rep(x2[j],2),c(lci[i1[j]],uci[i1[j]]),col="grey10", lwd=lwd)
        x3 <- which(names(y1)=="Mixedwood  80")
        lines(c(x2[1:5], x1[x3]),y1[c(i1, x3)],col="grey30", lty=1, lwd=lwd)
        points(x2[1:5],y1[i1],pch=18,cex=1,col="grey30")
        points(x2[1:5],y1[i1],pch=5,cex=0.7,col="grey10")

        par(op)

    invisible(pr2)
}

rank_fun <- function(x, l, u, n=1, col=NULL, lab=NULL) {
    p <- data.frame(x=x, l=l, u=u, n=n)
    rownames(p) <- lab
    if (is.null(col))
        col <- 1
    p$col <- col
    p <- p[!is.na(p$x),]
    p <- p[order(p$x),]
    cex <- log(p$n) / 10
    n <- nrow(p)

    op <- par(mar = c(6, 4, 2, 2)+0.1, las=2, cex.axis=0.6)
    plot(1:n, p$x, ann=FALSE, type="n", axes=FALSE, ylim=c(-100,100))
    abline(h=0, col="red4")
    segments(1:n, y0=p$l, y1=p$u, col=p$col)
    points(1:n, p$x, pch=19, cex=cex, col=p$col)
    axis(1, at=1:n, tick=FALSE, labels=rownames(p))
    axis(2)
    par(op)

    invisible(NULL)
}

pred_hf <-
function(est, Xn, n=200, fillin=0, remn=NULL)
{
    cc <- paste(c("THF", "Lin", "Nonlin", "Succ", "Alien", "Noncult", "Cult",
        "THF2", "Nonlin2", "Succ2", "Alien2", "Noncult2"),
        "KM", sep="_")
    if (!is.null(remn))
        #cc <- c(cc, paste(c("Remn", "Remn2"), "KM", sep="_"))
        cc <- c(cc, remn)
    ehf <- est
    ehf[] <- 0
    ehf[,cc] <- est[,cc]

        pr <- list()
        Types <- c("Cult", "UrbInd", "HFor", "SoftLin", "HardLin")
        if (!is.null(remn))
            Types <- c("Remn", Types)
        for (type in Types) {
            Range <- if (type %in% c("SoftLin", "HardLin"))
                c(0, 0.25) else c(0, 1)
            hf <- seq(Range[1], Range[2], len=n)
            Xhf <- Xn[1:n,] # matrix(0, n, length(cc))
            Xhf[] <- 0
            rownames(Xhf) <- NULL
            if (!is.null(remn) & type == "Remn") {
                Xhf[,remn] <- hf
                #Xhf[,"Remn_KM"] <- hf
                #Xhf[,"Remn2_KM"] <- hf^2
            }
            if (type == "Cult") {
                Xhf[,"THF_KM"] <- hf
                Xhf[,"THF2_KM"] <- hf^2
                Xhf[,"Alien_KM"] <- hf
                Xhf[,"Alien2_KM"] <- hf^2
                Xhf[,"Cult_KM"] <- hf
                Xhf[,"Nonlin_KM"] <- hf
                Xhf[,"Nonlin2_KM"] <- hf^2
                if (fillin) {
                    Xhf[,remn] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
                }
            }
            if (type == "UrbInd") {
                Xhf[,"THF_KM"] <- hf
                Xhf[,"THF2_KM"] <- hf^2
                Xhf[,"Alien_KM"] <- hf
                Xhf[,"Alien2_KM"] <- hf^2
                Xhf[,"Noncult_KM"] <- hf
                Xhf[,"Noncult2_KM"] <- hf^2
                Xhf[,"Nonlin_KM"] <- hf
                Xhf[,"Nonlin2_KM"] <- hf^2
                if (fillin) {
                    Xhf[,remn] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
                }
            }
            if (type == "HFor") {
                Xhf[,"THF_KM"] <- hf
                Xhf[,"THF2_KM"] <- hf^2
                Xhf[,"Succ_KM"] <- hf
                Xhf[,"Succ2_KM"] <- hf^2
                Xhf[,"Noncult_KM"] <- hf
                Xhf[,"Noncult2_KM"] <- hf^2
                Xhf[,"Nonlin_KM"] <- hf
                Xhf[,"Nonlin2_KM"] <- hf^2
                if (fillin) {
                    Xhf[,remn] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
                }
            }
            if (type == "SoftLin") {
                Xhf[,"THF_KM"] <- hf
                Xhf[,"THF2_KM"] <- hf^2
                Xhf[,"Succ_KM"] <- hf
                Xhf[,"Succ2_KM"] <- hf^2
                Xhf[,"Noncult_KM"] <- hf
                Xhf[,"Noncult2_KM"] <- hf^2
                Xhf[,"Lin_KM"] <- hf
                if (fillin) {
                    Xhf[,remn] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
                }
            }
            if (type == "HardLin") {
                Xhf[,"THF_KM"] <- hf
                Xhf[,"THF2_KM"] <- hf^2
                Xhf[,"Alien_KM"] <- hf
                Xhf[,"Alien2_KM"] <- hf^2
                Xhf[,"Noncult_KM"] <- hf
                Xhf[,"Noncult2_KM"] <- hf^2
                Xhf[,"Lin_KM"] <- hf
                if (fillin) {
                    Xhf[,remn] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    #Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
                }
            }
        #prst <- predStat(Xhf, ehf, level, n=0, ci=TRUE, raw=FALSE)[,c("Median", "CL1q","CL2q")]
        #prst[,1] <- prst[,1] / prst[1,1]
        #prst[,2] <- prst[,2] / prst[1,1]
        #prst[,2] <- prst[,2] / prst[1,1]
        prst <- t(exp(predStat(Xhf, ehf, raw=TRUE)))
        Mean <- colMeans(prst)
        prst <- t(apply(t(prst / prst[,1]), 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
        colnames(prst) <- c("Median", "CL1q","CL2q")
        pr[[type]] <- data.frame(hf=hf*100, prst, Mean=Mean)
    }
    pr
}

pred_any <-
function(what="WetWaterKM", est, Xn, n=200)
{

    ehf <- est
    ehf[] <- 0
    ehf[,what] <- est[,what]

    hf <- seq(0, 1, len=n)
    Xhf <- Xn[1:n,] # matrix(0, n, length(cc))
    Xhf[] <- 0
    rownames(Xhf) <- NULL
    Xhf[,what] <- hf

    prst <- t(exp(predStat(Xhf, ehf, raw=TRUE)))
    Mean <- colMeans(prst)
    prst <- t(apply(t(prst / prst[,1]), 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
    colnames(prst) <- c("Median", "CL1q","CL2q")

    data.frame(hf=hf*100, prst, Mean=Mean)
}

fig_hf_noremn <-
function(est, Xn, fillin=0,LAB="",remn=NULL)
{
    pr <- pred_hf(est, Xn, fillin=fillin,remn=remn)
    pr <- pr[c("HFor", "Cult", "UrbInd", "SoftLin", "HardLin")]
    ymax <- max(sapply(pr, function(z) max(z[,-1])))
    mmax <- min(max(max(sapply(pr, function(z) max(z[,2]))), 5), 10)
    ymax <- min(mmax, ymax)

    Colb <- c("#00FF0015", "#DD000010", "#D0D00025", "#80008015", "#0000DD10")
    Coll <- c(rgb(0.2,0.6,0.2), "red3", "orange3", "#A000A0", "blue3")

    #op <- par(mfrow=c(1,1))

    plot(pr[["Cult"]]$hf, pr[["Cult"]]$Median,
        xlab="", ylim=c(0,ymax),
        cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",typ="n",bty="n",
        ylab="Relative abundance",col.lab="grey50",cex.lab=1.2)
    axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80")
    axis(side=2,at=seq(0,ymax,0.5),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80",las=2)
    axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60")
    axis(side=2,at=seq(0,ymax,0.5),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60",las=2)
    mtext(side=1,at=50,adj=0.5,"Surrounding footprint (%)",line=2,col="grey50",cex=1.2)
    box(col="grey60")
    mtext(side=3,at=0,adj=0,LAB,col="grey30")

    for (i in 1:length(pr)) {
        polygon(c(pr[[i]]$hf, rev(pr[[i]]$hf)),
            c(pr[[i]][,"CL1q"], rev(pr[[i]][,"CL2q"])),
            col=Colb[i], border=NA)  # Transparent colours so they overlap
    }
    for (i in 1:length(pr))
        lines(pr[[i]]$hf, pr[[i]][,"Median"], col=Coll[i],lwd=2)
    text(rep(0.1, 5), c(0.9, 0.75, 0.6, 0.45, 0.3), names(pr), pos=4, col=Coll)

    #par(op)

    invisible(NULL)
}

fig_any <-
function(what, est, Xn, LAB="", xlab="Percent")
{
    pr <- pred_any(what, est, Xn)
    ymax <- max(pr[,-1])
    mmax <- min(max(max(max(pr[,2])), 5), 10)
    ymax <- min(mmax, ymax)

    Colb <- c("#00FF0015", "#DD000010", "#D0D00025", "#80008015", "#0000DD10")[1]
    Coll <- c(rgb(0.2,0.6,0.2), "red3", "orange3", "#A000A0", "blue3")[1]

    #op <- par(mfrow=c(1,1))

    plot(pr$hf, pr$Median,
        xlab="", ylim=c(0,ymax),
        cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",typ="n",bty="n",
        ylab="Relative abundance",col.lab="grey50",cex.lab=1.2)
    axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80")
    axis(side=2,at=seq(0,ymax,0.5),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80",las=2)
    axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60")
    axis(side=2,at=seq(0,ymax,0.5),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60",las=2)
    mtext(side=1,at=50,adj=0.5,xlab,line=2,col="grey50",cex=1.2)
    box(col="grey60")
    mtext(side=3,at=0,adj=0,LAB,col="grey30")

    polygon(c(pr$hf, rev(pr$hf)),
        c(pr[,"CL1q"], rev(pr[,"CL2q"])),
        col=Colb, border=NA)  # Transparent colours so they overlap
    lines(pr$hf, pr[,"Median"], col=Coll,lwd=2)
    #text(rep(0.1, 5), c(0.9, 0.75, 0.6, 0.45, 0.3), names(pr), pos=4, col=Coll)

    #par(op)

    invisible(NULL)
}

fig_any_black <-
function(what, est, Xn, LAB="", xlab="Percent")
{
    pr <- pred_any(what, est, Xn)
    ymax <- max(pr[,-1])
    mmax <- min(max(max(max(pr[,2])), 5), 10)
    ymax <- min(mmax, ymax)

    Colb <- c("#00FF0015", "#DD000010", "#D0D00025", "#80008015", "#0000DD10")[1]
    Coll <- c(rgb(0.2,0.6,0.2), "red3", "orange3", "#A000A0", "blue3")[1]

    #op <- par(mfrow=c(1,1))

    plot(pr$hf, pr$Median,
        xlab="", ylim=c(0,ymax),typ="n",
        #cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n",
        ylab="Relative abundance",col.lab=1,cex.lab=1.2)
#    axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1)
#    axis(side=2,at=seq(0,ymax,0.5),tck=1,cex.axis=1,las=2)
#    axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1)
#    axis(side=2,at=seq(0,ymax,0.5),tck=0.01,las=2)
    mtext(side=1,at=50,adj=0.5,xlab,line=2,cex=1.2)
    box()
    mtext(side=3,at=0,adj=0,LAB,col=1)

    polygon(c(pr$hf, rev(pr$hf)),
        c(pr[,"CL1q"], rev(pr[,"CL2q"])),
        col=Colb, border=NA)  # Transparent colours so they overlap
    lines(pr$hf, pr[,"Median"], col=Coll,lwd=2)
    #text(rep(0.1, 5), c(0.9, 0.75, 0.6, 0.45, 0.3), names(pr), pos=4, col=Coll)

    #par(op)

    invisible(NULL)
}

fig_linear <-
function(pr, LAB)
{
		p.mean <- pr[1]
		#p.softlin10 <- 0.9*pr[1] + 0.1*pr[2]
		p.softlin10 <- pr[1] * exp(0.1*log(pr[2]))
		p.hardlin10 <- pr[1]*pr[5]
		#p.hardlin10 <- 0.9*pr[1] + 0.1*pr[5]
		ymax1<-max(p.softlin10,p.hardlin10,2*p.mean)*1.03
		plot(c(1,1.95,2.05),c(p.mean,p.softlin10,p.hardlin10),pch=c(1,16,15),col=c("grey30","blue3","red4"),xlab="Human footprint",ylab="Relative abundance",xlim=c(0.8,2.8),ylim=c(0,ymax1),tck=0.01,yaxs="i",xaxt="n",yaxt="n",bty="l",cex=2,lwd=2,cex.lab=1.4,cex.axis=1.3,col.lab="grey40")
		axis(side=2,at=pretty(c(0,ymax1),n=5),cex.axis=1.3,tck=0.01,cex.axis=1.3,col.axis="grey40",col.ticks="grey40")
		axis(side=1,at=c(1,2),lab=c("None","10% linear"),tck=0.01,cex.axis=1.3,col.axis="grey40",col.ticks="grey40")
		box(bty="l",col="grey40")
		lines(c(1,1.95),c(p.mean,p.softlin10),col="blue3")
		lines(c(1,2.05),c(p.mean,p.hardlin10),col="red4")
		points(c(1,1.95,2.05),c(p.mean,p.softlin10,p.hardlin10),pch=c(1,16,15),col=c("grey30","blue3","red4"),cex=2,lwd=2)  # Put these back on top of the lines
		ly<-c(p.softlin10,p.hardlin10)  # Label y values - adjust so not overlapping
		if (abs(ly[2]-ly[1])<ymax1/20) ly<-c(mean(ly)+ymax1/40*sign(ly[1]-ly[2]),mean(ly)+ymax1/40*sign(ly[2]-ly[1]))
		text(c(2.15,2.15),ly,c("Soft linear","Hard linear"),col=c("blue3","red4"),cex=1.3,adj=0)
		mtext(side=3,at=0.8,adj=0,LAB,col="grey30",cex=1.3)
    invisible(NULL)
}

