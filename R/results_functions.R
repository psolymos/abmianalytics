loadSPP <- function(path) {
    if (!file.exists(path)) {
        cat("species results file does not exist\n")
        return(NULL)
    }
    e <- new.env()
    load(path, envir=e)
    e$res
}

fstat <- function(x) {
    c(Mean=mean(x), Median=median(x), quantile(x, c((1-level)/2, 1 - (1-level)/2)))
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
function(est, Xn)
{
    X <- Xn[1:11,colnames(est)]
    X[,-1] <- 0
    diag(X) <- 1
    rownames(X) <- c("Decid","Mixwood", "Conif", "Pine", "BSpr", "Larch",  
    "Wetland", "Shrub", "GrassHerb", "Cult", "UrbInd")

    X3 <- X[1,,drop=FALSE]
    rownames(X3) <- "Burn"
    X3[1,"hab1bBurn"] <- 1

    X <- rbind(predX_age_cc("Decid", CC=FALSE, Xn=Xn),
        predX_age_cc("Mixwood", CC=FALSE, Xn=Xn),
        predX_age_cc("Conif", CC=FALSE, Xn=Xn),
        predX_age_cc("Pine", CC=FALSE, Xn=Xn),
        predX_age_cc("BSpr", CC=FALSE, Xn=Xn),
        predX_age_cc("Larch", CC=FALSE, Xn=Xn),
        predX_age_cc("Decid", CC=TRUE, Xn=Xn),
        predX_age_cc("Mixwood", CC=TRUE, Xn=Xn),
        predX_age_cc("Conif", CC=TRUE, Xn=Xn),
        predX_age_cc("Pine", CC=TRUE, Xn=Xn),
        X[7:11,])
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

    pr <- predStat(X, est, level, n=0, ci=TRUE, raw=FALSE)
    out <- pr[,c(1,2,5,6)]
    ## Linear features
    MEAN <- mean(out[,"Median"])
    Soft <- MEAN * mean(exp(0.1*est[,"SoftLin_PC"]))
    Hard <- MEAN * mean(exp(est[,"ROAD01"]))
    attr(out, "linear") <- c(Baseline=MEAN, Soft=Soft, Hard=Hard)
    ## burn should not be shown when it is not selected (i.e. when sum == 0)
    ## REALLY: burn should be just part of young age class, and not being on its own
    attr(out, "burn") <- sum(abs(est[,"hab1bBurn"]))
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
#    Soft <- MEAN * mean(exp(0.1*est[,"SoftLin_PC"]))
    Hard <- MEAN * mean(exp(est[,"ROAD01"]))
    list(treed=pr1[,c(1,2,5,6)], nontreed=pr0[,c(1,2,5,6)],
        linear=c(Baseline=MEAN, Soft=MEAN, Hard=Hard))
}

fig_hab <- 
function(spp, burn=TRUE) 
{
    pr <- if (is.character(spp))
        pred_hab(spp) else spp

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

            "LarchFen  0", "LarchFen  10", "LarchFen  20", "LarchFen  40", "LarchFen  60", 
            "LarchFen  80", "LarchFen  100", "LarchFen  120", "LarchFen  140", 

            "Grassland", "Shrubland", "Marsh", "Swamp", "Bog", "Fen", 
            "Cult", "UrbInd", 

            "Deciduous CC 0", "Deciduous CC 10", "Deciduous CC 20", "Deciduous CC 40", "Deciduous CC 60", 
            "Mixedwood CC 0", "Mixedwood CC 10", "Mixedwood CC 20", "Mixedwood CC 40", "Mixedwood CC 60",
            "WhiteSpruce CC 0", "WhiteSpruce CC 10", "WhiteSpruce CC 20", "WhiteSpruce CC 40", "WhiteSpruce CC 60", 
            "Pine CC 0", "Pine CC 10", "Pine CC 20", "Pine CC 40", "Pine CC 60", 
            "Cult", "UrbInd", 
            "Burn")

        pr2 <- pr[labs,]
        lci <- pr2[,3]
        uci <- pr2[,4]
        y1 <- pr2[,2]
        ymax <- min(max(uci),2*max(y1))
        x <- c(rep(1:9,6)+rep(seq(0,50,10),each=9), 61,63,65,67,69,71, 74,76)
        space <- c(1,x[-1]-x[-length(x)])-0.99  # The spacing between bars
        col.r <- c(rep(0,9),seq(0.3,0.6,length.out=9),seq(0.5,1,length.out=9),
            seq(0.8,0.9,length.out=9),rep(0,9),rep(0,9),0.8,0.2,0,0,0,0,rep(0.2,3))  # The red part
        col.g <- c(seq(0.5,1,length.out=9),seq(0.4,0.8,length.out=9),seq(0.1,0.2,length.out=9),
            seq(0.4,0.8,length.out=9),seq(0.4,0.7,length.out=9),seq(0.15,0.5,length.out=9),
            0.8,0.8,0,0,0,0,rep(0.2,3))  # The green part
        col.b <- c(rep(0,9),rep(0,9),rep(0,9),seq(0.2,0.4,length.out=9),
            seq(0.2,0.6,length.out=9),seq(0.4,0.7,length.out=9),0,0,1,0.8,0.6,0.4,rep(0.2,3))  # The blue part
        idx <- 1:62
        x1 <- barplot(y1[idx],
            space=space,border="white",col=rgb(col.r,col.g,col.b),ylim=c(0,ymax),
            #xlim=c(-0.5,81.5),
            xlim=c(-0.5,77.5),
            xaxs="i",yaxt="n",
            ylab="Predicted population density (males/ha)", # ylab="Relative abundance",
            col.lab="grey50",
            cex.lab=1.2,axisnames=FALSE)[,1]
        ax <- axis(side=2,cex.axis=0.9,col.axis="grey50",col.ticks="grey50",las=2)
        abline(h=ax, col="grey80")
        #axis(side=2, tck=0.02, cex.axis=0.9,col.axis="grey50",col.ticks="grey50",
        #    las=2,at=seq(0,ymax,0.2))
        #abline(h=seq(0,ymax,0.2),col="grey80")
        x1 <- barplot(y1[idx],
            space=space,border="white",col=rgb(col.r,col.g,col.b),ylim=c(0,ymax),
            xlim=c(-0.5,83.5),xaxs="i",yaxt="n",
            #ylab="Relative abundance",
            col.lab="grey50",
            cex.lab=1.2,axisnames=FALSE, add=TRUE)[,1]
        box(bty="l",col="grey50")
        for (i in 1:length(x1)) {
            lines(rep(x1[i],2), c(lci[idx][i], y1[idx][i]),col="grey90")
            lines(rep(x1[i],2), c(uci[idx][i], y1[idx][i]),col=rgb(col.r[i],col.g[i],col.b[i]))
        }
        mtext(side=1,at=x1[c(5,14,23,32,41,50)],line=1.4,
            c("Upland Spruce","Pine","Deciduous","Mixedwood","Black Spruce","Larch Fen"),
            col=rgb(col.r[c(5,14,23,32,41,50)],col.g[c(5,14,23,32,41,50)],col.b[c(5,14,23,32,41,50)]),las=1)
        mtext(side=1,at=x1[c(55,56,57,58,59,60)],c("Grass","Shrub","Marsh", "Swamp","Open Bog", "Fen"),
            col=rgb(col.r[c(55,56,57,58,59,60)],col.g[c(55,56,57,58,59,60)],col.b[c(55,56,57,58,59,60)]),
            las=2,adj=1.1)
        mtext(side=1,at=x1[c(61,62)],c("Cultivated HF","Urban/Industry HF"),
            col=rgb(col.r[c(61,62)],col.g[c(61,62)],col.b[c(61,62)]),las=2,adj=1.1)
        at1<-rep(seq(1,9,2),6)+rep(c(0,9,18,27,36,45),each=5)
        mtext(side=1,at=x1[at1]-0.3,rep(c("0","20","60","100","140"),6),
            line=0.2,adj=0.5,cex=0.8,col=rgb(col.r[at1],col.g[at1],col.b[at1]))
        mtext(side=1,at=-0.25,adj=1,line=0.2,"Age:",col="grey40",cex=0.8)
        mtext(side=3,at=0,adj=0,as.character(taxa(mm)[spp,"English_Name"]),col="grey30")
    

        ## Add soft and hard linear 10% effects
        #lines(2*x1[63]-x1[62]+c(-0.2,0.2),  pr2[c(nrow(pr2)-1, nrow(pr2)-0),2],col="grey60")
        #lines(3*x1[63]-2*x1[62]+c(-0.2,0.2),pr2[c(nrow(pr2)-3, nrow(pr2)-2),2],col="grey60")
        #points(2*x1[63]-x1[62]+c(-0.2,0.2),pr2[c(nrow(pr2)-1, nrow(pr2)-0),2],pch="-",cex=2.1,col="grey50")
        #points(3*x1[63]-2*x1[62]+c(-0.2,0.2),pr2[c(nrow(pr2)-3, nrow(pr2)-2),2],pch="-",
        #    cex=2.1,col="grey50")
        #mtext(side=1,at=c(2*x1[63]-x1[62],3*x1[63]-2*x1[62]),line=0.7,
        #    c("10% Soft linear","10% Hard linear"),col="grey50",las=2)


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

        ## white zone HF
#        points(x1[names(y1)=="Cult"]+0.15,y1[names(y1)=="Cult W=1"],pch="-",cex=2.1,col="red4")
#        lines(rep(x1[names(y1)=="Cult"]+0.15, 2), pr2[names(y1)=="Cult W=1", 3:4], col="red1")
#        points(x1[names(y1)=="UrbInd"]+0.15,y1[names(y1)=="UrbInd W=1"],pch="-",cex=2.1,col="red4")
#        lines(rep(x1[names(y1)=="UrbInd"]+0.15, 2), pr2[names(y1)=="UrbInd W=1", 3:4], col="red1")

        ## burn
        if (burn && attr(pr, "burn") > 0) {
            iii <- grep("  0", names(y1))
            points(x1[iii]-0.25,rep(y1["Burn"], length(iii)),pch="x",cex=1,col="black")
        }

        par(op)

    invisible(pr2)
}
