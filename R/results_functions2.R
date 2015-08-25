## predict coefficients
predX_age_cc <- 
function(ht=c("Deciduous", "Mixedwood", "WhiteSpruce", "Pine", "BlackSpruce", "LarchFen"), 
ages=c(0,10,20,40,60,80,100,120,140), CC=FALSE) 
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

    fCC2 <- if (ht %in% c("Deciduous", "Mixedwood"))
        approxfun(age, decid)(wtage) else approxfun(age, conif)(wtage)
    fCC2[is.na(fCC2)] <- 0
    if (!CC)
        fCC2 <- 0

#    X[,"fCC1"] <- fCC1
    X[,"fCC2"] <- fCC2

#    if (ht=="Deciduous") {
#    }
    if (ht=="Mixedwood") {
        X[,"hab1Mixedwood"] <- 1
        X[,"isMix:wtAge05"] <- sqrt(wtage)
        X[,"isMix:wtAge"] <- wtage
        X[,"isMix:wtAge2"] <- wtage^2
    }
    if (ht=="WhiteSpruce") {
        X[,"hab1WhiteSpruce"] <- 1
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
    if (ht %in% c("LarchFen", "BlackSpruce")) {
        if (ht=="LarchFen")
            X[,"hab1LarchFen"] <- 1
        if (ht=="BlackSpruce")
            X[,"hab1BlackSpruce"] <- 1
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

pred_hab <- 
function(res)
{
    if (is.character(res))
        res <- loadSPP(res)
    stage <- 4
    est <- getEst(res, stage, na.out=FALSE)#[2,,drop=FALSE]
    #est <- est[,colSums(abs(est))>0]
    X <- Xn[1:14,colnames(est)]
    X[,-1] <- 0
    diag(X) <- 1
    rownames(X) <- levels(xn$hab1)

#    X2 <- X[11:14,]
#    rownames(X2) <- paste(rownames(X2), "W=1")
#    X2["Shrubland W=1","isOpenWZ"] <- 1
#    X2["Grassland W=1","isOpenWZ"] <- 1
#    X2["Cult W=1","isCultWZ"] <- 1
#    X2["UrbInd W=1","isUrbIndWZ"] <- 1

    X3 <- X[1,,drop=FALSE]
    rownames(X3) <- "Burn"
    X3[1,"hab1bBurn"] <- 1

    X <- rbind(predX_age_cc("Deciduous", CC=FALSE),
        predX_age_cc("Mixedwood", CC=FALSE),
        predX_age_cc("WhiteSpruce", CC=FALSE),
        predX_age_cc("Pine", CC=FALSE),
        predX_age_cc("BlackSpruce", CC=FALSE),
        predX_age_cc("LarchFen", CC=FALSE),
        predX_age_cc("Deciduous", CC=TRUE),
        predX_age_cc("Mixedwood", CC=TRUE),
        predX_age_cc("WhiteSpruce", CC=TRUE),
        predX_age_cc("Pine", CC=TRUE),
        X[7:14,])
    X[,"hab1bMixedwood"] <- X[,"hab1Mixedwood"]
    X[,"hab1bWhiteSpruce"] <- X[,"hab1WhiteSpruce"]
    X[,"hab1bPine"] <- X[,"hab1Pine"]
    X[,"hab1bBlackSpruce"] <- X[,"hab1BlackSpruce"]
    X[,"hab1bLarchFen"] <- X[,"hab1LarchFen"]
    X[,"hab1bBog"] <- X[,"hab1Bog"]
    X[,"hab1bFen"] <- X[,"hab1Fen"]
    X[,"hab1bSwamp"] <- X[,"hab1Swamp"]
    X[,"hab1bMarsh"] <- X[,"hab1Marsh"]
    X[,"hab1bShrubland"] <- X[,"hab1Shrubland"]
    X[,"hab1bGrassland"] <- X[,"hab1Grassland"]
    X[,"hab1bCult"] <- X[,"hab1Cult"]
    X[,"hab1bUrbInd"] <- X[,"hab1UrbInd"]
    X <- rbind(X, X3)

    pr <- predStat(X, est, level, n=0, ci=TRUE, raw=FALSE)
    #pr0 <- prR1 <- prS1 <- colMeans(pr)
    #prR1["Median"] <- median(exp(log(pr0["Median"]) + est[,"ROAD01"])) # 0/1
    #prS1["Median"] <- median(exp(log(pr0["Median"]) + 0.1*est[,"SoftLin_PC"])) # 10%
    #pr <- rbind(pr, MeanNoRoad=pr0, MeanRoad=prR1)
    #pr <- rbind(pr, MeanNoSofLin=pr0, MeanSofLin=prS1)
    out <- pr[,c(1,2,5,6)]
    ## burn should not be shown when it is not selected
    attr(out, "burn") <- sum(abs(est[,"hab1bBurn"]))
    out
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



pred_hf <- 
function(est, Xn, n=200, fillin=0, remn=FALSE) 
{
    cc <- paste(c("THF", "Lin", "Nonlin", "Succ", "Alien", "Noncult", "Cult", 
        "THF2", "Nonlin2", "Succ2", "Alien2", "Noncult2"),
        "KM", sep="_")
    if (remn)
        cc <- c(cc, paste(c("Remn", "Remn2"), "KM", sep="_"))
    ehf <- est
    ehf[] <- 0
    ehf[,cc] <- est[,cc]

        pr <- list()
        Types <- c("Cult", "UrbInd", "HFor", "SoftLin", "HardLin")
        if (remn)
            Types <- c("Remn", Types)
        for (type in Types) {
            Range <- if (type %in% c("SoftLin", "HardLin"))
                c(0, 0.25) else c(0, 1)
            hf <- seq(Range[1], Range[2], len=n)
            Xhf <- Xn[1:n,] # matrix(0, n, length(cc))
            Xhf[] <- 0
            rownames(Xhf) <- NULL
            if (remn & type == "Remn") {
                Xhf[,"Remn_KM"] <- hf
                Xhf[,"Remn2_KM"] <- hf^2
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
                    Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
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
                    Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
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
                    Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
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
                    Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
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
                    Xhf[,"Remn_KM"] <- pmin(fillin, 1-hf)
                    Xhf[,"Remn2_KM"] <- Xhf[,"Remn_KM"]^2
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

fig_hf_noremn <- 
function(spp, fillin=0) 
{
    pr <- pred_hf(spp, fillin=fillin)
    prR <- pr[c("Remn")]
    pr <- pr[c("HFor", "Cult", "UrbInd", "SoftLin", "HardLin")]
    ymax <- max(sapply(pr, function(z) max(z[,-1])))
    mmax <- min(max(max(sapply(pr, function(z) max(z[,2]))), 5), 10)
    ymax <- min(mmax, ymax)
    
    ymaxR <- max(prR[[1]][,-1])

    op <- par(mfrow=c(1,2))

    plot(prR[[1]]$hf, prR[[1]]$Median,
        xlab="", ylim=c(0,ymaxR),
        cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",typ="n",bty="n",
        ylab="Relative abundance",col.lab="grey50",cex.lab=1.2)
    axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80")
    axis(side=2,at=seq(0,ymaxR,0.5),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80",las=2)
    axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60")
    axis(side=2,at=seq(0,ymaxR,0.5),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60",las=2)
    mtext(side=1,at=50,adj=0.5,"Surrounding suitable habitat (%)",line=2,col="grey50",cex=1.2)
    box(col="grey60")
    mtext(side=3,at=0,adj=0,as.character(taxa(mm)[spp,"English_Name"]),col="grey30")

    Colb <- c("#00FF0015", "#DD000010", "#D0D00025", "#80008015", "#0000DD10")
    Coll <- c(rgb(0.2,0.6,0.2), "red3", "orange3", "#A000A0", "blue3")
    polygon(c(prR[[1]]$hf, rev(prR[[1]]$hf)),
        c(prR[[1]][,"CL1q"], rev(prR[[1]][,"CL2q"])),
        col=Colb[1], border=NA)  # Transparent colours so they overlap
    lines(prR[[1]]$hf, prR[[1]][,"Median"], col=Coll[1],lwd=2)

    plot(pr[[1]]$hf, pr[[1]]$Median,
        xlab="", ylim=c(0,ymax),
        cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",typ="n",bty="n",
        ylab="Relative abundance",col.lab="grey50",cex.lab=1.2)
    axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80")
    axis(side=2,at=seq(0,ymax,0.5),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80",las=2)
    axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60")
    axis(side=2,at=seq(0,ymax,0.5),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60",las=2)
    mtext(side=1,at=50,adj=0.5,"Surrounding human footprint (%)",line=2,col="grey50",cex=1.2)
    box(col="grey60")
#    mtext(side=3,at=0,adj=0,paste(as.character(taxa(mm)[spp,"English_Name"]),
#        "- All province"),col="grey30")

    for (i in 1:length(pr)) {
        polygon(c(pr[[i]]$hf, rev(pr[[i]]$hf)),
            c(pr[[i]][,"CL1q"], rev(pr[[i]][,"CL2q"])),
            col=Colb[i], border=NA)  # Transparent colours so they overlap
    }
    for (i in 1:length(pr))
        lines(pr[[i]]$hf, pr[[i]][,"Median"], col=Coll[i],lwd=2)
    text(rep(0.1, 5), c(0.9, 0.75, 0.6, 0.45, 0.3), names(pr), pos=4, col=Coll)
    
    par(op)
    
    invisible(NULL)
}

fig_hf2 <- 
function(spp) 
{
    prz <- pred_hf(spp, fillin=0)
    prf <- pred_hf(spp, fillin=1)
    prR <- prz[c("Remn")]
    prz <- prz[c("HFor", "Cult", "UrbInd", "SoftLin", "HardLin")]
    prf <- prf[c("HFor", "Cult", "UrbInd", "SoftLin", "HardLin")]
    ymax <- max(sapply(c(prz,prf), function(z) max(z[,-1])))
    mmax <- min(max(max(sapply(c(prz,prf), function(z) max(z[,2]))), 5), 10)
    ymax <- min(mmax, ymax)
    
    ymaxR <- max(prR[[1]][,-1])

    op <- par(mfrow=c(1,3))

    plot(prR[[1]]$hf, prR[[1]]$Median,
        xlab="", ylim=c(0,ymaxR),
        cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",typ="n",bty="n",
        ylab="Relative abundance",col.lab="grey50",cex.lab=1.2)
    axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80")
    axis(side=2,at=seq(0,ymaxR,0.5),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80",las=2)
    axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60")
    axis(side=2,at=seq(0,ymaxR,0.5),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60",las=2)
    mtext(side=1,at=50,adj=0.5,"Surrounding suitable habitat (%)",line=2,col="grey50",cex=1.2)
    mtext(side=1,at=50,adj=0.5,"Surrounding footprint = 0%",line=4,col="grey50",cex=1)
    box(col="grey60")
    mtext(side=3,at=0,adj=0,as.character(taxa(mm)[spp,"English_Name"]),col="grey30")

    Colb <- c("#00FF0015", "#DD000010", "#D0D00025", "#80008015", "#0000DD10")
    Coll <- c(rgb(0.2,0.6,0.2), "red3", "orange3", "#A000A0", "blue3")

    polygon(c(prR[[1]]$hf, rev(prR[[1]]$hf)),
        c(prR[[1]][,"CL1q"], rev(prR[[1]][,"CL2q"])),
        col=Colb[1], border=NA)  # Transparent colours so they overlap
    lines(prR[[1]]$hf, prR[[1]][,"Median"], col=Coll[1],lwd=2)

    for (iii in 1:2) {
        
        if (iii == 1) {
            pr <- prz
            sub <- "Surrounding suitable habitat = 0%"
        } else {
            pr <- prf
            sub <- "Surrounding suitable habitat = 100% - HF"
        }
        
        plot(pr[[1]]$hf, pr[[1]]$Median,
            xlab="", ylim=c(0,ymax),
            cex.axis=1.3,xaxs="i",yaxs="i",xaxt="n",yaxt="n",typ="n",bty="n",
            ylab="Relative abundance",col.lab="grey50",cex.lab=1.2)
        axis(side=1,at=seq(0,100,25),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80")
        axis(side=2,at=seq(0,ymax,0.5),tck=1,cex.axis=1,col.axis="grey60",col.ticks="grey80",las=2)
        axis(side=1,at=seq(0,100,25),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60")
        axis(side=2,at=seq(0,ymax,0.5),tck=0.01,cex.axis=1,col.axis="grey60",col.ticks="grey60",las=2)
        mtext(side=1,at=50,adj=0.5,"Surrounding human footprint (%)",line=2,col="grey50",cex=1.2)
        mtext(side=1,at=50,adj=0.5,sub,line=4,col="grey50",cex=1)
        box(col="grey60")
    #    mtext(side=3,at=0,adj=0,paste(as.character(taxa(mm)[spp,"English_Name"]),
    #        "- All province"),col="grey30")

        for (i in 1:length(pr)) {
            polygon(c(pr[[i]]$hf, rev(pr[[i]]$hf)),
                c(pr[[i]][,"CL1q"], rev(pr[[i]][,"CL2q"])),
                col=Colb[i], border=NA)  # Transparent colours so they overlap
        }
        for (i in 1:length(pr))
            lines(pr[[i]]$hf, pr[[i]][,"Median"], col=Coll[i],lwd=2)
        text(rep(0.1, 5), c(0.9, 0.75, 0.6, 0.45, 0.3), names(pr), pos=4, col=Coll)
    }
    
    par(op)
    
    invisible(NULL)
}

pred_yr <- 
function(spp, Wzone=0, trend=FALSE) 
{
    res <- loadSPP(spp)
    est0 <- getEst(res, length(mods), na.out=FALSE)
    est <- est0
    est[] <- 0
    est[,grepl("YR", colnames(est))] <- est0[,grepl("YR", colnames(est))]

    yr <- 1997:2014
    Xyr <- Xn[1:length(yr),]
    Xyr[] <- 0
    rownames(Xyr) <- NULL
    Xyr[,"YR"] <- yr - 1997
    Xyr[,"YR51"] <- ifelse(yr >= 2003 & yr < 2009, 1, 0)
    Xyr[,"YR52"] <- ifelse(yr >= 2009, 1, 0)
    if (trend) {
        pr <- exp(predStat(Xyr, est, level, raw=TRUE))
        yr0 <- 1997:2014 - 1997
        ## this is to fit a line to YR5 etc
        d <- numeric(ncol(pr))
        for (i in 1:ncol(pr)) {
            z <- pr[,i]
            tmp1 <- coef(lm(log(z+0.5) ~ yr0))
            tmp2 <- try(coef(nls(z ~ exp(a + b*yr0), start=list(a=tmp1[1], b=tmp1[2]))))
            d[i] <- if (inherits(tmp2, "try-error"))
                tmp1[2] else tmp2[2]
        }
        #d <- apply(pr, 2, function(z) coef(lm(log(z+1) ~ yr0))[2])
    } else {
        d <- data.frame(year=yr, predStat(Xyr, est, level, n=0, ci=TRUE, raw=FALSE)[,c("Median", "Mean", "CL1q","CL2q")])
    }
    d
}
