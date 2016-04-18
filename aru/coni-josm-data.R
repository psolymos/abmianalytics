library(mefa4)
library(pbapply)

ROOT <- "e:/peter/AB_data_v2016"

e <- new.env()
load("e:/peter/AB_data_v2016/out/birds/data/data-full-withrevisit.Rdata", envir=e)

load("e:/peter/AB_data_v2016/out/birds/data/data-full-north-josm.Rdata")

YY <- cBind(YY, CONI=e$YY[rownames(YY),"CONI"])
table(YY[,"CONI"])

#> exp(coni_rem_2[["0"]]$coef)
#[1] 0.5940393
#> exp(coni_dis_2[["0"]]$coef)
#[1] 1.139891
stopifnot(all(rownames(DAT)==rownames(OFF)))
OFF <- cbind(OFF, CONI=log(ifelse(is.finite(DAT$MAXDIS), 
    DAT$MAXDIS^2*pi * (1-exp(-0.5940393*DAT$MAXDUR)) * 
        (1.139891^2 * (1 - exp(-DAT$MAXDIS^2/1.139891^2))/DAT$MAXDIS^2),
    1.139891^2*pi * (1-exp(-0.5940393*DAT$MAXDUR)))))
OFF <- OFF[,"CONI",drop=FALSE]
summary(OFF[,"CONI"])
SPP <- "CONI"

save(list=c("BB", "DAT", "HSH", "mods", "OFF", "ROOT", "SPP", "TAX", "YY"),
    file="e:/peter/AB_data_v2016/out/birds/data/data-full-north-josm-coni.Rdata")

