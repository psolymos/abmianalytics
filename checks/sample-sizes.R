library(mefa4)

e1 <- new.env()
e2 <- new.env()
#load("e:/peter/AB_data_v2018/data/analysis/site/veg-hf_SiteCenter_v6verified.Rdata", envir=e3)
load("e:/peter/AB_data_v2018/data/analysis/site/veg-hf_BAM-BBS-BU_v6verified.Rdata", envir=e1)
load("e:/peter/AB_data_v2018/data/analysis/site/veg-hf_CameraARU_v6verified.Rdata", envir=e2)

x <- rbind(e1$dd_150m[[1]], e2$dd_150m[[1]])

summary(rowSums(x))
table(rowSums(x) / (150^2*pi) < 1.1)

x <- x[rowSums(x) / (150^2*pi) < 1.1,]

x <- x / rowSums(x)
m <- find_max(x)
d <- data.frame(table(m$index))
colnames(d) <- c("LCC", "n_dom")
rownames(d) <- d$LCC
d$avg_dom <- NA
d$min_dom <- NA
d$max_dom <- NA
d$sd_dom <- NA
for (i in rownames(d)) {
    tmp <- m[m$index==i,]
    if (nrow(tmp)) {
        d[i,"avg_dom"] <- mean(tmp$value)
        d[i,"min_dom"] <- min(tmp$value)
        d[i,"max_dom"] <- max(tmp$value)
        d[i,"sd_dom"] <- sd(tmp$value)
    }
}
write.csv(d, row.names=FALSE, file="HF-sample-sizes.csv")

summary(x[,"SeismicLineNarrow"]+x[,"SeismicLineWide"])
