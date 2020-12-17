## need sf to read in gdb
library(sf)
## stuff for veg processing
od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)

## gdb file
f <- "S:/GC_eric/FromEric/To_Brandon/Report_layers_3by7.gdb"

## layer names and years for 3x7
l <- st_layers(f)
n <- l$name
n <- sort(n[grepl("Veg_id", n)])
y <- as.integer(sapply(strsplit(n, "_"), "[[", 5))


## read in the i'th layer
i <- 1
for (i in 2:15) {
    gc()
    cat("---------------------\nTaking note of the time:", as.character(Sys.time()))
    YR <- y[i]
    cat("\nReading gdb layer for year", YR, "...\n")
    p <- st_read(f, n[i])
    d <- as.data.frame(p)
    d$Shape <- NULL
    attr(d, "sf_column") <- NULL
    attr(d, "agr") <- NULL
    if (is.null(d$YEAR))
        d$YEAR <- d$year
    u <- d$FEATURE_TY %in% c("HARVEST-AREA", "CUTBLOCK") & is.na(d$YEAR)
    if (sum(u) > 0) {
        cat("\n >>> found",
            sum(u), "rows with unknown CC age amounting to",
            round(sum(d$Shape_Area[u]/10^6), 3), "km^2 -> set to R")
        d$YEAR[u] <- YR
    }
    cat(" OK\nMaking veg/hf/soil info ...")
    dd <- make_vegHF_wide_v6(d,
        col.label="ABMI_ID",
        col.year=YR,
        col.HFyear="YEAR",
        col.HABIT="Combined_ChgByCWCS",
        col.SOIL="Soil_Type_1",
        HF_fine=TRUE, wide=TRUE)
    cat(" OK\nFixing age 0 ...")
    dx <- nonDuplicated(d[d$NSRNAME != "",], ABMI_ID, TRUE)[rownames(dd[[1]]),]
    dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)
#    xx <- dd[[1]]
#    colnames(xx)[endsWith(colnames(xx), "0")]
    cat(" OK\nSaving results for year", YR, "...")
    save(dd, file=paste0("s:/AB_data_v2020/data/raw/veghf/3x7/", YR, ".RData"))
    cat(" OK\n\n")
}


library(mefa4)
ff <- paste0("s:/AB_data_v2020/data/raw/veghf/3x7/", y, ".RData")
L <- list()
for (i in seq_along(ff)) {
    load(ff[i])
    L[[as.character(y[i])]] <- as.matrix(dd$veg_current)
}
str(L)

A3x7 <- array(0, c(nrow(L[[2]]), ncol(L[[2]]), length(ff)))
dimnames(A3x7) <- list(1:1656, colnames(L[[2]]), y)
dim(A3x7)
for (i in seq_along(ff)) {
    A3x7[rownames(L[[i]]),colnames(L[[i]]),i] <- L[[i]]
}
save(A3x7, file="s:/AB_data_v2020/data/raw/veghf/3x7/VegHF3x7_1999_2017_age0fix.RData")


D <- dim(A3x7)
rs <- sapply(1:D[3], function(i) rowSums(A3x7[,,i])/(21*10^6))
colnames(rs) <- y
boxplot(rs, ylab="x 21 km^2")

round(t(apply(rs, 2, summary)),2)


## check areas
i <- 1
p <- st_read(f, n[i])

d <- as.data.frame(p)
d$Shape <- NULL
attr(d, "sf_column") <- NULL
attr(d, "agr") <- NULL
if (is.null(d$YEAR))
    d$YEAR <- d$year
if (sum(u) > 0)
    d$YEAR[u] <- YR

(sum(d$Shape_Area)/1656)/10^6
(sum(dd[[1]])/1656)/10^6

z <- sum_by(d$Shape_Area/10^6, d$ABMI_ID)
range(z[,1])
s <- rowSums(dd[[1]])[rownames(z)]/10^6

plot(z[,1], s)
abline(0,1)

s2 <- rowSums(L[[i]])[rownames(z)]/10^6
plot(s2, s)
abline(0,1)

lapply(1:15, function(i) setdiff(colnames(L[[i]]), colnames(L[[1]])))

sapply(2:15, function(i) colSums(L[[i]][,c("Decid0", "Mixedwood0",
    "Pine0", "Spruce0", "TreedBog0",
    "TreedFen0", "TreedSwamp0", "CutBlocks")])) / 10^6


## summarizing 2017 HFI in 1km grid

od <- setwd("~/repos/recurring/veghf")
source("00-setup.R")
setwd(od)


f <- "s:/GC_eric/FromEric/BDQT_2019_veg61RefCond2017_HFI2017/20190705_Veg61RefCond2017HFI2017_grid1sqkm.sqlite"


db <- dbConnect(RSQLite::SQLite(), f)

dbListTables(db)
dbListFields(db, "Veg61RefCond2017HF2017_grid1sqkm")

d <- dbGetQuery(db,
    "SELECT
      GRID_LABEL,
      Origin_Year,
      Pct_of_Larch,
      NSRNAME,
      Soil_Type_1,
      FEATURE_TY,
      Combined_ChgByCWCS,
      YEAR,
      SHAPE_Area
    FROM
      `Veg61RefCond2017HF2017_grid1sqkm`;")

dbDisconnect(db)

d$Shape_Area <- d$SHAPE_Area
d$SHAPE_Area <- NULL
dd <- make_vegHF_wide_v6(d,
        col.label="GRID_LABEL",
        col.year=2017,
        col.HFyear="YEAR",
        col.HABIT="Combined_ChgByCWCS",
        col.SOIL="Soil_Type_1",
        HF_fine=TRUE, wide=TRUE)
dx <- nonDuplicated(d[d$NSRNAME != "",], GRID_LABEL, TRUE)[rownames(dd[[1]]),]
dd <- fill_in_0ages_v6(dd, dx$NSRNAME, ages_list)

summary(rowSums(dd[[1]])/10^6)

save(dd, file="s:/AB_data_v2020/data/raw/veghf/3x7/VegHFw2w_2017_age0fix.RData")

