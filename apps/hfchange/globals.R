#devtools::install_github("ramnathv/rCharts")

#DIR <- "e:/peter/sppweb2017/apps/hfchange"

#library(mefa4)
#library(rgdal)
#library(rgeos)
library(sp)
#library(gstat)
library(raster)
library(viridis)
#library(ggplot2)
#library(rCharts)
library(shiny)
library(googleVis)
library(mapview)

load("hfchange.rda")

get_data0 <- function(r, c, byregion=TRUE, which_region="nr") {
    if (length(c) < 1)
        stop("Select at least one industrial sector.")
    YR <- as.integer(dimnames(HF)[[3]])
    n <- table(GIS[[which_region]])
    NAM <- levels(GIS[[which_region]])
    if (byregion) {
        hf <- lapply(dimnames(HF)[[3]], function(z)
            groupMeans(HF[,,z], 1, GIS[[which_region]])[NAM[r],c,drop=FALSE])
        out <- matrix(NA, length(YR), length(r))
        rownames(out) <- dimnames(HF)[[3]]
        colnames(out) <- dimnames(hf[[1]])[[1]]
        for (i in 1:length(YR))
            out[i,] <- rowSums(hf[[i]])
    } else {
        hf <- lapply(dimnames(HF)[[3]], function(z)
            groupSums(HF[,,z], 1, GIS[[which_region]])[NAM[r],c,drop=FALSE])
        out <- matrix(NA, length(YR), length(c))
        rownames(out) <- dimnames(HF)[[3]]
        colnames(out) <- dimnames(hf[[1]])[[2]]
        for (i in 1:length(YR))
            out[i,] <- colSums(hf[[i]]) / sum(n[r])
    }
    list(x=YR, y=out)
}

get_mape <- function(r, c, which_region="nr") {

    if (length(c) < 1)
        stop("Select at least one industrial sector.")
    nr <- seq_len(nlevels(GIS[[which_region]]))
    names(nr) <- dimnames(HF[[which_region]])[[1]]
    Show <- nr %in% r

    Footprint <- RHF[[c[1]]]
    if (length(c) > 1)
        for (i in 2:length(c))
            Footprint <- RHF[[c[i]]] + Footprint
    Footprint[!(REG[[which_region]] %in% nr[r])] <- NA

    m <- mapview(Footprint, legend=TRUE,
        col.regions=rev(viridis::magma(256)),
        trim = FALSE,
        alpha.regions=0.7, legend.opacity = 0.7)
    m@map
}
get_gplot <- function(r, c, byregion=TRUE, which_region="nr") {
    d <- get_data0(r, c, byregion, which_region)
    gvisLineChart(data.frame(Year=d$x, d$y),
        options=list(gvis.editor="Edit", width="100%", height="450",
            hAxis="{title:'Year'}", vAxis="{title:'Percent'}"))
}

if (FALSE) {

get_data <- function(r, c, byregion=TRUE) {
    d <- get_data0(r, c, byregion)
    out <- data.frame(year=d$x, stack(data.frame(round(d$y, 2))))
    #WHAT <- if (byregion)
    #    "Region" else "Sector"
    colnames(out)[1:2] <- c("Year", "Percent")
    #out$int <- out[[WHAT]]
    out$yr <- as.numeric(as.POSIXct(paste0(d$x, "-01-01")))
    out
}

get_map <- function(r) {
    nr <- c(4, 5, 3, 1, 2, 6)
    names(nr) <- c("Grassland", "Parkland", "Foothills", "Boreal",
        "Canadian Shield", "Rocky Mountain")
    Col <- ifelse(nr %in% r, "#0000FF80", "#00000020")
    Col2 <- ifelse(pts@data$NATURAL_REGIONS %in% names(nr)[r], "#0000FFFF", "#00000060")
    op <- par(mar=c(1,1,1,1))
    on.exit(par(op))
    plot(ABnr, col=Col, border=Col)
    plot(pts, pch=".", col=Col2, add=TRUE)
    invisible(NULL)
}

get_plot <- function(r, c, byregion=TRUE, ...) {
    d <- get_data(r, c, byregion)
    pl <- ggplot(d, aes(x=Year, y=Percent, colour=ind)) +
        geom_line(aes(group=ind))
    pl
}

get_hplot <- function(r, c, byregion=TRUE, ...) {
    d <- get_data(r, c, byregion)
    pl <- hPlot(Percent ~ Year, group = "ind", data = d, type = "line",
        radius=5, title="Human Footprint",
        subtitle=if (byregion) "By Region" else "By Sector")
    pl$addParams(dom = 'myChart')
    pl
}

get_rplot <- function(r, c, byregion=TRUE, ...) {
    d <- get_data(r, c, byregion)
    pl <- Rickshaw$new()
    pl$layer(Percent ~ yr, group = "ind", data = d, type = "area", width = 560,
        title="Human Footprint",
        subtitle=if (byregion) "By Region" else "By Sector")
    # add a helpful slider this easily; other features TRUE as a default
    pl$set(slider = TRUE)
    pl$addParams(dom = 'myChart')
    pl
}

get_lmap <- function(r) {
    nr <- c(4, 5, 3, 1, 2, 6)
    names(nr) <- c("Grassland", "Parkland", "Foothills", "Boreal",
        "Canadian Shield", "Rocky Mountain")
    Show <- nr %in% r
    Show2 <- gis$NATURAL_REGIONS %in% names(nr)[r]

    leaflet() %>%
        addTiles() %>%
        addFeatures(as(ab[!Show], "sf"), color="#0000FF", opacity=0.2, stroke=TRUE) %>%
        addFeatures(as(ab[Show], "sf"), color="#FF0000", opacity=0.5, stroke=TRUE) %>%
        addFeatures(as(spTransform(pts[!Show2,], proj4string(ab)), "sf"), radius=2,
            color="#0000FF", opacity=0.4, stroke=FALSE) %>%
        addFeatures(as(spTransform(pts[Show2,], proj4string(ab)), "sf"), radius=2,
            color="#FF0000", opacity=0.8, stroke=FALSE)
}


get_rmap <- function(r, c) {

    if (length(c) < 1)
        stop("Select at least one industrial sector.")
    nr <- c(4, 5, 3, 1, 2, 6)
    names(nr) <- c("Grassland", "Parkland", "Foothills", "Boreal",
        "Canadian Shield", "Rocky Mountain")
    Show <- nr %in% r

    rrr <- rr[[c[1]]]
    if (length(c) > 1)
        for (i in 2:length(c))
            rrr <- rr[[c[i]]] + rrr
    msk <- rr[[1]]
    OK <- rnr %in% nr[r]
    msk[OK] <- NA
    rrr[!OK] <- NA
    op <- par(mar=c(1,1,1,1))
    on.exit(par(op))
    plot(rrr, col=rev(viridis::magma(100)), axes=FALSE, box=FALSE)
    plot(msk, col="lightgrey", add=TRUE, legend=FALSE)
    #plot(ABnr[!Show], add=TRUE, col="#808080", border=NA)
    plot(ABnr, add=TRUE, border="darkgrey")

    invisible(NULL)
}

}

