library(jsonlite)
library(whisker)

tmp <- readLines("pipeline/2020/template.txt")

X <- fromJSON("s:/reports_staging/2020/images/index.json")
TAXA <- X$taxa$id
fun <- function(x) {
    for (i in names(x))
        if (is.logical(x[[i]]))
            x[[i]] <- tolower(as.character(x[[i]]))
    x
}
OUT <- "~/tmp/nuxt/pages/%s/%s/"


for (taxon in TAXA) {
    SPP <- X$species$id[X$species$taxonid == taxon]
    for (spp in SPP) {
        cat(taxon, spp, "\n")
        flush.console()
        x <- fromJSON(
            sprintf(
                "s:/reports_staging/2020/images/%s/%s/index.json",
                taxon, spp)
        )
        w <- whisker.render(tmp, fun(x))
        DIR <- sprintf(OUT, taxon, spp)
        if (!dir.exists(DIR))
            dir.create(DIR)
        writeLines(w, paste0(DIR, "index.vue"))
    }
}
