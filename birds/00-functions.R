#' SQL tables represent character not factor
make_char2fact <- function(x) {
    if (is.null(dim(x)))
        if (is.character(x))
            return(as.factor(x))
    for (i in seq_len(ncol(x)))
        if (is.character(x[,i]))
            x[,i] <- as.factor(x[,i])
        x
}
#' This turns ('1', '11') into (' 1', '11') [same char width]
fill_char <- function(x, fill=" ") {
    x <- as.character(x)
    n <- nchar(x)
    m <- max(n)
    sapply(seq_len(length(x)), function(i) {
        paste0(c(rep(fill, m-n[i]), x[i]), collapse="")
    })
}
#' cleans up AOU codes
normalize_sppcode <- function(x) {
    f <- function(x) {
        x <- toupper(x)
        x <- gsub(" ", "", x)
        x
    }
    if (is.factor(x)) {
        levels(x) <- f(levels(x))
        return(x)
    }
    if (is.character(x)) {
        return(f(x))
    }
    stop("x must be character or factor")
}
