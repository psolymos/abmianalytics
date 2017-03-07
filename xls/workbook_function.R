#
# Wrapper to export R objects (dataframes, matrix, etc...) as a sheets
# within a workbook
# Created January 16, 2017

save.xlsx <- function (file, ...)

        {
        require(xlsx, quietly = TRUE)
        objects <- list(...)
        fargs <- as.list(match.call(expand.dots = TRUE))
        objnames <- as.character(fargs)[-c(1, 2)]
        nobjects <- length(objects)
        for (i in 1:nobjects) {
                if (i == 1)
                        write.xlsx(objects[[i]], file, sheetName = objnames[i],
                                   row.names = FALSE, showNA = FALSE)
                else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                                append = TRUE, row.names = FALSE, showNA = FALSE)
        }
        print(paste("Workbook", file, "has", nobjects, "worksheets."))
}