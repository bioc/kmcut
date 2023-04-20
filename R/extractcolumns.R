#' Extract a sub-set of columns (such as a sub-set of samples)
#' from a data table. All rows will be preserved.
#' @param fnamein character vector that specifies the name of
#' tab-delimited text file with the input data table.
#' @param fids character vector that specifies the name of
#' text file with column ids (such as sample ids).
#' The file must contain one column id per line, without
#' any trailing spaces or any other additional symbols.
#' @param fnameout character vector that specifies the name of
#' output file where the new data table will be saved.
#' @param wdir character vector that specifies the name of
#' the working directory for the input/output files (defaults to the
#' current R directory).
#' @return no return value
#'
#' @export
#' @examples
#'
#' # Example with built-in data files:
#'
#' library(data.table)
#' library(kmcut)
#'
#' # Load example gene expression data table for 2 genes
#' fdat <- system.file("extdata", "example_genes_295.txt", package = "kmcut")
#' # Load a list that contains column (sample) ids
#' idlist <- system.file("extdata", "columnids.txt", package = "kmcut")
#' # Run the function
#' extractcolumns(fnamein = fdat, fids = idlist, 
#' fnameout = "example_samples_subset.txt")
#'
#' # This will create in the current working directory a tab-delimited text file
#' # "example_samples_subset.txt"

extractcolumns<-function(
    # The name of tab-delimited text file with the input data table
    fnamein=NULL,
    # The name of text file with column ids (such as sample ids)
    fids=NULL,
    # The name of output file where the new data table will be saved
    fnameout=NULL,
    # Working directory for the input/output files
    wdir=getwd()
)
# begin function
{
if(!is.null(wdir)) setwd(wdir)

dat <- read.delim(fnamein, header = TRUE, row.names = 1,
                    stringsAsFactors = FALSE)

ids <- unique( unlist(read.delim(fids, header = FALSE,
                    stringsAsFactors = FALSE)) )

snames <- colnames(dat)

df <- dat[, snames %in% ids]
n <- length(colnames(df))

# Convert to table and convert row names into a column
setDT(df, keep.rownames=TRUE)
colnames(df)[1] <- "tracking_id"
write.table(df, file = fnameout, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
}
# end function

