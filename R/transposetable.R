#' Transpose a data table
#' 
#' Converts table rows to columns and columns to rows.
#' 
#' @param fnamein character vector that specifies the name of
#' tab-delimited text file with the input data table.
#' @param fnameout character vector that specifies the name of output file
#' where the transposed data table will be saved.
#' @param wdir character vector that specifies the name of the working
#' directory for the input/output files (defaults to the current R directory).
#' @return no return value
#'
#' @export
#' @examples
#'
#' # Example with data files included in the package:
#'
#' # Load example gene expression data table for 2 genes
#' fdat <- system.file("extdata", "example_genes.txt", package = "kmcut")
#'
#' transpose_table(fnamein = fdat,
#'                 fnameout = "example_genes_transposed.txt")
#'
#' # This will create in the current working directory a tab-delimited text
#' # file with the transposed table: "example_genes_transposed.txt"

transpose_table<-function(
    # The name of tab-delimited text file with the input data table
    fnamein,
    # The name of output file where the transposed data table will be saved
    fnameout,
    # Working directory for the input/output files
    wdir = getwd()
)
# begin function
{
  setwd(wdir)
  error <- character(0)
  
  if(file.access(fnamein, mode = 4) == -1)
  {
    error <- c(error, sprintf("Unable to open file <%s>\n", fnamein))
  }
  if(length(error) > 0) stop(error)
  
dat <- read.delim(fnamein, header = TRUE, row.names = 1,
            stringsAsFactors = FALSE)

df <- as.data.frame(t(dat))
colnames(df) <- row.names(dat)

# Convert row names into a column
df <- cbind(rownames(df), df)
colnames(df)[1] <- "tracking_id"
write.table(df, file = fnameout, quote = FALSE, row.names = FALSE,
        col.names = TRUE, sep = "\t")
}
# end function
