#' Extract a sub-set of columns
#' 
#' Extract a sub-set of columns (such as a sub-set of samples)
#' from a data table. All rows will be preserved.
#' 
#' @param fnamein a character string (character vector of length 1) that
#' specifies the name of tab-delimited text file with the input data table.
#' @param fids a character string (character vector of length 1) that
#' specifies the name of text file with column ids (such as sample ids).
#' The file must contain one column id per line, without
#' any trailing spaces or any other additional symbols.
#' @param fnameout a character string (character vector of length 1) that
#' specifies the name of output file where the new data table will be saved.
#' @param wdir a character string (character vector of length 1) that
#' specifies the name of the working directory for the input/output files 
#' (defaults to the current R directory).
#' 
#' @return no return value
#'
#' @export
#' 
#' @examples
#'
#' # Example with built-in data files:
#'
#' # Load example gene expression data table for 2 genes
#' fdat <- system.file("extdata", "example_genes.txt", package = "kmcut")
#' # Load a list that contains column (sample) ids
#' idlist <- system.file("extdata", "columnids.txt", package = "kmcut")
#' # Run the function
#' extract_columns(fnamein = fdat, fids = idlist, 
#'                 fnameout = "example_samples_subset.txt")
#'
#' # This will create in the current working directory a tab-delimited text file
#' # "example_samples_subset.txt"

extract_columns<-function(
    # The name of tab-delimited text file with the input data table
    fnamein,
    # The name of text file with column ids (such as sample ids)
    fids,
    # The name of output file where the new data table will be saved
    fnameout,
    # Working directory for the input/output files
    wdir=getwd()
)
# begin function
{
  setwd(wdir)
  error <- character(0)
  
  if(file.access(fnamein, mode = 4) == -1)
  {
    error <- c(error, sprintf("Unable to open file <%s>\n", fnamein))
  }
  if(file.access(fids, mode = 4) == -1)
  {
    error <- c(error, sprintf("Unable to open file <%s>\n", fids))
  }
  if(length(error) > 0) stop(error)

dat <- read.delim(fnamein, header = TRUE, row.names = 1,
                    stringsAsFactors = FALSE)

ids <- unique( unlist(read.delim(fids, header = FALSE,
                    stringsAsFactors = FALSE)) )
if(length(ids) == 0) stop("Input file with column IDs is blank\n")

snames <- colnames(dat)

df <- dat[, snames %in% ids]
if(ncol(df) == 0)
  stop("No matching sample IDs were found in the input data table\n")

# Convert row names into a column
df <- cbind(rownames(df), df)
colnames(df)[1] <- "tracking_id"
write.table(df, file = fnameout, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
}
# end function

