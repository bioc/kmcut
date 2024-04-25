#' Extract a sub-set of rows
#' 
#' Extract a sub-set of rows (such as a group of gene ids)
#' from a data table. All columns will be preserved.
#' 
#' @param fnamein a character string (character vector of length 1) that
#' specifies the name of tab-delimited text file with the input data table.
#' @param fids a character string (character vector of length 1) that
#' specifies the name of text file with row ids (such as gene ids).
#' The file must contain one row id per line, without any trailing spaces or
#' any other additional symbols.
#' @param fnameout a character string (character vector of length 1) that
#' specifies the name of output file where the new data table will be saved.
#' @param wdir a character string (character vector of length 1) that
#' specifies the name of the working 
#' directory for the input/output files (defaults to the current R directory).
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
#' # Load a list that contains one gene id (MYCN)
#' idlist <- system.file("extdata", "rowids.txt", package = "kmcut")
#' # Run the function
#' extract_rows(fnamein = fdat, fids = idlist,
#'              fnameout = "example_genes_subset.txt")
#'
#' # This will create in the current working directory a tab-delimited text file
#' # "example_genes_subset.txt" with one row "MYCN".

extract_rows<-function(
        # The name of tab-delimited text file with the input data table
    fnamein,
    # The name of text file with row ids (such as gene ids)
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
  
dat <- read.delim(fnamein, header = TRUE, stringsAsFactors = FALSE)
dat[,1] <- make.names(dat[,1], unique = TRUE)

# The gene list
ids <- unique(as.vector(unlist(read.delim(fids, header = FALSE))))

if(length(ids) == 0) stop("Input file with gene IDs is blank\n")

dat <- dat[dat[,1] %in% ids,]

if(nrow(dat) == 0)
  stop("No matching gene IDs were found in the input data table\n")

write.table(dat, file = fnameout, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
}
# end function

