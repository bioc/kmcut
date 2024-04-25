#' Create SummarizedExperiment object
#' 
#' Reads a file with expression data and a file with survival data. Then, uses 
#' the data to create a SummarizedExperiment object. 
#'
#' @param efile a character string (character vector of length 1) that
#' specifies the name of the file with expression data 
#' for each sample. The file must be tab-delimited,
#' where genes are in rows and samples are in columns. First column must
#' contain gene names. Column names must contain sample ids.
#' @param sfile a character string (character vector of length 1) that
#' specifies the name of the file with 
#' right-censored survival time data. The file must be tab-delimited,
#' where samples are in rows. First column must contain sample ids that match
#' those in 'efile'. The file must contain columns called 'stime' and 'scens',
#' with survival time and censoring variable (0 or 1), respectively.
#' @param wdir a character string (character vector of length 1) that 
#' specifies the name of the working directory for the input files 
#' (defaults to the current R directory).
#' 
#' @return a SummarizedExperiment object that contains expression matrix along
#' with survival data as column data.
#'
#' @export
#' 
#' @examples
#'
#' # Example with data files included in the package:
#'
#' # Load example gene expression data and survival data for 2 genes
#' # and 93 samples:
#' fdat <- system.file("extdata", "example_genes.txt", package = "kmcut")
#' sdat <- system.file("extdata", "survival_data.txt", package = "kmcut")
#' 
#' # Create SummarizedExperiment object
#' se <- create_se_object(efile = fdat, sfile = sdat)

create_se_object<-function(
  # The file with  expression data
  efile,
  # The file with survival data
  sfile,
  # Working directory
  wdir = getwd()
)
{
  setwd(wdir)
  error <- character(0)
  
  if(file.access(efile, mode = 4) == -1)
  {
    error <- c(error, sprintf("Unable to open file <%s>\n", efile))
  }
  if(file.access(sfile, mode = 4) == -1)
  {
    error <- c(error, sprintf("Unable to open file <%s>\n", sfile))
  }
  if(length(error) > 0) stop(error)
  
  edat <- read.delim(efile, header = TRUE, row.names = 1,
                     stringsAsFactors = FALSE)
  sdat <- read.delim(sfile, header = TRUE, stringsAsFactors = FALSE)
  
  if("sample_id" %in% colnames(sdat) == FALSE) 
    error <- c(error, "Column 'sample_id' is missing in survival data file\n")
  if("stime" %in% colnames(sdat) == FALSE) 
    error <- c(error, "Column 'stime' is missing in survival data file\n")
  if("scens" %in% colnames(sdat) == FALSE) 
    error <- c(error, "Column 'scens' is missing in survival data file\n")
  if(length(error) > 0) stop(error)

  row.has.na <- apply(edat, 1,
                      function(x){any(is.na(x) | is.nan(x) | is.infinite(x))})
  s <- sum(row.has.na)
  if(s > 0)
  {
    error <- c(error,
               sprintf("File <%s> has missing or non-numeric values\n", efile))
  }
  col.has.na <- apply(sdat[, c("stime","scens")], 2,
                      function(x){any(is.na(x) | is.nan(x) | is.infinite(x))})
  s <- sum(col.has.na)
  if(s > 0)
  {
    error <- c(error,
               sprintf("File <%s> has missing or non-numeric values\n", sfile))
  }
  if(length(error) > 0) stop(error)
  
  if(all(sdat["scens"] == 0 | sdat["scens"] == 1) == FALSE) 
     error <- c(error,"The 'scens' column contains values other than 0 or 1\n")
  if(any(sdat["stime"] < 0) == TRUE) 
    error <- c(error,"The 'stime' column contains negative values\n")
  if(length(error) > 0) stop(error)
  
  sdat <- sdat[c("sample_id", "stime", "scens")]
  
  colnames(edat) <- make.names(colnames(edat), unique = TRUE)
  sdat[,"sample_id"] <- make.names(sdat[,"sample_id"], unique = TRUE)
 
  # Merge edat and sdat by samples
  edat <- as.data.frame(t(edat))
  # Convert row names into first column
  rnames <- rownames(edat)
  edat <- cbind(rnames, edat)
  colnames(edat)[1] <- "sample_id"
  edat <- merge(edat, sdat, by.x = 1, by.y = 1, all = FALSE)
  # Check for failed merge
  if(dim(edat)[1] <= 1 || dim(edat)[2] <= 1 || is.null(dim(edat)) == TRUE) 
    stop("No common sample IDs found in expression and survival data files")
  rownames(edat) <- unlist(edat["sample_id"])
  edat["sample_id"] <- NULL
  sdat <- edat[ , c("stime","scens")]
  edat <- as.data.frame(t(edat))
  edat <- edat[-c(nrow(edat)-1, nrow(edat)) ,]
  assay <- as.matrix(edat)
  coldata <- DataFrame(sdat)
  obj <- SummarizedExperiment(list(expr=assay), colData=coldata)
  return(obj)
}
