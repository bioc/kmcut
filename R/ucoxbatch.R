#' Fit Cox regression models in batch mode
#' 
#' For each feature, fits a univariate Cox regression and performs the
#' likelihood ratio test.
#'
#' @param obj SummarizedExperiment object with expression-like data
#  and survival data.
#' @param bfname a character string (character vector of length 1) that 
#' specifies the base name used to create the output file name, which 
#' is created by adding '_ucoxbatch.txt' to 'bfname'.
#' @param wdir a character string (character vector of length 1) that
#' specifies the name of the working 
#' directory for the output file (defaults to the current R directory).
#' @param min_uval numeric value that specifies the minimal percentage
#' of unique values per feature (default is 50)
#' Features that have less than 'min_uval' percent unique values are
#' excluded from the analysis.
#' @param psort logical value whether to sort the output table by p-values
#' in increasing order (default is FALSE).
#' @param verbose logical value whether to print progress (default is TRUE).
#' @return no return value
#'
#' @export
#' @examples
#'
#' # Example with data files included in the package:
#'
#' # Load example gene expression data and survival data for 2 genes
#' # and 93 samples
#' fdat <- system.file("extdata", "example_genes.txt", package = "kmcut")
#' sdat <- system.file("extdata", "survival_data.txt", package = "kmcut")
#'
#' # Create SummarizedExperiment object
#' se <- create_se_object(efile = fdat, sfile = sdat)
#' 
#' ucox_batch(obj = se, bfname = "test")
#'
#' # This will create in the current working directory a tab-delimited text
#' # file with the results: "test_ucoxbatch.txt"

ucox_batch<-function(
    # SummarizedExperiment object with expression data and survival data
    obj,
    # Base name for the output files
    bfname,
    # Working directory for the input/output files
    wdir = getwd(),
    # Min percentage of unique values in ]0, 100] for each feature
    min_uval = 50,
    # Option to sort the output table by p-values in increasing order
    # (FALSE by default)
    psort = FALSE,
    # Print progress (TRUE by default)
    verbose = TRUE
)
# begin function
{
setwd(wdir)
error <- character(0)

if(min_uval <= 0 || min_uval > 100)
{
    error<-c(error, "min_uval must be in ]0, 100]\n")
}

if(length(error) > 0) stop(error)

bfname <- file_path_sans_ext(bfname)
if(length(bfname) == 0)
  stop("The base file name must be 1 or more characters long\n")

# Name of the output TXT file
txt_file <- sprintf("%s_ucoxpbatch.txt", bfname)

# The survival time data
sdat_init <- get_sdat(obj)

# Must have "status" and "time" variables
sdat <- data.frame(sdat_init$sample_id, sdat_init$scens, sdat_init$stime)
colnames(sdat) <- c("sample_id", "status", "time")

# The input data table
edat <- as.data.frame(assay(obj))
edat <- filter_unique(edat, min_uval)
edat <- as.data.frame(t(edat))

# Convert row names into first column
rnames <- rownames(edat)
edat <- cbind(rnames, edat)
colnames(edat)[1] <- "sample_id"

edat <- merge(edat, sdat, by.x = 1, by.y = 1, all = FALSE)
rownames(edat) <- unlist(edat["sample_id"])
edat["sample_id"] <- NULL

# Recalculate the number of features in the table
n_genes <- length(colnames(edat)) - 2

# Results table
results <- matrix(data = NA, nrow = n_genes, ncol = 4)
colnames(results) <- c("CC","HR","P","FDR_P")
rownames(results) <- colnames(edat)[seq_len(n_genes)]

# Go over each feature and calculate C-index
for(i in seq_len(n_genes))
{
    if(verbose == TRUE)
    {
        message("Processing ",i," of ",n_genes)
    }

    res <- coxph(Surv(time, status) ~ edat[,i], data=edat, model = FALSE)
    r <- summary(res)
    results[i,"HR"] <- r$coef[2]
    results[i,"CC"] <- unlist(r$concordance["C"])
    df <- as.data.frame(r["logtest"])
    results[i,"P"] <- df["pvalue",1]
}
results[,"FDR_P"]  <- p.adjust(results[,"P"], method = "fdr")

if(length(rownames(results)) > 1 && psort == TRUE)
{
    results <- results[order(results[,"P"]),]
}

# Convert row names into a column
df <- as.data.frame(results)
df <- cbind(rownames(df), df)
colnames(df)[1] <- "tracking_id"
write.table(df, file = txt_file, quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

}
# end function
