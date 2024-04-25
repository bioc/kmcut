#' Fit and validate Cox regression models
#' 
#' For each feature, fits a univariate Cox regression model on training data
#' and then uses the model to predict the risk score for test data.
#'
#' @param obj1 SummarizedExperiment object with training expression 
#  and survival data.
#' @param obj2 SummarizedExperiment object with test expression 
#  and survival data.
#' @param bfname a character string (character vector of length 1) that 
#' specifies the base name used to create the output file names, which 
#' are created by adding '_cox_train_sum.txt', '_train_scores.txt', and 
#' '_test_scores.txt' to 'bfname'.
#' @param wdir a character string (character vector of length 1) that
#' specifies the name of the working directory for the output files 
#' (defaults to the current R directory).
#' @param min_uval numeric value that specifies the minimal percentage of
#' unique values per feature (default is 50).
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
#' # Load training (fdat1) and test (fdat2) gene expression data
#' # files and survival data file (sdat).
#' fdat1 <- system.file("extdata", "expression_data_1.txt", package = "kmcut")
#' fdat2 <- system.file("extdata", "expression_data_2.txt", package = "kmcut")
#' sdat <- system.file("extdata", "survival_data.txt", package = "kmcut")
#'
#' # Create SummarizedExperiment object with training data
#' se1 <- create_se_object(efile = fdat1, sfile = sdat)
#' 
#' # Create SummarizedExperiment object with test data
#' se2 <- create_se_object(efile = fdat2, sfile = sdat)
#' 
#' # Fit Cox model on the training data and use it to calculate the risk
#' # scores for the test data.
#' ucox_pred(obj1 = se1, obj2 = se2, bfname = "demo", min_uval = 90)
#'
#' # This will create three output files in the current working directory:
#' # 1) Tab-delimited text file with Cox summary for the training data:
#' # "demo_cox_train_sum.txt"
#' # 2) Tab-delimited text file with the risk scores for training data:
#' # "demo_train_scores.txt"
#' # 3) Tab-delimited text file with the risk scores for test data:
#' # "demo_test_scores.txt"

ucox_pred<-function(
    # SummarizedExperiment object with training data
    obj1,
    # SummarizedExperiment object with test data
    obj2,
    # Base file name
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

# Training data
train_edat <- as.data.frame(assay(obj1))
# Test data
test_edat <- as.data.frame(assay(obj2))

# Name of the output file with COXPH summary for the training data
cox_out_train <- sprintf("%s_cox_train_sum.txt", bfname)
# Name of the output file with the risk scores for training data
risk_out_train <- sprintf("%s_train_scores.txt", bfname)
# Name of the output file with the risk scores for test data
risk_out_test <- sprintf("%s_test_scores.txt", bfname)

# The training survival time data
train_sdat_init <- get_sdat(obj1)
# Must have "status" and "time" variables
train_sdat <- data.frame(train_sdat_init$sample_id, train_sdat_init$scens,
                train_sdat_init$stime)
colnames(train_sdat) <- c("sample_id","status","time")

# The test survival time data
test_sdat_init <- get_sdat(obj2)
# Must have "status" and "time" variables
test_sdat <- data.frame(test_sdat_init$sample_id, test_sdat_init$scens,
                test_sdat_init$stime)
colnames(test_sdat) <- c("sample_id","status","time")

# Remove genes that have less than U% of unique values in the training dataset
train_edat <- filter_unique(train_edat, min_uval)
# Remove genes that have less than U% of unique values in the test dataset
test_edat <- filter_unique(test_edat, min_uval)

# Merge train expression data with survival data
train_edat <- as.data.frame(t(train_edat))

# Convert row names into first column
rnames <- rownames(train_edat)
train_edat <- cbind(rnames, train_edat)
colnames(train_edat)[1] <- "sample_id"

train_edat <- merge(train_edat, train_sdat, by.x = 1, by.y = 1, all = FALSE)
rownames(train_edat) <- unlist(train_edat["sample_id"])
train_edat["sample_id"] <- NULL

n_train_samples <- length(rownames(train_edat))

# Recalculate the number of features in the table
n_genes_train <- length(colnames(train_edat)) - 2

# Merge test expression data with survival data
test_edat <- as.data.frame(t(test_edat))
# Convert row names into first column
rnames <- rownames(test_edat)
test_edat <- cbind(rnames, test_edat)
colnames(test_edat)[1] <- "sample_id"

test_edat <- merge(test_edat, test_sdat, by.x = 1, by.y = 1, all = FALSE)
rownames(test_edat) <- unlist(test_edat["sample_id"])
test_edat["sample_id"] <- NULL

n_test_samples <- length(rownames(test_edat))

# Recalculate the number of features in the table
n_genes_test <- length(colnames(test_edat)) - 2
n_genes_test

# Table with coxph summary for the training data
cox_train <- matrix(data = NA, nrow = n_genes_train, ncol = 4)
colnames(cox_train) <- c("CC","HR","P","FDR_P")
rownames(cox_train) <- colnames(train_edat)[seq_len(n_genes_train)]

# Table with risk scores for training data
results_train <- matrix(data = NA, nrow = n_genes_train, ncol = n_train_samples)
colnames(results_train) <- rownames(train_edat)
rownames(results_train) <- colnames(train_edat)[seq_len(n_genes_train)]

# Table with risk scores for test data
results_test <- matrix(data = NA, nrow = n_genes_test, ncol = n_test_samples)
colnames(results_test) <- rownames(test_edat)
rownames(results_test) <- colnames(test_edat)[seq_len(n_genes_test)]

# Calculate the risk score for each feature in the training and test datasets
for(i in seq_len(n_genes_train))
{
    if(verbose == TRUE)
    {
        message("Processing ",i," of ",n_genes_train)
    }

    res <- coxph(as.formula(paste('Surv(time, status)~',
                colnames(train_edat)[i])),
                data=train_edat, model = TRUE)
    r <- summary(res)
    cox_train[i,"CC"] <- unlist(r$concordance["C"])
    cox_train[i,"HR"] <- r$coefficients[2]
    df <- as.data.frame(r["logtest"])
    cox_train[i,"P"] <- df["pvalue",1]

    pred_resubst <- predict(res, newdata = train_edat, type="risk")
    results_train[i,] <- pred_resubst

    crow <- which(colnames(train_edat) == colnames(test_edat)[i])
    if(length(crow) != 1) next

    pred_test <- predict(res, newdata = test_edat, type="risk")
    results_test[crow, ] <- pred_test
}

cox_train[,"FDR_P"]  <- p.adjust(cox_train[,"P"], method = "fdr")

# Convert row names into a column
df <- as.data.frame(cox_train)
df <- cbind(rownames(df), df)
colnames(df)[1] <- "tracking_id"
write.table(df, file = cox_out_train, quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

df <- as.data.frame(results_train)
df <- cbind(rownames(df), df)
colnames(df)[1] <- "tracking_id"
write.table(df, file = risk_out_train, quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

df <- as.data.frame(results_test)
df <- cbind(rownames(df), df)
colnames(df)[1] <- "tracking_id"
write.table(df, file = risk_out_test, quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

}
# end function
