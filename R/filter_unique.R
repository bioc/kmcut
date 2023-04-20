# Removes genes that have less than min_uval% of unique values
#' @keywords internal
filter_unique<-function(dat, min_uval)
{
    n_genes <- length(rownames(dat))
    N <- length(colnames(dat))
    rv <- logical(length = n_genes)
    for(i in seq_len(n_genes))
    {
    un <- 100 * length(unique(unlist(dat[i, ]))) / N
    if(un < min_uval)
    {
        rv[i] <- FALSE
    }else
    {
        rv[i] <- TRUE
    }
    }
    dat <- dat[rv, ]
    return(dat)
}
