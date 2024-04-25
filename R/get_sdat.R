#' @keywords internal

get_sdat<-function(obj)
{
  sdat <- as.data.frame(colData(obj))
  sdat <- cbind(rownames(colData(obj)), sdat)
  colnames(sdat)[1] <- "sample_id"
  return(sdat)
}
