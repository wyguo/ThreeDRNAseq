#' Convert list to data frame
#' @param x a list object with elements of vectors. The lengths of these vectors can be different.
#' @return a data frame.
list2dataframe <- function(x){
  n <- max(sapply(x, length))
  y <- lapply(x,function(i){
    c(i,rep(NA,n-length(i)))
  })
  z <- do.call(cbind,y)
  data.frame(z,check.names = F)
}