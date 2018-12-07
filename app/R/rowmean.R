#' Calculate column mean of a matrix or data frame based on a grouping variable
#'
#' Compute column means across rows of a numeric matrix-like object for each level
#' of a grouping variable.
#' @param x a matrix or data frame.
#' @param group a vector of factor giving grouping, with one element per row of x.
#' @param reorder if TRUE, then the result will be in order of \code{sort(unique(group))}.
#' @param na.rm logical (TRUE or FALSE). Should NA (including NaN) values be replaced by value 0?
#' @return \code{rowmean} returns a matrix or data frame containing the means. There
#' will be one row per unique value of group.
#'
#' @examples
#' x <- matrix(runif(50), ncol = 5)
#' group <- sample(1:4, 10, TRUE)
#' xmean <- rowmean(x, group)
#'
#' @export
#' @seealso \code{\link{rowsum}}, \code{\link{rowratio}}
rowmean<-function(x,group,reorder=F,na.rm=T){
  order.idx<-as.character(unique(group))
  if (reorder)
    order.idx<-gtools::mixedsort(order.idx)
  
  counts <- table(group)[order.idx]
  sums <- rowsum(x, group = group)[order.idx,]
  means <- (diag(1/counts)) %*% as.matrix(sums)
  rownames(means) <- order.idx
  if (na.rm)
    means[is.na(means)] <- 0
  return(means)
}



