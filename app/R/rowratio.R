#' Calculate column raito of a matrix or data frame based on a grouping variable
#'
#' Compute column ratio across rows of a numeric matrix-like object for each level
#' of a grouping variable.
#' @param x a matrix or data frame.
#' @param group a vector of factor giving grouping, with one element per row of x.
#' @param reorder if TRUE, then the result will be in order of row names of x. If row names of x is null, the
#' results is not reordered.
#' @param na.rm logical (TRUE or FALSE). Should NA (including NaN) values be replaced by value 0?
#' @return \code{rowratio} returns a matrix or data frame containing the ratios. There
#' will be one row per unique value of group.
#'
#' @examples
#' x <- matrix(runif(50), ncol = 5)
#' group <- sample(1:4, 10, TRUE)
#' xratio <- rowratio(x, group)
#'
#' @export
#'
#' @seealso \code{\link{rowsum}}, \code{\link{rowmean}}
#' 
rowratio <- function (x, group, reorder = F, na.rm = T){
  y <- rowsum(x, group = group)
  y <- y[group, ]
  ratio = x/y
  rownames(ratio) <- rownames(x)
  if (na.rm) 
    ratio[is.na(ratio)] <- 0
  if (reorder & !is.null(rownames(ratio))) 
    ratio <- ratio[gtools::mixedsort(rownames(ratio)), ]
  return(ratio)
}