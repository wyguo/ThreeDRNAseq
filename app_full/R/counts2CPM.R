#' Convert read counts to counts per million reads (CPM)
#' @param obj a read count matrix or a \code{\link{DGEList}} object.
#' @param lib.size a numeric vector of library size (sequencing depth) of each sample. If \code{NULL}, library
#' size of a sample is calculated from sum of all reads in this sample.
#' @param Log logcial, whether to take \eqn{\log_2} of CPM.
#' 
#' @details CPM is defined as:
#' \deqn{\frac{counts}{lib.size}\times 10^6}
#' To avoid taking \eqn{\log_2} of zero values, small offsets 0.5 and 1 were added to the numerator and 
#' denominator, respectively, i.e.
#' \deqn{\log_2 \frac{counts + 0.5}{lib.size+1}\times 10^6}
#' 
#' @return a numeric matrix of CPM or \eqn{\log_2}-CPM.
#' @export
#' @seealso \code{\link[edgeR]{cpm}}.
#' @examples 
#' y <- matrix(rnbinom(20,size=1,mu=10),5,4)
#' counts2CPM(y,Log = T)
#' counts2CPM(y,Log = F)

counts2CPM <- function (obj,lib.size = NULL,Log=F){
  if (is(obj, "DGEList")) {
    counts <- obj$counts
    if (is.null(lib.size)) 
      lib.size <- with(obj$samples, lib.size * norm.factors)
  } else {
    counts <- as.matrix(obj)
    if(is.null(lib.size))
      lib.size <- colSums(counts)
  }
  if(Log){
    t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  } else {
    t(t(counts)/(lib.size) * 1e+06)
  }
}