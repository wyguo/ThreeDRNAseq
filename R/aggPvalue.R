#' Aggregate multiple p-values to a single p-value
#' @param pvals a vector of p-values.
#' @param method either "Fisher" or "Simes" method can be used to aggregate p-values.
#' @param returnstat logical. If "TRUE" then return the statistic as well as the p-values.
#' @export
#' 
aggPvalue <- function(pvals,method=c('Fisher','Simes'), returnstat = FALSE){
  method <- match.arg(arg = method,c('Fisher','Simes'))
  pval <- switch(EXPR = method,
                 Simes={
                   r = rank(pvals)
                   x = min(length(pvals) * pvals/r)
                   if (returnstat) 
                     c(x, x)
                   else x
                 },
                 Fisher={
                   x <- -2*sum(log(pvals))
                   p = 1 - pchisq(x, df = 2 * length(pvals))
                   if (returnstat) 
                     c(p, x)
                   else p
                 })
  return(pval)
}