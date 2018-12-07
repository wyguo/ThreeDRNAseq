#' Estimate mean-variance trend of expression
#' @param obj a matrix of gene/transcript read counts, or a \code{\link[edgeR]{DGEList}} object. Read counts
#' must be non-negative and NAs are not permitted.
#' @param design design matrix with rows of samples and columns of conditions to fit a linear model.
#' @param lib.size a vector of library size for each sample. If \code{NULL}, the library size will be either extracted
#' from the \code{DGEList} object or calculated by summing up all the reads of genes/transcripts in each sample.
#' @param span a numeric value passed to \code{\link{lowess}} smoothing window as a proportion to fit the mean-variance trend.
#' @param ... additional arguments passed to \code{\link[limma]{lmFit}} function.
#' @return  a list object with elements:
#'     \itemize{
#'       \item{ fit: }{ the linear regression results of expression based on design matrix.}
#'       \item{ sx: }{ a numeric vector of mean values. }
#'       \item{ sy: }{ a numeric vector of variance values.}
#'       \item{ l: }{ the \code{lowess} fit result based on sx and sy.}
#'     }
#' @export
#' @seealso \code{\link[limma]{voom}}, \code{\link{lowess}} and \code{\link{condition2design}}.
#'
trend.mean.variance <- function(obj, design, lib.size = NULL, span = 0.5, ...){
  if (is(obj, "DGEList")) {
    counts <- obj$counts
    if (is.null(lib.size))
      lib.size <- with(obj$samples, lib.size * norm.factors)
  } else {
    counts <- as.matrix(obj)
    if(is.null(lib.size))
      lib.size <- colSums(counts)
  }
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean))
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  list(fit=fit,sx=sx,sy=sy,l=l)
}


#' Compare expression mean-variance trends before and after low expression filters.
#' @param counts.raw a matrix/data.frame of raw read counts before low expression filters.
#' @param counts.filtered a matrix/data.frame of read counts after low expression filters.
#' @param condition a vector of characters to distinguish conditions of samples (e.g. c('A','A','B','B')), which is used to make the design
#' matrix to fit linear regression model.
#' @param span a numeric value passed to \code{\link{lowess}} smoothing window as a proportion to fit the mean-variance trend.
#' @param ... additional arguments passed to \code{\link{trend.mean.variance}} function.
#' @return a list object with elements "fit.raw" and "fit.flitered", which are the \code{\link{trend.mean.variance}} fit results
#' for \code{counts.raw} and \code{counts.filtered}, respectively.
#' @export
#' @seealso \code{\link{trend.mean.variance}} and \code{\link[limma]{voom}}.
check.mean.variance <- function(counts.raw,
                                counts.filtered,
                                condition,
                                span = 0.5,...){
  ###---design matrix
  condition <- factor(condition,levels = unique(condition))
  design<-model.matrix(~0+condition)
  colnames(design) <- gsub('condition','',colnames(design))

  ###---generate DGEList object
  message('=> Generate DGEList object')
  dge.raw <- DGEList(counts=counts.raw)
  dge.raw <- calcNormFactors(dge.raw)
  dge.filtered<- DGEList(counts=counts.filtered)
  dge.filtered <- calcNormFactors(dge.filtered)

  ###---fit mean-variance trend
  message('=> Fit mean-variance trend')
  #before filter
  fit.raw <- trend.mean.variance(obj = dge.raw,design = design,span = span,...)
  #after filter
  fit.filtered <- trend.mean.variance(obj = dge.filtered,design = design,span = span,...)

  message('Done!!!')
  return(list(fit.raw=fit.raw,fit.filtered=fit.filtered))
}

#' Plot the mean-variance trend
#' @param x a numeric vector of mean value from the results of \code{\link{trend.mean.variance}}.
#' @param y a numeric vector of variance value from the results of \code{\link{trend.mean.variance}}.
#' @param x.lab,y.lab a string of x-axis label and y-axis label, respectively.
#' @param main plot title.
#' @param ... additional arguments passed to \code{\link{plot}}.
#' @return a plot
#' @export
plotMeanVariance <- function(x,y,l,fit.line.col='red',
                               x.lab = "log2( count size + 0.5 )",
                               y.lab = "Sqrt(standard deviation)",
                               main="Mean-variance trend",lwd=1.5,...){
  plot(x, y, pch = 16, cex = 0.25,xlab=x.lab,ylab=y.lab,main=main,...)
  lines(l,col=fit.line.col,lwd=lwd)
}

