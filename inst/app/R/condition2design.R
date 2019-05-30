#' Convert conditions to design matrix for linear regression
#' @param condtion a vector of conditions of samples corresponding to the expression dataset.
#' @param batch.effect batch effect terms estimated from \code{remove.batch}. Default is \code{NULL}.
#' @return A design matrix for linear regression.
#' @export
condition2design <- function(condition,batch.effect=NULL){
  design.data <- data.frame(condition=condition)
  if(!is.null(batch.effect)){
    colnames(batch.effect) <- paste0('batch',1:ncol(batch.effect))
    design.data <- data.frame(design.data,batch.effect)
  }
  design.fomula <- as.formula(paste0('~0+',paste0(colnames(design.data),collapse = '+')))
  design <- model.matrix(design.fomula, data=design.data)
  idx <- grep('condition',colnames(design))
  design.idx <- colnames(design)
  design.idx <- substr(design.idx,nchar('condition')+1,nchar(design.idx))
  colnames(design)[idx] <- design.idx[idx]
  design
}