#' Estimate batch effects between biological replicates
#' @param read.counts data.frame/matrix of read counts.
#' @param condition condition labels of the columns in read count data, e.g. c('A','A','A','B','B','B','C','C','C').
#' @param design design matrix relating to conditions to be preserved. If \code{NULL}, design matrix will beg
#' generated from provided \code{condition}.
#' @param contrast a vector of contrast groups for condition comparisons, e.g. c('B-A','C-A').
#' @param group a vector or factor giving the experimental group/condition for each sample/library. See \code{\link[edgeR]{DGEList}} for details.
#' @param method one of "RUVr", "RUVg" and "RUVs". See the RUVSeq R package for details \url{http://bioconductor.org/packages/release/bioc/html/RUVSeq.html}.
#' @param cIdx a vector of character, logical or numeric to indicate the subset of genes to 
#' be used as negative controls in the estimation of the factors of unwanted variation. Default is \code{NULL}.
#' @param k a integer number of factors of unwanted variation to be estimated from the data.
#' @param ... additional arguments passed to \code{\link[edgeR]{DGEList}}.
#' @details The provided arguments \code{condition}, \code{design}, \code{contrast} and \code{group} are used
#' to fit a GLM model in \code{\link{edgeR}}, which generates residuals for the \code{\link[RUVSeq]{RUVr}} method.
#' If the negative controls (not differential expressed genes/transcripts) \code{cIdx=NULL}, the targets with p-value > 0.1
#' from \code{\link[edgeR]{glmQLFit}} will be used as negative controls. See the packages \code{\link{edgeR}} and \code{\link{RUVSeq}}
#' for detials.
#' @return A list object with elements:
#'     \itemize{
#'       \item{\code{W:} the biological-replicates-by-factors matrix with the estiamted factors of batch effects, 
#'       which can be passed to the design matrix of linear regression as batch effect terms.}
#'       \item{\code{normalizedCounts:} a expression matrix the same dimension as the input read counts, in which
#'       the batch effects have been removed according to the estimated factors. }
#'       \item{\code{method:} the method used to estimate the batch effects (RUVr, RUVg or RUVs).}
#'     }
#' @export
#' 
#' @examples 
##------> sample information
#' samples <- exp.data$samples
#' 
#' ##------> sum sequencing replicates
#' trans_counts <- exp.data$trans_counts
#' idx <- paste0(samples$condition,'_',samples$brep)
#' y <- sumarrays(trans_counts,group = idx)
#' 
#' ##------> update sample information after sum 
#' samples <- samples[samples$srep==samples$srep[1],]
#' 
#' ##------> normalisation
#' dge <- DGEList(counts=y)
#' dge <- calcNormFactors(dge)
#' data2pca <- t(counts2CPM(obj = dge,Log = T))
#' 
#' ##------> plot PCA
#' groups <- samples$brep
#' plot.PCA.ind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = groups,ellipse.type='polygon')
#' 
#' ##------> remove batch effects
#' y.new <- remove.batch(read.counts = y,
#'                       condition = samples$condition,
#'                       method = 'RUVr')
#' 
#' ##------> normalisation and plot PCA again
#' dge <- DGEList(counts=genes_batch$normalizedCounts)
#' dge <- suppressWarnings(calcNormFactors(dge))
#' data2pca <- t(counts2CPM(obj = dge,Log = T))
#' plot.PCA.ind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = groups,ellipse.type='polygon')

remove.batch<-function(read.counts,
                       condition,
                       design=NULL,
                       contrast=NULL,
                       group=NULL,
                       method=c('RUVr','RUVg','RUVs'),
                       cIdx=NULL,
                       k=1,...){
  start.time <- Sys.time()
  method <- match.arg(method,c('RUVr','RUVg','RUVs'))
  read.counts<-as.matrix(round(read.counts,0))
  
  ##########################################################################
  #---get the negative control
  if(is.null(cIdx) | method=='RUVr'){
    if(is.null(design)){
      condition<-factor(condition,levels = unique(condition))
      design<-model.matrix(~0+condition)
      colnames(design)<-gsub('condition','',colnames(design))
    }
    if(is.null(contrast)){
      contrast <- unique(condition)
      contrast <- paste0(contrast[-1],'-',contrast[1])
    }
    
    contrast <- makeContrasts(contrasts = contrast, levels=design)
    message('Estimate norm factor...')
    y <- DGEList(counts=read.counts, group=group)
    y <- calcNormFactors(y)
    message('Estimate Common Dispersion for Negative Binomial GLMs ...')
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    message('Fit genewise Negative Binomial Generalized Linear Models and calculate residuals...')
    fit <- glmQLFit(y, design)
    
    lrt <- glmQLFTest(fit, contrast = contrast)
    top <- topTags(lrt, n=nrow(read.counts))$table
    empirical <- top[order(top$PValue,decreasing = T),]
    empirical <- rownames(empirical)[empirical$PValue>0.1]
  }
  
  switch(method, 
         RUVr = {
           empirical <- rep(TRUE,dim(read.counts)[1])
           message('Remove Unwanted Variation Using RUVr...')
           res <- residuals(fit, type="deviance")
           results <-  RUVr(x = read.counts,cIdx = empirical,residuals = res,k=k)
         },
         RUVg = {
           message('Remove Unwanted Variation Using RUVg...')
           results <- RUVg(x = read.counts,cIdx = empirical,k=k)
         },
         RUVs = {
           message('Remove Unwanted Variation Using RUVs...')
           results <- RUVs(x = read.counts,cIdx = empirical,k=k,scIdx = makeGroups(condition))
         }
  )
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
  message('Done!!!')
  results$method <- method
  return(results)
}


remove.batch.shiny<-function(read.counts,
                             condition,
                             design=NULL,
                             contrast=NULL,
                             group=NULL,
                             method=c('RUVr','RUVg','RUVs'),
                             cIdx=NULL,
                             k=1,...){
  withProgress(message = 'Estimate batch effects...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0.1)
                 start.time <- Sys.time()
                 method <- match.arg(method,c('RUVr','RUVg','RUVs'))
                 read.counts<-as.matrix(round(read.counts,0))
                 
                 ##########################################################################
                 #########################################################################
                 #---get the negative control
                 if(is.null(cIdx) | method=='RUVr'){
                   if(is.null(design)){
                     condition<-factor(condition,levels = unique(condition))
                     design<-model.matrix(~0+condition)
                     colnames(design)<-gsub('condition','',colnames(design))
                   }
                   if(is.null(contrast)){
                     contrast <- unique(condition)
                     contrast <- paste0(contrast[-1],'-',contrast[1])
                   }
                   
                   contrast <- makeContrasts(contrasts = contrast, levels=design)
                   
                   
                   message('Estimate norm factor...')
                   incProgress(0.2)
                   y <- DGEList(counts=read.counts, group=group)
                   y <- calcNormFactors(y)
                   
                   message('Estimate Common Dispersion for Negative Binomial GLMs ...')
                   incProgress(0.3)
                   y <- estimateGLMCommonDisp(y, design)
                   y <- estimateGLMTagwiseDisp(y, design)
                   
                   message('Fit genewise Negative Binomial Generalized Linear Models and calculate residuals...')
                   incProgress(0.4)
                   fit <- glmQLFit(y, design)
                   
                   
                   lrt <- glmQLFTest(fit, contrast = contrast)
                   top <- topTags(lrt, n=nrow(read.counts))$table
                   empirical <- top[order(top$PValue,decreasing = T),]
                   empirical <- rownames(empirical)[empirical$PValue>0.1]
                 }
                 
                 switch(method, 
                        RUVr = {
                          empirical <- rep(TRUE,dim(read.counts)[1])
                          message('Remove Unwanted Variation Using RUVr...')
                          res <- residuals(fit, type="deviance")
                          results <-  RUVr(x = read.counts,cIdx = empirical,residuals = res,k=k)
                        },
                        RUVg = {
                          message('Remove Unwanted Variation Using RUVg...')
                          results <- RUVg(x = read.counts,cIdx = empirical,k=k)
                        },
                        RUVs = {
                          message('Remove Unwanted Variation Using RUVs...')
                          results <- RUVs(x = read.counts,cIdx = empirical,k=k,scIdx = makeGroups(condition))
                        }
                 )
                 incProgress(0.7)
                 end.time <- Sys.time()
                 time.taken <- end.time - start.time
                 message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
                 showNotification("Done!!!")
                 message('Done!!!')
                 results$method <- method
                 return(results)
               })
}

