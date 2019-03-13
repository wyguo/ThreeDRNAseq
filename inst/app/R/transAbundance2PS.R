#' Calculate transcript PS (percent spliced) and \eqn{\Delta}PS from abundance
#' 
#' @details PS is defined as the ratio of transcript abundance divided by the gene abundance 
#' (sum of all transcript abundance of the same gene).
#' @param transAbundance matrix of transcript level abundance, e.g. TPM and read count exprssion.
#' @param PS PS value matrix. If a PS matrix provided, the \eqn{\Delta}PS values is directly calculated from the PS matrix 
#' based on the contrast groups, otherwise, new PS value will be generate. Default is \code{NULL}.
#' @param contrast a vector of contrast groups to calculate \eqn{\Delta}PS values, e.g. c('B-A','C-A'), which compares conditions 
#' "B" and "C" to condition "A".
#' @param condition a vector of condition labels which match the columns in the abundance matrix, 
#' e.g. c('A','A','A','B','B','B','C','C','C'), which means the abundance matrix corresponds to 3 conditions and each condition 
#' has 3 replicates.
#' @param mapping transcript-gene mapping matrix. First column is transcript list with column name "TXNAME" and 
#' second column is gene list with column name "GENEID".
#' @return a list of two data.frames: "PS" and "deltaPS".
#' @export
#' @examples 
#' set.seed(36)
#' transAbundance <- matrix(rnorm(36,mean=10,sd=1),ncol = 9,nrow = 3)
#' PS <- NULL
#' contrast <- c('B-A','C-A')
#' condition <- c('A','A','A','B','B','B','C','C','C')
#' mapping <- data.frame(TXNAME=c('trans1.1','trans1.2','trans2'),
#'                       GENEID=c('gene1','gene1','gene2'))
#' transAbundance2PS(transAbundance = transAbundance,PS = PS,
#'                   contrast = contrast,condition = condition,mapping = mapping)
#'                  

transAbundance2PS <- function(transAbundance=NULL,
                              PS=NULL,
                              contrast,
                              condition,
                              mapping){
  colnames(mapping) <- c('TXNAME','GENEID')
  if(!is.null(transAbundance))
    rownames(transAbundance) <- mapping$TXNAME
  
  if(is.null(PS)){
    transAbundance.mean <- t(rowmean(t(transAbundance),group = condition,reorder = F))
    ##---PS
    PS <- rowratio(x = transAbundance.mean[as.vector(mapping$TXNAME),],
                   group = as.vector(mapping$GENEID), reorder = F)
    PS <- PS[mapping$TXNAME,]
  }

  contrast.matrix <- makeContrasts(contrasts = contrast,levels = factor(condition,levels = unique(condition)))
  deltaPS <- data.frame(as.matrix(PS)%*%as.matrix(contrast.matrix),check.names = F)
  ##---deltaPSI
  # deltaPS <- lapply(strsplit(contrast,'-'),function(x){
  #   PS[,x[1]]-PS[,x[2]]
  # })
  # deltaPS <- do.call(cbind,deltaPS)
  colnames(deltaPS) <- contrast
  rownames(PS) <- rownames(deltaPS) <- mapping$TXNAME
  list(PS=PS,deltaPS=deltaPS)
}
