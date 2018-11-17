#' Filter low expression based on CPM (count per million reads)
#' @details An expressed target must have \eqn{\geq} \code{sample.n} with CPM \eqn{\geq} \code{cpm.cut}.
#' @param counts a data.frame or matrix of read counts.
#' @param sample.n number of samples
#' @param cpm.cut a numeric value of CPM cut-off
#' @return A vector of logcial values (TRUE/FALSE) to indicate the rows of the input data to keep/filter.
#' @export
#' @examples 
#' x <- rbind(matrix(rpois(50,lambda=4),5,10),matrix(rpois(50,lambda=0.1),5,10))
#' rownames(x) <- letters[1:10]
#' filter.idx <- CPM.filter(x,sample.n = 3,cpm.cut = 1)
#' ##The remaining
#' x[filter.idx,]
#' rownames(x)[filter.idx]
CPM.filter<-function(counts,sample.n=3,cpm.cut=1,...){
  rowSums(counts2CPM(counts,...)>=cpm.cut)>=sample.n
}

#' Filter low expression based on TPM (transcript per million reads)
#' @details An expressed target must have \eqn{\geq} \code{sample.n} with TPM \eqn{\geq} \code{cpm.cut}.
#' @param counts a data.frame or matrix of read counts.
#' @param sample.n number of samples
#' @param tpm.cut a numeric value of TPM cut-off
#' @return A vector of logcial values (TRUE/FALSE) to indicate the rows of the input data to keep/filter.
#' @export
TPM.filter<-function(TPM,sample.n=3,tpm.cut=1){
  rowSums(TPM>=tpm.cut)>=sample.n
}

#' Filter low expressed transcripts and genes based on read counts
#' @details Read counts are converted to CPM (count per million reads) to filter low expressed transcripts. 
#' 
#' 1) An expressed transcript must have \eqn{\geq} \code{sample.n} with CPM \eqn{\geq} \code{cpm.cut}.
#' 
#' 2) An expressed gene must have at least one expressed transcript.
#' 
#' @param counts read counts at transcript level.
#' @param mapping a data.frame of transcript-gene name mapping, with first column of transcript list and second column of gene list.
#' @param cpm.cut CPM cut-off. 
#' @param sample.n number of samples cut across all samples.
#' @param Log if TRUE, the log2-CPM is used for filters.
#' @return  A list of following elements:
#'     \itemize{
#'       \item{ trans_high: }{ the expressed transcripts}
#'       \item{ genes_high: }{ the expressed genes}
#'       \item{ mapping: }{ the remaining mapping data.frame after filter the low expressed transcripts}
#'     }
#' @export
#' 
counts.filter <- function(counts,mapping,cpm.cut=1,sample.n=3,Log=F){
  colnames(mapping) <- c('TXNAME','GENEID')
  rownames(mapping) <- mapping$TXNAME
  filter.idx<-CPM.filter(counts=counts,sample.n = sample.n,cpm.cut =cpm.cut,Log=Log)
  trans_high <- names(filter.idx)[filter.idx==T]
  genes_high <- unique(mapping[trans_high,]$GENEID)
  list(trans_high=trans_high,genes_high=genes_high,mapping_high=mapping[trans_high,])
}
