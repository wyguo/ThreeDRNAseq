#' Filtering the scores for isoform switch
#'
#' Filtering the scores output from \code{\link{iso.switch}}.
#'
#' Users can set cut-offs, such as for the probability/frequency of switch and sum of average differences, to further refine the switch results.
#'
#' @param scores the scores object output from \code{\link{iso.switch}}.
#' @param prob.cutoff,diff.cutoff,t.points.cutoff,pval.cutoff,FDR.cutoff,cor.cutoff the cut-offs corresponding to switch frequencies/probablities,
#' sum of average sample differences, p-value, adjusted p-value (FDR) and time points cut-offs for both intervals before and after switch and Pearson correlation.
#' @param data.exp,mapping the expression and gene-isoform mapping data.
#' @param sub.isoform.list a vector of isoforms to output the corresponding results (see \code{TSIS.data$sub.isoforms}).
#' @param sub.isoform logical, to output subset of the results(TRUE) or not (FALSE). If TRUE, \code{sub.isoform.list} must be provided.
#' @param max.ratio logical, to show maximum abundant isoform results(TRUE) or not (FALSE). If TRUE, data.exp and mapping data must be
#' provided to calculate the isoform ratios to the genes using  \code{\link{rowratio}}.
#' @param x.value.limit the region of x axis (time) for investigation.
#'
#' @export
#'
#' @return a table of scores after filtering.
#'
score.filter<-function(scores,
                       prob.cutoff=0.5,
                       diff.cutoff=1,
                       t.points.cutoff=2,
                       pval.cutoff=0.01,
                       FDR.cutoff=0.01,
                       cor.cutoff=0.5,
                       data.exp=NULL,
                       mapping=NULL,
                       sub.isoform.list=NULL,
                       sub.isoform=F,
                       max.ratio=F,
                       x.value.limit=c(9,17)){

  if((is.null(data.exp) | is.null(mapping)) & max.ratio){
    msg<-data.frame('Expression data and mapping data are not provided.',row.names = NULL)
    colnames(msg)<-'Warnings:'
    return(msg)
  } else if(is.null(sub.isoform.list) & sub.isoform){
    msg<-data.frame('Subset of isoform names is not provided.',row.names = NULL)
    colnames(msg)<-'Warnings:'
    return(msg)
  } else {
    if(cor.cutoff==0){
      scores=scores[which(scores$before.t.points>=t.points.cutoff
                          & scores$after.t.points >= t.points.cutoff
                          & scores$prob>prob.cutoff
                          & scores$diff >diff.cutoff
                          & scores$before.pval < pval.cutoff
                          & scores$after.pval< pval.cutoff
                          & scores$before.FDR < FDR.cutoff
                          & scores$after.FDR< FDR.cutoff
                          & scores$x.value >=x.value.limit[1]
                          & scores$x.value <=x.value.limit[2]),] 
    } else {
      scores=scores[which(scores$before.t.points>=t.points.cutoff
                          & scores$after.t.points >= t.points.cutoff
                          & scores$prob>prob.cutoff
                          & scores$diff >diff.cutoff
                          & scores$before.pval < pval.cutoff
                          & scores$after.pval< pval.cutoff
                          & scores$before.FDR < FDR.cutoff
                          & scores$after.FDR< FDR.cutoff
                          & abs(scores$cor)>cor.cutoff
                          & scores$x.value >=x.value.limit[1]
                          & scores$x.value <=x.value.limit[2]),]
    }

    scores<-scores[order(scores[,'prob'],decreasing = T),]
    rownames(scores)<-NULL

    ##subsets of the data

    if(sub.isoform & nrow(scores)>0)
      scores<-scores[which(scores$iso1 %in% sub.isoform.list | scores$iso2 %in% sub.isoform.list),]

    if(max.ratio & nrow(scores)>0){
      genes0<-mapping[,1]
      isoforms0<-mapping[,2]

      ##expression ratio aboundance
      data.exp.mean.ratio<-apply(data.exp,1,mean)
      data.exp.mean.ratio<-data.frame(usage=rowratio(x = data.exp.mean.ratio,group = genes0))
      max.ratio.idx<-by(data.exp.mean.ratio,INDICES = genes0,function(x){
        rownames(x)[x==max(x)]
      },simplify = F)
      max.ratio.idx<-do.call(cbind,max.ratio.idx)
      scores<-scores[which(scores$iso1 %in% max.ratio.idx | scores$iso2 %in% max.ratio.idx),]
    }
    return(scores)

  }


}


