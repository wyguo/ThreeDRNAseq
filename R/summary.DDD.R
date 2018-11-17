#' Summary DE genes and transcripts
#' @param stat a data.frame object with first column of "target", second column of "contrast",
#' third and fourth columns of DE gene/transcript statistics (i.e. "adj.pval" and "log2FC", respectively). 
#' @param cutoff a numeric vector of cut-offs to "adj.pval" and "log2FC", e.g. \code{cutoff=c(adj.pval=0.01,log2FC=1)}.
#' @return a data.frame object, which is a subset of input \code{stat} after applying the \code{cutoff}.
#' 
summary.DE.target <- function(stat,cutoff=c(adj.pval=0.01,log2FC=1)){
  names(cutoff) <- c('adj.pval','log2FC')
  stat$up.down <- 'up_regulate'
  stat$up.down[stat[,'log2FC']<0] <- 'down_regulate'
  idx <- (abs(stat[,'adj.pval'])<=cutoff['adj.pval']) & (abs(stat[,'log2FC'])>=cutoff['log2FC'])
  stat <- stat[idx,]
  # split(stat , f = stat$contrast)
}

#' Summary DAS genes and DTU transcripts
#' @param stat a data.frame object with first column of "target", second column of "contrast",
#' third and fourth columns of DAS gene/DTU transcript statistics (i.e. "adj.pval" and "maxdeltaPS"/"deltaPS", respectively). 
#' @param cutoff a numeric vector of cut-offs for the statistics.
#' @return a data.frame object, which is a subset of input \code{stat} after applying the \code{cutoff}.
#' 
summary.DAS.target <- function(stat,lfc,cutoff=c(adj.pval=0.01,deltaPS=0.1)){
  names(cutoff) <- c('adj.pval','deltaPS')
  lfc <- lfc[which(lfc$target %in% stat$target),]
  stat <- merge(stat,lfc)
  stat$up.down <- 'up_regulate'
  stat$up.down[stat[,'log2FC']<0] <- 'down_regulate'
  idx <- (abs(stat[,'adj.pval'])<=cutoff['adj.pval']) & (abs(stat[,grep('deltaPS',colnames(stat))])>=cutoff['deltaPS'])
  stat <- stat[idx,]
  # split(stat , f = stat$contrast)
}

#' Compare DE and DAS results
#' @param DE_genes,DAS_genes,DE_trans,DTU_trans data.frame objects of DE genes, DAS genes, DE transcripts and DTU transcripts,
#' which are the outputs of \code{\link{summary.DE.target}} and \code{\link{summary.DAS.target}}.
#' @param contrast a vector of contrast groups, e.g. \code{contrast = c('B-A','C-A')}, which compares condition B and C to 
#' condition A.
#' @details In each contrast group, \code{DEvsDAS} compares the DE and DAS genes; \code{DEvsDTU} compares the DE and DTU
#' transcripts; and \code{summary.DDD.number} summarises the DE/DAS/DTU gene/transcript numbers.
#' 
#' @return a data.frame with first column of contrast goups, second column of DE only genes/transcripts, third column of
#' DE&DAS genes or DE&DTU transcripts and fourth column of DAS only genes or DTU only transcripts.
#' @export
#' @rdname DEvsDAS
DEvsDAS <- function(DE_genes,DAS_genes,contrast){
  idx <- c('DEonly','DE&DAS','DASonly')
  x <-  split(DE_genes$target,DE_genes$contrast)
  y <-  split(DAS_genes$target,DAS_genes$contrast)
  
  num <- lapply(contrast,function(i){
    x <- set2(x[[i]],y[[i]])
    names(x) <- idx
    x
  })
  names(num) <- contrast
  
  n1 <- lapply(contrast,function(i){
    data.frame(contrast=i,t(sapply(num[[i]],length)))
  })
  n1 <- do.call(rbind,n1)
  colnames(n1) <- c('Contrast',idx)
  
  if(length(contrast)==1){
    n1
  } else {
    num <- unlist(num,recursive=FALSE)
    n2 <- sapply(idx,function(i){
      length(Reduce(intersect,num[grep(paste0('\\.',i,'$'),names(num))]))
    })
    n2 <- data.frame(Contrast='Intersection',t(n2))
    names(n2) <- c('Contrast',idx)
    rbind(n1,n2)
  }
}

#' @export
#' @rdname DEvsDAS
DEvsDTU <- function(DE_trans,DTU_trans,contrast){
  idx <-  c('DEonly','DE&DTU','DTUonly')
  x <-  split(DE_trans$target,DE_trans$contrast)
  y <-  split(DTU_trans$target,DTU_trans$contrast)
  
  num <- lapply(contrast,function(i){
    x <- set2(x[[i]],y[[i]])
    names(x) <- idx
    x
  })
  names(num) <- contrast
  
  n1 <- lapply(contrast,function(i){
    data.frame(contrast=i,t(sapply(num[[i]],length)))
  })
  n1 <- do.call(rbind,n1)
  colnames(n1) <- c('Contrast',idx)
  
  if(length(contrast)==1){
    n1
  } else {
    num <- unlist(num,recursive=FALSE)
    n2 <- sapply(idx,function(i){
      length(Reduce(intersect,num[grep(paste0('\\.',i,'$'),names(num))]))
    })
    n2 <- data.frame(Contrast='Intersection',t(n2))
    names(n2) <- c('Contrast',idx)
    rbind(n1,n2)
  }
}

#' @export
#' @rdname DEvsDAS
summary.DDD.number <- function(DE_genes,DAS_genes,DE_trans,DTU_trans,contrast){
  # n1 <- lapply(DE_genes,function(x) x$target)
  idx <- factor(DE_genes$contrast,levels = contrast)
  n1 <- split(DE_genes$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x<- data.frame(`DE genes`=c(sapply(n1,length),`Intersection`=n2))
  
  # n1 <- lapply(DAS_genes,function(x) x$target)
  idx <- factor(DAS_genes$contrast,levels = contrast)
  n1 <- split(DAS_genes$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DAS genes`=c(sapply(n1,length),`Intersection`=n2)))
  
  # n1 <- lapply(DE_trans,function(x) x$target)
  idx <- factor(DE_trans$contrast,levels = contrast)
  n1 <- split(DE_trans$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DE transcripts`=c(sapply(n1,length),`Intersection`=n2)))
  
  
  # n1 <- lapply(DTU_trans,function(x) x$target)
  idx <- factor(DTU_trans$contrast,levels = contrast)
  n1 <- split(DTU_trans$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DTU transcripts`=c(sapply(n1,length),`Intersection`=n2)))
  
  x <- data.frame(contrast=rownames(x),x,row.names = NULL)
  colnames(x) <- c('contrast','DE genes','DAS genes','DE transcripts','DTU transcripts')
  if(length(contrast)==1)
    x <- x[-nrow(x),,drop=F]
  x
}

