#######################################################################################################
###---> Install packages from Cran
packages = c('shiny','shinydashboard','rhandsontable','shinyFiles','shinyjs','DT',
             'plotly','ggplot2','eulerr','ggrepel',
             'gridExtra','fastcluster','Gmisc','rmarkdown')
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE,repos = 'https://cran.ma.imperial.ac.uk/')
    library(x, character.only = TRUE)
  }
})


###---> Install packages from Bioconductor
bioconductor.package.list <- c('tximport','edgeR','limma','RUVSeq','ComplexHeatmap','rhdf5')
for(i in bioconductor.package.list){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = 'https://cran.ma.imperial.ac.uk/')
  if(!(i %in% rownames(installed.packages()))){
    message('Installing package: ',i)
    BiocManager::install(i, version = "3.8")
  } else next
}

####Other packages may be required for the analysis, please install accordingly.

##Load R packages
library(shiny)
library(shinydashboard)
library(rhandsontable)
library(shinyFiles)
library(shinyjs)
library(DT)
library(tximport)
library(edgeR)
library(limma)
library(plotly)
library(ggplot2)
library(RUVSeq)
library(eulerr)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(fastcluster)
library(Gmisc)
options(stringsAsFactors=F)











##########################################################################
# ========================== define R functions ======================== #
#' Sum Over Replicate Arrays
#' 
#' @param x a matrix-like object.
#' @param group grouping of sample identifier.
#' @return A data object of the same class as x with a column for sums according grouping.
#' @examples 
#' set.seed(100)
#' x<- matrix(rnorm(8*4),8,4)
#' colnames(x) <- c("a","a","b","b")
#' sumarrays(x)
#' data.frame(a=x[,1]+x[,2],b=x[,3]+x[,4])
#' @export

sumarrays<-function(x,group=NULL){
  if(is.null(group))
    group<-colnames(x)
  colnames(x)<-group
  
  x<-rowsum(t(x),group = group)
  #order the columns as the input group odering
  x<-data.frame(t(x)[,unique(group)])
  colnames(x) <- unique(group)
  return(x)
}


##########################################################################
######################################################################################################################
#' Bar plot of up- and down-regulation
#' @param data2plot a matrix/data.frame with three columns: "contrast", "regulation" and "number". See examples for details.
#' @param contrast a vector of contrast groups, e.g. c('C-A','B-A').
#' @param plot.title the the plot titile.
#' @param angle the angle to rotate x-axis labels. 
#' @param hjust numeric value of horizontal adjust for x-axis labels.
#' @return a bar plot in \code{ggplot} format.
#' @export
#' @examples
#' data2plot <- data.frame(contrast = c('T10-T2','T10-T2','T19-T2','T19-T2'),
#'                         regulation = c('down_regulate','up_regulate','down_regulate','up_regulate'),
#'                         number = c(305,727,1062,1805)
#'                         )
#' plotUpdown(data2plot = data2plot,contrast = c('T10-T2','T19-T2'),plot.title = 'Up-down regulation')
#'
plotUpdown <- function(data2plot,contrast,plot.title=NULL,angle=0,hjust=0.5){
  data2plot$contrast <- factor(data2plot$contrast,levels = unique(data2plot$contrast))
  data2plot$regulation <- factor(data2plot$regulation,levels = c('up_regulate','down_regulate'))
  g <- ggplot(data2plot,aes(x=contrast,y=number,group=regulation,
                            fill=regulation,label=number))+
    geom_bar(stat='identity')+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme_bw()+
    labs(title=plot.title)
  g <- g+theme(axis.text.x = element_text(angle = angle, hjust = hjust))
  g
}

######################################################################################################################
#' Plot gene and transcript expression profiles
#' @param data.exp a data.frame of transcript level expression, e.g. read counts and TPM.
#' @param gene a gene name to plot.
#' @param mapping a data.frame of transcript-gene mapping (first column is transcript list and second column is gene list).
#' @param genes.ann a data.frame of gene annotation with first column "target" and second column "annotation". Default is \code{NULL}.
#' @param trans.ann a data.frame of transcript annotation with first column "target" and second column "annotation".
#' Default is \code{NULL}.
#' @param trans.expressed a vector of expressed transcripts. If not \code{NULL}, only the profiles of provided transcripts
#' will be shown in the plot.
#' @param reps a vecotr of replicate labels, which provide the grouping information to calculate the average expression in each conditions.
#' @param x.lab,y.lab characters for x-axis and y-axis labels, respectively.
#' @param marker.size a number passed to \code{geom_point} to control the point size in the plot.
#' @param legend.text.size a number to control the legend text size.
#' @param error.type error type to make the error bars. Options are: "stderr" for standard error and
#' "sd" for standard deviation.
#' @param error.bar logical, whether to show error bars on the profile plot.
#' @param show.annotation logical, whether to show the annotations provided in \code{genes.ann} and \code{trans.ann} on the plot.
#' @param plot.gene logical, whether to show the gene level expression on the plot. The gene expression is the total of all the
#' transcript expression.
#' @details \code{plotAbundance} is used to plot the gene and/or trnascript level abundance while \code{plotPS} is used to plot
#' the percent spliced (PS) of transcript. PS is defined as the ratio of transcript expression to the total of all transcript
#' expression (i.e. gene expression) in the same gene.
#'
#' @return a plot in \code{ggplot} format.
#' @export
#' @rdname plot.abundance
#' @examples
#' data(exp.data)
#' plotAbundance(data.exp = exp.data$trans_TPM,
#'               gene = 'AT1G01020',
#'               mapping = exp.data$mapping,
#'               genes.ann = exp.data$genes.ann,
#'               trans.ann = exp.data$trans.ann,
#'               reps = exp.data$samples$condition)
#' 
#' plotPS(data.exp = exp.data$trans_TPM,
#'               gene = 'AT1G01020',
#'               mapping = exp.data$mapping,
#'               genes.ann = exp.data$genes.ann,
#'               trans.ann = exp.data$trans.ann,
#'               reps = exp.data$samples$condition)
#' @seealso \code{\link{data.error}}
plotAbundance<- function(
  data.exp,
  gene,
  mapping,
  genes.ann=NULL,
  trans.ann=NULL,
  trans.expressed=NULL,
  reps,
  x.lab='Conditions',
  y.lab='TPM',
  marker.size=3,
  legend.text.size=11,
  error.type='stderr',
  error.bar=T,
  show.annotation=T,
  plot.gene=T
){
  #####################################################################
  ##prepare plot data
  #####################################################################
  rownames(mapping) <- mapping$TXNAME
  if(!is.null(trans.expressed)){
    data.exp <- data.exp[trans.expressed,]
    mapping <- mapping[trans.expressed,]
  }
  rownames(genes.ann) <- genes.ann$target
  rownames(trans.ann) <- trans.ann$target
  trans <- mapping$TXNAME[mapping$GENEID==gene]
  plot.title <- paste0('Gene: ',gene)
  # expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  # trans.idx <- trans.idx[expression.sum>0]
  
  if(length(trans)==0){
    message(paste0(gene, ' is not in the dataset'))
    return(NULL)
  }
  
  data2plot<-data.exp[trans,,drop=F]
  
  if(length(trans)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)
  
  genes.idx <- gene
  trans.idx <- trans
  if(show.annotation){
    if(gene %in% genes.ann$target){
      idx <- genes.ann[gene,'annotation']
      genes.idx <- paste0(gene,' (',idx,')')
    }
    
    if(any(trans %in% trans.ann$target)){
      trans.idx <- trans
      idx <- which(trans %in% trans.ann$target)
      trans.idx[idx] <- paste0(trans.idx[idx],' (',trans.ann[trans.idx[idx],'annotation'],')')
    }
  }
  rownames(data2plot) <- c(genes.idx,trans.idx)
  legend.ncol <- ceiling(length(trans.idx)/15)
  
  ##--mean and error
  mean2plot <- t(rowmean(t(data2plot),group = reps))
  sd2plot <- by(t(data2plot),INDICES = reps,FUN = function(x){
    apply(x,2,function(y){
      data.error(y,error.type = error.type)
    })
  })
  sd2plot <- do.call(cbind,sd2plot)
  sd2plot <- sd2plot[rownames(mean2plot),]
  
  mean2plot <- reshape2::melt(mean2plot)
  sd2plot <- reshape2::melt(sd2plot)
  colnames(mean2plot) <- c('Targets','Conditions','mean')
  colnames(sd2plot) <- c('Targets','Conditions','error')
  data2plot <- merge(mean2plot,sd2plot)
  
  ##--plot gene or not
  if(!plot.gene){
    data2plot <- data2plot[-which(data2plot$Targets %in% genes.idx),]
    plot.title <- paste0('Gene: ',genes.idx)
  }
  
  
  data2plot$Conditions <- factor(data2plot$Conditions,levels = unique(reps))
  
  profiles <- ggplot(data2plot,aes(x=Conditions,y=mean,group=Targets,color=Targets))+
    geom_line(size=1)+
    geom_point(size=marker.size,aes(fill=Targets,shape=Targets))+
    scale_shape_manual(name="Targets",values=c(25:0,25:0,rep(25,500)))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))
  if(error.bar)
    profiles <- profiles+
    geom_errorbar(aes(ymin=mean-error,ymax=mean+error),width=0.1,color='black',size=0.3)
  profiles
}

#' @rdname plot.abundance
#' @export
plotPS <- function(
  data.exp,
  gene,
  mapping,
  genes.ann=NULL,
  trans.ann=NULL,
  trans.expressed=NULL,
  reps,
  y.lab='TPM',
  x.lab='Conditions',
  marker.size=3,
  legend.text.size=11,
  show.annotation=T
){
  #####################################################################
  ##prepare plot data
  #####################################################################
  rownames(mapping) <- mapping$TXNAME
  if(!is.null(trans.expressed)){
    data.exp <- data.exp[trans.expressed,]
    mapping <- mapping[trans.expressed,]
  }
  rownames(genes.ann) <- genes.ann$target
  rownames(trans.ann) <- trans.ann$target
  trans <- mapping$TXNAME[mapping$GENEID==gene]
  # expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  # trans.idx <- trans.idx[expression.sum>0]
  
  if(length(trans)==0){
    message(paste0(gene, ' is not in the dataset'))
    return(NULL)
  }
  data2plot<-data.exp[trans,,drop=F]
  data2plot <- t(rowmean(t(data2plot),reps))
  
  if(length(trans)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)
  
  genes.idx <- gene
  trans.idx <- trans
  if(show.annotation){
    if(gene %in% genes.ann$target){
      idx <- genes.ann[gene,'annotation']
      genes.idx <- paste0(gene,' (',idx,')')
    }
    
    if(any(trans %in% trans.ann$target)){
      trans.idx <- trans
      idx <- which(trans %in% trans.ann$target)
      trans.idx[idx] <- paste0(trans.idx[idx],' (',trans.ann[trans.idx[idx],'annotation'],')')
    }
  }
  rownames(data2plot) <- c(genes.idx,trans.idx)
  legend.ncol <- ceiling(length(trans.idx)/15)
  plot.title <- paste0('Gene: ',genes.idx)
  ##PS value
  data2plot[trans.idx,] <- t(t(data2plot[trans.idx,])/data2plot[genes.idx,])
  data2plot <- data2plot[-which(rownames(data2plot) %in% genes.idx),,drop=F]
  data2plot <- reshape2::melt(data2plot)
  colnames(data2plot) <- c('Targets','Conditions','PS')
  data2plot$Conditions <- factor(data2plot$Conditions,levels = unique(reps))
  
  profiles <- ggplot(data2plot,aes(x=Conditions,y=PS,group=Targets,color=Targets))+
    geom_bar(stat='identity',aes(fill=Targets))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))
  profiles
}

######################################################################################################################
#' Bar plot of GO annotation
#' @param go.table a data.frame of GO annotation table, with first column "Category" (i.e. BP, CC and MF),
#' second column of Go "Term" and remaining columns with significant statistics, e.g. "Count", "-log10(FDR)", etc.
#' @param col.idx a character to indicate which column of statistics to use to report the significance, e.g.
#' \code{col.idx="-log10(FDR)"}. \code{col.idx} must match to the column names of \code{go.table}.
#' @param plot.title the titile to show on the plot.
#' @return a plot in \code{ggplot} format.
#' @export
#' @examples
#' data(go.table)
#' plotGO(go.table = go.table,col.idx = '-log10(FDR)',plot.title = 'GO annotation: DE genes')
#'
plotGO <- function(go.table,col.idx,plot.title='GO annotation'){
  class(go.table) <- c("GO", class(go.table))
  data2plot <- go.table[,c(1,2,which(colnames(go.table)==col.idx))]
  colnames(data2plot) <- c('Category','Term','Value')
  data2plot <- by(data2plot,data2plot$Category,function(x){
    x[order(x$Value,decreasing = T),]
  })
  data2plot <- do.call(rbind,data2plot)
  data2plot$Term <- factor(data2plot$Term,levels = rev(data2plot$Term))
  g <- ggplot(data2plot,aes(x=Term,y=Value))+
    geom_bar(stat='identity',aes(fill=Category))+
    coord_flip() +
    theme_bw()+
    theme(axis.title.y = element_blank())+
    facet_grid(Category~.,scales = 'free_y',space='free_y')+
    labs(y=col.idx,title=plot.title)
  g
}

######################################################################################################################
#' Make principal component analysis (PCA) plot of individual samples
#' @param data2pca a matrix/data.frame to make the PCA plot of which the rows are defined individuals to plot.
#' @param dim1,dim2 characters to indicate the PCs to plot on x-axis (dim1) and y-axis (dim2), e.g. \code{dim1='PC1'},
#' \code{dim2='PC2'}; \code{dim1='PC2'}, \code{dim2='PC3'}; etc.
#' @param groups a vector of characters to specify the groups of conditions. The scatter points will be coloured
#' according to provided group information.
#' @param plot.title the title to show on the plot.
#' @param ellipse.type add colour shadow to the sample groups. Options are: "none","ellipse" and "polygon".
#' @param add.label logical, whether to add sample labels to the scatter points.
#' @param adj.label logical, whether to adjust the position of sample labels to avoid overlapping.
#' @return a plot in \code{ggplot} format.
#' @export
#' @examples
#' data(exp.data)
#' library(edgeR)
#' library(RUVSeq)
#' library(ggplot2)
#' trans_counts <- exp.data$trans_counts
#' 
#' ##------> sample informationsamples <- exp.data$samples
#' samples <- exp.data$samples[,c('condition','brep','srep')]
#' 
#' ##------> sum sequencing replicates
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
#' plotPCAind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = groups,ellipse.type='polygon')
#' 
#' ##------> remove batch effects
#' trans_batch <- remove.batch(read.counts = y,
#'                       condition = samples$condition,
#'                       method = 'RUVr')
#' 
#' ##------> normalisation and plot PCA again
#' dge <- DGEList(counts=trans_batch$normalizedCounts)
#' dge <- suppressWarnings(calcNormFactors(dge))
#' data2pca <- t(counts2CPM(obj = dge,Log = T))
#' plotPCAind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = groups,ellipse.type='polygon')

plotPCAind <- function(data2pca,
                       dim1='PC1',dim2='PC2',
                       groups,
                       plot.title='PCA plot',
                       ellipse.type=c('none','ellipse','polygon'),
                       add.label=T,adj.label=F){
  
  ellipse.type <- match.arg(ellipse.type,c('none','ellipse','polygon'))
  dim1 <- toupper(dim1)
  dim2 <- toupper(dim2)
  
  fit <- prcomp(data2pca,scale = T)
  fit.stat <- summary(fit)$importance
  dim1.p <- round(fit.stat[2,dim1]*100,2)
  dim2.p <- round(fit.stat[2,dim2]*100,2)
  data2plot <- data.frame(groups=groups,dim1=fit$x[,dim1],dim2=fit$x[,dim2])
  data2plot$groups <- factor(data2plot$groups,levels = unique(data2plot$groups))
  data2plot$labels <- rownames(data2plot)
  
  g <- ggplot(data2plot,aes(x=dim1,y=dim2),frame=T)+
    geom_point(aes(colour=groups,shape=groups))+
    geom_hline(yintercept = 0,linetype=2)+
    geom_vline(xintercept=0,linetype=2)+
    theme_bw()+
    theme(panel.border = element_blank())+
    labs(x=paste0(dim1,' (',dim1.p,'%)'),
         y=paste0(dim2,' (',dim2.p,'%)'),
         title=plot.title)+
    scale_shape_manual(values = rep(c(18:0,19:25),3)[1:length(unique(groups))])+
    coord_cartesian(xlim = c(min(data2plot$dim1)*1.1,max(data2plot$dim1)*1.1),
                    ylim = c(min(data2plot$dim2)*1.1,max(data2plot$dim2)*1.1))
  
  if(ellipse.type=='ellipse')
    g <- g+stat_ellipse(geom = "polygon",aes(fill=groups,colour=groups),alpha=0.2)
  if(ellipse.type=='polygon'){
    hulls <- plyr::ddply(data2plot, "groups", function(x){
      x[chull(x$dim1,x$dim2),]
    })
    g <- g+geom_polygon(data = hulls,aes(x=dim1,y=dim2,fill=groups),alpha = 0.3)
  }
  if(add.label & !adj.label)
    g <- g+geom_text(aes(label=labels,colour=groups),hjust=0, vjust=0)
  if(add.label & adj.label)
    g <- g+ggrepel::geom_text_repel(label=rownames(data2pca),aes(colour=groups))
  g
  
}


#
#' Boxplot of expression distribution in each sample.
#' @param data.before a numeric matrix of expression before normalisation, e.g. raw read counts.
#' @param data.after a numeric matrix of expression after normalisation, e.g. \eqn{\log_2}-CPM.
#' @param condition a vector of conditions, which is used to distinguish and colour the conditions of
#' samples.
#' @param sample.name a vector of sample names, which is used to name and order the samples in the boxplot.
#' @return a list with two elements: "g1" and "g2", which are \code{ggplot} format boxplots of the data
#' before and after normalisation.
#' @export
#'
boxplot.normalised <- function(data.before,data.after,condition,sample.name){
  data2plot <- data.before
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g1 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x='samples',y='read counts',title='Data distribution before normalisation')
  
  data2plot <- data.after
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g2 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x='samples',y='log2(CPM)',title='Data distribution after normalisation')
  
  if(length(sample.name)>20){
    g1 <- g1+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
    g2 <- g2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
  }
  return(list(g1=g1,g2=g2))
}

#' Plot Euler diagram
#' @param x (1) a list of individuals in different sets for comparisons or (2) a vecotr of numbers of set relations.
#' @param fill a vector of colours to fill the diagram. 
#' @param shape geometric shape used in the diagram.
#' @param ... arguments passed down to the \code{\link{eulerr::euler}} function.
#' @return a Euler diagram
#' @examples 
#' x <- letters[1:10]
#' y <- letters[5:15]
#' z <- set2(x,y)
#' combo <- c(x=length(z$x.only),y=length(z$y.only),"x&y"=length(z$xy))
#' plotEulerDiagram(list(x=x,y=y))
#' plotEulerDiagram(combo)
#' @export
#' @seealso \code{\link{eulerr::euler}}.
#' 
plotEulerDiagram <- function(x,
                             fill = gg.color.hue(length(x)),
                             shape = c("ellipse", "circle"),scale.plot=1,
                             ...){
  shape <- match.arg(shape,c("ellipse", "circle"))
  fit <- eulerr::euler(x,...)
  # grid.newpage()
  g <- plot(fit,quantities = TRUE,
            labels = list(font =1),
            fill=fill,shape = shape)
  vp = viewport(height=unit(1*scale.plot, "npc"), 
                width=unit(1*scale.plot, "npc"))
  g <- arrangeGrob(g,vp=vp)
  g
  # grid.newpage()
  # grid.rect(vp=vp,gp=gpar(lty=1, col="white", lwd=0))
  # grid.draw(g)
}

##########################################################################
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
##########################################################################
#' Summary DE genes and transcripts
#' @param stat a data.frame object with first column of "target", second column of "contrast",
#' third and fourth columns of DE gene/transcript statistics (i.e. "adj.pval" and "log2FC", respectively). 
#' @param cutoff a numeric vector of cut-offs to "adj.pval" and "log2FC".
#' @return a data.frame object, which is a subset of input \code{stat} after applying the \code{cutoff}.
#' 
summaryDEtarget <- function(stat,cutoff=c(adj.pval=0.01,log2FC=1)){
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
summaryDAStarget <- function(stat,lfc,cutoff=c(adj.pval=0.01,deltaPS=0.1)){
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
#' which are the outputs of \code{\link{summaryDEtarget}} and \code{\link{summaryDAStarget}}.
#' @param contrast a vector of contrast groups, e.g. \code{contrast = c('B-A','C-A')}, which compares condition B and C to 
#' condition A.
#' @details In each contrast group, \code{DEvsDAS} compares the DE and DAS genes; \code{DEvsDTU} compares the DE and DTU
#' transcripts; and \code{summary3Dnumber} summarises the DE/DAS/DTU gene/transcript numbers.
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
summary3Dnumber <- function(DE_genes,DAS_genes,DE_trans,DTU_trans,contrast){
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


#' Merge two testing statistics
#' @param x,y data.frame object of statistics (e.g. p-values, log2-fold changes, deltaPS), with rows of targets and
#' columns statistics in contrast groups.
#' @param contrast a vector of contrast groups, e.g. \code{contrast = c('B-A','C-A')}. If \code{NULL}, the column names of
#' \code{x} are used as \code{contrast}.
#' @param stat.type a vector of statistic types of \code{x} and \code{y}, respectively, e.g.
#' \code{stat.type = c('adj.pval','lfc')}. \code{stat.type} is passed to third and fourth column names of the output data.frame.
#' @param srot.by a column name to sort the output data.frame.
#' @return \code{summaryStat} returns a data.frame object with first column of target, second column of contrast groups, third column name of
#' statistics \code{x} and fourth column of statistic \code{y}.
#' @export
summaryStat  <- function(x,y,target,
                         contrast = NULL,
                         stat.type = c('adj.pval','lfc'),
                         srot.by = stat.type[1]){
  x <- x[target,,drop=F]
  y <- y[target,,drop=F]
  if(is.null(contrast))
    contrast <- colnames(x)
  x <- x[,contrast]
  y <- y[,contrast]
  
  stat <- lapply(contrast,function(i){
    z <- data.frame(targets =target,contrast=i, x[,i],y[,i],row.names = NULL)
    colnames(z) <- c('target','contrast',stat.type)
    z <- z[order(z[,srot.by]),]
  })
  names(stat) <- contrast
  stat <- do.call(rbind,stat)
  rownames(stat) <- NULL
  stat
}


##########################################################################
#' Plot flow-chart graph of DE and DAS results
#' @description The flow-chart shows the results of DE vs DAS genes or DE vs DTU transcripts.
#' For example, the expressed genes are divided into DE genes and not DE genes; the DE genes
#' are divided into DE only (only transcription regulation) and DE&DAS (both transcription
#' and AS regulation) genes; the not DE genes are divided into DAS only (only AS regulation) and
#' no expression/AS change (no regulation). The same divisions to DE and/or DTU transcripts.
#' @param expressed a vector of expressed genes/transcripts.
#' @param x a vector of DE genes/transcripts.
#' @param y a vector of DAS genes/DTU transcripts.
#' @param type a character to indicate the provided list type. Options are "genes" and "transcripts".
#' @param pval.cutoff,lfc.cutoff,deltaPS.cutoff the cut-offs used to determine the significant
#' genes/transcripts.
#' @return a flow-chart plot.
#' @export
plotFlowChart <- function(expressed,x,y,
                          type=c('genes','transcripts'),
                          pval.cutoff=0.01,
                          lfc.cutoff=1,
                          deltaPS.cutoff=0.1){
  type <- match.arg(type,c('genes','transcripts'))
  n.expressed <- length(expressed)
  n.x <- length(x)
  n.notx <- n.expressed-n.x
  n.onlyx <- length(setdiff(x,y))
  n.xy <- length(intersect(x,y))
  n.onlyy <- length(setdiff(y,x))
  n.no<- n.notx-n.onlyy
  n.y <- length(y)
  
  grid.newpage()
  width <- 0.2
  ##---Expressed
  exp.bx <- boxGrob(
    label = paste0('Expressed\n',length(expressed)),
    x = 0.1,
    y = 0.5,
    width = 0.14,
    box_gp=gpar(fill = 'white'),
    just = 'center')
  
  ##---DE genes
  de.genes.bx <- boxGrob(
    label = paste0('DE ',type,'\n',n.x),
    x = 0.35,
    y = 0.7,
    width = 0.16,
    box_gp=gpar(fill = 'gold'),
    just = 'center')
  
  ##---not DE genes
  notde.genes.bx <- boxGrob(
    label = paste0('Not\n DE ',type,'\n',n.notx),
    x = 0.35,
    y = 0.3,
    width = 0.16,
    box_gp=gpar(fill = 'white'),
    just = 'center')
  
  ##---Only transcription regulation
  de.only.bx <- boxGrob(
    label = paste0('Only transcription\nregulation\n',n.onlyx),
    x = 0.65,
    y = 0.8,
    width = width,
    box_gp=gpar(fill = 'gold'),
    just = 'center')
  
  ##---Transcription and AS regulation
  de.das.bx <- boxGrob(
    label = paste0('Transcription\nand AS regulation\n',n.xy),
    x = 0.65,
    y = 0.6,
    width = width,
    box_gp=gpar(fill = 'seagreen3'),
    just = 'center')
  
  ##---Only AS regulation
  das.only.bx <- boxGrob(
    label = paste0('Only AS\nregulation\n',n.onlyy),
    x = 0.65,
    y = 0.4,
    width = width,
    box_gp=gpar(fill = 'deepskyblue3'),
    just = 'center')
  
  ##---No regulation
  no.regulation.bx <- boxGrob(
    label = paste0('No\nregulation\n',n.no),
    x = 0.65,
    y = 0.2,
    width = width,
    box_gp=gpar(fill = 'white'),
    just = 'center')
  
  ##---DAS genes
  das.bx <- boxGrob(
    label = ifelse(type=='genes',paste0('DAS genes\n',n.y),paste0('DTU transcripts\n',n.y)),
    x = 0.9,
    y = 0.5,
    width = 0.15,
    box_gp=gpar(fill = 'deepskyblue3'),
    just = 'center')
  
  
  arrow_obj <- getOption('connectGrobArrow',
                         default = arrow(angle = 0, ends = 'last', type = 'closed'))
  lty_gp <- getOption('connectGrob', default = gpar(fill = 'black',lwd=4))
  plot(connectGrob(exp.bx, de.genes.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(exp.bx, notde.genes.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(de.genes.bx, de.only.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(de.genes.bx, de.das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(notde.genes.bx, das.only.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(notde.genes.bx, no.regulation.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  
  lty_gp <- getOption('connectGrob', default = gpar(col = 'red',lwd=4))
  plot(connectGrob(de.das.bx, das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  plot(connectGrob(das.only.bx, das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
  
  plot(exp.bx)
  plot(de.genes.bx)
  plot(notde.genes.bx)
  plot(de.only.bx)
  plot(de.das.bx)
  plot(das.only.bx)
  plot(no.regulation.bx)
  plot(das.bx)
  grid.text(ifelse(type=='genes','DE and DAS genes','DE and DTU transcripts'), x=0.5, y=0.95,
            gp=gpar(fontsize=16, col='black'))
  grid.text(paste0('P-value cut-off: ',pval.cutoff,
                   '\nLog2 fold change cut-off: ',lfc.cutoff,
                   '\ndeltaPS cut-off: ', deltaPS.cutoff),
            x=0.05, y=0.1, just = 'left',
            gp=gpar(fontsize=12, col='black'))
}


##########################################################################
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


##########################################################################
#' Dfferential expression (DE) and alternative splicing (DAS) analysis
#' @description Two pipelines are provided to perform differential expression (DE), differential
#' alternative splicing (DAS) gene and differential transcript usage (DTU) transcript analysis.
#' @param dge a numeric matrix of read counts or a \code{\link[edgeR]{DGEList}} object.
#' Counts must be non-negative and NAs are not permitted.
#' @param method method to use in the \code{\link{edgeR.pipeline}}. Options are "glm" or "glmQL".
#' @param contrast a vector of contrast groups for expression comparisons,
#' e.g. c('B-A','C-A') compares conditions "B" and "C" to condition "A".
#' @param span width of the \code{\link{lowess}} smoothing window
#' as a proportion, which is passed to the \code{\link[limma]{voom}} function to estimate
#' expression mean-variance trend.
#' @param design design matrix with rows of samples and columns of conditions to fit a linear model.
#' @param deltaPS a matrix of \eqn{\Delta}PS values, with rows of transcripts and columns of contrast groups.
#' @param diffAS logical, whether to perform DAS analysis. If \code{TRUE}, the \code{\link[limma]{diffSplice}} (limma pipeline)
#' or \code{\link[edgeR]{diffSpliceDGE}} (edgeR pipeline) function is used to peform DAS gene and DTU transcript analysis.
#' @param adjust.method method the adjust p-values for multiple comparisons. Default is "BH" and full details can be found in \code{\link{p.adjust}.}
#' @details
#' \bold{Abbreviation}
#' \itemize{
#' \item{DE: differential expressed genes/transcripts.}
#' \item{DAS: differential alternative spliced genes.}
#' \item{DTU: differential transcript usage transcripts.}
#' }
#'
#' \bold{Statistics}
#' \itemize{
#'     \item{
#'     Contrast groups: refer to the conditions to compare, e.g. "B-A" and "C-A" compares the conditions
#'     "B" and "C" to condition "A".
#'     }
#'     \item{
#'     Adjusted p-value: in the expression comparitive analysis, each gene/transcript is fitted with a statistical model, reporting a p-value.
#'     However the testing is performed over a substantial number of genes/transcripts and each testing has a probability of making errors.
#'     Hence, all p-values in the analsysis are adjusted by controlling the false discovery rate (FDR). Please see the R function \code{\link{p.adjust}}
#'     for details.
#'     }
#'     \item{
#'     \eqn{L_2FC}: \eqn{\log_2} fold change of expression based on contrast groups.
#'    }
#'    \item{
#'    \eqn{\Delta}PS: Transcript level PS (percent of splice) is defined as the ratio of transcript average abundance of conditions divided by the
#'    average gene abundance. \eqn{\Delta}PS is the PS differences of conditions based on the contrast groups. The abundance to calculate PS values
#'    can be TPMs or read counts.
#'    }
#' }
#' \bold{DE gene/transcript analysis}
#'
#' A gene/transcript is identified as DE in a contrast group if \eqn{L_2FC} of expression \eqn{\geq} a cut-off and with adjusted p-value < a cut-off.
#' 
#' Two pipelines are provided for DE analysis: 
#' 
#' ``limma-voom'' and ``edgeR''. The ``edgeR'' pipeline includes two methods: ``glmQL'' (Genewise Negative Binomial Generalized Linear Models with 
#' Quasi-likelihood Tests) and ``glm'' (Genewise Negative Binomial Generalized Linear Models). In limma pipeline, L2FCs are calculated based count 
#' per million reads (CPMs) while in edgeR, L2FCs are based on read counts. From several real RNA-seq data analyses, high consensus is achieved between 
#' ``limma-voom'' and ``glmQL'' (>90%) and they have more stringent controls of false discovery rate than the ``glm'' method.
#'
#' \bold{DAS gene and DTU transcript analysis}
#'
#' To test differential expression, the gene level expression is compared to transcript level expression in the contrast groups.
#' A gene is DAS in a contrast group if adjusted p-value < cut-off and at least one transcript of the gene with \eqn{\Delta}PS \eqn{\geq} a cut-off.
#' A transcript is DTU if adjusted p-value < a cut-off and \eqn{\Delta}PS \eqn{\geq} a cut-off.
#'
#' Two pipelines for AS analysis:
#' \itemize{
#'     \item{
#'     limma pipeline: To identify DAS genes, the expression of each transcript is compared to the weighted average expression of all the transcripts
#'     for the same gene (weight on transcript expression variance; the average expression can be treated as gene level expression). P-value of each
#'     test is converted to genewise p-value by using F-test or Simes method. To identify DTU transcript, each transcript is compared to the weighted
#'     average of all the other transcripts for the same gene. See the \code{\link{diffSplice}} function in \code{\link{limma}} package for details.
#'     }
#'     \item{
#'     edgeR pipeline: To identify DTU transcripts, log-fold-change of each transcript is compared to log-fold-change of entire gene
#'     (sum of all transcripts). The DTU transcript test p-values are summarised to genewise p-value by using F-test or Simes method
#'     to report DAS genes. See \code{\link{diffSpliceDGE}} function in \code{\link{edgeR}} R package for details.
#'     }
#' }
#' @return a list object with elements:
#' \itemize{
#'    \item{
#'    \code{voom.object:} the results of \code{\link[limma]{voom}} function to estimate expression mean-variance trend.
#'    }
#'    \item{
#'    \code{fit.lmFit:} the results of \code{\link[limma]{lmFit}} to fit a linear regression model on conditions.
#'    }
#'    \item{
#'    \code{fit.glm:} the results of \code{\link[edgeR]{glmFit}} in the edgeR pipeline.
#'    }
#'    \item{
#'    \code{fit.glmQL:} the  results of \code{\link[edgeR]{glmQLFit}} in the edgeR pipeline.
#'    }
#'    \item{
#'    \code{fit.contrast:} the results of \code{\link[limma]{contrasts.fit}}, i.e. incorporate the contrast groups to the linear regression model.
#'    }
#'    \item{
#'    \code{fit.eBayes:} the results of \code{\link[limma]{eBayes}}.
#'    }
#'    \item{
#'    \code{DE.pval.list:} a list object of which each element is the result of applying \code{\link[limma]{topTable}} or \code{\link[edgeR]{topTags}}
#'    to individual contrast groups.
#'    }
#'    \item{
#'    \code{DE.pval:} a data.frame object of adjusted p-values, with rows of genes/transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DE.lfc:} a data.frame object of \eqn{L_2FC}, with rows of genes/transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DE.stat:} a data.frame object with first column of target (genes/transcripts), second column of contrast (contrast groups), third column of
#'    adj.pval and fourth column of log2FC.
#'    }
#'    \item{
#'    \code{DE.stat.overalltest:} a data.frame object of testing statistics over all contrast groups rather than testing in individual contrast groups.
#'    }
#' }
#' If alternative splicing analysis is performed (\code{diffAS=TRUE}), more results are returned:
#' \itemize{
#'    \item{
#'    \code{fit.splice:} the results of \code{\link[limma]{diffSplice}} or \code{\link[edgeR]{diffSpliceDGE}} for individual contrast groups.
#'    }
#'    \item{
#'    \code{DTU.pval.list:} a list object of which each element is the DTU transcript result of \code{\link[limma]{topSplice}} or
#'    \code{\link[edgeR]{topSpliceDGE}} for individual contrast groups.
#'    }
#'    \item{
#'    \code{DTU.pval:} a data.frame object of adjusted p-values, with rows of transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DTU.deltaPS:} a data.frame object of \eqn{\Delta}PS, with rows of transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DTU.stat:} a data.frame object with first column of target (transcripts), second column of contrast groups, third column of adj.pval and
#'     fourth column of deltaPS.
#'    }
#'    \item{
#'    \code{maxdeltaPS:} a data.frame object with rows of genes and columns of maximum \eqn{\Delta}PS of the transcripts in the same gene in different
#'    contrast groups.
#'    }
#'    \item{
#'    \code{DAS.pval.F.list:} a list object of which each element is the DAS gene result of "F-test" using \code{\link[limma]{topSplice}} or
#'    \code{\link[edgeR]{topSpliceDGE}} in individual contrast groups.
#'    }
#'    \item{
#'    \code{DAS.pval.F:} a data.frame object of adjusted p-value with rows of genes and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DAS.F.stat:} the DAS gene statistics using "F-test".
#'    }
#'    \item{
#'    \code{DAS.pval.simes.list:} a list object of "Simes-test" using \code{\link[limma]{topSplice}} or
#'    \code{\link[edgeR]{topSpliceDGE}} in individual contrast groups.
#'    }
#'    \item{
#'    \code{DAS.pval.simes:} a data.frame object of DAS gene adjusted p-values using "Simes" method.
#'    }
#'    \item{
#'    \code{DAS.simes.stat:} the DAS gene statistics using "Simes-test".
#'    }
#' }
#' @export
#' @rdname limma.pipeline
#' @seealso R packages: \code{\link{limma}} and \code{\link{edgeR}}.

limma.pipeline <- function(dge,
                           contrast,
                           span=0.5,
                           design,
                           deltaPS=NULL,
                           diffAS=F,
                           adjust.method='BH'){
  start.time <- Sys.time()
  results <- list()
  ##########################################################
  ##--limma voom
  message('Limma-voon to estimate mean-vriance trend ...')
  voom.object<-voom(dge,design,plot=F,span = span)
  results$voom.object <- voom.object
  targets <- rownames(voom.object$E)
  
  ##########################################################
  ##--Fit a basic linear model
  message('Fit a basic linear model ...')
  fit.lmFit <- lmFit(voom.object, design)
  results$fit.lmFit <- fit.lmFit
  
  ##########################################################
  ##--Fit a basic linear model
  message('Fit the contrast model ...')
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  print(paste0('Contrast groups: ',paste0(contrast,collapse = '; ')))
  fit.contrast<-contrasts.fit(fit.lmFit, contrast.matrix)
  results$fit.contrast<-fit.contrast
  
  ##########################################################
  ##--Fit a eBayes model
  message('Fit a eBayes model ...')
  fit.eBayes<-eBayes(fit.contrast)
  results$fit.eBayes<-fit.eBayes
  
  ##########################################################
  ##--Testing statistics for each contrast group
  message('Testing for each contrast group ...')
  DE.pval.list <- lapply(contrast,function(i){
    x <- topTable(fit.eBayes,
                  coef=i,
                  adjust.method =adjust.method,
                  number = Inf)
    x <- x[targets,]
    x
  })
  names(DE.pval.list) <- contrast
  results$DE.pval.list<-DE.pval.list
  
  ###---DE pval and lfc
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$adj.P.Val))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast
  
  DE.stat <- summaryStat (x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))
  
  # DE.stat <- cbind(DE.pval,DE.lfc)
  # colnames(DE.stat) <- c(paste0('pval:',contrast),paste0('lfc:',contrast))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  
  ##########################################################
  ##--Testing statistics for across all contrast groups
  message('Testing across all contrast groups ...')
  DE.stat.overalltest<-topTable(fit.eBayes,number = Inf,
                                coef = contrast,adjust.method =adjust.method )
  # DE.stat.overalltest <- DE.stat.overalltest[targets,]
  col.idx <- gsub('-','.',contrast)
  col.idx <- grep(paste0(col.idx,collapse = '|'),colnames(DE.stat.overalltest))
  DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                    contrast='overall',DE.stat.overalltest,row.names = NULL)
  colnames(DE.stat.overalltest)[col.idx] <- gsub('[.]','-',colnames(DE.stat.overalltest)[col.idx])
  results$DE.stat.overalltest<-DE.stat.overalltest
  
  if(diffAS){
    message('Fit a splicing model ...')
    if(is.null(deltaPS))
      stop('Please provide deltaPS for DAS analysis...')
    fit.splice<-diffSplice(fit.contrast, geneid = 'GENEID')
    results$fit.splice<-fit.splice
    
    # ##########################################################
    ##---DTU transcripts
    ##DTU transcript pval list
    genes.idx <- unique(fit.splice$genes$GENEID)
    trans.idx <- unique(fit.splice$genes$TXNAME)
    
    DTU.pval.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="t", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast
    
    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast
    
    ##DTU transcript deltaPS
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]
    
    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))
    
    # DTU.stat <- cbind(DTU.pval,DTU.deltaPS)
    # colnames(DTU.stat) <- c(paste0('pval:',contrast),paste0('deltaPS:',contrast))
    
    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat
    
    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice$genes$TXNAME,,drop=F],
                     INDICES = fit.splice$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS
    
    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="F", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.F.list) <- contrast
    
    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast
    
    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    
    # DAS.F.stat <- cbind(DAS.pval.F,maxdeltaPS)
    # colnames(DAS.F.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    #
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat
    
    ##---simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.simes.list) <- contrast
    
    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast
    
    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))
    #
    # DAS.simes.stat <- cbind(DAS.pval.simes,maxdeltaPS)
    # colnames(DAS.simes.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    
    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3)))
  message('Done!!! ')
  return(results)
}



#############################################################################################
##-----------------------------------------------------------------------------------------
#' @rdname limma.pipeline
#' @export
edgeR.pipeline <- function(dge,
                           method=c('glm','glmQL'),
                           design,
                           contrast,
                           diffAS=F,
                           deltaPS=NULL,
                           adjust.method='BH'){
  
  start.time <- Sys.time()
  results <- list()
  method <- match.arg(method,c('glm','glmQL'))
  targets <- rownames(dge$counts)
  message('Estimate dispersion ...')
  Disp <- estimateGLMCommonDisp(dge, design)
  Disp <- estimateGLMTagwiseDisp(Disp, design)
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  
  switch(method,
         glm = {
           message('Fit Genewise Negative Binomial Generalized Linear Models...')
           fit <- glmFit(Disp, design = design)
           results$fit.glm <- fit
           
           ###individual test
           message('Fit the contrast model ...')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmLRT(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })
           ###overall test
           x <- glmLRT(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         },
         glmQL = {
           message('Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests...')
           fit <- glmQLFit(Disp,design = design)
           results$fit.glmQL <- fit
           
           ###individual test
           message('Fit the contrast model ...')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmQLFTest(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })
           
           ###overall test
           x <- glmQLFTest(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         }
  )
  names(DE.pval.list) <- contrast
  results$DE.pval.list <- DE.pval.list
  
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$FDR))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast
  
  DE.stat <- summaryStat (x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  ##################################################################
  ###---DAS
  if(diffAS){
    fit.splice <- lapply(contrast,function(i){
      message(paste0('\ndiffSpliceDGE of contrast: ', i))
      diffSpliceDGE(fit, contrast=contrast.matrix[,i,drop=F], geneid="GENEID")
    })
    names(fit.splice) <- contrast
    results$fit.splice<-fit.splice
    
    genes.idx <- unique(fit.splice[[1]]$genes$GENEID)
    trans.idx <- unique(fit.splice[[1]]$genes$TXNAME)
    
    DTU.pval.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="exon", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast
    
    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast
    
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]
    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))
    
    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat
    
    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice[[1]]$genes$TXNAME,,drop=F],
                     INDICES = fit.splice[[1]]$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS
    
    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="gene", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.F.list) <- contrast
    
    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast
    
    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat
    
    ##---Simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="Simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.simes.list) <- contrast
    
    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast
    
    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))
    
    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3)))
  message('Done!!! ')
  return(results)
}



##########################################################################
#' Generate all possible logical relations between set x and y.
#' @param x,y vectors of set x and set y.
#' @return a list of three elements: "x.only", "xy" and "y.only".
#' @examples 
#' x <- letters[1:10]
#' y <- letters[5:15]
#' z <- set2(x,y)
#' combo <- c(x=length(z$x.only),y=length(z$y.only),"x&y"=length(z$xy))
#' plotEulerDiagram(combo)
#' 
#' @export
#' @seealso \code{\link{plotEulerDiagram}} and \code{\link{eulerr::euler}}
#' 
set2 <- function(x,y){
  x.only <- setdiff(x,y)
  xy <- intersect(x,y)
  y.only <- setdiff(y,x)
  results <- list(x.only=x.only,xy=xy,y.only=y.only)
  attributes(results) <- list(
    x.only='x-x&y',
    xy='x&y',
    y.only='y-x&y'
  )
  names(results) <- c('x.only','xy','y.only')
  return(results)
}



##########################################################################
#' Generate html, pdf and word report from Rmd (RMardkown) file
#' @param input.Rmd the Rmd file used to generate the report.
#' @param report.folder the folder to save the reports.
#' @param type the ouput document type. Options are "all", "html_document", "pdf_document" and "word_document".
#' @return the generated html, pdf or word report will be saved into the provided \code{report.folder}.
#' @export
#' 
generate.report <- function(input.Rmd,report.folder,
                            type=c('all','html_document','pdf_document','word_document')){
  if(!file.exists(report.folder))
    dir.create(report.folder,recursive = T)
  type <- match.arg(type,c('all','html_document','pdf_document','word_document'))
  if(type=='all'){
    try(rmarkdown::render(input = input.Rmd,
                          params = list(doc.type = "html_document"
                          ),
                          output_format = "html_document",
                          output_dir = report.folder))
    try(rmarkdown::render(input = input.Rmd, 
                          params = list(doc.type = "pdf_document"
                          ),
                          output_format = "pdf_document",output_dir = report.folder))
    
    try(rmarkdown::render(input = input.Rmd, 
                          params = list(doc.type = "word_document"
                          ),
                          output_format = "word_document",output_dir = report.folder))
  } else {
    try(rmarkdown::render(input = input.Rmd, 
                          params = list(doc.type = type
                          ),
                          output_format = type,output_dir = report.folder))
  }
}

##########################################################################
#' Filter low expression based on CPM (count per million reads)
#' @details An expressed target must have \eqn{\geq} \code{sample.n} with CPM \eqn{\geq} \code{cpm.cut}.
#' @param counts a matrix/data.frame or matrix of read counts.
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
#' @details An expressed target must have \eqn{\geq} \code{sample.n} with TPM \eqn{\geq} \code{tpm.cut}.
#' @param TPM a data.frame/matrix of TPM.
#' @param sample.n number of samples
#' @param tpm.cut a numeric value of TPM cut-off
#' @return A vector of logcial values (TRUE/FALSE) to indicate the rows of the input data to keep/filter.
#' @export
TPM.filter<-function(TPM,sample.n=3,tpm.cut=1){
  rowSums(TPM>=tpm.cut)>=sample.n
}

#' Filter low expressed transcripts and genes
#' @details If read counts are used, they are converted to CPM (count per million reads) to filter low expressed transcripts. 
#' 
#' 1) An expressed transcript must have \eqn{\geq} \code{sample.n} with CPM \eqn{\geq} \code{cpm.cut}.
#' 
#' 2) An expressed gene must have at least one expressed transcript.
#' 
#' @param abundance a data.frame/matrix of read counts/TPMs at transcript level.
#' @param mapping a data.frame of transcript-gene name mapping, with first column of transcript list and second column of gene list.
#' @param abundance.cut a numeric value of abundance cut-off. 
#' @param sample.n a number of samples cut across all samples.
#' @param Log logical. If TRUE, the log2-CPM is used for filters when read counts are used as abundance.
#' @param unit the abundance used to filter the low expressed genes/transcripts. Options are "counts" and "TPM".
#' @return  A list of following elements:
#'     \itemize{
#'       \item{ trans_high: }{ the expressed transcripts}
#'       \item{ genes_high: }{ the expressed genes}
#'       \item{ mapping: }{ the remaining mapping data.frame after filter the low expressed transcripts}
#'     }
#' @export
#' 
low.expression.filter <- function(abundance,
                                  mapping,
                                  abundance.cut=1,
                                  sample.n=3,
                                  Log=F,
                                  unit=c('counts','TPM')){
  unit <- match.arg(unit,c('counts','TPM'))
  colnames(mapping) <- c('TXNAME','GENEID')
  rownames(mapping) <- mapping$TXNAME
  if(unit == 'counts'){
    filter.idx<-CPM.filter(counts = abundance,sample.n = sample.n,cpm.cut = abundance.cut,Log = Log)
  } else {
    filter.idx<-TPM.filter(TPM = abundance,sample.n = sample.n,tpm.cut = abundance.cut)
  }
  trans_high <- names(filter.idx)[filter.idx==T]
  genes_high <- unique(mapping[trans_high,]$GENEID)
  list(trans_high=trans_high,genes_high=genes_high,mapping_high=mapping[trans_high,])
}

##########################################################################
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

##########################################################################
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

##########################################################################
#' Calculate standard deviation or standard errors
#'
#' Calculate standard deviation or standard errors from a numerical vector.
#'
#' @param x a numerical vector
#' @param error.type method used to calculate data error. Options include "sd" for standard deviation and "stderr" for standard error.
#'
#' @return data error
#' @seealso \code{\link{sd}}
#'
#' @examples
#'
#' set.seed(2000)
#' x=rnorm(1000,mean=0,sd=1)
#' ##standard error
#' data.error(x,error.type = 'stderr')
#' ##standard deviation
#' data.error(x,error.type = 'sd')
#'
#' @export
#'

data.error <- function (x, error.type = c("stderr","sd")){
  error.type <- match.arg(error.type,c("stderr","sd"))
  if (error.type == "stderr") 
    error <- sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))
  if (error.type == "sd") 
    error <- sd(na.omit(x))
  return(error)
}

##########################################################################
#' Generate HCL colours, which are the same to \code{ggplot} default colour palette
#' @param n number of colours to generate.
#' @return a vector of \code{n} colours.
#' @export
gg.color.hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate distinct colours
#' @param n number of colours to generate.
#' @return a vector of \code{n} colours.
distinct.color <- function(n){
  col.lib <- c("#C0D57E","#BF6CF4","#D435BB","#C2BBCB","#CAE678","#BEC3AD","#EC75F6","#F8DAD6","#D9FBD7","#46C3CF","#5FECDE","#AEE1F5","#AA67CC",
               "#81C1F1","#D6A381","#81DADB","#99AD6D","#B8ECF3","#828D6F","#F832A9","#8DD092","#7FC0D0","#6EDDD0","#D4BE95","#B8F7F2","#9752C4",
               "#AF9D7C","#468CE2","#9AA843","#4F9F7A","#8C7DB6","#CA44D8","#93D0B9","#4A718E","#F259C8","#E1C358","#BB924E","#5682F1","#B3E18E",
               "#59C1A6","#3BA7EB","#4ED06A","#A28FF0","#F4F2E7","#F7C2DA","#6921ED","#E8AAE6","#C2F2DF","#A3F7A1","#F0A85C","#BFF399","#E93B87",
               "#D7C7F0","#BCF7B6","#39B778","#DD714B","#B7D9A6","#E19892","#BF70C4","#B882F4","#74F38B","#D6E3C2","#9E23B2","#65F55F","#97F3EB",
               "#6BCD37","#8A7B8E","#DEEDA7","#F03E4F","#8FE77B","#ACF165","#8FCAC4","#CA1FF5","#E854EE","#43E11C","#6EEDA7","#B584B3","#BAAAF2",
               "#DE88C9","#415B96","#4EF4F3","#81F8CA","#C1B4AC","#BBA4CF","#AB55E5","#F789EC","#CEDFF8","#A9BEF0","#658DC6","#94ACEA","#91CD6F",
               "#BAA431","#8FB2C5","#F8F1A3","#D1CCCE","#F533D8","#844C56","#F3DA7D","#BDE5C0","#F1F69C","#F17C85","#C5E743","#B8C8ED","#3ED3EC",
               "#E6C5D9","#AD7967","#EAF4C3","#3F4DC1","#F4DDEC","#F1B136","#C247F3","#6FD3A1","#DCF7A3","#B8D1C2","#545CED","#2C7A29","#E4B1BB",
               "#785095","#8373F1","#D0F2C5","#B53757","#49F79F","#F0BE8B","#9630DD","#D6D89F","#79EED3","#EEEFD7","#7FA41F","#D1BF71","#6AEFF6",
               "#D661D8","#CFEEF8","#9ECB59","#67D382","#A8945E","#BB9B9E","#F3E4B9","#9964F3","#835BD4","#A5E2BF","#9CA6C8","#CD9BE9","#8FF038",
               "#B6C08A","#D19BB5","#8893EB","#51746E","#5920AC","#F3A1BB","#D7E1DA","#D668C3","#EC9EF2","#C2738A","#EDE979","#615AB4","#E9EF68",
               "#F3DECA","#90A6A0","#F2B8E7","#36DEA3","#EED89A","#5AADA5","#92EFAD","#A3288A","#ECC9F4","#91C9E7","#ECC2B8","#EB9BCC","#2735BB",
               "#5E7A2D","#C3BDF7","#DCF8F4","#52E0BC","#C2F8D6","#B3CFDA","#6E93A2","#9AECCA","#E3D8F0","#D97EE7","#9CBE9A","#60A858","#55F9D7",
               "#F47AA2","#DCF5E3","#B55398","#AC81D3","#644AF2","#546DEA","#E7EB37","#EA7E2A","#A1517E","#E131E4","#DE92F4","#EC66B1","#82DDF3",
               "#E99475","#E4EAF3","#E8D144","#40C0EE","#559ABA")
  if(n>length(col.lib))
    stop(paste0('n should be smaller than ',length(col.lib)))
  col.lib[1:n]
}

# gg_legend <- function(g){ 
#   tmp <- ggplot_gtable(ggplot_build(g)) 
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#   legend <- tmp$grobs[[leg]] 
#   return(legend)
#   } 


##########################################################################
reorder.clusters <- function(clusters,dat){
  targets <- rownames(dat)
  names(clusters) <- targets
  idx0 <- rowmean(dat,group = clusters,reorder = T)
  para <- lapply(1:ncol(idx0),function(i) idx0[,i])
  idx0 <- idx0[do.call(order,para),]
  idx0 <- idx0[rev(rownames(idx0)),]
  order.idx <- rownames(idx0)
  idx <- 1:nrow(idx0)
  names(idx) <- order.idx
  idx <- idx[as.character(clusters)]
  names(idx) <- names(clusters)
  idx
}

##########################################################################
#' Calculate column mean of a matrix or data frame based on a grouping variable
#'
#' Compute column means across rows of a numeric matrix-like object for each level
#' of a grouping variable.
#' @param x a matrix or data frame.
#' @param group a vector of factor giving grouping, with one element per row of x.
#' @param reorder if TRUE, then the result will be in order of \code{sort(unique(group))}.
#' @param na.rm logical (TRUE or FALSE). Should NA (including NaN) values be replaced by value 0?
#' @return \code{rowmean} returns a matrix or data frame containing the means. There
#' will be one row per unique value of group.
#'
#' @examples
#' x <- matrix(runif(50), ncol = 5)
#' group <- sample(1:4, 10, TRUE)
#' xmean <- rowmean(x, group)
#'
#' @export
#' @seealso \code{\link{rowsum}}, \code{\link{rowratio}}
rowmean<-function(x,group,reorder=F,na.rm=T){
  order.idx<-as.character(unique(group))
  if (reorder)
    order.idx<-gtools::mixedsort(order.idx)
  
  counts <- table(group)[order.idx]
  sums <- rowsum(x, group = group)[order.idx,]
  means <- (diag(1/counts)) %*% as.matrix(sums)
  rownames(means) <- order.idx
  if (na.rm)
    means[is.na(means)] <- 0
  return(means)
}


##########################################################################
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
##########################################################################
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
#' transAbundance <- matrix(rnorm(36),ncol = 9,nrow = 3)
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
  rownames(transAbundance) <- mapping$TXNAME
  
  if(is.null(PS)){
    transAbundance.mean <- t(rowmean(t(transAbundance),group = condition,reorder = F))
    ##---PS
    PS <- rowratio(x = transAbundance.mean[as.vector(mapping$TXNAME),],
                   group = as.vector(mapping$GENEID), reorder = F)
    PS <- data.frame(PS[mapping$TXNAME,])
  }
  
  ##---deltaPSI
  deltaPS <- lapply(strsplit(contrast,'-'),function(x){
    PS[,x[1]]-PS[,x[2]]
  })
  deltaPS <- do.call(cbind,deltaPS)
  colnames(deltaPS) <- contrast
  rownames(PS) <- rownames(deltaPS) <- mapping$TXNAME
  list(PS=PS,deltaPS=deltaPS)
}

##########################################################################
##########################################################################



##dashboardHeader######################################################
# ========================== dashboardHeader ======================== #
mainheader <- dashboardHeader(
  title = '3D RNA-seq App',
  titleWidth = 200,
  dropdownMenuOutput("messageMenu")
)

##dashboardSidebar#####################################################
# ========================= dashboardSidebar ======================== #
mainsidebar <- dashboardSidebar(
  width = 200,
  sidebarMenu(
    menuItem(text = 'Introduction',tabName = 'introduction',icon = icon("book"),selected = T),
    menuItem(text = 'Data generation',tabName = 'generation',icon = icon("database"),selected = F),
    menuItem(text = 'Data pre-processing',tabName = 'preprocessing',icon = icon("cogs"),startExpanded = F,selected = F
    ),
    menuItem(text = 'DE DAS and DTU',tabName = 'ddd',icon = icon("random"),selected = F),
    menuItem(text = 'Result summary',tabName = 'resultssummary',icon = icon("sitemap",lib='font-awesome'),selected = F),
    menuItem(text = 'Advanced plot',tabName = 'advancedanalysis',icon = icon("line-chart",lib='font-awesome'),selected = F),
    menuItem(text = 'Generate report',tabName = 'report',icon = icon("file",lib='font-awesome'),selected = F)
  )
)

##dashboardBody########################################################
# ========================== dashboardBody ========================== #
mainbody <- dashboardBody(
  tags$head(includeScript(system.file("google-analytics.js", package = "ThreeDRNAseq"))),
  withMathJax(),
  shinyjs::useShinyjs(),
  tabItems(
    tabItem("introduction",
            fluidRow(
              column(width = 12,
                     box(title='Introduction',
                         width = NULL,status = 'primary', solidHeader = T,
                         # shiny::uiOutput('page1'),
                         # DT::dataTableOutput("info")
                         HTML('
                              <center><img src="https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/fig/header.gif?raw=true" align="middle" style="height: 45%; width: 45%; object-fit: contain;"></center>
                              <h3>Description</h3>
                              <div align="justify">
                              <p><strong>ThreeDRNAseq (3D RNA-seq)</strong> R package provides an interactive graphical user interface (GUI) for RNA-seq data differential expression (DE), differential alternative splicing (DAS) and differential transcript usage (DTU) (3D) analyses based on two popular pipelines: <a href="https://bioconductor.org/packages/release/bioc/html/limma.html" target="_blank">limma</a> (Smyth et al. 2013) and <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">edgeR</a> (Robinson et al., 2010). The 3D RNA-seq GUI is based on R <a href="https://shiny.rstudio.com/" target="_blank">shiny</a> App and enables the RNA-seq analysis to be done only with in <strong>3 Days (3D)</strong>. To perform the 3D analysis,</p>
                              <ul>
                              <li>The first step is to generate transcript quantification from tools, such as <a href="https://combine-lab.github.io/salmon/" target="_blank">Salmon</a> (Patro et al., 2017) and <a href="https://pachterlab.github.io/kallisto/" target="_blank">Kallisto</a> (Bray et al., 2016). The transcript quantification procedure often takes <strong>1-2 Days</strong>.</li>
                              <li>Then users can do mouse click on the App to upload transcript read counts, perform 3D analysis, and make beautiful plots, e.g. expression mean-variance trend plots, PCA plots, heatmap, GO annotation plots, etc.. The results of significance and all the plots can be wrapped into html, pdf and/or word reports by one-button-click. The entire 3D mouse-click analysis only takes <strong>1 Day or even less</strong>.</li>
                              </ul>
                              <p>The 3D RNA-seq analysis pipeline has steps of proper data pre-processing, e.g. low expression filters based on expression mean-variance trend, visualise data variation, batch effect estimation, data normalisation, etc.. The optimal parameters for each step are determined by quality control plots, e.g. mean-variance trend plots, PCA plots, data distribution plots, etc.. The pipeline also has stringent controls of false positive, leading to robust 3D predictions. The pipeline has been successfully applied in different RNA-seq studies from Arabidopsis (Calixto et al., 2018), barley (Bull et al., 2017) and potato to identify novel condition-responsive genes/transcripts, especially those with significant alternative splicing changes. To use our pipeline in your work, please cite:</p>
                              <p><a href="http://www.plantcell.org/content/30/7/1424" target="_blank">Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.</a></p>
                              <h3 id="user-manuals">User manuals</h3>
                              <p>Two versions of step-by-step user manuals are provided:</p>
                              <ul>
                              <li>3D RNA-seq App manual (command-line free) for easy to use: <a href="https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md" target="_blank">https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md</a></li>
                              <li>3D RNA-seq command-line based manual for advanced R users: <a href="https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/Command_line_based_manul.md" target="_blank">https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/Command_line_based_manul.md</a></li>
                              </ul>
                              <h3 id="references">References</h3>
                              <ul>
                              <li><p>Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525–527.</p></li>
                              <li><p>Bull,H., Casao,M.C., Zwirek,M., Flavell,A.J., Thomas,W.T.B.B., Guo,W., Zhang,R., Rapazote-Flores,P., Kyriakidis,S., Russell,J., Druka,A., McKim,S.M., and Waugh,R. (2017) Barley SIX-ROWED SPIKE3 encodes a putative Jumonji C-type H3K9me2/me3 demethylase that represses lateral spikelet fertility. Nat. Commun., 8, 936.</p></li>
                              <li><p>Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell, tpc.00177.2018.</p></li>
                              <li><p>Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.</p></li>
                              <li><p>Robinson,M.D., McCarthy,D.J., and Smyth,G.K. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–40.</p></li>
                              <li>Smyth,G.K., Ritchie,M., Thorne,N., Wettenhall,J., and Shi,W. (2013) limma:Linear Models for Microarray Data User’s Guide(Now Including RNA-Seq Data Analysis). R Man., 1–123.</li>
                              </ul>
                              </div>
                              ')
                         )
                         )
                         )
            
                         ),
    
    #=======================>> Data generation panel <<============================
    tabItem("generation",
            tags$head(tags$style("#TxtOut {white-space: normal;}")),
            fluidRow(
              column(width = 12,
                     box(title = 'Choose working directory',
                         width=13,status = 'primary', solidHeader = T,
                         # HTML('<h4><strong>Choose working directory</strong></h4>'),
                         # wellPanel(
                         # column(width = 2,
                         shinyDirButton("wdir", "Chose directory", "Upload",
                                        buttonType ="primary"),
                         br(),
                         # ),
                         # column(width = 10,
                         verbatimTextOutput("wdir.info"),
                         # ),
                         actionButton(inputId = 'create.folders',label = 'Create folders to save results',
                                      icon("folder",lib="font-awesome"),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                         br(),
                         br(),
                         HTML('If "data", "result", "figure" and "report" folders are not in working directory,
                              please click this button to create.')
                         )
                     )
              ),
            fluidRow(
              column(width = 12,
                     box(width = 13,title = 'Generate/Load data',
                         status = 'primary',solidHeader = T,
                         # shiny::dataTableOutput('loaded.data.info')
                         verbatimTextOutput("DDD.data.size"),
                         tableOutput('loaded.data.info'),
                         h4('Load data for analysis'),
                         div(style="display:inline-block;vertical-align:middle",
                             actionButton('loaded.data.all','Load all',icon("upload",lib="font-awesome"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                             h4('OR Load data for individual steps'),
                             actionButton('loaded.data.preprocessing','Load Data pre-processing',icon("upload",lib="font-awesome"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                             actionButton('loaded.data.ddd','Load DE DAS and DTU',icon("upload",lib="font-awesome"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                             actionButton('loaded.data.summary','Load Result summary',icon("upload",lib="font-awesome"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                             actionButton('loaded.data.advanced','Load Advanced plot',icon("upload",lib="font-awesome"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         ),
                         hr(),
                         HTML('<p>The whole pipeline has a number of steps and each step is related to several intermediate datasets. This table summarise the generation and usage steps of these datasets and whether (True or False) they have been loaded in workspace ready for use. To avoid running the pipeline from very beginning every time the App get refreshed, the following options are provided:</p>
                              <ul>
                              <li>If all/part of the required datasets (in .RData format) do not exist in the “data” folder of working directory, users can generate these datasets from corresponding steps.</li>
                              <li>If all/part of the required datasets (in .RData format) already exist in the “data” folder, users can click “Load all” button to load them.</li>
                              <li>To perform analysis or visualise results in a specific step, users also can load the required datasets for individual steps.</li>
                              </ul>
                              ')
                         ))
                         ),
            fluidRow(
              #load sample tables####
              column(width = 6,
                     box(title='Step 1: Load quantification and experimental information',
                         width = 13,status = 'primary', solidHeader = T,
                         # HTML('<h4><strong>Step 1: load gene-transcript mapping</strong></h4>'),
                         fileInput("mapping.input", "Choose mapping.csv file",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                         ),
                         br(),
                         HTML('This file includes the gene-transcript mapping information. See examples in the manual.'),
                         hr(),
                         fileInput("sample.input", "Choose samples.csv file",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                         ),
                         br(),
                         HTML('This file includes the information of conditions for linear regression,
                              biological replicates, sequencing replicates and the complete paths of
                              transcript quantification.
                              See details in the manual.')
                         )
                         ),
              #generate gene expression####
              # #generate transcript expression####
              column(width = 6,
                     box(title='Step 2: Generate gene and transcript expression',
                         width = 13,status = 'primary', solidHeader = T,
                         HTML('<h4><strong>Generate gene expression</strong></h4>'),
                         selectInput(inputId = "tximport.quant.method",
                                     label = "Abundance generated from:",
                                     choices = c("salmon","kallisto","sailfish","rsem","stringtie","none"),
                                     selected = "salmon",
                                     width = 200
                                     
                         ),
                         div(style="display:inline-block;vertical-align:middle",
                             selectInput(inputId = "tximport.method",
                                         label = "Choose a tximport method:",
                                         choices = c("lengthScaledTPM", "scaledTPM", "no"),
                                         selected = "lengthScaledTPM",
                                         width = 200
                             )
                         ),
                         div(style="display:inline-block;vertical-align:middle;height:30px",
                             actionButton('run.txi.genes','Run',icon("send outline icon"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         ),
                         radioButtons('generate.new.genes.data',label = NULL,inline = T,
                                      choices = c('Generate new data','Load txi_genes.RData'),
                                      selected = 'Generate new data'
                         ),
                         HTML('If the <i>txi_genes.RData</i> file is already in the "data" folder, to save time, user can select "Load txi_genes.RData"
                              and click "Load" to load the data directly.'
                         ),
                         br(),
                         verbatimTextOutput("rungeneinfo"),
                         hr(),
                         HTML('<h4><strong>Generate transcript expression</strong></h4>'),
                         div(style="display:inline-block;vertical-align:middle;height:30px",
                             actionButton('run.txi.trans','Run',icon("send outline icon"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4"
                             )
                         ),
                         br(),
                         br(),
                         radioButtons('generate.new.trans.data',label = NULL,inline = T,
                                      choices = c('Generate new data','Load txi_trans.RData'),
                                      selected = 'Generate new data'
                                      
                         ),
                         HTML('If the <i>txi_trans.RData</i> file is already in the "data" folder, to save time, user can select "Load txi_trans.RData"
                              and click "Run" to load the data directly.'
                         ),
                         br(),
                         verbatimTextOutput("runtransinfo")
                         )
                     )
                         ),
            # output sample table####
            fluidRow(
              tabBox(width = 13,
                     # Title can include an icon
                     # title = tagList(shiny::icon("gear"), "Loaded tables"),
                     title = "Loaded tables",side = 'right',
                     tabPanel("Mapping",
                              "Part of gene-transcript mapping table:",
                              DT::dataTableOutput("mapping.output")
                     ),
                     tabPanel("Samples",
                              "The sample information:",
                              DT::dataTableOutput("sample.output")
                     )
              )
            )
                     ),
    
    #========================>> Data pre-processing panel <<=======================
    tabItem('preprocessing',
            ##----------Step 1: merge sequencing reps------------
            fluidRow(
              # column(width = 3,
              box(title = 'Step 1: Merge sequencing replicates',
                  width=6,status = 'primary', solidHeader = T,
                  verbatimTextOutput("seqrep.info"),
                  radioButtons(inputId = 'merge.or.skip',label = 'Merge options',
                               choices = c('Merge','Skip'),inline = T),
                  actionButton('run.seq.merge','Merge',icon("plus-square"),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                  hr(),
                  HTML('If there are no sequencing replicates, please select "Skip".')
              ),
              box(title = 'Step 2: Filter low expression',
                  width = 6,status = 'primary',solidHeader = T,
                  sliderInput(inputId = 'cpm.cut',label = 'CPM cut-off',
                              min = 0,max = 10,value = 1),
                  sliderInput(inputId = 'sample.n.cut',label = 'Sample number cut-off',
                              min = 0,max = 10,value = 1),
                  actionButton('run.filter','Filter',icon("filter"),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                  actionButton('mv.plot.button','Mean-variance plot',icon("pencil"),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                  br(),
                  br(),
                  verbatimTextOutput("low.filter.info"),
                  hr(),
                  HTML('Note: Click "Filter" and then click "Mean-variance plot".
                       <ul><li>An expressed transcript must have &#x2265 <i>n</i> samples &#x2265 <i>m</i> CPM (Count Per Million reads) expression.</li>
                       <li>An expressed gene must have at least one expressed transcript.</li></ul>')
                  )
                  ),
            ##----------Step 2: Filter low expression------------
            fluidRow(
              box(title = NULL,
                  width=12,status = 'primary', solidHeader = T,
                  HTML('RNA-seq data information: target numbers'),
                  DT::dataTableOutput("data.info"),
                  # tableOutput("data.info"),
                  br(),
                  actionButton(inputId = 'save.data.info',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
              ),
              ##---------- Mean-variance trend plot------------
              column(width = 12,
                     tabBox(title = 'Mean-variance trend plot',id = 'mv.plot', width = 13,
                            tabPanel(title = 'Transcript level',
                                     plotOutput(outputId = 'mv.trans.plot',click ='mv.trans.plot.click')
                            ),
                            tabPanel(title = 'Gene level',
                                     plotOutput(outputId = 'mv.genes.plot',click ='mv.genes.plot.click')
                            )
                     ),
                     actionButton(inputId = 'save.mv.plot',label = 'Save',
                                  icon = icon('download',lib = 'font-awesome'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px;"),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "mv.plot.res",label = 'PNG plot resolution',value = '150',
                                      min = 0,max = 600,step = 60,width = "100%")
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "mv.plot.height",label = 'Plot height (inch)',
                                      value = 4.5,width = "100%",step = 0.2)
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "mv.plot.width",label = 'Plot width (inch)',
                                      value = 11.5,width = "100%",step = 0.2)
                     )
              )
            ),
            br(),
            br(),
            br(),
            ##----------Step 3: Principal component analysis------------
            fluidRow(
              column(width = 3,
                     box(title = 'Step 3: Principal component analysis (PCA)',width = 13,
                         status = 'primary',solidHeader = T,
                         div(style="display:inline-block;vertical-align:middle;width:150px",
                             selectInput(inputId = 'dim1',label = 'X-axis',
                                         choices = c('PC1','PC2','PC3','PC4'),selected = 'PC1')
                         ),
                         div(style="display:inline-block;vertical-align:middle;width:150px",
                             selectInput(inputId = 'dim2',label = 'Y-axis',
                                         choices = c('PC1','PC2','PC3','PC4'), selected = 'PC2')
                         ),
                         div(style="display:inline-block;vertical-align:middle;height:30px",
                             actionButton('plot.pca.button','Plot',icon("pencil",lib = 'font-awesome'),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         ),
                         br(),
                         br(),
                         radioButtons(inputId = 'pca.plot.type',label = 'Plot type',
                                      choices = c('Bio-reps','Average expression'),
                                      selected = 'Bio-reps',inline = T),
                         # radioButtons(inputId = 'pca.color.option',label = 'Colour on',
                         #              choices = c('Bio-reps','Conditions'),
                         #              selected = 'Bio-reps',inline = T),
                         uiOutput('pca.color.option.ui'),
                         radioButtons(inputId = 'pca.add.circle',label = 'Add polygon to clusters',
                                      choices = c('none','ellipse','polygon'),
                                      selected = 'none',inline = T),
                         HTML('<ol>
                              <li style="margin-left: -20px;">Select PCs to visualise</li>
                              <li style="margin-left: -20px;">Select "Plot type", either show all Bio-reps or average expression of conditions</li>
                              <li style="margin-left: -20px;">Select colour according to columns in the samples.csv (multi-select allowed)</li>
                              <li style="margin-left: -20px;">Select "ellipse" or "polygon" to highlight the clusters</li>
                              </ol>')
                         )
                         ),
              column(width = 9,
                     tabBox(title = 'PCA plot',id = 'pca.plot',width = 13,
                            tabPanel(title = 'Transcript level',
                                     plotlyOutput('trans.pca.plot')
                            ),
                            tabPanel(title = 'Gene level',
                                     plotlyOutput('genes.pca.plot')
                            )
                     ),
                     actionButton(inputId = 'save.pca.plot',label = 'Save',
                                  icon = icon('download',lib = 'font-awesome'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "pca.plot.res",label = 'PNG plot resolution',value = '150',
                                      min = 0,max = 600,step = 60,width = "100%")
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "pca.plot.height",label = 'Plot height (inch)',
                                      value = 8,width = "100%",step = 0.2)
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "pca.plot.width",label = 'Plot width (inch)',
                                      value = 9,width = "100%",step = 0.2)
                     ),
                     br(),
                     br(),
                     br(),
                     br()
              )
                         ),
            fluidRow(column(width = 3,
                            box(title = 'Step 3-continued: Batch effect estimation',width = 13,
                                status = 'primary',solidHeader = T,
                                HTML('<h4><strong>If no distinct batch effects between bio-reps, please skip this step.</strong></h4>'),
                                hr(),
                                radioButtons(inputId = 'ruvseq.method',label = 'RUVSeq method',
                                             choices = c('RUVr','RUVs','RUVg'),selected = 'RUVr',
                                             inline = T),
                                sliderInput(inputId = 'ruvseq.k',label = 'Number of possible factors caused batch effects',
                                            min = 1,max = 6,value = 1,step = 1),
                                actionButton('remove.batch','Remove batch effects',icon("cut"),
                                             style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                                br(),
                                br(),
                                actionButton('update.pca.plot','Update PCA plot',icon("refresh",lib="font-awesome"),
                                             style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                                br(),
                                br(),
                                verbatimTextOutput("br.saved.info"),
                                br(),
                                HTML('The parameters to plot PCA are inherited from previous PCA panel.')
                            )
            ),
            column(width = 9,
                   tabBox(title = 'PCA plot after remove batch effects',id = 'pca.plot.batch.removed',width = 13,
                          tabPanel(title = 'Transcript level',
                                   plotlyOutput('trans.pca.br.plot')
                          ),
                          tabPanel(title = 'Gene level',
                                   plotlyOutput('genes.pca.br.plot')
                          )
                   ),
                   actionButton(inputId = 'save.pca.br.plot',label = 'Save',
                                icon = icon('download',lib = 'font-awesome'),
                                style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                   br(),
                   br(),
                   br(),
                   br()
            )
            ),
            
            ##----------Step 4: Normalization------------
            fluidRow(
              column(width = 3,
                     box(title = 'Step 4: Normalization',width = 13,
                         status = 'primary',solidHeader = T,
                         radioButtons(inputId = 'norm.method',label = 'Normalization methods',
                                      choices = c('TMM','RLE','upperquartile'),
                                      selected = 'TMM',inline = T),
                         actionButton('run.norm','Run',icon("send outline icon",lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                         actionButton('norm.plot.button','Plot distribution',icon("pencil"),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                         br(),
                         br(),
                         verbatimTextOutput("norm.info")
                     )
                     
              ),
              column(width = 9,
                     tabBox(title = 'Data distribution plot',id = 'dist.plot',width = 13,
                            tabPanel(title = 'Transcript level',
                                     plotlyOutput('trans.dist.plot',height = 850)
                            ),
                            tabPanel(title = 'Gene level',
                                     plotlyOutput('genes.dist.plot',height = 850)
                            )
                     ),
                     actionButton(inputId = 'save.dist.plot',label = 'Save',
                                  icon = icon('download',lib = 'font-awesome'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "dist.plot.res",label = 'PNG plot resolution',value = '150',
                                      min = 0,max = 600,step = 60,width = "100%")
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "dist.plot.height",label = 'Plot height (inch)',
                                      value = 7,width = "100%",step = 0.2)
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "dist.plot.width",label = 'Plot width (inch)',
                                      value = 8,width = "100%",step = 0.2)
                     ),
                     br(),
                     br(),
                     br(),
                     br()
              )
            )
              ),
    #========================>> 3D panel <<=======================
    tabItem('ddd',
            ##---------------Contrast group------------
            fluidRow(
              column(width = 4,
                     box(title='Step 1: Load contrast group',
                         width = 13,status = 'primary', solidHeader = T,
                         # HTML('<h4><strong>Step 1: load gene-transcript mapping</strong></h4>'),
                         fileInput("contrast.input", "Choose contrast.csv file",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                         ),
                         hr(),
                         HTML('<h4><strong>Contrast groups</strong></h4>'),
                         DT::dataTableOutput('contrast.table'),
                         br(),
                         HTML('The expression of "Treatment" conditions is compared to "Control" conditions.')
                     )
              ),
              column(width = 4,
                     box(title='Step 2: DE genes',
                         width = 13,status = 'primary', solidHeader = T,
                         HTML('<h4><strong>Choose a pipeline</strong></h4>'),
                         div(style="display:inline-block;vertical-align:middle;height:30px",
                             selectInput(inputId = "DE.pipeline", label = NULL,
                                         choices = list(`limma-voom:`=list('limma'),`edgeR:` = c("glmQL", "glm")),
                                         selected = 'limma',width = 160
                             )
                         ),
                         div(style="display:inline-block;vertical-align:middle;height:30px",
                             actionButton('run.DE','Run',icon("send outline icon"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         ),
                         br(),
                         br(),
                         HTML('<p align="justify">In practics, limma-voom and edgeR-glmQL have comparable performance, and both have
                              better control of false positives than edgeR-glm.</p>'),
                         hr(),
                         HTML('<h4><strong>Set thresholds</strong></h4>'),
                         selectInput(inputId = 'p.adjust.method',label = 'P-value adjust method',
                                     choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                 "fdr", "none"),selected = 'BH'),
                         sliderInput(inputId = 'pval.cutoff',label = 'Adjusted p-value',
                                     min = 0,max = 0.05,value = 0.01,step = 0.005),
                         sliderInput(inputId = 'lfc.cutoff',label = 'Absolute \\(\\log_2\\)FC',
                                     min = 0,max = 5,value = 1,step = 0.5),
                         verbatimTextOutput("run.DE.info")
                         )
              ),
              column(width = 4,
                     box(title='Step 3: DAS genes, DE and DTU transcripts',
                         width = 13,status = 'primary', solidHeader = T,
                         HTML('<h4><strong>Generate transcript deltaPS</strong></h4>'),
                         div(style="display:inline-block;vertical-align:middle",
                             selectInput(inputId = "PS.data.type",
                                         label = "Data for deltaPS:",
                                         choices = c("TPM", "Read counts"),
                                         selected = "TPM",
                                         width = 160
                             )
                         ),
                         div(style="display:inline-block;vertical-align:middle;height:30px",
                             actionButton('run.deltaPS','Run',icon("send outline icon"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         ),
                         br(),
                         verbatimTextOutput("rundeltaPSinfo"),
                         hr(),
                         HTML('<h4><strong>Set threshods</strong></h4>'),
                         selectInput(inputId = 'DAS.p.method',label = 'DAS gene test method',
                                     choices = c('F-test','Simes'),selected = 'F-test'),
                         sliderInput(inputId = 'deltaPS.cutoff',label = 'Absolute \\(\\Delta\\)PS',
                                     min = 0,max = 1,value = 0.1,step = 0.05),
                         br(),
                         HTML('<p align="justify">PS (percent of splice) is defined as the ratio of transcript average abundance of conditions divided by the average
                              gene abundance. \\(\\Delta\\)PS is the PS differences of conditions based on the contrast groups. The abundance to calculate
                              PS values can be TPMs or read counts.</p>'),
                         hr(),
                         HTML('<h4><strong>Proceed the analysis</strong></h4>'),
                         div(style="display:inline-block;vertical-align:middle;height:45px",
                             verbatimTextOutput("diffsplice.method")
                         ),
                         div(style="display:inline-block;vertical-align:middle;height:40px",
                             actionButton('run.diffsplice','Run',icon("send outline icon"),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         ),
                         br(),
                         br(),
                         verbatimTextOutput("run.DAS.info")
                         )
                     )
            ),
            
            ##-----------DE DAS DTU transcripts------------
            fluidRow(
              column(width = 12,
                     box(title='Part of deltaPS',
                         width = 13,status = 'primary', solidHeader = T,
                         DT::dataTableOutput('deltaPS.table'),
                         actionButton(inputId = 'save.deltaPS.table',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                     )
              )
            )
            ),
    #========================>> Results summary- panel <<=======================
    tabItem('resultssummary',
            fluidRow(
              column(width = 3,
                     box(title='3D target tables',
                         width = 13,status = 'primary', solidHeader = T,
                         HTML('<h4><strong>Visulise statistics</strong></h4>'),
                         selectInput(inputId = 'stat.type',label = 'Targets',
                                     choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts'),
                                     selected = 'DE genes'
                         ),
                         uiOutput('show.contrast.groups'),
                         sliderInput(inputId = 'top.stat',label = 'Top ranked statistics in each contrast',
                                     min = 10,max = 100,value = 10,step = 10),
                         hr(),
                         actionButton(inputId = 'save.stat',label = 'Save all tables',icon = icon('download',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                         br(),
                         br(),
                         HTML('<p align="justify">Up to top 100 ranked targets can be visulise in this App. To visulise the statistics of
                              full gene/transcript sets, please click the "Save all tables" button. The full tables of statistics and target numbers in this page will be saved
                              as csv files in the "result" folder.</p>')
                         )
                         ),
              column(width = 9,
                     uiOutput('show.top.stat.table'),
                     box(title='DE DAS DTU target numbers',
                         width = 13,status = 'primary', solidHeader = T,
                         tableOutput('DDD.numbers')
                     ),
                     fluidRow(
                       column(width = 6,
                              box(title='DE vs DAS genes',
                                  width = 13,status = 'primary', solidHeader = T,
                                  tableOutput('DE.vs.DAS')
                              )
                       ),
                       column(width = 6,
                              box(title='DE vs DTU transcripts',
                                  width = 13,status = 'primary', solidHeader = T,
                                  tableOutput('DE.vs.DTU')
                              )
                       )
                     )
              )
                     ),
            fluidRow(
              column(width = 12,
                     tabBox(width = 13,
                            
                            # Title can include an icon
                            # title = tagList(shiny::icon("gear"), "Loaded tables"),
                            title = "Up-down regulation",side = 'right',
                            tabPanel("DE genes",
                                     plotOutput("up.down.bar.plot.DE.genes")
                            ),
                            tabPanel("DAS genes",
                                     plotOutput("up.down.bar.plot.DAS.genes")
                            ),
                            tabPanel("DE transcripts",
                                     plotOutput("up.down.bar.plot.DE.trans")
                            ),
                            tabPanel("DTU transcripts",
                                     plotOutput("up.down.bar.plot.DTU.trans")
                            )
                     ),
                     actionButton(inputId = 'save.up.down.bar.plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4;float:right; margin-top: 25px;"),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "updown.x.rotate",label = 'Rotate x-labels',value = '0',
                                      min = 0,max = 360,step = 15,width = "100%")
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "updown.x.hjust",label = 'Horizontal adjust x-labels',value = '0.5',
                                      min = 0,max = 1,step = 0.1,width = "100%")
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "updown.plot.res",label = 'PNG plot resolution',value = '150',
                                      min = 0,max = 600,step = 60,width = "100%")
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "updown.plot.height",label = 'Plot height (inch)',
                                      value = 5,width = "100%",step = 0.2)
                     ),
                     div(style="float:right;margin-right: 25px;",
                         numericInput(inputId = "updown.plot.width",label = 'Plot width (inch)',
                                      value = 8,width = "100%",step = 0.2)
                     ),
                     br(),
                     br(),
                     br(),
                     br()
              ),
              column(width = 12,
                     box(title='Comparison across individual contrast groups',
                         width = 13,status = 'primary', solidHeader = T,
                         uiOutput('across.contrast.euler'),
                         splitLayout(cellWidths = c("50%", "50%"),
                                     plotOutput("DE.genes.euler"),
                                     plotOutput("DAS.genes.euler")
                         ),
                         br(),
                         splitLayout(cellWidths = c("50%", "50%"),
                                     plotOutput("DE.trans.euler"),
                                     plotOutput("DTU.trans.euler")
                         ),
                         actionButton(inputId = 'save.across.contrast.euler.plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.contrast.euler.plot.res",label = 'PNG plot resolution',value = '150',
                                          min = 0,max = 600,step = 60,width = "100%")
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.contrast.euler.plot.height",label = 'Plot height (inch)',
                                          value = 4.7,width = "100%",step = 0.2)
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.contrast.euler.plot.width",label = 'Plot width (inch)',
                                          value = 5,width = "100%",step = 0.2)
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.contrast.euler.plot.scale",label = 'Plot scales (0-1)',
                                          value = 1,min = 0,max = 1,step = 0.1,width = "100%")
                         )
                     ),
                     box(title='Comparison across targets in individual conatrast groups',
                         width = 13,status = 'primary', solidHeader = T,
                         uiOutput('across.target.euler'),
                         splitLayout(cellWidths = c("50%", "50%"),
                                     plotOutput("DEvsDAS.euler"),
                                     plotOutput("DEvsDTU.euler")
                         ),
                         actionButton(inputId = 'save.across.target.euler.plot',label = 'Save',
                                      icon = icon('download',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.target.euler.plot.res",label = 'PNG plot resolution',value = '150',
                                          min = 0,max = 600,step = 60,width = "100%")
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.target.euler.plot.height",label = 'Plot height (inch)',
                                          value = 4.7,width = "100%",step = 0.2)
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.target.euler.plot.width",label = 'Plot width (inch)',
                                          value = 5,width = "100%",step = 0.2)
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "across.target.euler.plot.scale",label = 'Plot scales (0-1)',
                                          value = 1,min = 0,max = 1,step = 0.1,width = "100%")
                         )
                     ),
                     box(title='Union set of all contrast groups',
                         width = 13,status = 'primary', solidHeader = T,
                         plotOutput("genes.flow.chart"),
                         hr(),
                         plotOutput("trans.flow.chart"),
                         actionButton(inputId = 'save.union.flow.chart',label = 'Save',
                                      icon = icon('download',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "flow.plot.res",label = 'PNG plot resolution',value = '150',
                                          min = 0,max = 600,step = 60,width = "100%")
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "flow.plot.height",label = 'Plot height (inch)',
                                          value = 5,width = "100%",step = 0.2)
                         ),
                         div(style="float:right;margin-right: 25px;",
                             numericInput(inputId = "flow.plot.width",label = 'Plot width (inch)',
                                          value = 8.7,width = "100%",step = 0.2)
                         )
                         
                     )
              )
            )
    ),
    
    #========================>> Advanced plot <<=======================
    tabItem('advancedanalysis',
            ##----------Heatmap------------
            fluidRow(
              # column(width = 12,
              box(title='Heatmap',
                  width = 3,status = 'primary', solidHeader = T,
                  
                  HTML('<h4><strong>Choose targets</strong></h4>'),
                  radioButtons(inputId = 'heatmap.select.or.upload',label = '',
                               choices = c('Select targets','Upload target list'),inline = T),
                  selectInput(inputId = 'heatmap.target.type',label = 'Select targets',
                              choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts')
                  ),
                  hr(),
                  fileInput("heatmap.target.list.input", "OR Choose target list csv file",
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv")
                  ),
                  radioButtons(inputId = 'heatmap.targetlist.type',label = 'Target list type',inline = T,
                               choices = c('Gene list','Transcript list')
                  ),
                  hr(),
                  HTML('<h4><strong>Clustering</strong></h4>'),
                  selectInput(inputId = 'dist.method',label = 'Distance measurement',
                              choices = c("euclidean","maximum","manhattan","canberra","binary","minkowski")),
                  selectInput(inputId = 'cluster.method',label = 'Clustering method',
                              choices = c("ward.D","ward.D2","average","complete","single","mcquitty","median","centroid"),
                              selected = 'complete'),
                  numericInput(inputId = 'cluster.number',label = 'Cluster number',value = 10,min = 1,max = 50),
                  HTML('See the R functions for hierarchical clustering: <a href="https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust" target="_blank">hclust</a> and
                       <a href="https://www.rdocumentation.org/packages/proxy/versions/0.4-22/topics/dist" target="_blank">dist</a>'),
                  br(),
                  br(),
                  actionButton(inputId = 'plot.heatmap',label = 'Plot heatmap',
                               icon = icon('pencil',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                  ),
              box(title=NULL,
                  width = 9,status = 'primary', solidHeader = T,
                  plotOutput('heatmap.plot.panel',height = 730),
                  actionButton(inputId = 'save.target.in.cluster',label = 'Save target list in clusters',
                               icon = icon('download',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px; margin-left: 25px"),
                  actionButton(inputId = 'save.heatmap',label = 'Save heatmap',
                               icon = icon('download',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                  div(style="float:right;margin-right: 25px;",
                      numericInput(inputId = "heatmap.plot.res",label = 'PNG plot resolution',value = '150',
                                   min = 0,max = 600,step = 60,width = "100%")
                  ),
                  div(style="float:right;margin-right: 25px;",
                      numericInput(inputId = "heatmap.plot.height",label = 'Plot height (inch)',
                                   value = 8,width = "100%",step = 0.2)
                  ),
                  div(style="float:right;margin-right: 25px;",
                      numericInput(inputId = "heatmap.plot.width",label = 'Plot width (inch)',
                                   value = 4.5,width = "100%",step = 0.2)
                  )
              )
            ),
            fluidRow(
              ##---------------multiple plots----------------
              # column(width = 12,
              box(title='Profiles',
                  width = 12,status = 'primary', solidHeader = T,
                  column(width = 6,
                         wellPanel(
                           HTML('<h4><strong>Visualise single gene plot</strong></h4>'),
                           textInput(inputId = 'gene.id',label = 'Input a gene',value = ''),
                           HTML('E.g. AT1G01060 in Arabidopsis (gene name must match to provided transcript-gene mapping).'),
                           br(),
                           radioButtons(inputId = 'profile.data.type',label = 'Select data type',
                                        choices = c('TPM','Read counts'),inline = T),
                           radioButtons(inputId = 'profile.filter.lowexpressed',label = 'Filter low expressed transcripts?',
                                        choices = c('Yes','No'),inline = T),
                           actionButton(inputId = 'make.profile.plot',label = 'Plot',
                                        icon = icon('pencil',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float"),
                           br()
                         )
                  ),
                  column(width = 6,
                         wellPanel(
                           HTML('<h4><strong>Make multiple gene plots</strong></h4>'),
                           fileInput("multiple.gene.input", "Choose gene list csv file",
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                           ),
                           radioButtons(inputId = 'multiple.plot.type',label = 'Select plot type',
                                        choices = c('Abundance','PS','Both'),inline = T),
                           radioButtons(inputId = 'multi.plot.format',label = 'Select format',
                                        choices = c('png','pdf','both'),inline = T),
                           actionButton(inputId = 'make.multiple.plot',label = 'Run',
                                        icon = icon('send outline icon',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float"),
                           br()
                         )
                  )
              ),
              # column(width = 12,
              box(title=NULL,
                  width = 12,status = 'primary', solidHeader = T,
                  plotlyOutput('profile.plot.panel'),
                  br(),
                  hr(),
                  br(),
                  plotlyOutput('ps.plot.panel'),
                  actionButton(inputId = 'save.ps.plot',label = 'Save PS plot',
                               icon = icon('download',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px; margin-left: 25px"),
                  actionButton(inputId = 'save.abundance.plot',label = 'Save abundance plot',
                               icon = icon('download',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
                  div(style="float:right;margin-right: 25px;",
                      numericInput(inputId = "ps.plot.res",label = 'PNG plot resolution',value = '150',
                                   min = 0,max = 600,step = 60,width = "100%")
                  ),
                  div(style="float:right;margin-right: 25px;",
                      numericInput(inputId = "ps.plot.height",label = 'Plot height (inch)',
                                   value = 4.5,width = "100%",step = 0.2)
                  ),
                  div(style="float:right;margin-right: 25px;",
                      numericInput(inputId = "ps.plot.width",label = 'Plot width (inch)',
                                   value = 7,width = "100%",step = 0.2)
                  )
              )
            ),
            fluidRow(
              ##---------------GO plots----------------
              column(width = 3,
                     box(title = 'GO annotation',
                         width = 13,status = 'primary', solidHeader = T,
                         HTML('<h4><strong>Generate target list</strong></h4>'),
                         actionButton(inputId = 'generate.target.list',label = 'Generate',
                                      icon = icon('send outline icon',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                         br(),
                         br(),
                         HTML('Generate lists of DE genes, DAS genes, DE transcripts and DTU transcripts, which are the union
                              set across all contrast groups.'),
                         hr(),
                         fileInput("go.annotation.input", "Choose GO annotation csv file",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                         ),
                         uiOutput('show.go.table.column'),
                         radioButtons(inputId = 'go.plot.label',label = 'Select gene list type',
                                      choices  = c('DE genes','DAS genes'),inline = T),
                         numericInput(inputId = 'go.text.cut',label = 'GO term text length',value = 50),
                         actionButton(inputId = 'make.go.plot',label = 'Plot',
                                      icon = icon('pencil',lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                         
                         )
              ),
              column(width = 9,
                     box(title = NULL,
                         width = 13,status = 'primary', solidHeader = T,
                         h4('Loaded GO annotation table'),
                         DT::dataTableOutput('go.table')
                     )
              )
    ),
    fluidRow(
      box(title = NULL,
          width = 12,status = 'primary', solidHeader = T,
          h4('GO annotation plot'),
          uiOutput('show.go.plot.panel'),
          actionButton(inputId = 'save.go.plot',label = 'Save',
                       icon = icon('download',lib = 'font-awesome'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4;float: right; margin-top: 25px"),
          div(style="float:right;margin-right: 25px;",
              numericInput(inputId = "go.plot.res",label = 'PNG plot resolution',value = '150',
                           min = 0,max = 600,step = 60,width = "100%")
          ),
          div(style="float:right;margin-right: 25px;",
              numericInput(inputId = "go.plot.height",label = 'Plot height (inch)',
                           value = 4.5,width = "100%",step = 0.2)
          ),
          div(style="float:right;margin-right: 25px;",
              numericInput(inputId = "go.plot.width",label = 'Plot width (inch)',
                           value = 7.8,width = "100%",step = 0.2)
          )
      )
    )
            ),
    
    #=======================>> Generate report panel <<============================
    tabItem("report",
            tags$head(tags$style("#TxtOut {white-space: normal;}")),
            box(title = 'Parameter summary',
                width = 12,status = 'primary', solidHeader = T,
                HTML('<h5><strong>Table: Paramter used for the analysis. </strong> If the parameters are incorrect, please
                     double click the the cell in the table to correct it.</h5>'),
                DT::DTOutput('para.summary'),
                br(),
                br(),
                actionButton(inputId = 'save.para.summary',label = 'Save',
                             icon = icon('download',lib = 'font-awesome'),
                             style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                ),
            actionButton(inputId = 'generate.report',label = 'Generate report',
                         icon = icon('send outline icon',lib = 'font-awesome'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                         font-size:200%; padding: 30px 60px; margin-left:40%")
            )
    )
            )

ui <- dashboardPage(
  mainheader,
  mainsidebar,
  mainbody
)



Logged = F
my_username <- "User"
my_password <- "3DRNAseq"

server <- function(input, output, session) {
  # set data size
  
  #Log in ##############################################################
  # ============================= password ============================ #
  values <- reactiveValues(authenticated = FALSE)
  
  # Return the UI for a modal dialog with data selection input. If 'failed'
  # is TRUE, then display a message that the previous value was invalid.
  dataModal <- function(failed = FALSE) {
    modalDialog(
      textInput("username", "Username:",value = 'User'),
      passwordInput("password", "Password:",value = '3DRNAseq'),
      footer = tagList(
        # modalButton("Cancel"),
        div(style="display:inline-block;vertical-align:middle;height:40px;float: left",
            verbatimTextOutput('log.info')
        ),
        div(style="display:inline-block;vertical-align:middle;height:40px",
            actionButton("ok", "OK",
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4")
        )
      ),size = 's'
    )
  }
  
  # Show modal when button is clicked.
  # This `observe` is suspended only whith right user credential
  
  obs1 <- observe({
    showModal(dataModal())
  })
  
  # When OK button is pressed, attempt to authenticate. If successful,
  # remove the modal.
  
  obs2 <- observe({
    req(input$ok)
    isolate({
      Username <- input$username
      Password <- input$password
    })
    Id.username <- which(my_username == Username)
    Id.password <- which(my_password == Password)
    if (length(Id.username) > 0 & length(Id.password) > 0) {
      if (Id.username == Id.password) {
        Logged <<- TRUE
        values$authenticated <- TRUE
        obs1$suspend()
        removeModal()
      } else {
        values$authenticated <- FALSE
      }
    }
  })
  
  observeEvent(input$ok,{
    output$log.info <- renderText({
      if (values$authenticated) "OK!!!!!"
      else "Wrong password!!!"
    })
  })
  # ============================= password ============================ #
  
  DDD.data <- reactiveValues(
    data.folder=NULL,
    figure.folder=NULL,
    result.folder=NULL,
    report.folder=NULL,
    path=getwd(),
    mapping=NULL,
    samples=NULL,
    samples_new=NULL,
    trans_counts=NULL,
    genes_counts=NULL,
    trans_TPM=NULL,
    genes_TPM=NULL,
    target_high=NULL,
    trans_batch=NULL,
    genes_batch=NULL,
    trans_dge=NULL,
    genes_dge=NULL,
    contrast=NULL,
    deltaPS=NULL,
    genes_log2FC=NULL,
    trans_log2FC=NULL,
    DE_genes=NULL,
    DE_genes_stat=NULL,
    DAS_genes=NULL,
    DE_trans=NULL,
    DTU_trans=NULL,
    loaded.data=NULL,
    genes.ann=NULL,
    trans.ann=NULL
    # info=list()
  )
  
  output$messageMenu <- renderMenu({
    dropdownMenu(type = "tasks", badgeStatus = "success",
                 taskItem(value = 90, color = "green",
                          "Documentation"
                 ),
                 taskItem(value = 17, color = "aqua",
                          "Project X"
                 ),
                 taskItem(value = 75, color = "yellow",
                          "Server deployment"
                 ),
                 taskItem(value = 80, color = "red",
                          "Overall project"
                 )
    )
  })
  
  
  ##------------->>   Data generation      <<--------------
  ##-------------  Loaded data in workspace   --------------
  
  observe({
    Data <- c(
      'mapping',
      'samples',
      'samples_new',
      'trans_counts',
      'genes_counts',
      'trans_TPM',
      'genes_TPM',
      'target_high',
      'trans_batch',
      'genes_batch',
      'trans_dge',
      'genes_dge',
      'contrast',
      'deltaPS',
      'DE_genes',
      'DAS_genes',
      'DE_trans',
      'DTU_trans'
    )
    if(is.null(DDD.data$data.folder)){
      Existed1 <- F
    } else {
      Existed1 <- sapply(Data,function(i) file.exists(paste0(DDD.data$data.folder,'/',i,'.RData')))
    }
    DDD.data.list <- reactiveValuesToList(DDD.data)
    Existed2 <- sapply(1:length(Data),function(i) !is.null(DDD.data.list[[Data[i]]]))
    
    Description=c('Gene-transcript mapping',
                  'Sample information',
                  'Sample information (Sequencing replicates are merged if exist)',
                  'Transcript read counts (sequencing replicates are merged if exist)',
                  'Gene read counts (sequencing replicates are merged if exist)',
                  'Transcript level TPM',
                  'Gene level TPM',
                  'Remaining transcripts and genes after filtering the low expressed',
                  'Estimated transcript batch effect only if exist',
                  'Estimated gene batch effect only if exist',
                  'Transcript level normalization factor',
                  'Gene level normalization factor',
                  'Contrast groups of conditions for comparisons',
                  'delta Percent spliced value based on contrast groups',
                  'DE genes in different contrast groups',
                  'DAS genes in different contrast groups',
                  'DE transcripts in different contrast groups',
                  'DTU transcripts in different contrast groups'
    )
    
    Generate=c('User provided',
               'User provided',
               'Data pre-processing',
               'Data pre-processing',
               'Data pre-processing',
               'Data generation',
               'Data generation',
               'Data pre-processing',
               'Data pre-processing',
               'Data pre-processing',
               'Data pre-processing',
               'Data pre-processing',
               'User provided',
               'DE DAS and DTU',
               'DE DAS and DTU',
               'DE DAS and DTU',
               'DE DAS and DTU',
               'DE DAS and DTU'
    )
    
    Use=c('All steps',
          'All steps',
          'All steps',
          'Data pre-processing and DE DAS and DTU',
          'Data pre-processing and DE DAS and DTU',
          'DE DAS and DTU and Advanced plot',
          'Advanced plot',
          'All steps',
          'Data pre-processing and DE DAS and DTU',
          'Data pre-processing and DE DAS and DTU',
          'Data pre-processing and DE DAS and DTU',
          'Data pre-processing and DE DAS and DTU',
          'DE DAS and DTU and Result summary',
          'DE DAS and DTU',
          'Result summary and Advanced plot',
          'Result summary and Advanced plot',
          'Result summary and Advanced plot',
          'Result summary and Advanced plot'
    )
    
    x <- data.frame(Data=Data,`In data folder`=Existed1,
                    `In workspace`=Existed2,
                    Description=Description,
                    `Generate from`=Generate,
                    `Use in`=Use)
    colnames(x) <- c('Data','In data folder','In workspace','Description','Generate from','Use in')
    DDD.data$loaded.data <- x
  })
  
  DDD.data.size <- reactive({
    x <- sum(sapply(reactiveValuesToList(DDD.data),object.size))
    utils:::format.object_size(x, "auto")
  })
  output$DDD.data.size <- renderText({
    paste0('Data in use: ',DDD.data.size())
  })
  
  output$loaded.data.info <- renderTable(
    DDD.data$loaded.data
  )
  
  observeEvent(input$loaded.data.all,{
    idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F & DDD.data$loaded.data$`In data folder`==T]
    withProgress(message = 'Loading exisiting data',
                 detail = 'This may take a while...', value = 0, {
                   for(i in idx){
                     incProgress(1/length(idx))
                     cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                     load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                     if(!exists(i))
                       next
                     DDD.data[[i]] <- get(i)
                   }
                 })
    message('Done!!!')
    showNotification("Done!!!")
    
    if(length(idx==0))
      return(NULL)
    ###disable the buttons
  })
  
  observeEvent(input$loaded.data.preprocessing,{
    idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                       DDD.data$loaded.data$`In data folder`==T &
                                       grepl(pattern = 'All steps|Data pre-processing',x = DDD.data$loaded.data$`Use in`)]
    withProgress(message = 'Loading exisiting data',
                 detail = 'This may take a while...', value = 0, {
                   for(i in idx){
                     incProgress(1/length(idx))
                     cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                     load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                     if(!exists(i))
                       next
                     DDD.data[[i]] <- get(i)
                   }
                 })
    message('Done!!!')
    showNotification("Done!!!")
    
    if(length(idx==0))
      return(NULL)
    ###disable the buttons
  })
  
  observeEvent(input$loaded.data.ddd,{
    idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                       DDD.data$loaded.data$`In data folder`==T &
                                       grepl(pattern = 'All steps|DE DAS and DTU',x = DDD.data$loaded.data$`Use in`)]
    withProgress(message = 'Loading exisiting data',
                 detail = 'This may take a while...', value = 0, {
                   for(i in idx){
                     incProgress(1/length(idx))
                     cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                     load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                     if(!exists(i))
                       next
                     DDD.data[[i]] <- get(i)
                   }
                 })
    message('Done!!!')
    showNotification("Done!!!")
    
    if(length(idx==0))
      return(NULL)
    ###disable the buttons
  })
  
  
  observeEvent(input$loaded.data.summary,{
    idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                       DDD.data$loaded.data$`In data folder`==T &
                                       grepl(pattern = 'All steps|Result summary',x = DDD.data$loaded.data$`Use in`)]
    withProgress(message = 'Loading exisiting data',
                 detail = 'This may take a while...', value = 0, {
                   for(i in idx){
                     incProgress(1/length(idx))
                     cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                     load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                     if(!exists(i))
                       next
                     DDD.data[[i]] <- get(i)
                   }
                 })
    message('Done!!!')
    showNotification("Done!!!")
    
    if(length(idx==0))
      return(NULL)
    ###disable the buttons
  })
  
  observeEvent(input$loaded.data.advanced,{
    idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                       DDD.data$loaded.data$`In data folder`==T &
                                       grepl(pattern = 'All steps|Advanced plot',x = DDD.data$loaded.data$`Use in`)]
    withProgress(message = 'Loading exisiting data',
                 detail = 'This may take a while...', value = 0, {
                   for(i in idx){
                     incProgress(1/length(idx))
                     cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                     load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                     if(!exists(i))
                       next
                     DDD.data[[i]] <- get(i)
                   }
                 })
    message('Done!!!')
    showNotification("Done!!!")
    
    if(length(idx==0))
      return(NULL)
    ###disable the buttons
  })
  
  
  ##---------------working directory-----------------
  shinyDirChoose(input, 'wdir', roots = c(getVolumes()()))
  wdir <- reactive({
    input$wdir
  })
  
  observeEvent(ignoreNULL = T,
               eventExpr = {
                 input$wdir
               },
               handlerExpr = {
                 if(is.integer(wdir())){
                   path.folder <- NULL
                 } else {
                   root <- getVolumes()()[wdir()$root]
                   root <- gsub('[\\]','',normalizePath(root))
                   current_wk <- unlist(wdir()$path[-1])
                   path.folder <- file.path(root, paste(current_wk, collapse = .Platform$file.sep))
                   setwd(path.folder)
                 }
                 DDD.data$path <- path.folder
               })
  
  output$wdir.info <- renderText({
    DDD.data$path
  })
  
  observe({
    DDD.data$data.folder <- paste0(DDD.data$path,'/data')
    DDD.data$figure.folder <- paste0(DDD.data$path,'/figure')
    DDD.data$result.folder <- paste0(DDD.data$path,'/result')
    DDD.data$report.folder <- paste0(DDD.data$path,'/report')
  })
  
  #####data save folder
  data.folder <- observeEvent(input$create.folders,{
    data.folder <- paste0(DDD.data$path,'/data')
    if(!file.exists(data.folder))
      dir.create(data.folder,recursive = T)
    DDD.data$data.folder <- data.folder
    showNotification('data folder is created.')
    message('data folder is created.')
  })
  
  #####figure save folder
  figure.folder <- observeEvent(input$create.folders,{
    figure.folder <- paste0(DDD.data$path,'/figure')
    if(!file.exists(figure.folder))
      dir.create(figure.folder,recursive = T)
    DDD.data$figure.folder <- figure.folder
    showNotification('figure folder is created.')
    message('figure folder is created.')
  })
  
  #####results save folder
  result.folder <- observeEvent(input$create.folders,{
    result.folder <- paste0(DDD.data$path,'/result')
    if(!file.exists(result.folder))
      dir.create(result.folder,recursive = T)
    DDD.data$result.folder <- result.folder
    showNotification('result folder is created.')
    message('result folder is created.')
  })
  
  #####reports save folder
  report.folder <- observeEvent(input$create.folders,{
    report.folder <- paste0(DDD.data$path,'/report')
    if(!file.exists(report.folder))
      dir.create(report.folder,recursive = T)
    DDD.data$report.folder <- report.folder
    showNotification('report folder is created.')
    message('report folder is created.')
  })
  
  
  
  ##mapping table####
  observe({
    inFile <- input$mapping.input
    if (is.null(inFile))
      return(NULL)
    # withProgress(message = 'Loading mapping information', value = 0, {
    mapping <- read.csv(inFile$datapath, header = T)
    colnames(mapping) <- c('TXNAME','GENEID')
    rownames(mapping) <- mapping$TXNAME
    save(mapping,file=paste0(DDD.data$data.folder,'/mapping.RData'))
    DDD.data$mapping <- mapping
    # })
  })
  
  output$mapping.output <- DT::renderDataTable({
    DDD.data$mapping[1:25,]
  },options = list(scrollX = TRUE,
                   columnDefs = list(list(className = 'dt-center', targets="_all"))
  ),rownames= FALSE)
  
  ##sample table####
  samples <- observe({
    inFile <- input$sample.input
    if (is.null(inFile))
      return(NULL)
    # withProgress(message = 'Loading sample information', value = 0, {
    samples <- read.csv(inFile$datapath, header = T)
    # })
    if('srep' %in% colnames(samples)){
      if(length(unique(samples$srep))==1){
        samples$sample.name <- paste0(samples$condition,'_',samples$brep)
      } else {
        samples$sample.name <- paste0(samples$condition,'_',samples$brep,'_',samples$srep)
      }
    } else {
      
    }
    save(samples,file=paste0(DDD.data$data.folder,'/samples.RData'))
    
    DDD.data$samples <- samples
    
  })
  
  output$sample.output <- DT::renderDataTable({
    if(!is.null(DDD.data$samples_new))
      DDD.data$samples_new else DDD.data$samples
  },options = list(scrollX = TRUE,
                   columnDefs = list(list(className = 'dt-left', targets="_all"))))
  
  ##------------------->> generate  expression  <<--------------------
  ##--------------------generate genes expression---------------------
  observe({
    if('abundance.h5' %in% DDD.data$samples$path)
      updateSelectInput(session,inputId = "tximport.quant.method",selected = "kallisto")
  })
  
  observe({
    if(input$generate.new.genes.data=='Load txi_genes.RData'){
      updateActionButton(session,inputId = 'run.txi.genes',label = 'Load',icon = icon('upload'))
    } else {
      updateActionButton(session,inputId = 'run.txi.genes',label = 'Run')
    }
  })
  
  observeEvent(input$run.txi.genes,{
    if (is.null(DDD.data$samples))
      return(NULL)
    if(input$generate.new.genes.data=='Load txi_genes.RData'){
      cat('\nLoading gene expression data: txi_genes.RData')
      load(paste0(DDD.data$data.folder,'/txi_genes.RData'))
      txi_genes <- txi_genes
      message(paste0('txi_genes.RData is loaded from folder: ',DDD.data$data.folder))
    } else {
      cat('\nRunning tximport to generate gene expression using method: ',
          input$tximport.method)
      
      data.files <-  as.character(DDD.data$samples$path)
      if(!all(file.exists(data.files))){
        warning('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
        showNotification('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
        return(NULL)
      }
      
      withProgress(message = 'Generate gene expression...', value = 0, {
        incProgress(0.1)
        txi_genes <- tximport(data.files,dropInfReps = T,
                              type = input$tximport.quant.method, tx2gene = DDD.data$mapping,
                              countsFromAbundance = input$tximport.method)
        incProgress(0.7)
        colnames(txi_genes$counts)<-DDD.data$samples$sample.name
        write.csv(txi_genes$counts,file=paste0(DDD.data$result.folder,'/counts_genes.csv'))
        incProgress(0.8)
        colnames(txi_genes$abundance)<-DDD.data$samples$sample.name
        write.csv(txi_genes$abundance,file=paste0(DDD.data$result.folder,'/TPM_genes.csv'))
        incProgress(0.9)
        colnames(txi_genes$length)<-DDD.data$samples$sample.name
        save(txi_genes,file=paste0(DDD.data$data.folder,'/txi_genes.RData'))
        genes_TPM <- txi_genes$abundance
        save(genes_TPM,file=paste0(DDD.data$data.folder,'/genes_TPM.RData'))
        message(paste0('txi_genes.RData and genes_TPM.RData are saved in folder: ',DDD.data$data.folder))
        incProgress(1)
      })
    }
    DDD.data$genes_counts=txi_genes$counts
    DDD.data$genes_TPM <- txi_genes$abundance
    message('Done!!!')
    showNotification("Done!!!")
  })
  
  rungeneinfo <- reactive({
    if(input$run.txi.genes>0){
      if(input$generate.new.genes.data=='Load txi_genes.RData'){
        paste0('txi_genes.RData is loaded from the data folder.')
      } else {
        paste0('txi_genes.RData is saved in the data folder.')
      }
    } else {
      if(file.exists(paste0(DDD.data$data.folder,'/txi_genes.RData'))){
        paste0('txi_genes.RData is already in the data folder.\nGenerate new?')
      } else {
        paste0('txi_genes.RData is not in the data folder.\nPlease generate new.')
      }
    }
  })
  
  output$rungeneinfo <- renderText({
    rungeneinfo()
  })
  
  ##-------------------- generate trans expression --------------------
  ###---update actionbutton if load data
  observe({
    if(input$generate.new.trans.data=='Load txi_trans.RData'){
      updateActionButton(session,inputId = 'run.txi.trans',label = 'Load',icon = icon('upload'))
    } else {
      updateActionButton(session,inputId = 'run.txi.trans',label = 'Run')
    }
  })
  
  observeEvent(input$run.txi.trans,ignoreNULL = T,{
    if (is.null(DDD.data$samples))
      return(NULL)
    
    if(input$generate.new.trans.data=='Load txi_trans.RData'){
      cat('\nLoading transcript expression data: txi_trans.RData')
      load(paste0(DDD.data$data.folder,'/txi_trans.RData'))
      txi_trans <- txi_trans
      message(paste0('txi_trans.RData is loaded from folder: ',DDD.data$data.folder))
    } else {
      cat(paste0('\nRunning tximport to generate transcript expression using method: ',
                 input$tximport.method))
      
      data.files <-  as.character(DDD.data$samples$path)
      if(!all(file.exists(data.files))){
        warning('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
        showNotification('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
        return(NULL)
      }
      
      # cat(paste0(data.files[1]))
      withProgress(message = 'Generate transcript expression...', value = 0, {
        incProgress(0.1)
        txi_trans<- tximport(data.files, type = input$tximport.quant.method, tx2gene = NULL,
                             countsFromAbundance = input$tximport.method,
                             txOut = T,dropInfReps = T)
        incProgress(0.7)
        colnames(txi_trans$counts)<-DDD.data$samples$sample.name
        write.csv(txi_trans$counts,file=paste0(DDD.data$result.folder,'/counts_trans.csv'))
        incProgress(0.8)
        colnames(txi_trans$abundance)<-DDD.data$samples$sample.name
        write.csv(txi_trans$abundance,file=paste0(DDD.data$result.folder,'/TPM_trans.csv'))
        incProgress(0.9)
        colnames(txi_trans$length)<-DDD.data$samples$sample.name
        save(txi_trans,file=paste0(DDD.data$data.folder,'/txi_trans.RData'))
        trans_TPM <- txi_trans$abundance
        save(trans_TPM,file=paste0(DDD.data$data.folder,'/trans_TPM.RData'))
        message(paste0('txi_trans.RData and trans_TPM.RData are saved in folder: ',DDD.data$data.folder))
        incProgress(1)
      })
    }
    DDD.data$trans_counts <- txi_trans$counts
    DDD.data$trans_TPM <- txi_trans$abundance
    message('Done!!!')
    showNotification("Done!!!")
  })
  
  
  runtransinfo <- reactive({
    if(input$run.txi.trans>0){
      if(input$generate.new.trans.data=='Load txi_trans.RData'){
        paste0('txi_trans.RData is loaded from the data folder.')
      } else {
        paste0('txi_trans.RData is saved in the data folder.')
      }
    } else {
      if(file.exists(paste0(DDD.data$data.folder,'/txi_trans.RData'))){
        paste0('txi_trans.RData is already in the data folder.\nGenerate new?')
      } else {
        paste0('txi_trans.RData is not in the data folder.\nPlease generate new.')
      }
    }
  })
  
  output$runtransinfo <- renderText({
    runtransinfo()
  })
  
  
  
  ##-------------->>    data pre-processing   <<--------------
  ##-------------  Step 1: merge sequencing reps -------------
  
  
  seqrepinfo <- reactive({
    if(is.null(DDD.data$samples)){
      info <- HTML('Load samples.csv to check seq-reps.')
    } else {
      if(all(DDD.data$samples$srep=='srep1') | !('srep' %in% colnames(DDD.data$samples))){
        info <- HTML('The data has no sequencing replicates. Please skip.')
      } else if(!is.null(DDD.data$samples_new)) {
        info <- HTML('The sequencing replicates have been merged.')
      } else {
        info <- HTML(paste0('The data has ', length(unique(DDD.data$samples$srep)),' sequencing replicates. Please merge.'))
      }
    }
    return(info)
  })
  output$seqrep.info <- renderText(seqrepinfo())
  
  observe({
    if(input$merge.or.skip=='Merge'){
      updateActionButton(session,inputId = 'run.seq.merge',label = 'Merge',icon = icon("plus-square"))
    } else {
      updateActionButton(session,inputId = 'run.seq.merge',label = 'Skip',icon = icon("angle-right"))
    }
  })
  
  
  
  ###---update the sample information after merge---
  observeEvent(input$run.seq.merge,{
    if(is.null(DDD.data$samples))
      return(NULL)
    if(input$merge.or.skip=='Merge'){
      samples_new <- DDD.data$samples[DDD.data$samples$srep==DDD.data$samples$srep[1],]
      samples_new$srep <- 'merged'
    } else {
      samples_new <- DDD.data$samples
    }
    save(samples_new,file=paste0(DDD.data$data.folder,'/samples_new.RData'))
    message(paste0('samples_new.RData is saved in folder: ',DDD.data$data.folder))
    DDD.data$samples_new <- samples_new
    
  })
  
  ###---updateRadioButtons
  observe({
    if(is.null(DDD.data$samples))
      return(NULL)
    if('srep' %in% names(DDD.data$samples)){
      x <- 'Merge'
    } else {
      x <- 'Skip'
    }
    updateRadioButtons(session,label = 'Merge options', "merge.or.skip",selected = x)
  })
  #
  observeEvent(input$run.seq.merge,{
    if(is.null(DDD.data$samples) | is.null(DDD.data$trans_counts))
      return(NULL)
    
    if(input$merge.or.skip=='Merge'){
      cat('\nMerge sequencing replicates')
      idx <- paste0(DDD.data$samples$condition,'_',DDD.data$samples$brep)
      if(length(idx)==ncol(DDD.data$trans_counts)){
        trans_counts <- sumarrays(DDD.data$trans_counts,group = idx)
      } else {
        trans_counts <- DDD.data$trans_counts
      }
    } else {
      trans_counts <- DDD.data$trans_counts
    }
    save(trans_counts,file=paste0(DDD.data$data.folder,'/trans_counts.RData'))
    message(paste0('trans_counts.RData is saved in folder: ',DDD.data$data.folder))
    
    DDD.data$trans_counts <- trans_counts
  })
  
  observeEvent(input$run.seq.merge,{
    if(is.null(DDD.data$samples) | is.null(DDD.data$genes_counts))
      return(NULL)
    
    if(input$merge.or.skip=='Merge'){
      idx <- paste0(DDD.data$samples$condition,'_',DDD.data$samples$brep)
      if(length(idx)==ncol(DDD.data$genes_counts)){
        genes_counts <- sumarrays(DDD.data$genes_counts,group = idx)
      } else {
        genes_counts <- DDD.data$genes_counts
      }
    } else {
      genes_counts <- DDD.data$genes_counts
    }
    save(genes_counts,file=paste0(DDD.data$data.folder,'/genes_counts.RData'))
    message(paste0('genes_counts.RData is saved in folder: ',DDD.data$data.folder))
    
    showNotification("Done!!!")
    DDD.data$genes_counts <- genes_counts
  })
  
  
  
  observeEvent(input$run.seq.merge,{
    if(is.null(DDD.data$samples_new)){
      max.n <- 10
    } else {
      max.n <- nrow(DDD.data$samples_new)
    }
    updateSliderInput(session, "sample.n.cut", max=max.n)
  })
  #
  ##----------Step 2: filter low expressed genes------------
  
  observeEvent(input$run.filter,{
    if(is.null(DDD.data$trans_counts) | is.null(DDD.data$mapping))
      return(NULL)
    cat('\nFilter low expression')
    target_high <- low.expression.filter(abundance = DDD.data$trans_counts,
                                         mapping = DDD.data$mapping,
                                         abundance.cut = input$cpm.cut,
                                         sample.n = input$sample.n.cut,
                                         unit = 'counts',
                                         Log=F)
    save(target_high,file=paste0(DDD.data$data.folder,'/target_high.RData'))
    message(paste0('target_high.RData is saved in folder: ',DDD.data$data.folder))
    
    
    message('Done!!!')
    showNotification("Done!!!")
    DDD.data$target_high <- target_high
  })
  
  low.filter.info <- reactive({
    if(is.null(DDD.data$data.folder))
      return(NULL)
    if(input$run.filter>0){
      paste0('target_high.RData is saved in the data folder.')
    } else {
      if(!is.null(DDD.data$target_high)){
        paste0('target_high.RData is already uploaded in the workspace.')
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/target_high.RData'))){
          paste0('target_high.RData is already in "data" folder.\nCan upload directly in "Data generation".')
        } else {
          paste0('target_high.RData is not in "data" folder.\nPlease generate new.')
        }
      }
    }
  })
  
  
  output$low.filter.info <- renderText({
    if(is.null(low.filter.info()))
      return(NULL)
    low.filter.info()
  })
  
  
  
  ##--------------mean-variance plots---------------
  ###--mean variance trans trend plot
  mv.trans <- eventReactive(input$mv.plot.button,{
    condition <- paste0(DDD.data$samples_new$condition)
    ###---transcript levels
    counts.raw = DDD.data$trans_counts[rowSums(DDD.data$trans_counts>0)>0,]
    counts.filtered = DDD.data$trans_counts[DDD.data$target_high$trans_high,]
    # graphics.off()
    check.mean.variance(counts.raw = counts.raw,
                        counts.filtered = counts.filtered,
                        condition = condition)
  })
  
  mv.trans.plot <- eventReactive(input$mv.plot.button,{
    mv.trans.plot <- function(){
      fit.raw <- mv.trans()$fit.raw
      fit.filtered <- mv.trans()$fit.filtered
      par(mfrow=c(1,2))
      plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                       l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
      title('\n\nRaw counts')
      plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                       l = fit.filtered$l,lwd=2,col='black')
      title(paste0('\n\nFiltered counts (cutoff: cpm=',input$cpm.cut,'; sample=',input$sample.n.cut,')'))
      lines(fit.raw$l, col = "gold",lty=4,lwd=2)
      legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
             legend = c('low-exp removed','low-exp kept'))
    }
    return(mv.trans.plot)
  })
  #
  output$mv.trans.plot <- renderPlot({
    mv.trans.plot()()
  })
  
  ###--mean variance genes trend plot
  mv.genes <- eventReactive(input$mv.plot.button,{
    condition <- paste0(DDD.data$samples_new$condition)
    ###---genes levels
    counts.raw = DDD.data$genes_counts[rowSums(DDD.data$genes_counts>0)>0,]
    counts.filtered = DDD.data$genes_counts[DDD.data$target_high$genes_high,]
    # graphics.off()
    check.mean.variance(counts.raw = counts.raw,
                        counts.filtered = counts.filtered,
                        condition = condition)
  })
  
  mv.genes.plot <- eventReactive(input$mv.plot.button,{
    mv.genes.plot <- function(){
      fit.raw <- mv.genes()$fit.raw
      fit.filtered <- mv.genes()$fit.filtered
      par(mfrow=c(1,2))
      plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                       l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
      title('\n\nRaw counts')
      plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                       l = fit.filtered$l,lwd=2,col='black')
      title(paste0('\n\nFiltered counts (cutoff: cpm=',input$cpm.cut,'; sample=',input$sample.n.cut,')'))
      lines(fit.raw$l, col = "gold",lty=4,lwd=2)
      legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
             legend = c('low-exp removed','low-exp kept'))
    }
    return(mv.genes.plot)
  })
  #
  output$mv.genes.plot <- renderPlot({
    mv.genes.plot()()
  })
  
  observeEvent(input$save.mv.plot,{
    ##transcript level
    png(filename = paste0(DDD.data$figure.folder,'/Transcript mean-variance trend.png'),
        width = input$mv.plot.width,height = input$mv.plot.height,res=input$mv.plot.res, 
        units = 'in')
    mv.trans.plot()()
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Transcript mean-variance trend.pdf'),
        width = input$mv.plot.width,height = input$mv.plot.height)
    mv.trans.plot()()
    dev.off()
    
    ###gene level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Gene mean-variance trend.png'),
        width = input$mv.plot.width,height = input$mv.plot.height,res=input$mv.plot.res, 
        units = 'in')
    mv.genes.plot()()
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Gene mean-variance trend.pdf'),
        width = input$mv.plot.width,height = input$mv.plot.height)
    mv.genes.plot()()
    dev.off()
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  })
  
  output$save.mv.genes.plot.text <- renderText({
    if(input$save.mv.genes.plot==0)
      return(NULL)
    else HTML(paste0('Figure is saved in: ',DDD.data$figure.folder))
  })
  
  
  ##--------------Step 3: PCA plot---------------
  ###update the pca selections in the selectinput
  
  observe({
    if(is.null(input$sample.input)){
      x <- 5
    } else {
      if(input$pca.plot.type=='Bio-reps'){
        x <- length(DDD.data$samples_new$condition)
      } else {
        x <- length(unique(DDD.data$samples_new$condition))
      }
    }
    
    x <- paste0('PC',1:x)
    updateSelectInput(session, "dim1",
                      label = 'x-axis',
                      choices = x,
                      selected = 'PC1'
    )
    updateSelectInput(session, "dim2",
                      label = 'y-axis',
                      choices = x,
                      selected = 'PC2'
    )
  })
  
  ###---pca color options
  output$pca.color.option.ui <- renderUI({
    span(
      selectInput(inputId = 'pca.color.option',
                  label = 'Colour on groups',
                  selected = 'brep',multiple = T,
                  choices = colnames(DDD.data$samples_new))
    )
  })
  
  
  ###---trans level
  pca.trans.g <- eventReactive(input$plot.pca.button,{
    data2pca <- DDD.data$trans_counts[DDD.data$target_high$trans_high,]
    dge <- DGEList(counts=data2pca)
    dge <- calcNormFactors(dge)
    data2pca <- t(counts2CPM(obj = dge,Log = T))
    dim1 <- input$dim1
    dim2 <- input$dim2
    ellipse.type <- input$pca.add.circle
    
    ##--Bio-reps
    if(input$pca.plot.type=='Bio-reps'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      # if(input$pca.color.option=='Bio-reps')
      #   groups <- DDD.data$samples_new$brep
      # if(input$pca.color.option=='Conditions')
      #   groups <- DDD.data$samples_new$condition
      groups <- interaction(DDD.data$samples_new[,input$pca.color.option])
      g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                      groups = groups,plot.title = 'Transcript PCA: bio-reps',
                      ellipse.type = ellipse.type,
                      add.label = T,adj.label = F)
    }
    
    ##--average expression plot
    if(input$pca.plot.type=='Average expression'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      groups <- DDD.data$samples_new$brep
      data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
      groups <- unique(DDD.data$samples_new$condition)
      g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                      groups = groups,plot.title = 'Transcript PCA: average expression',
                      ellipse.type = 'none',add.label = T,adj.label = T)
    }
    g
  })
  
  output$trans.pca.plot <- renderPlotly({
    ggplotly(pca.trans.g())
  })
  
  ###---genes level
  pca.genes.g <- eventReactive(input$plot.pca.button,{
    data2pca <- DDD.data$genes_counts[DDD.data$target_high$genes_high,]
    dge <- DGEList(counts=data2pca)
    dge <- calcNormFactors(dge)
    data2pca <- t(counts2CPM(obj = dge,Log = T))
    dim1 <- input$dim1
    dim2 <- input$dim2
    ellipse.type <- input$pca.add.circle
    
    ##--Bio-reps
    if(input$pca.plot.type=='Bio-reps'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      # if(input$pca.color.option=='Bio-reps')
      #   groups <- DDD.data$samples_new$brep
      # if(input$pca.color.option=='Conditions')
      #   groups <- DDD.data$samples_new$condition
      groups <- interaction(DDD.data$samples_new[,input$pca.color.option])
      g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                      groups = groups,plot.title = 'Gene PCA: bio-reps',
                      ellipse.type = ellipse.type,
                      add.label = T,adj.label = F)
    }
    
    ##--average expression plot
    if(input$pca.plot.type=='Average expression'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      groups <- DDD.data$samples_new$brep
      data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
      groups <- unique(DDD.data$samples_new$condition)
      g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                      groups = groups,plot.title = 'Gene PCA: average expression',
                      ellipse.type = 'none',add.label = T,adj.label = T)
    }
    g
  })
  
  output$genes.pca.plot <- renderPlotly({
    ggplotly(pca.genes.g())
  })
  
  observeEvent(input$save.pca.plot,{
    ##transcript level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA ',input$pca.plot.type,'.png'),
        width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
    print(pca.trans.g())
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA ',input$pca.plot.type,'.pdf'),
        width = input$pca.plot.width,height = input$pca.plot.height)
    print(pca.trans.g())
    dev.off()
    ###gene level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Gene PCA ',input$pca.plot.type,'.png'),
        width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
    print(pca.genes.g())
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA ',input$pca.plot.type,'.pdf'),
        width = input$pca.plot.width,height = input$pca.plot.height)
    print(pca.genes.g())
    dev.off()
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  })
  
  
  ##--------------remove batch effects-trans level---------------
  ###---update remove.batch action button
  
  ###---transcript level
  observe({
    if(is.null(DDD.data$samples_new)){
      max.n <- 6
    } else {
      max.n <- nrow(DDD.data$samples_new)
    }
    updateSliderInput(session, "ruvseq.k", max=max.n,min = 1,step = 1)
  })
  
  
  observeEvent(input$remove.batch,{
    cat('\nEstimate transcript level batch effects...')
    trans_batch <- remove.batch.shiny(read.counts = DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                                      condition = DDD.data$samples_new$condition,
                                      design = NULL,
                                      contrast=NULL,
                                      group = DDD.data$samples_new$condition,
                                      method = input$ruvseq.method,
                                      k = input$ruvseq.k)
    save(trans_batch,file=paste0(DDD.data$data.folder,'/trans_batch.RData'))
    message(paste0('trans_batch.RData is saved in folder: ',DDD.data$data.folder))
    DDD.data$trans_batch <- trans_batch
  })
  
  ###---gene level
  observeEvent(input$remove.batch,{
    cat('\nEstimate gene level batch effects')
    genes_batch <- remove.batch.shiny(read.counts = DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                                      condition = DDD.data$samples_new$condition,
                                      design = NULL,
                                      contrast=NULL,
                                      group = DDD.data$samples_new$condition,
                                      method = input$ruvseq.method,
                                      k = input$ruvseq.k)
    save(genes_batch,file=paste0(DDD.data$data.folder,'/genes_batch.RData'))
    message(paste0('genes_batch.RData is saved in folder: ',DDD.data$data.folder))
    DDD.data$genes_batch <- genes_batch
  })
  
  br.saved.info <- reactive({
    if(is.null(DDD.data$data.folder))
      return(NULL)
    if(input$remove.batch>0){
      paste0('trans(genes)_batch.RData are saved in the data folder.')
    } else {
      if(!is.null(DDD.data$genes_batch) & !is.null(DDD.data$trans_batch)){
        paste0('trans(genes)_batch.RData are already uploaded in the workspace.')
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/genes_batch.RData')) &
           file.exists(paste0(DDD.data$data.folder,'/trans_batch.RData'))){
          paste0('trans(genes)_batch.RData are already in "data" folder.\nCan upload directly in "Data generation".')
        } else {
          paste0('trans(genes)_batch.RData are not in "data" folder.\nPlease generate new.')
        }
      }
    }
  })
  
  output$br.saved.info <- renderText({
    if(is.null(br.saved.info()))
      return(NULL)
    br.saved.info()
  })
  #
  ##--------------Step 3-continued: PCA after remove batch---------------
  # ###---trans level
  pca.trans.br.g <- eventReactive(input$update.pca.plot,{
    message('Transcript level PCA--batch removed')
    data2pca <- DDD.data$trans_batch$normalizedCounts
    dge <- DGEList(counts=data2pca)
    dge <- suppressWarnings(calcNormFactors(dge))
    data2pca <- t(counts2CPM(obj = dge,Log = T))
    dim1 <- input$dim1
    dim2 <- input$dim2
    ellipse.type <- input$pca.add.circle
    
    ##--Bio-reps
    if(input$pca.plot.type=='Bio-reps'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      # if(input$pca.color.option=='Bio-reps')
      #   groups <- DDD.data$samples_new$brep
      # if(input$pca.color.option=='Conditions')
      #   groups <- DDD.data$samples_new$condition
      groups <- interaction(DDD.data$samples_new[,input$pca.color.option])
      g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                      groups = groups,plot.title = 'Transcript PCA: bio-reps',
                      ellipse.type = ellipse.type,
                      add.label = T,adj.label = F)
    }
    
    ##--average expression plot
    if(input$pca.plot.type=='Average expression'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      groups <- DDD.data$samples_new$brep
      data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
      groups <- unique(DDD.data$samples_new$condition)
      g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                      groups = groups,plot.title = 'Transcript PCA: average expression',
                      ellipse.type = 'none',add.label = T,adj.label = T)
    }
    g
  })
  #
  
  output$trans.pca.br.plot <- renderPlotly({
    ggplotly(pca.trans.br.g())
  })
  
  
  # ###---genes level
  pca.genes.br.g <- eventReactive(input$update.pca.plot,{
    message('Gene level PCA--batch removed')
    data2pca <- DDD.data$genes_batch$normalizedCounts
    dge <- DGEList(counts=data2pca)
    dge <- suppressWarnings(calcNormFactors(dge))
    data2pca <- t(counts2CPM(obj = dge,Log = T))
    dim1 <- input$dim1
    dim2 <- input$dim2
    ellipse.type <- input$pca.add.circle
    
    ##--Bio-reps
    if(input$pca.plot.type=='Bio-reps'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      # if(input$pca.color.option=='Bio-reps')
      #   groups <- DDD.data$samples_new$brep
      # if(input$pca.color.option=='Conditions')
      #   groups <- DDD.data$samples_new$condition
      groups <- interaction(DDD.data$samples_new[,input$pca.color.option])
      g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                      groups = groups,plot.title = 'Gene PCA: bio-reps',
                      ellipse.type = ellipse.type,
                      add.label = T,adj.label = F)
    }
    
    ##--average expression plot
    if(input$pca.plot.type=='Average expression'){
      rownames(data2pca) <- gsub('_','.',rownames(data2pca))
      groups <- DDD.data$samples_new$brep
      data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
      groups <- unique(DDD.data$samples_new$condition)
      g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                      groups = groups,plot.title = 'Gene PCA: average expression',
                      ellipse.type = 'none',add.label = T,adj.label = T)
    }
    g
  })
  #
  
  output$genes.pca.br.plot <- renderPlotly({
    ggplotly(pca.genes.br.g())
  })
  
  observeEvent(input$save.pca.br.plot,{
    ###transcript level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',input$pca.plot.type,'.png'),
        width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
    print(pca.trans.br.g())
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',input$pca.plot.type,'.pdf'),
        width = input$pca.plot.width,height = input$pca.plot.height)
    print(pca.trans.br.g())
    dev.off()
    
    ###gene level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',input$pca.plot.type,'.png'),
        width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
    print(pca.genes.br.g())
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',input$pca.plot.type,'.pdf'),
        width = input$pca.plot.width,height = input$pca.plot.height)
    print(pca.genes.br.g())
    dev.off()
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  })
  
  ##--------------Step 4-normalization--------------
  observeEvent(input$run.norm,{
    cat('\nNormalise transcript read counts')
    dge <- DGEList(counts=DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                   group = DDD.data$samples_new$condition,
                   genes = DDD.data$mapping[DDD.data$target_high$trans_high,])
    trans_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm.method))
    save(trans_dge,file=paste0(DDD.data$data.folder,'/trans_dge.RData'))
    message(paste0('trans_dge.RData is saved in folder: ',DDD.data$data.folder))
    message('Done!!!')
    DDD.data$trans_dge <- trans_dge
  })
  
  observeEvent(input$run.norm,{
    cat('\nNormalise gene read counts')
    dge <- DGEList(counts=DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                   group = DDD.data$samples_new$condition)
    genes_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm.method))
    save(genes_dge,file=paste0(DDD.data$data.folder,'/genes_dge.RData'))
    message(paste0('genes_dge.RData is saved in folder: ',DDD.data$data.folder))
    message('Done!!!')
    DDD.data$genes_dge <- genes_dge
  })
  
  
  norm.info <- reactive({
    if(is.null(DDD.data$data.folder))
      return(NULL)
    if(input$run.norm>0){
      paste0('trans(genes)_dge.RData are saved in the data folder.')
    } else {
      if(!is.null(DDD.data$genes_dge) & !is.null(DDD.data$trans_dge)){
        paste0('trans(genes)_dge.RData are already uploaded in the workspace.')
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/genes_dge.RData')) &
           file.exists(paste0(DDD.data$data.folder,'/trans_dge.RData'))){
          paste0('trans(genes)_dge.RData are already in "data" folder.\nCan upload directly in "Data generation".')
        } else {
          paste0('trans(genes)_dge.RData are not in "data" folder.\nPlease generate new.')
        }
      }
    }
  })
  
  
  output$norm.info <- renderText({
    if(is.null(norm.info()))
      return(NULL)
    norm.info()
  })
  
  
  ##--------------make the distribution plot-------------
  ###---transcript level
  trans.dist.g <- eventReactive(input$norm.plot.button,{
    sample.name <- paste0(DDD.data$samples_new$condition,'.',DDD.data$samples_new$brep)
    condition <- DDD.data$samples_new$condition
    data.before <- DDD.data$trans_counts[DDD.data$target_high$trans_high,]
    data.after <- counts2CPM(obj = DDD.data$trans_dge,Log = T)
    g <- boxplot.normalised(data.before = data.before,
                            data.after = data.after,
                            condition = condition,
                            sample.name = sample.name)
    g
  })
  
  output$trans.dist.plot <- renderPlotly({
    q <- subplot(style(trans.dist.g()$g1+theme(legend.position = 'none',
                                               axis.title.x = element_blank(),
                                               axis.text.x = element_blank()), showlegend = FALSE),
                 trans.dist.g()$g2+labs(title='Data distribution (before and after normalization)'),
                 nrows = 2,
                 titleX=T,
                 titleY = T,margin=0.04)
    q
  })
  
  
  ###---gene level
  genes.dist.g <- eventReactive(input$norm.plot.button,{
    sample.name <- paste0(DDD.data$samples_new$condition,'.',DDD.data$samples_new$brep)
    condition <- DDD.data$samples_new$condition
    data.before <- DDD.data$genes_counts[DDD.data$target_high$genes_high,]
    data.after <- counts2CPM(obj = DDD.data$genes_dge,Log = T)
    g <- boxplot.normalised(data.before = data.before,
                            data.after = data.after,
                            condition = condition,
                            sample.name = sample.name)
    g
  })
  
  output$genes.dist.plot <- renderPlotly({
    q <- subplot(style(genes.dist.g()$g1+theme(legend.position = 'none',
                                               axis.title.x = element_blank(),
                                               axis.text.x = element_blank()), showlegend = FALSE),
                 genes.dist.g()$g2+labs(title='Data distribution (before and after normalization)'),
                 nrows = 2,
                 titleX=T,
                 titleY = T,margin=0.04)
    q
  })
  
  observeEvent(input$save.dist.plot,{
    ##transcript level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Transcript data distribution.png'),
        width = input$dist.plot.width,height = input$dist.plot.height,res=input$dist.plot.res, units = 'in')
    gridExtra::grid.arrange(trans.dist.g()$g1,trans.dist.g()$g2,ncol=1)
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Transcript data distribution.pdf'),
        width = input$dist.plot.width,height = input$dist.plot.height)
    gridExtra::grid.arrange(trans.dist.g()$g1,trans.dist.g()$g2,ncol=1)
    dev.off()
    
    ##gene level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Gene data distribution.png'),
        width = input$dist.plot.width,height = input$dist.plot.height,res=input$dist.plot.res, units = 'in')
    gridExtra::grid.arrange(genes.dist.g()$g1,genes.dist.g()$g2,ncol=1)
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Gene data distribution.pdf'),
        width = input$dist.plot.width,height = input$dist.plot.height)
    gridExtra::grid.arrange(genes.dist.g()$g1,genes.dist.g()$g2,ncol=1)
    dev.off()
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  })
  
  ##---------------data information-----------------
  
  data.info <- observe({
    if(is.null(DDD.data$samples) | is.null(DDD.data$trans_counts) | is.null(DDD.data$genes_counts))
      return(NULL)
    data.info <- data.frame(
      raw.number=c(length(DDD.data$mapping$TXNAME),length(unique(DDD.data$mapping$GENEID))),
      sample=rep(nrow(DDD.data$samples),2),
      condition=rep(length(unique(DDD.data$samples$condition)),2),
      brep=rep(length(unique(DDD.data$samples$brep)),2),
      srep=rep(length(unique(DDD.data$samples$srep)),2)
    )
    
    if(!is.null(DDD.data$trans_counts))
      data.info$merged.sample <- c(ncol(DDD.data$trans_counts),ncol(DDD.data$genes_counts))
    
    if(!is.null(DDD.data$target_high)){
      data.info$cpm.cut <- c(input$cpm.cut,'--')
      data.info$sample.n.cut <- c(input$sample.n.cut,'--')
      data.info$after.filter <- c(length(DDD.data$target_high$trans_high),length(DDD.data$target_high$genes_high))
    }
    rownames(data.info) <- c('Transcript','Gene')
    DDD.data$data.info <- data.info
  })
  
  output$data.info <- DT::renderDataTable({
    DDD.data$data.info
  },
  options = list(
    scrollX=T,
    searching = FALSE,
    paging=F,
    info = FALSE,
    columnDefs = list(list(className = 'dt-center', targets="_all")))
  )
  
  observeEvent(input$save.data.info,{
    write.csv(DDD.data$data.info,file=paste0(DDD.data$result.folder,'/data.info.csv'),row.names = T)
    message(paste0('Table is saved in folder: ',DDD.data$result.folder))
    showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
  })
  
  ##--------------->> DE DAS and DTU <<---------------
  ##---------------- Load contrast ----------------
  observe({
    inFile <- input$contrast.input
    if (is.null(inFile))
      return(NULL)
    contrast <- read.csv(inFile$datapath, header = T)
    contrast <- paste0(contrast[,1],'-',contrast[,2])
    save(contrast,file=paste0(DDD.data$data.folder,'/contrast.RData'))
    DDD.data$contrast <- contrast
    
  })
  
  output$contrast.table<- DT::renderDataTable({
    if(is.null(DDD.data$contrast))
      return(NULL)
    x <- do.call(rbind,strsplit(DDD.data$contrast,'-'))
    colnames(x) <- c('Treatmeant','Control')
    rownames(x) <- paste0('Contrast group',1:nrow(x))
    x
  },options = list(
    columnDefs = list(list(className = 'dt-center', targets = "_all")))
  )
  
  ##--------------- DE genes -----------------
  observeEvent(input$run.DE,{
    cat('\nDE gene analysis')
    withProgress(message = 'DE gene analysis...', detail = 'This may take a while...',
                 value = 0, {
                   if(is.null(DDD.data$genes_batch)){
                     batch.effect <- NULL
                   } else {
                     batch.effect <- DDD.data$genes_batch$W}
                   incProgress(0.1)
                   design <- condition2design(condition = DDD.data$samples_new$condition,
                                              batch.effect = batch.effect)
                   DE.pipeline <- input$DE.pipeline
                   incProgress(0.2)
                   switch(DE.pipeline,
                          limma={
                            genes_3D_stat <- limma.pipeline(dge = DDD.data$genes_dge,
                                                            design = design,
                                                            deltaPS = DDD.data$deltaPS,
                                                            contrast = DDD.data$contrast,
                                                            diffAS = F,
                                                            adjust.method = input$p.adjust.method)
                          },
                          glmQL={
                            genes_3D_stat <- edgeR.pipeline(dge = DDD.data$genes_dge,
                                                            design = design,
                                                            deltaPS = DDD.data$deltaPS,
                                                            contrast = DDD.data$contrast,
                                                            diffAS = F,
                                                            method = 'glmQL',
                                                            adjust.method = input$p.adjust.method)
                          },
                          glm={
                            genes_3D_stat <- edgeR.pipeline(dge = DDD.data$genes_dge,
                                                            design = design,
                                                            deltaPS = NULL,
                                                            contrast = DDD.data$contrast,
                                                            diffAS = F,
                                                            method = 'glm',
                                                            adjust.method = input$p.adjust.method)
                          }
                   )
                   incProgress(0.7)
                   DE_genes <- summaryDEtarget(stat = genes_3D_stat$DE.stat,
                                               cutoff = c(adj.pval=input$pval.cutoff,
                                                          log2FC=input$lfc.cutoff))
                   
                   DDD.data$DE_genes <- DE_genes
                   DDD.data$genes_log2FC <- genes_3D_stat$DE.lfc
                   save(DE_genes,file=paste0(DDD.data$data.folder,'/DE_genes.RData'))
                   incProgress(0.8)
                   save(genes_3D_stat,file=paste0(DDD.data$data.folder,'/genes_3D_stat.RData'))
                   message(paste0('genes_3D_stat.RData is saved in folder: ',DDD.data$data.folder))
                 })
  })
  
  
  ###run DE information
  run.DE.info <- reactive({
    if(is.null(DDD.data$data.folder))
      return(NULL)
    if(input$run.DE>0){
      paste0('genes_3D_stat.RData is saved in the data folder.')
    } else {
      if(file.exists(paste0(DDD.data$data.folder,'/genes_3D_stat.RData'))){
        paste0('genes_3D_stat.RData is already in "data" folder.\nGenerate new?.')
      } else {
        paste0('genes_3D_stat.RData is not in "data" folder.\nPlease generate new.')
      }
    }
  })
  
  output$run.DE.info <- renderText({
    if(is.null(run.DE.info()))
      return(NULL)
    run.DE.info()
  })
  #
  
  ##--------------- generate deltaPS -----------------
  observeEvent(input$run.deltaPS,{
    if(is.null(DDD.data$PS)){
      ##load require data for analysis
      if(is.null(DDD.data$target_high)){
        load(paste0(DDD.data$data.folder,'/target_high.RData'))
        DDD.data$target_high <- target_high
      }
      
      if(is.null(DDD.data$mapping)){
        load(paste0(DDD.data$data.folder,'/mapping.RData'))
        DDD.data$mapping <- mapping
      }
      
      if(input$PS.data.type=='TPM'){
        if(is.null(DDD.data$samples)){
          load(paste0(DDD.data$data.folder,'/samples.RData'))
          DDD.data$samples <- samples
        }
        
        if(is.null(DDD.data$trans_TPM)){
          load(paste0(DDD.data$data.folder,'/trans_TPM.RData'))
          DDD.data$trans_TPM <- trans_TPM
        }
        transAbundance <- DDD.data$trans_TPM
        condition = DDD.data$samples$condition
      }
      if(input$PS.data.type=='Read counts'){
        if(is.null(DDD.data$samples_new)){
          load(paste0(DDD.data$data.folder,'/samples_new.RData'))
          DDD.data$samples_new <- samples_new
        }
        if(is.null(DDD.data$trans_counts)){
          load(paste0(DDD.data$data.folder,'/trans_counts.RData'))
          DDD.data$trans_counts <- trans_counts
        }
        transAbundance <- DDD.data$trans_counts
        condition = DDD.data$samples_new$condition
      }
      transAbundance <- transAbundance[DDD.data$target_high$trans_high,]
    } else {
      transAbundance <- NULL
      condition <- NULL
    }
    
    cat('\nGenerating deltaPS')
    deltaPS <- transAbundance2PS(transAbundance = transAbundance,
                                 PS = DDD.data$PS,
                                 contrast = DDD.data$contrast,
                                 condition = condition,
                                 mapping = DDD.data$mapping[DDD.data$target_high$trans_high,])
    PS <- deltaPS$PS
    deltaPS <- deltaPS$deltaPS
    save(deltaPS,file=paste0(DDD.data$data.folder,'/deltaPS.RData'))
    save(PS,file=paste0(DDD.data$data.folder,'/PS.RData'))
    message(paste0('deltaPS.RData are saved in folder: ',DDD.data$data.folder))
    showNotification(paste0('deltaPS.RData are saved in folder: ',DDD.data$data.folder))
    message('Done!!!')
    DDD.data$PS <- PS
    DDD.data$deltaPS <- deltaPS
  })
  
  output$deltaPS.table <- DT::renderDataTable({
    if(is.null(DDD.data$deltaPS))
      return(NULL)
    round(DDD.data$deltaPS[1:20,,drop=F],5)
  },options = list(
    columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  
  observeEvent(input$save.deltaPS.table,{
    write.csv(DDD.data$deltaPS,file=paste0(DDD.data$result.folder,'/deltaPS.csv'),row.names = T)
    message(paste0('Table is saved in folder: ',DDD.data$result.folder))
    showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
  })
  
  ### message of run deltaPS
  rundeltaPSinfo <- reactive({
    if(is.null(DDD.data$data.folder))
      return(NULL)
    if(input$run.deltaPS>0){
      paste0('deltaPS.RData is saved in the data folder.')
    } else {
      if(!is.null(DDD.data$deltaPS)){
        paste0('deltaPS.RData is already uploaded in the workspace.')
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/deltaPS.RData'))){
          paste0('deltaPS.RData is already in "data" folder.\nCan upload directly in "Data generation".')
        } else {
          paste0('deltaPS.RData is not in "data" folder.\nPlease generate new.')
        }
      }
    }
  })
  
  output$rundeltaPSinfo <- renderText({
    if(is.null(rundeltaPSinfo()))
      return(NULL)
    rundeltaPSinfo()
  })
  
  ##--------------- DE DAS DTU transcripts ---------------
  observeEvent(input$run.diffsplice,{
    if(is.null(DDD.data$deltaPS)){
      message('Please generate deltaPS values.')
      showNotification('Please generate deltaPS values.')
      return(NULL)
    }
    withProgress(message = 'DE gene analysis...', detail = 'This may take a while...',
                 value = 0, {
                   cat('\nDAS gene, DE and DTU transcript analysis')
                   if(is.null(DDD.data$trans_batch)){
                     batch.effect <- NULL
                   } else {
                     batch.effect <- DDD.data$trans_batch$W
                   }
                   incProgress(0.1)
                   design <- condition2design(condition = DDD.data$samples_new$condition,
                                              batch.effect = batch.effect)
                   DE.pipeline <- input$DE.pipeline
                   incProgress(0.2)
                   switch(DE.pipeline,
                          limma={
                            trans_3D_stat <- limma.pipeline(dge = DDD.data$trans_dge,
                                                            design = design,
                                                            deltaPS = DDD.data$deltaPS,
                                                            contrast = DDD.data$contrast,
                                                            diffAS = T,
                                                            adjust.method = input$p.adjust.method)
                          },
                          glmQL={
                            trans_3D_stat <- edgeR.pipeline(dge = DDD.data$trans_dge,
                                                            design = design,
                                                            deltaPS = DDD.data$deltaPS,
                                                            contrast = DDD.data$contrast,
                                                            diffAS = T,
                                                            method = 'glmQL',
                                                            adjust.method = input$p.adjust.method)
                          },
                          glm={
                            trans_3D_stat <- edgeR.pipeline(dge = DDD.data$trans_dge,
                                                            design = design,
                                                            deltaPS = DDD.data$deltaPS,
                                                            contrast = DDD.data$contrast,
                                                            diffAS = T,
                                                            method = 'glm',
                                                            adjust.method = input$p.adjust.method)
                          }
                   )
                   incProgress(0.7)
                   ##DE trans
                   DE_trans <- summaryDEtarget(stat = trans_3D_stat$DE.stat,
                                               cutoff = c(adj.pval=input$pval.cutoff,
                                                          log2FC=input$lfc.cutoff))
                   DDD.data$DE_trans <- DE_trans
                   DDD.data$trans_log2FC <- trans_3D_stat$DE.lfc
                   save(DE_trans,file=paste0(DDD.data$data.folder,'/DE_trans.RData'))
                   
                   ##DAS genes
                   if(input$DAS.p.method=='F-test')
                     DAS.stat <- trans_3D_stat$DAS.F.stat else DAS.stat <- trans_3D_stat$DAS.Simes.stat
                   
                   lfc <- DDD.data$genes_log2FC
                   lfc <- reshape2::melt(as.matrix(lfc))
                   colnames(lfc) <- c('target','contrast','log2FC')
                   DAS_genes <- summaryDAStarget(stat = DAS.stat,
                                                 lfc = lfc,
                                                 cutoff=c(input$pval.cutoff,input$deltaPS.cutoff))
                   DDD.data$DAS_genes <- DAS_genes
                   save(DAS_genes,file=paste0(DDD.data$data.folder,'/DAS_genes.RData'))
                   
                   ##DTU  trans
                   lfc <- DDD.data$trans_log2FC
                   lfc <- reshape2::melt(as.matrix(lfc))
                   colnames(lfc) <- c('target','contrast','log2FC')
                   DTU_trans <- summaryDAStarget(stat = trans_3D_stat$DTU.stat,
                                                 lfc = lfc,cutoff = c(adj.pval=input$pval.cutoff,
                                                                      deltaPS=input$deltaPS.cutoff))
                   DDD.data$DTU_trans <- DTU_trans
                   save(DTU_trans,file=paste0(DDD.data$data.folder,'/DTU_trans.RData'))
                   incProgress(0.8)
                   
                   save(trans_3D_stat,file=paste0(DDD.data$data.folder,'/trans_3D_stat.RData'))
                   DDD.data$trans_3D_stat <- trans_3D_stat
                   message(paste0('trans_3D_stat.RData is saved in folder: ',DDD.data$data.folder))
                 })
  })
  
  diffsplice.method <- reactive({
    if(input$DE.pipeline=='limma')
      'Use limma::diffSplice' else 'Use edgeR::diffSpliceDGE'
  })
  
  output$diffsplice.method <- renderText({
    diffsplice.method()
  })
  
  ###run DAS information
  run.DAS.info <- reactive({
    if(is.null(DDD.data$data.folder))
      return(NULL)
    if(input$run.DE>0){
      paste0('trans_3D_stat.RData is saved in the data folder.')
    } else {
      if(file.exists(paste0(DDD.data$data.folder,'/trans_3D_stat.RData'))){
        paste0('trans_3D_stat.RData is already in "data" folder.\nGenerate new?.')
      } else {
        paste0('trans_3D_stat.RData is not in "data" folder.\nPlease generate new.')
      }
    }
  })
  
  output$run.DAS.info <- renderText({
    if(is.null(run.DAS.info()))
      return(NULL)
    run.DAS.info()
  })
  
  ##---------------->>  Result summary   <<-----------------
  ##---------------- DE DAS DTU statistics -----------------
  
  output$show.contrast.groups <- renderUI({
    selectInput(inputId = 'contrast.groups',label = 'Contrast groups',
                choices = DDD.data$contrast)
  })
  
  output$show.top.stat.table <- renderUI({
    box(title=paste0('Top ',input$top.stat, ' statistics of ',input$stat.type,'; Contrast: ',input$contrast.groups),
        width = 13,status = 'primary', solidHeader = T,
        DT::dataTableOutput('top.stat.table')
    )
  })
  
  
  ###--show the table
  DDD.stat <- reactive({
    if(is.null(input$contrast.groups))
      return(NULL)
    
    if(input$stat.type=='DE genes'){
      if(is.null(DDD.data$DE_genes))
        return(NULL)
      x <- DDD.data$DE_genes
    }
    
    if(input$stat.type=='DAS genes'){
      if(is.null(DDD.data$DAS_genes))
        return(NULL)
      x <- DDD.data$DAS_genes
    }
    
    if(input$stat.type=='DE transcripts'){
      if(is.null(DDD.data$DE_trans))
        return(NULL)
      x <- DDD.data$DE_trans
    }
    
    if(input$stat.type=='DTU transcripts'){
      if(is.null(DDD.data$DTU_trans))
        return(NULL)
      x <- DDD.data$DTU_trans
    }
    
    x <- split(x,x$contrast)
    x <- lapply(x,function(i) i[order(i$adj.pval,decreasing = F),])
    
    stat <- x[[input$contrast.groups]][1:input$top.stat,]
    stat <- data.frame(lapply(stat, function(y) if(is.numeric(y)) format(y,digits=5) else y))
    return(stat)
  })
  
  
  output$top.stat.table <- DT::renderDataTable({
    if(is.null(DDD.stat()))
      return(NULL)
    DDD.stat()
  },options = list(
    scrollX=T,
    columnDefs = list(list(className = 'dt-center',
                           targets = "_all"))))
  
  ##---------------DDD bar plot-----------------
  # output$up.down.bar.plot <- renderPlotly({
  g.updown <- reactive({
    if(is.null((DDD.data$DE_genes)) | is.null((DDD.data$DAS_genes)) | is.null((DDD.data$DE_trans)) | is.null((DDD.data$DTU_trans)))
      return(NULL)
    idx <- c('DE_genes','DAS_genes','DE_trans','DTU_trans')
    title.idx <- gsub('_',' ',idx)
    title.idx <- gsub('trans','transcripts',title.idx)
    names(title.idx) <- idx
    g <- lapply(idx,function(i){
      idx <- factor(DDD.data[[i]]$contrast,levels = DDD.data$contrast)
      targets <-  split(DDD.data[[i]],idx)
      data2plot <- lapply(DDD.data$contrast,function(i){
        if(nrow(targets[[i]])==0){
          x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
        } else {
          x <- data.frame(contrast=i,table(targets[[i]]$up.down))
          colnames(x) <- c('contrast','regulation','number')
        }
        x
      })
      data2plot <- do.call(rbind,data2plot)
      plotUpdown(data2plot,plot.title = title.idx[i],
                 contrast = DDD.data$contrast,
                 angle = input$updown.x.rotate,hjust = input$updown.x.hjust)
    })
    names(g) <- title.idx
    g
  })
  
  output$up.down.bar.plot.DE.genes <- renderPlot({
    if(is.null((g.updown())))
      return(NULL)
    print(g.updown()[[1]])
  })
  
  output$up.down.bar.plot.DAS.genes <- renderPlot({
    if(is.null((g.updown())))
      return(NULL)
    print(g.updown()[[2]])
  })
  
  output$up.down.bar.plot.DE.trans <- renderPlot({
    if(is.null((g.updown())))
      return(NULL)
    print(g.updown()[[3]])
  })
  
  output$up.down.bar.plot.DTU.trans <- renderPlot({
    if(is.null((g.updown())))
      return(NULL)
    print(g.updown()[[4]])
  })
  
  observeEvent(input$save.up.down.bar.plot,{
    lapply(names(g.updown()),function(i){
      png(paste0(DDD.data$figure.folder,'/',i,' updown regulation numbers.png'),
          width = input$updown.plot.width,height = input$updown.plot.height,units = 'in',
          res = input$updown.plot.res)
      print(g.updown()[[i]])
      dev.off()
      
      pdf(paste0(DDD.data$figure.folder,'/',i,' updown regulation numbers.pdf'),
          width = input$updown.plot.width,height = input$updown.plot.height)
      print(g.updown()[[i]])
      dev.off()
    })
    
    message(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
  })
  
  ##---------------Comparisons between contrast groups-----------------
  output$across.contrast.euler <- renderUI({
    span(
      selectInput(inputId = 'across.contrast.group',
                  label = 'Contrast groups (multiple selection allowed)',
                  selected = DDD.data$contrast[1:2],multiple = T,
                  choices = DDD.data$contrast)
    )
  })
  
  output$across.target.euler <- renderUI({
    span(
      selectInput(inputId = 'across.target',
                  label = 'Contrast groups',
                  selected = DDD.data$contrast[1],multiple = F,
                  choices = DDD.data$contrast)
    )
  })
  
  
  ##---------------flow chat-----------------
  genes.flow.chart <- reactive({
    DE.genes <- unique(DDD.data$DE_genes$target)
    DAS.genes <- unique(DDD.data$DAS_genes$target)
    genes.flow.chart <- function(){
      plotFlowChart(expressed = DDD.data$target_high$genes_high,
                    x = DE.genes,
                    y = DAS.genes,
                    type = 'genes',
                    pval.cutoff = input$pval.cutoff,lfc.cutoff = input$lfc.cutoff,deltaPS.cutoff = input$deltaPS.cutoff)
    }
    genes.flow.chart
  })
  output$genes.flow.chart <- renderPlot({
    genes.flow.chart()()
  })
  
  trans.flow.chart <- reactive({
    DE.trans<- unique(DDD.data$DE_trans$target)
    DTU.trans <- unique(DDD.data$DTU_trans$target)
    trans.flow.chart <- function(){
      plotFlowChart(expressed = DDD.data$target_high$trans_high,
                    x = DE.trans,
                    y = DTU.trans,
                    type = 'transcripts',
                    pval.cutoff = input$pval.cutoff,lfc.cutoff = input$lfc.cutoff,deltaPS.cutoff = input$deltaPS.cutoff)
    }
    trans.flow.chart
  })
  output$trans.flow.chart <- renderPlot({
    trans.flow.chart()()
  })
  
  observeEvent(input$save.union.flow.chart,{
    ##transcript level
    png(filename = paste0(DDD.data$figure.folder,'/Union set DE genes vs DAS genes.png'),
        width = input$flow.plot.width,height = input$flow.plot.height,
        res=input$flow.plot.res, units = 'in')
    genes.flow.chart()()
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Union set DE genes vs DAS genes.pdf'),
        width = input$flow.plot.width,height = input$flow.plot.height)
    genes.flow.chart()()
    dev.off()
    
    ###gene level
    # graphics.off()
    png(filename = paste0(DDD.data$figure.folder,'/Union set DE transcripts vs DTU transcripts.png'),
        width = input$flow.plot.width,height = input$flow.plot.height,
        res=input$flow.plot.res, units = 'in')
    trans.flow.chart()()
    dev.off()
    
    pdf(file = paste0(DDD.data$figure.folder,'/Union set DE transcripts vs DTU transcripts.pdf'),
        width = input$flow.plot.width,height = input$flow.plot.height)
    trans.flow.chart()()
    dev.off()
    
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  })
  
  
  
  ###---DE genes
  g.across.contrast <- reactive({
    if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans) | is.null(DDD.data$DTU_trans) | is.null(input$across.contrast.group))
      return(NULL)
    
    idx <- c('DE_genes','DAS_genes','DE_trans','DTU_trans')
    title.idx <- gsub('_',' ',idx)
    title.idx <- gsub('trans','transcripts',title.idx)
    g.list <- lapply(idx,function(i){
      targets <- lapply(input$across.contrast.group,function(j){
        subset(DDD.data[[i]],contrast==j)$target
      })
      names(targets) <- input$across.contrast.group
      g <- plotEulerDiagram(x = targets,scale.plot = input$across.contrast.euler.plot.scale)
      g
    })
    names(g.list) <- title.idx
    g.list
  })
  
  ###---DE genes
  output$DE.genes.euler <- renderPlot({
    if(is.null(g.across.contrast()))
      return(NULL)
    grid.arrange(g.across.contrast()[[1]],top=textGrob('DE genes', gp=gpar(cex=1.2)))
  })
  
  ###---DAS genes
  output$DAS.genes.euler <- renderPlot({
    if(is.null(g.across.contrast()))
      return(NULL)
    grid.arrange(g.across.contrast()[[2]],top=textGrob('DAS genes', gp=gpar(cex=1.2)))
  })
  
  ###---DE transcripts
  output$DE.trans.euler <- renderPlot({
    if(is.null(g.across.contrast()))
      return(NULL)
    grid.arrange(g.across.contrast()[[3]],top=textGrob('DAS genes', gp=gpar(cex=1.2)))
  })
  
  ###---DTU transcripts
  output$DTU.trans.euler <- renderPlot({
    if(is.null(g.across.contrast()))
      return(NULL)
    grid.arrange(g.across.contrast()[[4]],top=textGrob('DAS genes', gp=gpar(cex=1.2)))
  })
  
  observeEvent(input$save.across.contrast.euler.plot,{
    ###save DE.genes
    
    lapply(names(g.across.contrast()),function(i){
      figure.name <- paste0(DDD.data$figure.folder,'/',i,' euler plot across contrast')
      png(paste0(figure.name,'.png'),
          width = input$across.contrast.euler.plot.width,
          height = input$across.contrast.euler.plot.height,
          res=input$across.contrast.euler.plot.res, units = 'in')
      grid.arrange(g.across.contrast()[[i]],top=textGrob('DE genes', gp=gpar(cex=1.2)))
      dev.off()
      pdf(paste0(figure.name,'.pdf'),
          width = input$across.contrast.euler.plot.width,height = input$across.contrast.euler.plot.height)
      grid.arrange(g.across.contrast()[[i]],top=textGrob('DE genes', gp=gpar(cex=1.2)))
      dev.off()
    })
    
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    
  })
  
  ###target numbers
  DDD.numbers <- reactive({
    if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans)| is.null(DDD.data$DTU_trans))
      return(NULL)
    summary3Dnumber(DE_genes = DDD.data$DE_genes,
                    DAS_genes = DDD.data$DAS_genes,
                    DE_trans = DDD.data$DE_trans,
                    DTU_trans=DDD.data$DTU_trans,contrast = DDD.data$contrast)
  })
  
  output$DDD.numbers <- renderTable({
    DDD.numbers()
  })
  
  ##---------------DE vs DAS-----------------
  DEvsDAS.results <- reactive({
    if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes))
      return(NULL)
    # summary.3D.vs(x = split(DDD.data$DE_genes$target,DDD.data$DE_genes$contrast),
    #                y = split(DDD.data$DAS_genes$target,DDD.data$DAS_genes$contrast),
    #                contrast = DDD.data$contrast,
    #                idx = c('DEonly','DE&DAS','DASonly'))
    DEvsDAS(DE_genes = DDD.data$DE_genes,
            DAS_genes = DDD.data$DAS_genes,
            contrast = DDD.data$contrast)
    
  })
  
  output$DE.vs.DAS <- renderTable({
    if(is.null(DEvsDAS.results()))
      return(NULL)
    DEvsDAS.results()
  },align='c')
  
  ##euler plot
  ###---DE vs DAS
  g.across.target1 <- reactive({
    if(is.null(DEvsDAS.results())| is.null(input$across.target))
      return(NULL)
    x <- unlist(DEvsDAS.results()[DEvsDAS.results()$Contrast==input$across.target,-1])
    if(length(x)==0)
      return(NULL)
    names(x) <- c('DE','DE&DAS','DAS')
    g <- plotEulerDiagram(x = x,fill = gg.color.hue(2),scale.plot = input$across.target.euler.plot.scale)
    g
  })
  
  output$DEvsDAS.euler <- renderPlot({
    if(is.null(g.across.target1()))
      return(NULL)
    grid.arrange(g.across.target1(),top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
  })
  
  
  ##---------------DE vs DTU-----------------
  DEvsDTU.results <- reactive({
    if(is.null(DDD.data$DE_trans) | is.null(DDD.data$DTU_trans))
      return(NULL)
    # summary.3D.vs(x = split(DDD.data$DE_trans$target,DDD.data$DE_trans$contrast),
    #                y = split(DDD.data$DTU_trans$target,DDD.data$DTU_trans$contrast),
    #                contrast = DDD.data$contrast,
    #                idx = c('DEonly','DE&DTU','DTUonly'))
    DEvsDTU(DE_trans = DDD.data$DE_trans,
            DTU_trans = DDD.data$DTU_trans,
            contrast = DDD.data$contrast)
  })
  
  output$DE.vs.DTU <- renderTable({
    if(is.null(DEvsDTU.results()))
      return(NULL)
    DEvsDTU.results()
  },align='c')
  
  ###---DE vs DTU
  g.across.target2 <- reactive({
    if(is.null(DEvsDTU.results())| is.null(input$across.target))
      return(NULL)
    x <- unlist(DEvsDTU.results()[DEvsDTU.results()$Contrast==input$across.target,-1])
    if(length(x)==0)
      return(NULL)
    names(x) <- c('DE','DE&DTU','DTU')
    g <- plotEulerDiagram(x = x,fill = gg.color.hue(2),scale.plot = input$across.target.euler.plot.scale)
    g
  })
  
  output$DEvsDTU.euler <- renderPlot({
    if(is.null(g.across.target2()))
      return(NULL)
    grid.arrange(g.across.target2(),top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
  })
  
  observeEvent(input$save.across.target.euler.plot,{
    ##DE vs DAS genes
    figure.name <- paste0(DDD.data$figure.folder,'/DE vs DAS gene euler plot in contrast ',input$across.target)
    png(paste0(figure.name,'.png'),
        width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height,
        res=input$across.target.euler.plot.res, units = 'in')
    grid.arrange(g.across.target1(),top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
    dev.off()
    pdf(paste0(figure.name,'.pdf'),
        width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height)
    grid.arrange(g.across.target1(),top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
    dev.off()
    
    ##DE vs DTU transcript
    figure.name <- paste0(DDD.data$figure.folder,'/DE vs DTU transcript euler plot in contrast ',input$across.target)
    png(paste0(figure.name,'.png'),
        width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height,
        res=input$across.target.euler.plot.res, units = 'in')
    grid.arrange(g.across.target2(),top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
    dev.off()
    
    pdf(paste0(figure.name,'.pdf'),
        width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height)
    grid.arrange(g.across.target2(),top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
    dev.off()
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  })
  
  
  ###---save stat tables
  observeEvent(input$save.stat,{
    withProgress(message = 'Save test statistics',
                 detail = 'This may take a while...', value = 0, {
                   incProgress(0.1)
                   write.csv(DDD.data$DE_genes,
                             file=paste0(DDD.data$result.folder,'/DE genes.csv'),row.names = F)
                   incProgress(0.3)
                   write.csv(DDD.data$DAS_genes,
                             file=paste0(DDD.data$result.folder,'/DAS genes.csv'),row.names = F)
                   incProgress(0.5)
                   write.csv(DDD.data$DE_trans,
                             file=paste0(DDD.data$result.folder,'/DE transcripts.csv'),row.names = F)
                   incProgress(0.7)
                   write.csv(DDD.data$DTU_trans,
                             file=paste0(DDD.data$result.folder,'/DTU transcripts.csv'),row.names = F)
                   incProgress(0.8)
                   write.csv(DDD.numbers(),file=paste0(DDD.data$result.folder,'/DE DAS DTU numbers.csv'),
                             row.names = F)
                   incProgress(0.9)
                   write.csv(DEvsDAS.results(),file=paste0(DDD.data$result.folder,'/DE vs DAS gene number.csv'),
                             row.names = F)
                   incProgress(1)
                   write.csv(DEvsDTU.results(),file=paste0(DDD.data$result.folder,'/DE vs DTU transcript number.csv'),
                             row.names = F)
                 })
    message(paste0('All tables of test statistics and target numbers are saved in folder: ',DDD.data$result.folder))
    showNotification(paste0('All tables of test statistics are saved in folder: ',DDD.data$result.folder))
  })
  
  
  
  ##---------------->>    Advanced plot   <<-----------------
  ##----------------        heatmap        ------------------
  observe({
    shinyjs::toggleState(id = 'heatmap.targetlist.type',
                         condition = input$heatmap.select.or.upload=='Upload target list')
    shinyjs::toggleState(id = 'heatmap.target.list.input',
                         condition = input$heatmap.select.or.upload=='Upload target list')
    shinyjs::toggleState(id = 'heatmap.target.type',
                         condition = input$heatmap.select.or.upload=='Select targets')
  })
  
  heatmap.targets <- eventReactive(input$plot.heatmap,{
    if(is.null(DDD.data$genes_TPM)){
      load(paste0(DDD.data$data.folder,'/genes_TPM.RData'))
      DDD.data$genes_TPM <- genes_TPM
    }
    if(is.null(DDD.data$trans_TPM)){
      load(paste0(DDD.data$data.folder,'/trans_TPM.RData'))
      DDD.data$trans_TPM <- trans_TPM
    }
    
    if(input$heatmap.select.or.upload=='Select targets'){
      if(input$heatmap.target.type=='DE genes'){
        targets <- unique(DDD.data$DE_genes$target)
        data2heatmap <- DDD.data$genes_TPM[targets,]
        column_title <- paste0(length(targets),' DE genes')
      }
      
      if(input$heatmap.target.type=='DAS genes'){
        targets <- intersect(unique(DDD.data$DAS_genes$target),unique(DDD.data$DE_genes$target))
        data2heatmap <- DDD.data$genes_TPM[targets,]
        column_title <- paste0(length(targets),' DE&DAS genes (DAS only are filtered)')
      }
      
      if(input$heatmap.target.type=='DE transcripts'){
        targets <- unique(DDD.data$DE_trans$target)
        data2heatmap <- DDD.data$trans_TPM[targets,]
        column_title <- paste0(length(targets),' DE transcripts')
      }
      
      if(input$heatmap.target.type=='DTU transcripts'){
        targets <- intersect(unique(DDD.data$DTU_trans$target),unique(DDD.data$DE_trans$target))
        data2heatmap <- DDD.data$trans_TPM[targets,]
        column_title <- paste0(length(targets),' DE&DTU transcripts (DTU only are filtered)')
      }
    }
    
    if(input$heatmap.select.or.upload=='Upload target list'){
      inFile <- input$heatmap.target.list.input
      if (is.null(inFile))
        return(NULL)
      targets <- read.csv(inFile$datapath, header = F)
      targets <- as.vector(t(targets))
      if(input$heatmap.targetlist.type=='Gene list'){
        data2heatmap <- DDD.data$genes_TPM[targets,]
        column_title <- paste0(length(targets),' Genes')
      } else {
        data2heatmap <- DDD.data$trans_TPM[targets,]
        column_title <- paste0(length(targets),' Transcripts')
      }
    }
    list(data2heatmap=data2heatmap,column_title=column_title)
  })
  
  
  heatmap.g <- eventReactive(input$plot.heatmap,{
    withProgress(message = 'Generate gene expression...', value = 0, {
      incProgress(0.1)
      data2plot <- rowmean(x = t(heatmap.targets()$data2heatmap),
                           group = DDD.data$samples$condition,
                           reorder = F)
      data2plot <- t(scale(data2plot))
      incProgress(0.2)
      hc.dist <- dist(data2plot,method = input$dist.method)
      hc <- fastcluster::hclust(hc.dist,method = input$cluster.method)
      clusters <- cutree(hc, k = input$cluster.number)
      clusters <- reorder.clusters(clusters = clusters,dat = data2plot)
      incProgress(0.3)
      g <- Heatmap(as.matrix(data2plot), name = 'Z-scores',
                   cluster_rows = TRUE,
                   clustering_method_rows=input$cluster.method,
                   row_dend_reorder = T,
                   show_row_names = FALSE,
                   show_column_names = ifelse(ncol(data2plot)>10,F,T),
                   cluster_columns = FALSE,
                   split=clusters,
                   column_title= heatmap.targets()$column_title)
      incProgress(0.8)
    })
    list(g=g,clusters=clusters)
  })
  
  output$heatmap.plot.panel <- renderPlot({
    draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
  })
  
  observeEvent(input$save.target.in.cluster,{
    clusters <- heatmap.g()$clusters
    x <- split(names(clusters),clusters)
    x <- lapply(names(x),function(i){
      data.frame(Clusters=i,Targets=x[[i]])
    })
    x <- do.call(rbind,x)
    colnames(x) <- c('Clusters','Targets')
    write.csv(x,file=paste0(DDD.data$result.folder,'/Target in each cluster heatmap ', heatmap.targets()$column_title,'.csv'),row.names = F)
    message(paste0('Table is saved in folder: ',DDD.data$result.folder))
    showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
  })
  
  ##update heigt and width
  observe({
    updateNumericInput(session,inputId = 'heatmap.plot.width',
                       value = round(pmax(10,1*length(unique(DDD.data$samples$condition)))/2.54,1))
  })
  
  observeEvent(input$save.heatmap,{
    withProgress(message = 'Saving heatmap...', detail = 'This may take a while...',
                 value = 0.5, {
                   ##DE vs DTU transcript
                   figure.name <- paste0(DDD.data$figure.folder,'/Heatmap ', input$heatmap.target.type)
                   png(paste0(figure.name,'.png'),
                       width = input$heatmap.plot.width,height = input$heatmap.plot.height,
                       res=input$heatmap.plot.res, units = 'in')
                   draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                   dev.off()
                   
                   pdf(paste0(figure.name,'.pdf'),
                       width = input$heatmap.plot.width,height = input$heatmap.plot.height)
                   draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                   dev.off()
                 })
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    
  })
  
  ##---------------profile plot-----------------
  observe({
    if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans) |is.null(DDD.data$DTU_trans))
      return(NULL)
    DE.genes <- unique(DDD.data$DE_genes$target)
    DAS.genes <- unique(DDD.data$DAS_genes$target)
    DE.trans <- unique(DDD.data$DE_trans$target)
    DTU.trans <- unique(DDD.data$DTU_trans$target)
    
    genes.ann <- set2(DE.genes,DAS.genes)
    names(genes.ann) <- c('DEonly','DE&DAS','DASonly')
    genes.ann <- plyr::ldply(genes.ann,cbind)[,c(2,1)]
    colnames(genes.ann) <- c('target','annotation')
    DDD.data$genes.ann <- genes.ann
    
    trans.ann <- set2(DE.trans,DTU.trans)
    names(trans.ann) <- c('DEonly','DE&DTU','DTUonly')
    trans.ann <- plyr::ldply(trans.ann,cbind)[,c(2,1)]
    colnames(trans.ann) <- c('target','annotation')
    DDD.data$trans.ann <- trans.ann
  })
  
  g.profiles <- eventReactive(input$make.profile.plot,{
    if(is.null(input$gene.id) | input$gene.id=="")
      return(NULL)
    gene <- trimws(input$gene.id)
    
    mapping <- DDD.data$mapping
    rownames(mapping) <- mapping$TXNAME
    
    if(is.null(DDD.data$trans_TPM)){
      load(paste0(DDD.data$data.folder,'/trans_TPM.RData'))
      DDD.data$trans_TPM <- trans_TPM
    }
    
    if(is.null(DDD.data$target_high)){
      load(paste0(DDD.data$data.folder,'/target_high.RData'))
      DDD.data$target_high <- target_high
    }
    
    if(input$profile.data.type=='TPM'){
      data.exp <- DDD.data$trans_TPM
      reps <- DDD.data$samples$condition
    }
    
    if(input$profile.data.type=='Read counts'){
      data.exp <- DDD.data$trans_counts
      reps <- DDD.data$samples_new$condition
    }
    
    if(input$profile.filter.lowexpressed=='Yes'){
      data.exp <- data.exp[DDD.data$target_high$trans_high,]
      mapping <- mapping[DDD.data$target_high$trans_high,]
    }
    
    g.pr <- plotAbundance(data.exp = data.exp,
                          gene = gene,
                          mapping = mapping,
                          genes.ann = DDD.data$genes.ann,
                          trans.ann = DDD.data$trans.ann,
                          trans.expressed = NULL,
                          reps = reps,
                          y.lab = input$profile.data.type)
    
    g.ps <- plotPS(data.exp = data.exp,
                   gene = gene,
                   mapping = mapping,
                   genes.ann = DDD.data$genes.ann,
                   trans.ann = DDD.data$trans.ann,
                   trans.expressed = NULL,
                   reps = reps,
                   y.lab = 'PS')
    
    
    
    
    return(list(g.pr=g.pr,g.ps=g.ps,gene=gene))
  })
  
  output$profile.plot.panel <- renderPlotly({
    if(is.null(g.profiles()$g.pr))
      return(NULL)
    ggplotly(g.profiles()$g.pr)
  })
  
  output$ps.plot.panel <- renderPlotly({
    if(is.null(g.profiles()$g.ps))
      return(NULL)
    ggplotly(g.profiles()$g.ps)
  })
  
  ##update heigt and width
  observe({
    updateNumericInput(session,inputId = 'ps.plot.width',
                       value = round((16+floor((length(unique(g.profiles()$g.pr$data))-1)/15))/2.54,1))
  })
  
  observeEvent(input$save.abundance.plot,{
    n <- length(unique(g.profiles()$g.pr$data))-1
    gene <- g.profiles()$gene
    folder2save <- paste0(DDD.data$figure.folder,'/Abundance plot')
    if(!file.exists(folder2save))
      dir.create(folder2save,recursive = T)
    
    png(paste0(folder2save,'/Abundance ',gene,'.png'),
        width = input$ps.plot.width,height = input$ps.plot.height,
        res=input$ps.plot.res, units = 'in')
    print(g.profiles()$g.pr)
    dev.off()
    
    pdf(paste0(folder2save,'/Abundance ',gene,'.pdf'),
        width = input$ps.plot.width,height = input$ps.plot.height)
    print(g.profiles()$g.pr)
    dev.off()
    
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    
  })
  
  
  observeEvent(input$save.ps.plot,{
    n <- length(unique(g.profiles()$g.ps$data))
    gene <- g.profiles()$gene
    folder2save <- paste0(DDD.data$figure.folder,'/PS plot')
    if(!file.exists(folder2save))
      dir.create(folder2save,recursive = T)
    
    png(paste0(folder2save,'/PS ',gene,'.png'),
        width = input$ps.plot.width,height = input$ps.plot.height,
        res=input$ps.plot.res, units = 'in')
    print(g.profiles()$g.ps)
    dev.off()
    
    pdf(paste0(folder2save,'/PS ',gene,'.pdf'),
        width = input$ps.plot.width,height = input$ps.plot.height)
    print(g.profiles()$g.ps)
    dev.off()
    
    message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    
  })
  
  
  ##---------------multiple plots----------------
  observeEvent(input$make.multiple.plot,{
    inFile <- input$multiple.gene.input
    if (is.null(inFile))
      return(NULL)
    genes <- read.csv(inFile$datapath, header = T)
    genes <- trimws(as.vector(t(genes)))
    mapping <- DDD.data$mapping
    rownames(mapping) <- mapping$TXNAME
    
    if(input$profile.data.type=='TPM'){
      data.exp <- DDD.data$trans_TPM
      reps <- DDD.data$samples$condition
    }
    
    if(input$profile.data.type=='Read counts'){
      data.exp <- DDD.data$trans_counts
      reps <- DDD.data$samples_new$condition
    }
    
    if(input$profile.filter.lowexpressed=='Yes'){
      data.exp <- data.exp[DDD.data$target_high$trans_high,]
      mapping <- mapping[DDD.data$target_high$trans_high,]
    }
    
    ###plot abundance
    folder2Abundance <- paste0(DDD.data$figure.folder,'/Abundance plot')
    if(!file.exists(folder2Abundance))
      dir.create(folder2Abundance,recursive = T)
    
    folder2PS <- paste0(DDD.data$figure.folder,'/PS plot')
    if(!file.exists(folder2PS))
      dir.create(folder2PS,recursive = T)
    
    withProgress(message = paste0('Plotting ',length(genes),' genes'),
                 detail = 'This may take a while...', value = 0, {
                   for(gene in genes){
                     incProgress(1/length(genes))
                     if(input$multiple.plot.type %in% c('Abundance','both')){
                       g.pr <- plotAbundance(data.exp = data.exp,
                                             gene = gene,
                                             mapping = mapping,
                                             genes.ann = DDD.data$genes.ann,
                                             trans.ann = DDD.data$trans.ann,
                                             trans.expressed = NULL,
                                             reps = reps,
                                             y.lab = input$profile.data.type)
                       n <- length(unique(g.pr$data))-1
                       if(input$multi.plot.format %in% c('png','both')){
                         png(paste0(folder2Abundance,'/Abundance ',gene,'.png'),
                             width = (16+floor(n/15))/2.54,
                             height = 10/2.54,
                             units = 'in',res = 300)
                         print(g.pr)
                         dev.off()
                       }
                       
                       if(input$multi.plot.format %in% c('pdf','both')){
                         pdf(paste0(folder2Abundance,'/Abundance ',gene,'.pdf'),
                             width = (16+floor(n/15))/2.54,
                             height = 10/2.54)
                         print(g.pr)
                         dev.off()
                       }
                     }
                     
                     if(input$multiple.plot.type %in% c('PS','both')){
                       g.ps <- plotPS(data.exp = data.exp,
                                      gene = gene,
                                      mapping = mapping,
                                      genes.ann = DDD.data$genes.ann,
                                      trans.ann = DDD.data$trans.ann,
                                      trans.expressed = NULL,
                                      reps = reps,
                                      y.lab = 'PS')
                       n <- length(unique(g.ps$data))
                       if(input$multi.plot.format %in% c('png','both')){
                         png(paste0(folder2PS,'/PS ',gene,'.png'),
                             width = (16+floor(n/15))/2.54,
                             height = 10/2.54,
                             units = 'in',res = 300)
                         print(g.ps)
                         dev.off()
                       }
                       if(input$multi.plot.format %in% c('pdf','both')){
                         pdf(paste0(folder2PS,'/PS ',gene,'.pdf'),
                             width = (16+floor(n/15))/2.54,
                             height = 10/2.54)
                         print(g.ps)
                         dev.off()
                       }
                     }
                     
                   }
                 })
    
    
    
    
    if(input$multiple.plot.type %in% c('PS','both')){
      g.ps <- plotPS(data.exp = data.exp,
                     gene = gene,
                     mapping = mapping,
                     genes.ann = DDD.data$genes.ann,
                     trans.ann = DDD.data$trans.ann,
                     trans.expressed = NULL,
                     reps = reps,
                     y.lab = 'PS')
      
    }
  })
  
  ##---------------GO plot----------------
  observeEvent(input$generate.target.list,{
    write.csv(data.frame(`DE genes`=unique(DDD.data$DE_genes$target)),
              file = paste0(DDD.data$result.folder,'/DE gene list union set of all contrast group.csv'),
              row.names = F)
    write.csv(data.frame(`DAS genes`=unique(DDD.data$DAS_genes$target)),
              file = paste0(DDD.data$result.folder,'/DAS gene list union set of all contrast group.csv'),
              row.names = F)
    write.csv(data.frame(`DE transcripts`=unique(DDD.data$DE_trans$target)),
              file = paste0(DDD.data$result.folder,'/DE transcript list union set of all contrast group.csv'),
              row.names = F)
    write.csv(data.frame(`DTU transcripts`=unique(DDD.data$DTU_trans$target)),
              file = paste0(DDD.data$result.folder,'/DTU transcript list union set of all contrast group.csv'),
              row.names = F)
    message(paste0('Tables of target list save in: ',DDD.data$result.folder))
    showNotification(paste0('Tables of target list save in: ',DDD.data$result.folder))
  })
  
  go.table <- reactive({
    inFile <- input$go.annotation.input
    if (is.null(inFile))
      return(NULL)
    go.table <- suppressMessages(readr::read_csv(inFile$datapath))
    go.table
  })
  
  output$go.table <- DT::renderDataTable({
    if(is.null(go.table()))
      return(NULL)
    go.table()
  },options = list(
    scrollX=T,
    columnDefs = list(list(className = 'dt-left',
                           targets = "_all"))))
  
  
  
  output$show.go.table.column <- renderUI({
    if(is.null(go.table()))
      return(NULL)
    span(
      selectInput(inputId = 'select.go.table.column',label = 'Select column of statistics',
                  choices = colnames(go.table())[-c(1:2)])
    )
  })
  
  data2goplot <- eventReactive(input$make.go.plot,{
    col.idx<- input$select.go.table.column
    data2plot <- go.table()[,c(1,2,which(colnames(go.table())==col.idx))]
    colnames(data2plot) <- c('Category','Term','Value')
    ####
    data2plot <- by(data2plot,data2plot$Category,function(x){
      x[order(x$Value,decreasing = T),]
    })
    data2plot <- do.call(rbind,data2plot)
    idx<-sapply(data2plot$Term,function(x){
      x<-as.vector(x)
      if(nchar(x)>input$go.text.cut)
        x<-paste0(substr(x,1,input$go.text.cut),'...')
      x
    })
    data2plot$Term <- factor(idx,levels = rev(idx))
    data2plot
  })
  
  row.height <- eventReactive(input$make.go.plot,{
    if(is.null(go.table()))
      return(NULL)
    nrow(go.table())
  })
  
  observeEvent(input$make.go.plot,{
    output$show.go.plot.panel <- renderUI({
      plotlyOutput('go.plot',height = pmax(row.height()*20,350))
    })
    
  })
  
  output$go.plot <- renderPlotly({
    if(is.null(data2goplot()))
      return(NULL)
    Category.idx <- unique(data2goplot()$Category)
    color.idx <- gg.color.hue(length(Category.idx))
    g.list <- lapply(1:length(Category.idx),function(i){
      x <- subset(data2goplot(),Category==Category.idx[i])
      x$Term <- factor(x$Term,levels = rev(x$Term))
      g <- ggplot(x,aes(x=Term,y=Value))+
        geom_bar(stat='identity',aes(fill=Category))+
        scale_fill_manual(values = color.idx[i])+
        coord_flip(ylim = c(0,max(data2goplot()$Value))) +
        theme_bw()+
        theme(axis.title.y = element_blank())+
        labs(y=input$select.go.table.column,title=paste0('GO annotation ',input$go.plot.label))
      if(i>1)
        g <- g+guides(fill=guide_legend(title=NULL))
      if(i!=length(Category.idx))
        g <- g+theme(axis.title.x = element_blank())
      ggplotly(g)
    })
    heights <- sapply(1:length(Category.idx),function(i) length(g.list[[i]]$x$data[[1]]$y))
    heights <- heights/sum(heights)
    subplot(g.list,nrows=length(g.list),heights = heights, titleX = TRUE)
    
  })
  
  ## update height
  observe({
    updateNumericInput(session,inputId = 'go.plot.height',
                       value = round(pmax(10,row.height()*0.4)/2.54,1))
  })
  
  
  observeEvent(input$save.go.plot,{
    if(is.null(data2goplot()))
      return(NULL)
    col.idx<- input$select.go.table.column
    data2plot <- go.table()[,c(1,2,which(colnames(go.table())==col.idx))]
    colnames(data2plot) <- c('Category','Term','Value')
    ####
    data2plot <- by(data2plot,data2plot$Category,function(x){
      x[order(x$Value,decreasing = T),]
    })
    data2plot <- do.call(rbind,data2plot)
    idx<-sapply(data2plot$Term,function(x){
      x<-as.vector(x)
      if(nchar(x)>input$go.text.cut)
        x<-paste0(substr(x,1,input$go.text.cut),'...')
      x
    })
    data2plot$Term <- factor(idx,levels = rev(idx))
    
    g <- ggplot(data2plot,aes(x=Term,y=Value))+
      geom_bar(stat='identity',aes(fill=Category))+
      coord_flip() +
      theme_bw()+
      theme(axis.title.y = element_blank())+
      facet_grid(Category~.,scales = 'free_y',space='free_y')+
      labs(y=col.idx,title=paste0('GO annotation ',input$go.plot.label))
    png(paste0(DDD.data$figure.folder,'/',input$go.plot.label,' GO annotation plot.png'),
        width = input$go.plot.width,height = input$go.plot.height,
        res=input$go.plot.res, units = 'in')
    print(g)
    dev.off()
    pdf(paste0(DDD.data$figure.folder,'/',input$go.plot.label,' GO annotation plot.pdf'),
        width = input$go.plot.width,height = input$go.plot.height)
    print(g)
    dev.off()
    message(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
    showNotification(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
  })
  
  ###---------------->>        Generate report      <<-----------------
  data.para <- reactive({
    if(is.null(DDD.data$samples) | is.null(DDD.data$samples_new)){
      srep.merge <- 'No information'
    } else {
      srep.merge <- ifelse(nrow(DDD.data$samples)>nrow(DDD.data$samples_new),'Yes','No')
    }
    
    if(is.null(DDD.data$samples)){
      srep.n <- 'No information'
    } else {
      srep.n <- length(unique(DDD.data$samples$srep))
    }
    
    if(is.null(DDD.data$path)){
      path.idx <- getwd()
    } else {
      path.idx <- DDD.data$path
    }
    
    if(is.null(DDD.data$data.folder)){
      batch.idx <- 'No information'
    } else {
      batch.idx <- ifelse(length(list.files(DDD.data$data.folder,'*batch.RData'))>0,'Yes','No')
    }
    
    if(is.null(DDD.data$data.folder)){
      data.folder.idx <- 'No information'
    } else {
      data.folder.idx <- DDD.data$data.folder
    }
    
    if(is.null(DDD.data$figure.folder)){
      figure.folder.idx <- 'No information'
    } else {
      figure.folder.idx <- DDD.data$figure.folder
    }
    
    if(is.null(DDD.data$result.folder)){
      result.folder.idx <- 'No information'
    } else {
      result.folder.idx <- DDD.data$result.folder
    }
    
    if(is.null(DDD.data$report.folder)){
      report.folder.idx <- 'No information'
    } else {
      report.folder.idx <- DDD.data$report.folder
    }
    
    x <- rbind(
      c('Folders','Directory',path.idx),#
      c('Folders','Data',data.folder.idx),#
      c('','Results',result.folder.idx),#
      c('','Figures',figure.folder.idx),#
      c('','Reports',report.folder.idx),#
      c('Data generation','tximport method',input$tximport.method),
      c('Data pre-processing','Sequencing replicates',srep.n),
      c('','Sequencing replicate merged',srep.merge),
      c('','Low expression CPM cut-off',input$cpm.cut),
      c('','Sample number for CPM cut-off', input$sample.n.cut),
      c('','Batch effect estimation',batch.idx),#
      c('','Batch effect estimation method',input$ruvseq.method),#
      c('','Normalization method',input$norm.method),
      c('DE DAS and DTU','Pipeline',input$DE.pipeline),
      c('','AS function',ifelse(input$DE.pipeline=='limma',
                                'limma::diffSplice','edgeR::diffSpliceDGE')),
      c('','P-value adjust method',input$p.adjust.method),
      c('','Adjusted p-value cut-off',input$pval.cutoff),
      c('','Log2 fold change cut-off',input$lfc.cutoff),
      c('','delta PS cut-off',input$deltaPS.cutoff),
      c('Heatmap','Distance method',input$dist.method),
      c('','Cluster method',input$cluster.method),
      c('','Number of clusters',input$cluster.number)
    )
    
    x <- data.frame(x)
    colnames(x) <- c('Step','Description','Parameter (Double click to edit if not correct)')
    x
  })
  
  x = reactiveValues(df = NULL)
  
  observe({
    df <- data.para()
    x$df <- df
  })
  
  
  output$para.summary <- DT::renderDT({
    x$df
  },editable = T,options = list(pageLength = 50))
  
  
  proxy = DT::dataTableProxy('para.summary')
  observeEvent(input$para.summary_cell_edit, {
    info = input$para.summary_cell_edit
    i = info$row
    j = info$col
    v = info$value
    
    # problem starts here
    x$df[i, j] <- isolate(suppressWarnings(DT::coerceValue(v, x$df[i, j])))
  })
  
  observeEvent(input$save.para.summary,{
    write.csv(x$df,file=paste0(DDD.data$result.folder,'/Parameter summary.csv'),row.names = F)
    message(paste0('Table of parameters is save in: ',DDD.data$result.folder))
    showNotification(paste0('Table of parameters is save in: ',DDD.data$result.folder))
  })
  
  observeEvent(input$generate.report,{
    ###download report rmarkdown
    # download.file(url = 'https://raw.githubusercontent.com/wyguo/ThreeDRNAseq/master/vignettes/report.Rmd',
    #               destfile = 'report.Rmd')
    
    if(!file.exists(paste0(DDD.data$result.folder,'/Parameter summary.csv')))
      write.csv(x$df,file=paste0(DDD.data$result.folder,'/Parameter summary.csv'),row.names = F)
    
    if(!file.exists(paste0(DDD.data$result.folder,'/data.info.csv')))
      write.csv(DDD.data$data.info,file=paste0(DDD.data$result.folder,'/data.info.csv'),row.names = T)
    
    if(!file.exists(paste0(DDD.data$result.folder,'/DE DAS DTU numbers.csv')))
      write.csv(DDD.numbers(),file=paste0(DDD.data$result.folder,'/DE DAS DTU numbers.csv'),
                row.names = F)
    if(!file.exists(paste0(DDD.data$result.folder,'/DE vs DAS gene number.csv')))
      write.csv(DEvsDAS.results(),file=paste0(DDD.data$result.folder,'/DE vs DAS gene number.csv'),
                row.names = F)
    if(!file.exists(paste0(DDD.data$result.folder,'/DE vs DTU transcript number.csv')))
      write.csv(DEvsDTU.results(),file=paste0(DDD.data$result.folder,'/DE vs DTU transcript number.csv'))
    
    
    withProgress(message = 'Generating report',
                 detail = 'This may take a while...', value = 0, {
                   incProgress(0.3)
                   ##down load report
                   report.file <- paste0(DDD.data$path,'/report.Rmd')
                   if(!file.exists(paste0(DDD.data$path,'/report.Rmd')))
                     download.file(url = 'https://raw.githubusercontent.com/wyguo/ThreeDRNAseq/master/vignettes/report.Rmd',
                                   destfile = report.file)
                   
                   ##generate report
                   generate.report(input.Rmd = paste0(DDD.data$path,"/report.Rmd"),
                                   report.folder = DDD.data$report.folder)
                   message(paste0('Reports in html, pdf and word format are saved in : ',paste0(DDD.data$path,'/report')))
                   showNotification(paste0('Reports in html, pdf and word format are saved in : ',paste0(DDD.data$path,'/report')))
                 })
  })
}
shinyApp(ui, server)

