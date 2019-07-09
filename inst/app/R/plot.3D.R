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
plotUpdown <- function(data2plot,contrast,plot.title=NULL,angle=0,
                       hjust=0.5,vjust=0.5,
                       col.up=NULL,
                       col.down=NULL){
  data2plot$contrast <- factor(data2plot$contrast,levels = unique(data2plot$contrast))
  data2plot$regulation <- factor(data2plot$regulation,levels = c('up-regulated','down-regulated'))
  g <- ggplot(data2plot,aes(x=contrast,y=number,group=regulation,
                            fill=regulation,label=number))+
    geom_bar(stat='identity')+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme_bw()+
    labs(title=plot.title)
  g <- g+theme(axis.text.x = element_text(angle = angle, hjust = hjust,vjust = vjust))
  if(!is.null(col.up) & !is.null(col.down))
    g <- g+scale_fill_manual(values=c(col.up,col.down))
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
#' @param trans.highlight a vector of transcripts. If not \code{NULL}, all the other transcripts which are not in this vector will be hidden in the plot.
#' @param reps a vecotr of replicate labels, which provide the grouping information of biolocial replciates to 
#' calculate the average expression in each conditions.
#' @param groups a vecotr of grouping labels for \code{reps}, which provide the grouping information to slice the profile plot into different segments, e.g. samples in "treatment"
#' vs "Control", samples in different development stages, etc. The vector length of \code{groups} is the same to the length of \code{reps}.
#' @param sliceProfile logical, whether to slice profile plot according to \code{groups} labels.
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
  trans.highlight=NULL,
  reps,
  groups=NULL,
  sliceProfile=F,
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
  if(!is.null(groups)){
    group.idx <- data.frame(reps=reps,groups=groups)
    idx2check <- as.vector(unlist(by(data = group.idx,INDICES = group.idx$reps,FUN = function(x){
      length(unique(x$groups))
    })))
    if(any(idx2check>1)){
      groups <- NULL
    } else {
      group.idx <- group.idx[!duplicated(group.idx$reps),]
      if(length(unique(group.idx$groups))==1){
        groups <- NULL
      } else {
        rownames(group.idx) <- group.idx$reps
        data2plot$Group <- factor(group.idx[data2plot$Conditions,'groups'])
      }
    }
  }
  
  if(!is.null(trans.highlight)){
    if(plot.gene){
      data2plot <- droplevels(data2plot[data2plot$Targets %in% c(genes.idx,trans.highlight),])
    } else {
      data2plot <- droplevels(data2plot[data2plot$Targets %in% trans.highlight,])
    }
  }
  
  legend.ncol <- ceiling(length(unique(data2plot$Targets))/15)
  
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
  if(sliceProfile & !is.null(groups))
    profiles <- profiles+facet_grid(.~Group,scales = 'free_x')
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
  groups=NULL,
  sliceProfile=F,
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
  
  if(!is.null(groups)){
    group.idx <- data.frame(reps=reps,groups=groups)
    idx2check <- as.vector(unlist(by(data = group.idx,INDICES = group.idx$reps,FUN = function(x){
      length(unique(x$groups))
    })))
    if(any(idx2check>1)){
      groups <- NULL
    } else {
      group.idx <- group.idx[!duplicated(group.idx$reps),]
      if(length(unique(group.idx$groups))==1){
        groups <- NULL
      } else {
        rownames(group.idx) <- group.idx$reps
        data2plot$Group <- factor(group.idx[data2plot$Conditions,'groups'])
      }
    }
  }
  profiles <- ggplot(data2plot,aes(x=Conditions,y=PS,group=Targets,color=Targets))+
    geom_bar(stat='identity',aes(fill=Targets))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))
  
  if(sliceProfile & !is.null(groups))
    profiles <- profiles+facet_grid(.~Group,scales = 'free_x')
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
plotGO <- function(go.table,col.idx,plot.title='GO annotation',go.text.cut=50,legend.position='right'){
  class(go.table) <- c("GO", class(go.table))
  data2plot <- go.table[,c('Category','Term',colnames(go.table)[which(colnames(go.table)==col.idx)])]
  colnames(data2plot) <- c('Category','Term','Value')
  data2plot <- by(data2plot,data2plot$Category,function(x){
    x[order(x$Value,decreasing = T),]
  })
  data2plot <- do.call(rbind,data2plot)
  idx<-sapply(data2plot$Term,function(x){
    x<-as.vector(x)
    if(nchar(x)>go.text.cut)
      x<-paste0(substr(x,1,go.text.cut),'...')
    x
  })
  data2plot$Term <- factor(idx,levels = rev(idx))
  g <- ggplot(data2plot,aes(x=Term,y=Value))+
    geom_bar(stat='identity',aes(fill=Category))+
    coord_flip() +
    theme_bw()+
    theme(axis.title.y = element_blank(),legend.position=legend.position)+
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
boxplotNormalised <- function(data.before,data.after,condition,sample.name){
  data2plot <- data.before
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g1 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5))+
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
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5))+
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
                             fill = NULL,
                             shape = c("ellipse", "circle"),scale.plot=1,
                             ...){
  if(any(is.null(fill)))
    fill <- gg.color.hue(length(x))
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

#' Volcano plot
#' @param data2plot a dataframe with a column of "target" names, a column of x coordinates, a column of y coordinates and
#' a column of "Significant" or "Not significant" to distinguish significance.
#' @param xlab the x-axis label.
#' @param ylab the y-axis label.
#' @param title the title of the plot.
#' @param col0 colour of not significant points.
#' @param col1 colour of significant points.
#' @param size the point size on the plot.
#' @param top.n show target labels of top.n distance to (0,0) coordinate.
#' @return a Volcano diagram
#' @export

#' 
plotVolcano <- function(data2plot,
                        xlab,
                        ylab,
                        title='Volcano plot',
                        col0='black',col1='red',
                        size=1,
                        top.n=10){
  data2plot$distance <- sqrt((data2plot$x)^2+(data2plot$y)^2)
  if('Significant' %in% data2plot$significance){
    data2label <- data2plot[data2plot$significance=='Significant',]
    data2label <- droplevels(data2label[order(data2label$distance,decreasing = T),][1:top.n,])
    g <- ggplot(data = data2plot,aes(x=x,y=y))+
      geom_point(aes(colour=significance),size=size)+
      scale_color_manual(values=c('Not significant'=col0,'Significant'=col1))+
      theme_bw()+
      labs(x=xlab,y=ylab,title = title)+
      scale_x_continuous(breaks = pretty(data2plot$x, n = 10))
    q <- g+
      geom_label_repel(data=data2label,
                       aes(x=x,y=y,label=target,fill=contrast),
                       alpha=0.6,
                       nudge_y      = -0.35,
                       direction    = "both",
                       segment.color = "black",
                       hjust        = 1,
                       segment.size = 0.2)
  } else {
    g <- ggplot(data = data2plot,aes(x=x,y=y))+
      geom_point(aes(colour=significance),size=size)+
      scale_color_manual(values=c('Not significant'=col0))+
      theme_bw()+
      labs(x=xlab,y=ylab,title = title)+
      scale_x_continuous(breaks = pretty(data2plot$x, n = 10))
    q <- g
  }
  q
}
