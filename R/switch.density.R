#' Density plot
#'
#' @param x a numeric vector to plot th density.
#' @param time.points a numeric vector of the time points of time-series, e.g. 1,2,3,...
#' @param make.plotly logical, to plot \code{\link{plotly}} format figures (TRUE)
#' or general plot (FALSE)?. See details in \code{\link{ggplotly}} in \code{\link{plotly}} R pacakge.
#' @param plot.type the plot types. Options are "density" for density bar plot and "frequency" for frequency bar plot.
#' @param show.line logical, to show density or frequency line or not?
#' @param titile the title of the plot.
#' @param ... additional parameters pass to \code{\link{plotly::ggplotly}}
#'
#' @return density plot in ggplot2 format or \code{\link{plotly}} format if \code{make.plotly=T}
#'
#' @export
#'

switch.density<-function(x,time.points,
                         make.plotly=T,
                         plot.type='density',
                         show.line=T,
                         title="Density of switch points",...){

  ##plot type
  plot.type<-match.arg(plot.type,c('density','frequency'))

  data2plot<-data.frame(x=x+0.001)

  if(plot.type=='density'){
    ##density plot
    g<-ggplot(data2plot,aes(x,fill=1))+geom_histogram(aes(y=..density..),
                                                      breaks=time.points,closed='left')+theme_bw()+theme(legend.position='none')
    if(show.line)
      g<-g+stat_density(geom='line',size=1,color='red')
  } else {
    ##frequency plot
    g<-ggplot(data2plot,aes(x,fill=1))+geom_histogram(aes(y=..count..),
                                                      breaks=time.points,closed='left')+theme_bw()+theme(legend.position='none')
    if(show.line)
      g<-g+geom_freqpoly(binwidth=1,color='red',size=1)
  }

  g<-g+coord_cartesian(xlim=c(min(time.points),max(time.points)))+labs(x='Switch time points',title=title)

  if(make.plotly)
    plotly::ggplotly(g,...) else g
}

#' Barplot of switch numbers in contrast groups
#' @param x a data frame with first column of contrast groups and second column of isoform switch numbers in the contrast groups
#' @param xlab.angle rotate angle of x-axis labels.
#' @param xlab.hjust horizontal adjust of x-axis labels.
#' @param fill colour to fill the bars
#' @export
#' @examples 
#' x <- data.frame(contrast=c('B-A','C-A'),number=c(500,1000))
#' plotPWISnumber(x)
plotPWISnumber <- function(x,angle = 0,hjust = 0.5,vjust = 0.5,fill='dodgerblue3'){
  data2plot <- x
  colnames(data2plot) <- c('contrast','number')
  data2plot$contrast <- factor(data2plot$contrast,levels = unique(data2plot$contrast))
  g <- ggplot(data2plot,aes(x=contrast,y=number))+geom_bar(stat='identity',fill= fill)+
    geom_text(aes(label=number), vjust=0)+
    theme_bw()+
    theme(legend.position = 'none',axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust))+
    labs(title='Isoform switch number in contrast groups')
  return(g)
}

#' Barplot of switch numbers across time-points
#' @param x a numeric vector of isoform switch x-axis coordinates
#' @param time.points either a numeric or character vector of time-points, where the switches occur.
#' @param xlab.angle rotate angle of x-axis labels.
#' @param xlab.hjust horizontal adjust of x-axis labels.
#' @param fill colour to fill the bars
#' @export
#' @examples 
#' set.seed(1000)
#' x <- runif(100,min = 1,max = 3)
#' plotTSISnumber(x = x,time.points = c('T2','T10','T19'))
#' plotTSISnumber(x = x,time.points = 1:3)
#' 
plotTSISnumber <- function(x,time.points,angle = 0,hjust = 0.5,vjust = 0.5,fill='dodgerblue3'){
  if(!is.numeric(time.points))
    breaks <- as.numeric(factor(time.points,levels = unique(time.points))) else breaks <- time.points
  data2plot <- data.frame(x=x+0.0001)
  g<-ggplot(data2plot,aes(x))+
    geom_histogram(aes(y=..count..),breaks=breaks,closed='left',color='white',fill = fill)+
    stat_bin(geom="text",aes(label=..count..),breaks=breaks,closed='left',vjust=0) +
    theme_bw()+
    theme(legend.position='none',axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust))+
    scale_x_continuous(breaks = breaks,labels = time.points)+
    labs(x='time-points',y='number',title='Isoform switch numbers across time-points')
  return(g)
}

