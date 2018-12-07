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
#' @param pval.cutoff,lfc.cutoff,deltaPS.cutoff the cut-offs used to determine the significant genes/transcripts.
#' @param shadow.dist a numeric value to adjust the shadow distance to the rectangles.
#' @param width.adj a numeric value to adjust the width of rectangles.
#' @param height.adj a numeric value to adjust the height of rectangles.
#' @param text.cex a numeric value of text size.
#' @param line.lwd a numeric value of line width.
#' @param border.lwd a numeric value of rectangle border line width.
#' @param plot.title plot titile
#' @return a flow-chart plot.
#' @seealso \code{\link[graphics]{rect}}
#' @export

plotFlowChart <- function(expressed,x,y,
                          type = c('genes','transcripts'),
                          pval.cutoff = 0.01,
                          lfc.cutoff = 1,
                          deltaPS.cutoff = 0.1,
                          shadow.dist = 0.7,
                          width.adj = 12,
                          height.adj = 8,
                          text.cex = 1,
                          line.lwd = 2,
                          border.lwd = 1.5,
                          plot.title = NULL
){
  if(is.null(plot.title))
    plot.title <- ifelse(type=='genes','DE and DAS genes','DE and DTU transcripts')
  
  type <- match.arg(type,c('genes','transcripts'))
  #############################################################
  ##----> generate numbers
  ns <- list()
  ns$n.expressed <- length(expressed)
  ns$n.x <- length(x)
  ns$n.notx <- ns$n.expressed-ns$n.x
  ns$n.onlyx <- length(setdiff(x,y))
  ns$n.xy <- length(intersect(x,y))
  ns$n.onlyy <- length(setdiff(y,x))
  ns$n.no<- ns$n.notx-ns$n.onlyy
  ns$n.y <- length(y)
  
  ns <- lapply(ns,function(i) prettyNum(i,big.mark = ','))
  
  ##----> plot lines
  par(mar=c(0,0,2,0))
  plot(c(0,120),c(0,80),type="n",yaxt="n",xaxt="n",xlab = NA,ylab = NA,axes = 0,main=plot.title)
  xcenter <- c(15,45,45,75,75,75,75,105)
  ycenter <- c(40,60,20,70,50,30,10,40)
  x0 <- c(15,30,30,60,60,60,60,90)
  x1 <- c(30,60,60,75,90,90,75,105)
  y1 <- y0 <- ycenter
  segments(x0 = x0,y0 = y0,x1 = x1,y1 = y1,lwd = line.lwd)
  
  x1 <- x0 <- c(30,60,60,90)
  y0 <- c(20,10,50,30)
  y1 <- c(60,30,70,50)
  segments(x0 = x0,y0 = y0,x1 = x1,y1 = y1,lwd = line.lwd)
  
  x0 <- c(75,75,90,90)
  x1 <- c(90,90,105,90)
  y0 <- c(50,30,40,30)
  y1 <- c(50,30,40,50)
  segments(x0 = x0,y0 = y0,x1 = x1,y1 = y1,col = 'red',lwd=line.lwd+1.5)
  
  ##----> add box
  xleft <- xcenter-width.adj
  xright <- xcenter+width.adj
  ybottom <- ycenter-height.adj
  ytop <- ycenter+height.adj
  rect.col <- rep('white',8)
  rect.col[c(2,4)] <- 'gold'
  rect.col[5] <- 'seagreen3'
  rect.col[c(6,8)] <- 'deepskyblue3'
  rect(xleft = xleft+shadow.dist,ybottom = ybottom-shadow.dist,xright = xright+shadow.dist,ytop = ytop-shadow.dist,
       col = 'gray',border = "transparent")
  rect(xleft = xleft,ybottom = ybottom,xright = xright,ytop = ytop,col = rect.col,lwd = border.lwd)
  
  
  ##----> add text
  text.labels <- c(
    paste0('Expressed\n',ns$n.expressed),
    paste0('DE ',type,'\n',ns$n.x),
    paste0('Not\n DE ',type,'\n',ns$n.notx),
    paste0('Only transcription\nregulation\n',ns$n.onlyx),
    paste0('Transcription\nand AS regulation\n',ns$n.xy),
    paste0('Only AS\nregulation\n',ns$n.onlyy),
    paste0('No\nregulation\n',ns$n.no),
    ifelse(type=='genes',paste0('DAS genes\n',ns$n.y),paste0('DTU transcripts\n',ns$n.y))
  )
  text(x = xcenter,y = ycenter,labels = text.labels,cex = text.cex)
  
  ##----> add statistics
  text(paste0('P-value cut-off: ',pval.cutoff,
              '\nLog2 fold change cut-off: ',lfc.cutoff,
              '\ndeltaPS cut-off: ', deltaPS.cutoff),x=0, y=10,pos=4)
  
}









#' #' Plot flow-chart graph of DE and DAS results
#' #' @description The flow-chart shows the results of DE vs DAS genes or DE vs DTU transcripts.
#' #' For example, the expressed genes are divided into DE genes and not DE genes; the DE genes
#' #' are divided into DE only (only transcription regulation) and DE&DAS (both transcription
#' #' and AS regulation) genes; the not DE genes are divided into DAS only (only AS regulation) and
#' #' no expression/AS change (no regulation). The same divisions to DE and/or DTU transcripts.
#' #' @param expressed a vector of expressed genes/transcripts.
#' #' @param x a vector of DE genes/transcripts.
#' #' @param y a vector of DAS genes/DTU transcripts.
#' #' @param type a character to indicate the provided list type. Options are "genes" and "transcripts".
#' #' @param pval.cutoff,lfc.cutoff,deltaPS.cutoff the cut-offs used to determine the significant
#' #' genes/transcripts.
#' #' @return a flow-chart plot.
#' #' @export
#' plotFlowChart <- function(expressed,x,y,
#'                             type=c('genes','transcripts'),
#'                             pval.cutoff=0.01,
#'                             lfc.cutoff=1,
#'                             deltaPS.cutoff=0.1){
#'   type <- match.arg(type,c('genes','transcripts'))
#'   n.expressed <- length(expressed)
#'   n.x <- length(x)
#'   n.notx <- n.expressed-n.x
#'   n.onlyx <- length(setdiff(x,y))
#'   n.xy <- length(intersect(x,y))
#'   n.onlyy <- length(setdiff(y,x))
#'   n.no<- n.notx-n.onlyy
#'   n.y <- length(y)
#' 
#'   grid.newpage()
#'   width <- 0.2
#'   ##---Expressed
#'   exp.bx <- boxGrob(
#'     label = paste0('Expressed\n',length(expressed)),
#'     x = 0.1,
#'     y = 0.5,
#'     width = 0.14,
#'     box_gp=gpar(fill = 'white'),
#'     just = 'center')
#' 
#'   ##---DE genes
#'   de.genes.bx <- boxGrob(
#'     label = paste0('DE ',type,'\n',n.x),
#'     x = 0.35,
#'     y = 0.7,
#'     width = 0.16,
#'     box_gp=gpar(fill = 'gold'),
#'     just = 'center')
#' 
#'   ##---not DE genes
#'   notde.genes.bx <- boxGrob(
#'     label = paste0('Not\n DE ',type,'\n',n.notx),
#'     x = 0.35,
#'     y = 0.3,
#'     width = 0.16,
#'     box_gp=gpar(fill = 'white'),
#'     just = 'center')
#' 
#'   ##---Only transcription regulation
#'   de.only.bx <- boxGrob(
#'     label = paste0('Only transcription\nregulation\n',n.onlyx),
#'     x = 0.65,
#'     y = 0.8,
#'     width = width,
#'     box_gp=gpar(fill = 'gold'),
#'     just = 'center')
#' 
#'   ##---Transcription and AS regulation
#'   de.das.bx <- boxGrob(
#'     label = paste0('Transcription\nand AS regulation\n',n.xy),
#'     x = 0.65,
#'     y = 0.6,
#'     width = width,
#'     box_gp=gpar(fill = 'seagreen3'),
#'     just = 'center')
#' 
#'   ##---Only AS regulation
#'   das.only.bx <- boxGrob(
#'     label = paste0('Only AS\nregulation\n',n.onlyy),
#'     x = 0.65,
#'     y = 0.4,
#'     width = width,
#'     box_gp=gpar(fill = 'deepskyblue3'),
#'     just = 'center')
#' 
#'   ##---No regulation
#'   no.regulation.bx <- boxGrob(
#'     label = paste0('No\nregulation\n',n.no),
#'     x = 0.65,
#'     y = 0.2,
#'     width = width,
#'     box_gp=gpar(fill = 'white'),
#'     just = 'center')
#' 
#'   ##---DAS genes
#'   das.bx <- boxGrob(
#'     label = ifelse(type=='genes',paste0('DAS genes\n',n.y),paste0('DTU transcripts\n',n.y)),
#'     x = 0.9,
#'     y = 0.5,
#'     width = 0.15,
#'     box_gp=gpar(fill = 'deepskyblue3'),
#'     just = 'center')
#' 
#' 
#'   arrow_obj <- getOption('connectGrobArrow',
#'                          default = arrow(angle = 0, ends = 'last', type = 'closed'))
#'   lty_gp <- getOption('connectGrob', default = gpar(fill = 'black',lwd=4))
#'   plot(connectGrob(exp.bx, de.genes.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#'   plot(connectGrob(exp.bx, notde.genes.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#'   plot(connectGrob(de.genes.bx, de.only.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#'   plot(connectGrob(de.genes.bx, de.das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#'   plot(connectGrob(notde.genes.bx, das.only.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#'   plot(connectGrob(notde.genes.bx, no.regulation.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#' 
#'   lty_gp <- getOption('connectGrob', default = gpar(col = 'red',lwd=4))
#'   plot(connectGrob(de.das.bx, das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#'   plot(connectGrob(das.only.bx, das.bx, 'Z', 'l',arrow_obj = arrow_obj,lty_gp = lty_gp))
#' 
#'   plot(exp.bx)
#'   plot(de.genes.bx)
#'   plot(notde.genes.bx)
#'   plot(de.only.bx)
#'   plot(de.das.bx)
#'   plot(das.only.bx)
#'   plot(no.regulation.bx)
#'   plot(das.bx)
#'   grid.text(ifelse(type=='genes','DE and DAS genes','DE and DTU transcripts'), x=0.5, y=0.95,
#'             gp=gpar(fontsize=16, col='black'))
#'   grid.text(paste0('P-value cut-off: ',pval.cutoff,
#'                    '\nLog2 fold change cut-off: ',lfc.cutoff,
#'                    '\ndeltaPS cut-off: ', deltaPS.cutoff),
#'             x=0.05, y=0.1, just = 'left',
#'             gp=gpar(fontsize=12, col='black'))
#' }
#' 
#' 
#' 
