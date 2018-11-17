#' Plot flow-chart graph of DE and DAS results
#' @description The flow-chart shows the results of DE vs DAS genes or DE vs DTU transcripts.
#' For example, the expressed genes are divided into DE genes and not DE genes; the DE genes
#' are divided into DE only (only transcription regulation) and DE&DAS (both transcription
#' and AS regulation) genes; the not DE genes are divided into DAS only (only AS regulation) and
#' no expression/AS change (no regulation). The same divisions to DE and/or DTU transcripts.
#' @param expressed a vector of expressed genes/transcripts.
#' @param x a vector of DE genes/transcripts.
#' @param x a vector of DAS genes/DTU transcripts.
#' @param type a character to indicate the provided list type. Options are "genes" and "transcripts".
#' @param pval.cutoff,lfc.cutoff,deltaPS.cutoff the cut-offs used to determine the significant
#' genes/transcripts.
#' @return a flow-chart plot.
#' @export
plot.flow.chart <- function(expressed,x,y,
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



