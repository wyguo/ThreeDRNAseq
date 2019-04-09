##----------Step 1: Filter low expression------------
##update parameters
observe({
  if(is.null(DDD.data$params_list$cpm_cut))
    return(NULL)
  updateSliderInput(session = session,inputId = 'cpm_cut',value = DDD.data$params_list$cpm_cut)
})

observe({
  if(is.null(DDD.data$params_list$cpm_samples_n) | is.null(DDD.data$samples_new))
    return(NULL)
  updateSliderInput(session = session,inputId = 'cpm_samples_n',
                    value = DDD.data$params_list$cpm_samples_n,min = 0,max = nrow(DDD.data$samples_new))
})

observeEvent(input$run_filter,{
  if(is.null(DDD.data$trans_counts) | is.null(DDD.data$mapping))
    return(NULL)
  cat('\nFilter low expression\n')
  target_high <- low.expression.filter(abundance = DDD.data$trans_counts,
                                       mapping = DDD.data$mapping,
                                       abundance.cut = input$cpm_cut,
                                       sample.n = input$cpm_samples_n,
                                       unit = 'counts',
                                       Log=F)
  # save(target_high,file=paste0(DDD.data$data.folder,'/target_high.RData'))
  # message(paste0('target_high.RData is saved in folder: ',DDD.data$data.folder))
  message('Done!!!')
  showNotification("Done!!!")
  DDD.data$target_high <- target_high
})

observe({
  x <- data.frame(
    Description=c('Raw transcripts',
                  'Raw genes',
                  'Samples',
                  'Samples after merging seq-reps',
                  'Condition of interest',
                  'CPM cut-off',
                  'Min samples to CPM cut-off',
                  'Expressed transcripts',
                  'Expressed genes'),
    Number=c(ifelse(is.null(DDD.data$mapping),'--',length(DDD.data$mapping$TXNAME)),
             ifelse(is.null(DDD.data$mapping),'--',length(unique(DDD.data$mapping$GENEID))),
             ifelse(is.null(DDD.data$samples),'--',nrow(DDD.data$samples)),
             ifelse(is.null(DDD.data$samples_new),'--',nrow(DDD.data$samples_new)),
             ifelse(is.null(DDD.data$samples_new),'--',length(unique(DDD.data$samples[,'condition']))),
             ifelse(is.null(DDD.data$params_list$cpm_cut),input$cpm_cut,DDD.data$params_list$cpm_cut),
             ifelse(is.null(DDD.data$params_list$cpm_samples_n),input$cpm_samples_n,DDD.data$params_list$cpm_samples_n),
             ifelse(is.null(DDD.data$target_high),'--',length(DDD.data$target_high$trans_high)),
             ifelse(is.null(DDD.data$target_high),'--',length(DDD.data$target_high$genes_high))
             )
  )
  DDD.data$RNAseq_info <- x
})

output$RNAseq_datainfo <- DT::renderDataTable({
  if(is.null(DDD.data$RNAseq_info))
    return(NULL)
  datatable(DDD.data$RNAseq_info,options = list(pageLength = 20))
})

##-------------- > mean-variance plots---------------
###--mean variance trans trend plot
mv.trans <- eventReactive(input$run_filter,{
  if(is.null(DDD.data$samples_new) | is.null(DDD.data$trans_counts) | is.null(DDD.data$target_high))
    return(NULL)
  condition <- paste0(DDD.data$samples_new$condition)
  ###---transcript levels
  counts.raw = DDD.data$trans_counts[rowSums(DDD.data$trans_counts>0)>0,]
  counts.filtered = DDD.data$trans_counts[DDD.data$target_high$trans_high,]
  # graphics.off()
  check.mean.variance(counts.raw = counts.raw,
                      counts.filtered = counts.filtered,
                      condition = condition)
})

mv.trans.plot <- eventReactive(input$run_filter,{
  if(is.null(mv.trans()))
    return(NULL)
  mv.trans.plot <- function(){
    fit.raw <- mv.trans()$fit.raw
    fit.filtered <- mv.trans()$fit.filtered
    par(mfrow=c(1,2))
    plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
    title('\n\nRaw counts')
    plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     l = fit.filtered$l,lwd=2,col='black')
    title(paste0('\n\nFiltered counts (cutoff: cpm=',input$cpm_cut,'; sample=',input$cpm_samples_n,')'))
    lines(fit.raw$l, col = "gold",lty=4,lwd=2)
    legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
           legend = c('low-exp removed','low-exp kept'))
  }
  cat('\nMake mean-variance trend plot\n')
  return(mv.trans.plot)
})
#
output$mv.trans.plot <- renderPlot({
  if(is.null(mv.trans.plot()))
    return(NULL)
  # cat('\nMake mean-variance trend plot\n')
  mv.trans.plot()()
})

###--mean variance genes trend plot
mv.genes <- eventReactive(input$run_filter,{
  if(is.null(DDD.data$samples_new) | is.null(DDD.data$genes_counts) | is.null(DDD.data$target_high))
    return(NULL)
  condition <- paste0(DDD.data$samples_new$condition)
  ###---genes levels
  counts.raw = DDD.data$genes_counts[rowSums(DDD.data$genes_counts>0)>0,]
  counts.filtered = DDD.data$genes_counts[DDD.data$target_high$genes_high,]
  # graphics.off()
  check.mean.variance(counts.raw = counts.raw,
                      counts.filtered = counts.filtered,
                      condition = condition)
})

mv.genes.plot <- eventReactive(input$run_filter,{
  if(is.null(mv.genes()))
    return(NULL)
  mv.genes.plot <- function(){
    fit.raw <- mv.genes()$fit.raw
    fit.filtered <- mv.genes()$fit.filtered
    par(mfrow=c(1,2))
    plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
    title('\n\nRaw counts')
    plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     l = fit.filtered$l,lwd=2,col='black')
    title(paste0('\n\nFiltered counts (cutoff: cpm=',input$cpm_cut,'; sample=',input$cpm_samples_n,')'))
    lines(fit.raw$l, col = "gold",lty=4,lwd=2)
    legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
           legend = c('low-exp removed','low-exp kept'))
  }
  return(mv.genes.plot)
})
#
output$mv.genes.plot <- renderPlot({
  if(is.null(mv.genes.plot()))
    return(NULL)
  mv.genes.plot()()
})

observeEvent(input$save.mv.plot,{
  if(is.null(mv.trans.plot()) | is.null(mv.genes.plot()))
    return(NULL)
 
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


##--------------Step 2: PCA plot---------------
###update the pca selections in the selectinput
observe({
  if(is.null(input$sample.input)){
    x <- 5
  } else {
    if(input$pca.plot.type=='all samples'){
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
observe({
  if(is.null(DDD.data$samples_new))
    return(NULL)
  updateSelectInput(session = session,inputId = 'pca_color_option',
                    label = 'Colour samples based on',
                    choices = colnames(DDD.data$samples_new))
})

###---trans level
pca.trans.g <- eventReactive(input$plot.pca.button,{
  # if(is.null(input$pca_color_option))
  #   return(warning('Please select a factor/factors to colo'))
  data2pca <- DDD.data$trans_counts[DDD.data$target_high$trans_high,]
  dge <- DGEList(counts=data2pca)
  dge <- calcNormFactors(dge)
  data2pca <- t(counts2CPM(obj = dge,Log = T))
  dim1 <- input$dim1
  dim2 <- input$dim2
  ellipse.type <- input$pca.add.circle
  
  ##--Bio-reps
  if(input$pca.plot.type=='all samples'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    groups <- as.vector(interaction(DDD.data$samples_new[,input$pca_color_option]))
    g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                    groups = groups,plot.title = 'Transcript PCA: all samples',
                    ellipse.type = ellipse.type,
                    add.label = T,adj.label = F)
  }
  
  ##--average expression plot
  if(input$pca.plot.type=='average expression'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
    groups <- unique(DDD.data$samples_new$condition)
    g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                    groups = groups,plot.title = 'Transcript PCA: average expression',
                    ellipse.type = 'none',add.label = T,adj.label = F)
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
  if(input$pca.plot.type=='all samples'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    # if(input$pca_color_option=='Bio-reps')
    #   groups <- DDD.data$samples_new$brep
    # if(input$pca_color_option=='Conditions')
    #   groups <- DDD.data$samples_new$condition
    groups <- as.vector(interaction(DDD.data$samples_new[,input$pca_color_option]))
    g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                    groups = groups,plot.title = 'Gene PCA: all samples',
                    ellipse.type = ellipse.type,
                    add.label = T,adj.label = F)
  }
  
  ##--average expression plot
  if(input$pca.plot.type=='average expression'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
    groups <- unique(DDD.data$samples_new$condition)
    g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                    groups = groups,plot.title = 'Gene PCA: average expression',
                    ellipse.type = 'none',add.label = T,adj.label = F)
  }
  g
})

output$genes.pca.plot <- renderPlotly({
  ggplotly(pca.genes.g())
})

observeEvent(input$save.pca.plot,{
  ##transcript level
  # graphics.off()
  if(input$pca.plot.type=='all samples'){
    figure.label <- paste0('all samples ',paste0(input$pca_color_option,collapse = '_'))
  } else {
    figure.label <- paste0('average expression')
  }
  
  png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  suppressWarnings(print(pca.trans.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper))))
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  suppressWarnings(print(pca.trans.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper))))
  dev.off()
  ###gene level
  # graphics.off()
  png(filename = paste0(DDD.data$figure.folder,'/Gene PCA ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  suppressWarnings(print(pca.genes.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper))))
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  suppressWarnings(print(pca.genes.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper))))
  dev.off()
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})



##--------------Step 3: Remove batch effects-trans level---------------
## update does the data has batch effects
observe({
  if(is.null(DDD.data$params_list$has_batcheffect))
    return(NULL)
  updateRadioButtons(session = session,inputId = 'has_batcheffect',selected = DDD.data$params_list$has_batcheffect)
})

observe({
  if(is.null(DDD.data$params_list$RUVseq_method))
    return(NULL)
  updateRadioButtons(session = session,inputId = 'RUVseq_method',selected = DDD.data$params_list$RUVseq_method)
})

###---update remove_batch action button
observe({
  if(input$has_batcheffect=='Yes'){
    updateButton(session, "remove_batch", disabled = F,icon = icon('cut'),
                 style="primary")
    updateButton(session, "update_pca_plot", disabled = F,icon = icon('refresh'),
                 style="primary")
  } else {
    updateButton(session, "remove_batch", disabled = T,icon = icon('ban'))
    updateButton(session, "update_pca_plot", disabled = T,icon = icon('ban'))
  }
})

###---transcript level
observeEvent(input$remove_batch,{
  cat('\nEstimate transcript level batch effects...\n')
  withProgress(message = 'Estimate transcript level batch effects...', value = 0, {
    trans_batch <- remove.batch.shiny(read.counts = DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                                      condition = DDD.data$samples_new$condition,
                                      design = NULL,
                                      contrast=NULL,
                                      group = DDD.data$samples_new$condition,
                                      method = input$RUVseq_method,
                                      k = 1)
    incProgress(0.7)
    # save(trans_batch,file=paste0(DDD.data$data.folder,'/trans_batch.RData'))
    # message(paste0('trans_batch.RData is saved in folder: ',DDD.data$data.folder))
    DDD.data$trans_batch <- trans_batch
    incProgress(1)
  })
})

###---gene level
observeEvent(input$remove_batch,{
  cat('\nEstimate gene level batch effects...\n')
  withProgress(message = 'Estimate gene level batch effects...', value = 0, {
    genes_batch <- remove.batch.shiny(read.counts = DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                                      condition = DDD.data$samples_new$condition,
                                      design = NULL,
                                      contrast=NULL,
                                      group = DDD.data$samples_new$condition,
                                      method = input$RUVseq_method,
                                      k = 1)
    incProgress(0.7)
    # save(genes_batch,file=paste0(DDD.data$data.folder,'/genes_batch.RData'))
    # message(paste0('genes_batch.RData is saved in folder: ',DDD.data$data.folder))
    DDD.data$genes_batch <- genes_batch
    incProgress(1)
  })
})

# ###---trans level
pca.trans.br.g <- eventReactive(input$update_pca_plot,{
  message('Transcript level PCA--batch removed')
  data2pca <- DDD.data$trans_batch$normalizedCounts
  dge <- DGEList(counts=data2pca)
  dge <- suppressWarnings(calcNormFactors(dge))
  data2pca <- t(counts2CPM(obj = dge,Log = T))
  dim1 <- input$dim1
  dim2 <- input$dim2
  ellipse.type <- input$pca.add.circle
  
  ##--Bio-reps
  if(input$pca.plot.type=='all samples'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    groups <- as.vector(interaction(DDD.data$samples_new[,input$pca_color_option]))
    g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                    groups = groups,plot.title = 'Transcript PCA: all samples (batch effect removed)',
                    ellipse.type = ellipse.type,
                    add.label = T,adj.label = F)
  }
  
  ##--average expression plot
  if(input$pca.plot.type=='average expression'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    groups <- DDD.data$samples_new$brep
    data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
    groups <- unique(DDD.data$samples_new$condition)
    g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                    groups = groups,plot.title = 'Transcript PCA: average expression (batch effect removed)',
                    ellipse.type = 'none',add.label = T,adj.label = F)
  }
  g
})
#

output$trans.pca.br.plot <- renderPlotly({
  ggplotly(pca.trans.br.g())
})


# ###---genes level
pca.genes.br.g <- eventReactive(input$update_pca_plot,{
  message('Gene level PCA--batch removed')
  data2pca <- DDD.data$genes_batch$normalizedCounts
  dge <- DGEList(counts=data2pca)
  dge <- suppressWarnings(calcNormFactors(dge))
  data2pca <- t(counts2CPM(obj = dge,Log = T))
  dim1 <- input$dim1
  dim2 <- input$dim2
  ellipse.type <- input$pca.add.circle
  
  ##--Bio-reps
  if(input$pca.plot.type=='all samples'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    groups <- as.vector(interaction(DDD.data$samples_new[,input$pca_color_option]))
    g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                    groups = groups,plot.title = 'Gene PCA: all samples (batch effect removed)',
                    ellipse.type = ellipse.type,
                    add.label = T,adj.label = F)
  }
  
  ##--average expression plot
  if(input$pca.plot.type=='average expression'){
    rownames(data2pca) <- gsub('_','.',rownames(data2pca))
    groups <- DDD.data$samples_new$brep
    data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
    groups <- unique(DDD.data$samples_new$condition)
    g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                    groups = groups,plot.title = 'Gene PCA: average expression (batch effect removed)',
                    ellipse.type = 'none',add.label = T,adj.label = F)
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
  if(input$pca.plot.type=='all samples'){
    figure.label <- paste0('all samples ',paste0(input$pca_color_option,collapse = '_'))
  } else {
    figure.label <- paste0('average expression')
  }
  
  png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  print(pca.trans.br.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  print(pca.trans.br.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  dev.off()
  
  ###gene level
  # graphics.off()
  png(filename = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  print(pca.genes.br.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  print(pca.genes.br.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  dev.off()
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})

##--------------Step 4-normalization--------------
##update normalisation method
observe({
  if(is.null(DDD.data$params_list$norm_method))
    return(NULL)
  updateRadioButtons(session = session,inputId = 'norm_method',selected = DDD.data$params_list$norm_method)
})

observeEvent(input$run_norm,{
  cat('\nNormalise transcript read counts\n')
  withProgress(message = 'Normalise transcript read counts...', value = 0, {
    dge <- DGEList(counts=DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                   group = DDD.data$samples_new$condition,
                   genes = DDD.data$mapping[DDD.data$target_high$trans_high,])
    incProgress(0.7)
    trans_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm_method))
    # save(trans_dge,file=paste0(DDD.data$data.folder,'/trans_dge.RData'))
    # message(paste0('trans_dge.RData is saved in folder: ',DDD.data$data.folder))
    message('Done!!!')
    DDD.data$trans_dge <- trans_dge
    incProgress(1)
  })
})

observeEvent(input$run_norm,{
  cat('\nNormalise gene read counts\n')
  withProgress(message = 'Normalise gene read counts...', value = 0, {
    dge <- DGEList(counts=DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                   group = DDD.data$samples_new$condition)
    incProgress(0.7)
    genes_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm_method))
    # save(genes_dge,file=paste0(DDD.data$data.folder,'/genes_dge.RData'))
    # message(paste0('genes_dge.RData is saved in folder: ',DDD.data$data.folder))
    message('Done!!!')
    DDD.data$genes_dge <- genes_dge
    incProgress(1)
  })
})

##-------------> make the distribution plot-------------
###---transcript level
trans.dist.g <- eventReactive(input$norm.plot.button,{
  # sample.name <- paste0(DDD.data$samples_new$condition,'.',DDD.data$samples_new$brep)
  sample.name <- DDD.data$samples_new$sample.name
  condition <- DDD.data$samples_new$condition
  data.before <- DDD.data$trans_counts[DDD.data$target_high$trans_high,]
  data.after <- counts2CPM(obj = DDD.data$trans_dge,Log = T)
  g <- boxplotNormalised(data.before = data.before,
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
  sample.name <- DDD.data$samples_new$sample.name
  condition <- DDD.data$samples_new$condition
  data.before <- DDD.data$genes_counts[DDD.data$target_high$genes_high,]
  data.after <- counts2CPM(obj = DDD.data$genes_dge,Log = T)
  g <- boxplotNormalised(data.before = data.before,
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
  png(filename = paste0(DDD.data$figure.folder,'/Transcript expression distribution.png'),
      width = input$dist.plot.width,height = input$dist.plot.height,res=input$dist.plot.res, units = 'in')
  gridExtra::grid.arrange(trans.dist.g()$g1,trans.dist.g()$g2,ncol=1)
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Transcript expression distribution.pdf'),
      width = input$dist.plot.width,height = input$dist.plot.height)
  gridExtra::grid.arrange(trans.dist.g()$g1,trans.dist.g()$g2,ncol=1)
  dev.off()
  
  ##gene level
  # graphics.off()
  png(filename = paste0(DDD.data$figure.folder,'/Gene expression distribution.png'),
      width = input$dist.plot.width,height = input$dist.plot.height,res=input$dist.plot.res, units = 'in')
  gridExtra::grid.arrange(genes.dist.g()$g1,genes.dist.g()$g2,ncol=1)
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Gene expression distribution.pdf'),
      width = input$dist.plot.width,height = input$dist.plot.height)
  gridExtra::grid.arrange(genes.dist.g()$g1,genes.dist.g()$g2,ncol=1)
  dev.off()
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})

observeEvent(input$page_before_preprocessing, {
  newtab <- switch(input$tabs, "preprocessing" = "generation","generation" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
})

observeEvent(input$page_after_preprocessing, {
  newtab <- switch(input$tabs, "preprocessing" = "ddd","ddd" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
})

