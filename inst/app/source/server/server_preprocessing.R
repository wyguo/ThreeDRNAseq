##----------Step 1: Filter low expression------------
observeEvent(input$tabs,{
  if(input$tabs=='preprocessing'){
    text2show <- 'Page: Data pre-processing'
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
    showmessage(text = text2show,showNoteify = F)
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
  }
})

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

observe({
  if(is.null(DDD.data$brep_column))
    return(NULL)
  updateSelectInput(session =session,inputId = 'pca_color_option',selected = DDD.data$brep_column)
})

observeEvent(input$run_filter,{
  if(is.null(DDD.data$trans_counts) | is.null(DDD.data$mapping)){
    showmessage('Please finish the "Data generation" step before data pre-processing.',
                duration = NULL,type = 'error')
    return(NULL)
  }
    

  startmessage(text = 'Filter low expression')
  target_high <- low.expression.filter(abundance = DDD.data$trans_counts,
                                       mapping = DDD.data$mapping,
                                       abundance.cut = input$cpm_cut,
                                       sample.n = input$cpm_samples_n,
                                       unit = 'counts',
                                       Log=F)
  # save(target_high,file=paste0(DDD.data$data.folder,'/target_high.RData'))
  # message(paste0('target_high.RData is saved in folder: ',DDD.data$data.folder))
  DDD.data$target_high <- target_high
  endmessage(text = 'Filter low expression')
})

observeEvent(input$run_filter,{
  DDD.data$params_list$cpm_cut <- input$cpm_cut
  DDD.data$params_list$cpm_samples_n <- input$cpm_samples_n
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
mv.trans.plot <- eventReactive(input$run_filter,{
  if(is.null(DDD.data$samples_new) | is.null(DDD.data$trans_counts) | is.null(DDD.data$target_high))
    return(NULL)
  startmessage(text = 'Estimate transcript mean-variance')
  condition <- paste0(DDD.data$samples_new$condition)
  ###---transcript levels
  counts.raw = DDD.data$trans_counts[rowSums(DDD.data$trans_counts>0)>0,]
  counts.filtered = DDD.data$trans_counts[DDD.data$target_high$trans_high,]
  # graphics.off()
  mv.trans <- check.mean.variance(counts.raw = counts.raw,
                      counts.filtered = counts.filtered,
                      condition = condition)
  endmessage(text = 'Estimate transcript mean-variance')
  
  mv.trans.plot <- function(){
    fit.raw <- mv.trans$fit.raw
    fit.filtered <- mv.trans$fit.filtered
    par(mfrow=c(1,2))
    plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     xlim=c(min(fit.raw$sx),max(fit.raw$sx)),
                     ylim=c(min(fit.raw$sy),max(fit.raw$sy)),
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
    title('\n\nRaw counts')
    plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     xlim=c(min(fit.raw$sx),max(fit.raw$sx)),
                     ylim=c(min(fit.raw$sy),max(fit.raw$sy)),
                     l = fit.filtered$l,lwd=2,col='black')
    title(paste0('\n\nFiltered counts (cutoff: cpm=',DDD.data$params_list$cpm_cut,
                 '; sample=',DDD.data$params_list$cpm_samples_n,')'))
    lines(fit.raw$l, col = "gold",lty=4,lwd=2)
    legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
           legend = c('low-exp removed','low-exp kept'))
  }
  return(mv.trans.plot)
})
#
output$mv.trans.plot <- renderPlot({
  if(is.null(mv.trans.plot()))
    return(NULL)
  startmessage(text = 'Plot transcript mean-variance trend')
  mv.trans.plot()()
  removeNotification(id = 'message_id')
  endmessage(text = 'Plot transcript mean-variance trend')
})

###--mean variance genes trend plot
mv.genes.plot <- eventReactive(input$run_filter,{
  if(is.null(DDD.data$samples_new) | is.null(DDD.data$genes_counts) | is.null(DDD.data$target_high))
    return(NULL)
  startmessage(text = 'Estimate gene mean-variance trend')
  condition <- paste0(DDD.data$samples_new$condition)
  ###---genes levels
  counts.raw = DDD.data$genes_counts[rowSums(DDD.data$genes_counts>0)>0,]
  counts.filtered = DDD.data$genes_counts[DDD.data$target_high$genes_high,]
  # graphics.off()
  mv.genes <- check.mean.variance(counts.raw = counts.raw,
                      counts.filtered = counts.filtered,
                      condition = condition)
  
  endmessage(text = 'Estimate gene mean-variance trend')
  
  mv.genes.plot <- function(){
    fit.raw <- mv.genes$fit.raw
    fit.filtered <- mv.genes$fit.filtered
    par(mfrow=c(1,2))
    plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     xlim=c(min(fit.raw$sx),max(fit.raw$sx)),
                     ylim=c(min(fit.raw$sy),max(fit.raw$sy)),
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
    title('\n\nRaw counts')
    plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     xlim=c(min(fit.raw$sx),max(fit.raw$sx)),
                     ylim=c(min(fit.raw$sy),max(fit.raw$sy)),
                     l = fit.filtered$l,lwd=2,col='black')
    title(paste0('\n\nFiltered counts (cutoff: cpm=',DDD.data$params_list$cpm_cut,
                 '; sample=',DDD.data$params_list$cpm_samples_n,')'))
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
  startmessage(text = 'Plot gene mean-variance trend')
  mv.genes.plot()()
  endmessage(text = 'Plot gene mean-variance trend')
})

observeEvent(input$save_mv_plot,{
  if(is.null(mv.trans.plot()) | is.null(mv.genes.plot()))
    return(NULL)
  startmessage(text = 'Save mean-variance plots')
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
  endmessage(text = 'Save mean-variance plots')
  showmessage(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})

##*****************preview saved plot******************
observeEvent(input$save_mv_plot_view, {
  file2save1 <- paste0(DDD.data$figure.folder,'/Transcript mean-variance trend.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/Gene mean-variance trend.png')
  
  if(!file.exists(file2save1) | !file.exists(file2save2)){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      p('Transcript level mean-variance trend plot'),
      tags$img(src=base64enc::dataURI(file = file2save1),
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
             border:1px solid #000000"),
      p('Gene level mean-variance trend plot'),
      tags$img(src=base64enc::dataURI(file = file2save2),
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
             border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  }

})

##---------------> Transcripts per gene number -----------
tpg_g <- eventReactive(input$tpg_plot_button,{
  if(is.null(DDD.data$mapping) | is.null(DDD.data$target_high))
    return(NULL)
  
  ##before
  x <- DDD.data$mapping$GENEID
  y <- table(table(x))
  n.idx <- as.numeric(names(y))
  y <- y[order(n.idx,decreasing = F)]
  z <- c(y[n.idx<=10],sum(y[n.idx > 10 & n.idx <=20]),sum(y[n.idx > 20 & n.idx <=30]),sum(y[n.idx > 30]))
  names(z) <- c(names(y[n.idx<=10]),'11-20','21-30','>30')
  data2plot <- data.frame(trans=names(z),genes=z)
  data2plot$trans <- factor(data2plot$trans,levels = data2plot$trans)
  data2plot$Group <- 'Before filter'
  data2plot.before <- data2plot
  
  ##after
  x <- DDD.data$mapping[DDD.data$target_high$trans_high,]$GENEID
  y <- table(table(x))
  n.idx <- as.numeric(names(y))
  y <- y[order(n.idx,decreasing = F)]
  z <- c(y[n.idx<=10],sum(y[n.idx > 10 & n.idx <=20]),sum(y[n.idx > 20 & n.idx <=30]),sum(y[n.idx > 30]))
  names(z) <- c(names(y[n.idx<=10]),'11-20','21-30','>30')
  data2plot <- data.frame(trans=names(z),genes=z)
  data2plot$trans <- factor(data2plot$trans,levels = data2plot$trans)
  data2plot$Group <- 'After filter'
  data2plot.after <- data2plot
  
  data2plot <- rbind(data2plot.before,data2plot.after)
  data2plot$Group <- factor(data2plot$Group,levels = c('Before filter','After filter'))
  g <- ggplot(data = data2plot,aes(x=trans,y=genes,fill=Group))+
    geom_bar(stat='identity',position = "dodge")+
    labs(x='Number of transcript per gene',y='Number of genes',
         title='Distribution of the number of transcripts per gene')+
    # geom_text(size = 3, position = position_dodge(0.9))+
    geom_text(aes(label = genes),size=2, position = position_dodge(0.9),vjust = 0)+
    scale_fill_manual(values = c(input$tpg_color_before,input$tpg_color_after))+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = input$tpg_number_x_rotate,
                                     hjust = input$tpg_number_x_hjust,
                                     vjust = input$tpg_number_x_vjust))
  return(g)
})

output$tpg_plot <- renderPlot({
  if(is.null(tpg_g()))
    return(NULL)
  startmessage(text = 'Plot number of transcript per gene')
  print(tpg_g())
  endmessage(text = 'Plot number of transcript per gene')
})

##save plot
observeEvent(input$tpg_save_button,{
  if(is.null(tpg_g()))
    return(NULL)
  startmessage(text = 'Save number of transcript per gene plot')
  create.folders(wd = DDD.data$folder)
  folder2save <- paste0(DDD.data$folder,'/figure')
  png(paste0(folder2save,'/Distribution of the number of transcripts per gene.png'),
      width = input$tpg_number_width,height = input$tpg_number_height,res = input$tpg_number_res,units = 'in')
  print(tpg_g())
  dev.off()
  
  pdf(paste0(folder2save,'/Distribution of the number of transcripts per gene.pdf'),
      width = input$tpg_number_width,height = input$tpg_number_height)
  print(tpg_g())
  dev.off()
  endmessage(text = 'Save number of transcript per gene plot')
  showmessage(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})


##*****************preview saved plot******************
observeEvent(input$tpg_save_button_view, {
  file2save <- paste0(DDD.data$figure.folder,'/Distribution of the number of transcripts per gene.png')
  
  if(!file.exists(file2save)){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      tags$img(src=base64enc::dataURI(file = file2save),
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
      ))
  }
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
                    choices = colnames(DDD.data$samples_new),
                    selected = ifelse(is.null(input$brep_column),'',input$brep_column))
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
    if(is.null(input$pca_color_option)){
      showmessage(text = 'Please select a factor to group and colour the samples on PCA plots.')
      return(NULL)
    }
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

observeEvent(input$plot.pca.button,{
  startmessage(text = 'Plot transcript level PCA')
  output$trans.pca.plot <- renderPlotly({
    if(is.null(pca.trans.g()))
      return(NULL)
    ggplotly(pca.trans.g())
    
  })
  endmessage(text = 'Plot transcript level PCA')
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
    if(is.null(input$pca_color_option)){
      showmessage(text = 'Please select a factor to group and colour the samples on PCA plots.')
      return(NULL)
    }
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

observeEvent(input$plot.pca.button,{
  startmessage(text = 'Plot gene level PCA')
  output$genes.pca.plot <- renderPlotly({
    if(is.null(pca.genes.g()))
      return(NULL)
    ggplotly(pca.genes.g())
  })
  endmessage(text = 'Plot gene level PCA')
})

observeEvent(input$save_pca_plot,{
  ##transcript level
  # graphics.off()
  startmessage(text = 'Save transcript and gene level PCA plots')
  if(input$pca.plot.type=='all samples'){
    figure.label <- paste0('all samples ',paste0(input$pca_color_option,collapse = '_'))
  } else {
    figure.label <- paste0('average expression')
  }
  
  pca.trans.g <- suppressMessages(pca.trans.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  print(pca.trans.g)
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  print(pca.trans.g)
  dev.off()
  
  ###gene level
  # graphics.off()
  pca.genes.g <- suppressMessages(pca.genes.g()+coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  png(filename = paste0(DDD.data$figure.folder,'/Gene PCA ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  print(pca.genes.g)
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  print(pca.genes.g)
  dev.off()
  endmessage(text = 'Save transcript and gene level PCA plots')
  
  showmessage(paste0('Figures are saved in folder: ',DDD.data$figure.folder))

})


##***************** preview saved plots *****************
observeEvent(input$save_pca_plot_view, {
  if(input$pca.plot.type=='all samples'){
    figure.label <- paste0('all samples ',paste0(input$pca_color_option,collapse = '_'))
  } else {
    figure.label <- paste0('average expression')
  }
  
  gene <- trimws(input$gene.id)
  file2save1 <- paste0(DDD.data$figure.folder,'/Transcript PCA ',figure.label,'.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/Gene PCA ',figure.label,'.png')
  
  if(!file.exists(file2save1) | !file.exists(file2save2)){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      tags$img(src=base64enc::dataURI(file = file2save1),
               style="width: 80%; display: block; margin-left: auto; margin-right: auto; 
               border:1px solid #000000"),
      p(),
      tags$img(src=base64enc::dataURI(file = file2save2),
               style="width: 80%; display: block; margin-left: auto; margin-right: auto; 
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
      )
    )
  }
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
    DDD.data$genes_batch <- NULL
    DDD.data$trans_batch <- NULL
    updateButton(session, "remove_batch", disabled = T,icon = icon('ban'))
    updateButton(session, "update_pca_plot", disabled = T,icon = icon('ban'))
  }
})

###---transcript level
observeEvent(input$remove_batch,{
  startmessage(text = 'Estimate transcript level batch effects')
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
  endmessage(text = 'Estimate transcript level batch effects')
})

###---gene level
observeEvent(input$remove_batch,{
  startmessage(text = 'Estimate gene level batch effects')
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
  endmessage(text = 'Estimate gene level batch effects')
})

# ###---trans level
pca.trans.br.g <- eventReactive(input$update_pca_plot,{
  if(is.null(DDD.data$trans_batch))
    return(NULL)
  showmessage('Transcript level PCA--batch removed',showNoteify = F)
  data2pca <- DDD.data$trans_batch$normalizedCounts
  dge <- DGEList(counts=data2pca)
  dge <- suppressWarnings(calcNormFactors(dge))
  data2pca <- t(counts2CPM(obj = dge,Log = T))
  dim1 <- input$dim1
  dim2 <- input$dim2
  ellipse.type <- input$pca.add.circle
  
  ##--Bio-reps
  if(input$pca.plot.type=='all samples'){
    if(is.null(input$pca_color_option)){
      showmessage(text = 'Please select a factor to group and colour the samples on PCA plots.')
      return(NULL)
    }
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

observeEvent(input$update_pca_plot,{
  startmessage(text = 'Plot transcript level PCA after batch effect removed')
  output$trans.pca.br.plot <- renderPlotly({
    if(is.null(pca.trans.br.g()))
      return(NULL)
    ggplotly(pca.trans.br.g())
  })
  endmessage(text = 'Plot transcript level PCA after batch effect removed')
})


# ###---genes level
pca.genes.br.g <- eventReactive(input$update_pca_plot,{
  if(is.null(DDD.data$trans_batch))
    return(NULL)
  showmessage('Gene level PCA--batch removed',showNoteify = F)
  data2pca <- DDD.data$genes_batch$normalizedCounts
  dge <- DGEList(counts=data2pca)
  dge <- suppressWarnings(calcNormFactors(dge))
  data2pca <- t(counts2CPM(obj = dge,Log = T))
  dim1 <- input$dim1
  dim2 <- input$dim2
  ellipse.type <- input$pca.add.circle
  
  ##--Bio-reps
  if(input$pca.plot.type=='all samples'){
    if(is.null(input$pca_color_option)){
      showmessage(text = 'Please select a factor to group and colour the samples on PCA plots.')
      return(NULL)
    }
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

observeEvent(input$update_pca_plot,{
  startmessage(text = 'Plot gene level PCA after batch effect removed')
  output$genes.pca.br.plot <- renderPlotly({
    if(is.null(pca.genes.br.g()))
      return(NULL)
    ggplotly(pca.genes.br.g())
  })
  endmessage(text = 'Plot gene level PCA after batch effect removed')
})

observeEvent(input$save_pca_br_plot,{
  ###transcript level
  # graphics.off()
  startmessage(text = 'Save PCA plots after batch effect removed')
  if(input$pca.plot.type=='all samples'){
    figure.label <- paste0('all samples ',paste0(input$pca_color_option,collapse = '_'))
  } else {
    figure.label <- paste0('average expression')
  }
  
  pca.trans.br.g <- suppressMessages(pca.trans.br.g()+
                                       coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  print(pca.trans.br.g)
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  print(pca.trans.br.g)
  dev.off()
  
  ###gene level
  # graphics.off()
  pca.genes.br.g <- suppressMessages(pca.genes.br.g()+
                                       coord_cartesian(xlim=c(input$pca.plot.x.limit.lower,input$pca.plot.x.limit.upper)))
  png(filename = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',figure.label,'.png'),
      width = input$pca.plot.width,height = input$pca.plot.height,res=input$pca.plot.res, units = 'in')
  print(pca.genes.br.g)
  dev.off()
  
  pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',figure.label,'.pdf'),
      width = input$pca.plot.width,height = input$pca.plot.height)
  print(pca.genes.br.g)
  dev.off()
  endmessage(text = 'Save PCA plots after batch effect removed')
})

##***************** preview saved plots *****************
observeEvent(input$save_pca_br_plot_view, {
  if(input$pca.plot.type=='all samples'){
    figure.label <- paste0('all samples ',paste0(input$pca_color_option,collapse = '_'))
  } else {
    figure.label <- paste0('average expression')
  }
  
  gene <- trimws(input$gene.id)
  file2save1 <- paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed ',figure.label,'.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed ',figure.label,'.png')
  
  if(!file.exists(file2save1) | !file.exists(file2save2)){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      tags$img(src=base64enc::dataURI(file = file2save1),
               style="width: 80%; display: block; margin-left: auto; margin-right: auto; 
               border:1px solid #000000"),
      p(),
      tags$img(src=base64enc::dataURI(file = file2save2),
               style="width: 80%; display: block; margin-left: auto; margin-right: auto; 
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
      ))
  }
})


##--------------Step 4-normalization--------------
##update normalisation method
observe({
  if(is.null(DDD.data$params_list$norm_method))
    return(NULL)
  updateRadioButtons(session = session,inputId = 'norm_method',selected = DDD.data$params_list$norm_method)
})

observeEvent(input$run_norm,{
  
  startmessage(text = 'Normalise transcript read counts')
  withProgress(message = 'Normalise transcript read counts...', value = 0, {
    dge <- DGEList(counts=DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                   group = DDD.data$samples_new$condition,
                   genes = DDD.data$mapping[DDD.data$target_high$trans_high,])
    incProgress(0.7)
    trans_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm_method))
    DDD.data$trans_dge <- trans_dge
    incProgress(1)
  })
  endmessage(text = 'Normalise transcript read counts')
})

observeEvent(input$run_norm,{
  startmessage(text = 'Normalise gene read counts')
  withProgress(message = 'Normalise gene read counts...', value = 0, {
    dge <- DGEList(counts=DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                   group = DDD.data$samples_new$condition)
    incProgress(0.7)
    genes_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm_method))
    DDD.data$genes_dge <- genes_dge
    incProgress(1)
  })
  endmessage(text = 'Normalise gene read counts')
})

##-------------> make the distribution plot-------------
###---transcript level
trans.dist.g <- eventReactive(input$run_norm,{
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

observeEvent(input$run_norm,{
  startmessage(text = 'Make transcript level distribution plot')
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
  endmessage(text = 'Make transcript level distribution plot')
})



###---gene level
genes.dist.g <- eventReactive(input$run_norm,{
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

observeEvent(input$run_norm,{
  startmessage(text = 'Make gene level distribution plot')
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
  endmessage(text = 'Make gene level distribution plot')
})


observeEvent(input$save_dist_plot,{
  ##transcript level
  # graphics.off()
  startmessage(text = 'Save distribution plots')
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
  endmessage(text = 'Save distribution plots')
  
  showmessage(text = paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})

##*****************preview saved plot******************
observeEvent(input$save_dist_plot_view, {
  file2save1 <- paste0(DDD.data$figure.folder,'/Transcript expression distribution.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/Gene expression distribution.png')
  
  if(!file.exists(file2save1) | !file.exists(file2save2)){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      p('Transcript level distribution plot'),
      tags$img(src=base64enc::dataURI(file = file2save1),
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
               border:1px solid #000000"),
      p(),
      p('Gene level distribution plot'),
      tags$img(src=base64enc::dataURI(file = file2save2),
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
      ))
  }
})

##***************** previous and next page ******************
observeEvent(input$page_before_preprocessing, {
  newtab <- switch(input$tabs, "preprocessing" = "generation","generation" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_preprocessing, {
  newtab <- switch(input$tabs, "preprocessing" = "ddd","ddd" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

