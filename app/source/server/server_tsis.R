# output$TSIS_time_point_cut_btn <- renderUI({
#   if(input$TSISorisokTSP=='isokTSP')
#     return(NULL)
#   selectInput('TSIS_time_point_cut',label='Min time in interval:',choices = 1:10,selected = 2,multiple = F)
# })

observe({
  if(is.null(DDD.data$samples0))
    return(NULL)
  if(input$TSISorisokTSP=='isokTSP'){
    updateSelectInput(session = session,inputId = 'TSIS_time_point_cut',choices = 1,selected = 1)
  } else {
    updateSelectInput(session = session,inputId = 'TSIS_time_point_cut',choices = 1:floor((nrow(DDD.data$samples0)-1)/2),
                      selected = 1)
  }
})

# output$tsis.sample.table <- DT::renderDataTable({
#   if(is.null(DDD.data$samples))
#     return(NULL)
#   x <- DDD.data$samples[,c(DDD.data$factor_column,'condition')]
#   x$condition <- factor(x$condition,levels = unique(x$condition))
#   x$`Numeric coordinates`<- as.numeric(x$condition)
#   
#   colnames(x)[colnames(x)=='condition'] <- 'Condition of interest'
#   values <- unique(x[,'Condition of interest'])
#   colour.n <- colorRampPalette(c("pink", "red"))(length(values))
#   x <- datatable(x)
#   formatStyle(x,'Condition of interest',
#               backgroundColor = styleEqual(levels = values,values = colour.n))
# })

output$show.spline.df <- renderUI({
  if(is.null(input$method.intersection))
    return(NULL)
  if(input$method.intersection == 'Mean' | input$TSISorisokTSP=='isokTSP')
    return(NULL)
  numericInput('spline.df',label='Spline degree:',value = ceiling(length(unique(DDD.data$samples$condition))*2/3))
})

output$show.method.intersection <- renderUI({
  if(input$TSISorisokTSP=='isokTSP'){
    selectInput('method.intersection','Search intersections',c('Mean'))
  } else {
    selectInput('method.intersection','Search intersections',c('Mean','Spline'))
  }
})



##---------->isokTSP or TSIS scoring------------
observeEvent(input$scoring,{
  contrast <- DDD.data$contrast0
  condition <- DDD.data$samples$condition
  samples <- DDD.data$samples
  DAS_genes <- DDD.data$DAS_genes
  TSIS.mapping.full <- DDD.data$mapping
  rownames(TSIS.mapping.full) <- TSIS.mapping.full$TXNAME
  TSIS.mapping.full <- TSIS.mapping.full[DDD.data$target_high$trans_high,]
  trans_TPM <- DDD.data$trans_TPM
  
  if(input$TSISorisokTSP == 'isokTSP'){
    message('Peform isokTSP analysis ...')
    if(is.null(contrast) | is.null(condition) | is.null(samples) | is.null(DAS_genes) | is.null(TSIS.mapping.full) | is.null(trans_TPM))
      return(NULL)
    
    scores.results <- list()
    for(i in seq_along(contrast)){
      contrast.idx <- contrast[i]
      samples.idx <- unlist(strsplit(contrast.idx,'-'))
      if(any((samples.idx %in% condition)==F)){
        message(paste0('IS analysis of contrast group: ',contrast.idx, ' is not applied'))
        next
      }

      message(paste0('IS analysis of contrast group: ',contrast.idx))
      col.idx <- which(condition %in% samples.idx)
      DAS.genes <- unique(DAS_genes$target[DAS_genes$contrast==contrast.idx])
      TSIS.mapping <- TSIS.mapping.full[TSIS.mapping.full$GENEID %in% DAS.genes,]
      TSIS.mapping <- TSIS.mapping[,c('GENEID','TXNAME')]
      TSIS.tpm <- trans_TPM[TSIS.mapping$TXNAME,col.idx]
      times <- samples$condition[col.idx]
      
      scores<-iso.switch.pairwise(data.exp=TSIS.tpm,
                                  mapping = TSIS.mapping,
                                  times=times,
                                  rank=F,
                                  min.difference = 1,
                                  spline = F,
                                  verbose = T,shiny = T)
      
      scores <- data.frame(contrast=contrast.idx,scores,row.names = NULL,check.names = F)
      scores.results <- c(scores.results,setNames(list(scores),contrast.idx))
    }
    scores.results <- do.call(rbind,scores.results)
    rownames(scores.results) <- NULL
    DDD.data$scores <- scores.results
  } else {
    DAS.genes <- unique(DAS_genes$target)
    TSIS.mapping <- TSIS.mapping.full[TSIS.mapping.full$GENEID %in% DAS.genes,]
    TSIS.mapping <- TSIS.mapping[,c('GENEID','TXNAME')]
    TSIS.tpm <- trans_TPM[TSIS.mapping$TXNAME,]
    times <- samples$condition
    scores.results <- iso.switch.shiny(data.exp = TSIS.tpm,
                                       mapping = TSIS.mapping,
                                       times=times,
                                       rank=F,
                                       min.t.points = 1,
                                       min.difference = 1,
                                       spline = (input$method.intersection == 'Spline'),
                                       spline.df = input$spline.df,
                                       verbose = T)
    # scores.results <- roundDF(scores,digits = 6)
    
    rownames(scores.results) <- NULL
    DDD.data$scores <- scores.results
  }
})


observeEvent(input$filtering,{
  if(is.null(DDD.data$scores))
    return(NULL)
  scores.filtered <-score.filter(
    scores = DDD.data$scores,
    prob.cutoff = input$TSIS_prob_cut,
    diff.cutoff = input$TSIS_diff_cut,
    t.points.cutoff = ifelse(input$TSISorisokTSP == 'isokTSP',1,input$TSIS_time_point_cut),
    pval.cutoff = 0.05,
    FDR.cutoff = input$TSIS_adj_pval_cut,
    cor.cutoff = input$TSIS_cor_cut,
    data.exp = NULL,
    mapping = NULL,
    sub.isoform.list = NULL,
    sub.isoform = F,
    max.ratio = F,
    x.value.limit = c(1,length(unique(DDD.data$samples$condition)))
  )
  DDD.data$scores_filtered <- scores.filtered
})

scores <- reactive({
  if(is.null(DDD.data$scores) & is.null(DDD.data$scores_filtered))
    return(NULL)
  if(!is.null(DDD.data$scores) & is.null(DDD.data$scores_filtered)){
    scores <- roundDF(DDD.data$scores,digits = 6)
  } else {
    scores <- roundDF(DDD.data$scores_filtered,digits = 6)
  }
  return(scores)
})

output$score.table <- DT::renderDataTable({
  if(is.null(scores()))
    return(NULL)
  x <- datatable(scores())
  formatStyle(x,c('x.value','before.FDR','after.FDR','before.t.points','after.t.points','prob','diff'),
              backgroundColor = 'red')
})

##---------- plot switch number ------------
ISnumber.g <- eventReactive(input$plotISnumber,{
  if(is.null(DDD.data$scores_filtered) & is.null(DDD.data$scores))
    return(NULL)
  if(is.null(DDD.data$scores_filtered)){
    scores2plot <- DDD.data$scores
  } else {
    scores2plot <- DDD.data$scores_filtered
  }
  
  if(input$TSISorisokTSP=='TSIS'){
    x <- scores2plot$x.value
    time.points <- unique(DDD.data$samples$condition)
    g <- plotTSISnumber(x = x,time.points = time.points,
                        fill = input$switch.density.color.pick,
                        angle = input$is_number_x_rotate,
                        hjust = input$is_number_x_hjust,
                        vjust = input$is_number_x_vjust
    )
  } else {
    if(!('contrast' %in% colnames(scores2plot))){
      message('Please select correct type of isoform switch.')
      return(NULL)
    }
    
    x <- data.frame(table(scores2plot$contrast))
    g <- plotPWISnumber(x,
                        fill = input$switch.density.color.pick,
                        angle = input$is_number_x_rotate,
                        hjust = input$is_number_x_hjust,
                        vjust = input$is_number_x_vjust
    )
  }
  return(g)
})


output$switch.number <- renderPlot({
  if(is.null(ISnumber.g()))
    return(NULL)
  print(ISnumber.g())
  
})

observeEvent(input$plotISnumber.save,{
  create.folders(wd = DDD.data$folder)
  folder2save <- paste0(DDD.data$folder,'/figure')
  png(paste0(folder2save,'/Isoform switch number.png'),
      width = input$is.number.width,height = input$is.number.height,res = input$is.number.res,units = 'in')
  print(ISnumber.g())
  dev.off()
  
  pdf(paste0(folder2save,'/Isoform switch number.pdf'),
      width = input$is.number.width,height = input$is.number.height)
  print(ISnumber.g())
  dev.off()
  showmessage()
})

##---------- plot isoform switch ------------
output$show.contrst.button <- renderUI({
  if(input$TSISorisokTSP=='TSIS')
    return(NULL)
  conditionalPanel(condition = "input.TSISorisokTSP == 'isokTSP'",
                   selectInput(inputId = 'select.contrast',label = 'Select a contrast group',
                               choices = DDD.data$contrast,
                               multiple = F)
  )
})

# output$show.contrst.button <- renderUI({
#   if(is.null(DDD.data$contrast) | input$TSISorisokTSP=='TSIS')
#     return(NULL)
#   selectInput(inputId = 'select.contrast',label = 'Select a contrast group',
#               choices = DDD.data$contrast,
#               multiple = F)
# })


IS.g <- eventReactive(input$plot.IS,{
  if(is.null(input$iso1) | is.null(input$iso2) | is.null(DDD.data$samples))
    return(NULL)
  if(is.null(DDD.data$scores_filtered) & is.null(DDD.data$scores))
    return(NULL)
  if(is.null(DDD.data$scores_filtered)){
    scores2plot <- DDD.data$scores
  } else {
    scores2plot <- DDD.data$scores_filtered
  }
  
  iso1 <- input$iso1
  iso2 <- input$iso2
  # iso1 <- 'AT5G65060_ID6'
  # iso2 <- 'AT5G65060_s1'
  times0 <- DDD.data$samples$condition
  data2plot <- DDD.data$trans_TPM[c(iso1,iso2),]
  
  if(input$TSISorisokTSP == 'isokTSP'){
    if(!('contrast' %in% colnames(scores2plot))){
      message('Please select correct type of isoform switch.')
      return(NULL)
    }
    contrast.idx <- input$select.contrast
    contrast.idx <- unlist(strsplit(contrast.idx,'-'))[2:1]
    times.idx <- which(times0 %in% contrast.idx)
    times0 <- times0[times.idx]
    data2plot <- data2plot[,times.idx]
    scores2plot <- scores2plot[scores2plot$contrast==input$select.contrast,]
  }
  
  if(!is.numeric(times0)){
    times <- as.numeric(factor(times0,levels = unique(times0)))
    xlabs <- paste0(unique(times0),' (x=',unique(times),')')
  } else {
    times <- times0
    xlabs <- unique(times0)
  }

  g <- suppressMessages(plotTSIS(data2plot = data2plot,
                                 scores = scores2plot,
                                 iso1 = iso1,
                                 ribbon.plot = (input$ribbon.plot=='Ribbon'),
                                 iso2 = iso2,
                                 y.lab = 'TPM',
                                 spline = (input$method.intersection == 'Spline'),
                                 spline.df = input$spline.df,
                                 times = times,
                                 errorbar.width = 0.03,
                                 x.lower.boundary = ifelse(is.numeric(times),min(times),1),
                                 x.upper.boundary = ifelse(is.numeric(times),max(times),length(unique(times))),
                                 show.scores = input$show.scores,
                                 show.region = F)+
                          scale_x_continuous(breaks = unique(times),labels = xlabs)+
                          scale_color_manual(values = c(input$color.iso1,input$color.iso2,'black'))+
                          scale_fill_manual(values = c(input$color.iso1,input$color.iso2,'black'))
  )
  return(g)
})

output$IS.plot <- renderPlot({
  if(is.null(IS.g()))
    return(NULL)
  print(IS.g())
})

observeEvent(input$save_tsis_plot,{
  create.folders(wd = DDD.data$folder)
  
  folder2save <- paste0(DDD.data$folder,'/figure/IS plots')
  if(!file.exists(folder2save))
    dir.create(path = folder2save,recursive = T)
  
  
  
  if(input$TSISorisokTSP == 'isokTSP'){
    folder2save <- paste0(DDD.data$folder,'/figure/IS plots/',input$select.contrast)
    if(!file.exists(folder2save))
      dir.create(path = folder2save,recursive = T)
  } 
  
  file2save <- paste0(folder2save,'/',paste0(input$iso1,'_vs_',input$iso2))
  png(filename = paste0(file2save,'.png'),width = input$tsis.plot.width,
      height = input$tsis.plot.height,units = 'in',res = input$tsis.plot.res)
  print(IS.g())
  dev.off()
  
  pdf(file = paste0(file2save,'.pdf'),width = input$tsis.plot.width,
      height = input$tsis.plot.height)
  print(IS.g())
  dev.off()
  showmessage('Done! Plots in pdf and png formats are saved into figure folder')
})

output$number_of_IS <- renderText({
  if(is.null(DDD.data$scores_filtered))
    return(NULL)
  paste0('There are ',nrow(DDD.data$scores_filtered), ' signficant isoform switches in totoal.')
})

# 
observeEvent(input$plot_all_IS_btn,{
  if(is.null(DDD.data$scores_filtered))
    return(NULL)
  
  create.folders(wd = DDD.data$folder)
  folder2save0 <- paste0(DDD.data$folder,'/figure/IS plots')
  if(!file.exists(folder2save0))
    dir.create(path = folder2save0,recursive = T)
  
  scores2plot <- DDD.data$scores_filtered
  withProgress(message = 'Plot isoform switch...', value = 0, {
    for(i in 1:nrow(scores2plot)){
      message(paste0(i,' of ',nrow(scores2plot)))
      incProgress(1/nrow(scores2plot))
      iso1 <- scores2plot$iso1[i]
      iso2 <- scores2plot$iso2[i]
      times0 <- DDD.data$samples$condition
      data2plot <- DDD.data$trans_TPM[c(iso1,iso2),]
      
      folder2save <- folder2save0
      if(input$TSISorisokTSP == 'isokTSP'){
        contrast.idx <- scores2plot$contrast[i]
        ###
        folder2save <- paste0(folder2save0,'/',contrast.idx)
        if(!file.exists(folder2save))
          dir.create(path = folder2save,recursive = T)
        
        contrast.idx <- unlist(strsplit(contrast.idx,'-'))[2:1]
        times.idx <- which(times0 %in% contrast.idx)
        times0 <- times0[times.idx]
        data2plot <- data2plot[,times.idx]
      }
      
      if(!is.numeric(times0)){
        times <- as.numeric(factor(times0,levels = unique(times0)))
        xlabs <- paste0(unique(times0),' (x=',unique(times),')')
      } else {
        times <- times0
        xlabs <- unique(times0)
      }
      
      g <- suppressMessages(plotTSIS(data2plot = data2plot,
                                     scores = scores2plot[i,,drop=F],
                                     iso1 = iso1,
                                     ribbon.plot = (input$ribbon.plot=='Ribbon'),
                                     iso2 = iso2,
                                     y.lab = 'TPM',
                                     spline = (input$method.intersection == 'Spline'),
                                     spline.df = input$spline.df,
                                     times = times,
                                     errorbar.width = 0.03,
                                     x.lower.boundary = ifelse(is.numeric(times),min(times),1),
                                     x.upper.boundary = ifelse(is.numeric(times),max(times),length(unique(times))),
                                     show.scores = input$show.scores,
                                     show.region = F)+
                              scale_x_continuous(breaks = unique(times),labels = xlabs)+
                              scale_color_manual(values = c(input$color.iso1,input$color.iso2,'black'))+
                              scale_fill_manual(values = c(input$color.iso1,input$color.iso2,'black'))
                            )
      
      file2save <- paste0(folder2save,'/',paste0(iso1,'_vs_',iso2))
      
      if(input$tsis.multiple.plot.format=='all' | input$tsis.multiple.plot.format=='png'){
        png(filename = paste0(file2save,'.png'),width = input$tsis.all.plot.width,
            height = input$tsis.all.plot.height,units = 'in',res = input$tsis.all.plot.res)
        print(g)
        dev.off()
      }
      
      if(input$tsis.multiple.plot.format=='all' | input$tsis.multiple.plot.format=='pdf'){
        pdf(file = paste0(file2save,'.pdf'),width = input$tsis.all.plot.width,
            height = input$tsis.all.plot.height)
        print(g)
        dev.off()
      }
    }
  })
  showmessage('Done!')
  
})

observeEvent(input$page_before_tsis, {
  newtab <- switch(input$tabs, "TSIS" = "function","function" = "TSIS")
  updateTabItems(session, "tabs", newtab)
})

observeEvent(input$page_after_tsis, {
  newtab <- switch(input$tabs, "report" = "TSIS","TSIS" = "report")
  updateTabItems(session, "tabs", newtab)
})

