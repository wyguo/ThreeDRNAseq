##--------------- >>Step 1: Select factors <<-----------------
# output$correct_condition_select_btn <- renderUI({
#   if(input$correct_condition=='Yes')
#     return(NULL)
#   selectInput(inputId = 'factor_column_new',
#               label = 'Select the correct factors of interest',
#               selected = NULL,multiple = T,width = '500px',
#               choices = colnames(DDD.data$samples0))
# })
# 
# output$refresh_condition_select_ui <- renderUI({
#   if(input$correct_condition=='Yes')
#     return(NULL)
#   actionButton('refresh_condition_select','Refresh',icon("refresh"),
#                style="color: #fff; background-color: #428bca; border-color: #2e6da4")
# })
# 
# observeEvent(input$refresh_condition_select,{
#   if(is.null(DDD.data$samples0))
#     return(NULL)
#   if(input$correct_condition=='No'){
#     DDD.data$factor_column <- input$factor_column_new
#     DDD.data$samples$condition <- as.vector(interaction(DDD.data$samples[,input$factor_column_new]))
#     DDD.data$samples_new$condition <- as.vector(interaction(DDD.data$samples_new[,input$factor_column_new]))
#   } else {
#     return(NULL)
#   }
# })

selected.samples <- reactive({
  if(is.null(DDD.data$samples_new) | all(is.null(DDD.data$factor_column)))
    return(NULL)
  x <- DDD.data$samples_new[,DDD.data$factor_column,drop=F]
  if(length(DDD.data$factor_column)>1){
    x$Labels <- make.names(interaction(x))
  } else {
    x$Labels  <- make.names(as.vector(t(x)))
  }
  x[sapply(x,is.character)] <- lapply(x[sapply(x,is.character)],as.factor)
  # x <- x[order(x$`Condition to contrast`),]
  return(x)
})

observe({
  if(is.null(DDD.data$samples_new) | is.null(DDD.data$factor_column))
    return(NULL)
  DDD.data$conditions <- as.vector(t(selected.samples()$Labels))
})

output$factor_table <- DT::renderDataTable({
  if(is.null(selected.samples()))
    return(NULL)
  x <- selected.samples()
  column.idx <- colnames(x)
  n <- length(column.idx)
  colour.from <- gg.color.hue(n)
  values <- lapply(1:ncol(x),function(i) unique(x[,i]))
  values <- lapply(x,unique)
  values.n <- sapply(values,length)
  colour.n <- lapply(1:length(colour.from),function(i){
    sub.colour <- colorRampPalette(c(colour.from[i],'white'))(values.n[i]+1)[1:values.n[i]]
    sub.colour <- paste0(sub.colour,'70')
  })
  names(colour.n) <- column.idx
  # DT::datatable(x)
  datatable(x,options = list(pageLength = 100)) %>% formatStyle(
    names(x),
    backgroundColor = styleEqual(levels = unlist(values),values = unlist(colour.n))
  )
})

##--------------- >>Step 2: Set contrast groups <<-----------------
observe({
  callModule(module = selectorServer, id = 1, 
             thisList = DDD.data$conditions)
  updateSelectInput(session = session,inputId = NS(1)('Treatment'),
                    choices = DDD.data$conditions,selected = NULL)
  updateSelectInput(session = session,inputId = NS(1)('Control'),
                    choices = DDD.data$conditions,selected = NULL)
})

params <- reactiveValues(btn = 1)
contrastIdx <- reactiveValues(n=1)
# 
observeEvent(input$generate_contrast,{
  x <- lapply(contrastIdx$n,function(i){
    x1 <- input[[NS(i, "Treatment")]]
    if(length(x1)>1)
      x1 <- paste0('(',paste0(x1,collapse = '+'),')/',length(x1))
    x2 <- input[[NS(i, "Control")]]
    if(length(x2)>1)
      x2 <- paste0('(',paste0(x2,collapse = '+'),')/',length(x2))
    data.frame(Treatment=x1,
               Control=x2)
  })
  x <- data.frame(do.call(rbind,x))
  colnames(x) <- c('Treatment','Control')
  x <- x[!duplicated(x),]
  x <- na.omit(x)
  idx <- which(as.vector(t(x[,1]))==as.vector(t(x[,2])))
  if(length(idx)>0)
    x <- x[-idx,,drop=F]
  if(nrow(x)==0)
    return(NULL)
  contrast <- paste0(x$Treatment,'-',x$Control)
  DDD.data$contrast0 <- contrast
})

observeEvent(input$insertParamBtn, {
  n <- contrastIdx$n
  n <- c(n,length(n)+1)
  contrastIdx$n <- n
  
  params$btn <- params$btn + 1
  callModule(selectorServer, params$btn,
             thisList = DDD.data$conditions)
  insertUI(
    selector = '#placeholder',
    ui = selectorUI(params$btn)
  )
})

observeEvent(input$removeParamBtn, {
  removeUI(
    ## pass in appropriate div id
    selector = paste0('#param', params$btn)
  )
  params$btn <- params$btn - 1
  n <- contrastIdx$n
  n <- n[-length(n)]
  contrastIdx$n <- n
})

observeEvent(input$refresh_contrast,{
  if(params$btn==1)
    return(NULL)
  for(i in params$btn:2){
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#param', i)
    )
    params$btn <- params$btn - 1
    n <- contrastIdx$n
    n <- n[-length(n)]
    contrastIdx$n <- n
  }

})

## contrast of contrast groups is disabled
# ##--------------->contrast between contrast groups-----------------
# observe({
#   if(input$contrast_of_contrast=='Yes'){
#     updateButton(session, "insertParamBtn_cofc", disabled = F,icon = icon('plus'),
#                  style="primary")
#     updateButton(session, "removeParamBtn_cofc", disabled = F,icon = icon('minus'),
#                  style="primary")
#     updateButton(session, "refresh_contrast_cofc", disabled = F,icon = icon('refresh'),
#                  style="primary")
#   } else {
#     updateButton(session, "insertParamBtn_cofc", disabled = T,icon = icon('ban'))
#     updateButton(session, "removeParamBtn_cofc", disabled = T,icon = icon('ban'))
#     updateButton(session, "refresh_contrast_cofc", disabled = T,icon = icon('ban'))
#   }
# })
# 
# observe({
#   if(input$contrast_of_contrast=='No')
#     return(NULL)
#   callModule(module = selectorServer2, id = 1,
#              thisList = DDD.data$contrast0)
#   updateSelectInput(session = session,inputId = NS(1)('Contrast2'),
#                     choices = DDD.data$contrast0,selected = NULL)
#   updateSelectInput(session = session,inputId = NS(1)('Contrast1'),
#                     choices = DDD.data$contrast0,selected = NULL)
# })
# 
# params_cofc <- reactiveValues(btn = 1)
# contrastIdx_cofc <- reactiveValues(n=1)
# #
# observeEvent(input$generate_contrast_cofc,{
#   x <- lapply(contrastIdx_cofc$n,function(i){
#     x1 <- input[[NS(i, "Contrast2")]]
#     x2 <- input[[NS(i, "Contrast1")]]
#     data.frame(Contrast1=x1,
#                Contrast2=x2)
#   })
#   x <- data.frame(do.call(rbind,x))
#   colnames(x) <- c('Treatment','Control')
#   x <- x[!duplicated(x),]
#   x <- na.omit(x)
#   idx <- which(as.vector(t(x[,1]))==as.vector(t(x[,2])))
#   if(length(idx)>0)
#     x <- x[-idx,,drop=F]
#   if(nrow(x)==0)
#     return(NULL)
#   contrast <- paste0('(',x$Treatment,')-(',x$Control,')')
#   DDD.data$contrast_cofc <- contrast
# })
# 
# 
# observeEvent(input$insertParamBtn_cofc, {
#   n <- contrastIdx_cofc$n
#   n <- c(n,length(n)+1)
#   contrastIdx_cofc$n <- n
# 
#   params_cofc$btn <- params_cofc$btn + 1
#   callModule(selectorServer2, params_cofc$btn,
#              thisList = DDD.data$contrast0)
#   insertUI(
#     selector = '#placeholder_cofc',
#     ui = selectorUI2(params_cofc$btn)
#   )
# })
# 
# observeEvent(input$removeParamBtn_cofc, {
#   removeUI(
#     ## pass in appropriate div id
#     selector = paste0('#param_cofc', params_cofc$btn)
#   )
#   params_cofc$btn <- params_cofc$btn - 1
#   n <- contrastIdx_cofc$n
#   n <- n[-length(n)]
#   contrastIdx_cofc$n <- n
# })
# 
# observeEvent(input$refresh_contrast,{
#   if(params_cofc$btn==1)
#     return(NULL)
#   for(i in params_cofc$btn:2){
#     removeUI(
#       ## pass in appropriate div id
#       selector = paste0('#param_cofc', i)
#     )
#     params_cofc$btn <- params_cofc$btn - 1
#     n <- contrastIdx_cofc$n
#     n <- n[-length(n)]
#     contrastIdx_cofc$n <- n
#   }
# })
# 
# observeEvent(input$refresh_contrast,{
#   DDD.data$contrast <- NULL
#   DDD.data$contrast_cofc <- NULL
#   DDD.data$contrast0 <- NULL
# })
# 
# observe({
#   DDD.data$contrast <- c(DDD.data$contrast0,DDD.data$contrast_cofc)
# })

observe({
  DDD.data$contrast <- c(DDD.data$contrast0)
})

## show contrast table
output$view_contrast_table <- DT::renderDataTable({
  x <- DDD.data$contrast
  if(is.null(x))
    return(NULL)
  x <- data.frame(Contrast=paste0('Contrast group',1:length(x)),
                  `Treatment-Control`= x,check.names = F)
  datatable(x,options=list(searching=F,info=F))
})



##--------------->>Step 3: 3D analysis<<-----------------

observeEvent(input$run.DE,{
  cat('\nDE gene analysis\n')
  ##---------------> DE genes -----------------
  withProgress(message = 'DE gene analysis...', detail = 'This may take a while...',
               value = 0, {
                 if(is.null(DDD.data$genes_batch)){
                   batch.effect <- NULL
                 } else {
                   batch.effect <- DDD.data$genes_batch$W
                 }
                 incProgress(0.1)
               
                 
                 design <- condition2design(condition = factor(DDD.data$conditions,
                                                               levels = unique(DDD.data$conditions)),
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
                                                          adjust.method = input$pval_adj_method,
                                                          block = DDD.data$block)
                        },
                        glmQL={
                          genes_3D_stat <- edgeR.pipeline(dge = DDD.data$genes_dge,
                                                          design = design,
                                                          deltaPS = DDD.data$deltaPS,
                                                          contrast = DDD.data$contrast,
                                                          diffAS = F,
                                                          method = 'glmQL',
                                                          adjust.method = input$pval_adj_method)
                        },
                        glm={
                          genes_3D_stat <- edgeR.pipeline(dge = DDD.data$genes_dge,
                                                          design = design,
                                                          deltaPS = NULL,
                                                          contrast = DDD.data$contrast,
                                                          diffAS = F,
                                                          method = 'glm',
                                                          adjust.method = input$pval_adj_method)
                        }
                 )
                 DDD.data$genes_3D_stat <- genes_3D_stat
                 
                 incProgress(0.7)
                 DE_genes <- summaryDEtarget(stat = genes_3D_stat$DE.stat,
                                             cutoff = c(adj.pval=input$pval_cut,
                                                        log2FC=input$l2fc_cut))
                 
                 DDD.data$DE_genes <- DE_genes
                 DDD.data$genes_log2FC <- genes_3D_stat$DE.lfc
                 # save(DE_genes,file=paste0(DDD.data$data.folder,'/DE_genes.RData'))
                 incProgress(1)
                 # save(genes_3D_stat,file=paste0(DDD.data$data.folder,'/genes_3D_stat.RData'))
                 # message(paste0('genes_3D_stat.RData is saved in folder: ',DDD.data$data.folder))
               })
  
  ##---------------> Generating deltaPS -----------------
  message('Generating deltaPS...')
  condition <- as.vector(interaction(DDD.data$samples0[,DDD.data$factor_column]))
  mapping <- DDD.data$mapping
  rownames(mapping) <- mapping$TXNAME
  deltaPS <- transAbundance2PS(transAbundance = DDD.data$trans_TPM[DDD.data$target_high$trans_high,],
                               PS = DDD.data$PS,
                               contrast = DDD.data$contrast,
                               condition = condition,
                               mapping = mapping[DDD.data$target_high$trans_high,])
  DDD.data$PS <- deltaPS$PS
  DDD.data$deltaPS <- deltaPS$deltaPS
  
  ##---------------> DAS genes, DE/DTU transcripts -----------------
  withProgress(message = 'DE gene analysis...', detail = 'This may take a while...',
               value = 0, {
                 cat('\nDAS gene, DE and DTU transcript analysis\n')
                 if(is.null(DDD.data$trans_batch)){
                   batch.effect <- NULL
                 } else {
                   batch.effect <- DDD.data$trans_batch$W
                 }
                 incProgress(0.1)
                 design <- condition2design(condition = factor(DDD.data$conditions,levels = unique(DDD.data$conditions)),
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
                                                          adjust.method = input$pval_adj_method,
                                                          block = DDD.data$block)
                        },
                        glmQL={
                          trans_3D_stat <- edgeR.pipeline(dge = DDD.data$trans_dge,
                                                          design = design,
                                                          deltaPS = DDD.data$deltaPS,
                                                          contrast = DDD.data$contrast,
                                                          diffAS = T,
                                                          method = 'glmQL',
                                                          adjust.method = input$pval_adj_method)
                        },
                        glm={
                          trans_3D_stat <- edgeR.pipeline(dge = DDD.data$trans_dge,
                                                          design = design,
                                                          deltaPS = DDD.data$deltaPS,
                                                          contrast = DDD.data$contrast,
                                                          diffAS = T,
                                                          method = 'glm',
                                                          adjust.method = input$pval_adj_method)
                        }
                 )
                 DDD.data$trans_3D_stat <- trans_3D_stat
                
                 incProgress(0.7)
                 ##DE trans
                 DE_trans <- summaryDEtarget(stat = trans_3D_stat$DE.stat,
                                             cutoff = c(adj.pval=input$pval_cut,
                                                        log2FC=input$l2fc_cut))
                 DDD.data$DE_trans <- DE_trans
                 DDD.data$trans_log2FC <- trans_3D_stat$DE.lfc
                 # save(DE_trans,file=paste0(DDD.data$data.folder,'/DE_trans.RData'))
                 
                 ##DAS genes
                 if(input$DAS_pval_method=='F-test')
                   DAS.stat <- trans_3D_stat$DAS.F.stat else DAS.stat <- trans_3D_stat$DAS.Simes.stat
                 
                 lfc <- DDD.data$genes_log2FC
                 lfc <- reshape2::melt(as.matrix(lfc))
                 colnames(lfc) <- c('target','contrast','log2FC')
                 DAS_genes <- summaryDAStarget(stat = DAS.stat,
                                               lfc = lfc,
                                               cutoff=c(input$pval_cut,input$deltaPS_cut))
                 DDD.data$DAS_genes <- DAS_genes
                 # save(DAS_genes,file=paste0(DDD.data$data.folder,'/DAS_genes.RData'))
                 
                 ##DTU  trans
                 lfc <- DDD.data$trans_log2FC
                 lfc <- reshape2::melt(as.matrix(lfc))
                 colnames(lfc) <- c('target','contrast','log2FC')
                 DTU_trans <- summaryDAStarget(stat = trans_3D_stat$DTU.stat,
                                               lfc = lfc,cutoff = c(adj.pval=input$pval_cut,
                                                                    deltaPS=input$deltaPS_cut))
                 DDD.data$DTU_trans <- DTU_trans
                 # save(DTU_trans,file=paste0(DDD.data$data.folder,'/DTU_trans.RData'))
                 incProgress(1)
                 # 
                 # save(trans_3D_stat,file=paste0(DDD.data$data.folder,'/trans_3D_stat.RData'))
                 # 
                 # message(paste0('trans_3D_stat.RData is saved in folder: ',DDD.data$data.folder))
               })
  
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

##--------------->DDD up- and down-regulation-----------------
# output$up.down.bar.plot <- renderPlotly({
g.updown <- reactive({
  if(is.null((DDD.data$DE_genes)) | 
     is.null((DDD.data$DE_trans)) | 
     is.null((DDD.data$contrast)))
    return(NULL)
  idx <- c('DE_genes','DE_trans')
  title.idx <- gsub('_',' ',idx)
  title.idx <- gsub('trans','transcripts',title.idx)
  names(title.idx) <- idx
  g <- lapply(idx,function(i){
    subidx <- factor(DDD.data[[i]]$contrast,levels = DDD.data$contrast)
    targets <-  split(DDD.data[[i]],subidx)
    data2plot <- lapply(DDD.data$contrast,function(i){
      if(nrow(targets[[i]])==0){
        x <- data.frame(contrast=i,regulation=c('down-regulated','up-regulated'),number=0,check.names = F)
      } else {
        x <- data.frame(contrast=i,table(targets[[i]]$up.down))
        colnames(x) <- c('contrast','regulation','number')
      }
      x
    })
    data2plot <- do.call(rbind,data2plot)
    plotUpdown(data2plot,plot.title = title.idx[i],
               col.up = input$color4up_regulate,
               col.down = input$color4down_regulate,
               contrast = DDD.data$contrast,
               angle = input$updown.x.rotate,
               hjust = input$updown.x.hjust,
               vjust = input$updown.x.vjust)
  })
  names(g) <- title.idx
  g
})

output$up.down.bar.plot.DE.genes <- renderPlot({
  if(is.null((g.updown())))
    return(NULL)
  print(g.updown()[[1]])
})

# output$up.down.bar.plot.DAS.genes <- renderPlot({
#   if(is.null((g.updown())))
#     return(NULL)
#   print(g.updown()[[2]])
# })

output$up.down.bar.plot.DE.trans <- renderPlot({
  if(is.null((g.updown())))
    return(NULL)
  print(g.updown()[[2]])
})

# output$up.down.bar.plot.DTU.trans <- renderPlot({
#   if(is.null((g.updown())))
#     return(NULL)
#   print(g.updown()[[4]])
# })

observeEvent(input$save.up.down.bar.plot,{
  withProgress(message = 'Saving plots ...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0)
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
                 incProgress(0.8)
                 message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
                 showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
                 incProgress(1)
               })
})

##---------------Comparisons between contrast groups-----------------
observe({
  updateSelectInput(
    session = session,
    inputId = 'across.contrast.group',
    label = 'Contrast groups (multiple selection allowed)',
    selected = DDD.data$contrast[1:2],
    choices = DDD.data$contrast
  )
})

observe({
  updateSelectInput(
    session = session,
    inputId = 'across.target',
    label = 'Contrast groups',
    selected = DDD.data$contrast[1],
    choices = DDD.data$contrast
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
                  pval.cutoff = input$pval_cut,lfc.cutoff = input$l2fc_cut,deltaPS.cutoff = input$deltaPS_cut)
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
                  pval.cutoff = input$pval_cut,lfc.cutoff = input$l2fc_cut,
                  deltaPS.cutoff = input$deltaPS_cut)
  }
  trans.flow.chart
})
output$trans.flow.chart <- renderPlot({
  trans.flow.chart()()
})

observeEvent(input$save.union.flow.chart,{
  ##transcript level
  withProgress(message = 'Saving plots ...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0)
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
                 incProgress(0.5)
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
                 incProgress(0.9)
                 message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
                 showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
                 incProgress(1)
               })
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
    col2fill <- colorRampPalette(c(input$color4_3d_number_from,input$color4_3d_number_to))(length(targets))
    g <- plotEulerDiagram(x = targets,
                          scale.plot = input$across.contrast.euler.plot.scale,
                          fill = col2fill)
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
  grid.arrange(g.across.contrast()[[3]],top=textGrob('DE transcripts', gp=gpar(cex=1.2)))
})

###---DTU transcripts
output$DTU.trans.euler <- renderPlot({
  if(is.null(g.across.contrast()))
    return(NULL)
  grid.arrange(g.across.contrast()[[4]],top=textGrob('DTU transcripts', gp=gpar(cex=1.2)))
})

observeEvent(input$save.across.contrast.euler.plot,{
  ###save DE.genes
  lapply(names(g.across.contrast()),function(i){
    figure.name <- paste0(DDD.data$figure.folder,'/',i,' euler plot across contrast ',paste0(input$across.contrast.group,collapse = ' vs '))
    png(paste0(figure.name,'.png'),
        width = input$across.contrast.euler.plot.width,
        height = input$across.contrast.euler.plot.height,
        res=input$across.contrast.euler.plot.res, units = 'in')
    grid.arrange(g.across.contrast()[[i]],top=textGrob(i, gp=gpar(cex=1.2)))
    dev.off()
    pdf(paste0(figure.name,'.pdf'),
        width = input$across.contrast.euler.plot.width,height = input$across.contrast.euler.plot.height)
    grid.arrange(g.across.contrast()[[i]],top=textGrob(i, gp=gpar(cex=1.2)))
    dev.off()
  })
  
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  
})

###target numbers
observe({
  if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans)| is.null(DDD.data$DTU_trans))
    return(NULL)
  DDD.data$DDD_numbers <- summary3Dnumber(DE_genes = DDD.data$DE_genes,
                                            DAS_genes = DDD.data$DAS_genes,
                                            DE_trans = DDD.data$DE_trans,
                                            DTU_trans = DDD.data$DTU_trans,contrast = DDD.data$contrast)
  
})

output$DDD.numbers <- renderTable({
  DDD.data$DDD_numbers
},align='c')

##---------------DE vs DAS-----------------
observe({
  if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$contrast))
    return(NULL)
  # summary.3D.vs(x = split(DDD.data$DE_genes$target,DDD.data$DE_genes$contrast),
  #                y = split(DDD.data$DAS_genes$target,DDD.data$DAS_genes$contrast),
  #                contrast = DDD.data$contrast,
  #                idx = c('DEonly','DE&DAS','DASonly'))
  DDD.data$DEvsDAS_results <- DEvsDAS(DE_genes = DDD.data$DE_genes,
                                        DAS_genes = DDD.data$DAS_genes,
                                        contrast = DDD.data$contrast)
  
})

output$DE.vs.DAS <- renderTable({
  DDD.data$DEvsDAS_results
},align='c')

##euler plot
###---DE vs DAS
g.across.target1 <- reactive({
  if(is.null(DDD.data$DEvsDAS_results)| is.null(input$across.target))
    return(NULL)
  x <- unlist(DDD.data$DEvsDAS_results[DDD.data$DEvsDAS_results$Contrast==input$across.target,-1])
  if(length(x)==0)
    return(NULL)
  names(x) <- c('DE','DE&DAS','DAS')
  g <- plotEulerDiagram(x = x,
                        fill = c(input$color4_dedas_number_from,input$color4_dedas_number_to),
                        scale.plot = input$across.target.euler.plot.scale)
  g
})

output$DEvsDAS.euler <- renderPlot({
  if(is.null(g.across.target1()))
    return(NULL)
  grid.arrange(g.across.target1(),top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',input$across.target), 
                                               gp=gpar(cex=1.2)))
})


##---------------DE vs DTU-----------------
observe({
  if(is.null(DDD.data$DE_trans) | is.null(DDD.data$DTU_trans) | is.null(DDD.data$contrast))
    return(NULL)
  DDD.data$DEvsDTU_results <- DEvsDTU(DE_trans = DDD.data$DE_trans,
                                        DTU_trans = DDD.data$DTU_trans,
                                        contrast = DDD.data$contrast)
})

output$DE.vs.DTU <- renderTable({
  DDD.data$DEvsDTU_results
},align='c')

###---DE vs DTU
g.across.target2 <- reactive({
  if(is.null(DDD.data$DEvsDTU_results)| is.null(input$across.target))
    return(NULL)
  x <- unlist(DDD.data$DEvsDTU_results[DDD.data$DEvsDTU_results$Contrast==input$across.target,-1])
  if(length(x)==0)
    return(NULL)
  names(x) <- c('DE','DE&DTU','DTU')
  g <- plotEulerDiagram(x = x,fill = gg.color.hue(2),scale.plot = input$across.target.euler.plot.scale)
  g
})

output$DEvsDTU.euler <- renderPlot({
  if(is.null(g.across.target2()))
    return(NULL)
  grid.arrange(g.across.target2(),top=textGrob(paste0('DE vs DTU transcripts\nin contrast group: ',input$across.target), 
                                               gp=gpar(cex=1.2)))
})

observeEvent(input$save.across.target.euler.plot,{
  
  g.across.target1 <- reactive({
    if(is.null(DDD.data$DEvsDAS_results)| is.null(input$across.target))
      return(NULL)
    
  })
  if(is.null(DDD.data$contrast) | is.null(DDD.data$DEvsDAS_results)| is.null(input$across.target))
    return(NULL)
  
  withProgress(message = paste0('Saving plots ...'),
               detail = 'This may take a while...', value = 0, {
                 for(i in DDD.data$contrast){
                   ##DE vs DAS genes
                   incProgress(1/length(DDD.data$contrast))
                   x <- unlist(DDD.data$DEvsDAS_results[DDD.data$DEvsDAS_results$Contrast==i,-1])
                   if(length(x)==0)
                     next
                   names(x) <- c('DE','DE&DAS','DAS')
                   g.across.target1 <- plotEulerDiagram(x = x,
                                                        fill = c(input$color4_dedas_number_from,input$color4_dedas_number_to),
                                                        scale.plot = input$across.target.euler.plot.scale)
                   figure.name <- paste0(DDD.data$figure.folder,'/DE vs DAS gene euler plot in contrast ',gsub('/','over',i))
                   png(paste0(figure.name,'.png'),
                       width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height,
                       res=input$across.target.euler.plot.res, units = 'in')
                   grid.arrange(g.across.target1,top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',i), 
                                                              gp=gpar(cex=1.2)))
                   dev.off()
                   pdf(paste0(figure.name,'.pdf'),
                       width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height)
                   grid.arrange(g.across.target1,top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',i), 
                                                              gp=gpar(cex=1.2)))
                   dev.off()
                 }
               })
    
    ##DE vs DTU transcript
  withProgress(message = paste0('Plotting venn diagrams...'),
               detail = 'This may take a while...', value = 0, {
                 for(i in DDD.data$contrast){ 
                   incProgress(1/length(DDD.data$contrast))
                   x <- unlist(DDD.data$DEvsDTU_results[DDD.data$DEvsDTU_results$Contrast==i,-1])
                   if(length(x)==0)
                     next
                   names(x) <- c('DE','DE&DTU','DTU')
                   g.across.target2 <- plotEulerDiagram(x = x,fill = gg.color.hue(2),
                                                        scale.plot = input$across.target.euler.plot.scale)
                   figure.name <- paste0(DDD.data$figure.folder,'/DE vs DTU transcript euler plot in contrast ', gsub('/','over',i))
                   png(paste0(figure.name,'.png'),
                       width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height,
                       res=input$across.target.euler.plot.res, units = 'in')
                   grid.arrange(g.across.target2,top=textGrob(paste0('DE vs DTU transcripts\nin contrast group: ',i), 
                                                              gp=gpar(cex=1.2)))
                   dev.off()
                   
                   pdf(paste0(figure.name,'.pdf'),
                       width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height)
                   grid.arrange(g.across.target2,top=textGrob(paste0('DE vs DTU transcripts\nin contrast group: ',i), 
                                                              gp=gpar(cex=1.2)))
                   dev.off()
                 }
               })
  
  ##DE vs DAS genes
  figure.name <- paste0(DDD.data$figure.folder,'/DE vs DAS gene euler plot in contrast ',input$across.target)
  png(paste0(figure.name,'.png'),
      width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height,
      res=input$across.target.euler.plot.res, units = 'in')
  grid.arrange(g.across.target1(),top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',input$across.target), 
                                               gp=gpar(cex=1.2)))
  dev.off()
  pdf(paste0(figure.name,'.pdf'),
      width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height)
  grid.arrange(g.across.target1(),top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',input$across.target), 
                                               gp=gpar(cex=1.2)))
  dev.off()
  
  
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})


observeEvent(input$page_before_ddd, {
  newtab <- switch(input$tabs, "preprocessing" = "ddd","ddd" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
})

observeEvent(input$page_after_ddd, {
  newtab <- switch(input$tabs, "function" = "ddd","ddd" = "function")
  updateTabItems(session, "tabs", newtab)
})
