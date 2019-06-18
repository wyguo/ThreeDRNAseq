##--------------- >>Step 1: Select factors <<-----------------
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
  rownames(x) <- NULL
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
##--------- + Difference of pair-wise----
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

##insert a line
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

##remove a line
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

# 
observeEvent(input$generate_contrast,{
  DDD.data$contrast_pw <- NULL
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
  DDD.data$contrast_pw <- unique(paste0(x$Treatment,'-',x$Control))
})

##--------- + Difference of multiple group mean----
observe({
  if('Difference of group mean' %in% input$contrast_type)
    thisList = DDD.data$conditions else thisList = NULL
  callModule(module = selectorServer_mean, id = 1, 
             thisList = thisList)
  updateSelectInput(session = session,inputId = NS(1)('Treatment_mean'),
                    choices = thisList,selected = NULL)
  updateSelectInput(session = session,inputId = NS(1)('Control_mean'),
                    choices = thisList,selected = NULL)
})

params_mean <- reactiveValues(btn = 1)
contrastIdx_mean <- reactiveValues(n=1)

##insert a line
observe({
  if('Difference of group mean' %in% input$contrast_type){
    updateButton(session, "insertParamBtn_mean", disabled = F,icon = icon('plus'),
                 style="primary")
    updateButton(session, "removeParamBtn_mean", disabled = F,icon = icon('minus'),
                 style="primary")
  } else {
    updateButton(session, "insertParamBtn_mean", disabled = T,icon = icon('ban'))
    updateButton(session, "removeParamBtn_mean", disabled = T,icon = icon('ban'))
  }
})

observeEvent(input$insertParamBtn_mean, {
  n <- contrastIdx_mean$n
  n <- c(n,length(n)+1)
  contrastIdx_mean$n <- n
  params_mean$btn <- params_mean$btn + 1
  callModule(selectorServer_mean, params_mean$btn,
             thisList = DDD.data$conditions)
  insertUI(
    selector = '#placeholder_mean',
    ui = selectorUI_mean(params_mean$btn)
  )
})

##remove a line
observeEvent(input$removeParamBtn_mean, {
  removeUI(
    ## pass in appropriate div id
    selector = paste0('#param_mean', params_mean$btn)
  )
  params_mean$btn <- params_mean$btn - 1
  n <- contrastIdx_mean$n
  n <- n[-length(n)]
  contrastIdx_mean$n <- n
})

# 
observeEvent(input$generate_contrast,{
  DDD.data$contrast_mean <- NULL
  if(!('Difference of group mean' %in% input$contrast_type))
    return(NULL)
  
  x <- lapply(contrastIdx_mean$n,function(i){
    x1 <- input[[NS(i, "Treatment_mean")]]
    x2 <- input[[NS(i, "Control_mean")]]
    if(length(x1)==0 | length(x2)==0){
      NA
    } else {
      if(length(x1)>1)
        x1 <- paste0('(',paste0(x1,collapse = '+'),')/',length(x1))
      if(length(x2)>1)
        x2 <- paste0('(',paste0(x2,collapse = '+'),')/',length(x2))
      data.frame(Treatment=x1, Control=x2)
    }
  })
  x <- x[!is.na(x)]
  if(length(x)==0)
    return(NULL)
  x <- data.frame(do.call(rbind,x))
  colnames(x) <- c('Treatment','Control')
  x <- x[!duplicated(x),]
  x <- na.omit(x)
  idx <- which(as.vector(t(x[,1]))==as.vector(t(x[,2])))
  if(length(idx)>0)
    x <- x[-idx,,drop=F]
  if(nrow(x)==0)
    return(NULL)
  DDD.data$contrast_mean <- unique(paste0(x$Treatment,'-',x$Control))
})

##--------- + Difference of pair-wise group diffence----
observe({
  if('Difference of group difference' %in% input$contrast_type)
    thisList = DDD.data$conditions else thisList = NULL
    
  callModule(module = selectorServer_pgdiff, id = 1, 
             thisList = thisList)
  updateSelectInput(session = session,inputId = NS(1)('Treatment_pgdiff'),
                    choices = thisList,selected = NULL)
  updateSelectInput(session = session,inputId = NS(1)('Control_pgdiff'),
                    choices =thisList,selected = NULL)
})

params_pgdiff <- reactiveValues(btn = 1)
contrastIdx_pgdiff <- reactiveValues(n=1)

##insert a line
observe({
  if('Difference of group difference' %in% input$contrast_type){
    updateButton(session, "insertParamBtn_pgdiff", disabled = F,icon = icon('plus'),
                 style="primary")
    updateButton(session, "removeParamBtn_pgdiff", disabled = F,icon = icon('minus'),
                 style="primary")
  } else {
    updateButton(session, "insertParamBtn_pgdiff", disabled = T,icon = icon('ban'))
    updateButton(session, "removeParamBtn_pgdiff", disabled = T,icon = icon('ban'))
  }
})

observeEvent(input$insertParamBtn_pgdiff, {
  n <- contrastIdx_pgdiff$n
  n <- c(n,length(n)+1)
  contrastIdx_pgdiff$n <- n
  params_pgdiff$btn <- params_pgdiff$btn + 1
  callModule(module = selectorServer_pgdiff, id = params_pgdiff$btn,
             thisList = DDD.data$conditions)
  insertUI(
    selector = '#placeholder_pgdiff',
    ui = selectorUI_pgdiff(params_pgdiff$btn)
  )
})

##remove a line
observeEvent(input$removeParamBtn_pgdiff, {
  removeUI(
    ## pass in appropriate div id
    selector = paste0('#param_pgdiff', params_pgdiff$btn)
  )
  params_pgdiff$btn <- params_pgdiff$btn - 1
  n <- contrastIdx_pgdiff$n
  n <- n[-length(n)]
  contrastIdx_pgdiff$n <- n
})

# 
observeEvent(input$generate_contrast,{
  DDD.data$contrast_pgdiff <- NULL
  if(!('Difference of group difference' %in% input$contrast_type))
    return(NULL)
  
  x <- lapply(contrastIdx_pgdiff$n,function(i){
    x1 <- input[[NS(i, "Treatment_pgdiff")]]
    x2 <- input[[NS(i, "Control_pgdiff")]]
    if(length(x1)!=2 | length(x2)!=2) {
      NA
    } else {
      x1 <- paste0('(',paste0(x1,collapse = '-'),')')
      x2 <- paste0('(',paste0(x2,collapse = '-'),')')
      data.frame(Treatment=x1, Control=x2)
    }
  })
  x <- x[!is.na(x)]
  if(length(x)==0)
    return(NULL)
  x <- data.frame(do.call(rbind,x))
  colnames(x) <- c('Treatment','Control')
  x <- x[!duplicated(x),]
  x <- na.omit(x)
  idx <- which(as.vector(t(x[,1]))==as.vector(t(x[,2])))
  if(length(idx)>0)
    x <- x[-idx,,drop=F]
  if(nrow(x)==0)
    return(NULL)
  DDD.data$contrast_pgdiff <- unique(paste0(x$Treatment,'-',x$Control))
})

##--------- + Refresh contrast group----
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

observe({
  if(('Difference of group mean' %in% input$contrast_type))
    return(NULL)
  if(params_mean$btn==1)
    return(NULL)
  for(i in params_mean$btn:2){
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#param_mean', i)
    )
    params_mean$btn <- params_mean$btn - 1
    n <- contrastIdx_mean$n
    n <- n[-length(n)]
    contrastIdx_mean$n <- n
  }
})

observeEvent(input$refresh_contrast,{
  if(params_mean$btn==1)
    return(NULL)
  for(i in params_mean$btn:2){
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#param_mean', i)
    )
    params_mean$btn <- params_mean$btn - 1
    n <- contrastIdx_mean$n
    n <- n[-length(n)]
    contrastIdx_mean$n <- n
  }
})


observe({
  if(('Difference of group difference' %in% input$contrast_type))
    return(NULL)
  if(params_pgdiff$btn==1)
    return(NULL)
  for(i in params_pgdiff$btn:2){
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#param_pgdiff', i)
    )
    params_pgdiff$btn <- params_pgdiff$btn - 1
    n <- contrastIdx_pgdiff$n
    n <- n[-length(n)]
    contrastIdx_pgdiff$n <- n
  }
})

observeEvent(input$refresh_contrast,{
  if(params_pgdiff$btn==1)
    return(NULL)
  for(i in params_pgdiff$btn:2){
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#param_pgdiff', i)
    )
    params_pgdiff$btn <- params_pgdiff$btn - 1
    n <- contrastIdx_pgdiff$n
    n <- n[-length(n)]
    contrastIdx_pgdiff$n <- n
  }
})

observeEvent(input$generate_contrast,{
  DDD.data$contrast <- unique(c(DDD.data$contrast_pw,DDD.data$contrast_mean,DDD.data$contrast_pgdiff))
  message(paste0('Contrast groups: ',paste0(DDD.data$contrast,collapse = '; ')))
})

## show contrast table
output$view_contrast_table <- DT::renderDataTable({
  if(is.null(DDD.data$contrast))
    return(NULL)
  x <- DDD.data$contrast
  x <- data.frame(Contrast=paste0('Contrast group',1:length(x)),
                  `Treatment-Control`= x,check.names = F)
  datatable(x,options=list(searching=F,info=F))
})

##--------------->>Step 3: 3D analysis<<-----------------
##update parameters
observe({
  if(is.null(DDD.data$params_list$DE_pipeline))
    return(NULL)
  updateSelectInput(session = session,inputId = 'DE_pipeline',selected = DDD.data$params_list$DE_pipeline)
})

observe({
  if(is.null(DDD.data$params_list$pval_adj_method))
    return(NULL)
  updateSelectInput(session = session,inputId = 'pval_adj_method',selected = DDD.data$params_list$pval_adj_method)
})

observe({
  if(is.null(DDD.data$params_list$DAS_pval_method))
    return(NULL)
  updateSelectInput(session = session,inputId = 'DAS_pval_method',selected = DDD.data$params_list$DAS_pval_method)
})

observe({
  if(is.null(DDD.data$params_list$pval_cut))
    return(NULL)
  updateSliderInput(session = session,
                    inputId = 'pval_cut',value = DDD.data$params_list$pval_cut)
})

observe({
  if(is.null(DDD.data$params_list$l2fc_cut))
    return(NULL)
  updateSliderInput(session = session,
                    inputId = 'l2fc_cut',value = DDD.data$params_list$l2fc_cut)
})

observe({
  if(is.null(DDD.data$params_list$deltaPS_cut))
    return(NULL)
  updateSliderInput(session = session,
                    inputId = 'deltaPS_cut',value = DDD.data$params_list$deltaPS_cut)
})


observeEvent(input$run.DE,{
  if(is.null(DDD.data$contrast)){
    showmessage('Please generate contrast groups.')
    return(NULL)
  }
   
  cat('\nDE gene analysis\n')
  ##---------------+ DE genes -----------------
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
                 DE_pipeline <- input$DE_pipeline
                 incProgress(0.2)
                 switch(DE_pipeline,
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
  
  ##---------------+ Generating deltaPS -----------------
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
  
  ##---------------+ DAS genes, DE/DTU transcripts -----------------
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
                 DE_pipeline <- input$DE_pipeline
                 incProgress(0.2)
                 switch(DE_pipeline,
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

##---------------->>  Step 4: 3D testing statistics   <<-----------------
##----------------+ DE DAS DTU statistics -----------------

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
  
  stat <- x[[input$contrast.groups]][1:min(nrow(x[[input$contrast.groups]]),input$top.stat),]
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

##---------------->>  Step 5: Significant 3D numbers   <<-----------------
##---------------+ Transcripts per gene number -----------
tpg_g <- eventReactive(input$tpg_plot_button,{
  if(is.null(DDD.data$mapping) | is.null(DDD.data$target_high))
    return(NULL)
  x <- DDD.data$mapping[DDD.data$target_high$trans_high,]$GENEID
  y <- table(table(x))
  z <- c(y[1:10],sum(y[11:20]),sum(y[21:length(y)]))
  data2plot <- data.frame(trans=c(1:10,'11-20',paste0('21-',length(y))),genes=z)
  data2plot$trans <- factor(data2plot$trans,levels = data2plot$trans)
  g <- ggplot(data = data2plot,aes(x=trans,y=genes,label=genes))+
    geom_bar(stat='identity',fill=input$tpg_color)+
    labs(x='Number of transcript per gene',y='Number of genes',
         title='Distribution of the number of transcripts per gene')+
    geom_text(size = 3, position=position_dodge(width=0.9), vjust=-0.25)+
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
  print(tpg_g())
})

##save plot
observeEvent(input$tpg_save_button,{
  if(is.null(tpg_g()))
    return(NULL)
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
  showmessage()
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

##---------------+ Euler plot between contrast groups-----------------
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

observeEvent(input$save_3D_euler_plot,{
  ###save DE.genes
  lapply(names(g.across.contrast()),function(i){
    figure.name <- paste0(DDD.data$figure.folder,'/',i,' euler plot across contrast ',
                          paste0(gsub('/','over',input$across.contrast.group),collapse = ' vs '))
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

#*****************preview saved plot******************
observeEvent(input$save_3D_euler_plot_view, {
  title.idx <- c('DE genes','DAS genes','DE transcripts','DTU transcripts')
  file2save <- lapply(title.idx,function(i){
    paste0(DDD.data$figure.folder,'/',i,' euler plot across contrast ',
           paste0(gsub('/','over',input$across.contrast.group),collapse = ' vs '),'.png')
  })
  if(!all(sapply(file2save,file.exists))){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    file2save.div <- lapply(file2save,function(i){
      div(tags$img(src=base64enc::dataURI(file = i),
                   style="height: 100%; width: 100%; display: block; margin-left: auto;
                   margin-right: auto;border:1px solid #000000"),p())
      })
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      file2save.div,
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='m'
    ))
  }
})


###target numbers
observe({
  if(is.null(DDD.data$DE_genes) | 
     is.null(DDD.data$DAS_genes) | 
     is.null(DDD.data$DE_trans)| 
     is.null(DDD.data$DTU_trans))
    return(NULL)
  DDD.data$DDD_numbers <- summary3Dnumber(DE_genes = DDD.data$DE_genes,
                                          DAS_genes = DDD.data$DAS_genes,
                                          DE_trans = DDD.data$DE_trans,
                                          DTU_trans = DDD.data$DTU_trans,
                                          contrast = DDD.data$contrast)
  
})

output$DDD.numbers <- renderTable({
  DDD.data$DDD_numbers
},align='c')

##---------------+ DDD up- and down-regulation-----------------
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

output$up.down.bar.plot.DE.trans <- renderPlot({
  if(is.null((g.updown())))
    return(NULL)
  print(g.updown()[[2]])
})

observeEvent(input$save_up_down_bar_plot,{
  if(is.null(g.updown()))
    return(NULL)
  withProgress(message = 'Saving plots ...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0)
                 lapply(names(g.updown()),function(i){
                   png(paste0(DDD.data$figure.folder,'/',i,' up and down regulation numbers.png'),
                       width = input$updown.plot.width,height = input$updown.plot.height,units = 'in',
                       res = input$updown.plot.res)
                   print(g.updown()[[i]])
                   dev.off()
                   
                   pdf(paste0(DDD.data$figure.folder,'/',i,' up and down regulation numbers.pdf'),
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

##*****************preview saved plot******************
observeEvent(input$save_up_down_bar_plot_view, {
  file2save1 <- paste0(DDD.data$figure.folder,'/DE genes up and down regulation numbers.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/DE transcripts up and down regulation numbers.png')
  
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
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
               border:1px solid #000000"),
      p(),
      tags$img(src=base64enc::dataURI(file = file2save2),
               style="height: 100%; width: 100%; display: block; margin-left: auto; margin-right: auto;
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
      ))
  }
})

##---------------+ Volcano plot-----------------
g.volcano <- eventReactive(input$Volcano_plot_button,{
  if(is.null((DDD.data$DE_genes)) | 
     is.null((DDD.data$DAS_genes)) | 
     is.null((DDD.data$DE_trans)) | 
     is.null((DDD.data$DTU_trans)))
    return(NULL)
  top.n <- input$Volcano_top_label
  size <- input$Volcano_point_size
  col0 <- input$color4Volcano0
  col1 <- input$color4Volcano1
  
  idx <- c('DE genes','DAS genes','DE transcripts','DTU transcripts')
  title.idx <- c(paste0("Volcano plot: DE genes (Low expression filtered; \nAdjusted p<",
                        input$pval_cut,'; |L2FC|>=',input$l2fc_cut,'; Labels: top ',top.n,' distance to (0,0))'),
                 paste0("Volcano plot: DAS genes (Low expression filtered; \nAdjusted p<",
                        input$pval_cut,'; |MaxdeltaPS|>=',input$deltaPS_cut,'; Labels: top ',top.n,' distance to (0,0))'),
                 paste0("Volcano plot: DE transcripts (Low expression filtered; \nAdjusted p<",
                        input$pval_cut,'; |L2FC|>=',input$l2fc_cut,'; Labels: top ',top.n,' distance to (0,0))'),
                 paste0("Volcano plot: DTU transcripts (Low expression filtered; \nAdjusted p<",
                        input$pval_cut,'; |deltaPS|>=',input$deltaPS_cut,'; Labels: top ',top.n,' distance to (0,0))')
                 )
  names(title.idx) <- idx
  g <- lapply(idx,function(i){
    if(i=='DE genes'){
      DDD.stat <- DDD.data$genes_3D_stat$DE.stat
      data2plot <- data.frame(target=DDD.stat$target,
                              contrast=DDD.stat$contrast,
                              x=DDD.stat$log2FC,
                              y=-log10(DDD.stat$adj.pval))
      data2plot$significance <- 'Not significant'
      data2plot$significance[DDD.stat$adj.pval< input$pval_cut & abs(DDD.stat$log2FC)>=input$l2fc_cut] <- 'Significant'
      q <- plotVolcano(data2plot = data2plot,
                       xlab = 'log2FC of genes',
                       ylab='-log10(FDR)',title = title.idx[i],
                       col0 = col0,col1 = col1,size = size,top.n = top.n)
    }
    
    if(i=='DAS genes'){
      if(input$DAS_pval_method=='F-test')
        DDD.stat <- DDD.data$trans_3D_stat$DAS.F.stat else DDD.stat <- DDD.data$trans_3D_stat$DAS.simes.stat
      data2plot <- data.frame(target=DDD.stat$target,
                              contrast=DDD.stat$contrast,
                              x=DDD.stat$maxdeltaPS,
                              y=-log10(DDD.stat$adj.pval))
      data2plot$significance <- 'Not significant'
      data2plot$significance[DDD.stat$adj.pval< input$pval_cut & abs(DDD.stat$maxdeltaPS)>=input$deltaPS_cut] <- 'Significant'
      q <- plotVolcano(data2plot = data2plot,xlab = 'MaxdeltaPS of transcripts in each gene',ylab='-log10(FDR)',
                       title = title.idx[i],col0 = col0,col1 = col1,size = size,top.n = top.n)
    }
    
    if(i=='DE transcripts'){
      DDD.stat <- DDD.data$trans_3D_stat$DE.stat
      data2plot <- data.frame(target=DDD.stat$target,
                              contrast=DDD.stat$contrast,
                              x=DDD.stat$log2FC,
                              y=-log10(DDD.stat$adj.pval))
      data2plot$significance <- 'Not significant'
      data2plot$significance[DDD.stat$adj.pval< input$pval_cut & abs(DDD.stat$log2FC)>=input$l2fc_cut] <- 'Significant'
      q <- plotVolcano(data2plot = data2plot,xlab = 'log2FC of transcripts',ylab='-log10(FDR)',title = title.idx[i],
                       col0 = col0,col1 = col1,size = size,top.n = top.n)
    }
    
    if(i=='DTU transcripts'){
      DDD.stat <- DDD.data$trans_3D_stat$DTU.stat 
      data2plot <- data.frame(target=DDD.stat$target,
                              contrast=DDD.stat$contrast,
                              x=DDD.stat$deltaPS,
                              y=-log10(DDD.stat$adj.pval))
      data2plot$significance <- 'Not significant'
      data2plot$significance[DDD.stat$adj.pval< input$pval_cut & abs(DDD.stat$deltaPS)>=input$deltaPS_cut] <- 'Significant'
      q <- plotVolcano(data2plot = data2plot,xlab = 'deltaPS of transcripts',ylab='-log10(FDR)',title = title.idx[i],
                       col0 = col0,col1 = col1,size = size,top.n = top.n)
    }
    q
  })
  names(g) <- idx
  g
})

##make plots
output$DE_genes_Volcano <- renderPlot({
  if(is.null((g.volcano())))
    return(NULL)
  showNotification('Plotting volcano plot, please wait ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'volcano_de_genes_message')
  print(g.volcano()[['DE genes']])
  removeNotification(id = 'volcano_de_genes_message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})

output$DAS_genes_Volcano <- renderPlot({
  if(is.null((g.volcano())))
    return(NULL)
  showNotification('Plotting volcano plot, please wait ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'volcano_das_genes_message')
  print(g.volcano()[['DAS genes']])
  removeNotification(id = 'volcano_das_genes_message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})

output$DE_trans_Volcano <- renderPlot({
  if(is.null((g.volcano())))
    return(NULL)
  showNotification('Plotting volcano plot, please wait ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'volcano_de_trans_message')
  print(g.volcano()[['DE transcripts']])
  removeNotification(id = 'volcano_de_trans_message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})

output$DTU_trans_Volcano <- renderPlot({
  if(is.null((g.volcano())))
    return(NULL)
  showNotification('Plotting volcano plot, please wait ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'volcano_dtu_trans_message')
  print(g.volcano()[['DTU transcripts']])
  removeNotification(id = 'volcano_dtu_trans_message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})

##save plots
observeEvent(input$save_Volcano_plot,{
  if(is.null(g.volcano()))
    return(NULL)
  showNotification('Saving volcano plot, please wait ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'save_volcano__message')
  withProgress(message = 'Saving plots ...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0)
                 
                 lapply(names(g.volcano()),function(i){
                   incProgress(1/4)
                   png(paste0(DDD.data$figure.folder,'/',i,' volcano plot.png'),
                       width = input$Volcano_width,height = input$Volcano_height,units = 'in',
                       res = input$Volcano_res)
                   print(g.volcano()[[i]])
                   dev.off()
                   
                   pdf(paste0(DDD.data$figure.folder,'/',i,' volcano plot.pdf'),
                       width = input$Volcano_width,height = input$Volcano_height)
                   print(g.volcano()[[i]])
                   dev.off()
                 })
                 message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
                 showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
                 incProgress(1)
               })
  removeNotification(id = 'save_volcano__message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})

##*****************preview saved plot******************
observeEvent(input$save_Volcano_plot_view, {
  title.idx <- c('DE genes','DAS genes','DE transcripts','DTU transcripts')
  file2save <- lapply(title.idx,function(i){
    paste0(DDD.data$figure.folder,'/',i,' volcano plot.png')
  })
  if(!all(sapply(file2save,file.exists))){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    file2save.div <- lapply(file2save,function(i){
      div(tags$img(src=base64enc::dataURI(file = i),
                   style="height: 100%; width: 100%; display: block; margin-left: auto;
                   margin-right: auto;border:1px solid #000000"),p())
    })
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      file2save.div,
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  }
})

##---------------->>  Step 6: Transcriptional vs alternative splicing regulation   <<-----------------

##---------------+ flow chat-----------------
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

observeEvent(input$save_union_flow_chart,{
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

##*****************preview saved plot******************
observeEvent(input$save_union_flow_chart_view, {
  file2save1 <- paste0(DDD.data$figure.folder,'/Union set DE genes vs DAS genes.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/Union set DE transcripts vs DTU transcripts.png')
  
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



##---------------+ DE vs DAS-----------------
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


##---------------+ DE vs DTU-----------------

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
  g <- plotEulerDiagram(x = x,fill = c(input$color4_dedas_number_from,input$color4_dedas_number_to),
                        scale.plot = input$across.target.euler.plot.scale)
  g
})

output$DEvsDTU.euler <- renderPlot({
  if(is.null(g.across.target2()))
    return(NULL)
  grid.arrange(g.across.target2(),top=textGrob(paste0('DE vs DTU transcripts\nin contrast group: ',input$across.target), 
                                               gp=gpar(cex=1.2)))
})

observeEvent(input$save_transvsAS_euler_plot,{
  
  # g.across.target1 <- reactive({
  #   if(is.null(DDD.data$DEvsDAS_results)| is.null(input$across.target))
  #     return(NULL)
  #   
  # })
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
  
  # ##DE vs DAS genes
  # figure.name <- paste0(DDD.data$figure.folder,'/DE vs DAS gene euler plot in contrast ',input$across.target)
  # png(paste0(figure.name,'.png'),
  #     width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height,
  #     res=input$across.target.euler.plot.res, units = 'in')
  # grid.arrange(g.across.target1(),top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',input$across.target), 
  #                                              gp=gpar(cex=1.2)))
  # dev.off()
  # pdf(paste0(figure.name,'.pdf'),
  #     width = input$across.target.euler.plot.width,height = input$across.target.euler.plot.height)
  # grid.arrange(g.across.target1(),top=textGrob(paste0('DE vs DAS genes\nin contrast group: ',input$across.target), 
  #                                              gp=gpar(cex=1.2)))
  # dev.off()
  # 
  
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})

#*****************preview saved plot******************
observeEvent(input$save_transvsAS_euler_plot_view, {
  
  file2save1 <- lapply(DDD.data$contrast,function(i){
    paste0(DDD.data$figure.folder,'/DE vs DAS gene euler plot in contrast ', gsub('/','over',i),'.png')
  })
  file2save2 <- lapply(DDD.data$contrast,function(i){
    paste0(DDD.data$figure.folder,'/DE vs DTU transcript euler plot in contrast ', gsub('/','over',i),'.png')
  })
  if(!all(sapply(file2save1,file.exists)) | !all(sapply(file2save2,file.exists)) | is.null(DDD.data$contrast)){
    showModal(modalDialog(
      h4('Figures are not found in the directory. Please double check if the figures are saved.'),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  } else {
    file2save.div <- lapply(1:length(DDD.data$contrast),function(i){
      div(
      div(tags$img(src=base64enc::dataURI(file = file2save1[[i]]),
                   style="width: 45%; display: inline-block; margin-left: auto;
                   margin-right: auto;border:1px solid #000000; float: auto;"),
          tags$img(src=base64enc::dataURI(file = file2save2[[i]]),
                   style="width: 45%; display: inline-block; margin-left: auto;
                   margin-right: auto;border:1px solid #000000; float: auto")),
      p())
    })
    showModal(modalDialog(
      h4('Saved figure. If the labels are not properly placed, please correct width, height, ect.'),
      file2save.div,
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
    ))
  }
})

##


observeEvent(input$page_before_ddd, {
  newtab <- switch(input$tabs, "preprocessing" = "ddd","ddd" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_ddd, {
  newtab <- switch(input$tabs, "function" = "ddd","ddd" = "function")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
