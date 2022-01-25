#--------- >> Step 1: Time-series experimental design----
observeEvent(input$tabs,{
  if(input$tabs=='tstrend'){
    text2show <- 'Page: Time-series trend'
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
    showmessage(text = text2show,showNoteify = F)
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
  }
})

observe({
  if(is.null(DDD.data$params_list$whether_tstrend))
    return(NULL)
  updateRadioButtons(session = session,inputId = 'whether_tstrend',selected = DDD.data$params_list$whether_tstrend)
})

observe({
  if(is.null(DDD.data$samples_new))
    return(NULL)
  if(is.null(DDD.data$params_list$tstrend_time_group)){
    updateSelectInput(session = session,inputId = 'tstrend_time_group',
                      choices = colnames(DDD.data$samples_new),selected = '')
  } else {
    updateSelectInput(session = session,inputId = 'tstrend_time_group',
                      choices = colnames(DDD.data$samples_new),selected = DDD.data$params_list$tstrend_time_group)
  }
  
  if(is.null(DDD.data$params_list$tstrend_time)){
    updateSelectInput(session = session,inputId = 'tstrend_time',
                      choices = colnames(DDD.data$samples_new),selected = '')
  } else {
    updateSelectInput(session = session,inputId = 'tstrend_time',
                      choices = colnames(DDD.data$samples_new),
                      selected = DDD.data$params_list$tstrend_time)
  }
})

observe({
  if(is.null(DDD.data$tstrend_spline))
    return(NULL)
  updateRadioButtons(session = session,inputId = 'tstrend_spline',selected = DDD.data$tstrend_spline)
})

observe({
  if(is.null(DDD.data$params_list$tstrend_time))
    return(NULL)
  updateSelectInput(session = session,inputId = 'tstrend_time',selected = DDD.data$params_list$tstrend_time)
})

observe({
  if(is.null(DDD.data$params_list$tstrend_time_group))
    return(NULL)
  updateSelectInput(session = session,inputId = 'tstrend_time_group',selected = DDD.data$params_list$tstrend_time_group)
})



observe({
  if(!is.null(DDD.data$tstrend_spline_df)){
    updateSliderInput(session = session,inputId = 'tstrend_spline_df',value = DDD.data$tstrend_spline_df)
  } else {
    if(is.null(DDD.data$tstrend_timePoints))
      return(NULL)
    updateSliderInput(session = session,inputId = 'tstrend_spline_df',
                      max = length(unique(DDD.data$tstrend_timePoints))-1)
  }
})

observe({
  if(is.null(DDD.data$params_list$mean_variance_tstrend_method))
    return(NULL)
  updateSelectInput(session = session,inputId = 'mean_variance_tstrend_method',
                    selected = DDD.data$params_list$mean_variance_tstrend_method)
})

observe({
  if(is.null(DDD.data$params_list$pval_adj_method_tstrend))
    return(NULL)
  updateSelectInput(session = session,inputId = 'pval_adj_method_tstrend',
                    selected = DDD.data$params_list$pval_adj_method_tstrend)
})

observe({
  if(is.null(DDD.data$params_list$DAS_pval_method_tstrend))
    return(NULL)
  updateSelectInput(session = session,inputId = 'DAS_pval_method_tstrend',
                    selected = DDD.data$params_list$DAS_pval_method_tstrend)
})

observe({
  if(is.null(DDD.data$params_list$DAS_pval_method_across_group))
    return(NULL)
  updateSelectInput(session = session,inputId = 'DAS_pval_method_across_group',
                    selected = DDD.data$params_list$DAS_pval_method_across_group)
})


observe({
  if(is.null(DDD.data$params_list$pval_cut_tstrend))
    return(NULL)
  updateSliderInput(session = session,inputId = 'pval_cut_tstrend',value = DDD.data$params_list$pval_cut_tstrend)
})


observeEvent(input$tstrend_time,{
  if(is.null(DDD.data$samples_new) | is.null(input$tstrend_time) | input$tstrend_time=='')
    return(NULL)
  x <- DDD.data$samples_new[,input$tstrend_time]
  x <- factor(x,levels = unique(x))
  DDD.data$tstrend_timePoints <- x
},ignoreInit = T)

# observeEvent(input$tstrend_time,{
#   if(!is.null(DDD.data$tstrend_timePoints) & length(unique(DDD.data$tstrend_timePoints))<=2){
#     updateRadioButtons(session = session,inputId = 'tstrend_spline',selected = 'No')
#     showmessage('Time-points are too few to use spline method')
#   }
# })

observeEvent(input$tstrend_time_group,{
  if(is.null(DDD.data$samples_new) | is.null(input$tstrend_time_group) | input$tstrend_time_group=='')
    return(NULL)
  x <- DDD.data$samples_new[,input$tstrend_time_group]
  x <- factor(x,levels = unique(x))
  DDD.data$tstrend_timeGroups <- x
},ignoreInit = T)

# observeEvent(input$tstrend_time_group,{
#   if(!is.null(DDD.data$tstrend_timePoints) & length(unique(DDD.data$tstrend_timePoints))<=2){
#     updateRadioButtons(session = session,inputId = 'tstrend_spline',selected = 'No')
#     showmessage('Time-points are too few to use spline method')
#   }
# })


output$tstrend_design_table <- DT::renderDataTable({
  if(is.null(DDD.data$tstrend_timePoints) | is.null(DDD.data$tstrend_timeGroups))
    return(NULL)
  x <- data.frame(Group=DDD.data$tstrend_timeGroups,Time=DDD.data$tstrend_timePoints)
  if(input$tstrend_spline=='Yes')
    x$`Numeric time-coordinates for Spline` <- as.numeric(factor(DDD.data$tstrend_timePoints))
  x
})


output$tstrend_example_plot <- renderPlot({
  if(is.null(DDD.data$tstrend_timePoints) | is.null(DDD.data$tstrend_timeGroups))
    return(NULL)
  # Time0 <- Time.points()
  # Time <- as.numeric(factor(Time0))
  Time.points <- DDD.data$tstrend_timePoints
  Time.groups <- DDD.data$tstrend_timeGroups
  Time <- as.numeric(factor(Time.points))
  set.seed(100)
  Expression <- sin(Time*pi/6)+rnorm(length(Time),1,0.5)
  Expression <- Expression*as.numeric(factor(Time.groups))
  
  if(input$tstrend_spline == 'Yes'){
    Time.ns <- splines::ns(Time,df=input$tstrend_spline_df)
    knots <- sort(c(attr(Time.ns,'knots'),attr(Time.ns,'Boundary.knots')),decreasing = F)
    data2fit <- data.frame(Time=Time,Expression=Expression,Group=Time.groups)
    data2line <- by(data = data2fit,INDICES = data2fit$Group,function(x){
      fit <- tryCatch({lm(Expression ~ splines::ns(Time,df=input$tstrend_spline_df),data = x)},
                      error=function(e){})
      if(is.null(fit))
        return(NULL)
      pred.value <- predict(fit,data.frame(Time=knots))
      names(pred.value) <- knots
      pred.value
    })
    data2line <- do.call(rbind,data2line)
    if(is.null(data2line)){
      showmessage(text = 'Please select correct columns of time-points and time-points groups.')
      return(NULL)
    }
    
    data2line <- reshape2::melt(data2line)
    colnames(data2line) <- c('Group','Time','Expression')
    
    data2plot <- data.frame(Time=Time,Label=Time.points,Group=Time.groups,Expression=Expression)
    # data2line <- data.frame(Time=knots,Expression=pred.value)
    g <- ggplot(data = data2plot,aes(x=Time,y=Expression))+
      geom_point(aes(pch='Samples reps'),colour='black')+
      scale_shape_manual(values = c('Samples reps'=1))+
      theme_bw()+
      labs(title='Spline time-series with Natural Cubic method',shape='Data points')
    q <- g+geom_line(data = data2line,aes(x=Time,y=Expression,colour='red'))+
      geom_point(data = data2line,pch=16,size=3,aes(colour='blue'))+facet_grid(.~Group)+
      scale_color_identity(name = "Spline fit",
                           breaks = c( "red", "blue"),
                           labels = c("Spline fitted line", "Spline knots"),
                           guide = "legend")
  } else {
    data2plot <- data.frame(Time=Time.points,Expression=Expression,Group=Time.groups)
    q <- ggplot(data2plot, aes(x=Time, y=Expression)) + 
      geom_point(aes(pch='Samples reps'),colour='black') +
      stat_summary(aes(y = Expression,group=1,colour="Mean value"), fun.y=mean, geom="line",group=1)+
      facet_grid(.~Group)+
      theme_bw()+
      labs(title='Time-series trend',shape='Data points',colour='Fitted line')+
      scale_shape_manual(values = c('Samples reps'=1))
  }
  print(q)
})

##--------- >> Step 2: Compare trend across time groups----
observe({
  if(is.null(DDD.data$tstrend_timeGroups))
    return(NULL)
  Time.groups <- DDD.data$tstrend_timeGroups
  if(!is.factor(Time.groups))
    Time.groups <- factor(Time.groups,levels = unique(Time.groups))
  thisList <- levels(Time.groups)
  
  callModule(module = selectorServer_tstrend, id = 1, 
             thisList = thisList,selected = thisList)
  updateSelectInput(session = session,inputId = NS(1)('tstrend_groups'),
                    choices = thisList,selected = thisList)
})

params_tstrend <- reactiveValues(btn = 1)
contrastIdx_tstrend <- reactiveValues(n=1)

##----->insert lines
observeEvent(input$insertParamBtn_tstrend, {
  if(is.null(DDD.data$tstrend_timeGroups))
    return(NULL)
  Time.groups <- DDD.data$tstrend_timeGroups
  if(!is.factor(Time.groups))
    Time.groups <- factor(Time.groups,levels = unique(Time.groups))
  thisList <- levels(Time.groups)
  
  n <- contrastIdx_tstrend$n
  n <- c(n,length(n)+1)
  contrastIdx_tstrend$n <- n
  params_tstrend$btn <- params_tstrend$btn + 1
  callModule(module = selectorServer_tstrend, id = params_tstrend$btn,
             thisList = thisList)
  insertUI(
    selector = '#placeholder_tstrend',
    ui = selectorUI_tstrend(params_tstrend$btn)
  )
})

##----->remove lines
observeEvent(input$removeParamBtn_tstrend, {
  removeUI(
    ## pass in appropriate div id
    selector = paste0('#param_tstrend', params_tstrend$btn)
  )
  params_tstrend$btn <- params_tstrend$btn - 1
  n <- contrastIdx_tstrend$n
  n <- n[-length(n)]
  contrastIdx_tstrend$n <- n
})

observeEvent(input$generate_tstrend_groups,{
  startmessage(text = 'Generate contrast groups for time-series trend analysis')
  x <- lapply(contrastIdx_tstrend$n,function(i){
   input[[NS(i, "tstrend_groups")]]
  })
  n.x <- length(x)
  x <- x[!unlist(sapply(x,is.na))]
  x <- x[!sapply(x,is.null)]
  x <- x[!duplicated(x)]
  x <- x[sapply(x,length)>1]
  
  if(length(x)!=n.x)
    showmessage(text = 'Incorrect contrast groups are removed.')
  
  if(length(x)==0)
    return(NULL)
  
  # message(x)
  # save(x,file='x.RData')
  names(x) <- sapply(x, function(i) paste0(i,collapse = '='))
  DDD.data$tstrend_Contrastgroups <- x
  endmessage(text = 'Generate contrast groups for time-series trend analysis')
},ignoreInit = T )


# observeEvent(input$generate_tstrend_groups,{
#   if(!is.null(DDD.data$tstrend_timePoints) & length(unique(DDD.data$tstrend_timePoints))<=2){
#     updateRadioButtons(session = session,inputId = 'tstrend_spline',selected = 'No')
#     showmessage('Time-points are too few to use spline method')
#   }
# })

output$tstrend_group_table <- DT::renderDataTable({
  if(is.null(DDD.data$tstrend_Contrastgroups))
    return(NULL)
  x <- sapply(DDD.data$tstrend_Contrastgroups,function(i){
    paste(i,collapse = ' = ')
  })
  x <- paste0('Time-series trend: ',x)
  x <- data.frame(`Time-series groups`=paste0('Time-series group',1:length(x)),
                  `Null hypothesis`= x,`Alternative hypothesis`='Not all equal',check.names = F)
  DDD.data$DDD_tstrend_design <- x
  datatable(x,options=list(searching=F,info=F))
})


output$tstrend_test_number_table <- renderTable({
  if(is.null(DDD.data$DE_tstrend_genes) | 
     is.null(DDD.data$DAS_tstrend_genes) | 
     is.null(DDD.data$DE_tstrend_trans)| 
     is.null(DDD.data$DTU_tstrend_trans))
    return(NULL)
  idx <- c('DE_tstrend_genes','DAS_tstrend_genes','DE_tstrend_trans','DTU_tstrend_trans')
  o.idx <- names(DDD.data$tstrend_Contrastgroups)
  x0 <- rep(0,length(o.idx))
  names(x0) <- o.idx
  
  x <- lapply(idx, function(i){
    split.idx <- factor(DDD.data[[i]]$contrast,levels = o.idx)
    x <- sapply(split(DDD.data[[i]]$targets,split.idx),length)
    x[o.idx]
  })
  x <- t(do.call(rbind,x))
  x <- data.frame(contrast=rownames(x),x,row.names = NULL)
  colnames(x) <- c('Contrast','DE TS trendgenes','DAS TS trend genes',
                   'DE TS trend transcripts','DTU TS trend transcripts')
  DDD.data$DDD_tstrend_numbers <- x
  x
})


##--------- >> Step 3: Time-series 3D trend analysis----
observeEvent(input$run_3d_tstrend,{
  # if(length(unique(DDD.data$tstrend_timePoints)) < 2 & input$tstrend_spline=='Yes'){
  #   showmessage('Time-points are too few to use spline method. Please select "No" for spline method in Step 1.')
  #   return(NULL)
  # }
  
  if(is.null(DDD.data$tstrend_timePoints) | 
     is.null(DDD.data$tstrend_timeGroups) )
    return(NULL)
  if(is.null(DDD.data$tstrend_Contrastgroups)){
    showmessage('Please generate time groups for comparisons.')
    return(NULL)
  }
  
  ###====>data prepare
  time.points <- DDD.data$tstrend_timePoints
  time.groups <- DDD.data$tstrend_timeGroups
  if(!is.factor(time.points))
    time.points <- factor(time.points,levels = unique(time.points))
  if(!is.factor(time.groups))
    time.groups <- factor(time.groups,levels = unique(time.groups))


  if(input$tstrend_spline == 'Yes'){
    time.points <- as.numeric(factor(time.points,levels = unique(time.points)))
    Time <- splines::ns(time.points, df=pmax(2,df=input$tstrend_spline_df))
    colnames(Time) <- paste0('ns',1:ncol(Time))
    DDD.data$tstrend_spline_df <- pmax(2,df=input$tstrend_spline_df)
  } else {
    Time <- time.points
  }
  
  Group <- time.groups
  
  ######Gene level analysis
  startmessage(text = 'Gene level 3D trend analysis')
  withProgress(message = 'DE gene time-series trend analysis...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0.1)
                 batch.effect <- DDD.data$genes_batch$W
                 dge <- DDD.data$genes_dge
                 if(is.null(batch.effect)){
                   design <- model.matrix(~Group*Time)
                 } else {
                   design <- model.matrix(~Group*Time+batch.effect)
                 }
                 
                 if(input$tstrend_spline == 'Yes' && ncol(Time)==1){
                   design.col <- gsub('Time$','Timens1',colnames(design))
                   colnames(design) <- design.col
                 }
                 
                 colnames(design) <- gsub(pattern = 'Group|Time','',colnames(design))
                 colnames(design) <- gsub(pattern = '[:]','.',colnames(design))
                 colnames(design)[colnames(design)=='(Intercept)'] <- 'Intercept'
                 # DDD.data$tstrend_Contrastgroups <- list(c('Day0','Day1','Day4'),c('Day0','Day4'),c('Day1','Day4'),c('Day0','Day1')) 
                 if(input$tstrend_spline == 'Yes'){
                   contrast.list <- lapply(DDD.data$tstrend_Contrastgroups,function(i){
                     TSgroup2contrast(Time = colnames(Time),Group = i,Intercept = list(Time=NULL,Group=levels(Group)[1]))
                   }) 
                 } else {
                   contrast.list <- lapply(DDD.data$tstrend_Contrastgroups,function(i){
                     TSgroup2contrast(Time = levels(Time),Group = i,Intercept = list(Time=levels(Time)[1],Group=levels(Group)[1]))
                   })
                 }
                 names(contrast.list) <- names(DDD.data$tstrend_Contrastgroups)
                 incProgress(0.5)
                 ##trend analysis
                 genes_3D_tstrend_stat <- TStrend.pipeline(dge = dge,
                                                           design = design,
                                                           span = 0.5,
                                                           contrast.list = contrast.list,
                                                           adjust.method = input$pval_adj_method_tstrend,
                                                           diffAS = F,
                                                           aggPvalue.method = input$DAS_pval_method_across_group,
                                                           voomWeights = input$mean_variance_tstrend_method=='voomWeights')
                 DDD.data$genes_3D_tstrend_stat <- genes_3D_tstrend_stat
                 
                 incProgress(0.8)
                 ##DE genes
                 stat <- genes_3D_tstrend_stat$DE.pval
                 DE_tstrend_genes <- lapply(colnames(stat),function(i){
                   sub.stat <- stat[,i]
                   x <- data.frame(targets=rownames(stat),contrast=i,adj.pval=sub.stat)
                   x <- x[x$adj.pval < input$pval_cut_tstrend,]
                   x
                 })
                 DE_tstrend_genes <- do.call(rbind,DE_tstrend_genes)
                 rownames(DE_tstrend_genes) <- NULL
                 DDD.data$DE_tstrend_genes <- DE_tstrend_genes
                 incProgress(1)
               })
  endmessage(text = 'Gene level 3D trend analysis')
  
  startmessage(text = 'AS level 3D trend analysis')
  ######Transcript level analysis
  withProgress(message = 'DE transcript time-series trend analysis...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0.1)
                 batch.effect <- DDD.data$trans_batch$W
                 dge <- DDD.data$trans_dge
                 if(is.null(batch.effect)){
                   design <- model.matrix(~Group*Time)
                 } else {
                   design <- model.matrix(~Group*Time+batch.effect)
                 }
                 
                 if(input$tstrend_spline == 'Yes' && ncol(Time)==1){
                   design.col <- gsub('Time$','Timens1',colnames(design))
                   colnames(design) <- design.col
                 }
                 
                 colnames(design) <- gsub(pattern = 'Group|Time','',colnames(design))
                 colnames(design) <- gsub(pattern = '[:]','.',colnames(design))
                 colnames(design)[colnames(design)=='(Intercept)'] <- 'Intercept'
                 # DDD.data$tstrend_Contrastgroups <- list(c('Day0','Day1','Day4'),c('Day0','Day4'),c('Day1','Day4'),c('Day0','Day1')) 
                 if(input$tstrend_spline == 'Yes'){
                   contrast.list <- lapply(DDD.data$tstrend_Contrastgroups,function(i){
                     TSgroup2contrast(Time = colnames(Time),Group = i,Intercept = list(Time=NULL,Group=levels(Group)[1]))
                   }) 
                 } else {
                   contrast.list <- lapply(DDD.data$tstrend_Contrastgroups,function(i){
                     TSgroup2contrast(Time = levels(Time),Group = i,Intercept = list(Time=levels(Time)[1],Group=levels(Group)[1]))
                   })
                 }
                 names(contrast.list) <- names(DDD.data$tstrend_Contrastgroups)
                 incProgress(0.5)
                 ##trend analysis
                 trans_3D_tstrend_stat <- TStrend.pipeline(dge = dge,
                                                           design = design,
                                                           span = 0.5,
                                                           contrast.list = contrast.list,
                                                           adjust.method = input$pval_adj_method_tstrend,
                                                           diffAS = T,
                                                           aggPvalue.method = input$DAS_pval_method_across_group,
                                                           voomWeights = input$mean_variance_tstrend_method=='voomWeights')
                 DDD.data$trans_3D_tstrend_stat <- trans_3D_tstrend_stat
                 
                 incProgress(0.6)
                 ##DE trans
                 stat <- trans_3D_tstrend_stat$DE.pval
                 DE_tstrend_trans <- lapply(colnames(stat),function(i){
                   sub.stat <- stat[,i]
                   x <- data.frame(targets=rownames(stat),contrast=i,adj.pval=sub.stat)
                   x <- x[x$adj.pval < input$pval_cut_tstrend,]
                   x
                 })
                 DE_tstrend_trans <- do.call(rbind,DE_tstrend_trans)
                 rownames(DE_tstrend_trans) <- NULL
                 DDD.data$DE_tstrend_trans <- DE_tstrend_trans
                 
                 incProgress(0.7)
                 ##DAS genes
                 
                 # DDD.data$DAS_pval_method_tstrend <- input$DAS_pval_method_tstrend
                 
                 if(input$DAS_pval_method_tstrend=='F-test')
                   stat <- trans_3D_tstrend_stat$DAS.pval.F else stat <- trans_3D_tstrend_stat$DAS.pval.simes
                 
                 
                 DAS_tstrend_genes <- lapply(colnames(stat),function(i){
                   sub.stat <- stat[,i]
                   x <- data.frame(targets=rownames(stat),contrast=i,adj.pval=sub.stat)
                   x <- x[x$adj.pval < input$pval_cut_tstrend,]
                   x
                 })
                 DAS_tstrend_genes <- do.call(rbind,DAS_tstrend_genes)
                 rownames(DAS_tstrend_genes) <- NULL
                 DDD.data$DAS_tstrend_genes <- DAS_tstrend_genes
                 
                 incProgress(0.8)
                 ##DTU trans
                 stat <- trans_3D_tstrend_stat$DTU.pval
                 DTU_tstrend_trans <- lapply(colnames(stat),function(i){
                   sub.stat <- stat[,i]
                   x <- data.frame(targets=rownames(stat),contrast=i,adj.pval=sub.stat)
                   x <- x[x$adj.pval < input$pval_cut_tstrend,]
                   x
                 })
                 DTU_tstrend_trans <- do.call(rbind,DTU_tstrend_trans)
                 rownames(DTU_tstrend_trans) <- NULL
                 DDD.data$DTU_tstrend_trans <- DTU_tstrend_trans
                 incProgress(1)
               })

  endmessage(text = 'AS level 3D trend analysis')

})

output$tstrend_test_stat_table <- DT::renderDataTable({
  if(input$tstrend_test_stat_type=='DE trend genes'){
    x <- DDD.data$DE_tstrend_genes
  } else if(input$tstrend_test_stat_type=='DE trend transcripts'){
    x <- DDD.data$DE_tstrend_trans
  } else if(input$tstrend_test_stat_type=='DAS trend genes'){
    x <- DDD.data$DAS_tstrend_genes
  } else if(input$tstrend_test_stat_type=='DTU trend transcripts'){
    x <- DDD.data$DTU_tstrend_trans
  }
  
  if(is.null(x))
    return(NULL)
  x <- x[order(x$adj.pval,decreasing = F),]
  x <- x[1:pmin(input$top_n_tstrend_stat,nrow(x)),]
  x
})

##--------------->> Step 4: Profile plot-----------------
##---------------+ profile plot-----------------
output$ddd_tstrend_profile_plot_folder_text <- renderText({
  if(is.null(DDD.data$figure.folder))
    return(NULL)
  paste0('Figures are saved in:\n', DDD.data$figure.folder)
})



observe({
  if(is.null(DDD.data$samples))
    return(NULL)
  updateSelectInput(session = session,inputId = 'tstrend_profile_slice_group',choices = c('None',colnames(DDD.data$samples)),
                    selected = 'None')
})

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

tstrend_g_profiles <- eventReactive(input$tstrend_make_profile_plot,{
  if(is.null(input$tstrend_gene_id) | 
     input$tstrend_gene_id=="" |
     is.null(DDD.data$trans_TPM) | 
     is.null(DDD.data$mapping) |
     is.null(DDD.data$target_high))
    return(NULL)
  
  gene <- trimws(input$tstrend_gene_id)
  mapping <- DDD.data$mapping
  rownames(mapping) <- mapping$TXNAME
  
  if(input$tstrend_profile_data_type=='TPM'){
    data.exp <- DDD.data$trans_TPM
    reps <- DDD.data$samples$condition
  }
  
  if(input$tstrend_profile_data_type=='Read counts'){
    data.exp <- DDD.data$trans_counts
    reps <- DDD.data$samples_new$condition
  }
  
  if(input$tstrend_profile_filter_lowexpressed=='Yes'){
    data.exp <- data.exp[DDD.data$target_high$trans_high,]
    mapping <- mapping[DDD.data$target_high$trans_high,]
  }
  
  if(!(gene %in% mapping$GENEID)){
    showmessage(paste0('Gene ',gene,' is not in the dataset.'),duration = NULL)
    return(NULL)
  }
  
  if(input$tstrend_profile_slice_group=='None'){
    groups <- NULL 
  } else if(input$tstrend_profile_data_type=='TPM') {
    groups <-DDD.data$samples[,input$tstrend_profile_slice_group]
  } else {
    groups <-DDD.data$samples_new[,input$tstrend_profile_slice_group]
  }
  
  g.pr <- plotAbundance(data.exp = data.exp,
                        sliceProfile = T,
                        gene = gene,
                        groups = groups,
                        mapping = mapping,
                        genes.ann = DDD.data$genes.ann,
                        trans.ann = DDD.data$trans.ann,
                        trans.expressed = NULL,
                        reps = reps,
                        y.lab = input$tstrend_profile_data_type)+
    theme(axis.text.x = element_text(angle = input$tstrend.profile.x.rotate, 
                                     hjust = input$tstrend.profile.x.hjust,
                                     vjust = input$tstrend.profile.x.vjust))
  
  g.ps <- plotPS(data.exp = data.exp,
                 sliceProfile = T,
                 gene = gene,
                 mapping = mapping,
                 groups = groups,
                 genes.ann = DDD.data$genes.ann,
                 trans.ann = DDD.data$trans.ann,
                 trans.expressed = NULL,
                 reps = reps,
                 y.lab = 'PS')+
    theme(axis.text.x = element_text(angle = input$tstrend.profile.x.rotate, 
                                     hjust = input$tstrend.profile.x.hjust,
                                     vjust = input$tstrend.profile.x.vjust))
  return(list(g.pr=g.pr,g.ps=g.ps,gene=gene))
})

observeEvent(input$tstrend_make_profile_plot,{
  startmessage(paste0('Plot profile of gene: ',trimws(input$tstrend_gene_id)))
  output$tstrend_profile_plot_panel <- renderPlotly({
    if(is.null(tstrend_g_profiles()$g.pr))
      return(NULL)
    ggplotly(tstrend_g_profiles()$g.pr)
  })
  endmessage(paste0('Plot profile of gene: ',trimws(input$tstrend_gene_id)))
})


observeEvent(input$tstrend_make_profile_plot,{
  startmessage(paste0('Plot PS of gene: ',trimws(input$tstrend_gene_id)))
  output$tstrend_ps_plot_panel <- renderPlotly({
    if(is.null(tstrend_g_profiles()$g.ps))
      return(NULL)
    ggplotly(tstrend_g_profiles()$g.ps)
  })
  endmessage(paste0('Plot PS of gene: ',trimws(input$tstrend_gene_id)))
})

##update heigt and width
observe({
  if(is.null(tstrend_g_profiles()$g.pr$data))
    return(NULL)
  updateNumericInput(session,inputId = 'tstrend_ps_plot_width',
                     value = round((16+floor((length(unique(tstrend_g_profiles()$g.pr$data))-1)/15))/2.54,1))
})

observeEvent(input$tstrend_save_abundance_plot,{
  startmessage('Save profile plot')
  n <- length(unique(tstrend_g_profiles()$g.pr$data))-1
  gene <- tstrend_g_profiles()$gene
  folder2save <- paste0(DDD.data$figure.folder,'/Abundance plot')
  if(!file.exists(folder2save))
    dir.create(folder2save,recursive = T)
  
  png(paste0(folder2save,'/Abundance ',gene,'.png'),
      width = input$tstrend_ps_plot_width,height = input$tstrend_ps_plot_height,
      res=input$tstrend_ps_plot_res, units = 'in')
  print(tstrend_g_profiles()$g.pr)
  dev.off()
  
  pdf(paste0(folder2save,'/Abundance ',gene,'.pdf'),
      width = input$tstrend_ps_plot_width,height = input$tstrend_ps_plot_height)
  print(tstrend_g_profiles()$g.pr)
  dev.off()
  endmessage(paste0('Profile plot is saved in folder: ',DDD.data$figure.folder))
  
})


observeEvent(input$tstrend_save_abundance_plot,{
  startmessage('Save PS plot')
  n <- length(unique(tstrend_g_profiles()$g.ps$data))
  gene <- tstrend_g_profiles()$gene
  folder2save <- paste0(DDD.data$figure.folder,'/PS plot')
  if(!file.exists(folder2save))
    dir.create(folder2save,recursive = T)
  
  png(paste0(folder2save,'/PS ',gene,'.png'),
      width = input$tstrend_ps_plot_width,height = input$tstrend_ps_plot_height,
      res=input$tstrend_ps_plot_res, units = 'in')
  print(tstrend_g_profiles()$g.ps)
  dev.off()
  
  pdf(paste0(folder2save,'/PS ',gene,'.pdf'),
      width = input$tstrend_ps_plot_width,height = input$tstrend_ps_plot_height)
  print(tstrend_g_profiles()$g.ps)
  dev.off()
  
  endmessage(paste0('PS plot is saved in folder: ',DDD.data$figure.folder))
  
})


##preview saved plots
observeEvent(input$tstrend_save_abundance_plot_view, {
  gene <- trimws(input$tstrend_gene_id)
  file2save1 <- paste0(DDD.data$figure.folder,'/Abundance plot/Abundance ',gene,'.png')
  file2save2 <- paste0(DDD.data$figure.folder,'/PS plot/PS ',gene,'.png')
  
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

##---------------+ multiple plots----------------
observeEvent(input$tstrend_make_multiple_plot,{
  inFile <- input$tstrend_multiple_gene_input
  if (is.null(inFile))
    return(NULL)
  genes <- read.csv(inFile$datapath, header = T)
  genes <- trimws(as.vector(t(genes)))
  
  if(is.null(DDD.data$trans_TPM) | 
     is.null(DDD.data$mapping) |
     is.null(DDD.data$target_high))
    return(NULL)
  
  mapping <- DDD.data$mapping
  rownames(mapping) <- mapping$TXNAME
  
  if(input$tstrend_profile_data_type=='TPM'){
    data.exp <- DDD.data$trans_TPM
    reps <- DDD.data$samples$condition
  }
  
  if(input$tstrend_profile_data_type=='Read counts'){
    data.exp <- DDD.data$trans_counts
    reps <- DDD.data$samples_new$condition
  }
  
  if(input$tstrend_profile_filter_lowexpressed=='Yes'){
    data.exp <- data.exp[DDD.data$target_high$trans_high,]
    mapping <- mapping[DDD.data$target_high$trans_high,]
  }
  
  ###plot abundance
  folder2Abundance <- paste0(DDD.data$figure.folder,'/Profile plots/Abundance')
  if(!file.exists(folder2Abundance))
    dir.create(folder2Abundance,recursive = T)
  
  folder2PS <- paste0(DDD.data$figure.folder,'/Profile plots/PS')
  if(!file.exists(folder2PS))
    dir.create(folder2PS,recursive = T)
  
  startmessage(paste0('Plot profile of ',length(genes),' genes'))
  start.time <- Sys.time()
  withProgress(message = paste0('Plotting ',length(genes),' genes'),
               detail = 'This may take a while...', value = 0, {
                 graphics.off()
                 if(F){ 
                   incProgress(1/length(genes))
                   mclapply(genes,function(gene){
                     if(input$tstrend_multiple_plot_type %in% c('Abundance','Both')){
                       # message(paste0('Plotting abundance prfiles => Gene: ', gene, ' (',which(genes %in% gene), ' of ', length(genes),')'))
                       g.pr <- plotAbundance(data.exp = data.exp,
                                             gene = gene,
                                             mapping = mapping,
                                             genes.ann = DDD.data$genes.ann,
                                             trans.ann = DDD.data$trans.ann,
                                             trans.expressed = NULL,
                                             reps = reps,
                                             y.lab = input$tstrend_profile_data_type)+
                         theme(axis.text.x = element_text(angle = input$tstrend.profile.x.rotate, 
                                                          hjust = input$tstrend.profile.x.hjust,
                                                          vjust = input$tstrend.profile.x.vjust))
                       n <- length(unique(g.pr$data))-1
                       if(input$tstrend_multi_plot_format %in% c('png','both')){
                         png(paste0(folder2Abundance,'/Abundance ',gene,'.png'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height,
                             units = 'in',res = input$tstrend_ps_multiplot_res)
                         print(g.pr)
                         dev.off()
                       }
                       
                       if(input$tstrend_multi_plot_format %in% c('pdf','both')){
                         pdf(paste0(folder2Abundance,'/Abundance ',gene,'.pdf'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height)
                         print(g.pr)
                         dev.off()
                       }
                     }
                     
                     if(input$tstrend_multiple_plot_type %in% c('PS','Both')){
                       # message(paste0('Plotting PS prfiles => Gene: ', gene, ' (',which(genes %in% gene), ' of ', length(genes),')'))
                       g.ps <- plotPS(data.exp = data.exp,
                                      gene = gene,
                                      mapping = mapping,
                                      genes.ann = DDD.data$genes.ann,
                                      trans.ann = DDD.data$trans.ann,
                                      trans.expressed = NULL,
                                      reps = reps,
                                      y.lab = 'PS')+
                         theme(axis.text.x = element_text(angle = input$tstrend.profile.x.rotate, 
                                                          hjust = input$tstrend.profile.x.hjust,
                                                          vjust = input$tstrend.profile.x.vjust))
                       n <- length(unique(g.ps$data))
                       if(input$tstrend_multi_plot_format %in% c('png','both')){
                         png(paste0(folder2PS,'/PS ',gene,'.png'),
                             width = input$ps.multiplot.width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height,
                             units = 'in',res = input$tstrend_ps_multiplot_res)
                         print(g.ps)
                         dev.off()
                       }
                       if(input$multi.plot.format %in% c('pdf','both')){
                         pdf(paste0(folder2PS,'/PS ',gene,'.pdf'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height)
                         print(g.ps)
                         dev.off()
                       }
                     }
                   },mc.cores = DDD.data$num_cores)
                 } else {
                   lapply(genes,function(gene){
                     incProgress(1/length(genes))
                     if(input$tstrend_multiple_plot_type %in% c('Abundance','Both')){
                       # message(paste0('Plotting abundance prfiles => Gene: ', gene, ' (',which(genes %in% gene), ' of ', length(genes),')'))
                       g.pr <- plotAbundance(data.exp = data.exp,
                                             gene = gene,
                                             mapping = mapping,
                                             genes.ann = DDD.data$genes.ann,
                                             trans.ann = DDD.data$trans.ann,
                                             trans.expressed = NULL,
                                             reps = reps,
                                             y.lab = input$tstrend_profile_data_type)+
                         theme(axis.text.x = element_text(angle = input$tstrend.profile.x.rotate, 
                                                          hjust = input$tstrend.profile.x.hjust,
                                                          vjust = input$tstrend.profile.x.vjust))
                       n <- length(unique(g.pr$data))-1
                       if(input$tstrend_multi_plot_format %in% c('png','both')){
                         png(paste0(folder2Abundance,'/Abundance ',gene,'.png'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height,
                             units = 'in',res = input$tstrend_ps_multiplot_res)
                         print(g.pr)
                         dev.off()
                       }
                       
                       if(input$tstrend_multi_plot_format %in% c('pdf','both')){
                         pdf(paste0(folder2Abundance,'/Abundance ',gene,'.pdf'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height)
                         print(g.pr)
                         dev.off()
                       }
                     }
                     
                     if(input$tstrend_multiple_plot_type %in% c('PS','Both')){
                       # message(paste0('Plotting PS prfiles => Gene: ', gene, ' (',which(genes %in% gene), ' of ', length(genes),')'))
                       g.ps <- plotPS(data.exp = data.exp,
                                      gene = gene,
                                      mapping = mapping,
                                      genes.ann = DDD.data$genes.ann,
                                      trans.ann = DDD.data$trans.ann,
                                      trans.expressed = NULL,
                                      reps = reps,
                                      y.lab = 'PS')+
                         theme(axis.text.x = element_text(angle = input$tstrend.profile.x.rotate, 
                                                          hjust = input$tstrend.profile.x.hjust,
                                                          vjust = input$tstrend.profile.x.vjust))
                       n <- length(unique(g.ps$data))
                       if(input$tstrend_multi_plot_format %in% c('png','both')){
                         png(paste0(folder2PS,'/PS ',gene,'.png'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height,
                             units = 'in',res = input$tstrend_ps_multiplot_res)
                         print(g.ps)
                         dev.off()
                       }
                       if(input$multi.plot.format %in% c('pdf','both')){
                         pdf(paste0(folder2PS,'/PS ',gene,'.pdf'),
                             width = input$tstrend_ps_multiplot_width+floor(n/15)/2.54,
                             height = input$tstrend_ps_multiplot_height)
                         print(g.ps)
                         dev.off()
                       }
                     }
                   })
                 }
                 
                 ###
               })
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  showmessage(format(time.taken),showNoteify = F)
  endmessage(paste0('Plot profile of ',length(genes),' genes'))

})

# target2save <- reactive({
#   if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans) | is.null(DDD.data$DTU_trans))
#     return(NULL)
#   DEgenes <- unique(DDD.data$DE_genes$target)
#   DASgenes <- unique(DDD.data$DAS_genes$target)
#   DEtrans <- unique(DDD.data$DE_trans$target)
#   genesDEtrans <- unique(DDD.data$mapping[DEtrans,'GENEID'])
#   DTUtrans <- unique(DDD.data$DTU_trans$target)
#   genesDTUtrans <- unique(DDD.data$mapping[DTUtrans,'GENEID'])
#   x <- list(DEgenes=DEgenes,DASgenes=DASgenes,genesDEtrans=genesDEtrans,genesDTUtrans=genesDTUtrans)
#   n <- max(sapply(x, length))
#   y <- lapply(x,function(i){
#     c(i,rep(NA,n-length(i)))
#   })
#   z <- do.call(cbind,y)
#   colnames(z) <- c('DE genes','DAS genes','Genes of DE transcripts','Genes of DTU transcripts')
#   return(z)
# })

observeEvent(input$page_before_tstrend, {
  newtab <- switch(input$tabs, "ddd" = "tstrend","tstrend" = "ddd")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_tstrend, {
  newtab <- switch(input$tabs, "function" = "tstrend","tstrend" = "function")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

