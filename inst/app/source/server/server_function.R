observeEvent(input$tabs,{
  if(input$tabs=='function'){
    text2show <- 'Page: Functional plot'
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
    showmessage(text = text2show,showNoteify = F)
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
  }
})

##----------------        >> Heatmap        ------------------
observe({
  if(is.null(DDD.data$params_list$dist_method))
    return(NULL)
  updateSelectInput(session = session,
                    inputId = 'dist_method',selected = DDD.data$params_list$dist_method)
})

observe({
  if(is.null(DDD.data$params_list$cluster_method))
    return(NULL)
  updateSelectInput(session = session,
                    inputId = 'cluster_method',selected = DDD.data$params_list$cluster_method)
})

observe({
  if(is.null(DDD.data$params_list$cluster_number))
    return(NULL)
  updateNumericInput(session = session,
                    inputId = 'cluster_number',value = DDD.data$params_list$cluster_number)
})

# observe({
#   if(input$whether_tstrend=='No')
#     return(NULL)
#   updateSelectInput(session = session,
#                     inputId = 'heatmap.target.type',
#                     choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts',
#                                 'DE tstrend genes','DAS tstrend genes','DE tstrend transcripts','DTU tstrend transcripts'))
# })

heatmap.targets <- eventReactive(input$plot.heatmap,{

  # if(input$heatmap.select.or.upload=='Select targets'){
    if(input$heatmap.target.type=='DE genes'){
      targets <- unique(DDD.data$DE_genes$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$genes_TPM[targets,]
      column_title <- paste0(length(targets),' DE genes')
    }
    
    if(input$heatmap.target.type=='DAS genes'){
      targets <- intersect(unique(DDD.data$DAS_genes$target),unique(DDD.data$DE_genes$target))
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$genes_TPM[targets,]
      column_title <- paste0(length(targets),' DE+DAS genes (DAS only are filtered)')
    }
    
    if(input$heatmap.target.type=='DE transcripts'){
      targets <- unique(DDD.data$DE_trans$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$trans_TPM[targets,]
      column_title <- paste0(length(targets),' DE transcripts')
    }
    
    if(input$heatmap.target.type=='DTU transcripts'){
      targets <- intersect(unique(DDD.data$DTU_trans$target),unique(DDD.data$DE_trans$target))
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$trans_TPM[targets,]
      column_title <- paste0(length(targets),' DE+DTU transcripts (DTU only are filtered)')
    }
  
  ###############ts trend
  
  # if(input$heatmap.select.or.upload=='Select targets'){
  if(input$heatmap.target.type=='DE tstrend genes'){
    targets <- unique(DDD.data$DE_tstrend_genes$target)
    if(length(targets)<5){
      showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
      return(NULL)
    }
    data2heatmap <- DDD.data$genes_TPM[targets,]
    column_title <- paste0(length(targets),' DE tstrend genes')
  }
  
  if(input$heatmap.target.type=='DAS tstrend genes'){
    targets <- intersect(unique(DDD.data$DAS_tstrend_genes$target),unique(DDD.data$DE_genes$target))
    if(length(targets)<5){
      showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
      return(NULL)
    }
    data2heatmap <- DDD.data$genes_TPM[targets,]
    column_title <- paste0(length(targets),' DE+DAS tstrend genes (DAS only are filtered)')
  }
  
  if(input$heatmap.target.type=='DE tstrend transcripts'){
    targets <- unique(DDD.data$DE_tstrend_trans$target)
    if(length(targets)<5){
      showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
      return(NULL)
    }
    data2heatmap <- DDD.data$trans_TPM[targets,]
    column_title <- paste0(length(targets),' DE tstrend transcripts')
  }
  
  if(input$heatmap.target.type=='DTU tstrend transcripts'){
    targets <- intersect(unique(DDD.data$DTU_tstrend_trans$target),unique(DDD.data$DE_trans$target))
    if(length(targets)<5){
      showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
      return(NULL)
    }
    data2heatmap <- DDD.data$trans_TPM[targets,]
    column_title <- paste0(length(targets),' DE+DTU tstrend transcripts (DTU only are filtered)')
  }
  
  list(data2heatmap=data2heatmap,column_title=column_title)
})


heatmap.g <- eventReactive(input$plot.heatmap,{
  if(is.null(heatmap.targets()))
    return(NULL)
  startmessage(text = 'Plot heatmap')
  withProgress(message = 'Making heat-map...', value = 0, {
    incProgress(0.1)
    data2plot <- rowmean(x = t(heatmap.targets()$data2heatmap),
                         group = DDD.data$samples$condition,
                         reorder = F)
    # if(input$TPM_or_Zscore=='TPM'){
    #   data2plot <- t(data2plot)
    # } else {
    data2plot <- t(scale(data2plot))
    # }
    incProgress(0.2)
    hc.dist <- dist(data2plot,method = input$dist_method)
    hc <- fastcluster::hclust(hc.dist,method = input$cluster_method)
    clusters <- cutree(hc, k = pmin(input$cluster_number,nrow(hc$merge)+1))
    clusters <- reorderClusters(clusters = clusters,dat = data2plot)
    incProgress(0.3)
    g <- Heatmap(as.matrix(data2plot), 
                 col =  colorRampPalette(c(input$color4_heatmap_neg,input$color4_heatmap_0,input$color4_heatmap_pos))(100),
                 name = 'Z-scores',
                 cluster_rows = TRUE,
                 clustering_method_rows=input$cluster_method,
                 row_dend_reorder = T,
                 show_row_names = ifelse(input$show_row_labels=='Yes',T,F),
                 show_column_names = ifelse(input$show_column_labels=='Yes',T,F),
                 # show_column_names = ifelse(ncol(data2plot)>10,F,T),
                 cluster_columns = FALSE,
                 split=clusters,
                 column_title= heatmap.targets()$column_title)
    incProgress(0.8)
  })
  list(g=g,clusters=clusters)
})

output$heatmap.plot.panel <- renderPlot({
  if(is.null(heatmap.g()))
    return(NULL)
  draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
  endmessage(text = 'Plot heatmap')
})

##update heigt and width
observe({
  if(is.null(DDD.data$samples$condition))
    return(NULL)
  updateNumericInput(session,inputId = 'heatmap.plot.width',
                     value = round(pmax(10,1.5*length(unique(DDD.data$samples$condition)))/2.54,1))
})



##save plot to working directory
observeEvent(input$save_heatmap,{
  # if(save_heatmap_flag$download_flag==0)
  #   return(NULL)
  startmessage(text = 'Save heatmap plot')
  withProgress(message = 'Saving heatmap...', detail = 'This may take a while...',
               value = 0.5, {
                 ##DE vs DTU transcript
                 figure.name <- paste0(DDD.data$figure.folder,'/Heatmap ', input$heatmap.target.type)
                 png(paste0(figure.name,'.png'),
                     width = input$heatmap.plot.width,height = input$heatmap.plot.height,
                     res=input$heatmap.plot.res, units = 'in')
                 draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                 dev.off()
                 
                 pdf(paste0(figure.name,'.pdf'),
                     width = input$heatmap.plot.width,height = input$heatmap.plot.height)
                 draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                 dev.off()
                 
                 ###save target list
                 clusters <- heatmap.g()$clusters
                 x <- split(names(clusters),clusters)
                 x <- lapply(names(x),function(i){
                   data.frame(Clusters=i,Targets=x[[i]])
                 })
                 x <- do.call(rbind,x)
                 colnames(x) <- c('Clusters','Targets')
                 write.csv(x,file=paste0(DDD.data$result.folder,'/Target in each cluster heatmap ', 
                                         heatmap.targets()$column_title,'.csv'),row.names = F)
                 showmessage(paste0('Table is saved in folder: ',DDD.data$result.folder))
               })
  endmessage(text = 'Save heatmap plot')
  showmessage(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})


###preview saved plot
observeEvent(input$save_heatmap_view, {
  file2save <- paste0(DDD.data$figure.folder,'/Heatmap ', input$heatmap.target.type,'.png')
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
               style="width: 80%; display: block; margin-left: auto; margin-right: auto; 
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='m'
      ))
  }

})

target2save <- reactive({
  if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans) | is.null(DDD.data$DTU_trans))
    return(NULL)
  DEgenes <- unique(DDD.data$DE_genes$target)
  DASgenes <- unique(DDD.data$DAS_genes$target)
  DEtrans <- unique(DDD.data$DE_trans$target)
  genesDEtrans <- unique(DDD.data$mapping[DEtrans,'GENEID'])
  DTUtrans <- unique(DDD.data$DTU_trans$target)
  genesDTUtrans <- unique(DDD.data$mapping[DTUtrans,'GENEID'])
  x <- list(DEgenes=DEgenes,DASgenes=DASgenes,genesDEtrans=genesDEtrans,genesDTUtrans=genesDTUtrans)
  n <- max(sapply(x, length))
  y <- lapply(x,function(i){
    c(i,rep(NA,n-length(i)))
  })
  z <- do.call(cbind,y)
  colnames(z) <- c('DE genes','DAS genes','Genes of DE transcripts','Genes of DTU transcripts')
  return(z)
})



##----------------      >> TS trend heatmap        ------------------
observe({
  if(is.null(DDD.data$params_list$tstrend_dist_method))
    return(NULL)
  updateSelectInput(session = session,
                    inputId = 'tstrend_dist_method',selected = DDD.data$params_list$tstrend_dist_method)
})

observe({
  if(is.null(DDD.data$params_list$tstrend_cluster_method))
    return(NULL)
  updateSelectInput(session = session,
                    inputId = 'tstrend_cluster_method',selected = DDD.data$params_list$tstrend_cluster_method)
})

observe({
  if(is.null(DDD.data$params_list$tstrend_cluster_number))
    return(NULL)
  updateNumericInput(session = session,
                     inputId = 'tstrend_cluster_number',value = DDD.data$params_list$tstrend_cluster_number)
})

observe({
  if(is.null(DDD.data$tstrend_Contrastgroups))
    return(NULL)
  updateSelectInput(session = session,inputId = 'tstrend_heatmap_contrast_groups',
                    choices = c('All',names(DDD.data$tstrend_Contrastgroups)),
                    selected = 'All')
})

tstrend.heatmap.targets <- eventReactive(input$plot.tstrend.heatmap,{
  if(is.null(input$tstrend_time_group) | 
     is.null(input$tstrend_heatmap_contrast_groups) | 
     is.null(DDD.data$samples) |
     input$tstrend_time_group==''
     ) 
    return(NULL)
  
  ### DE tstrend genes
  if(input$tstrend.heatmap.target.type=='DE tstrend genes'){
    if(input$tstrend_heatmap_contrast_groups=='All'){
      targets <- unique(DDD.data$DE_tstrend_genes$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$genes_TPM[targets,]
      group <- DDD.data$samples$condition
    } else {
      dat <- subset(DDD.data$DE_tstrend_genes,contrast==input$tstrend_heatmap_contrast_groups)
      targets <- unique(dat$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      idx1 <- DDD.data$samples[,input$tstrend_time_group]
      idx2 <- unlist(strsplit(input$tstrend_heatmap_contrast_groups,'='))
      data2heatmap <- DDD.data$genes_TPM[targets,idx1 %in% idx2]
      group <- DDD.data$samples$condition[idx1 %in% idx2]
    }
    column_title <- paste0(length(targets),' TS DE trend genes in contrast ',input$tstrend_heatmap_contrast_groups)
  }
  
  ### DAS tstrend genes
  if(input$tstrend.heatmap.target.type=='DAS tstrend genes'){
    if(input$tstrend_heatmap_contrast_groups=='All'){
      targets <- unique(DDD.data$DAS_tstrend_genes$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$genes_TPM[targets,]
      group <- DDD.data$samples$condition
    } else {
      dat <- subset(DDD.data$DAS_tstrend_genes,contrast==input$tstrend_heatmap_contrast_groups)
      targets <- unique(dat$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      idx1 <- DDD.data$samples[,input$tstrend_time_group]
      idx2 <- unlist(strsplit(input$tstrend_heatmap_contrast_groups,'='))
      data2heatmap <- DDD.data$genes_TPM[targets,idx1 %in% idx2]
      group <- DDD.data$samples$condition[idx1 %in% idx2]
    }
    column_title <- paste0(length(targets),' TS DAS trend genes in contrast ',input$tstrend_heatmap_contrast_groups)
  }
  
  ### DE tstrend transcripts
  if(input$tstrend.heatmap.target.type=='DE tstrend transcripts'){
    if(input$tstrend_heatmap_contrast_groups=='All'){
      targets <- unique(DDD.data$DE_tstrend_trans$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$trans_TPM[targets,]
      group <- DDD.data$samples$condition
    } else {
      dat <- subset(DDD.data$DE_tstrend_trans,contrast==input$tstrend_heatmap_contrast_groups)
      targets <- unique(dat$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      idx1 <- DDD.data$samples[,input$tstrend_time_group]
      idx2 <- unlist(strsplit(input$tstrend_heatmap_contrast_groups,'='))
      data2heatmap <- DDD.data$trans_TPM[targets,idx1 %in% idx2]
      group <- DDD.data$samples$condition[idx1 %in% idx2]
    }
    column_title <- paste0(length(targets),' TS DE trend transcripts in contrast ',input$tstrend_heatmap_contrast_groups)
  }
  
  ### DTU tstrend transcripts
  if(input$tstrend.heatmap.target.type=='DTU tstrend transcripts'){
    if(input$tstrend_heatmap_contrast_groups=='All'){
      targets <- unique(DDD.data$DTU_tstrend_trans$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      data2heatmap <- DDD.data$trans_TPM[targets,]
      group <- DDD.data$samples$condition
    } else {
      dat <- subset(DDD.data$DTU_tstrend_trans,contrast==input$tstrend_heatmap_contrast_groups)
      targets <- unique(dat$target)
      if(length(targets)<5){
        showmessage('Too few targets to make heatmap!',duration = 20,type = 'error')
        return(NULL)
      }
      idx1 <- DDD.data$samples[,input$tstrend_time_group]
      idx2 <- unlist(strsplit(input$tstrend_heatmap_contrast_groups,'='))
      data2heatmap <- DDD.data$trans_TPM[targets,idx1 %in% idx2]
      group <- DDD.data$samples$condition[idx1 %in% idx2]
    }
    column_title <- paste0(length(targets),' TS DTU trend transcripts in contrast ',input$tstrend_heatmap_contrast_groups)
  }
  list(data2heatmap=data2heatmap,group=group,column_title=column_title)
})


tstrend.heatmap.g <- eventReactive(input$plot.tstrend.heatmap,{
  if(is.null(tstrend.heatmap.targets()))
    return(NULL)
  
  startmessage(text = 'Plot heatmap')
  withProgress(message = 'Making heat-map...', value = 0, {
    incProgress(0.1)
    data2plot <- rowmean(x = t(tstrend.heatmap.targets()$data2heatmap),
                         group = tstrend.heatmap.targets()$group,
                         reorder = F)
    
    
    # if(input$TPM_or_Zscore=='TPM'){
    #   data2plot <- t(data2plot)
    # } else {
    data2plot <- t(scale(data2plot))
    # }
    incProgress(0.2)
    hc.dist <- dist(data2plot,method = input$tstrend_dist_method)
    hc <- fastcluster::hclust(hc.dist,method = input$tstrend_cluster_method)
    clusters <- cutree(hc, k = pmin(input$tstrend_cluster_number,nrow(hc$merge)+1))
    clusters <- reorderClusters(clusters = clusters,dat = data2plot)
    incProgress(0.3)
    
    g <- Heatmap(as.matrix(data2plot), 
                 col =  colorRampPalette(c(input$color4_tstrend_heatmap_neg,
                                           input$color4_tstrend_heatmap_0,
                                           input$color4_tstrend_heatmap_pos))(100),
                 name = 'Z-scores',
                 cluster_rows = TRUE,
                 clustering_method_rows=input$tstrend_cluster_method,
                 row_dend_reorder = T,
                 show_row_names = ifelse(input$tstrend_show_row_labels=='Yes',T,F),
                 show_column_names = ifelse(input$tstrend_show_column_labels=='Yes',T,F),
                 # show_column_names = ifelse(ncol(data2plot)>10,F,T),
                 cluster_columns = FALSE,
                 split=clusters,
                 column_title= tstrend.heatmap.targets()$column_title)
    incProgress(0.8)
  })
  list(g=g,clusters=clusters)
})

output$tstrend.heatmap.plot.panel <- renderPlot({
  if(is.null(tstrend.heatmap.g()))
    return(NULL)
  draw(tstrend.heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
  endmessage(text = 'Plot heatmap')
})

##update heigt and width
observe({
  if(is.null(tstrend.heatmap.targets()))
    return(NULL)
  updateNumericInput(session,inputId = 'tstrend.heatmap.plot.width',
                     value = round(pmax(10,1.5*length(unique(tstrend.heatmap.targets()$group)))/2.54,1))
})


##save plot to working directory
observeEvent(input$save_tstrend_heatmap,{
  if(is.null(tstrend.heatmap.targets()))
    return(NULL)
  # if(save_heatmap_flag$download_flag==0)
  #   return(NULL)
  startmessage(text = 'Save heatmap plot')
  withProgress(message = 'Saving heatmap...', detail = 'This may take a while...',
               value = 0.5, {
                 ##DE vs DTU transcript
                 figure.name <- paste0(DDD.data$figure.folder,'/Heatmap ', tstrend.heatmap.targets()$column_title)
                 png(paste0(figure.name,'.png'),
                     width = input$tstrend.heatmap.plot.width,height = input$tstrend.heatmap.plot.height,
                     res=input$tstrend.heatmap.plot.res, units = 'in')
                 draw(tstrend.heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                 dev.off()
                 
                 pdf(paste0(figure.name,'.pdf'),
                     width = input$tstrend.heatmap.plot.width,height = input$tstrend.heatmap.plot.height)
                 draw(tstrend.heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                 dev.off()
                 
                 ###save target list
                 clusters <- tstrend.heatmap.g()$clusters
                 x <- split(names(clusters),clusters)
                 x <- lapply(names(x),function(i){
                   data.frame(Clusters=i,Targets=x[[i]])
                 })
                 x <- do.call(rbind,x)
                 colnames(x) <- c('Clusters','Targets')
                 write.csv(x,file=paste0(DDD.data$result.folder,'/Target in each cluster heatmap ', 
                                         tstrend.heatmap.targets()$column_title,'.csv'),row.names = F)
                 showmessage(paste0('Table is saved in folder: ',DDD.data$result.folder))
               })
  endmessage(text = 'Save heatmap plot')
  showmessage(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
})



###preview saved plot
observeEvent(input$save_tstrend_heatmap_view, {
  if(is.null(tstrend.heatmap.targets()))
    return(NULL)
  file2save <- paste0(DDD.data$figure.folder,'/Heatmap ', tstrend.heatmap.targets()$column_title,'.png')
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
               style="width: 80%; display: block; margin-left: auto; margin-right: auto; 
               border:1px solid #000000"),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='m'
      ))
  }
  
})

##--------------->> GO plot----------------
output$download_target_list <- downloadHandler(
  filename=function(){
    'Gene list of 3D genes and transcripts.csv'
  },
  content=function(file){
    write.csv(target2save(),file,row.names = F,quote = F,na = '')
  }
)


go.table <- reactive({
  inFile <- input$go.annotation.input
  if (is.null(inFile))
    return(NULL)
  go.table <- suppressMessages(readr::read_csv(inFile$datapath))
  go.table
})

output$go.table <- DT::renderDataTable({
  if(is.null(go.table()))
    return(NULL)
  showmessage('Show the uploaded GO annotation table',showNoteify = F)
  go.table()
},options = list(
  scrollX=T,
  columnDefs = list(list(className = 'dt-left',
                         targets = "_all"))))


observe({
  if(is.null(go.table()))
    return(NULL)
  updateSelectInput(session = session,inputId = 'select.go.table.column',label = 'Select column of statistics',
                    choices = colnames(go.table())[-c(1:2)])
})

# output$show.go.table.column <- renderUI({
#   if(is.null(go.table()))
#     return(NULL)
#   span(
#     selectInput(inputId = 'select.go.table.column',label = 'Select column of statistics',
#                 choices = colnames(go.table())[-c(1:2)])
#   )
# })

data2goplot <- eventReactive(input$make.go.plot,{
  col.idx<- input$select.go.table.column
  data2plot <- go.table()[,c(1,2,which(colnames(go.table())==col.idx))]
  colnames(data2plot) <- c('Category','Term','Value')
  ####
  data2plot <- by(data2plot,data2plot$Category,function(x){
    x[order(x$Value,decreasing = T),]
  })
  data2plot <- do.call(rbind,data2plot)
  idx<-sapply(data2plot$Term,function(x){
    x<-as.vector(x)
    if(nchar(x)>input$go.text.cut)
      x<-paste0(substr(x,1,input$go.text.cut),'...')
    x
  })
  data2plot$Term <- factor(idx,levels = rev(idx))
  data2plot
})

row.height <- eventReactive(input$make.go.plot,{
  if(is.null(go.table()))
    return(NULL)
  nrow(go.table())
})

observeEvent(input$make.go.plot,{
  output$show.go.plot.panel <- renderUI({
    plotlyOutput('go.plot',height = pmax(row.height()*20,350))
  })
  
})

observeEvent(input$make.go.plot,{
  startmessage('Plot GO annotation')
  output$go.plot <- renderPlotly({
    if(is.null(data2goplot()))
      return(NULL)
    Category.idx <- unique(data2goplot()$Category)
    color.idx <- gg.color.hue(length(Category.idx))
    g.list <- lapply(1:length(Category.idx),function(i){
      x <- subset(data2goplot(),Category==Category.idx[i])
      x$Term <- factor(x$Term,levels = rev(x$Term))
      g <- ggplot(x,aes(x=Term,y=Value))+
        geom_bar(stat='identity',aes(fill=Category))+
        scale_fill_manual(values = color.idx[i])+
        coord_flip(ylim = c(0,max(data2goplot()$Value))) +
        theme_bw()+
        theme(axis.title.y = element_blank())+
        labs(y=input$select.go.table.column,title=paste0('GO annotation: ',input$go_plot_label))
      if(i>1)
        g <- g+guides(fill=guide_legend(title=NULL))
      if(i!=length(Category.idx))
        g <- g+theme(axis.title.x = element_blank())
      ggplotly(g)
    })
    heights <- sapply(1:length(Category.idx),function(i) length(g.list[[i]]$x$data[[1]]$y))
    heights <- heights/sum(heights)
    subplot(g.list,nrows=length(g.list),heights = heights, titleX = TRUE)
  })
  endmessage('Plot GO annotation')
})

## update height
observe({
  if(is.null(row.height()))
    return(NULL)
  updateNumericInput(session,inputId = 'go.plot.height',
                     value = round(pmax(10,row.height()*0.4)/2.54,1))
})


observeEvent(input$save_go_plot,{
  if(is.null(data2goplot()))
    return(NULL)
  startmessage('Save GO annotation plot')
  col.idx<- input$select.go.table.column
  data2plot <- go.table()[,c(1,2,which(colnames(go.table())==col.idx))]
  colnames(data2plot) <- c('Category','Term','Value')
  ####
  data2plot <- by(data2plot,data2plot$Category,function(x){
    x[order(x$Value,decreasing = T),]
  })
  data2plot <- do.call(rbind,data2plot)
  idx<-sapply(data2plot$Term,function(x){
    x<-as.vector(x)
    if(nchar(x)>input$go.text.cut)
      x<-paste0(substr(x,1,input$go.text.cut),'...')
    x
  })
  data2plot$Term <- factor(idx,levels = rev(idx))
  
  g <- ggplot(data2plot,aes(x=Term,y=Value))+
    geom_bar(stat='identity',aes(fill=Category))+
    coord_flip() +
    theme_bw()+
    theme(axis.title.y = element_blank())+
    facet_grid(Category~.,scales = 'free_y',space='free_y')+
    labs(y=col.idx,title=paste0('GO annotation ',input$go_plot_label))
  png(paste0(DDD.data$figure.folder,'/',input$go_plot_label,' GO annotation plot.png'),
      width = input$go.plot.width,height = input$go.plot.height,
      res=input$go.plot.res, units = 'in')
  print(g)
  dev.off()
  pdf(paste0(DDD.data$figure.folder,'/',input$go_plot_label,' GO annotation plot.pdf'),
      width = input$go.plot.width,height = input$go.plot.height)
  print(g)
  dev.off()
  endmessage('Save GO annotation plot')
  showmessage(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
})

### Preview saved plot
observeEvent(input$save_go_plot_view, {
  file2save <- paste0(DDD.data$figure.folder,'/',input$go_plot_label,' GO annotation plot.png')
  
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
               border:1px solid #000000" ),
      easyClose = TRUE,
      footer = modalButton("Cancel"),
      size='l'
      ))
  }
  
})

observeEvent(input$page_before_function, {
  newtab <- switch(input$tabs, "function" = "tstrend","tstrend" = "function")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_function, {
  newtab <- switch(input$tabs, "function" = "TSIS","TSIS" = "function")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

# observeEvent(input$page_before_function, {
#   newtab <- switch(input$tabs, "function" = "ddd","ddd" = "function")
#   updateTabItems(session, "tabs", newtab)
#   shinyjs::runjs("window.scrollTo(0, 50)")
# })
# 
# observeEvent(input$page_after_function, {
#   newtab <- switch(input$tabs, "function" = "TSIS","TSIS" = "function")
#   updateTabItems(session, "tabs", newtab)
#   shinyjs::runjs("window.scrollTo(0, 50)")
# })

