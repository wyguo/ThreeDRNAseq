##---------------->>    Advanced plot   <<-----------------
##----------------        heatmap        ------------------

heatmap.targets <- eventReactive(input$plot.heatmap,{

  # if(input$heatmap.select.or.upload=='Select targets'){
    if(input$heatmap.target.type=='DE genes'){
      targets <- unique(DDD.data$DE_genes$target)
      data2heatmap <- DDD.data$genes_TPM[targets,]
      column_title <- paste0(length(targets),' DE genes')
    }
    
    if(input$heatmap.target.type=='DAS genes'){
      targets <- intersect(unique(DDD.data$DAS_genes$target),unique(DDD.data$DE_genes$target))
      data2heatmap <- DDD.data$genes_TPM[targets,]
      column_title <- paste0(length(targets),' DE+DAS genes (DAS only are filtered)')
    }
    
    if(input$heatmap.target.type=='DE transcripts'){
      targets <- unique(DDD.data$DE_trans$target)
      data2heatmap <- DDD.data$trans_TPM[targets,]
      column_title <- paste0(length(targets),' DE transcripts')
    }
    
    if(input$heatmap.target.type=='DTU transcripts'){
      targets <- intersect(unique(DDD.data$DTU_trans$target),unique(DDD.data$DE_trans$target))
      data2heatmap <- DDD.data$trans_TPM[targets,]
      column_title <- paste0(length(targets),' DE+DTU transcripts (DTU only are filtered)')
    }
  # }
  # if(input$heatmap.select.or.upload=='Upload target list'){
  #   inFile <- input$heatmap.target.list.input
  #   if (is.null(inFile))
  #     return(NULL)
  #   targets <- read.csv(inFile$datapath, header = F)
  #   targets <- as.vector(t(targets))
  #   if(input$heatmap.targetlist.type=='Gene list'){
  #     data2heatmap <- DDD.data$genes_TPM[targets,]
  #     column_title <- paste0(length(targets),' Genes')
  #   } else {
  #     data2heatmap <- DDD.data$trans_TPM[targets,]
  #     column_title <- paste0(length(targets),' Transcripts')
  #   }
  # }
  list(data2heatmap=data2heatmap,column_title=column_title)
})


heatmap.g <- eventReactive(input$plot.heatmap,{
  if(is.null(heatmap.targets()))
    return(NULL)
  withProgress(message = 'Generate gene expression...', value = 0, {
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
    hc.dist <- dist(data2plot,method = input$dist.method)
    hc <- fastcluster::hclust(hc.dist,method = input$cluster.method)
    clusters <- cutree(hc, k = input$cluster.number)
    clusters <- reorderClusters(clusters = clusters,dat = data2plot)
    incProgress(0.3)
    g <- Heatmap(as.matrix(data2plot), 
                 col =  colorRampPalette(c(input$color4_heatmap_neg,input$color4_heatmap_0,input$color4_heatmap_pos))(100),
                 name = 'Z-scores',
                 cluster_rows = TRUE,
                 clustering_method_rows=input$cluster.method,
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
})

##update heigt and width
observe({
  updateNumericInput(session,inputId = 'heatmap.plot.width',
                     value = round(pmax(10,1*length(unique(DDD.data$samples$condition)))/2.54,1))
})

observeEvent(input$save.heatmap,{
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
                 message(paste0('Table is saved in folder: ',DDD.data$result.folder))
                 showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
               })
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  
})

##---------------profile plot-----------------
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

g.profiles <- eventReactive(input$make.profile.plot,{
  if(is.null(input$gene.id) | input$gene.id=="")
    return(NULL)
  gene <- trimws(input$gene.id)
  
  mapping <- DDD.data$mapping
  rownames(mapping) <- mapping$TXNAME
  
  if(is.null(DDD.data$trans_TPM)){
    load(paste0(DDD.data$data.folder,'/trans_TPM.RData'))
    DDD.data$trans_TPM <- trans_TPM
  }
  
  if(is.null(DDD.data$target_high)){
    load(paste0(DDD.data$data.folder,'/target_high.RData'))
    DDD.data$target_high <- target_high
  }
  
  if(input$profile.data.type=='TPM'){
    data.exp <- DDD.data$trans_TPM
    reps <- DDD.data$samples$condition
  }
  
  if(input$profile.data.type=='Read counts'){
    data.exp <- DDD.data$trans_counts
    reps <- DDD.data$samples_new$condition
  }
  
  if(input$profile.filter.lowexpressed=='Yes'){
    data.exp <- data.exp[DDD.data$target_high$trans_high,]
    mapping <- mapping[DDD.data$target_high$trans_high,]
  }
  
  g.pr <- plotAbundance(data.exp = data.exp,
                        gene = gene,
                        mapping = mapping,
                        genes.ann = DDD.data$genes.ann,
                        trans.ann = DDD.data$trans.ann,
                        trans.expressed = NULL,
                        reps = reps,
                        y.lab = input$profile.data.type)
  
  g.ps <- plotPS(data.exp = data.exp,
                 gene = gene,
                 mapping = mapping,
                 genes.ann = DDD.data$genes.ann,
                 trans.ann = DDD.data$trans.ann,
                 trans.expressed = NULL,
                 reps = reps,
                 y.lab = 'PS')
  
  
  
  
  return(list(g.pr=g.pr,g.ps=g.ps,gene=gene))
})

output$profile.plot.panel <- renderPlotly({
  if(is.null(g.profiles()$g.pr))
    return(NULL)
  ggplotly(g.profiles()$g.pr)
})

output$ps.plot.panel <- renderPlotly({
  if(is.null(g.profiles()$g.ps))
    return(NULL)
  ggplotly(g.profiles()$g.ps)
})

##update heigt and width
observe({
  updateNumericInput(session,inputId = 'ps.plot.width',
                     value = round((16+floor((length(unique(g.profiles()$g.pr$data))-1)/15))/2.54,1))
})

observeEvent(input$save.abundance.plot,{
  n <- length(unique(g.profiles()$g.pr$data))-1
  gene <- g.profiles()$gene
  folder2save <- paste0(DDD.data$figure.folder,'/Abundance plot')
  if(!file.exists(folder2save))
    dir.create(folder2save,recursive = T)
  
  png(paste0(folder2save,'/Abundance ',gene,'.png'),
      width = input$ps.plot.width,height = input$ps.plot.height,
      res=input$ps.plot.res, units = 'in')
  print(g.profiles()$g.pr)
  dev.off()
  
  pdf(paste0(folder2save,'/Abundance ',gene,'.pdf'),
      width = input$ps.plot.width,height = input$ps.plot.height)
  print(g.profiles()$g.pr)
  dev.off()
  
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  
})


observeEvent(input$save.ps.plot,{
  n <- length(unique(g.profiles()$g.ps$data))
  gene <- g.profiles()$gene
  folder2save <- paste0(DDD.data$figure.folder,'/PS plot')
  if(!file.exists(folder2save))
    dir.create(folder2save,recursive = T)
  
  png(paste0(folder2save,'/PS ',gene,'.png'),
      width = input$ps.plot.width,height = input$ps.plot.height,
      res=input$ps.plot.res, units = 'in')
  print(g.profiles()$g.ps)
  dev.off()
  
  pdf(paste0(folder2save,'/PS ',gene,'.pdf'),
      width = input$ps.plot.width,height = input$ps.plot.height)
  print(g.profiles()$g.ps)
  dev.off()
  
  message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
  
})


##---------------multiple plots----------------
observeEvent(input$make.multiple.plot,{
  inFile <- input$multiple.gene.input
  if (is.null(inFile))
    return(NULL)
  genes <- read.csv(inFile$datapath, header = T)
  genes <- trimws(as.vector(t(genes)))
  mapping <- DDD.data$mapping
  rownames(mapping) <- mapping$TXNAME
  
  if(input$profile.data.type=='TPM'){
    data.exp <- DDD.data$trans_TPM
    reps <- DDD.data$samples$condition
  }
  
  if(input$profile.data.type=='Read counts'){
    data.exp <- DDD.data$trans_counts
    reps <- DDD.data$samples_new$condition
  }
  
  if(input$profile.filter.lowexpressed=='Yes'){
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
  
  withProgress(message = paste0('Plotting ',length(genes),' genes'),
               detail = 'This may take a while...', value = 0, {
                 graphics.off()
                 for(gene in genes){
                   incProgress(1/length(genes))
                   if(input$multiple.plot.type %in% c('Abundance','Both')){
                     message(paste0('Plotting abundance prfiles => Gene: ', gene, ' (',which(genes %in% gene), ' of ', length(genes),')'))
                     g.pr <- plotAbundance(data.exp = data.exp,
                                           gene = gene,
                                           mapping = mapping,
                                           genes.ann = DDD.data$genes.ann,
                                           trans.ann = DDD.data$trans.ann,
                                           trans.expressed = NULL,
                                           reps = reps,
                                           y.lab = input$profile.data.type)
                     n <- length(unique(g.pr$data))-1
                     if(input$multi.plot.format %in% c('png','both')){
                       png(paste0(folder2Abundance,'/Abundance ',gene,'.png'),
                           width = input$ps.multiplot.width+floor(n/15)/2.54,
                           height = input$ps.multiplot.height,
                           units = 'in',res = input$ps.multiplot.res)
                       print(g.pr)
                       dev.off()
                     }
                     
                     if(input$multi.plot.format %in% c('pdf','both')){
                       pdf(paste0(folder2Abundance,'/Abundance ',gene,'.pdf'),
                           width = input$ps.multiplot.width+floor(n/15)/2.54,
                           height = input$ps.multiplot.height)
                       print(g.pr)
                       dev.off()
                     }
                   }
                   
                   if(input$multiple.plot.type %in% c('PS','Both')){
                     message(paste0('Plotting PS prfiles => Gene: ', gene, ' (',which(genes %in% gene), ' of ', length(genes),')'))
                     g.ps <- plotPS(data.exp = data.exp,
                                    gene = gene,
                                    mapping = mapping,
                                    genes.ann = DDD.data$genes.ann,
                                    trans.ann = DDD.data$trans.ann,
                                    trans.expressed = NULL,
                                    reps = reps,
                                    y.lab = 'PS')
                     n <- length(unique(g.ps$data))
                     if(input$multi.plot.format %in% c('png','both')){
                       png(paste0(folder2PS,'/PS ',gene,'.png'),
                           width = input$ps.multiplot.width+floor(n/15)/2.54,
                           height = input$ps.multiplot.height,
                           units = 'in',res = input$ps.multiplot.res)
                       print(g.ps)
                       dev.off()
                     }
                     if(input$multi.plot.format %in% c('pdf','both')){
                       pdf(paste0(folder2PS,'/PS ',gene,'.pdf'),
                           width = input$ps.multiplot.width+floor(n/15)/2.54,
                           height = input$ps.multiplot.height)
                       print(g.ps)
                       dev.off()
                     }
                   }
                   
                 }
               })
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

##---------------GO plot----------------
output$download_target_list <- downloadHandler(
  filename=function(){
      'Gene list of 3D genes and transcripts.csv'
  },
  content=function(file){
      write.csv(target2save(),file,row.names = F,quote = F,na = '')
  }
)
# output$download_target_list <- downloadHandler(
#   filename=function(){
#     if(input$list_to_download=='DE genes'){
#       'Gene list of DE genes.csv'
#     }
#     if(input$list_to_download=='DAS genes'){
#       'Gene list of DAS genes.csv'
#     }
#     if(input$list_to_download=='DE transcripts'){
#       'Gene list of DE transcripts.csv'
#     }
#     if(input$list_to_download=='DTU transcripts'){
#       'Gene list of DTU transcripts.csv'
#     }
#   },
#   content=function(file){
#     if(input$list_to_download=='DE genes'){
#       write.csv(data.frame(`DE genes`=unique(DDD.data$DE_genes$target),check.names = F),file,row.names = F)
#     }
#     if(input$list_to_download=='DAS genes'){
#       write.csv(data.frame(`DAS genes`=unique(DDD.data$DAS_genes$target),check.names = F),file,row.names = F)
#     }
#   }
# )


# observeEvent(input$generate.target.list,{
#   write.csv(data.frame(`DE genes`=unique(DDD.data$DE_genes$target)),
#             file = paste0(DDD.data$result.folder,'/DE gene list union set of all contrast group.csv'),
#             row.names = F)
#   write.csv(data.frame(`DAS genes`=unique(DDD.data$DAS_genes$target)),
#             file = paste0(DDD.data$result.folder,'/DAS gene list union set of all contrast group.csv'),
#             row.names = F)
#   write.csv(data.frame(`DE transcripts`=unique(DDD.data$DE_trans$target)),
#             file = paste0(DDD.data$result.folder,'/DE transcript list union set of all contrast group.csv'),
#             row.names = F)
#   write.csv(data.frame(`DTU transcripts`=unique(DDD.data$DTU_trans$target)),
#             file = paste0(DDD.data$result.folder,'/DTU transcript list union set of all contrast group.csv'),
#             row.names = F)
#   message(paste0('Tables of target list save in: ',DDD.data$result.folder))
#   showNotification(paste0('Tables of target list save in: ',DDD.data$result.folder))
# })

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

## update height
observe({
  updateNumericInput(session,inputId = 'go.plot.height',
                     value = round(pmax(10,row.height()*0.4)/2.54,1))
})


observeEvent(input$save.go.plot,{
  if(is.null(data2goplot()))
    return(NULL)
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
  message(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
  showNotification(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
})

observeEvent(input$page_before_function, {
  newtab <- switch(input$tabs, "function" = "ddd","ddd" = "function")
  updateTabItems(session, "tabs", newtab)
})

observeEvent(input$page_after_function, {
  newtab <- switch(input$tabs, "function" = "TSIS","TSIS" = "function")
  updateTabItems(session, "tabs", newtab)
})
