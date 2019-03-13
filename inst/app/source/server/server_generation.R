#=======================>> Data generation panel <<============================


####upload intermediate data is disabled
# observeEvent(input$load_intermediate_data_btn,{
#   file.idx <- 'data/txi_genes.RData'
#   if(!file.exists(file.idx))
#     return(NULL)
#   load(file.idx)
# 
#   file.idx <- 'data/txi_trans.RData'
#   if(!file.exists(file.idx))
#     return(NULL)
#   load(file.idx)
# 
#   file.idx <- 'data/intermediate_data.RData'
#   if(!file.exists(file.idx))
#     return(NULL)
#   load(file.idx)
# 
#   for(i in names(intermediate_data)){
#     DDD.data[[i]] <- intermediate_data[[i]]
#   }
# })
# 
# ##----------Step 1: Select folders to save results------------
# output$data_folder_button_ui <- renderUI({
#   if(input$use_diff_folder=='No')
#     return(NULL)
#   shinyDirButton(id = 'data_folder_button',label = 'Select a folder',
#                  title = 'Select a folder to save intermidate data',
#                  buttonType = 'primary',
#                  icon = icon('folder'))
# })
# 
# ##----------Step 2: How to generate data?------------
# output$upload_data_ui <- renderUI({
#   if(input$generate_new_data=='Generate new data')
#     return(NULL)
#   fileInput("upload_data", "Select intermediate_data.RData",
#             accept = c(
#               ".RData")
#   )
# })
# 
# observe({
#   file_path <- input$upload_data$datapath
#   if (is.null(file_path) | length(file_path)==0)
#     return(NULL)
#   withProgress(message = 'Loading intermediate_data.RData...', value = 0, {
#     suppressWarnings(load(file_path))
#     for(i in names(intermediate_data)){
#       DDD.data[[i]] <- intermediate_data[[i]]
#     }
#     incProgress(0.7)
#   })
# })

# volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
# volumes <- c(wd='D:/PhD project/R projects/test round 2018/ThreeDRNAseq_improve/')
volumes <- c(wd=getwd())
shinyFileChoose(input, "mapping_file_button", roots = volumes, session = session,filetypes='csv')
shinyFileChoose(input, "sample_file_button", roots = volumes, session = session,filetypes='csv')
shinyDirChoose(input = input,id = 'data_folder_button',roots=volumes,session = session,restrictions = system.file(package = "base"))
shinyDirChoose(input = input,id = 'quant_folder_button',roots=volumes,session = session,restrictions = system.file(package = "base"))

# observe({
#   path.idx <- parseDirPath(roots = volumes,selection = input$data_folder_button)
#   if(!identical(x = path.idx,y = character(0)))
#     DDD.data$folder <- path.idx
# })

observe({
  DDD.data$folder <- getwd()
  DDD.data$data.folder <- 'data'
  DDD.data$figure.folder <- 'figure'
  DDD.data$result.folder <- 'result'
  DDD.data$report.folder <- 'report'
})

output$data_folder_path_text <- renderText({
  if(is.null(DDD.data$data.folder))
    return(NULL)
  DDD.data$folder
})

##----------Step 3: Inputs of 3D analysis------------
####mapping 
observe({
  file_path <- input$mapping_file_button$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  withProgress(message = 'Loading transcript-gene mapping...', value = 0, {
  mapping <- read.csv(file_path,header = T)
  rownames(mapping) <- mapping$TXNAME
  DDD.data$mapping <- mapping
  incProgress(1)
  })
})

####sample 
observe({
  file_path <- input$sample_file_button$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  withProgress(message = 'Loading transcript-gene mapping...', value = 0, {
    samples0 <- read.csv(file_path,header = T)
    DDD.data$samples0 <- samples0
    incProgress(1)
  })
})

####quant zip file
observe({
  file_path <- input$quant_zip_button$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  withProgress(message = 'Unzipping transcript quantification file...', value = 0, {
    incProgress(0.1)
    incProgress(0.5)
    unzip(file_path, list = F, exdir = getwd())
    incProgress(1)
    showmessage('Done!!!')
  })
})


output$show_srep_column <- renderUI({
  if(input$has_srep=='No')
    return(NULL)
  selectInput(inputId = 'srep_column',label = '>Select sequencing replicate column',
              choices = colnames(DDD.data$samples0),selected = NULL,multiple = F,width = NULL)
  
})

observe({
  if(is.null(DDD.data$samples0))
    return(NULL)
  updateSelectInput(session = session,inputId = 'factor_column',
                    choices = colnames(DDD.data$samples0),selected = NULL)
  updateSelectInput(session = session,inputId = 'brep_column',
                    choices = colnames(DDD.data$samples0),selected = NULL)
  updateSelectInput(session = session,inputId = 'quant_column',
                    choices = colnames(DDD.data$samples0),selected = NULL)
})

output$qaunt_folder_path_text <- renderText({
  path.idx <- parseDirPath(roots = volumes,selection = input$quant_folder_button)
  if(identical(x = path.idx,y = character(0)))
    return(NULL)
  paste0('Transcript quantification in directory:\n', path.idx)
  
})

observeEvent(input$samples_update,{
  path.idx <- parseDirPath(roots = volumes,selection = input$quant_folder_button)
  if(any(is.null(input$factor_column)) | 
     any(is.null(input$brep_column)) | 
     any(is.null(input$quant_column)) | 
     identical(x = path.idx,y = character(0)))
    return(NULL)
  
  DDD.data$brep_column <- input$brep_column
  DDD.data$srep_column <- input$srep_column
  DDD.data$factor_column <- input$factor_column
  
  idx <- paste0(path.idx,'/',DDD.data$samples0[,input$quant_column],'/',ifelse(input$quant_method=='salmon','quant.sf','abundance.h5'))
  if(input$has_srep=='Yes'){
    samples <- data.frame(DDD.data$samples0,
                          quant.path=idx)
    samples$sample.name <- as.vector(interaction(DDD.data$samples0[,c(input$factor_column,input$brep_column,input$srep_column)]))
    samples$condition <- make.names(as.vector(interaction(DDD.data$samples0[,c(input$factor_column)])))
    DDD.data$samples <- samples
    idx2keep <- as.vector(interaction(DDD.data$samples0[,c(input$factor_column,input$brep_column)]))
    idx2keep <- which(!duplicated(idx2keep))
    DDD.data$samples_new <- samples[idx2keep,]
  } else {
    samples <- data.frame(DDD.data$samples0,
                          quant.path=idx)
    samples$sample.name <- as.vector(interaction(DDD.data$samples0[,c(input$factor_column,input$brep_column)]))
    samples$condition <- make.names(interaction(DDD.data$samples0[,c(input$factor_column)]))
    DDD.data$samples_new <- DDD.data$samples <- samples
  }
})

output$mapping_table <- DT::renderDataTable({
  if(is.null(DDD.data$mapping))
    return(NULL)
  DDD.data$mapping[1:20,]
})

output$sample_table <- DT::renderDataTable({
  if(is.null(DDD.data$samples0))
    return(NULL)
  if(is.null(DDD.data$samples))
    DDD.data$samples0 else DDD.data$samples
})

##----------Step 4: Generate gene and transcript read counts and TPMs------------

observeEvent(input$tximport_run,{
  if(is.null(DDD.data$samples) | is.null(DDD.data$mapping) | is.null(DDD.data$folder))
    return(NULL)
  
  data.files <-  as.character(DDD.data$samples$quant.path)
  if(!all(file.exists(data.files))){
    showmessage('The provided quantification folder information is incorrect. 
  Please double check the sample information csv spreadsheet or the selected directory of the quantification files.',duration = 20)
    return(NULL)
  }
  
  sample.name <- as.character(t(DDD.data$samples$sample.name))
  create.folders(wd = DDD.data$folder)
  
  withProgress(message = 'Generate gene expression...', value = 0, {
    incProgress(0.1)
    txi_genes <- tximport(data.files,dropInfReps = T,
                          type = input$quant_method, 
                          tx2gene = DDD.data$mapping,
                          countsFromAbundance = input$tximport_method)
    colnames(txi_genes$counts) <-
      colnames(txi_genes$abundance) <- 
      colnames(txi_genes$length)<-sample.name
    incProgress(0.7)
    
    ##save all
    save(txi_genes,file=paste0(DDD.data$data.folder,'/txi_genes.RData'))
    
    ##save TPMs
    genes_TPM <- txi_genes$abundance
    write.csv(genes_TPM,file=paste0(DDD.data$result.folder,'/Gene TPM.csv'))
    
    ##save counts
    idx <- as.vector(interaction(DDD.data$samples[,c(input$factor_column,input$brep_column)]))
    genes_counts <- sumarrays(x = txi_genes$counts,group = idx)
    write.csv(genes_counts,file=paste0(DDD.data$result.folder,'/Gene read counts.csv'))
    
    ####
    DDD.data$genes_TPM <- genes_TPM
    DDD.data$genes_counts <- genes_counts
    
    incProgress(1)
    showNotification('Gene read counts and TPMs are generated.')
    message('Gene read counts and TPMs are generated.')
    
  })
  
  withProgress(message = 'Generate transcript expression...', value = 0, {
    incProgress(0.1)
    txi_trans<- tximport(data.files, 
                         type = input$quant_method, 
                         tx2gene = NULL,
                         countsFromAbundance = input$tximport_method,
                         txOut = T,
                         dropInfReps = T)
    colnames(txi_trans$counts) <-
      colnames(txi_trans$abundance) <- 
      colnames(txi_trans$length)<-sample.name
    incProgress(0.7)
    
    ##save all
    save(txi_trans,file=paste0(DDD.data$data.folder,'/txi_trans.RData'))
    
    ##save TPMs
    trans_TPM <- txi_trans$abundance
    write.csv(trans_TPM,file=paste0(DDD.data$result.folder,'/Transcript TPM.csv'))
    
    ##save counts
    idx <- as.vector(interaction(DDD.data$samples[,c(input$factor_column,input$brep_column)]))
    trans_counts <- sumarrays(txi_trans$counts,group = idx)
    write.csv(trans_counts,file=paste0(DDD.data$result.folder,'/Transcript read counts.csv'))
    
    ####
    DDD.data$trans_TPM <- trans_TPM
    DDD.data$trans_counts <- trans_counts
    incProgress(1)
    showNotification('Transcript read counts and TPMs are generated.')
    message('Transcript read counts and TPMs are generated.')
  })
  # showNotification('Read count and TPM generation are done.',
  #                  action = a(href="#shiny-tab-preprocessing","Next page"))
})

observeEvent(input$page_before_generation, {
  newtab <- switch(input$tabs, "generation" = "introduction","introduction" = "generation")
  updateTabItems(session, "tabs", newtab)
})

observeEvent(input$page_after_generation, {
  newtab <- switch(input$tabs, "preprocessing" = "generation","generation" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
})
