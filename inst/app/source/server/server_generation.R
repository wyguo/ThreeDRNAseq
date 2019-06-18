#=======================>> Data generation panel <<============================
##----------Step 2: How to generate data?------------
##show how to generate data
output$show_how_to_generate_data <- renderUI({
  # if(!DDD.data$docker_image){
  if(!DDD.data$docker_image){
    tagList(
      fluidRow(
        ##----------Step 2: How to generate data of analysis?------------
        box(title = 'How to generate data of analysis?',
            width=6,status = 'primary', solidHeader = T,
            column(12,
                   # wellPanel(
                   # style = "background-color: #ffffff;",
                   radioButtons(inputId = 'generate_new_data',label = 'How to generate data?',
                                choices = c('Generate new data','Upload intermediate data'),selected = 'Generate new data',inline = T),
                   uiOutput('upload_data_ui'),
                   HTML('<h5 align="justify">If you have already run the 3D RNA-seq analysis once and the intermediate data has been saved as
                        "intermediate_data.RData" object, which can be uploaded here. <font color="red"><strong>By doing this, you can skip the following
                        data genertion steps and visualise/redo the analysis directly from "Data pre-processing" until to "Generate report" step.</strong></font></h4>')
                   # )
                   )
            )
        )
        )
  }
})

output$upload_data_ui <- renderUI({
  if(input$generate_new_data=='Generate new data')
    return(NULL)
  fileInput("upload_data", "Select intermediate_data.RData",
            accept = c(
              ".RData")
  )
})

observe({
  file_path <- input$upload_data$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  
  withProgress(message = 'Loading intermediate_data.RData...', value = 0, {
    if(!grepl('.RData',file_path))
      return(showmessage('Please select a file of required format.'))
    suppressWarnings(load(file_path))
    if(!exists('intermediate_data')){
      showmessage('Please select the correct data object: intermediate_data.RData!')
      return(NULL)
    }
    for(i in names(intermediate_data)){
      DDD.data[[i]] <- intermediate_data[[i]]
    }
    incProgress(0.7)
  })
})

observe({
  if(DDD.data$docker_image){
    ##make temp folder
    tempfolder <- tempfile(tmpdir=Sys.getenv("TMP"))
    while(file.exists(tempfolder))
      tempfolder <- tempfile(tmpdir=Sys.getenv("TMP"))
    
    ##for JHI docker
    DDD.data$upload_folder <- paste0('/tmp',tempfolder)
    dir.create(file.path('/tmp', tempfolder))
    
    ##for local
    # dir.create(file.path(tempfolder))
    # DDD.data$upload_folder <- tempfolder
    
    DDD.data$folder <- DDD.data$upload_folder
    DDD.data$data.folder <- file.path(DDD.data$folder,'data')
    DDD.data$figure.folder <- file.path(DDD.data$folder,'figure')
    DDD.data$result.folder <- file.path(DDD.data$folder,'result')
    DDD.data$report.folder <- file.path(DDD.data$folder,'report')
  } else {
    DDD.data$folder <- getwd()
    DDD.data$data.folder <- file.path(DDD.data$folder,'data')
    DDD.data$figure.folder <- file.path(DDD.data$folder,'figure')
    DDD.data$result.folder <- file.path(DDD.data$folder,'result')
    DDD.data$report.folder <- file.path(DDD.data$folder,'report')
  }
  create.folders(wd = DDD.data$folder)
  message(DDD.data$upload_folder)
  message(DDD.data$data.folder)
  message(DDD.data$figure.folder)
  message(DDD.data$result.folder)
  message(DDD.data$report.folder)
})

##----------Step 3: Inputs of 3D analysis------------
####mapping 
observeEvent(input$mapping_file_button,{
  file_path <- input$mapping_file_button$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  
  withProgress(message = 'Loading transcript-gene mapping...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0)
                 if(input$mapping_file_type=='gtf'){
                   if(!grepl('.gtf',file_path))
                     return(showmessage('Please select a file of required format.'))
                   incProgress(0.2)
                   mapping <- gtf2mapping(gtfile = file_path)
                   incProgress(0.7)
                 } else if(input$mapping_file_type=='fa'){
                   if(!grepl('.fa',file_path))
                     return(showmessage('Please select a file of required format.'))
                   incProgress(0.2)
                   mapping <- fa2mapping(pafile = file_path)
                   incProgress(0.7)
                 } else {
                   incProgress(0.2)
                   if(!grepl('.csv',file_path))
                     return(showmessage('Please select a file of required format.'))
                   mapping <- read.csv(file = file_path,header = T)
                   if(ncol(mapping)!=2)
                     return(showmessage('Please provide transcript-gene association table in correct format.'))
                   incProgress(0.7)
                 }
                 colnames(mapping) <- c('TXNAME','GENEID')
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
    if(!grepl('.csv',file_path))
      return(showmessage('Please select a file of required format.'))
    samples0 <- read.csv(file_path,header = T)
    DDD.data$samples0 <- samples0
    incProgress(1)
  })
})

####quant zip file
output$show_quant_load <- renderUI({
  if(DDD.data$docker_image){
    tagList(
      HTML('<strong> (3) Select transcript quantification "zip" file.</strong> '),
      p(),
      radioButtons(inputId = 'quant_method',label = '>Transcripts are quantified by:',
                   choices = c('salmon','kallisto'),
                   selected = 'salmon',inline = T),
      HTML('<strong> > Select the zipped quantification file:</strong> '),
      fileInput("quant_zip_button", "", accept = c(".zip")),
      HTML('<strong>Note:</strong> Please zip the quantification folder before upload.')
    )
  } else {
    tagList(
      HTML('<strong> (3) Select transcript quantification folder from local directory.</strong> '),
      p(),
      radioButtons(inputId = 'quant_method',label = '>Transcripts are quantified by:',
                   choices = c('salmon','kallisto'),
                   selected = 'salmon',inline = T),
      HTML('<strong> > Select the quantification folder:</strong> '),
      p(),
      shinyDirButton(id = 'quant_folder_button',label = 'Select a folder',title = 'Select the folder of quantification',
                     buttonType = 'primary',
                     icon = icon('folder'))
    )
  }
})

## quantification directory 
volumes_quant <- getVolumes()
shinyDirChoose(input = input,id = 'quant_folder_button',roots=volumes_quant,
               session = session,restrictions = system.file(package = "base"))

observe({
  if(DDD.data$docker_image){
    file_path <- input$quant_zip_button$datapath
    if (is.null(file_path) | length(file_path)==0)
      return(NULL)
    showNotification('Uploading quantification zipped file, please wait until this step done ...',
                     # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                     action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                     duration = NULL,id = 'upload_zip_quant_message')
    
    withProgress(message = 'Unzipping transcript quantification file...', value = 0, {
      incProgress(0.1)
      incProgress(0.5)
      if(!grepl('.zip',file_path))
        return(showmessage('Please select a zipped file of required format "xxx.zip".'))
      ###--->quantification file names
      fileNames0 <- unzip(file_path,list = T)$Name
      
      if(input$quant_method=='salmon'){
        fileNames <- fileNames0[grep('.sf',fileNames0)]
      } else {
        fileNames <- fileNames0[grep('.h5',fileNames0)]
        if(length(fileNames)==0)
          fileNames <- fileNames0[grep('.tsv',fileNames0)]
      }
      
      if(length(fileNames)==0){
        showmessage('Please check (1) whether the quantification folder includes transcript quantification files ("quant.sf" for samlom and "abundance.h5/abundance.tsv" for kallisto).(2) whether you have selected the correct quantification method "salmon/kallisto".',
                    duration = 20,id = 'check_filenames_message')
        removeNotification(id = 'upload_zip_quant_message')
        return(NULL)
      }
      
      fileNames.idx <- basename(dirname(fileNames))
      names(fileNames) <- fileNames.idx
      DDD.data$quant_fileNames <- fileNames
      unzip(file_path, list = F, exdir = DDD.data$upload_folder)
      
      removeNotification(id = 'check_filenames_message')
      
      incProgress(1)
      removeNotification(id = 'upload_zip_quant_message')
      showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
      # cat(DDD.data$upload_folder)
    })
  } else {
    if(is.null(input$quant_folder_button))
      return(NULL)
    upload_folder <- parseDirPath(roots = volumes_quant,selection = input$quant_folder_button)
    if(identical(upload_folder,character(0)))
      return(NULL)
    DDD.data$upload_folder <- upload_folder 
    fileNames0 <- list.files(upload_folder,recursive = T)
    
    if(input$quant_method=='salmon'){
      fileNames <- fileNames0[grep('.sf',fileNames0)]
    } else {
      fileNames <- fileNames0[grep('.h5',fileNames0)]
      if(length(fileNames)==0)
        fileNames <- fileNames0[grep('.tsv',fileNames0)]
    }
    
    if(length(fileNames)==0){
      showmessage('Please check (1) whether the quantification folder includes transcript quantification files ("quant.sf" for samlom and "abundance.h5/abundance.tsv" for kallisto).(2) whether you have selected the correct quantification method "salmon/kallisto".',
                  duration = 20,id = 'check_filenames_message')
      return(NULL)
    }
    
    fileNames.idx <- basename(dirname(fileNames))
    names(fileNames) <- fileNames.idx
    DDD.data$quant_fileNames <- fileNames
  }
})


output$show_srep_column <- renderUI({
  if(input$has_srep=='No')
    return(NULL)
  selectInput(inputId = 'srep_column',label = '>Select sequencing replicate column',
              choices = colnames(DDD.data$samples0),
              selected = DDD.data$srep_column,
              multiple = F,width = NULL)
})



observe({
  if(is.null(DDD.data$samples0))
    return(NULL)
  updateSelectInput(session = session,inputId = 'factor_column',
                    choices = colnames(DDD.data$samples0),
                    selected = DDD.data$factor_column)
  updateSelectInput(session = session,inputId = 'brep_column',
                    choices = colnames(DDD.data$samples0),
                    selected = DDD.data$brep_column)
  updateSelectInput(session = session,inputId = 'quant_column',
                    choices = colnames(DDD.data$samples0),
                    selected = DDD.data$quant_column)
})


output$qaunt_folder_path_text <- renderText({
  if(is.null(DDD.data$upload_folder))
    return(NULL)
  paste0('Transcript quantification in directory:\n', DDD.data$upload_folder)
})

# observeEvent(input$samples_update,{
observeEvent(input$samples_update,{
  if(any(is.null(input$factor_column)) | 
     any(is.null(input$brep_column)) | 
     any(is.null(input$quant_column)) | 
     is.null(DDD.data$upload_folder) |
     is.null(DDD.data$samples0))
    return(NULL)
  
  DDD.data$brep_column <- input$brep_column
  DDD.data$srep_column <- input$srep_column
  DDD.data$factor_column <- input$factor_column
  DDD.data$quant_column <- input$quant_column
  
  ##make.name
  make.name.idx <- c(input$brep_column,input$srep_column,input$factor_column)
  make.name.idx <- make.name.idx[!is.null(make.name.idx)]
  for(i in make.name.idx){
    DDD.data$samples0[,i] <- make.names(DDD.data$samples0[,i])
  }
  
  
  quantfile2check <- DDD.data$samples0[,input$quant_column]
  if(all(quantfile2check %in% names(DDD.data$quant_fileNames))){
    quant_fileNames <- paste0(DDD.data$upload_folder,'/',DDD.data$quant_fileNames[quantfile2check])
  } else {
    showmessage('The quantification file names in the sample meta-table do not match to the quantification files. 
                Please double check the sample information csv spreadsheet or the selected directory of the quantification files.',
                duration = 10)
    return(NULL)
  }
  
  
  # ##new add
  # if(input$quant_method=='kallisto'){
  #   idx <- sapply(idx, function(i){
  #     if(!file.exists(i))
  #       gsub('abundance.h5','abundance.tsv',i)
  #     else i
  #   })
  #   idx <- as.vector(idx)
  # }
  
  if(input$has_srep=='Yes'){
    samples <- data.frame(DDD.data$samples0,
                          quant.path=quant_fileNames)
    samples$sample.name <- as.vector(interaction(DDD.data$samples0[,c(input$factor_column,input$brep_column,input$srep_column)]))
    samples$condition <- make.names(as.vector(interaction(DDD.data$samples0[,c(input$factor_column)])))
    DDD.data$samples <- samples
    idx2keep <- as.vector(interaction(DDD.data$samples0[,c(input$factor_column,input$brep_column)]))
    idx2keep <- which(!duplicated(idx2keep))
    DDD.data$samples_new <- samples[idx2keep,]
  } else {
    samples <- data.frame(DDD.data$samples0,
                          quant.path=quant_fileNames)
    samples$sample.name <- as.vector(interaction(DDD.data$samples0[,c(input$factor_column,input$brep_column)]))
    samples$condition <- make.names(interaction(DDD.data$samples0[,c(input$factor_column)]))
    DDD.data$samples_new <- DDD.data$samples <- samples
  }
  save(samples,file='samples.RData')
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
##update tximport method
observe({
  if(is.null(DDD.data$params_list$tximport_method))
    return(NULL)
  updateSelectInput(inputId = "tximport_method",session = session,
                    selected = DDD.data$params_list$tximport_method
  )
})

##update has srep
observe({
  if(is.null(DDD.data$params_list$has_srep))
    return(NULL)
  updateRadioButtons(inputId = 'has_srep',session = session,label = '>Does the data have sequencing replicates?',
                     choices = c('No','Yes'),selected = DDD.data$params_list$has_srep,inline = T)
})

##update quant method
observe({
  if(is.null(DDD.data$params_list$quant_method))
    return(NULL)
  updateRadioButtons(inputId = 'quant_method',session = session,label = '>Transcripts are quantified by:',
                     choices = c('salmon','kallisto'),
                     selected = DDD.data$params_list$quant_method,
                     inline = T)
})


observeEvent(input$tximport_run,{
  # if(is.null(DDD.data$samples) | is.null(DDD.data$mapping) | is.null(DDD.data$folder))
  #   return(NULL)
  path.idx <- parseDirPath(roots = volumes_quant,selection = input$quant_folder_button)
  if(identical(x = path.idx,y = character(0))){
    showmessage(text = 'Please select the quatification folder in Step 2.',
                duration = 20)
    return(NULL)
  }
  
  if(is.null(DDD.data$samples)){
    showmessage(text = 'Please select factor labels in Step 2 and CLICK the button to add selected information for the analysis.',
                duration = 20)
    return(NULL)
  }
  
  if(is.null(DDD.data$mapping)){
    showmessage(text = 'Please select or generate transcript-gene mapping table in Step 1 for the analysis',
                duration = 20)
    return(NULL)
  }
  
  data.files <-  as.character(DDD.data$samples$quant.path)
  if(!all(file.exists(data.files))){
    showmessage('The quantification file names in the sample meta-table do not match to the quantification files. 
                Please double check the sample information csv spreadsheet or the selected directory of the quantification files.',
                duration = 10)
    return(NULL)
  }
  
  sample.name <- as.character(t(DDD.data$samples$sample.name))
  
  showNotification('Generating expression, please wait until this step done ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'generate_expression_message')
  
  ###check if the transcript names in quantification file matches to the mapping table
  showNotification('Verifying transcript names...',
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'message_id')
  file2check <- data.files[1]
  txi_trans2check<- tximport(files = file2check, 
                             type = input$quant_method, 
                             tx2gene = NULL,
                             countsFromAbundance = input$tximport_method,
                             txOut = T,
                             dropInfReps = T)
  trans2check <- rownames(txi_trans2check$abundance)
  idx2keep <- intersect(trans2check,DDD.data$mapping$TXNAME)
  removeNotification(id = 'message_id')
  
  if(length(idx2keep)==0){
    removeNotification(id = 'generate_expression_message')
    showmessage('The transcript names in transctipt-gene association table do not match to the transcript names in the quantification files.',
                duration = 20)
    return(NULL)
  }
  
  DDD.data$mapping <- DDD.data$mapping[idx2keep,]
  
  withProgress(message = 'Generate gene expression...', value = 0, {
    incProgress(0.1)
    txi_genes <- tximport(data.files,dropInfReps = T,
                          type = input$quant_method, 
                          tx2gene = DDD.data$mapping,
                          countsFromAbundance = input$tximport_method)
    
    idx2keep <- intersect(DDD.data$mapping$GENEID,rownames(txi_genes$counts))
    txi_genes$counts <- txi_genes$counts[idx2keep,]
    txi_genes$abundance <- txi_genes$abundance[idx2keep,]
    txi_genes$length <- txi_genes$length[idx2keep,]
    
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
    
    
  })
  
  withProgress(message = 'Generate transcript expression...', value = 0, {
    incProgress(0.1)
    txi_trans<- tximport(data.files, 
                         type = input$quant_method, 
                         tx2gene = NULL,
                         countsFromAbundance = input$tximport_method,
                         txOut = T,
                         dropInfReps = T)
    
    idx2keep <- intersect(DDD.data$mapping$TXNAME,rownames(txi_trans$counts))
    txi_trans$counts <- txi_trans$counts[idx2keep,]
    txi_trans$abundance <- txi_trans$abundance[idx2keep,]
    txi_trans$length <- txi_trans$length[idx2keep,]
    
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
    removeNotification(id = 'generate_expression_message')
    showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
  })
  # showNotification('Read count and TPM generation are done.',
  #                  action = a(href="#shiny-tab-preprocessing","Next page"))
})

observeEvent(input$page_before_generation, {
  newtab <- switch(input$tabs, "generation" = "introduction","introduction" = "generation")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

observeEvent(input$page_after_generation, {
  newtab <- switch(input$tabs, "preprocessing" = "generation","generation" = "preprocessing")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})
