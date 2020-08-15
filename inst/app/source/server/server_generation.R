#=======================>> Data generation panel <<============================
##----------Step 2: How to generate data?------------
##show how to generate data
observeEvent(input$tabs,{
  if(input$tabs=='generation'){
    text2show <- 'Page: Data generation'
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
    showmessage(text = text2show,showNoteify = F)
    showmessage(text = '############################################################',showNoteify = F,showTime = F)
  }
})

output$show_how_to_generate_data <- renderUI({
  # if(!DDD.data$docker_image){
  if(F){
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
                        data genertion steps and visualise/redo the analysis directly from where you saved the data.</strong></font></h4>')
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
      return(showmessage('Please select a file of required format.',type = 'error'))
    suppressWarnings(load(file_path))
    if(!exists('intermediate_data')){
      showmessage('Please select the correct data object: intermediate_data.RData!',type = 'error')
      return(NULL)
    }
    for(i in names(intermediate_data)){
      DDD.data[[i]] <- intermediate_data[[i]]
    }
    DDD.data$folder <-
      DDD.data$data.folder <-
      DDD.data$figure.folder <-
      DDD.data$result.folder <-
      DDD.data$report.folder <- NULL
    incProgress(0.7)
  })
})

observe({
  if(DDD.data$docker_image){
    ##make temp folder
    tempfolder <- tempfile(tmpdir=Sys.getenv("TMP"))
    while(file.exists(tempfolder))
      tempfolder <- tempfile(tmpdir=Sys.getenv("TMP"))
    tempfolder <- gsub('\\\\','/',tempfolder)

    ##for JHI docker
    #=================================================
    #start from here
    # DDD.data$upload_folder <- paste0('/tmp',tempfolder)
    # dir.create(file.path('/tmp', tempfolder))
    #end to here
    #=================================================

    ##for local
    #=================================================
    #start from here
    if(grepl('Win',Sys.info()[1])){
      DDD.data$upload_folder <- tempfolder
      dir.create(file.path(tempfolder))
    } else {
      DDD.data$upload_folder <- paste0('/tmp',tempfolder)
      dir.create(file.path('/tmp', tempfolder))
    }
    #end to here
    #=================================================

    # DDD.data$folder <- getwd()##for test
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
  ##create log

  createLog(logfile = 'log.txt',path = DDD.data$folder)

  # start.time.3D <- Sys.time()
  message(paste0('3D RNA-seq log information (',Sys.time(),')'))
  message('##########################################################\n')
  message('Many thanks for using our 3D RNA-seq shiny App o^_^o!')
  message('If you have questions to raise or are experiencing difficulties using the 3D RNA-seq, please use:','\n',
      'Github Issues: https://github.com/wyguo/ThreeDRNAseq/issues','\n',
      'Group group: https://groups.google.com/forum/#!forum/3d-rna-seq-app-user-group','\n')
  message('Or contact us:','\n',
      'Dr. Wenbin Guo: wenbin.guo@hutton.ac.uk','\n',
      'Dr. Runxuan Zhang: runxuan.zhang@hutton.ac.uk','\n')


  cat('Outputs are saved in:','\n')
  cat(DDD.data$upload_folder,'\n')
  cat(DDD.data$data.folder,'\n')
  cat(DDD.data$figure.folder,'\n')
  cat(DDD.data$result.folder,'\n')
  cat(DDD.data$report.folder,'\n')
  cat('\n')
})

##----------Step 3: Inputs of 3D analysis------------
####mapping
observeEvent(input$mapping_file_button,{
  file_path <- input$mapping_file_button$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  showmessage('Loading transcript-gene mapping...',
              action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
              duration = NULL,id = 'message_id')
  withProgress(message = 'Loading transcript-gene mapping...', detail = 'This may take a while...',
               value = 0, {
                 incProgress(0)
                 if(input$mapping_file_type=='gtf'){
                   if(!grepl('.gtf',file_path))
                     return(showmessage('Please select a file of required format.',type = 'error'))
                   incProgress(0.2)
                   mapping <- gtf2mapping(gtfile = file_path)
                   incProgress(0.7)
                 } else if(input$mapping_file_type=='fa'){
                   if(!grepl('.fa|.fasta',file_path))
                     return(showmessage('Please select a file of required format.',type = 'error'))
                   incProgress(0.2)
                   mapping <- fa2mapping(pafile = file_path)
                   incProgress(0.7)
                 } else {
                   incProgress(0.2)
                   if(!grepl('.csv',file_path))
                     return(showmessage('Please select a file of required format.',type = 'error'))
                   mapping <- read.csv(file = file_path,header = T,fileEncoding="UTF-8-BOM")
                   if(ncol(mapping)!=2)
                     return(showmessage('Please provide transcript-gene association table in correct format.',type = 'error'))
                   incProgress(0.7)
                 }
                 colnames(mapping) <- c('TXNAME','GENEID')
                 rownames(mapping) <- mapping$TXNAME
                 DDD.data$mapping <- mapping
                 incProgress(1)
               })
  showmessage('Loading transcript-gene mapping -> Done!!!',
              action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
              duration = NULL,showNoteify = F)
  cat('\n')
  removeNotification(id = 'message_id')
})

####sample
observe({
  file_path <- input$sample_file_button$datapath
  if (is.null(file_path) | length(file_path)==0)
    return(NULL)
  withProgress(message = 'Loading metadata table...', value = 0, {
    if(!grepl('.csv',file_path))
      return(showmessage('Please select a file of required format.',type = 'error'))
    showmessage('Loading metadata table ...',
                action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
                duration = NULL,id = 'message_id')
    samples0 <- read.csv(file_path,header = T,fileEncoding="UTF-8-BOM")
    DDD.data$samples0 <- samples0
    incProgress(1)
    showmessage('Loading metadata table -> Done!!!',
                action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
                duration = NULL,showNoteify = F)
    cat('\n')
    removeNotification(id = 'message_id')
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
      selectInput(inputId = 'quant_platform',label = 'with',width = '50%',
                  choices = c('Command-line','Galaxy'),
                  selected = 'Command-line',multiple = F)
      %>% helper(icon = "question-circle",
                 colour = NULL,
                 type = 'markdown',
                 title = 'Prepare quantification input files for 3D RNA-seq App',
                 size = 'l',
                 content = 'galaxy',
                 style="font-size: 2.0rem;margin-right:50px;"),
      HTML('<strong> > Select the zipped quantification file:</strong> '),
      fileInput("quant_zip_button", "", accept = c(".zip",".tar",".tgz")),
      HTML('<strong>Note:</strong> Please zip the quantification folder before upload.
           <ul style="padding-left: 20px">
           <li>It is recommended to compress the quantification folder to "zip" format.
              If it is in other formats, special symbols in the compressed folder name may cause errors in the unzipping process.</li>
           <li>If transcripts are quantified by using Salmon command line, 3D RNA-seq App will read "quant.sf" quantification files from zipped folder.</li>
           <li>If transcripts are quantified by using Kallisto command line, 3D RNA-seq App will read "abundance.tsv" quantification files (size is smaller than ".h5" file) from zipped folder.</li>
           <li>If transcripts are quantified by using Salmon/Kallisto on Galaxy, 3D RNA-seq App will read "xxx.tabular" quantification files with ".tabular" extension from zipped folder for both tools.</li>
           <li>Click the above question mark for more details.</li></ul>')
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

##---------> Bookmark: zipped file -----------
observe({
  if(DDD.data$docker_image){
    file_path <- input$quant_zip_button$datapath
    if (is.null(file_path) | length(file_path)==0)
      return(NULL)

    # if(!grepl('.zip',file_path)){
    #   showmessage('Please select a zipped file of required format "xxx.zip".')
    #   return(NULL)
    # }

    ##unzip function
    pos <- regexpr("\\.([[:alnum:]]+)$", file_path)
    zipp <- ifelse(pos > -1L, substring(file_path, pos + 1L), "")
    if(zipp == 'zip'){
      unzipfun <- unzip
      fileNames0 <- unzip(file_path,list = T)$Name
    }


    if(zipp %in% c('gz' ,'tar','tgz')){
      unzipfun <- untar
      fileNames0 <- untar(file_path,list = T)
    }

    if(!grepl(dirname(fileNames0)[1],file_path)){
      showmessage('The compressed quantification file name or its sub-file names include special symbols. Please uncompressed it and then compressed again to "zip" format.',
                  duration = NULL,type = 'error')
      return(NULL)
    }


    ###--->quantification file names

    if(input$quant_platform=='Command-line'){
      if(input$quant_method=='salmon'){
        idx <- grep(pattern = '/quant.sf$',fileNames0)
        if(length(idx)==0){
          showmessage('Please double check whether: (1) quant.sf files are in the zipped file and (2) you have selected the correct quantification method in Step 1: (3).',
                      duration = NULL,type = 'error')
          return(NULL)
        }
        fileNames0 <- fileNames0[idx]
        name.idx <- basename(dirname(fileNames0))
        names(fileNames0) <- name.idx
      } else {
        idx <- grep(pattern = '/abundance.tsv$',fileNames0)
        if(length(idx)==0){
          showmessage('Please double check whether: (1) abundance.tsv files are in the zipped file and (2) you have selected the correct quantification method in Step 1: (3).',
                      duration = NULL,type = 'error')
          return(NULL)
        }
        fileNames0 <- fileNames0[idx]
        name.idx <- basename(dirname(fileNames0))
        names(fileNames0) <- name.idx
      }
    } else {
      idx <- grep(pattern = '[.]tabular',fileNames0)
      if(length(idx)==0){
        showmessage('Please double check whether: (1) .tabular files are in the zipped file and (2) you have selected the correct quantification method in Step 1: (3).',
                    duration = NULL,type = 'error')
        return(NULL)
      }
      fileNames0 <- fileNames0[idx]
      name.idx <- gsub('[.]tabular','',basename(fileNames0))
      names(fileNames0) <- name.idx
    }
    withProgress(message = 'Unzipping transcript quantification file...', value = 0, {
      incProgress(0.1)
      incProgress(0.5)

      startmessage(text = 'Unzip quantification zipped file')
      DDD.data$quant_fileNames <- fileNames0
      unzipfun(file_path, list = F, exdir = DDD.data$upload_folder)
      incProgress(1)
      endmessage(text = 'Unzip quantification zipped file')
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

    if(input$quant_platform=='Command-line'){
      if(input$quant_method=='salmon'){
        idx <- grep(pattern = '/quant.sf$',fileNames0)
        if(length(idx)==0){
          showmessage('Please double check whether: (1) quant.sf files are in the quantification file and (2) you have selected the correct quantification method in Step 1: (3).',
                      duration = NULL,type = 'error')
          return(NULL)
        }
        fileNames0 <- fileNames0[idx]
        name.idx <- basename(dirname(fileNames0))
        names(fileNames0) <- name.idx
      } else {
        idx <- grep(pattern = '/abundance.tsv$',fileNames0)
        if(length(idx)==0){
          showmessage('Please double check whether: (1) abundance.tsv files are in the quantification file and (2) you have selected the correct quantification method in Step 1: (3).',
                      duration = NULL,type = 'error')
          return(NULL)
        }
        fileNames0 <- fileNames0[idx]
        name.idx <- basename(dirname(fileNames0))
        names(fileNames0) <- name.idx
      }
    } else {
      idx <- grep(pattern = '[.]tabular$',fileNames0)
      if(length(idx)==0){
        showmessage('Please double check whether: (1) .tabular files are in the quantification file and (2) you have selected the correct quantification method in Step 1: (3).',
                    duration = NULL,type = 'error')
        return(NULL)
      }
      fileNames0 <- fileNames0[idx]
      name.idx <- gsub('[.]tabular$','',basename(fileNames0))
      names(fileNames0) <- name.idx
    }
    DDD.data$quant_fileNames <- fileNames0
  }

})


output$show_srep_column <- renderUI({
  if(input$has_srep=='No')
    return(NULL)
  selectInput(inputId = 'srep_column',label = '>Select sequencing (technical) replicate column',
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
  paste0('Outputs are saved in:\n', DDD.data$upload_folder)
})

# observeEvent(input$samples_update,{
observeEvent(input$samples_update,{
  if(any(is.null(input$factor_column)) |
     any(is.null(input$brep_column)) |
     any(is.null(input$quant_column)) |
     is.null(DDD.data$upload_folder) |
     is.null(DDD.data$samples0) )
    return(NULL)

  if(is.null(DDD.data$quant_fileNames)){
    showmessage('The quantification file names in the sample meta-table do not match to the quantification files.
                Please double check whether (1) the quantification file names in the metadata table match to the files
                in the quantification folder or (2) you have selected the correct quantification method in Step 1: (3).',
                duration = NULL,type = 'error')
    return(NULL)
  }

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

  ##check quant files names
  quantfile2check <- gsub('.tabular','',DDD.data$samples0[,input$quant_column])
  idx <- grepl(pattern = paste0(quantfile2check,collapse = '|'),DDD.data$quant_fileNames)
  if(!all(idx==T)){
    showmessage('The quantification file names in the sample meta-table do not match to the quantification files.
                Please double check whether (1) the quantification file names in the metadata table match to the files
            in the quantification folder or (2) you have selected the correct quantification method in Step 1: (3).',
                duration = NULL,type = 'error')
    return(NULL)
  }

  quant_fileNames <- file.path(DDD.data$upload_folder,DDD.data$quant_fileNames[quantfile2check])
  if(any(file.exists(quant_fileNames)==F)){
    showmessage('Transcript quantification files are not found. Please double check: (1) whether you have uploaded the quantification compressed folder and
                (2) whether the compressed folder name includes special symbols, e.g. extra white space.',
                duration = NULL,type = 'error')
    return(NULL)
  }

  ## check if select the correct tool
  x2check <- as.character(read.table(quant_fileNames[1],nrows = 1)[1])

  if(input$quant_method=='salmon' & x2check!='Name'){
    showmessage('The quantification file does not match to Salmon output format. Please double check whetehr you have selected the correct quantification tool.',
                duration = NULL,type = 'error')
    return(NULL)
  }

  if(input$quant_method=='kallisto' & x2check!='target_id'){
    showmessage('The quantification file does not match to Kallisto output format.
                Please double check whetehr you have selected the correct quantification tool.',
                duration = NULL,type = 'error')
    return(NULL)
  }

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

  ##check sequencing replicates
  condition2check <- paste0(DDD.data$samples_new$condition)
  idx2check <- table(condition2check)
  if(any(idx2check==1)){
    text.idx <-paste0(paste0(names(idx2check)[idx2check==1][1:3],collapse = ', '),', ... ')
    showmessage(paste0('Condition ',text.idx,' have only one biological replicate. Please double the experimental design.'),
                duration = NULL,id = 'message_id')
    DDD.data$samples_new <- NULL
    return(NULL)
  }

  breps <- DDD.data$samples_new[,input$brep_column]
  idx2check <- by(data = breps,INDICES = condition2check,FUN = function(x){
    table(x)
  })
  name.idx <- names(idx2check)
  idx2check <- as.vector(unlist(idx2check))
  names(idx2check) <- name.idx

  if(any(idx2check > 1 )){
    text.idx <-paste0(paste0(names(idx2check)[idx2check>1][1:3],collapse = ', '),', ... ')
    showmessage(paste0('Duplicated biological replicate labels: ',text.idx,'Please double check (1) the replicate labels in your metadata table or (2) whether sequencing or technical replcates were selected.'),
                duration = NULL,id = 'message_id')
    DDD.data$samples_new <- NULL
    return(NULL)
  }
  removeNotification(id = 'message_id')
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
                duration = NULL,type = 'error')
    return(NULL)
  }

  if(is.null(DDD.data$samples)){
    showmessage(text = 'Please select factor labels in Step 2 and CLICK the button to add selected information for the analysis.',
                duration = NULL,type = 'error')
    return(NULL)
  }

  if(is.null(DDD.data$mapping)){
    showmessage(text = 'Please select or generate correct transcript-gene mapping table in Step 1 for the analysis',
                duration = NULL,type = 'error')
    return(NULL)
  }

  data.files <-  as.character(DDD.data$samples$quant.path)
  if(!all(file.exists(data.files))){
    showmessage('The quantification file names in the sample meta-table do not match to the quantification files.
                Please double check the sample information csv spreadsheet or the selected directory of the quantification files.',
                duration = NULL,type = 'error')
    return(NULL)
  }

  sample.name <- as.character(t(DDD.data$samples$sample.name))

  showmessage('Generating expression, please wait until this step done ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-coffee'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
                   duration = NULL,id = 'generate_expression_message')

  ###check if the transcript names in quantification file matches to the mapping table
  file2check <- data.files[1]
  showmessage('Verifying transcript names...',
                   action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
                   duration = NULL,id = 'message_id')

  if(grepl('.h5',file2check) & DDD.data$docker_image){
    showmessage('Please use the kallisto quantification files "abundance.tsv" as input files for 3D Anlaysis.',type = 'error',
                duration = 20)
    return(NULL)
  }

  txi_trans2check<- tximport(files = file2check,
                             type = input$quant_method,
                             tx2gene = NULL,
                             countsFromAbundance = input$tximport_method,
                             txOut = T)


  trans2check <- rownames(txi_trans2check$abundance)
  idx2keep <- intersect(trans2check,DDD.data$mapping$TXNAME)
  showmessage('Verifying transcript names -> Done!!!',
              action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
  cat('\n')
  removeNotification(id = 'message_id')

  if(length(idx2keep)==0){
    removeNotification(id = 'generate_expression_message')
    showmessage('The transcript names in transctipt-gene association table do not
                match to the transcript names in the quantification files.',type = 'error',
                duration = 20)
    return(NULL)
  }

  DDD.data$mapping <- DDD.data$mapping[idx2keep,]

  showmessage('Generate gene expression...',
              action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
              duration = NULL,id = 'message_id')

  withProgress(message = 'Generate gene expression...', value = 0, {
    incProgress(0.1)
    txi_genes <- tximport(data.files,
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

  showmessage('Generate gene expression -> Done!!!',
              action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
  cat('\n')
  removeNotification(id = 'message_id')


  showmessage('Generate transcript expression...',
              action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
              duration = NULL,id = 'message_id')
  withProgress(message = 'Generate transcript expression...', value = 0, {
    incProgress(0.1)
    txi_trans<- tximport(data.files,
                         type = input$quant_method,
                         tx2gene = NULL,
                         countsFromAbundance = input$tximport_method,
                         txOut = T)

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
  })
  showmessage('Generate transcript expression -> Done!!!',
              action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
  cat('\n')
  removeNotification(id = 'message_id')

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

