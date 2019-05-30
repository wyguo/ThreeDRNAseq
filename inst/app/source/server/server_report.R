###parameters
observeEvent(input$save_ddd_data_button,{
  # if(is.null(DDD.data$params_list)){
  #   DDD.data$params_list <- list()
    DDD.data$params_list$condition_n = length(unique((DDD.data$samples_new$condition)))
    DDD.data$params_list$brep_n = length(unique(DDD.data$samples0[,DDD.data$brep_column]))
    DDD.data$params_list$srep_n = length(unique(DDD.data$samples0[,DDD.data$srep_column]))
    DDD.data$params_list$samples_n = nrow(DDD.data$samples_new)
    DDD.data$params_list$has_srep = input$has_srep
    DDD.data$params_list$quant_method = input$quant_method
    DDD.data$params_list$tximport_method = input$tximport_method
    DDD.data$params_list$cpm_cut = input$cpm_cut
    DDD.data$params_list$cpm_samples_n = input$cpm_samples_n
    DDD.data$params_list$norm_method = input$norm_method
    DDD.data$params_list$has_batcheffect = input$has_batcheffect
    DDD.data$params_list$RUVseq_method = input$RUVseq_method
    DDD.data$params_list$contrast = DDD.data$contrast
    DDD.data$params_list$DE_pipeline = input$DE_pipeline
    DDD.data$params_list$pval_adj_method = input$pval_adj_method
    DDD.data$params_list$pval_cut = input$pval_cut
    DDD.data$params_list$l2fc_cut = input$l2fc_cut
    DDD.data$params_list$deltaPS_cut = input$deltaPS_cut
    DDD.data$params_list$DAS_pval_method = input$DAS_pval_method
    
    ##heatmap
    DDD.data$params_list$dist_method <- input$dist_method
    DDD.data$params_list$cluster_method <- input$cluster_method
    DDD.data$params_list$cluster_number <- input$cluster_number
    
    ##TSIS
    DDD.data$params_list$TSISorisokTSP <- input$TSISorisokTSP
    DDD.data$params_list$TSIS_method_intersection <- input$method_intersection
    DDD.data$params_list$TSIS_spline_df <- input$spline_df
    DDD.data$params_list$TSIS_prob_cut <- input$TSIS_prob_cut
    DDD.data$params_list$TSIS_diff_cut <- input$TSIS_diff_cut
    DDD.data$params_list$TSIS_adj_pval_cut <- input$TSIS_adj_pval_cut
    DDD.data$params_list$TSIS_time_point_cut <- ifelse(input$TSISorisokTSP == 'isokTSP',1,input$TSIS_time_point_cut)
    DDD.data$params_list$TSIS_cor_cut <- input$TSIS_cor_cut
    
    x <- DDD.data$params_list
    x <- lapply(x,function(i){paste0(i,collapse = '; ')})
    x <- data.frame(Description=names(x),Parameter=unlist(x),row.names = NULL)
    x$Description <- gsub('_',' ',x$Description)
    x$Description <- gsub('TSIS','IS',x$Description)
    
    x$Description[x$Description=='IS method intersection'] <- 'Values to identify ISs'
    DDD.data$params_table <- x
  # }
  
})

########parameter tables
# output$params_table_panel <- DT::renderDataTable({
#   x <- DDD.data$params_table
#   rownames(x) <- NULL
#   x
# },editable = T,options = list(pageLength = 50))
# 
# proxy = DT::dataTableProxy('params_table_panel')
# observeEvent(input$params_table_cell_edit, {
#   info = input$params_table_cell_edit
#   i = info$row
#   j = info$col
#   v = info$value
#   # problem starts here
#   DDD.data$params_table[i, j] <- isolate(suppressWarnings(DT::coerceValue(v, x$df[i, j])))
# })


##----------Step 2: Save data and results ------------
##---------->Generate report ------------
observeEvent(input$generate_report,{
  showNotification('Generating report ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'generate_report_message')
  withProgress(message = 'Generating report',
               detail = 'This may take a while...', value = 0, {
                 incProgress(0.3)
                 for (i in c('html_document','word_document','pdf_document')) {
                   tryCatch({
                     rmarkdown::render(input = '3D_report.Rmd',
                                       output_format = i,
                                       output_dir = DDD.data$report.folder,
                                       intermediates_dir = DDD.data$report.folder,
                                       knit_root_dir = DDD.data$report.folder,
                                       params = c(DDD.data=list(isolate(reactiveValuesToList(DDD.data))))
                     )
                   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                 }
                 
                 ####
                 # folder2copy <- paste0(DDD.data$folder,'/www')
                 # if(!file.exists(folder2copy))
                 #   dir.create(folder2copy,recursive = T)
                 # file.copy(from = paste0(DDD.data$folder,'/report'), 
                 #           to = folder2copy, recursive=TRUE,overwrite = T)
                 showmessage('Done!!!')
                 incProgress(1)
               })
  removeNotification(id = 'generate_report_message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})


######
observeEvent(input$save_ddd_data_button,{
  showNotification('Saving results ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'save_result_message')
  ####save intermediate data
  withProgress(message = 'Saving results...',
               detail = 'This may take a while...', value = 0, {
                 incProgress(0)
                 intermediate_data <- isolate(reactiveValuesToList(DDD.data))
                 save(intermediate_data,file=paste0(DDD.data$data.folder,'/intermediate_data.RData'))
                 incProgress(0.3)
                 ####save results 
                 idx <- c('DE_genes','DAS_genes','DE_trans','DTU_trans','samples','contrast',
                          'DDD_numbers','DEvsDAS_results','DEvsDTU_results','RNAseq_info')
                 idx.names <-gsub('_',' ',idx)
                 idx.names <- gsub('trans','transcripts',idx.names)
                 idx.names[1:4] <- paste0('Significant ',idx.names[1:4],' list and statistics')
                 
                 idx <- c(idx,'scores','scores_filtered')
                 idx.names <- c(idx.names,'Raw isoform switch scores','Significant isoform switch scores')
                 incProgress(0.5)
                 for(i in seq_along(idx)){
                   if(is.null(DDD.data[[idx[i]]]))
                     next
                   write.csv(x = DDD.data[[idx[i]]],file = paste0(DDD.data$result.folder,'/',idx.names[i],'.csv'),row.names = F)
                 }
                 ### save transcript-gene mapping
                 write.csv(x = DDD.data$mapping,
                           file = paste0(DDD.data$result.folder,'/Transcript and gene mapping.csv'),
                           row.names = F,na = '')
                 
                 ### save 3d list
                 threeD.list <-lapply(idx[1:4],function(i){
                   unique(DDD.data[[i]]$target)
                 })
                 
                 if(!any(sapply(threeD.list,is.null))){
                   n <- max(sapply(threeD.list, length))
                   threeD.list <- lapply(threeD.list, function(x){
                     y <- rep(NA,n)
                     y[1:length(x)] <- x
                     y
                   })
                   names(threeD.list) <- idx[1:4]
                   threeD <- do.call(cbind,threeD.list)
                   write.csv(x = threeD,file = paste0(DDD.data$result.folder,'/DDD genes and transcript lists across all contrast groups.csv'),
                             row.names = F,na = '')
                 }
                 
                 ##save all gene/transcript statistics
                 incProgress(0.8)
                 write.csv(x = DDD.data$genes_3D_stat$DE.stat,
                           file = paste0(DDD.data$result.folder,'/DE gene testing statistics.csv'),
                           row.names = F,na = '')
                 
                 write.csv(x = DDD.data$trans_3D_stat$DE.stat,
                           file = paste0(DDD.data$result.folder,'/DE transcripts testing statistics.csv'),
                           row.names = F,na = '')
                 
                 if(DDD.data$params_list$DAS_pval_method=='F-test'){
                   write.csv(x = DDD.data$trans_3D_stat$DAS.F.stat,
                             file = paste0(DDD.data$result.folder,'/DAS genes testing statistics.csv'),
                             row.names = F,na = '')
                 }
                 
                 if(DDD.data$params_list$DAS_pval_method=='Simes'){
                   write.csv(x = DDD.data$trans_3D_stat$DAS.simes.stat,
                             file = paste0(DDD.data$result.folder,'/DAS genes testing statistics.csv'),
                             row.names = F,na = '')
                 }
                 
                 write.csv(x = DDD.data$trans_3D_stat$DTU.stat,
                           file = paste0(DDD.data$result.folder,'/DTU transcripts testing statistics.csv'),
                           row.names = F,na = '')
                 incProgress(1)
                 
               })
  removeNotification(id = 'save_result_message')
  showmessage('Done!!!',action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"))
})
# 

observeEvent(input$act_download_zip_results, {
  showNotification('Starting downloading, it may take quite a while ...',
                   # action = HTML("<span style='font-size:50px;'>&#9786;</span> <i style='font-size:25px;' class='fas fa-dizzy'></i>"),
                   action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
                   duration = NULL,id = 'start_download_message')
  shinyjs::runjs("document.getElementById('download_zip_results').click();")
  removeNotification(id = 'start_download_message')
  showmessage('Download is started, please wait ...',
              action = HTML("<i style='font-size:35px;' class='fas fa-smile-wink'></i>"),
              duration = NULL)
})

output$download_zip_results <- downloadHandler(
  filename=function(){
    '3D output.zip'
  },
  content=function(file){
    zip(zipfile = file,files = c(DDD.data$data.folder,DDD.data$figure.folder,
                                 DDD.data$result.folder,DDD.data$report.folder))
  },contentType = "application/zip"
)


observeEvent(input$page_before_report, {
  newtab <- switch(input$tabs, "TSIS" = "report","report" = "TSIS")
  updateTabItems(session, "tabs", newtab)
  shinyjs::runjs("window.scrollTo(0, 50)")
})

# observeEvent(input$page_after_report, {
#   newtab <- switch(input$tabs, "report" = "contact","contact" = "report")
#   updateTabItems(session, "tabs", newtab)
# })
