sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = '*.R')) {
    #if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    #if(trace) cat("/n")
  }
}
sourceDir(path = 'R',encoding = 'UTF-8')
# library(ThreeDRNAseq)
#######################################################################################################
#######################################################################################################
# data.size.max <- 500
cat('##########################################################\n')
cat('Many thanks for using our 3D RNA-seq shiny App o^_^o!','\n')
cat('If you have questions to raise or are experiencing difficulties using the 3D RNA-seq, please use:','\n',
    'Github Issues: https://github.com/wyguo/ThreeDRNAseq/issues','\n',
    'Group group: https://groups.google.com/forum/#!forum/3d-rna-seq-app-user-group','\n')
cat('Or contact us:','\n',
    'Dr. Wenbin Guo: wenbin.guo@hutton.ac.uk','\n',
    'Dr. Runxuan Zhang: runxuan.zhang@hutton.ac.uk','\n')



## app.R ##
library(shiny)
library(shinydashboard)
library(rhandsontable)
library(shinyFiles)
library(shinyjs)
library(shinyBS)
library(shinyhelper)
library(shinyWidgets)
library(magrittr)
library(DT)
library(plotly)
library(ggplot2)
library(eulerr)
library(gridExtra)
library(grid)
library(fastcluster)
library(rmarkdown)

library(tximport)
library(edgeR)
library(limma)
library(RUVSeq)
library(ComplexHeatmap)
###new added package
# library(parallel)
library(base64enc)
library(ggrepel)

# library(TSIS)

options(stringsAsFactors=F)
# setwd("D:/PhD project/R projects/test round 2019/ThreeDRNAseq")

##dashboardHeader######################################################
# ========================== dashboardHeader ======================== #
mainheader <- dashboardHeader(
  title = '3D RNA-seq App',
  titleWidth = 190,
  dropdownMenuOutput("messageMenu")
)

##dashboardSidebar#####################################################
# ========================= dashboardSidebar ======================== #
mainsidebar <- dashboardSidebar(
  width = 190,
  sidebarMenu(id = "tabs",
    menuItem(text = 'Introduction',tabName = 'introduction',icon = icon("user"),selected = T),
    menuItem(text = 'Data generation',tabName = 'generation',icon = icon("table"),selected = F),
    menuItem(text = 'Data pre-processing',tabName = 'preprocessing',icon = icon("filter"),startExpanded = F,selected = F),
    menuItem(text = '3D analysis',tabName = 'ddd',icon = icon("cubes",lib='font-awesome'),selected = F),
    menuItem(text = 'Time-series trend',tabName = 'tstrend',icon = icon("chart-line",lib='font-awesome'),selected = F),
    menuItem(text = 'Functional plot',tabName = 'function',icon = icon("cogs",lib='font-awesome'),selected = F),
    menuItem(text = 'Isoform switch',tabName = 'TSIS',icon = icon("random"),selected = F),
    menuItem(text = 'Generate report',tabName = 'report',icon = icon("file",lib='font-awesome'),selected = F)
    # menuItem(text = 'Contact us',tabName = 'contact',icon = icon("comment",lib='font-awesome'),selected = F)
  )
)

##dashboardBody########################################################
# ========================== dashboardBody ========================== #
###timeout if no mouse click
inactivity <- "function idleTimer() {
var t = setTimeout(logout, 3600000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 3600000);  // time is in milliseconds (10000 is 10 second; 3600000 is a hour)
}
}
idleTimer();"

mainbody <- dashboardBody(
  # tags$head(includeScript(system.file("google-analytics.js", package = "ThreeDRNAseq"))),
  # tags$head(includeScript("google-analytics.js")),
  tags$script(src = "google-analytics.js"),
  # tags$div(includeScript("https://raw.githubusercontent.com/wyguo/ThreeDRNAseq/master/inst/google-analytics.js")),
  
  # tags$head(tags$style(HTML('.box {margin-left: -15px;}'))),
  withMathJax(),
  shinyjs::useShinyjs(),
  # tags$script(inactivity),### not close browser
  tabItems(
    #=======================>> Introduction panel <<==========================
    source(file.path("source/ui", "ui_introduction.R"),local = TRUE)$value,
    
    #=======================>> Data generation panel <<==========================
    source(file.path("source/ui", "ui_generation.R"),local = TRUE)$value,
    
    #=======================>> Data pre-processing <<============================
    source(file.path("source/ui", "ui_preprocessing.R"),local = TRUE)$value,
    
    #=======================>> DE/DAS/DTU <<============================
    source(file.path("source/ui", "ui_ddd.R"),local = TRUE)$value,
    
    #=======================>> TS trend analysis <<============================
    source(file.path("source/ui", "ui_tstrend.R"),local = TRUE)$value,
    
    
    #=======================>> Functional interpretation <<============================
    source(file.path("source/ui", "ui_function.R"),local = TRUE)$value,
    
    #=======================>> Functional interpretation <<============================
    source(file.path("source/ui", "ui_tsis.R"),local = TRUE)$value,
    
    #=======================>> Generate report <<============================
    source(file.path("source/ui", "ui_report.R"),local = TRUE)$value
    
    #=======================>> Contact us <<============================
    # source(file.path("ui", "ui_contact.R"),local = TRUE)$value
  )
)



ui <- dashboardPage(
  mainheader,
  mainsidebar,
  mainbody
)


server <- function(input, output, session) {
  
  # set data size
  options(shiny.maxRequestSize = 10*1024^3,shiny.sanitize.errors = TRUE)
  options(DT.options = list(scrollX = TRUE,scrollCollapse=TRUE,autoWidth=F,
                            columnDefs = list(list(className = 'dt-left', targets="_all"))))
  
  observe_helpers(withMathJax = TRUE)
  DDD.data <- reactiveValues(
    docker_image=T,
    run_linux=grepl("unix|linux",.Platform$OS.type),
    num_cores=ceiling(parallel::detectCores()/4),
    dataClass='intermediate_data',
    folder=getwd(),
    ##added by GS-->
    upload_folder=getwd(),
    quant_folder=getwd(),
    quant_fileNames = NULL,
    ##------------->
    data.folder=NULL,
    result.folder=NULL,
    figure.folder=NULL,
    report.folder=NULL,
    mapping=NULL,
    samples0=NULL,
    samples=NULL,
    samples_new=NULL,
    brep_column=NULL,
    srep_column=NULL,
    block=NULL,
    factor_column=NULL,
    conditions=NULL,
    trans_counts=NULL,
    genes_counts=NULL,
    trans_TPM=NULL,
    target_high=NULL,
    trans_batch=NULL,
    genes_batch=NULL,
    contrast0=NULL,
    contrast=NULL,
    contrast_pw=NULL,
    contrast_mean=NULL,
    contrast_pgdiff=NULL,
    PS=NULL,
    deltaPS=NULL,
    genes_log2FC=NULL,
    trans_log2FC=NULL,
    DE_genes=NULL,
    genes_3D_stat=NULL,
    trans_3D_stat=NULL,
    DAS_genes=NULL,
    DE_trans=NULL,
    DTU_trans=NULL,
    trans.ann=NULL,
    genes.ann=NULL,
    scores=NULL,
    scores_filtered=NULL,
    RNAseq_info=NULL,
    DDD_numbers=NULL,
    DEvsDAS_results=NULL,
    DEvsDTU_results=NULL,
    params_list=NULL,
    params_table=NULL,
    ##tstrend
    ##tstrend
    DAS_pval_method_tstrend=NULL,
    whether_tstrend=NULL,
    tstrend_spline=NULL,
    tstrend_spline_df=NULL,
    tstrend_Contrastgroups=NULL,
    tstrend_timePoints=NULL,
    tstrend_timeGroups=NULL,
    genes_3D_tstrend_stat=NULL,
    trans_3D_tstrend_stat=NULL,
    DE_tstrend_genes=NULL,
    DAS_tstrend_genes=NULL,
    DE_tstrend_trans=NULL,
    DTU_tstrend_trans=NULL,
    DDD_tstrend_design=NULL,
    DDD_tstrend_numbers=NULL
  )
  # 
 
  #=======================>> Introduction <<============================
  observeEvent(input$page_after_introduction, {
    newtab <- switch(input$tabs, "generation" = "introduction","introduction" = "generation")
    updateTabItems(session, "tabs", newtab)
    shinyjs::runjs("window.scrollTo(0, 50)")
  })
  #Important information
  observe({
    if(DDD.data$docker_image){
      showModal(modalDialog(
        title = HTML("<font color='red'><strong>Important message </strong></font><i class='fa fa-bell'></i>"),
        HTML("<div align='justify'>Many thanks for using our 3D RNA-seq App for your RNA-seq data analysis.
             <ul>
             <li>Due to our server capacity, it may take a while for the web browser to response when multi-users are running the analysis at the same time. Please be patience and only click the button once and wait until one step done. </li>
             <li>The App will be disconnected from our server after a hour if there is no mouse action on the web browser.</li>
             <li>We have fixed most of the bugs in the 3D RNA-seq App. But disconnection from server may also happen due to unexpected bugs or server maintaince. In such case, please double check your input data format or re-try later.</li> 
             <li>If you have questions to raise or are experiencing difficulties using the 3D RNA-seq, please use the 3D RNA-seq user group: <a href='https://groups.google.com/forum/#!forum/3d-rna-seq-app-user-group' target='_blank'>https://groups.google.com/forum/#!forum/3d-rna-seq-app-user-group</a>.</li> 
             <li>The 3D RNA-seq App can also be run through RStudio on a local computer. Please go to the Github page for details: <a href='https://github.com/wyguo/ThreeDRNAseq' target='_blank'>https://github.com/wyguo/ThreeDRNAseq</a>. </li>
             </ul></div>"),
        HTML('<div align="justify"><font color="red"><strong>Update: 11/12/2019 </strong></font>
             <ul>
             <li>Fix the bug in the "Time-series trend" page.</li>
             <li>Improve some parts of the GUI for easy to use. </li> 
             </ul>
             Update history can be found in Github Wiki page: <a href="https://github.com/wyguo/ThreeDRNAseq/wiki" target="_blank">
             https://github.com/wyguo/ThreeDRNAseq/wiki</a>
             </div>')
      ))
    }
  })
  
  #=======================>> Data generation panel <<============================
  source(file.path("source/server", "server_generation.R"),local = TRUE)$value
  
  #=======================>> Data pre-processing <<==============================
  source(file.path("source/server", "server_preprocessing.R"),local = TRUE)$value
  
  
  #=======================>> 3D analysis <<============================
  source(file.path("source/server", "server_ddd.R"),local = TRUE)$value
  
  
  #=======================>> TS trend analysis <<============================
  source(file.path("source/server", "server_tstrend.R"),local = TRUE)$value
  
  #=======================>> Functional interpretation <<============================
  source(file.path("source/server", "server_function.R"),local = TRUE)$value
  
  #=======================>> Functional interpretation <<============================
  source(file.path("source/server", "server_tsis.R"),local = TRUE)$value
  
  #=======================>> Generate report <<============================
  source(file.path("source/server", "server_report.R"),local = TRUE)$value
  
  #=======================>> Contact us <<============================
  # observeEvent(input$page_before_contact, {
  #   newtab <- switch(input$tabs, "report" = "contact","contact" = "report")
  #   updateTabItems(session, "tabs", newtab)
  # })
  
  closeAllConnections()
  session$allowReconnect(TRUE)
  
  
  ##stop the session
  # session$onSessionEnded(function() {
  #   stopApp()
  # })
  
}
shinyApp(ui, server,options=list(launch.browser=T))
