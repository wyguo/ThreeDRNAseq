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
data.size.max <- 500
message('Many thanks for using our 3D RNA-seq shiny App o^_^o!')
message('The step-by-step user manual of this app can be found in:\n https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md')
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
    menuItem(text = '3D analysis',tabName = 'ddd',icon = icon("bar-chart-o",lib='font-awesome'),selected = F),
    menuItem(text = 'Functional plot',tabName = 'function',icon = icon("chart-area",lib='font-awesome'),selected = F),
    menuItem(text = 'Isoform switch',tabName = 'TSIS',icon = icon("random"),selected = F),
    menuItem(text = 'Generate report',tabName = 'report',icon = icon("file",lib='font-awesome'),selected = F)
    # menuItem(text = 'Contact us',tabName = 'contact',icon = icon("comment",lib='font-awesome'),selected = F)
  )
)

##dashboardBody########################################################
# ========================== dashboardBody ========================== #
mainbody <- dashboardBody(
  # tags$head(includeScript(system.file("google-analytics.js", package = "ThreeDRNAseq"))),
  # tags$head(includeScript("google-analytics.js")),
  tags$script(src = "google-analytics.js"),
  # tags$div(includeScript("https://raw.githubusercontent.com/wyguo/ThreeDRNAseq/master/inst/google-analytics.js")),
  
  # tags$head(tags$style(HTML('.box {margin-left: -15px;}'))),
  withMathJax(),
  shinyjs::useShinyjs(),
  tabItems(
    #=======================>> Introduction panel <<==========================
    source(file.path("source/ui", "ui_introduction.R"),local = TRUE)$value,
    
    #=======================>> Data generation panel <<==========================
    source(file.path("source/ui", "ui_generation.R"),local = TRUE)$value,
    
    #=======================>> Data pre-processing <<============================
    source(file.path("source/ui", "ui_preprocessing.R"),local = TRUE)$value,
    
    #=======================>> DE/DAS/DTU <<============================
    source(file.path("source/ui", "ui_ddd.R"),local = TRUE)$value,
    
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
  options(shiny.maxRequestSize = data.size.max*1024^2)
  options(DT.options = list(scrollX = TRUE,scrollCollapse=TRUE,autoWidth=F,
                            columnDefs = list(list(className = 'dt-left', targets="_all"))))
  
  observe_helpers(withMathJax = TRUE)
  DDD.data <- reactiveValues(
    dataClass='intermediate_data',
    folder=getwd(),
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
    params_table=NULL
  )
  # 
 
  #=======================>> Introduction <<============================
  observeEvent(input$page_after_introduction, {
    newtab <- switch(input$tabs, "generation" = "introduction","introduction" = "generation")
    updateTabItems(session, "tabs", newtab)
  })
  
  #=======================>> Data generation panel <<============================
  source(file.path("source/server", "server_generation.R"),local = TRUE)$value
  
  #=======================>> Data pre-processing <<==============================
  source(file.path("source/server", "server_preprocessing.R"),local = TRUE)$value
  
  
  #=======================>> 3D analysis <<============================
  source(file.path("source/server", "server_ddd.R"),local = TRUE)$value
  
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
  
  session$allowReconnect(TRUE)
  
}
shinyApp(ui, server,options=list(launch.browser=T))
