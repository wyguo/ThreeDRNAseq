#' The shiny app for DE, DAS and DTU analysis
#'
#' @param data.size.max maximum size limit for unload files in Shiny. Default is 300 (MB).
#'
#' @return the Shiny App.
#'
#' @examples RNAseq3D.app()
#'
#' @seealso \code{\link{shiny}}
#'
#' @export

ThreeDRNAseq.app <- function(data.size.max=300) {
  message('Many thanks for using our 3DRNAseq shiny app o^_^o.')
  message('The step-by-step user manual of this app can be found in: xxxx ')
  ## app.R ##
  library(shiny)
  library(shinydashboard)
  library(rhandsontable)
  library(shinyFiles)
  library(tximport)
  library(edgeR)
  library(limma)
  library(plotly)
  library(RUVSeq)
  library(eulerr)
  library(gridExtra)
  library(grid)
  library(ComplexHeatmap)
  library(Gmisc)
  options(stringsAsFactors=F)

  ##dashboardHeader######################################################
  # ========================== dashboardHeader ======================== #
  mainheader <- dashboardHeader(
    title = 'DDD-GUI',
    titleWidth = 200,
    dropdownMenuOutput("messageMenu")
  )

  ##dashboardSidebar#####################################################
  # ========================= dashboardSidebar ======================== #
  mainsidebar <- dashboardSidebar(
    width = 200,
    sidebarMenu(
      menuItem(text = 'Introduction',tabName = 'introduction',icon = icon("book"),selected = F),
      menuItem(text = 'Data generation',tabName = 'generation',icon = icon("database"),selected = F),
      menuItem(text = 'Data pre-processing',tabName = 'preprocessing',icon = icon("cogs"),startExpanded = F,selected = F
      ),
      menuItem(text = 'DE DAS and DTU',tabName = 'ddd',icon = icon("random"),selected = F),
      menuItem(text = 'Result summary',tabName = 'resultssummary',icon = icon("sitemap",lib='font-awesome'),selected = F),
      menuItem(text = 'Advanced plot',tabName = 'advancedanalysis',icon = icon("line-chart",lib='font-awesome'),selected = F),
      menuItem(text = 'Generate report',tabName = 'report',icon = icon("file",lib='font-awesome'),selected = T)
    )
  )

  ##dashboardBody########################################################
  # ========================== dashboardBody ========================== #
  mainbody <- dashboardBody(
    tags$head(includeScript(system.file("google-analytics.js", package = "RNAseq3D"))),
    withMathJax(),
    shinyjs::useShinyjs(),
    tabItems(
      tabItem("introduction",
              fluidRow(
                column(width = 12,
                       box(title='Introduction',
                           width = NULL,status = 'primary', solidHeader = T,
                           # shiny::uiOutput('page1'),
                           # DT::dataTableOutput("info")
                           HTML('<h5 align="justify">The DDD-GUI is used to perform differential expression (DE) gene and transcript, differental
                                alternaive spliced (DAS) gene and differential transcript usage (DTU) transcript analysis.
                                To use our pipeline in your work, please cite:<br><br>
                                Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A.,
                                Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R.,
                                and Brown,J.W.S. (2018) Rapid and dynamic alternative
                                splicing impacts the Arabidopsis cold response transcriptome.
                                Plant Cell. <a href="http://www.plantcell.org/content/30/7/1424" target="_blank"> http://www.plantcell.org/content/30/7/1424 </a></h5>')
                           )
                           )
                           )

                           ),
      ###########################################################################
      #======================= Data generation panel ============================
      tabItem("generation",
              tags$head(tags$style("#TxtOut {white-space: normal;}")),
              fluidRow(
                column(width = 12,
                       box(title = 'Choose working directory',
                           width=13,status = 'primary', solidHeader = T,
                           # HTML('<h4><strong>Choose working directory</strong></h4>'),
                           # wellPanel(
                           column(width = 2,
                                  shinyDirButton("wdir", "Chose directory", "Upload",
                                                 buttonType ="primary")
                           ),
                           column(width = 10,
                                  verbatimTextOutput("wdir.info")
                           )
                       )
                )
              ),
              fluidRow(
                column(width = 12,
                       box(width = 13,title = 'Generate/Load data',
                           status = 'primary',solidHeader = T,
                           # shiny::dataTableOutput('loaded.data.info')
                           verbatimTextOutput("DDD.data.size"),
                           tableOutput('loaded.data.info'),
                           h4('Load data for analysis'),
                           div(style="display:inline-block;vertical-align:middle",
                               actionButton('loaded.data.all','Load all',icon("upload",lib="font-awesome"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                               h4('OR Load data for individual steps'),
                               actionButton('loaded.data.preprocessing','Load Data pre-processing',icon("upload",lib="font-awesome"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                               actionButton('loaded.data.ddd','Load DE DAS and DTU',icon("upload",lib="font-awesome"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                               actionButton('loaded.data.summary','Load Result summary',icon("upload",lib="font-awesome"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                               actionButton('loaded.data.advanced','Load Advanced plot',icon("upload",lib="font-awesome"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           ),
                           hr(),
                           HTML('Notes: The whole analysis has a number of steps and each step will generate several datasets for downstream analysis.
                                These datasets will be saved as ".RData" format to "data" folder in the workding directory and R users can access
                                them by R command lines. This table summarises the steps to generate and use these intermediate datasets, and whether
                                the datasets are already saved in the "data" folder. To generate/upload the datasets, the following options are provided:
                                <ul>
                                <li> If the datasets (in .RData format) do not exist in the "data" folder, users can generate these datasets step-by-step. </li>
                                <li> If all the required datasets are generated and saved in the "data" folder, users can upload by clicking "Load all" button. </li>
                                <li> To perform analysis or visualise results in a specific step, users also can unload the required datasets for individual steps. </li>
                                </ul>')
                           ))
                           ),
              fluidRow(
                #load sample tables####
                column(width = 6,
                       box(title='Step 1: Load quantification and experimental information',
                           width = 13,status = 'primary', solidHeader = T,
                           # HTML('<h4><strong>Step 1: load gene-transcript mapping</strong></h4>'),
                           fileInput("mapping.input", "Choose mapping.csv file",
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                           ),
                           br(),
                           HTML('This file includes the gene-transcript mapping information. See examples in the manual.'),
                           hr(),
                           fileInput("sample.input", "Choose samples.csv file",
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                           ),
                           br(),
                           HTML('This file includes the information of conditions for linear regression,
                                biological replicates, sequencing replicates (if exist) and the complete paths of
                                <a href="https://combine-lab.github.io/salmon/" target="_blank">salmon</a> quantification files "quant.sf".
                                See examples in the manual.')
                           )
                           ),
                #generate gene expression####
                # #generate transcript expression####
                column(width = 6,
                       box(title='Step 2: Generate gene and transcript expression',
                           width = 13,status = 'primary', solidHeader = T,
                           HTML('<h4><strong>Generate gene expression</strong></h4>'),
                           selectInput(inputId = "tximport.quant.method",
                                       label = "Abundance generated from:",
                                       choices = c("salmon","kallisto","sailfish","rsem","stringtie","none"),
                                       selected = "salmon",
                                       width = 200

                           ),
                           div(style="display:inline-block;vertical-align:middle",
                               selectInput(inputId = "tximport.method",
                                           label = "Choose a tximport method:",
                                           choices = c("lengthScaledTPM", "scaledTPM", "no"),
                                           selected = "lengthScaledTPM",
                                           width = 200
                               )
                           ),
                           div(style="display:inline-block;vertical-align:middle;height:30px",
                               actionButton('run.txi.genes','Run',icon("send outline icon"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           ),
                           radioButtons('generate.new.genes.data',label = NULL,inline = T,
                                        choices = c('Generate new data','Load txi_genes.RData'),
                                        selected = 'Generate new data'
                           ),
                           HTML('If the <i>txi_genes.RData</i> file is already in the "data" folder, user can select "Load txi_genes.RData"
                                and click "Load" to load the data directly.'
                           ),
                           br(),
                           verbatimTextOutput("rungeneinfo"),
                           hr(),
                           HTML('<h4><strong>Generate transcript expression</strong></h4>'),
                           div(style="display:inline-block;vertical-align:middle;height:30px",
                               actionButton('run.txi.trans','Run',icon("send outline icon"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4"
                               )
                           ),
                           br(),
                           br(),
                           radioButtons('generate.new.trans.data',label = NULL,inline = T,
                                        choices = c('Generate new data','Load txi_trans.RData'),
                                        selected = 'Generate new data'

                           ),
                           HTML('If the <i>txi_trans.RData</i> file is already in the "data" folder, user can select "Load txi_trans.RData"
                                and click "Run" to load the data directly.'
                           ),
                           br(),
                           verbatimTextOutput("runtransinfo")
                           )
                       )
                           ),
              # output sample table####
              fluidRow(
                tabBox(width = 13,
                       # Title can include an icon
                       # title = tagList(shiny::icon("gear"), "Loaded tables"),
                       title = "Loaded tables",side = 'right',
                       tabPanel("Mapping",
                                "Part of gene-transcript mapping table:",
                                DT::dataTableOutput("mapping.output")
                       ),
                       tabPanel("Samples",
                                "The sample information:",
                                DT::dataTableOutput("sample.output")
                       )
                )
              )
                       ),
      ###########################################################################
      #======================== Data pre-processing panel =======================
      tabItem('preprocessing',
              ##----------Step 1: merge sequencing reps------------
              fluidRow(
                # column(width = 3,
                box(title = 'Step 1: Merge sequencing replicates',
                    width=6,status = 'primary', solidHeader = T,
                    verbatimTextOutput("seqrep.info"),
                    radioButtons(inputId = 'merge.or.skip',label = 'Merge options',
                                 choices = c('Merge','Skip'),inline = T),
                    actionButton('run.seq.merge','Merge',icon("plus-square"),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                    hr(),
                    HTML('If there are no sequencing replicates, please select "Skip".')
                ),
                box(title = 'Step 2: Filter low expression',
                    width = 6,status = 'primary',solidHeader = T,
                    sliderInput(inputId = 'cpm.cut',label = 'CPM cut-off',
                                min = 0,max = 10,value = 1),
                    sliderInput(inputId = 'sample.n.cut',label = 'Sample number cut-off',
                                min = 0,max = 10,value = 1),
                    actionButton('run.filter','Filter',icon("filter"),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                    actionButton('mv.plot.button','Mean-variance plot',icon("pencil"),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                    br(),
                    br(),
                    verbatimTextOutput("low.filter.info"),
                    hr(),
                    HTML('Note: Click "Filter" and then click "Mean-variance plot".
                         <ul><li>An expressed transcript must have &#x2265 <i>n</i> samples &#x2265 <i>m</i> CPM (Count Per Million reads) expression.</li>
                         <li>An expressed gene must have at least one expressed transcript.</li></ul>')
                    )
                    ),
              ##----------Step 2: Filter low expression------------
              fluidRow(
                box(title = NULL,
                    width=12,status = 'primary', solidHeader = T,
                    HTML('RNA-seq data information: target numbers'),
                    DT::dataTableOutput("data.info"),
                    # tableOutput("data.info"),
                    br(),
                    actionButton(inputId = 'save.data.info',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                ),
                ##---------- Mean-variance trend plot------------
                column(width = 12,
                       tabBox(title = 'Mean-variance trend plot',id = 'mv.plot', width = 13,
                              tabPanel(title = 'Transcript level',
                                       plotOutput(outputId = 'mv.trans.plot',click ='mv.trans.plot.click')
                              ),
                              tabPanel(title = 'Gene level',
                                       plotOutput(outputId = 'mv.genes.plot',click ='mv.genes.plot.click')
                              )
                       ),
                       actionButton(inputId = 'save.mv.plot',label = 'Save',
                                    icon = icon('download',lib = 'font-awesome'),
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                )
              ),
              br(),
              br(),
              br(),
              ##----------Step 3: Principal component analysis------------
              fluidRow(
                column(width = 3,
                       box(title = 'Step 3: Principal component analysis (PCA)',width = 13,
                           status = 'primary',solidHeader = T,
                           div(style="display:inline-block;vertical-align:middle;width:150px",
                               selectInput(inputId = 'dim1',label = 'X-axis',
                                           choices = c('PC1','PC2','PC3','PC4'),selected = 'PC1')
                           ),
                           div(style="display:inline-block;vertical-align:middle;width:150px",
                               selectInput(inputId = 'dim2',label = 'Y-axis',
                                           choices = c('PC1','PC2','PC3','PC4'), selected = 'PC2')
                           ),
                           div(style="display:inline-block;vertical-align:middle;height:30px",
                               actionButton('plot.pca.button','Plot',icon("pencil",lib = 'font-awesome'),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           ),
                           br(),
                           br(),
                           radioButtons(inputId = 'pca.plot.type',label = 'Plot type',
                                        choices = c('Bio-reps','Average expression'),
                                        selected = 'Bio-reps',inline = T),
                           radioButtons(inputId = 'pca.color.option',label = 'Colour on',
                                        choices = c('Bio-reps','Conditions'),
                                        selected = 'Bio-reps',inline = T),
                           radioButtons(inputId = 'pca.add.circle',label = 'Add polygon to clusters',
                                        choices = c('none','ellipse','polygon'),
                                        selected = 'none',inline = T),
                           HTML('<ol>
                                <li style="margin-left: -20px;">Select PCs to visualise</li>
                                <li style="margin-left: -20px;">Select "Plot type", either show all Bio-reps or average expression of conditions</li>
                                <li style="margin-left: -20px;">Select colour the Bio-reps or the Condtions</li>
                                <li style="margin-left: -20px;">Select "ellipse" or "polygon" to highlight the clusters</li>
                                </ol>')
                           )
                           ),
                column(width = 9,
                       tabBox(title = 'PCA plot',id = 'pca.plot',width = 13,
                              tabPanel(title = 'Transcript level',
                                       plotlyOutput('trans.pca.plot')
                              ),
                              tabPanel(title = 'Gene level',
                                       plotlyOutput('genes.pca.plot')
                              )
                       ),
                       actionButton(inputId = 'save.pca.plot',label = 'Save',
                                    icon = icon('download',lib = 'font-awesome'),
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                       br(),
                       br(),
                       br(),
                       br()
                )
                           ),
              fluidRow(column(width = 3,
                              box(title = 'Step 3-continued: Batch effect estimation',width = 13,
                                  status = 'primary',solidHeader = T,
                                  radioButtons(inputId = 'ruvseq.method',label = 'RUVSeq method',
                                               choices = c('RUVr','RUVs','RUVg'),selected = 'RUVr',
                                               inline = T),
                                  actionButton('remove.batch','Remove batch effects',icon("cut"),
                                               style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                                  br(),
                                  br(),
                                  actionButton('update.pca.plot','Update PCA plot',icon("refresh",lib="font-awesome"),
                                               style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                                  br(),
                                  br(),
                                  verbatimTextOutput("br.saved.info"),
                                  br(),
                                  HTML('The parameters to plot PCA are inherited from previous PCA panel.
                                       If no distinct batch effects between bio-reps, <strong>please skip this step</strong>.')
                                  )
                              ),
                       column(width = 9,
                              tabBox(title = 'PCA plot after remove batch effects',id = 'pca.plot.batch.removed',width = 13,
                                     tabPanel(title = 'Transcript level',
                                              plotlyOutput('trans.pca.br.plot')
                                     ),
                                     tabPanel(title = 'Gene level',
                                              plotlyOutput('genes.pca.br.plot')
                                     )
                              ),
                              actionButton(inputId = 'save.pca.br.plot',label = 'Save',
                                           icon = icon('download',lib = 'font-awesome'),
                                           style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                              br(),
                              br(),
                              br(),
                              br()
                       )
              ),

              ##----------Step 4: Normalisation------------
              fluidRow(
                column(width = 3,
                       box(title = 'Step 4: Normalisation',width = 13,
                           status = 'primary',solidHeader = T,
                           radioButtons(inputId = 'norm.method',label = 'Normalisation methods',
                                        choices = c('TMM','RLE','upperquartile'),
                                        selected = 'TMM',inline = T),
                           actionButton('run.norm','Run',icon("send outline icon",lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                           actionButton('norm.plot.button','Plot distribution',icon("pencil"),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                           br(),
                           br(),
                           verbatimTextOutput("norm.info")
                       )

                ),
                column(width = 9,
                       tabBox(title = 'Data distribution plot',id = 'dist.plot',width = 13,
                              tabPanel(title = 'Transcript level',
                                       plotlyOutput('trans.dist.plot',height = 850)
                              ),
                              tabPanel(title = 'Gene level',
                                       plotlyOutput('genes.dist.plot',height = 850)
                              )
                       ),
                       actionButton(inputId = 'save.dist.plot',label = 'Save',
                                    icon = icon('download',lib = 'font-awesome'),
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                       br(),
                       br(),
                       br(),
                       br()
                )
              )
              ),
      ##----------DE, DAS and DTU------------
      tabItem('ddd',
              ##-contrast group------------
              fluidRow(
                column(width = 4,
                       box(title='Step 1: Load contrast group',
                           width = 13,status = 'primary', solidHeader = T,
                           # HTML('<h4><strong>Step 1: load gene-transcript mapping</strong></h4>'),
                           fileInput("contrast.input", "Choose contrast.csv file",
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                           ),
                           hr(),
                           HTML('<h4><strong>Contrast groups</strong></h4>'),
                           DT::dataTableOutput('contrast.table'),
                           br(),
                           HTML('The expression of "Treatment" conditions is compared to "Control" conditions.')
                       )
                ),
                column(width = 4,
                       box(title='Step 2: DE genes',
                           width = 13,status = 'primary', solidHeader = T,
                           HTML('<h4><strong>Choose a pipeline</strong></h4>'),
                           div(style="display:inline-block;vertical-align:middle;height:30px",
                               selectInput(inputId = "DE.pipeline", label = NULL,
                                           choices = list(`limma-voom:`=list('limma'),`edgeR:` = c("glmQL", "glm")),
                                           selected = 'limma',width = 160
                               )
                           ),
                           div(style="display:inline-block;vertical-align:middle;height:30px",
                               actionButton('run.DE','Run',icon("send outline icon"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           ),
                           br(),
                           br(),
                           HTML('<p align="justify">In practics, limma-voom and edgeR-glmQL have comparable performance, and both have
                                better control of false positives than edgeR-glm.</p>'),
                           hr(),
                           HTML('<h4><strong>Set thresholds</strong></h4>'),
                           selectInput(inputId = 'p.adjust.method',label = 'P-value adjust method',
                                       choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                   "fdr", "none"),selected = 'BH'),
                           sliderInput(inputId = 'pval.cutoff',label = 'Adjusted p-value',
                                       min = 0,max = 0.05,value = 0.01,step = 0.005),
                           sliderInput(inputId = 'lfc.cutoff',label = 'Absolute \\(\\log_2\\)FC',
                                       min = 0,max = 5,value = 1,step = 0.5),
                           verbatimTextOutput("run.DE.info")
                           )
                ),
                column(width = 4,
                       box(title='Step 3: DAS genes, DE and DTU transcripts',
                           width = 13,status = 'primary', solidHeader = T,
                           HTML('<h4><strong>Generate transcript deltaPS</strong></h4>'),
                           div(style="display:inline-block;vertical-align:middle",
                               selectInput(inputId = "PS.data.type",
                                           label = "Data for deltaPS:",
                                           choices = c("TPM", "Read counts"),
                                           selected = "TPM",
                                           width = 160
                               )
                           ),
                           div(style="display:inline-block;vertical-align:middle;height:30px",
                               actionButton('run.deltaPS','Run',icon("send outline icon"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           ),
                           br(),
                           verbatimTextOutput("rundeltaPSinfo"),
                           hr(),
                           HTML('<h4><strong>Set threshods</strong></h4>'),
                           selectInput(inputId = 'DAS.p.method',label = 'DAS gene test method',
                                       choices = c('F-test','Simes'),selected = 'F-test'),
                           sliderInput(inputId = 'deltaPS.cutoff',label = 'Absolute \\(\\Delta\\)PS',
                                       min = 0,max = 1,value = 0.1,step = 0.05),
                           br(),
                           HTML('<p align="justify">PS (percent of splice) is defined as the ratio of transcript average abundance of conditions divided by the average
                                gene abundance. \\(\\Delta\\)PS is the PS differences of conditions based on the contrast groups. The abundance to calculate
                                PS values can be TPMs or read counts.</p>'),
                           hr(),
                           HTML('<h4><strong>Proceed the analysis</strong></h4>'),
                           div(style="display:inline-block;vertical-align:middle;height:45px",
                               verbatimTextOutput("diffsplice.method")
                           ),
                           div(style="display:inline-block;vertical-align:middle;height:40px",
                               actionButton('run.diffsplice','Run',icon("send outline icon"),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           ),
                           br(),
                           br(),
                           verbatimTextOutput("run.DAS.info")
                           )
                       )
              ),

              ##-DE DAS DTU transcripts------------
              fluidRow(
                column(width = 12,
                       box(title='Part of deltaPS',
                           width = 13,status = 'primary', solidHeader = T,
                           DT::dataTableOutput('deltaPS.table'),
                           actionButton(inputId = 'save.deltaPS.table',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                       )
                )
              )
      ),
      ##----------Results summary------------
      tabItem('resultssummary',
              fluidRow(
                column(width = 3,
                       box(title='DDD target tables',
                           width = 13,status = 'primary', solidHeader = T,
                           HTML('<h4><strong>Visulise statistics</strong></h4>'),
                           selectInput(inputId = 'stat.type',label = 'Targets',
                                       choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts'),
                                       selected = 'DE genes'
                           ),
                           uiOutput('show.contrast.groups'),
                           sliderInput(inputId = 'top.stat',label = 'Top ranked statistics in each contrast',
                                       min = 10,max = 100,value = 10,step = 10),
                           hr(),
                           actionButton(inputId = 'save.stat',label = 'Save all tables',icon = icon('download',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                           br(),
                           br(),
                           HTML('<p align="justify">Up to top 100 ranked targets can be visulise in this App. To visulise the statistics of
                                full gene/transcript sets, please click the "Save all tables" button. The full tables of statistics and target numbers in this page will be saved
                                as csv files in the "result" folder.</p>')
                           )
                           ),
                column(width = 9,
                       uiOutput('show.top.stat.table'),
                       box(title='DE DAS DTU target numbers',
                           width = 13,status = 'primary', solidHeader = T,
                           tableOutput('DDD.numbers')
                       ),
                       fluidRow(
                         column(width = 6,
                                box(title='DE vs DAS genes',
                                    width = 13,status = 'primary', solidHeader = T,
                                    tableOutput('DE.vs.DAS')
                                )
                         ),
                         column(width = 6,
                                box(title='DE vs DTU transcripts',
                                    width = 13,status = 'primary', solidHeader = T,
                                    tableOutput('DE.vs.DTU')
                                )
                         )
                       )
                )
                       ),
              fluidRow(
                column(width = 12,
                       tabBox(width = 13,

                              # Title can include an icon
                              # title = tagList(shiny::icon("gear"), "Loaded tables"),
                              title = "Up-down regulation",side = 'right',
                              tabPanel("DE genes",
                                       plotOutput("up.down.bar.plot.DE.genes")
                              ),
                              tabPanel("DAS genes",
                                       plotOutput("up.down.bar.plot.DAS.genes")
                              ),
                              tabPanel("DE transcripts",
                                       plotOutput("up.down.bar.plot.DE.trans")
                              ),
                              tabPanel("DTU transcripts",
                                       plotOutput("up.down.bar.plot.DTU.trans")
                              )
                       ),
                       actionButton(inputId = 'save.up.down.bar.plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                       br(),
                       br(),
                       br(),
                       br()
                ),
                column(width = 12,
                       box(title='Comparison across individual contrast groups',
                           width = 13,status = 'primary', solidHeader = T,
                           uiOutput('across.contrast.euler'),
                           splitLayout(cellWidths = c("50%", "50%"),
                                       plotOutput("DE.genes.euler"),
                                       plotOutput("DAS.genes.euler")
                           ),
                           br(),
                           splitLayout(cellWidths = c("50%", "50%"),
                                       plotOutput("DE.trans.euler"),
                                       plotOutput("DTU.trans.euler")
                           ),
                           actionButton(inputId = 'save.across.contrast.euler.plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                       ),
                       box(title='Comparison across targets in individual conatrast groups',
                           width = 13,status = 'primary', solidHeader = T,
                           uiOutput('across.target.euler'),
                           splitLayout(cellWidths = c("50%", "50%"),
                                       plotOutput("DEvsDAS.euler"),
                                       plotOutput("DEvsDTU.euler")
                           ),
                           actionButton(inputId = 'save.across.target.euler.plot',label = 'Save',
                                        icon = icon('download',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                       ),
                       box(title='Union set of all contrast groups',
                           width = 13,status = 'primary', solidHeader = T,
                           plotOutput("genes.flow.chart"),
                           hr(),
                           plotOutput("trans.flow.chart"),
                           actionButton(inputId = 'save.union.flow.chart',label = 'Save',
                                        icon = icon('download',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                       )
                )
              )
      ),
      ##----------Advanced plot------------
      tabItem('advancedanalysis',
              ##----------Heatmap------------
              fluidRow(
                # column(width = 12,
                box(title='Heatmap',
                    width = 3,status = 'primary', solidHeader = T,

                    HTML('<h4><strong>Choose targets</strong></h4>'),
                    radioButtons(inputId = 'heatmap.select.or.upload',label = '',
                                 choices = c('Select targets','Upload target list'),inline = T),
                    selectInput(inputId = 'heatmap.target.type',label = 'Select targets',
                                choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts')
                    ),
                    hr(),
                    fileInput("heatmap.target.list.input", "OR Choose target list csv file",
                              accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")
                    ),
                    radioButtons(inputId = 'heatmap.targetlist.type',label = 'Target list type',inline = T,
                                 choices = c('Gene list','Transcript list')
                    ),
                    hr(),
                    HTML('<h4><strong>Clustering</strong></h4>'),
                    selectInput(inputId = 'dist.method',label = 'Distance measurement',
                                choices = c("euclidean","maximum","manhattan","canberra","binary","minkowski")),
                    selectInput(inputId = 'cluster.method',label = 'Clustering method',
                                choices = c("ward.D","ward.D2","average","complete","single","mcquitty","median","centroid"),
                                selected = 'complete'),
                    numericInput(inputId = 'cluster.number',label = 'Cluster number',value = 10,min = 1,max = 50),
                    HTML('See the R functions for hierarchical clustering: <a href="https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust" target="_blank">hclust</a> and
                         <a href="https://www.rdocumentation.org/packages/proxy/versions/0.4-22/topics/dist" target="_blank">dist</a>'),
                    br(),
                    br(),
                    actionButton(inputId = 'plot.heatmap',label = 'Plot heatmap',
                                 icon = icon('pencil',lib = 'font-awesome'),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                    ),
                box(title=NULL,
                    width = 9,status = 'primary', solidHeader = T,
                    plotOutput('heatmap.plot.panel',height = 730),
                    actionButton(inputId = 'save.target.in.cluster',label = 'Save target list in clusters',
                                 icon = icon('download',lib = 'font-awesome'),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                    actionButton(inputId = 'save.heatmap',label = 'Save heatmap',
                                 icon = icon('download',lib = 'font-awesome'),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")

                )
              ),
              fluidRow(
                ##---------------multiple plots----------------
                # column(width = 12,
                box(title='Profiles',
                    width = 12,status = 'primary', solidHeader = T,
                    column(width = 6,
                           wellPanel(
                             HTML('<h4><strong>Visualise single gene plot</strong></h4>'),
                             textInput(inputId = 'gene.id',label = 'Input a gene',value = 'AT1G01060'),
                             br(),
                             radioButtons(inputId = 'profile.data.type',label = 'Select data type',
                                          choices = c('TPM','Read counts'),inline = T),
                             radioButtons(inputId = 'profile.filter.lowexpressed',label = 'Filter low expressed transcripts?',
                                          choices = c('Yes','No'),inline = T),
                             actionButton(inputId = 'make.profile.plot',label = 'Plot',
                                          icon = icon('pencil',lib = 'font-awesome'),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4; float"),
                             br()
                           )
                    ),
                    column(width = 6,
                           wellPanel(
                             HTML('<h4><strong>Make multiple gene plots</strong></h4>'),
                             fileInput("multiple.gene.input", "Choose gene list csv file",
                                       accept = c(
                                         "text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")
                             ),
                             radioButtons(inputId = 'multiple.plot.type',label = 'Select plot type',
                                          choices = c('Abundance','PS','Both'),inline = T),
                             radioButtons(inputId = 'multi.plot.format',label = 'Select format',
                                          choices = c('png','pdf','both'),inline = T),
                             actionButton(inputId = 'make.multiple.plot',label = 'Run',
                                          icon = icon('send outline icon',lib = 'font-awesome'),
                                          style="color: #fff; background-color: #428bca; border-color: #2e6da4; float"),
                             br()
                           )
                    )
                ),
                # column(width = 12,
                box(title=NULL,
                    width = 12,status = 'primary', solidHeader = T,
                    plotlyOutput('profile.plot.panel'),
                    br(),
                    hr(),
                    br(),
                    plotlyOutput('ps.plot.panel'),
                    actionButton(inputId = 'save.ps.plot',label = 'Save PS plot',
                                 icon = icon('download',lib = 'font-awesome'),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right"),
                    actionButton(inputId = 'save.abundance.plot',label = 'Save abundance plot',
                                 icon = icon('download',lib = 'font-awesome'),
                                 style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                )
              ),
              fluidRow(
                ##---------------GO plots----------------
                column(width = 3,
                       box(title = 'GO annotation',
                           width = 13,status = 'primary', solidHeader = T,
                           HTML('<h4><strong>Generate target list</strong></h4>'),
                           actionButton(inputId = 'generate.target.list',label = 'Generate',
                                        icon = icon('send outline icon',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                           br(),
                           br(),
                           HTML('Generate lists of DE genes, DAS genes, DE transcripts and DTU transcripts, which are the union
                                set across all contrast groups.'),
                           hr(),
                           fileInput("go.annotation.input", "Choose GO annotation csv file",
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                           ),
                           uiOutput('show.go.table.column'),
                           radioButtons(inputId = 'go.plot.label',label = 'Select gene list type',
                                        choices  = c('DE genes','DAS genes'),inline = T),
                           actionButton(inputId = 'make.go.plot',label = 'Plot',
                                        icon = icon('pencil',lib = 'font-awesome'),
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")

                           )
                ),
                column(width = 9,
                       box(title = NULL,
                           width = 13,status = 'primary', solidHeader = T,
                           h4('Loaded GO annotation table'),
                           DT::dataTableOutput('go.table')
                       )
                )
      ),
      fluidRow(
        box(title = NULL,
            width = 12,status = 'primary', solidHeader = T,
            h4('GO annotation plot'),
            uiOutput('show.go.plot.panel'),
            actionButton(inputId = 'save.go.plot',label = 'Save',
                         icon = icon('download',lib = 'font-awesome'),
                         style="color: #fff; background-color: #428bca; border-color: #2e6da4;float: right")
        )
      )
      ),
      ###########################################################################
      #======================= Generate report ============================
      tabItem("report",
              tags$head(tags$style("#TxtOut {white-space: normal;}")),
              box(title = 'Parameter summary',
                  width = 12,status = 'primary', solidHeader = T,
                  HTML('<h5><strong>Table: Paramter used for the analysis. </strong> If the parameters are incorrect, please
                       double click the the cell in the table to correct it.</h5>'),
                  DT::DTOutput('para.summary'),
                  br(),
                  br(),
                  actionButton(inputId = 'save.para.summary',label = 'Save',
                               icon = icon('download',lib = 'font-awesome'),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
                  ),
              actionButton(inputId = 'generate.report',label = 'Generate report',
                           icon = icon('send outline icon',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                           font-size:200%; padding: 30px 60px; margin-left:40%")
              )
      )
              )

  ui <- dashboardPage(
    mainheader,
    mainsidebar,
    mainbody
  )



  Logged = F
  my_username <- "User"
  my_password <- "test"

  server <- function(input, output, session) {
    # set data size
    options(shiny.maxRequestSize = data.size.max*1024^2)

    #Log in ##############################################################
    # ============================= password ============================ #
    values <- reactiveValues(authenticated = FALSE)

    # Return the UI for a modal dialog with data selection input. If 'failed'
    # is TRUE, then display a message that the previous value was invalid.
    dataModal <- function(failed = FALSE) {
      modalDialog(
        textInput("username", "Username:",value = 'User'),
        passwordInput("password", "Password:",value = 'test'),
        footer = tagList(
          # modalButton("Cancel"),
          div(style="display:inline-block;vertical-align:middle;height:40px;float: left",
              verbatimTextOutput('log.info')
          ),
          div(style="display:inline-block;vertical-align:middle;height:40px",
              actionButton("ok", "OK",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4")
          )
        ),size = 's'
      )
    }

    # Show modal when button is clicked.
    # This `observe` is suspended only whith right user credential

    obs1 <- observe({
      showModal(dataModal())
    })

    # When OK button is pressed, attempt to authenticate. If successful,
    # remove the modal.

    obs2 <- observe({
      req(input$ok)
      isolate({
        Username <- input$username
        Password <- input$password
      })
      Id.username <- which(my_username == Username)
      Id.password <- which(my_password == Password)
      if (length(Id.username) > 0 & length(Id.password) > 0) {
        if (Id.username == Id.password) {
          Logged <<- TRUE
          values$authenticated <- TRUE
          obs1$suspend()
          removeModal()
        } else {
          values$authenticated <- FALSE
        }
      }
    })

    observeEvent(input$ok,{
      output$log.info <- renderText({
        if (values$authenticated) "OK!!!!!"
        else "Wrong password!!!"
      })
    })
    # ============================= password ============================ #

    DDD.data <- reactiveValues(
      data.folder=NULL,
      figure.folder=NULL,
      result.folder=NULL,
      path=getwd(),
      mapping=NULL,
      samples=NULL,
      samples_new=NULL,
      trans_counts=NULL,
      genes_counts=NULL,
      trans_TPM=NULL,
      genes_TPM=NULL,
      target_high=NULL,
      trans_batch=NULL,
      genes_batch=NULL,
      trans_dge=NULL,
      genes_dge=NULL,
      contrast=NULL,
      deltaPS=NULL,
      genes_log2FC=NULL,
      trans_log2FC=NULL,
      DE_genes=NULL,
      DE_genes_stat=NULL,
      DAS_genes=NULL,
      DE_trans=NULL,
      DTU_trans=NULL,
      loaded.data=NULL,
      genes.ann=NULL,
      trans.ann=NULL
      # info=list()
    )

    output$messageMenu <- renderMenu({
      dropdownMenu(type = "tasks", badgeStatus = "success",
                   taskItem(value = 90, color = "green",
                            "Documentation"
                   ),
                   taskItem(value = 17, color = "aqua",
                            "Project X"
                   ),
                   taskItem(value = 75, color = "yellow",
                            "Server deployment"
                   ),
                   taskItem(value = 80, color = "red",
                            "Overall project"
                   )
      )
    })


    ##############################################################
    ######                    Data generation                 ####
    ##############################################################

    ##------------->>  Loaded data in workspace   <<--------------

    observe({
      Data <- c(
        'mapping',
        'samples',
        'samples_new',
        'trans_counts',
        'genes_counts',
        'trans_TPM',
        'genes_TPM',
        'target_high',
        'trans_batch',
        'genes_batch',
        'trans_dge',
        'genes_dge',
        'contrast',
        'deltaPS',
        'DE_genes',
        'DAS_genes',
        'DE_trans',
        'DTU_trans'
      )
      if(is.null(DDD.data$data.folder)){
        Existed1 <- F
      } else {
        Existed1 <- sapply(Data,function(i) file.exists(paste0(DDD.data$data.folder,'/',i,'.RData')))
      }
      DDD.data.list <- reactiveValuesToList(DDD.data)
      Existed2 <- sapply(1:length(Data),function(i) !is.null(DDD.data.list[[Data[i]]]))

      Description=c('Gene-transcript mapping',
                    'Sample information',
                    'Sample information (Sequencing replicates are merged if exist)',
                    'Transcript read counts (sequencing replicates are merged if exist)',
                    'Gene read counts (sequencing replicates are merged if exist)',
                    'Transcript level TPM',
                    'Gene level TPM',
                    'Remaining transcripts and genes after filtering the low expressed',
                    'Estimated transcript batch effect only if exist',
                    'Estimated gene batch effect only if exist',
                    'Transcript level normalisation factor',
                    'Gene level normalisation factor',
                    'Contrast groups of conditions for comparisons',
                    'delta Percent spliced value based on contrast groups',
                    'DE genes in different contrast groups',
                    'DAS genes in different contrast groups',
                    'DE transcripts in different contrast groups',
                    'DTU transcripts in different contrast groups'
      )

      Generate=c('User provided',
                 'User provided',
                 'Data pre-processing',
                 'Data pre-processing',
                 'Data pre-processing',
                 'Data generation',
                 'Data generation',
                 'Data pre-processing',
                 'Data pre-processing',
                 'Data pre-processing',
                 'Data pre-processing',
                 'Data pre-processing',
                 'User provided',
                 'DE DAS and DTU',
                 'DE DAS and DTU',
                 'DE DAS and DTU',
                 'DE DAS and DTU',
                 'DE DAS and DTU'
      )

      Use=c('All steps',
            'All steps',
            'All steps',
            'Data pre-processing and DE DAS and DTU',
            'Data pre-processing and DE DAS and DTU',
            'DE DAS and DTU and Advanced plot',
            'Advanced plot',
            'All steps',
            'Data pre-processing and DE DAS and DTU',
            'Data pre-processing and DE DAS and DTU',
            'Data pre-processing and DE DAS and DTU',
            'Data pre-processing and DE DAS and DTU',
            'DE DAS and DTU and Result summary',
            'DE DAS and DTU',
            'Result summary and Advanced plot',
            'Result summary and Advanced plot',
            'Result summary and Advanced plot',
            'Result summary and Advanced plot'
      )

      x <- data.frame(Data=Data,`In data folder`=Existed1,
                      `In workspace`=Existed2,
                      Description=Description,
                      `Generate from`=Generate,
                      `Use in`=Use)
      colnames(x) <- c('Data','In data folder','In workspace','Description','Generate from','Use in')
      DDD.data$loaded.data <- x
    })

    DDD.data.size <- reactive({
      x <- sum(sapply(reactiveValuesToList(DDD.data),object.size))
      utils:::format.object_size(x, "auto")
    })
    output$DDD.data.size <- renderText({
      paste0('Data in use: ',DDD.data.size())
    })

    output$loaded.data.info <- renderTable(
      DDD.data$loaded.data
    )

    observeEvent(input$loaded.data.all,{
      idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F & DDD.data$loaded.data$`In data folder`==T]
      withProgress(message = 'Loading exisiting data',
                   detail = 'This may take a while...', value = 0, {
                     for(i in idx){
                       incProgress(1/length(idx))
                       cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                       load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                       if(!exists(i))
                         next
                       DDD.data[[i]] <- get(i)
                     }
                   })
      message('Done!!!')
      showNotification("Done!!!")

      if(length(idx==0))
        return(NULL)
      ###disable the buttons
    })

    observeEvent(input$loaded.data.preprocessing,{
      idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                         DDD.data$loaded.data$`In data folder`==T &
                                         grepl(pattern = 'All steps|Data pre-processing',x = DDD.data$loaded.data$`Use in`)]
      withProgress(message = 'Loading exisiting data',
                   detail = 'This may take a while...', value = 0, {
                     for(i in idx){
                       incProgress(1/length(idx))
                       cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                       load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                       if(!exists(i))
                         next
                       DDD.data[[i]] <- get(i)
                     }
                   })
      message('Done!!!')
      showNotification("Done!!!")

      if(length(idx==0))
        return(NULL)
      ###disable the buttons
    })

    observeEvent(input$loaded.data.ddd,{
      idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                         DDD.data$loaded.data$`In data folder`==T &
                                         grepl(pattern = 'All steps|DE DAS and DTU',x = DDD.data$loaded.data$`Use in`)]
      withProgress(message = 'Loading exisiting data',
                   detail = 'This may take a while...', value = 0, {
                     for(i in idx){
                       incProgress(1/length(idx))
                       cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                       load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                       if(!exists(i))
                         next
                       DDD.data[[i]] <- get(i)
                     }
                   })
      message('Done!!!')
      showNotification("Done!!!")

      if(length(idx==0))
        return(NULL)
      ###disable the buttons
    })


    observeEvent(input$loaded.data.summary,{
      idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                         DDD.data$loaded.data$`In data folder`==T &
                                         grepl(pattern = 'All steps|Result summary',x = DDD.data$loaded.data$`Use in`)]
      withProgress(message = 'Loading exisiting data',
                   detail = 'This may take a while...', value = 0, {
                     for(i in idx){
                       incProgress(1/length(idx))
                       cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                       load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                       if(!exists(i))
                         next
                       DDD.data[[i]] <- get(i)
                     }
                   })
      message('Done!!!')
      showNotification("Done!!!")

      if(length(idx==0))
        return(NULL)
      ###disable the buttons
    })

    observeEvent(input$loaded.data.advanced,{
      idx <- DDD.data$loaded.data$Data[DDD.data$loaded.data$`In workspace`==F &
                                         DDD.data$loaded.data$`In data folder`==T &
                                         grepl(pattern = 'All steps|Advanced plot',x = DDD.data$loaded.data$`Use in`)]
      withProgress(message = 'Loading exisiting data',
                   detail = 'This may take a while...', value = 0, {
                     for(i in idx){
                       incProgress(1/length(idx))
                       cat(paste0('\nLoading: ',DDD.data$data.folder,'/',i,'.RData'))
                       load(paste0(DDD.data$data.folder,'/',i,'.RData'))
                       if(!exists(i))
                         next
                       DDD.data[[i]] <- get(i)
                     }
                   })
      message('Done!!!')
      showNotification("Done!!!")

      if(length(idx==0))
        return(NULL)
      ###disable the buttons
    })


    ##---------------working directory-----------------
    shinyDirChoose(input, 'wdir', roots = c(getVolumes()()))
    wdir <- reactive({
      input$wdir
    })

    observeEvent(ignoreNULL = T,
                 eventExpr = {
                   input$wdir
                 },
                 handlerExpr = {
                   if(is.integer(wdir())){
                     path.folder <- NULL
                   } else {
                     root <- getVolumes()()[wdir()$root]
                     root <- gsub('[\\]','',normalizePath(root))
                     current_wk <- unlist(wdir()$path[-1])
                     path.folder <- file.path(root, paste(current_wk, collapse = .Platform$file.sep))
                     setwd(path.folder)
                   }
                   DDD.data$path <- path.folder

                 })

    output$wdir.info <- renderText({
      DDD.data$path
    })

    #####data save folder
    data.folder <- observe({
      data.folder <- paste0(DDD.data$path,'/data')
      if(!file.exists(data.folder))
        dir.create(data.folder,recursive = T)
      DDD.data$data.folder <- data.folder
    })

    #####figure save folder
    figure.folder <- observe({
      figure.folder <- paste0(DDD.data$path,'/figure')
      if(!file.exists(figure.folder))
        dir.create(figure.folder,recursive = T)
      DDD.data$figure.folder <- figure.folder
    })

    #####results save folder
    result.folder <- observe({
      result.folder <- paste0(DDD.data$path,'/result')
      if(!file.exists(result.folder))
        dir.create(result.folder,recursive = T)
      DDD.data$result.folder <- result.folder
    })

    ##mapping table####
    observe({
      inFile <- input$mapping.input
      if (is.null(inFile))
        return(NULL)
      # withProgress(message = 'Loading mapping information', value = 0, {
      mapping <- read.csv(inFile$datapath, header = T)
      rownames(mapping) <- mapping$TXNAME
      save(mapping,file=paste0(DDD.data$data.folder,'/mapping.RData'))
      DDD.data$mapping <- mapping
      # })
    })

    output$mapping.output <- DT::renderDataTable({
      DDD.data$mapping[1:25,]
    },options = list(scrollX = TRUE,
                     columnDefs = list(list(className = 'dt-center', targets="_all"))
    ),rownames= FALSE)

    ##sample table####
    samples <- observe({
      inFile <- input$sample.input
      if (is.null(inFile))
        return(NULL)
      # withProgress(message = 'Loading sample information', value = 0, {
      samples <- read.csv(inFile$datapath, header = T)
      # })
      if('srep' %in% colnames(samples)){
        samples$sample.name <- paste0(samples$condition,'_',samples$brep,'_',samples$srep)
      } else {
        samples$sample.name <- paste0(samples$condition,'_',samples$brep)
      }
      save(samples,file=paste0(DDD.data$data.folder,'/samples.RData'))

      DDD.data$samples <- samples

    })

    output$sample.output <- DT::renderDataTable({
      if(!is.null(DDD.data$samples_new))
        DDD.data$samples_new else DDD.data$samples
    },options = list(scrollX = TRUE,
                     columnDefs = list(list(className = 'dt-left', targets="_all"))))

    ##generate genes expression####
    observe({
      if(input$generate.new.genes.data=='Load txi_genes.RData'){
        updateActionButton(session,inputId = 'run.txi.genes',label = 'Load',icon = icon('upload'))
      } else {
        updateActionButton(session,inputId = 'run.txi.genes',label = 'Run')
      }
    })

    observeEvent(input$run.txi.genes,{
      if (is.null(DDD.data$samples))
        return(NULL)
      if(input$generate.new.genes.data=='Load txi_genes.RData'){
        cat('\nLoading gene expression data: txi_genes.RData')
        load(paste0(DDD.data$data.folder,'/txi_genes.RData'))
        txi_genes <- txi_genes
        message(paste0('txi_genes.RData is loaded from folder: ',DDD.data$data.folder))
      } else {
        cat('\nRunning tximport to generate gene expression using method: ',
            input$tximport.method)

        data.files <-  as.character(DDD.data$samples$path)
        if(!all(file.exists(data.files))){
          warning('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
          showNotification('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
          return(NULL)
        }

        withProgress(message = 'Generate gene expression...', value = 0, {
          incProgress(0.1)
          txi_genes <- tximport(data.files,dropInfReps = T,
                                type = input$tximport.quant.method, tx2gene = DDD.data$mapping,
                                countsFromAbundance = input$tximport.method)
          incProgress(0.7)
          colnames(txi_genes$counts)<-DDD.data$samples$sample.name
          write.csv(txi_genes$counts,file=paste0(DDD.data$result.folder,'/counts_genes.csv'))
          incProgress(0.8)
          colnames(txi_genes$abundance)<-DDD.data$samples$sample.name
          write.csv(txi_genes$abundance,file=paste0(DDD.data$result.folder,'/TPM_genes.csv'))
          incProgress(0.9)
          colnames(txi_genes$length)<-DDD.data$samples$sample.name
          save(txi_genes,file=paste0(DDD.data$data.folder,'/txi_genes.RData'))
          genes_TPM <- txi_genes$abundance
          save(genes_TPM,file=paste0(DDD.data$data.folder,'/genes_TPM.RData'))
          message(paste0('txi_genes.RData and genes_TPM.RData are saved in folder: ',DDD.data$data.folder))
          incProgress(1)
        })
      }
      DDD.data$genes_counts=txi_genes$counts
      DDD.data$genes_TPM <- txi_genes$abundance
      message('Done!!!')
      showNotification("Done!!!")
    })

    rungeneinfo <- reactive({
      if(input$run.txi.genes>0){
        if(input$generate.new.genes.data=='Load txi_genes.RData'){
          paste0('txi_genes.RData is loaded from the data folder.')
        } else {
          paste0('txi_genes.RData is saved in the data folder.')
        }
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/txi_genes.RData'))){
          paste0('txi_genes.RData is already in the data folder.\nGenerate new?')
        } else {
          paste0('txi_genes.RData is not in the data folder.\nPlease generate new.')
        }
      }
    })

    output$rungeneinfo <- renderText({
      rungeneinfo()
    })

    ##generate trans expression####
    ###---update actionbutton if load data
    observe({
      if(input$generate.new.trans.data=='Load txi_trans.RData'){
        updateActionButton(session,inputId = 'run.txi.trans',label = 'Load',icon = icon('upload'))
      } else {
        updateActionButton(session,inputId = 'run.txi.trans',label = 'Run')
      }
    })

    observeEvent(input$run.txi.trans,ignoreNULL = T,{
      if (is.null(DDD.data$samples))
        return(NULL)

      if(input$generate.new.trans.data=='Load txi_trans.RData'){
        cat('\nLoading transcript expression data: txi_trans.RData')
        load(paste0(DDD.data$data.folder,'/txi_trans.RData'))
        txi_trans <- txi_trans
        message(paste0('txi_trans.RData is loaded from folder: ',DDD.data$data.folder))
      } else {
        cat(paste0('\nRunning tximport to generate transcript expression using method: ',
                   input$tximport.method))

        data.files <-  as.character(DDD.data$samples$path)
        if(!all(file.exists(data.files))){
          warning('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
          showNotification('The paths of quant.sf in the samples.csv file are incorrect. Please double check.')
          return(NULL)
        }

        # cat(paste0(data.files[1]))
        withProgress(message = 'Generate transcript expression...', value = 0, {
          incProgress(0.1)
          txi_trans<- tximport(data.files, type = input$tximport.quant.method, tx2gene = NULL,
                               countsFromAbundance = input$tximport.method,
                               txOut = T,dropInfReps = T)
          incProgress(0.7)
          colnames(txi_trans$counts)<-DDD.data$samples$sample.name
          write.csv(txi_trans$counts,file=paste0(DDD.data$result.folder,'/counts_trans.csv'))
          incProgress(0.8)
          colnames(txi_trans$abundance)<-DDD.data$samples$sample.name
          write.csv(txi_trans$abundance,file=paste0(DDD.data$result.folder,'/TPM_trans.csv'))
          incProgress(0.9)
          colnames(txi_trans$length)<-DDD.data$samples$sample.name
          save(txi_trans,file=paste0(DDD.data$data.folder,'/txi_trans.RData'))
          trans_TPM <- txi_trans$abundance
          save(trans_TPM,file=paste0(DDD.data$data.folder,'/trans_TPM.RData'))
          message(paste0('txi_trans.RData and trans_TPM.RData are saved in folder: ',DDD.data$data.folder))
          incProgress(1)
        })
      }
      DDD.data$trans_counts <- txi_trans$counts
      DDD.data$trans_TPM <- txi_trans$abundance
      message('Done!!!')
      showNotification("Done!!!")
    })


    runtransinfo <- reactive({
      if(input$run.txi.trans>0){
        if(input$generate.new.trans.data=='Load txi_trans.RData'){
          paste0('txi_trans.RData is loaded from the data folder.')
        } else {
          paste0('txi_trans.RData is saved in the data folder.')
        }
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/txi_trans.RData'))){
          paste0('txi_trans.RData is already in the data folder.\nGenerate new?')
        } else {
          paste0('txi_trans.RData is not in the data folder.\nPlease generate new.')
        }
      }
    })

    output$runtransinfo <- renderText({
      runtransinfo()
    })


    ##############################################################
    ######                data pre-processing
    ##############################################################

    ##------------->>  Step 1: merge sequencing reps


    seqrepinfo <- reactive({
      if(is.null(DDD.data$samples)){
        info <- HTML('Load samples.csv to check seq-reps.')
      } else {
        if(all(DDD.data$samples$srep=='srep1')){
          info <- HTML('The data has no sequencing replicates. Please skip.')
        } else if(!is.null(DDD.data$samples_new)) {
          info <- HTML('The sequencing replicates have been merged.')
        } else {
          info <- HTML(paste0('The data has ', length(unique(DDD.data$samples$srep)),' sequencing replicates. Please merge.'))
        }
      }
      return(info)
    })
    output$seqrep.info <- renderText(seqrepinfo())

    observe({
      if(input$merge.or.skip=='Merge'){
        updateActionButton(session,inputId = 'run.seq.merge',label = 'Merge',icon = icon("plus-square"))
      } else {
        updateActionButton(session,inputId = 'run.seq.merge',label = 'Skip',icon = icon("angle-right"))
      }
    })



    ###---update the sample information after merge---
    observeEvent(input$run.seq.merge,{
      if(is.null(DDD.data$samples))
        return(NULL)
      if(input$merge.or.skip=='Merge'){
        samples_new <- DDD.data$samples[DDD.data$samples$srep==DDD.data$samples$srep[1],]
        samples_new$srep <- 'merged'
      } else {
        samples_new <- DDD.data$samples
      }
      save(samples_new,file=paste0(DDD.data$data.folder,'/samples_new.RData'))
      message(paste0('samples_new.RData is saved in folder: ',DDD.data$data.folder))
      DDD.data$samples_new <- samples_new

    })

    ###---updateRadioButtons
    observe({
      if(is.null(DDD.data$samples))
        return(NULL)
      if('srep' %in% names(DDD.data$samples)){
        x <- 'Merge'
      } else {
        x <- 'Skip'
      }
      updateRadioButtons(session,label = 'Merge options', "merge.or.skip",selected = x)
    })
    #
    observeEvent(input$run.seq.merge,{
      if(is.null(DDD.data$samples) | is.null(DDD.data$trans_counts))
        return(NULL)

      if(input$merge.or.skip=='Merge'){
        cat('\nMerge sequencing replicates')
        idx <- paste0(DDD.data$samples$condition,'_',DDD.data$samples$brep)
        if(length(idx)==ncol(DDD.data$trans_counts)){
          trans_counts <- sumarrays(DDD.data$trans_counts,group = idx)
        } else {
          trans_counts <- DDD.data$trans_counts
        }
      } else {
        trans_counts <- DDD.data$trans_counts
      }
      save(trans_counts,file=paste0(DDD.data$data.folder,'/trans_counts.RData'))
      message(paste0('trans_counts.RData is saved in folder: ',DDD.data$data.folder))

      DDD.data$trans_counts <- trans_counts
    })

    observeEvent(input$run.seq.merge,{
      if(is.null(DDD.data$samples) | is.null(DDD.data$genes_counts))
        return(NULL)

      if(input$merge.or.skip=='Merge'){
        idx <- paste0(DDD.data$samples$condition,'_',DDD.data$samples$brep)
        if(length(idx)==ncol(DDD.data$genes_counts)){
          genes_counts <- sumarrays(DDD.data$genes_counts,group = idx)
        } else {
          genes_counts <- DDD.data$genes_counts
        }
      } else {
        genes_counts <- DDD.data$genes_counts
      }
      save(genes_counts,file=paste0(DDD.data$data.folder,'/genes_counts.RData'))
      message(paste0('genes_counts.RData is saved in folder: ',DDD.data$data.folder))

      showNotification("Done!!!")
      DDD.data$genes_counts <- genes_counts
    })



    observeEvent(input$run.seq.merge,{
      if(is.null(DDD.data$samples_new)){
        max.n <- 10
      } else {
        max.n <- nrow(DDD.data$samples_new)
      }
      updateSliderInput(session, "sample.n.cut", max=max.n)
    })
    #
    ##----------Step 2: filter low expressed genes------------

    observeEvent(input$run.filter,{
      if(is.null(DDD.data$trans_counts) | is.null(DDD.data$mapping))
        return(NULL)
      cat('\nFilter low expression')
      target_high <- counts.filter(counts = DDD.data$trans_counts,
                                   mapping = DDD.data$mapping,
                                   cpm.cut = input$cpm.cut,
                                   sample.n = input$sample.n.cut,
                                   Log=F)
      save(target_high,file=paste0(DDD.data$data.folder,'/target_high.RData'))
      message(paste0('target_high.RData is saved in folder: ',DDD.data$data.folder))


      message('Done!!!')
      showNotification("Done!!!")
      DDD.data$target_high <- target_high
    })

    low.filter.info <- reactive({
      if(is.null(DDD.data$data.folder))
        return(NULL)
      if(input$run.filter>0){
        paste0('target_high.RData is saved in the data folder.')
      } else {
        if(!is.null(DDD.data$target_high)){
          paste0('target_high.RData is already uploaded in the workspace.')
        } else {
          if(file.exists(paste0(DDD.data$data.folder,'/target_high.RData'))){
            paste0('target_high.RData is already in "data" folder.\nCan upload directly in "Data generation".')
          } else {
            paste0('target_high.RData is not in "data" folder.\nPlease generate new.')
          }
        }
      }
    })


    output$low.filter.info <- renderText({
      if(is.null(low.filter.info()))
        return(NULL)
      low.filter.info()
    })



    ##--------------mean-variance plots---------------
    ###--mean variance trans trend plot
    mv.trans <- eventReactive(input$mv.plot.button,{
      condition <- paste0(DDD.data$samples_new$condition)
      ###---transcript levels
      counts.raw = DDD.data$trans_counts[rowSums(DDD.data$trans_counts>0)>0,]
      counts.filtered = DDD.data$trans_counts[DDD.data$target_high$trans_high,]
      # graphics.off()
      check.mean.variance(counts.raw = counts.raw,
                          counts.filtered = counts.filtered,
                          condition = condition)
    })

    mv.trans.plot <- eventReactive(input$mv.plot.button,{
      mv.trans.plot <- function(){
        fit.raw <- mv.trans()$fit.raw
        fit.filtered <- mv.trans()$fit.filtered
        par(mfrow=c(1,2))
        plot.mean.variance(x = fit.raw$sx,y = fit.raw$sy,
                           l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
        title('\n\nRaw counts')
        plot.mean.variance(x = fit.filtered$sx,y = fit.filtered$sy,
                           l = fit.filtered$l,lwd=2,col='black')
        title('\n\nFiltered counts')
        lines(fit.raw$l, col = "gold",lty=4,lwd=2)
        legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
               legend = c('low-exp removed','low-exp kept'))
      }
      return(mv.trans.plot)
    })
    #
    output$mv.trans.plot <- renderPlot({
      mv.trans.plot()()
    })

    ###--mean variance genes trend plot
    mv.genes <- eventReactive(input$mv.plot.button,{
      condition <- paste0(DDD.data$samples_new$condition)
      ###---genes levels
      counts.raw = DDD.data$genes_counts[rowSums(DDD.data$genes_counts>0)>0,]
      counts.filtered = DDD.data$genes_counts[DDD.data$target_high$genes_high,]
      # graphics.off()
      check.mean.variance(counts.raw = counts.raw,
                          counts.filtered = counts.filtered,
                          condition = condition)
    })

    mv.genes.plot <- eventReactive(input$mv.plot.button,{
      mv.genes.plot <- function(){
        fit.raw <- mv.genes()$fit.raw
        fit.filtered <- mv.genes()$fit.filtered
        par(mfrow=c(1,2))
        plot.mean.variance(x = fit.raw$sx,y = fit.raw$sy,
                           l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
        title('\n\nRaw counts')
        plot.mean.variance(x = fit.filtered$sx,y = fit.filtered$sy,
                           l = fit.filtered$l,lwd=2,col='black')
        title('\n\nFiltered counts')
        lines(fit.raw$l, col = "gold",lty=4,lwd=2)
        legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
               legend = c('low-exp removed','low-exp kept'))
      }
      return(mv.genes.plot)
    })
    #
    output$mv.genes.plot <- renderPlot({
      mv.genes.plot()()
    })

    observeEvent(input$save.mv.plot,{
      ##transcript level
      png(filename = paste0(DDD.data$figure.folder,'/Transcript mean-variance trend.png'),
          width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
      mv.trans.plot()()
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Transcript mean-variance trend.pdf'),
          width = 25/2.54,height = 12/2.54)
      mv.trans.plot()()
      dev.off()

      ###gene level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Gene mean-variance trend.png'),
          width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
      mv.genes.plot()()
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Gene mean-variance trend.pdf'),
          width = 25/2.54,height = 12/2.54)
      mv.genes.plot()()
      dev.off()
      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    })

    output$save.mv.genes.plot.text <- renderText({
      if(input$save.mv.genes.plot==0)
        return(NULL)
      else HTML(paste0('Figure is saved in: ',DDD.data$figure.folder))
    })


    ##--------------Step 3: PCA plot---------------
    ###update the pca selections in the selectinput

    observe({
      if(is.null(input$sample.input)){
        x <- 5
      } else {
        if(input$pca.plot.type=='Bio-reps'){
          x <- length(DDD.data$samples_new$condition)
        } else {
          x <- length(unique(DDD.data$samples_new$condition))
        }
      }

      x <- paste0('PC',1:x)
      updateSelectInput(session, "dim1",
                        label = 'x-axis',
                        choices = x,
                        selected = 'PC1'
      )
      updateSelectInput(session, "dim2",
                        label = 'y-axis',
                        choices = x,
                        selected = 'PC2'
      )
    })

    ###---trans level
    pca.trans.g <- eventReactive(input$plot.pca.button,{
      data2pca <- DDD.data$trans_counts[DDD.data$target_high$trans_high,]
      dge <- DGEList(counts=data2pca)
      dge <- calcNormFactors(dge)
      data2pca <- t(counts2CPM(obj = dge,Log = T))
      dim1 <- input$dim1
      dim2 <- input$dim2
      ellipse.type <- input$pca.add.circle

      ##--Bio-reps
      if(input$pca.plot.type=='Bio-reps'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        if(input$pca.color.option=='Bio-reps')
          groups <- DDD.data$samples_new$brep
        if(input$pca.color.option=='Conditions')
          groups <- DDD.data$samples_new$condition
        g <- plot.PCA.ind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                          groups = groups,plot.title = 'Transcript PCA: bio-reps',
                          ellipse.type = ellipse.type,
                          add.label = T,adj.label = F)
      }

      ##--average expression plot
      if(input$pca.plot.type=='Average expression'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        groups <- DDD.data$samples_new$brep
        data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
        groups <- unique(DDD.data$samples_new$condition)
        g <- plot.PCA.ind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                          groups = groups,plot.title = 'Transcript PCA: average expression',
                          ellipse.type = 'none',add.label = T,adj.label = F)
      }
      g
    })

    output$trans.pca.plot <- renderPlotly({
      ggplotly(pca.trans.g())
    })

    ###---genes level
    pca.genes.g <- eventReactive(input$plot.pca.button,{
      data2pca <- DDD.data$genes_counts[DDD.data$target_high$genes_high,]
      dge <- DGEList(counts=data2pca)
      dge <- calcNormFactors(dge)
      data2pca <- t(counts2CPM(obj = dge,Log = T))
      dim1 <- input$dim1
      dim2 <- input$dim2
      ellipse.type <- input$pca.add.circle

      ##--Bio-reps
      if(input$pca.plot.type=='Bio-reps'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        if(input$pca.color.option=='Bio-reps')
          groups <- DDD.data$samples_new$brep
        if(input$pca.color.option=='Conditions')
          groups <- DDD.data$samples_new$condition
        g <- plot.PCA.ind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                          groups = groups,plot.title = 'Gene PCA: bio-reps',
                          ellipse.type = ellipse.type,
                          add.label = T,adj.label = F)
      }

      ##--average expression plot
      if(input$pca.plot.type=='Average expression'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        groups <- DDD.data$samples_new$brep
        data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
        groups <- unique(DDD.data$samples_new$condition)
        g <- plot.PCA.ind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                          groups = groups,plot.title = 'Gene PCA: average expression',
                          ellipse.type = 'none',add.label = T,adj.label = F)
      }
      g
    })

    output$genes.pca.plot <- renderPlotly({
      ggplotly(pca.genes.g())
    })

    observeEvent(input$save.pca.plot,{
      ##transcript level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA_',input$pca.plot.type,'.png'),
          width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
      print(pca.trans.g())
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA_',input$pca.plot.type,'.pdf'),
          width = 15/2.54,height = 13/2.54)
      print(pca.trans.g())
      dev.off()
      ###gene level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Gene PCA_',input$pca.plot.type,'.png'),
          width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
      print(pca.genes.g())
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA_',input$pca.plot.type,'.pdf'),
          width = 15/2.54,height = 13/2.54)
      print(pca.genes.g())
      dev.off()
      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    })




    ##--------------remove batch effects-trans level---------------
    ###---update remove.batch action button

    ###---transcript level
    observeEvent(input$remove.batch,{
      cat('\nEstimate transcript level batch effects...')
      trans_batch <- remove.batch.shiny(read.counts = DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                                        condition = DDD.data$samples_new$condition,
                                        design = NULL,
                                        contrast=NULL,
                                        group = DDD.data$samples_new$condition,
                                        method = input$ruvseq.method)
      save(trans_batch,file=paste0(DDD.data$data.folder,'/trans_batch.RData'))
      message(paste0('trans_batch.RData is saved in folder: ',DDD.data$data.folder))
      DDD.data$trans_batch <- trans_batch
    })

    ###---gene level
    observeEvent(input$remove.batch,{
      cat('\nEstimate gene level batch effects')
      genes_batch <- remove.batch.shiny(read.counts = DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                                        condition = DDD.data$samples_new$condition,
                                        design = NULL,
                                        contrast=NULL,
                                        group = DDD.data$samples_new$condition,
                                        method = input$ruvseq.method)
      save(genes_batch,file=paste0(DDD.data$data.folder,'/genes_batch.RData'))
      message(paste0('genes_batch.RData is saved in folder: ',DDD.data$data.folder))
      DDD.data$genes_batch <- genes_batch
    })

    br.saved.info <- reactive({
      if(is.null(DDD.data$data.folder))
        return(NULL)
      if(input$remove.batch>0){
        paste0('trans(genes)_batch.RData are saved in the data folder.')
      } else {
        if(!is.null(DDD.data$genes_batch) & !is.null(DDD.data$trans_batch)){
          paste0('trans(genes)_batch.RData are already uploaded in the workspace.')
        } else {
          if(file.exists(paste0(DDD.data$data.folder,'/genes_batch.RData')) &
             file.exists(paste0(DDD.data$data.folder,'/trans_batch.RData'))){
            paste0('trans(genes)_batch.RData are already in "data" folder.\nCan upload directly in "Data generation".')
          } else {
            paste0('trans(genes)_batch.RData are not in "data" folder.\nPlease generate new.')
          }
        }
      }
    })

    output$br.saved.info <- renderText({
      if(is.null(br.saved.info()))
        return(NULL)
      br.saved.info()
    })
    #
    ##--------------Step 3-continued: PCA after remove batch---------------
    # ###---trans level
    pca.trans.br.g <- eventReactive(input$update.pca.plot,{
      message('Transcript level PCA--batch removed')
      data2pca <- DDD.data$trans_batch$normalizedCounts
      dge <- DGEList(counts=data2pca)
      dge <- suppressWarnings(calcNormFactors(dge))
      data2pca <- t(counts2CPM(obj = dge,Log = T))
      dim1 <- input$dim1
      dim2 <- input$dim2
      ellipse.type <- input$pca.add.circle

      ##--Bio-reps
      if(input$pca.plot.type=='Bio-reps'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        if(input$pca.color.option=='Bio-reps')
          groups <- DDD.data$samples_new$brep
        if(input$pca.color.option=='Conditions')
          groups <- DDD.data$samples_new$condition
        g <- plot.PCA.ind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                          groups = groups,plot.title = 'Transcript PCA: bio-reps',
                          ellipse.type = ellipse.type,
                          add.label = T,adj.label = F)
      }

      ##--average expression plot
      if(input$pca.plot.type=='Average expression'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        groups <- DDD.data$samples_new$brep
        data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
        groups <- unique(DDD.data$samples_new$condition)
        g <- plot.PCA.ind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                          groups = groups,plot.title = 'Transcript PCA: average expression',
                          ellipse.type = 'none',add.label = T,adj.label = F)
      }
      g
    })
    #

    output$trans.pca.br.plot <- renderPlotly({
      ggplotly(pca.trans.br.g())
    })


    # ###---genes level
    pca.genes.br.g <- eventReactive(input$update.pca.plot,{
      message('Gene level PCA--batch removed')
      data2pca <- DDD.data$genes_batch$normalizedCounts
      dge <- DGEList(counts=data2pca)
      dge <- suppressWarnings(calcNormFactors(dge))
      data2pca <- t(counts2CPM(obj = dge,Log = T))
      dim1 <- input$dim1
      dim2 <- input$dim2
      ellipse.type <- input$pca.add.circle

      ##--Bio-reps
      if(input$pca.plot.type=='Bio-reps'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        if(input$pca.color.option=='Bio-reps')
          groups <- DDD.data$samples_new$brep
        if(input$pca.color.option=='Conditions')
          groups <- DDD.data$samples_new$condition
        g <- plot.PCA.ind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                          groups = groups,plot.title = 'Gene PCA: bio-reps',
                          ellipse.type = ellipse.type,
                          add.label = T,adj.label = F)
      }

      ##--average expression plot
      if(input$pca.plot.type=='Average expression'){
        rownames(data2pca) <- gsub('_','.',rownames(data2pca))
        groups <- DDD.data$samples_new$brep
        data2pca.ave <- rowmean(data2pca,DDD.data$samples_new$condition,reorder = F)
        groups <- unique(DDD.data$samples_new$condition)
        g <- plot.PCA.ind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                          groups = groups,plot.title = 'genescript PCA: average expression',
                          ellipse.type = 'none',add.label = T,adj.label = F)
      }
      g
    })
    #

    output$genes.pca.br.plot <- renderPlotly({
      ggplotly(pca.genes.br.g())
    })

    observeEvent(input$save.pca.br.plot,{
      ###transcript level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed_',input$pca.plot.type,'.png'),
          width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
      print(pca.trans.br.g())
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Transcript PCA batch effect removed_',input$pca.plot.type,'.pdf'),
          width = 15/2.54,height = 13/2.54)
      print(pca.trans.br.g())
      dev.off()

      ###gene level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed_',input$pca.plot.type,'.png'),
          width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
      print(pca.genes.br.g())
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Gene PCA batch effect removed_',input$pca.plot.type,'.pdf'),
          width = 15/2.54,height = 13/2.54)
      print(pca.genes.br.g())
      dev.off()
      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    })

    ##--------------Step 4-normalisation--------------
    observeEvent(input$run.norm,{
      cat('\nNormalise transcript read counts')
      dge <- DGEList(counts=DDD.data$trans_counts[DDD.data$target_high$trans_high,],
                     group = DDD.data$samples_new$condition,
                     genes = DDD.data$mapping[DDD.data$target_high$trans_high,])
      trans_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm.method))
      save(trans_dge,file=paste0(DDD.data$data.folder,'/trans_dge.RData'))
      message(paste0('trans_dge.RData is saved in folder: ',DDD.data$data.folder))
      message('Done!!!')
      DDD.data$trans_dge <- trans_dge
    })

    observeEvent(input$run.norm,{
      cat('\nNormalise gene read counts')
      dge <- DGEList(counts=DDD.data$genes_counts[DDD.data$target_high$genes_high,],
                     group = DDD.data$samples_new$condition)
      genes_dge <- suppressWarnings(calcNormFactors(dge,method = input$norm.method))
      save(genes_dge,file=paste0(DDD.data$data.folder,'/genes_dge.RData'))
      message(paste0('genes_dge.RData is saved in folder: ',DDD.data$data.folder))
      message('Done!!!')
      DDD.data$genes_dge <- genes_dge
    })


    norm.info <- reactive({
      if(is.null(DDD.data$data.folder))
        return(NULL)
      if(input$run.norm>0){
        paste0('trans(genes)_dge.RData are saved in the data folder.')
      } else {
        if(!is.null(DDD.data$genes_dge) & !is.null(DDD.data$trans_dge)){
          paste0('trans(genes)_dge.RData are already uploaded in the workspace.')
        } else {
          if(file.exists(paste0(DDD.data$data.folder,'/genes_dge.RData')) &
             file.exists(paste0(DDD.data$data.folder,'/trans_dge.RData'))){
            paste0('trans(genes)_dge.RData are already in "data" folder.\nCan upload directly in "Data generation".')
          } else {
            paste0('trans(genes)_dge.RData are not in "data" folder.\nPlease generate new.')
          }
        }
      }
    })


    output$norm.info <- renderText({
      if(is.null(norm.info()))
        return(NULL)
      norm.info()
    })


    ##--------------make the distribution plot-------------
    ###---transcript level
    trans.dist.g <- eventReactive(input$norm.plot.button,{
      sample.name <- paste0(DDD.data$samples_new$condition,'.',DDD.data$samples_new$brep)
      condition <- DDD.data$samples_new$condition
      data.before <- DDD.data$trans_counts[DDD.data$target_high$trans_high,]
      data.after <- counts2CPM(obj = DDD.data$trans_dge,Log = T)
      g <- boxplot.normalised(data.before = data.before,
                              data.after = data.after,
                              condition = condition,
                              sample.name = sample.name)
      g
    })

    output$trans.dist.plot <- renderPlotly({
      q <- subplot(style(trans.dist.g()$g1+theme(legend.position = 'none',
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()), showlegend = FALSE),
                   trans.dist.g()$g2+labs(title='Data distribution (before and after normalisation)'),
                   nrows = 2,
                   titleX=T,
                   titleY = T,margin=0.04)
      q
    })


    ###---gene level
    genes.dist.g <- eventReactive(input$norm.plot.button,{
      sample.name <- paste0(DDD.data$samples_new$condition,'.',DDD.data$samples_new$brep)
      condition <- DDD.data$samples_new$condition
      data.before <- DDD.data$genes_counts[DDD.data$target_high$genes_high,]
      data.after <- counts2CPM(obj = DDD.data$genes_dge,Log = T)
      g <- boxplot.normalised(data.before = data.before,
                              data.after = data.after,
                              condition = condition,
                              sample.name = sample.name)
      g
    })

    output$genes.dist.plot <- renderPlotly({
      q <- subplot(style(genes.dist.g()$g1+theme(legend.position = 'none',
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()), showlegend = FALSE),
                   genes.dist.g()$g2+labs(title='Data distribution (before and after normalisation)'),
                   nrows = 2,
                   titleX=T,
                   titleY = T,margin=0.04)
      q
    })

    observeEvent(input$save.dist.plot,{
      ##transcript level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Transcript data distribution.png'),
          width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
      gridExtra::grid.arrange(trans.dist.g()$g1,trans.dist.g()$g2,ncol=1)
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Transcript data distribution.pdf'),
          width = 20/2.54,height = 20/2.54)
      gridExtra::grid.arrange(trans.dist.g()$g1,trans.dist.g()$g2,ncol=1)
      dev.off()

      ##gene level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Gene data distribution.png'),
          width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
      gridExtra::grid.arrange(genes.dist.g()$g1,genes.dist.g()$g2,ncol=1)
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Gene data distribution.pdf'),
          width = 20/2.54,height = 20/2.54)
      gridExtra::grid.arrange(genes.dist.g()$g1,genes.dist.g()$g2,ncol=1)
      dev.off()
      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    })

    ##---------------data information-----------------

    data.info <- observe({
      if(is.null(DDD.data$samples) | is.null(DDD.data$trans_counts) | is.null(DDD.data$genes_counts))
        return(NULL)
      data.info <- data.frame(
        raw.number=c(length(DDD.data$mapping$TXNAME),length(unique(DDD.data$mapping$GENEID))),
        sample=rep(nrow(DDD.data$samples),2),
        condition=rep(length(unique(DDD.data$samples$condition)),2),
        brep=rep(length(unique(DDD.data$samples$brep)),2),
        srep=rep(length(unique(DDD.data$samples$srep)),2)
      )

      if(!is.null(DDD.data$trans_counts))
        data.info$merged.sample <- c(ncol(DDD.data$trans_counts),ncol(DDD.data$genes_counts))

      if(!is.null(DDD.data$target_high)){
        data.info$cpm.cut <- c(input$cpm.cut,'--')
        data.info$sample.n.cut <- c(input$sample.n.cut,'--')
        data.info$after.filter <- c(length(DDD.data$target_high$trans_high),length(DDD.data$target_high$genes_high))
      }
      rownames(data.info) <- c('Transcript','Gene')
      DDD.data$data.info <- data.info
    })

    output$data.info <- DT::renderDataTable({
      DDD.data$data.info
    },
    options = list(
      scrollX=T,
      searching = FALSE,
      paging=F,
      info = FALSE,
      columnDefs = list(list(className = 'dt-center', targets="_all")))
    )

    observeEvent(input$save.data.info,{
      write.csv(DDD.data$data.info,file=paste0(DDD.data$result.folder,'/data.info.csv'),row.names = T)
      message(paste0('Table is saved in folder: ',DDD.data$result.folder))
      showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
    })

    ###########################################################################################
    ##DE DAS and DTU----
    ###########################################################################################

    ##Load contrast----------------
    observe({
      inFile <- input$contrast.input
      if (is.null(inFile))
        return(NULL)
      contrast <- read.csv(inFile$datapath, header = T)
      contrast <- paste0(contrast[,1],'-',contrast[,2])
      save(contrast,file=paste0(DDD.data$data.folder,'/contrast.RData'))
      DDD.data$contrast <- contrast

    })

    output$contrast.table<- DT::renderDataTable({
      if(is.null(DDD.data$contrast))
        return(NULL)
      x <- do.call(rbind,strsplit(DDD.data$contrast,'-'))
      colnames(x) <- c('Treatmeant','Control')
      rownames(x) <- paste0('Contrast group',1:nrow(x))
      x
    },options = list(
      columnDefs = list(list(className = 'dt-center', targets = "_all")))
    )

    ##---------------DE genes-----------------
    observeEvent(input$run.DE,{
      cat('\nDE gene analysis')
      withProgress(message = 'DE gene analysis...', detail = 'This may take a while...',
                   value = 0, {
                     if(is.null(DDD.data$genes_batch)){
                       batch.effect <- NULL
                     } else {
                       batch.effect <- DDD.data$genes_batch$W}
                     incProgress(0.1)
                     design <- condition2design(condition = DDD.data$samples_new$condition,
                                                batch.effect = batch.effect)
                     DE.pipeline <- input$DE.pipeline
                     incProgress(0.2)
                     switch(DE.pipeline,
                            limma={
                              genes_DDD_stat <- limma.pipeline(dge = DDD.data$genes_dge,
                                                               design = design,
                                                               deltaPS = DDD.data$deltaPS,
                                                               contrast = DDD.data$contrast,
                                                               diffAS = F,
                                                               adjust.method = input$p.adjust.method)
                            },
                            glmQL={
                              genes_DDD_stat <- edgeR.pipeline(dge = DDD.data$genes_dge,
                                                               design = design,
                                                               deltaPS = DDD.data$deltaPS,
                                                               contrast = DDD.data$contrast,
                                                               diffAS = F,
                                                               method = 'glmQL',
                                                               adjust.method = input$p.adjust.method)
                            },
                            glm={
                              genes_DDD_stat <- edgeR.pipeline(dge = DDD.data$genes_dge,
                                                               design = design,
                                                               deltaPS = NULL,
                                                               contrast = DDD.data$contrast,
                                                               diffAS = F,
                                                               method = 'glm',
                                                               adjust.method = input$p.adjust.method)
                            }
                     )
                     incProgress(0.7)
                     DE_genes <- summary.DE.target(stat = genes_DDD_stat$DE.stat,
                                                   cutoff = c(adj.pval=input$pval.cutoff,
                                                              log2FC=input$lfc.cutoff))

                     DDD.data$DE_genes <- DE_genes
                     DDD.data$genes_log2FC <- genes_DDD_stat$DE.lfc
                     save(DE_genes,file=paste0(DDD.data$data.folder,'/DE_genes.RData'))
                     incProgress(0.8)
                     save(genes_DDD_stat,file=paste0(DDD.data$data.folder,'/genes_DDD_stat.RData'))
                     message(paste0('genes_DDD_stat.RData is saved in folder: ',DDD.data$data.folder))
                   })
    })


    ###run DE information
    run.DE.info <- reactive({
      if(is.null(DDD.data$data.folder))
        return(NULL)
      if(input$run.DE>0){
        paste0('genes_DDD_stat.RData is saved in the data folder.')
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/genes_DDD_stat.RData'))){
          paste0('genes_DDD_stat.RData is already in "data" folder.\nGenerate new?.')
        } else {
          paste0('genes_DDD_stat.RData is not in "data" folder.\nPlease generate new.')
        }
      }
    })

    output$run.DE.info <- renderText({
      if(is.null(run.DE.info()))
        return(NULL)
      run.DE.info()
    })
    #

    ##---------------generate deltaPS-----------------
    observeEvent(input$run.deltaPS,{
      if(is.null(DDD.data$PS)){
        ##load require data for analysis
        if(is.null(DDD.data$target_high)){
          load(paste0(DDD.data$data.folder,'/target_high.RData'))
          DDD.data$target_high <- target_high
        }

        if(is.null(DDD.data$mapping)){
          load(paste0(DDD.data$data.folder,'/mapping.RData'))
          DDD.data$mapping <- mapping
        }

        if(input$PS.data.type=='TPM'){
          if(is.null(DDD.data$samples)){
            load(paste0(DDD.data$data.folder,'/samples.RData'))
            DDD.data$samples <- samples
          }

          if(is.null(DDD.data$trans_TPM)){
            load(paste0(DDD.data$data.folder,'/trans_TPM.RData'))
            DDD.data$trans_TPM <- trans_TPM
          }
          transAbundance <- DDD.data$trans_TPM
          condition = DDD.data$samples$condition
        }
        if(input$PS.data.type=='Read counts'){
          if(is.null(DDD.data$samples_new)){
            load(paste0(DDD.data$data.folder,'/samples_new.RData'))
            DDD.data$samples_new <- samples_new
          }
          if(is.null(DDD.data$trans_counts)){
            load(paste0(DDD.data$data.folder,'/trans_counts.RData'))
            DDD.data$trans_counts <- trans_counts
          }
          transAbundance <- DDD.data$trans_counts
          condition = DDD.data$samples_new$condition
        }
        transAbundance <- transAbundance[DDD.data$target_high$trans_high,]
      } else {
        transAbundance <- NULL
        condition <- NULL
      }

      cat('\nGenerating deltaPS')
      deltaPS <- transAbundance2PS(transAbundance = transAbundance,
                                   PS = DDD.data$PS,
                                   contrast = DDD.data$contrast,
                                   condition = condition,
                                   mapping = DDD.data$mapping[DDD.data$target_high$trans_high,])
      PS <- deltaPS$PS
      deltaPS <- deltaPS$deltaPS
      save(deltaPS,file=paste0(DDD.data$data.folder,'/deltaPS.RData'))
      save(PS,file=paste0(DDD.data$data.folder,'/PS.RData'))
      message(paste0('deltaPS.RData are saved in folder: ',DDD.data$data.folder))
      showNotification(paste0('deltaPS.RData are saved in folder: ',DDD.data$data.folder))
      message('Done!!!')
      DDD.data$PS <- PS
      DDD.data$deltaPS <- deltaPS
    })

    output$deltaPS.table <- DT::renderDataTable({
      if(is.null(DDD.data$deltaPS))
        return(NULL)
      round(DDD.data$deltaPS[1:20,,drop=F],5)
    },options = list(
      columnDefs = list(list(className = 'dt-center', targets = "_all"))))

    observeEvent(input$save.deltaPS.table,{
      write.csv(DDD.data$deltaPS,file=paste0(DDD.data$result.folder,'/deltaPS.csv'),row.names = T)
      message(paste0('Table is saved in folder: ',DDD.data$result.folder))
      showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
    })

    ### message of run deltaPS
    rundeltaPSinfo <- reactive({
      if(is.null(DDD.data$data.folder))
        return(NULL)
      if(input$run.deltaPS>0){
        paste0('deltaPS.RData is saved in the data folder.')
      } else {
        if(!is.null(DDD.data$deltaPS)){
          paste0('deltaPS.RData is already uploaded in the workspace.')
        } else {
          if(file.exists(paste0(DDD.data$data.folder,'/deltaPS.RData'))){
            paste0('deltaPS.RData is already in "data" folder.\nCan upload directly in "Data generation".')
          } else {
            paste0('deltaPS.RData is not in "data" folder.\nPlease generate new.')
          }
        }
      }
    })

    output$rundeltaPSinfo <- renderText({
      if(is.null(rundeltaPSinfo()))
        return(NULL)
      rundeltaPSinfo()
    })

    ##--------------->> DE DAS DTU transcripts <<---------------
    observeEvent(input$run.diffsplice,{
      if(is.null(DDD.data$deltaPS)){
        message('Please generate deltaPS values.')
        showNotification('Please generate deltaPS values.')
        return(NULL)
      }
      withProgress(message = 'DE gene analysis...', detail = 'This may take a while...',
                   value = 0, {
                     cat('\nDAS gene, DE and DTU transcript analysis')
                     if(is.null(DDD.data$trans_batch)){
                       batch.effect <- NULL
                     } else {
                       batch.effect <- DDD.data$trans_batch$W
                     }
                     incProgress(0.1)
                     design <- condition2design(condition = DDD.data$samples_new$condition,
                                                batch.effect = batch.effect)
                     DE.pipeline <- input$DE.pipeline
                     incProgress(0.2)
                     switch(DE.pipeline,
                            limma={
                              trans_DDD_stat <- limma.pipeline(dge = DDD.data$trans_dge,
                                                               design = design,
                                                               deltaPS = DDD.data$deltaPS,
                                                               contrast = DDD.data$contrast,
                                                               diffAS = T,
                                                               adjust.method = input$p.adjust.method)
                            },
                            glmQL={
                              trans_DDD_stat <- edgeR.pipeline(dge = DDD.data$trans_dge,
                                                               design = design,
                                                               deltaPS = DDD.data$deltaPS,
                                                               contrast = DDD.data$contrast,
                                                               diffAS = T,
                                                               method = 'glmQL',
                                                               adjust.method = input$p.adjust.method)
                            },
                            glm={
                              trans_DDD_stat <- edgeR.pipeline(dge = DDD.data$trans_dge,
                                                               design = design,
                                                               deltaPS = DDD.data$deltaPS,
                                                               contrast = DDD.data$contrast,
                                                               diffAS = T,
                                                               method = 'glm',
                                                               adjust.method = input$p.adjust.method)
                            }
                     )
                     incProgress(0.7)
                     ##DE trans
                     DE_trans <- summary.DE.target(stat = trans_DDD_stat$DE.stat,
                                                   cutoff = c(adj.pval=input$pval.cutoff,
                                                              log2FC=input$lfc.cutoff))
                     DDD.data$DE_trans <- DE_trans
                     DDD.data$trans_log2FC <- trans_DDD_stat$DE.lfc
                     save(DE_trans,file=paste0(DDD.data$data.folder,'/DE_trans.RData'))

                     ##DAS genes
                     if(input$DAS.p.method=='F-test')
                       DAS.stat <- trans_DDD_stat$DAS.F.stat else DAS.stat <- trans_DDD_stat$DAS.Simes.stat

                     lfc <- DDD.data$genes_log2FC
                     lfc <- reshape2::melt(as.matrix(lfc))
                     colnames(lfc) <- c('target','contrast','log2FC')
                     DAS_genes <- summary.DAS.target(stat = DAS.stat,
                                                     lfc = lfc,
                                                     cutoff=c(input$pval.cutoff,input$deltaPS.cutoff))
                     DDD.data$DAS_genes <- DAS_genes
                     save(DAS_genes,file=paste0(DDD.data$data.folder,'/DAS_genes.RData'))

                     ##DTU  trans
                     lfc <- DDD.data$trans_log2FC
                     lfc <- reshape2::melt(as.matrix(lfc))
                     colnames(lfc) <- c('target','contrast','log2FC')
                     DTU_trans <- summary.DAS.target(stat = trans_DDD_stat$DTU.stat,
                                                     lfc = lfc,cutoff = c(adj.pval=input$pval.cutoff,
                                                                          deltaPS=input$deltaPS.cutoff))
                     DDD.data$DTU_trans <- DTU_trans
                     save(DTU_trans,file=paste0(DDD.data$data.folder,'/DTU_trans.RData'))
                     incProgress(0.8)

                     save(trans_DDD_stat,file=paste0(DDD.data$data.folder,'/trans_DDD_stat.RData'))
                     DDD.data$trans_DDD_stat <- trans_DDD_stat
                     message(paste0('trans_DDD_stat.RData is saved in folder: ',DDD.data$data.folder))

                     ######

                   })
    })

    diffsplice.method <- reactive({
      if(input$DE.pipeline=='limma')
        'Use limma::diffSplice' else 'Use edgeR::diffSpliceDGE'
    })

    output$diffsplice.method <- renderText({
      diffsplice.method()
    })

    ###run DAS information
    run.DAS.info <- reactive({
      if(is.null(DDD.data$data.folder))
        return(NULL)
      if(input$run.DE>0){
        paste0('trans_DDD_stat.RData is saved in the data folder.')
      } else {
        if(file.exists(paste0(DDD.data$data.folder,'/trans_DDD_stat.RData'))){
          paste0('trans_DDD_stat.RData is already in "data" folder.\nGenerate new?.')
        } else {
          paste0('trans_DDD_stat.RData is not in "data" folder.\nPlease generate new.')
        }
      }
    })

    output$run.DAS.info <- renderText({
      if(is.null(run.DAS.info()))
        return(NULL)
      run.DAS.info()
    })

    ##############################################################
    ####                       Result summary                 ####
    ##############################################################
    ##---------------->> DE DAS DTU statistics <<-----------------
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

      stat <- x[[input$contrast.groups]][1:input$top.stat,]
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

    ##---------------DDD bar plot-----------------
    # output$up.down.bar.plot <- renderPlotly({
    g.updown <- reactive({
      if(is.null((DDD.data$DE_genes)) | is.null((DDD.data$DAS_genes)) | is.null((DDD.data$DE_trans)) | is.null((DDD.data$DTU_trans)))
        return(NULL)
      idx <- c('DE_genes','DAS_genes','DE_trans','DTU_trans')
      title.idx <- gsub('_',' ',idx)
      title.idx <- gsub('trans','transcripts',title.idx)
      names(title.idx) <- idx
      g <- lapply(idx,function(i){
        idx <- factor(DDD.data[[i]]$contrast,levels = DDD.data$contrast)
        targets <-  split(DDD.data[[i]],idx)
        data2plot <- lapply(DDD.data$contrast,function(i){
          if(nrow(targets[[i]])==0){
            x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
          } else {
            x <- data.frame(contrast=i,table(targets[[i]]$up.down))
            colnames(x) <- c('contrast','regulation','number')
          }
          x
        })
        data2plot <- do.call(rbind,data2plot)
        plot.updown(data2plot,plot.title = title.idx[i],contrast = DDD.data$contrast)
      })
      names(g) <- title.idx
      g
    })

    output$up.down.bar.plot.DE.genes <- renderPlot({
      if(is.null((g.updown())))
        return(NULL)
      print(g.updown()[[1]])
    })

    output$up.down.bar.plot.DAS.genes <- renderPlot({
      if(is.null((g.updown())))
        return(NULL)
      print(g.updown()[[2]])
    })

    output$up.down.bar.plot.DE.trans <- renderPlot({
      if(is.null((g.updown())))
        return(NULL)
      print(g.updown()[[3]])
    })

    output$up.down.bar.plot.DTU.trans <- renderPlot({
      if(is.null((g.updown())))
        return(NULL)
      print(g.updown()[[4]])
    })

    observeEvent(input$save.up.down.bar.plot,{
      lapply(names(g.updown()),function(i){
        png(paste0(DDD.data$figure.folder,'/',i,' updown regulation numbers.png'),
            width = length(DDD.data$contrast)*5/2.54,10/2.54,units = 'in',res = 300)
        print(g.updown()[[i]])
        dev.off()

        pdf(paste0(DDD.data$figure.folder,'/',i,' updown regulation numbers.pdf'),
            width = length(DDD.data$contrast)*5/2.54,10/2.54)
        print(g.updown()[[i]])
        dev.off()
      })

      message(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
    })

    ##---------------consensus across contrast-----------------
    output$across.contrast.euler <- renderUI({
      span(
        selectInput(inputId = 'across.contrast.group',
                    label = 'Contrast groups (multiple selection allowed)',
                    selected = DDD.data$contrast[1:2],multiple = T,
                    choices = DDD.data$contrast)
      )
    })

    output$across.target.euler <- renderUI({
      span(
        selectInput(inputId = 'across.target',
                    label = 'Contrast groups',
                    selected = DDD.data$contrast[1],multiple = F,
                    choices = DDD.data$contrast)
      )
    })


    ##---------------flow chat-----------------
    genes.flow.chart <- reactive({
      DE.genes <- unique(DDD.data$DE_genes$target)
      DAS.genes <- unique(DDD.data$DAS_genes$target)
      genes.flow.chart <- function(){
        plot.flow.chart(expressed = DDD.data$target_high$genes_high,
                        x = DE.genes,
                        y = DAS.genes,
                        type = 'genes',
                        pval.cutoff = input$pval.cutoff,lfc.cutoff = input$lfc.cutoff,deltaPS.cutoff = input$deltaPS.cutoff)
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
        plot.flow.chart(expressed = DDD.data$target_high$trans_high,
                        x = DE.trans,
                        y = DTU.trans,
                        type = 'transcripts',
                        pval.cutoff = input$pval.cutoff,lfc.cutoff = input$lfc.cutoff,deltaPS.cutoff = input$deltaPS.cutoff)
      }
      trans.flow.chart
    })
    output$trans.flow.chart <- renderPlot({
      trans.flow.chart()()
    })

    observeEvent(input$save.union.flow.chart,{
      ##transcript level
      png(filename = paste0(DDD.data$figure.folder,'/Union set DE genes vs DAS genes.png'),
          width = 22/2.54,height = 13/2.54,units = 'in',res = 300)
      genes.flow.chart()()
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Union set DE genes vs DAS genes.pdf'),
          width = 22/2.54,height = 13/2.54)
      genes.flow.chart()()
      dev.off()

      ###gene level
      # graphics.off()
      png(filename = paste0(DDD.data$figure.folder,'/Union set DE transcripts vs DTU transcripts.png'),
          width = 22/2.54,height = 13/2.54,units = 'in',res = 300)
      trans.flow.chart()()
      dev.off()

      pdf(file = paste0(DDD.data$figure.folder,'/Union set DE transcripts vs DTU transcripts.pdf'),
          width = 22/2.54,height = 13/2.54)
      trans.flow.chart()()
      dev.off()

      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
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
        g <- plot.euler.diagram(x = targets)
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
      grid.arrange(g.across.contrast()[[3]],top=textGrob('DAS genes', gp=gpar(cex=1.2)))
    })

    ###---DTU transcripts
    output$DTU.trans.euler <- renderPlot({
      if(is.null(g.across.contrast()))
        return(NULL)
      grid.arrange(g.across.contrast()[[4]],top=textGrob('DAS genes', gp=gpar(cex=1.2)))
    })

    observeEvent(input$save.across.contrast.euler.plot,{
      ###save DE.genes

      lapply(names(g.across.contrast()),function(i){
        figure.name <- paste0(DDD.data$figure.folder,'/',i,' euler plot across contrast ',
                              paste0(input$across.contrast.group,collapse = '&'))
        png(paste0(figure.name,'.png'),
            width = 15/2.54,height = 12/2.54,units = 'in',res = 300)
        grid.arrange(g.across.contrast()[[i]],top=textGrob('DE genes', gp=gpar(cex=1.2)))
        dev.off()
        pdf(paste0(figure.name,'.pdf'),
            width = 15/2.54,height = 12/2.54)
        grid.arrange(g.across.contrast()[[i]],top=textGrob('DE genes', gp=gpar(cex=1.2)))
        dev.off()
      })

      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))

    })

    ###target numbers
    DDD.numbers <- reactive({
      if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes) | is.null(DDD.data$DE_trans)| is.null(DDD.data$DTU_trans))
        return(NULL)
      summary.DDD.number(DE_genes = DDD.data$DE_genes,
                         DAS_genes = DDD.data$DAS_genes,
                         DE_trans = DDD.data$DE_trans,
                         DTU_trans=DDD.data$DTU_trans,contrast = DDD.data$contrast)
    })

    output$DDD.numbers <- renderTable({
      DDD.numbers()
    })

    ##---------------DE vs DAS-----------------
    DEvsDAS <- reactive({
      if(is.null(DDD.data$DE_genes) | is.null(DDD.data$DAS_genes))
        return(NULL)
      # summary.DDD.vs(x = split(DDD.data$DE_genes$target,DDD.data$DE_genes$contrast),
      #                y = split(DDD.data$DAS_genes$target,DDD.data$DAS_genes$contrast),
      #                contrast = DDD.data$contrast,
      #                idx = c('DEonly','DE&DAS','DASonly'))
      DEvsDAS(DE_genes = DDD.data$DE_genes,
              DAS_genes = DDD.data$DAS_genes,
              contrast = DDD.data$contrast)

    })

    output$DE.vs.DAS <- renderTable({
      if(is.null(DEvsDAS()))
        return(NULL)
      DEvsDAS()
    },align='c')

    ##euler plot
    ###---DE vs DAS
    g.across.target1 <- reactive({
      if(is.null(DEvsDAS())| is.null(input$across.target))
        return(NULL)
      x <- unlist(DEvsDAS()[DEvsDAS()$Contrast==input$across.target,-1])
      if(length(x)==0)
        return(NULL)
      names(x) <- c('DE','DE&DAS','DAS')
      g <- plot.euler.diagram(x = x,fill = gg.color.hue(2))
      g
    })

    output$DEvsDAS.euler <- renderPlot({
      if(is.null(g.across.target1()))
        return(NULL)
      grid.arrange(g.across.target1(),top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
    })


    ##---------------DE vs DTU-----------------
    DEvsDTU <- reactive({
      if(is.null(DDD.data$DE_trans) | is.null(DDD.data$DTU_trans))
        return(NULL)
      # summary.DDD.vs(x = split(DDD.data$DE_trans$target,DDD.data$DE_trans$contrast),
      #                y = split(DDD.data$DTU_trans$target,DDD.data$DTU_trans$contrast),
      #                contrast = DDD.data$contrast,
      #                idx = c('DEonly','DE&DTU','DTUonly'))
      DEvsDTU(DE_trans = DDD.data$DE_trans,
              DTU_trans = DDD.data$DTU_trans,
              contrast = DDD.data$contrast)
    })

    output$DE.vs.DTU <- renderTable({
      if(is.null(DEvsDTU()))
        return(NULL)
      DEvsDTU()
    },align='c')

    ###---DE vs DTU
    g.across.target2 <- reactive({
      if(is.null(DEvsDTU())| is.null(input$across.target))
        return(NULL)
      x <- unlist(DEvsDTU()[DEvsDTU()$Contrast==input$across.target,-1])
      if(length(x)==0)
        return(NULL)
      names(x) <- c('DE','DE&DTU','DTU')
      g <- plot.euler.diagram(x = x,fill = gg.color.hue(2))
      g
    })

    output$DEvsDTU.euler <- renderPlot({
      if(is.null(g.across.target2()))
        return(NULL)
      grid.arrange(g.across.target2(),top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
    })

    observeEvent(input$save.across.target.euler.plot,{
      ##DE vs DAS genes
      figure.name <- paste0(DDD.data$figure.folder,'/DE vs DAS gene euler plot in contrast ',input$across.target)
      png(paste0(figure.name,'.png'),
          width = 15/2.54,height = 12/2.54,units = 'in',res = 300)
      grid.arrange(g.across.target1(),top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
      dev.off()
      pdf(paste0(figure.name,'.pdf'),
          width = 15/2.54,height = 12/2.54)
      grid.arrange(g.across.target1(),top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
      dev.off()

      ##DE vs DTU transcript
      figure.name <- paste0(DDD.data$figure.folder,'/DE vs DTU transcript euler plot in contrast ',input$across.target)
      png(paste0(figure.name,'.png'),
          width = 15/2.54,height = 12/2.54,units = 'in',res = 300)
      grid.arrange(g.across.target2(),top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
      dev.off()

      pdf(paste0(figure.name,'.pdf'),
          width = 15/2.54,height = 12/2.54)
      grid.arrange(g.across.target2(),top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
      dev.off()
      message(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figures are saved in folder: ',DDD.data$figure.folder))
    })


    ###---save stat tables
    observeEvent(input$save.stat,{
      withProgress(message = 'Save test statistics',
                   detail = 'This may take a while...', value = 0, {
                     incProgress(0.1)
                     write.csv(DDD.data$DE_genes,
                               file=paste0(DDD.data$result.folder,'/DE genes.csv'),row.names = F)
                     incProgress(0.3)
                     write.csv(DDD.data$DAS_genes,
                               file=paste0(DDD.data$result.folder,'/DAS genes.csv'),row.names = F)
                     incProgress(0.5)
                     write.csv(DDD.data$DE_trans,
                               file=paste0(DDD.data$result.folder,'/DE transcripts.csv'),row.names = F)
                     incProgress(0.7)
                     write.csv(DDD.data$DTU_trans,
                               file=paste0(DDD.data$result.folder,'/DTU transcripts.csv'),row.names = F)
                     incProgress(0.8)
                     write.csv(DDD.numbers(),file=paste0(DDD.data$result.folder,'/DE DAS DTU numbers.csv'),
                               row.names = F)
                     incProgress(0.9)
                     write.csv(DEvsDAS(),file=paste0(DDD.data$result.folder,'/DE vs DAS gene number.csv'),
                               row.names = F)
                     incProgress(1)
                     write.csv(DEvsDTU(),file=paste0(DDD.data$result.folder,'/DE vs DTU transcript number.csv'),
                               row.names = F)
                   })
      message(paste0('All tables of test statistics and target numbers are saved in folder: ',DDD.data$result.folder))
      showNotification(paste0('All tables of test statistics are saved in folder: ',DDD.data$result.folder))
    })


    ##############################################################
    ######                    Advanced plot                   ####
    ##############################################################

    ##---------------->>        heatmap        <<-----------------
    observe({
      shinyjs::toggleState(id = 'heatmap.targetlist.type',
                           condition = input$heatmap.select.or.upload=='Upload target list')
      shinyjs::toggleState(id = 'heatmap.target.list.input',
                           condition = input$heatmap.select.or.upload=='Upload target list')
      shinyjs::toggleState(id = 'heatmap.target.type',
                           condition = input$heatmap.select.or.upload=='Select targets')
    })

    heatmap.targets <- eventReactive(input$plot.heatmap,{
      if(is.null(DDD.data$genes_TPM)){
        load(paste0(DDD.data$data.folder,'/genes_TPM.RData'))
        DDD.data$genes_TPM <- genes_TPM
      }
      if(is.null(DDD.data$trans_TPM)){
        load(paste0(DDD.data$data.folder,'/trans_TPM.RData'))
        DDD.data$trans_TPM <- trans_TPM
      }

      if(input$heatmap.select.or.upload=='Select targets'){
        if(input$heatmap.target.type=='DE genes'){
          targets <- unique(DDD.data$DE_genes$target)
          data2heatmap <- DDD.data$genes_TPM[targets,]
          column_title <- paste0(length(targets),' DE genes')
        }

        if(input$heatmap.target.type=='DAS genes'){
          targets <- intersect(unique(DDD.data$DAS_genes$target),unique(DDD.data$DE_genes$target))
          data2heatmap <- DDD.data$genes_TPM[targets,]
          column_title <- paste0(length(targets),' DE&DAS genes (DAS only are filtered)')
        }

        if(input$heatmap.target.type=='DE transcripts'){
          targets <- unique(DDD.data$DE_trans$target)
          data2heatmap <- DDD.data$trans_TPM[targets,]
          column_title <- paste0(length(targets),' DE transcripts')
        }

        if(input$heatmap.target.type=='DTU transcripts'){
          targets <- intersect(unique(DDD.data$DTU_trans$target),unique(DDD.data$DE_trans$target))
          data2heatmap <- DDD.data$trans_TPM[targets,]
          column_title <- paste0(length(targets),' DE&DTU transcripts (DTU only are filtered)')
        }
      }

      if(input$heatmap.select.or.upload=='Upload target list'){
        inFile <- input$heatmap.target.list.input
        if (is.null(inFile))
          return(NULL)
        targets <- read.csv(inFile$datapath, header = F)
        targets <- as.vector(t(targets))
        if(input$heatmap.targetlist.type=='Gene list'){
          data2heatmap <- DDD.data$genes_TPM[targets,]
          column_title <- paste0(length(targets),' Genes')
        } else {
          data2heatmap <- DDD.data$trans_TPM[targets,]
          column_title <- paste0(length(targets),' Transcripts')
        }
      }
      list(data2heatmap=data2heatmap,column_title=column_title)
    })


    heatmap.g <- eventReactive(input$plot.heatmap,{
      withProgress(message = 'Generate gene expression...', value = 0, {
        incProgress(0.1)
        data2plot <- rowmean(x = t(heatmap.targets()$data2heatmap),
                             group = DDD.data$samples$condition,
                             reorder = F)
        data2plot <- t(scale(data2plot))
        incProgress(0.2)
        hc.dist <- dist(data2plot,method = input$dist.method)
        hc <- fastcluster::hclust(hc.dist,method = input$cluster.method)
        clusters <- cutree(hc, k = input$cluster.number)
        clusters <- reorder.clusters(clusters = clusters,dat = data2plot)
        incProgress(0.3)
        g <- Heatmap(as.matrix(data2plot), name = 'Z-scores',
                     cluster_rows = TRUE,
                     clustering_method_rows=input$cluster.method,
                     row_dend_reorder = T,
                     show_row_names = FALSE,
                     show_column_names = ifelse(ncol(data2plot)>10,F,T),
                     cluster_columns = FALSE,
                     split=clusters,
                     column_title= heatmap.targets()$column_title)
        incProgress(0.8)
      })
      list(g=g,clusters=clusters)
    })

    output$heatmap.plot.panel <- renderPlot({
      draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
    })

    observeEvent(input$save.target.in.cluster,{
      clusters <- heatmap.g()$clusters
      x <- split(names(clusters),clusters)
      x <- lapply(names(x),function(i){
        data.frame(Clusters=i,Targets=x[[i]])
      })
      x <- do.call(rbind,x)
      colnames(x) <- c('Clusters','Targets')
      write.csv(x,file=paste0(DDD.data$result.folder,'/Target in each cluster heatmap ', heatmap.targets()$column_title,'.csv'),row.names = F)
      message(paste0('Table is saved in folder: ',DDD.data$result.folder))
      showNotification(paste0('Table is saved in folder: ',DDD.data$result.folder))
    })

    observeEvent(input$save.heatmap,{
      withProgress(message = 'Saving heatmap...', detail = 'This may take a while...',
                   value = 0.5, {
                     ##DE vs DTU transcript
                     figure.name <- paste0(DDD.data$figure.folder,'/Heatmap ', input$heatmap.target.type)
                     png(paste0(figure.name,'.png'),
                         width = pmax(10,1*length(unique(DDD.data$samples$condition)))/2.54,height = 20/2.54,units = 'in',res = 300)
                     draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                     dev.off()

                     pdf(paste0(figure.name,'.pdf'),
                         width = pmax(10,1*length(unique(DDD.data$samples$condition)))/2.54,height = 20/2.54)
                     draw(heatmap.g()$g,column_title='Conditions',column_title_side = "bottom")
                     dev.off()
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

      g.pr <- plot.abundance(data.exp = data.exp,
                             gene = gene,
                             mapping = mapping,
                             genes.ann = DDD.data$genes.ann,
                             trans.ann = DDD.data$trans.ann,
                             trans.expressed = NULL,
                             reps = reps,
                             y.lab = input$profile.data.type)

      g.ps <- plot.PS(data.exp = data.exp,
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

    observeEvent(input$save.abundance.plot,{
      n <- length(unique(g.profiles()$g.pr$data))-1
      gene <- g.profiles()$gene
      folder2save <- paste0(DDD.data$figure.folder,'/Abundance plot')
      if(!file.exists(folder2save))
        dir.create(folder2save,recursive = T)

      png(paste0(folder2save,'/Abundance ',gene,'.png'),
          width = (16+floor(n/15))/2.54,
          height = 10/2.54,
          units = 'in',res = 300)
      print(g.profiles()$g.pr)
      dev.off()

      pdf(paste0(folder2save,'/Abundance ',gene,'.pdf'),
          width = (16+floor(n/15))/2.54,
          height = 10/2.54)
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
          width = (16+floor(n/15))/2.54,
          height = 10/2.54,
          units = 'in',res = 300)
      print(g.profiles()$g.ps)
      dev.off()

      pdf(paste0(folder2save,'/PS ',gene,'.pdf'),
          width = (16+floor(n/15))/2.54,
          height = 10/2.54)
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
      folder2Abundance <- paste0(DDD.data$figure.folder,'/Abundance plot')
      if(!file.exists(folder2Abundance))
        dir.create(folder2Abundance,recursive = T)

      folder2PS <- paste0(DDD.data$figure.folder,'/PS plot')
      if(!file.exists(folder2PS))
        dir.create(folder2PS,recursive = T)

      withProgress(message = paste0('Plotting ',length(genes),' genes'),
                   detail = 'This may take a while...', value = 0, {
                     for(gene in genes){
                       incProgress(1/length(genes))
                       if(input$multiple.plot.type %in% c('Abundance','both')){
                         g.pr <- plot.abundance(data.exp = data.exp,
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
                               width = (16+floor(n/15))/2.54,
                               height = 10/2.54,
                               units = 'in',res = 300)
                           print(g.pr)
                           dev.off()
                         }

                         if(input$multi.plot.format %in% c('pdf','both')){
                           pdf(paste0(folder2Abundance,'/Abundance ',gene,'.pdf'),
                               width = (16+floor(n/15))/2.54,
                               height = 10/2.54)
                           print(g.pr)
                           dev.off()
                         }
                       }

                       if(input$multiple.plot.type %in% c('PS','both')){
                         g.ps <- plot.PS(data.exp = data.exp,
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
                               width = (16+floor(n/15))/2.54,
                               height = 10/2.54,
                               units = 'in',res = 300)
                           print(g.ps)
                           dev.off()
                         }
                         if(input$multi.plot.format %in% c('pdf','both')){
                           pdf(paste0(folder2PS,'/PS ',gene,'.pdf'),
                               width = (16+floor(n/15))/2.54,
                               height = 10/2.54)
                           print(g.ps)
                           dev.off()
                         }
                       }

                     }
                   })




      if(input$multiple.plot.type %in% c('PS','both')){
        g.ps <- plot.PS(data.exp = data.exp,
                        gene = gene,
                        mapping = mapping,
                        genes.ann = DDD.data$genes.ann,
                        trans.ann = DDD.data$trans.ann,
                        trans.expressed = NULL,
                        reps = reps,
                        y.lab = 'PS')

      }
    })

    ##---------------GO plot----------------
    observeEvent(input$generate.target.list,{
      write.csv(data.frame(`DE genes`=unique(DDD.data$DE_genes$target)),
                file = paste0(DDD.data$result.folder,'/DE gene list union set of all contrast group.csv'),
                row.names = F)
      write.csv(data.frame(`DAS genes`=unique(DDD.data$DAS_genes$target)),
                file = paste0(DDD.data$result.folder,'/DAS gene list union set of all contrast group.csv'),
                row.names = F)
      write.csv(data.frame(`DE transcripts`=unique(DDD.data$DE_trans$target)),
                file = paste0(DDD.data$result.folder,'/DE transcript list union set of all contrast group.csv'),
                row.names = F)
      write.csv(data.frame(`DTU transcripts`=unique(DDD.data$DTU_trans$target)),
                file = paste0(DDD.data$result.folder,'/DTU transcript list union set of all contrast group.csv'),
                row.names = F)
      message(paste0('Tables of target list save in: ',DDD.data$result.folder))
      showNotification(paste0('Tables of target list save in: ',DDD.data$result.folder))
    })

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



    output$show.go.table.column <- renderUI({
      if(is.null(go.table()))
        return(NULL)
      span(
        selectInput(inputId = 'select.go.table.column',label = 'Select column of statistics',
                    choices = colnames(go.table())[-c(1:2)])
      )
    })

    data2goplot <- eventReactive(input$make.go.plot,{
      col.idx<- input$select.go.table.column
      data2plot <- go.table()[,c(1,2,which(colnames(go.table())==col.idx))]
      colnames(data2plot) <- c('Category','Term','Value')
      data2plot <- by(data2plot,data2plot$Category,function(x){
        x[order(x$Value,decreasing = T),]
      })
      data2plot <- do.call(rbind,data2plot)
      idx<-sapply(data2plot$Term,function(x){
        x<-as.vector(x)
        if(nchar(x)>100)
          x<-paste0(substr(x,1,100),'...')
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
          labs(y=input$select.go.table.column,title=paste0('GO annotation ',input$go.plot.label))
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

    observeEvent(input$save.go.plot,{
      if(is.null(data2goplot()))
        return(NULL)
      col.idx <- input$select.go.table.column
      data2plot <- go.table()[,c(1,2,which(colnames(go.table())==col.idx))]
      colnames(data2plot) <- c('Category','Term','Value')
      data2plot <- by(data2plot,data2plot$Category,function(x){
        x[order(x$Value,decreasing = T),]
      })
      data2plot <- do.call(rbind,data2plot)
      data2plot$Term <- factor(data2plot$Term,levels = rev(data2plot$Term))

      g <- ggplot(data2plot,aes(x=Term,y=Value))+
        geom_bar(stat='identity',aes(fill=Category))+
        coord_flip() +
        theme_bw()+
        theme(axis.title.y = element_blank())+
        facet_grid(Category~.,scales = 'free_y',space='free_y')+
        labs(y=col.idx,title=paste0('GO annotation ',input$go.plot.label))
      png(paste0(DDD.data$figure.folder,'/',input$go.plot.label,' GO annotation plot.png'),
          width = 20/2.54,height = pmax(10,row.height()*0.4)/2.54,
          units = 'in',res = 300)
      print(g)
      dev.off()
      pdf(paste0(DDD.data$figure.folder,'/',input$go.plot.label,' GO annotation plot.pdf'),
          width = 20/2.54,height = pmax(10,row.height()*0.4)/2.54)
      print(g)
      dev.off()
      message(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
      showNotification(paste0('Figure is saved in folder: ',DDD.data$figure.folder))
    })

    ###########################################################################
    #======================= Generate report ============================


    data.para <- reactive({
      x <- rbind(
        c('Folders','Data',ifelse(is.null(DDD.data$data.folder),'No information',DDD.data$data.folder)),
        c('','Results',ifelse(is.null(DDD.data$result.folder),'No information',DDD.data$result.folder)),
        c('','Figures',ifelse(is.null(DDD.data$figure.folder),'No information',DDD.data$figure.folder)),
        c('Data generation','tximport method',input$tximport.method),
        c('Data pre-processing','Merge sequencing replicates',
          ifelse(is.null(DDD.data$data.info),'No information',
                 ifelse(DDD.data$data.info$sample[1] > DDD.data$data.info$merged.sample[1],'Yes','No'))),
        c('','Low expression CPM cut-off',input$cpm.cut),
        c('','Sample number for CPM cut-off', input$sample.n.cut),
        c('','Batch effect estimation',ifelse(is.null(DDD.data$data.folder),'No information',
                                              ifelse(length(list.files(DDD.data$data.folder,'*batch.RData'))>0,'Yes','No'))),
        # ifelse(length(list.files(DDD.data$data.folder,'*batch.RData'))>0,'Yes','No'))),
        c('','Batch effect estimation method',input$ruvseq.method),
        c('','Normalisation method',input$norm.method),
        c('DE DAS and DTU','Pipeline',input$DE.pipeline),
        c('','AS function',ifelse(input$DE.pipeline=='limma',
                                  'limma::diffSplice','edgeR::diffSpliceDGE')),
        c('','P-value adjust method',input$p.adjust.method),
        c('','Adjusted p-value cut-off',input$pval.cutoff),
        c('','Log2 fold change cut-off',input$lfc.cutoff),
        c('','delta PS cut-off',input$deltaPS.cutoff)
      )
      x <- data.frame(x)
      colnames(x) <- c('Step','Description','Parameter (Double click to edit if not correct)')
      x
    })

    x = reactiveValues(df = NULL)

    observe({
      df <- data.para()
      x$df <- df
    })


    output$para.summary <- DT::renderDT({
      x$df
    },editable = T,options = list(pageLength = 50))


    proxy = DT::dataTableProxy('para.summary')
    observeEvent(input$para.summary_cell_edit, {
      info = input$para.summary_cell_edit
      i = info$row
      j = info$col
      v = info$value

      # problem starts here
      x$df[i, j] <- isolate(suppressWarnings(DT::coerceValue(v, x$df[i, j])))
    })

    observeEvent(input$save.para.summary,{
      write.csv(x$df,file=paste0(DDD.data$result.folder,'/Parameter summary.csv'),row.names = F)
      message(paste0('Table of parameters is save in: ',DDD.data$result.folder))
      showNotification(paste0('Table of parameters is save in: ',DDD.data$result.folder))
    })

    observeEvent(input$generate.report,{
      ###download report rmarkdown
      download.file(url = 'https://github.com/wyguo/AtRTD2-DE-DAS-DTU-pipeline/blob/master/fig/DE%20DAS%20and%20DTU%20analysis%20report.Rmd',
                    destfile = 'DE DAS and DTU analysis report.Rmd')

      if(!file.exists(paste0(DDD.data$result.folder,'/Parameter summary.csv')))
        write.csv(x$df,file=paste0(DDD.data$result.folder,'/Parameter summary.csv'),row.names = F)

      if(!file.exists(paste0(DDD.data$result.folder,'/data.info.csv')))
        write.csv(DDD.data$data.info,file=paste0(DDD.data$result.folder,'/data.info.csv'),row.names = T)

      if(!file.exists(paste0(DDD.data$result.folder,'/DE DAS DTU numbers.csv')))
        write.csv(DDD.numbers(),file=paste0(DDD.data$result.folder,'/DE DAS DTU numbers.csv'),
                  row.names = F)
      if(!file.exists(paste0(DDD.data$result.folder,'/DE vs DAS gene number.csv')))
        write.csv(DEvsDAS(),file=paste0(DDD.data$result.folder,'/DE vs DAS gene number.csv'),
                  row.names = F)
      if(!file.exists(paste0(DDD.data$result.folder,'/DE vs DTU transcript number.csv')))
        write.csv(DEvsDTU(),file=paste0(DDD.data$result.folder,'/DE vs DTU transcript number.csv'))

      report.folder <- paste0(DDD.data$path,'/report')
      if(!file.exists(report.folder))
        dir.create(report.folder,recursive = T)
      withProgress(message = 'Generating report',
                   detail = 'This may take a while...', value = 0, {
                     incProgress(0.3)
                     ##generate report
                     generate.report(input.Rmd = paste0(DDD.data$path,"/report.Rmd"),
                                     report.folder = report.folder)
                     # generate.report(input.Rmd = paste0("test.Rmd"),
                     #                 report.folder = report.folder)
                     message(paste0('Reports in html, pdf and word format are saved in : ',paste0(DDD.data$path,'/report')))
                     showNotification(paste0('Reports in html, pdf and word format are saved in : ',paste0(DDD.data$path,'/report')))
                   })
    })
  }
  shinyApp(ui, server,options=list(launch.browser=T))
}
