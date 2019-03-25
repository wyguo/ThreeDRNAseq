
#========================>> Functional interpretation <<=======================
tabItem('function',
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 70%; display: block; margin-left: auto; margin-right: auto;" 
                   src="function.png"/>')
          )
        ),
        ##----------Heatmap------------
        fluidRow(
          # column(width = 12,
          box(title='Heatmap',
              width = 3,status = 'primary', solidHeader = T,
              HTML('<h4><strong>Select targets</strong></h4>'),
              # radioButtons(inputId = 'heatmap.select.or.upload',label = '',
              #              choices = c('Select targets','Upload target list'),inline = T),
              selectInput(inputId = 'heatmap.target.type',label = 'Select a 3D list',
                          choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts')
              ),
              hr(),
              # fileInput("heatmap.target.list.input", "OR Choose target list csv file",
              #           accept = c(
              #             "text/csv",
              #             "text/comma-separated-values,text/plain",
              #             ".csv")
              # ),
              # radioButtons(inputId = 'heatmap.targetlist.type',label = 'Target list type',inline = T,
              #              choices = c('Gene list','Transcript list')
              # ),
              # hr(),
              HTML('<h4><strong>Clustering</strong></h4>'),
              # radioButtons(inputId = 'TPM_or_Zscore',label = 'Expression value',choices = c('Z-scores of TPM','TPM'),
              #              selected = 'Z-scores of TPM',inline = T),
              selectInput(inputId = 'dist.method',label = 'Distance measurement',
                          choices = c("euclidean","maximum","manhattan","canberra","binary","minkowski")),
              selectInput(inputId = 'cluster.method',label = 'Clustering method',
                          choices = c("ward.D","ward.D2","average","complete","single","mcquitty","median","centroid"),
                          selected = 'ward.D'),
              numericInput(inputId = 'cluster.number',label = 'Cluster number',value = 10,min = 1,max = 50),
              HTML('See the R functions for hierarchical clustering: <a href="https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust" target="_blank">hclust</a> and
                   <a href="https://www.rdocumentation.org/packages/proxy/versions/0.4-22/topics/dist" target="_blank">dist</a>'),
              hr(),
              HTML('<h4><strong>Options</strong></h4>'),
              radioButtons(inputId = 'show_row_labels',label = 'Show row labels',choices = c('No','Yes'),
                           selected = 'No',inline = T),
              radioButtons(inputId = 'show_column_labels',label = 'Show column labels',choices = c('Yes','No'),
                           selected = 'Yes',inline = T),
              HTML('<h5><strong>Select colours of values</strong></h5>'),
              div(style="display: inline-block; margin-right: 5px;",
                  spectrumInput(
                    inputId = "color4_heatmap_neg",
                    label = "  min",
                    choices = unique(c('blue',distinct.color(50)[1:50])),
                    options = list(`toggle-palette-more-text` = "Show more"),
                    width = "70px"
                  )
              ),
              div(style="display: inline-block; margin-right: 5px;",
                  spectrumInput(
                    inputId = "color4_heatmap_0",
                    label = HTML(" median"),
                    choices = unique(c('white',distinct.color(50)[1:50])),
                    options = list(`toggle-palette-more-text` = "Show more"),
                    width = "70px"
                  )
              ),
              div(style="display: inline-block; margin-right: 5px;",
                  spectrumInput(
                    inputId = "color4_heatmap_pos",
                    label = "  max",
                    choices = unique(c('red',distinct.color(50)[1:50])),
                    options = list(`toggle-palette-more-text` = "Show more"),
                    width = "70px"
                  )
              ),
              actionButton(inputId = 'plot.heatmap',label = 'Plot heatmap',
                           icon = icon('pencil',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4")
              ),
          box(title=NULL,
              width = 9,status = 'primary', solidHeader = T,
              plotOutput('heatmap.plot.panel',height = 730),
              actionButton(inputId = 'save.target.in.cluster',label = 'Save target list in clusters',
                           icon = icon('download',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px; margin-left: 25px"),
              actionButton(inputId = 'save.heatmap',label = 'Save heatmap',
                           icon = icon('download',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
              div(style="float:right;margin-right: 25px;",
                  numericInput(inputId = "heatmap.plot.res",label = 'PNG plot resolution',value = '150',
                               min = 0,max = 600,step = 60,width = "100%")
              ),
              div(style="float:right;margin-right: 25px;",
                  numericInput(inputId = "heatmap.plot.height",label = 'Plot height (inch)',
                               value = 8,width = "100%",step = 0.2)
              ),
              div(style="float:right;margin-right: 25px;",
                  numericInput(inputId = "heatmap.plot.width",label = 'Plot width (inch)',
                               value = 4.5,width = "100%",step = 0.2)
              )
          )
        ),
        fluidRow(
          ##---------------multiple plots----------------
          # column(width = 12,
          box(title='Profiles',
              width = 12,status = 'primary', solidHeader = T,
              column(width = 6,
                     wellPanel(
                       HTML('<h4><strong>Visualise plot of a single gene</strong></h4>'),
                       textInput(inputId = 'gene.id',label = 'Input a gene',value = ''),
                       HTML('E.g. AT1G01060 in Arabidopsis (gene name must match to provided transcript-gene mapping).'),
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
                       HTML('<h4><strong>Make plots of multiple genes</strong></h4>'),
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
                       div(style="display: inline-block; margin-right: 25px;",
                           numericInput(inputId = "ps.multiplot.res",label = 'PNG plot resolution',value = '150',
                                        min = 0,max = 600,step = 60,width = "100%")
                       ),
                       div(style="display: inline-block; margin-right: 25px;",
                           numericInput(inputId = "ps.multiplot.height",label = 'Plot height (inch)',
                                        value = 4.5,width = "100%",step = 0.2)
                       ),
                       div(style="display: inline-block; margin-right: 25px;",
                           numericInput(inputId = "ps.multiplot.width",label = 'Plot width (inch)',
                                        value = 7,width = "100%",step = 0.2)
                           
                       ),
                       # textInput(inputId = 'multiplot_folder_namae',
                       #           label = 'Figures are saved to',value = 'figure',width = "100%",
                       #           placeholder = 'New folder with this name will be created in the figure folder'),
                       actionButton(inputId = 'make.multiple.plot',label = 'Run',
                                    icon = icon('send outline icon',lib = 'font-awesome'),
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                       br()
                     )
              )
          ),
          # column(width = 12,
          tabBox(title = 'Profile plot',width = 12,side = 'left',
                 tabPanel(title = 'Abundance',
                          plotlyOutput('profile.plot.panel')
                 ),
                 tabPanel(title = 'Percent spliced',
                          plotlyOutput('ps.plot.panel')
                 )
          ),
          actionButton(inputId = 'save.ps.plot',label = 'Save PS plot',
                       icon = icon('download',lib = 'font-awesome'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px; margin-left: 25px"),
          actionButton(inputId = 'save.abundance.plot',label = 'Save abundance plot',
                       icon = icon('download',lib = 'font-awesome'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px"),
          div(style="float:right;margin-right: 25px;",
              numericInput(inputId = "ps.plot.res",label = 'PNG plot resolution',value = '150',
                           min = 0,max = 600,step = 60,width = "100%")
          ),
          div(style="float:right;margin-right: 25px;",
              numericInput(inputId = "ps.plot.height",label = 'Plot height (inch)',
                           value = 4.5,width = "100%",step = 0.2)
          ),
          div(style="float:right;margin-right: 25px;",
              numericInput(inputId = "ps.plot.width",label = 'Plot width (inch)',
                           value = 7,width = "100%",step = 0.2)
          )
        ),
        HTML('&nbsp;'),
        fluidRow(
          ##---------------GO plots----------------
          column(width = 3,
                 box(title = 'GO annotation at gene level',
                     width = 13,status = 'primary', solidHeader = T,
                     p() %>%
                       helper(icon = "question-circle",
                              colour = NULL,
                              type = 'markdown',
                              title = 'GO annotation plot',
                              size = 'l',
                              content = 'go',
                              style="font-size: 2.0rem;margin-right:50px;"),
                     HTML('&nbsp'),
                     downloadButton('download_target_list', 'Download gene list',class="btn btn-primary",
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4;float:left"),
                     br(),
                     br(),
                     br(),
                     HTML('<div align="justify">Gene lists of DE genes, DAS genes, DE transcripts and DTU transcripts will be downloaded and saved in scv file. Users can
                          use these lists and GO analysis tools/websites to create GO annotation table for plot.</div>'),
                     hr(),
                     fileInput("go.annotation.input", "Choose GO annotation csv file",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     ),
                     # uiOutput('show.go.table.column'),
                     selectInput(inputId = 'select.go.table.column',label = 'Select column of statistics',
                                 choices = NULL),
                     selectInput(inputId = 'go_plot_label',label = 'Plot title',
                                 choices = c('DE genes','DAS genes','Genes of DE transcripts','Genes of DTU transcripts'),
                                 selected = 'DE genes'),
                     bsTooltip(id = "go_plot_label", 
                               title = "Plot title to distinguish different plots",
                               placement = "bottom", options = list(container = "body")),
                     numericInput(inputId = 'go.text.cut',label = 'GO term text length',value = 50),
                     actionButton(inputId = 'make.go.plot',label = 'Plot',
                                  icon = icon('pencil',lib = 'font-awesome'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left")
                     
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
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;float: right; margin-top: 25px"),
              div(style="float:right;margin-right: 25px;",
                  numericInput(inputId = "go.plot.res",label = 'PNG plot resolution',value = '150',
                               min = 0,max = 600,step = 60,width = "100%")
              ),
              div(style="float:right;margin-right: 25px;",
                  numericInput(inputId = "go.plot.height",label = 'Plot height (inch)',
                               value = 4.5,width = "100%",step = 0.2)
              ),
              div(style="float:right;margin-right: 25px;",
                  numericInput(inputId = "go.plot.width",label = 'Plot width (inch)',
                               value = 7.8,width = "100%",step = 0.2)
              )
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_function',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>3D analysis</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>TSIS</i>'),
              actionButton(inputId = 'page_after_function',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )