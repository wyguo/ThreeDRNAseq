
#========================>> Functional interpretation <<=======================
tabItem('function',
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 50%; display: block; margin-left: auto; margin-right: auto;" 
                   src="function.png"/>')
              )
          ),
        ##---------->> Heatmap------------
        fluidRow(
          # column(width = 12,
          box(title='Step 1: Heatmap',
              width = 3,status = 'primary', solidHeader = T,
              HTML('<h4><strong>Select targets</strong></h4>'),
              # radioButtons(inputId = 'heatmap.select.or.upload',label = '',
              #              choices = c('Select targets','Upload target list'),inline = T),
              selectInput(inputId = 'heatmap.target.type',label = 'Select a 3D list',
                          choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts')
              ),
              hr(),
              # hr(),
              HTML('<h4><strong>Clustering</strong></h4>'),
              # radioButtons(inputId = 'TPM_or_Zscore',label = 'Expression value',choices = c('Z-scores of TPM','TPM'),
              #              selected = 'Z-scores of TPM',inline = T),
              selectInput(inputId = 'dist_method',label = 'Distance measurement',
                          choices = c("euclidean","maximum","manhattan","canberra","binary","minkowski")),
              selectInput(inputId = 'cluster_method',label = 'Clustering method',
                          choices = c("ward.D","ward.D2","average","complete","single","mcquitty","median","centroid"),
                          selected = 'ward.D'),
              numericInput(inputId = 'cluster_number',label = 'Cluster number',value = 10,min = 1,max = 50),
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
              div(style="float:left;margin-left: 25px;",
                  numericInput(inputId = "heatmap.plot.res",label = 'PNG plot resolution',value = '150',
                               min = 0,max = 600,step = 60,width = "100%")
              ),
              div(style="float:left;margin-left: 25px;",
                  numericInput(inputId = "heatmap.plot.height",label = 'Plot height (inch)',
                               value = 8,width = "100%",step = 0.2)
              ),
              div(style="float:left;margin-left: 25px;",
                  numericInput(inputId = "heatmap.plot.width",label = 'Plot width (inch)',
                               value = 4.5,width = "100%",step = 0.2)
              ),
              # downloadButton(outputId = 'save_heatmap',label = 'Save',
              #                class="btn btn-primary",
              #                style="color: #fff; background-color: #428bca; border-color: #2e6da4;float:left;margin-top: 25px; 
              #                margin-left: 25px")
              actionButton(inputId = 'save_heatmap',label = 'Save heatmap',
                           icon = icon('download',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                           float: left; margin-top: 25px;margin-left: 25px;"),
              bsTooltip(id = "save_heatmap", 
                        title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                        placement = "bottom", options = list(container = "body")),
              actionButton(inputId = 'save_heatmap_view',label = 'View saved plot',
                           icon = icon('eye',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                           float: left; margin-top: 25px; margin-left: 25px"),
              bsTooltip(id = "save_heatmap_view",
                        title = "Preview saved plots to correct width and height.",
                        placement = "bottom", options = list(container = "body"))
              )
        ),
        ##---------->> Time-series trend Heatmap------------
        conditionalPanel(condition = "input.whether_tstrend == 'Yes'",
                         fluidRow(
                           # column(width = 12,
                           box(title='Time-series trend heatmap',
                               width = 3,status = 'primary', solidHeader = T,
                               HTML('<h4><strong>Select targets</strong></h4>'),
                               # radioButtons(inputId = 'heatmap.select.or.upload',label = '',
                               #              choices = c('Select targets','Upload target list'),inline = T),
                               selectInput(inputId = 'tstrend.heatmap.target.type',label = 'Select a 3D list',
                                           choices = c('DE tstrend genes','DAS tstrend genes','DE tstrend transcripts','DTU tstrend transcripts')
                               ),
                               selectInput(inputId = 'tstrend_heatmap_contrast_groups',label = 'Select a contrast group',
                                           choices = 'All',multiple = F,selected = 'All'),
                               HTML('<strong>Note: </strong> If select "All", union set of all contrast groups is used.'),
                               hr(),
                               HTML('<h4><strong>Clustering</strong></h4>'),
                               # radioButtons(inputId = 'TPM_or_Zscore',label = 'Expression value',choices = c('Z-scores of TPM','TPM'),
                               #              selected = 'Z-scores of TPM',inline = T),
                               selectInput(inputId = 'tstrend_dist_method',label = 'Distance measurement',
                                           choices = c("euclidean","maximum","manhattan","canberra","binary","minkowski")),
                               selectInput(inputId = 'tstrend_cluster_method',label = 'Clustering method',
                                           choices = c("ward.D","ward.D2","average","complete","single","mcquitty","median","centroid"),
                                           selected = 'ward.D'),
                               numericInput(inputId = 'tstrend_cluster_number',label = 'Cluster number',value = 10,min = 1,max = 50),
                               HTML('See the R functions for hierarchical clustering: <a href="https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust" target="_blank">hclust</a> and
                                    <a href="https://www.rdocumentation.org/packages/proxy/versions/0.4-22/topics/dist" target="_blank">dist</a>'),
                               hr(),
                               HTML('<h4><strong>Options</strong></h4>'),
                               radioButtons(inputId = 'tstrend_show_row_labels',label = 'Show row labels',choices = c('No','Yes'),
                                            selected = 'No',inline = T),
                               radioButtons(inputId = 'tstrend_show_column_labels',label = 'Show column labels',choices = c('Yes','No'),
                                            selected = 'Yes',inline = T),
                               HTML('<h5><strong>Select colours of values</strong></h5>'),
                               div(style="display: inline-block; margin-right: 5px;",
                                   spectrumInput(
                                     inputId = "color4_tstrend_heatmap_neg",
                                     label = "  min",
                                     choices = unique(c('blue',distinct.color(50)[1:50])),
                                     options = list(`toggle-palette-more-text` = "Show more"),
                                     width = "70px"
                                   )
                               ),
                               div(style="display: inline-block; margin-right: 5px;",
                                   spectrumInput(
                                     inputId = "color4_tstrend_heatmap_0",
                                     label = HTML(" median"),
                                     choices = unique(c('white',distinct.color(50)[1:50])),
                                     options = list(`toggle-palette-more-text` = "Show more"),
                                     width = "70px"
                                   )
                               ),
                               div(style="display: inline-block; margin-right: 5px;",
                                   spectrumInput(
                                     inputId = "color4_tstrend_heatmap_pos",
                                     label = "  max",
                                     choices = unique(c('red',distinct.color(50)[1:50])),
                                     options = list(`toggle-palette-more-text` = "Show more"),
                                     width = "70px"
                                   )
                               ),
                               actionButton(inputId = 'plot.tstrend.heatmap',label = 'Plot heatmap',
                                            icon = icon('pencil',lib = 'font-awesome'),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                               ),
                           box(title=NULL,
                               width = 9,status = 'primary', solidHeader = T,
                               plotOutput('tstrend.heatmap.plot.panel',height = 730),
                               div(style="float:left;margin-left: 25px;",
                                   numericInput(inputId = "tstrend.heatmap.plot.res",label = 'PNG plot resolution',value = '150',
                                                min = 0,max = 600,step = 60,width = "100%")
                               ),
                               div(style="float:left;margin-left: 25px;",
                                   numericInput(inputId = "tstrend.heatmap.plot.height",label = 'Plot height (inch)',
                                                value = 8,width = "100%",step = 0.2)
                               ),
                               div(style="float:left;margin-left: 25px;",
                                   numericInput(inputId = "tstrend.heatmap.plot.width",label = 'Plot width (inch)',
                                                value = 4.5,width = "100%",step = 0.2)
                               ),
                               actionButton(inputId = 'save_tstrend_heatmap',label = 'Save heatmap',
                                            icon = icon('download',lib = 'font-awesome'),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                                            float: left; margin-top: 25px;margin-left: 25px;"),
                               bsTooltip(id = "save_tstrend_heatmap", 
                                         title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                                         placement = "bottom", options = list(container = "body")),
                               actionButton(inputId = 'save_tstrend_heatmap_view',label = 'View saved plot',
                                            icon = icon('eye',lib = 'font-awesome'),
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                                            float: left; margin-top: 25px; margin-left: 25px"),
                               bsTooltip(id = "save_tstrend_heatmap_view",
                                         title = "Preview saved plots to correct width and height.",
                                         placement = "bottom", options = list(container = "body"))
                               )
                         )
                         ),
        HTML('&nbsp;'),
        fluidRow(
          ##--------------->> GO plots----------------
          column(width = 3,
                 box(title = 'Step 2: GO annotation at gene level',
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
              div(style="float:left;margin-left: 0px;",
                  numericInput(inputId = "go.plot.res",label = 'PNG plot resolution',value = '150',
                               min = 0,max = 600,step = 60,width = "100%")
              ),
              div(style="float:left;margin-left: 25px;",
                  numericInput(inputId = "go.plot.height",label = 'Plot height (inch)',
                               value = 4.5,width = "100%",step = 0.2)
              ),
              div(style="float:left;margin-left: 25px;",
                  numericInput(inputId = "go.plot.width",label = 'Plot width (inch)',
                               value = 7.8,width = "100%",step = 0.2)
              ),
              actionButton(inputId = 'save_go_plot',label = 'Save',
                           icon = icon('download',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: 
                   #2e6da4;float: left; margin-top: 25px;margin-left: 25px;"),
              bsTooltip(id = "save_go_plot", 
                        title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                        placement = "bottom", options = list(container = "body")),
              actionButton(inputId = 'save_go_plot_view',label = 'View saved plot',
                           icon = icon('eye',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                   float: left; margin-top: 25px; margin-left: 25px"),
              bsTooltip(id = "save_go_plot_view",
                        title = "Preview saved plots to correct width and height.",
                        placement = "bottom", options = list(container = "body"))
              
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_function',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Time-series trend</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>TSIS</i>'),
              actionButton(inputId = 'page_after_function',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )
