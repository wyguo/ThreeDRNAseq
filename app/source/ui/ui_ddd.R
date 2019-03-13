#========================>> 3D analysis <<=======================
tabItem('ddd',
        tags$head(
          tags$style(type='text/css', 
                     ".nav-tabs-custom .nav-tabs li.active {font-weight: bold} ")),
        ##---------------Step 1: Select factors of interest------------
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 60%; display: block; margin-left: auto; margin-right: auto;"
                   src="DDD.png"/>')
              )
          ),
        fluidRow(
          box(title = 'Step 1: Visualise experimental design',
              width=12,status = 'primary', solidHeader = T,
              HTML('If the conditions to compare are incorrect, please go to "Data generation" page to check the selections of experimental design in the sample table.'),
              br(),
              DT::dataTableOutput('factor_table')
              
          ),
          ##---------------Step 2: Set contrast groups------------
          box(title = 'Step 2: Set contrast groups',
              width=6,status = 'primary', solidHeader = T,
              selectorUI(1),
              tags$div(id = 'placeholder'),
              actionButton(inputId = 'insertParamBtn', label = "Add a line",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                           icon = icon('plus')),
              actionButton(inputId = 'removeParamBtn', label = "Remove a line",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                           icon = icon('minus')),
              actionButton(inputId = 'generate_contrast', label = "Generate contrast groups",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                           icon = icon('hand-point-right')),
              HTML('<p>&nbsp;</p>'),
              actionButton(inputId = 'refresh_contrast', label = "Refresh contrast groups",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                           icon = icon('refresh')),
              bsTooltip(id = "refresh_contrast", 
                        title = "If the contrast groups are in a mess, click to refresh",
                        placement = "bottom", options = list(container = "body"))
              ## contrast of contrast groups is disabled
              # hr(),
              # radioButtons(inputId = 'contrast_of_contrast',label = 'Compare difference of pair-wise contrast groups?',
              #              choices = c('No','Yes'),selected = 'No',inline = T),
              # selectorUI2(1),
              # tags$div(id = 'placeholder_cofc'),
              # bsButton('insertParamBtn_cofc','Add a line',icon("ban"),
              #          style="primary"),
              # bsButton('removeParamBtn_cofc','Remove a line',icon("ban"),
              #          style="primary"),
              # actionButton(inputId = 'generate_contrast_cofc', label = "Generate contrast groups",
              #              style="color: #fff; background-color: #428bca; border-color: #2e6da4",
              #              icon = icon('hand-point-right')),
              # br(),
              # br(),
              # bsButton('refresh_contrast_cofc','Refresh contrast groups',icon("ban"),
              #          style="primary")
          ),
          box(title = 'Visualise contrast groups',status = 'primary',
              width=6, solidHeader = F,
              DT::dataTableOutput("view_contrast_table"),
              textOutput('test.contrast.text')
          )
        ),
        ##---------------Step 3: 3D analysis------------
        fluidRow(
          box(title='Step 3: 3D analysis',
              width = 6,status = 'primary', solidHeader = T,
              HTML('<h4><strong>Choose a pipeline</strong></h4>'),
              div(style="display:inline-block;vertical-align:middle;height:30px",
                  selectInput(inputId = "DE.pipeline", label = NULL,
                              choices = c("limma","glmQL"),
                              selected = 'limma',width = 160
                  )
              ),
              div(style="display:inline-block;vertical-align:middle;height:30px",
                  actionButton('run.DE','Run',icon("send outline icon"),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4")
              ),
              br(),
              br(),
              HTML('<p align="justify">In practics, limma and edgeR-glmQL have comparable performance, and both have
                   good control of false positives.</p>'),
              hr(),
              HTML('<h4><strong>Set thresholds</strong></h4>'),
              selectInput(inputId = 'pval_adj_method',label = 'P-value adjust method',
                          choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                      "fdr", "none"),selected = 'BH'),
              selectInput(inputId = 'DAS_pval_method',label = 'Summarise to DAS gene level p-value*',
                          choices = c('F-test','Simes'),selected = 'F-test'),
              sliderInput(inputId = 'pval_cut',label = 'Adjusted p-value',
                          min = 0,max = 0.05,value = 0.01,step = 0.005),
              sliderInput(inputId = 'l2fc_cut',label = 'Absolute \\(\\log_2\\)FC',
                          min = 0,max = 5,value = 1,step = 0.5),
              sliderInput(inputId = 'deltaPS_cut',label = 'Absolute \\(\\Delta\\)PS',
                          min = 0,max = 1,value = 0.1,step = 0.05),              
              HTML('<p align="justify">PS (percent of splice) is defined as the ratio of transcript average TPMs of conditions divided by the average
                   gene abundance. \\(\\Delta\\)PS is the PS differences of conditions based on the contrast groups.</p>')
              ),
          box(title='Conditions of interest for contrast groups',
              width = 6,status = 'primary', solidHeader = F,
              HTML('<img style="width: 95%; display: block; margin-left: auto; margin-right: auto;" 
                   src="contrast_design.png"/>')
              ),
          box(title='DE/DAS/DTU testing statistics',
              width = 6,status = 'primary', solidHeader = F,
              HTML('<img style="width: 100%; display: block; margin-left: auto; margin-right: auto;" 
                   src="DDD_test.png"/>')
              )
          ),
        ##---------------Step 4: 3D testing statistics------------
        fluidRow(
          box(title='Step 4: 3D testing statistics',
              width = 3,status = 'primary', solidHeader = T,
              HTML('<h4><strong>Visulise statistics</strong></h4>'),
              selectInput(inputId = 'stat.type',label = 'Targets',
                          choices = c('DE genes','DAS genes','DE transcripts','DTU transcripts'),
                          selected = 'DE genes'
              ),
              uiOutput('show.contrast.groups'),
              sliderInput(inputId = 'top.stat',label = 'Top ranked statistics in each contrast',
                          min = 10,max = 100,value = 10,step = 10),
              hr(),
              # actionButton(inputId = 'save.stat',label = 'Save all tables',icon = icon('download',lib = 'font-awesome'),
              #              style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
              # br(),
              # br(),
              HTML('<p align="justify">Up to top 100 ranked targets can be visulise in this App. The statistics of
                   full gene/transcript lists can be saved to local folder in the "Generate report" step.</p>')
              ),
          box(title='Top 100 3D targets',
              width = 9,status = 'primary', solidHeader = F,
              DT::dataTableOutput('top.stat.table')
          )
          
              ),
        fluidRow(
          ##---------------Step 5: Significant 3D numbers------------
          # uiOutput('show.top.stat.table'),
          box(title='Step 5: Significant 3D numbers',
              width = 12,status = 'primary', solidHeader = T,
              h4('Number of 3D genes/transcripts in each contrast group'),
              tableOutput('DDD.numbers')
          ),
          column(width = 3,
                 box(title='Venn diagrams of 3D targets',
                     width = 13,status = 'primary', solidHeader = F,
                     selectInput(
                       inputId = 'across.contrast.group',
                       label = 'Contrast groups (multiple selection allowed)',
                       selected = NULL,multiple = T,
                       choices = NULL
                     ),
                     spectrumInput(
                       inputId = "color4_3d_number_from",
                       label = "Gradient fill from:",
                       choices = distinct.color(50)[1:50],
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "180px"
                     ),
                     spectrumInput(
                       inputId = "color4_3d_number_to",
                       label = "Gradient fill to:",
                       choices = distinct.color(50)[2:50],
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "180px"
                     )
                 )
          ),
          column(width = 9,
                 tabBox(title = '3D',width = 13,
                        tabPanel(title = 'DE genes',
                                 plotOutput("DE.genes.euler")
                        ),
                        tabPanel(title = 'DAS genes',
                                 plotOutput("DAS.genes.euler")
                        ),
                        tabPanel(title = 'DE transcripts',
                                 plotOutput("DE.trans.euler")
                        ),
                        tabPanel(title = 'DTU transcripts',
                                 plotOutput("DTU.trans.euler")
                        )
                 ),
                 div(style="float:left;margin-left: 0px;",
                     numericInput(inputId = "across.contrast.euler.plot.scale",label = 'Plot scales (0-1)',
                                  value = 1,min = 0,max = 1,step = 0.1,width = "100%")
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "across.contrast.euler.plot.width",label = 'Plot width (inch)',
                                  value = 5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "across.contrast.euler.plot.height",label = 'Plot height (inch)',
                                  value = 4.7,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "across.contrast.euler.plot.res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 actionButton(inputId = 'save.across.contrast.euler.plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float:left;margin-top:25px;margin-left:25px")
          )
        ),
        HTML('<p>&nbsp;</p>'),
        fluidRow(
          column(width = 12,
                 tabBox(width = 13,
                        # Title can include an icon
                        # title = tagList(shiny::icon("gear"), "Loaded tables"),
                        title = "Number of up-down regulated targets",side = 'right',
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
                 )
          ),
          column(width = 12,
                 div(style="float:left;margin-left: 0px;",
                     spectrumInput(
                       inputId = "color4up_regulate",
                       label = "Color for up regulation:",
                       choices = distinct.color(50)[1:50],
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "100%"
                     )
                 ),
                 div(style="float:left;margin-left: 25px;",
                     spectrumInput(
                       inputId = "color4down_regulate",
                       label = "Color for down regulation:",
                       choices = distinct.color(50)[2:50],
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "100%"
                     )
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "updown.x.hjust",label = 'Horizontal-just x-labs',value = '0.5',
                                  min = 0,max = 1,step = 0.1,width = "100%")
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "updown.x.vjust",label = 'Vertical-just x-lab',value = '0.5',
                                  min = 0,max = 1,step = 0.1,width = "100%")
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "updown.x.rotate",label = 'Rotate x-lab',value = '0',
                                  min = 0,max = 360,step = 15,width = "100%")
                 )
          ),
          column(width = 12,
                 div(style="float:left;margin-left: 0px;",
                     numericInput(inputId = "updown.plot.width",label = 'Plot width (inch)',
                                  value = 8,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "updown.plot.height",label = 'Plot height (inch)',
                                  value = 5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "updown.plot.res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 actionButton(inputId = 'save.up.down.bar.plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; margin-top: 25px; margin-left: 25px")
          ),
          HTML("<p>&nbsp;</p>"),
          column(width = 12,
                 tabBox(width = 13,title = 'Number of union set of all contrast groups',side = 'right',
                        tabPanel(title = 'DE and DAS genes',
                                 plotOutput("genes.flow.chart")
                        ),
                        tabPanel(title = 'DE and DTU transcripts',
                                 plotOutput("trans.flow.chart")
                        )
                 ),
                 div(style="float:left;margin-left: 0px;",
                     numericInput(inputId = "flow.plot.width",label = 'Plot width (inch)',
                                  value = 8.7,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "flow.plot.height",label = 'Plot height (inch)',
                                  value = 5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "flow.plot.res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 actionButton(inputId = 'save.union.flow.chart',label = 'Save',
                              icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; margin-top: 25px; margin-left: 25px")
          )
        ),
        HTML('<p>&nbsp;</p>'),
        ##---------------Step 6: Transcriptional vs alternative splicing regulation------------
        fluidRow(
          box(title='Step 6: Transcriptional vs alternative splicing regulation',
              width = 12,status = 'primary', solidHeader = T,
              column(width = 6,
                     h4('DE vs DAS genes'),
                     wellPanel(style = "background-color: #ffffff;",
                               tableOutput('DE.vs.DAS')
                     )
              ),
              column(width = 6,
                     h4('DE vs DTU transcripts'),
                     wellPanel(style = "background-color: #ffffff;",
                               tableOutput('DE.vs.DTU') 
                     )
              )
              
          ),
          column(width = 3,
                 box(title='Transcriptional vs AS regulation',
                     width = 13,status = 'primary', solidHeader = F,
                     selectInput(inputId = 'across.target',
                                 label = 'Contrast groups',
                                 selected = NULL,
                                 multiple = F,
                                 choices = NULL),
                     spectrumInput(
                       inputId = "color4_dedas_number_from",
                       label = "First circle colour:",
                       choices = distinct.color(50)[1:50],
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "180px"
                     ),
                     spectrumInput(
                       inputId = "color4_dedas_number_to",
                       label = "Second circle colour:",
                       choices = distinct.color(50)[2:50],
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "180px"
                     )
                 )
          ),
          column(width = 9,
                 tabBox(width = 13,
                        title = "",side = 'right',
                        tabPanel("DE vs DAS genes",
                                 plotOutput("DEvsDAS.euler")
                        ),
                        tabPanel("DE vs DTU transcripts",
                                 plotOutput("DEvsDTU.euler")
                        )
                 ),
                 div(style="float:left;margin-left: 0px;",
                     numericInput(inputId = "across.target.euler.plot.scale",label = 'Plot scales (0-1)',
                                  value = 1,min = 0,max = 1,step = 0.1,width = "100%")
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "across.target.euler.plot.width",label = 'Plot width (inch)',
                                  value = 5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "across.target.euler.plot.height",label = 'Plot height (inch)',
                                  value = 4.7,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "across.target.euler.plot.res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 actionButton(inputId = 'save.across.target.euler.plot',label = 'Save',
                              icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float:left; margin-top:25px;margin-left: 25px")
                 
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_ddd',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Data pre-processing</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Functional plot</i>'),
              actionButton(inputId = 'page_after_ddd',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )


