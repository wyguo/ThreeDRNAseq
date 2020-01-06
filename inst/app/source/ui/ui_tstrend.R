#========================>> 3D analysis <<=======================
tabItem('tstrend',
        tags$head(
          tags$style(type='text/css', 
                     ".nav-tabs-custom .nav-tabs li.active {font-weight: bold} ")),
        ##--------- >> Step 1: Time-series experimental design----
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 80%; display: block; margin-left: auto; margin-right: auto;"
                   src="tspipeline.png"/>')
              # actionButton(inputId = 'load_example_data',label = 'load')
              )
        ),
        fluidRow(
          box(title = 'Description',
              width=12,status = 'primary', solidHeader = T,
              HTML("Time-series trend analysis is aim to study an experiment with many time-points in each group. 
                  When there are many time-points, it is recommended to analyse on smoothed expression changes over 
                  time instead of discrete comparisons of individual time-points in each group. This page is for DE, DAS and DTU 
                  trend analysis of group-wise (the same set experiments are performed in different groups/blocks):
                   <ul>
                   <li>Time-series</li>
                   <li>Developmental series</li>
                   <li>Consecutive conditions with a sequential order</li>
                   </ul>
                   "),
              hr(),
              radioButtons(inputId = 'whether_tstrend',label = 'Proceed trend analysis?',
                           choices = c('No','Yes'),inline = T),
              HTML('<strong>Note:</strong> If not applicable, please skip this page to the next analysis; 
                   Otherwise, please click "Yes" to expand analysis panels.'),
              HTML('<img style="width: 55%; display: block; margin-left: auto; margin-right: auto;"
                   src="tstrend.png"/>')
              )
              ),
        conditionalPanel(condition = "input.whether_tstrend == 'Yes'",
                         fluidRow(
                           box(title = 'Step 1: Time-series experimental design',
                               width=4,status = 'primary', solidHeader = T,
                               selectInput(inputId = 'tstrend_time',label = 'Select time-points/consecutive conditions column',
                                           choices = NULL,selected = NULL,width = NULL),
                               selectInput(inputId = 'tstrend_time_group',label = 'Select group column of time-points',
                                           choices = NULL,selected = NULL,width = NULL),
                               radioButtons(inputId = 'tstrend_spline',label = 'Using spline to smooth time-series',
                                            choices = c('Yes','No'),selected = 'Yes',inline = T),
                               conditionalPanel(
                                 condition = "input.tstrend_spline == 'Yes'",
                                 sliderInput(inputId = 'tstrend_spline_df',label = 'Spline degree of freedom',
                                             min = 1,max = 3,value = 1)
                               ),
                               HTML('<div align="justify"><strong>Note</strong>: (1) Development series or consecutive conditions with sequential order can be treated as
                                    "time-points" for trend analysis; (2) Spline method requires 3 or more time-points in each group; (3) If "No" is selected, trend will 
                                    be compared in the way of discrete jumping from one to another time-point in each group.</div>')
                               ),
                           box(title = NULL,
                               width=8,status = 'primary', solidHeader = F,
                               DT::dataTableOutput('tstrend_design_table')
                               
                           ),
                           box(title = NULL,
                               width=12,status = 'primary', solidHeader = F,
                               plotOutput('tstrend_example_plot')
                               
                           )
                           ),
                         ##--------- >> Step 2: Compare trend across time groupse----
                         fluidRow(
                           box(title = 'Step 2: Set groups to compare',
                               width=4,status = 'primary', solidHeader = T,
                               selectorUI_tstrend(1),
                               tags$div(id = 'placeholder_tstrend'),
                               actionButton(inputId = 'insertParamBtn_tstrend', label = "Add a line",
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                                            icon = icon('plus')),
                               actionButton(inputId = 'removeParamBtn_tstrend', label = "Remove a line",
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                                            icon = icon('minus')),
                               # br(style="clear:both"),
                               p(),
                               actionButton(inputId = 'generate_tstrend_groups', label = "Generate time groups",
                                            style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                                            icon = icon('hand-point-right'))
                           ),
                           box(title = NULL,
                               width=8,status = 'primary', solidHeader = F,
                               DT::dataTableOutput('tstrend_group_table')
                           )
                         ),
                         ##--------- >> Step 3: Time-series 3D trend analysis----
                         fluidRow(
                           column(4,
                                  box(title = 'Step 3: Time-series 3D trend analysis',
                                      width=13,status = 'primary', solidHeader = T,
                                      HTML('<strong>Choose a pipeline</strong>'),
                                      p(),
                                      div(style="display:inline-block;vertical-align:middle;height:30px",
                                          selectInput(inputId = "DE_tstrend_pipeline", label = NULL,
                                                      choices = 'limma',
                                                      selected = 'limma',width = 140
                                          )
                                      ),
                                      div(style="display:inline-block;vertical-align:middle;height:30px",
                                          selectInput(inputId = "mean_variance_tstrend_method", label = NULL,
                                                      choices = c('voom','voomWeights'),
                                                      selected = 'voom',width = 140
                                          )
                                      ),
                                      p(),
                                      selectInput(inputId = 'pval_adj_method_tstrend',label = 'P-value adjust method',
                                                  choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                              "fdr", "none"),selected = 'BH'),
                                      selectInput(inputId = 'DAS_pval_method_tstrend',label = 'Summarise DTU transcript test to DAS gene p-value',
                                                  choices = c('F-test','Simes'),selected = 'F-test'),
                                      selectInput(inputId = 'DAS_pval_method_across_group',label = 'Summarise DAS/DTU p-values across time-groups',
                                                  choices = c('Simes','Fisher'),selected = 'Simes'),
                                      HTML('<div align="justify"><strong>Note:</strong> "Simes" method is more stringent than "Fisher" method. Each contrast group involves several time groups for trend comparisons,
                                           e.g. "Null hypothesis: Day1=Day2=Day3". These methods are used to generate overall p-values across these time groups.</div>'),
                                      p(),
                                      sliderInput(inputId = 'pval_cut_tstrend',label = 'Adjusted p-value',
                                                  min = 0,max = 0.05,value = 0.01,step = 0.005),
                                      div(style="display:inline-block;vertical-align:middle;float:right; padding:4px;",
                                          actionButton('run_3d_tstrend','Run',icon("send outline icon"),
                                                       style="color: #fff; background-color: #428bca; border-color: #2e6da4; font-size:150%; width:100%; height:150%; ")
                                      )
                                  )
                           ),
                           column(8,
                                  box(title = 'Testing statistics',
                                      width=13,status = 'primary', solidHeader = F,
                                      div(style="display:inline-block;vertical-align:top;margin-right:15px; float: left;",
                                          selectInput(inputId = 'tstrend_test_stat_type',label = 'Targets',
                                                      choices = c('DE trend genes','DAS trend genes','DE trend transcripts','DTU trend transcripts'),
                                                      selected = 'DE trend genes',width = '250px')
                                      ),
                                      div(style="display:inline-block;vertical-align:top;margin-left:15px; float: left;",
                                          sliderInput(inputId = 'top_n_tstrend_stat',label = 'Top n statistics of significance',
                                                      min = 10,max = 100,value = 10,step = 10,width = '250px')
                                      ),
                                      DT::dataTableOutput('tstrend_test_stat_table')
                                  ),
                                  box(title = 'TS trend 3D number',
                                      width=13,status = 'primary', solidHeader = F,
                                      tableOutput('tstrend_test_number_table')
                                  )
                           )
                         ),
                         fluidRow(
                           ##--------------->> Step 4: Profile plot----------------
                           box(title='Step 4: Profile plot',
                               width = 12,status = 'primary', solidHeader = T,
                               column(width = 6,
                                      wellPanel(
                                        HTML('<h4><strong>Visualise plot of a single gene</strong></h4>'),
                                        textInput(inputId = 'tstrend_gene_id',label = 'Input a gene',value = ''),
                                        HTML('E.g. AT1G01060 in Arabidopsis (gene name must match to provided transcript-gene mapping).'),
                                        br(),
                                        radioButtons(inputId = 'tstrend_profile_data_type',label = 'Select data type',
                                                     choices = c('TPM','Read counts'),inline = T),
                                        radioButtons(inputId = 'tstrend_profile_filter_lowexpressed',label = 'Filter low expressed transcripts?',
                                                     choices = c('Yes','No'),inline = T),
                                        selectInput(inputId = 'tstrend_profile_slice_group',label = 'Slice profile plot on group',
                                                    choices = NULL,
                                                    selected = NULL),
                                        
                                        div(style="float:left; ",
                                            numericInput(inputId = "tstrend.profile.x.hjust",label = 'Horizontal-just x-labs',value = '0.5',
                                                         min = 0,max = 1,step = 0.1,width = "100%")
                                        ),
                                        div(style="float:left;margin-left: 25px;",
                                            numericInput(inputId = "tstrend.profile.x.vjust",label = 'Vertical-just x-lab',value = '0.5',
                                                         min = 0,max = 1,step = 0.1,width = "100%")
                                        ),
                                        div(style="float:left;margin-left: 25px;",
                                            numericInput(inputId = "tstrend.profile.x.rotate",label = 'Rotate x-lab',value = '0',
                                                         min = 0,max = 360,step = 15,width = "100%")
                                        ),
                                        br(style="clear:both"),
                                        div(style="float:left;",
                                            actionButton(inputId = 'tstrend_make_profile_plot',label = 'Plot',
                                                         icon = icon('pencil',lib = 'font-awesome'),
                                                         style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                                        ),
                                        br()
                                      )
                               ),
                               column(width = 6,
                                      wellPanel(
                                        HTML('<h4><strong>Make plots of multiple genes</strong></h4>'),
                                        fileInput("tstrend_multiple_gene_input", "Choose gene list csv file",
                                                  accept = c(
                                                    "text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv")
                                        ),
                                        HTML('<strong>Note:</strong> First line of the csv file is treated as header.'),
                                        p(),
                                        radioButtons(inputId = 'tstrend_multiple_plot_type',label = 'Select plot type',
                                                     choices = c('Abundance','PS','Both'),inline = T),
                                        radioButtons(inputId = 'tstrend_multi_plot_format',label = 'Select format',
                                                     choices = c('png','pdf','both'),inline = T),
                                        div(style="display: inline-block; margin-right: 25px;",
                                            numericInput(inputId = "tstrend_ps_multiplot_res",label = 'PNG plot resolution',value = '150',
                                                         min = 0,max = 600,step = 60,width = "100%")
                                        ),
                                        div(style="display: inline-block; margin-right: 25px;",
                                            numericInput(inputId = "tstrend_ps_multiplot_height",label = 'Plot height (inch)',
                                                         value = 4.5,width = "100%",step = 0.2)
                                        ),
                                        div(style="display: inline-block; margin-right: 25px;",
                                            numericInput(inputId = "tstrend_ps_multiplot_width",label = 'Plot width (inch)',
                                                         value = 7,width = "100%",step = 0.2)
                                            
                                        ),
                                        # textInput(inputId = 'multiplot_folder_namae',
                                        #           label = 'Figures are saved to',value = 'figure',width = "100%",
                                        #           placeholder = 'New folder with this name will be created in the figure folder'),
                                        actionButton(inputId = 'tstrend_make_multiple_plot',label = 'Run',
                                                     icon = icon('send outline icon',lib = 'font-awesome'),
                                                     style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                                        br()
                                      )
                               ),
                               column(12,
                                      verbatimTextOutput('ddd_tstrend_profile_plot_folder_text'),
                                      HTML('<strong>Note:</strong> If the App is not running locally, these plots can be downloaded in the final step "Generate report".')
                               )
                           ),
                           # column(width = 12,
                           tabBox(title = 'Profile plot',width = 12,side = 'left',
                                  tabPanel(title = 'Abundance',
                                           plotlyOutput('tstrend_profile_plot_panel')
                                  ),
                                  tabPanel(title = 'Percent spliced',
                                           plotlyOutput('tstrend_ps_plot_panel')
                                  )
                           ),
                           column(width = 12,
                                  div(style="float:left;margin-left: 25px;",
                                      numericInput(inputId = "tstrend_ps_plot_width",label = 'Plot width (inch)',
                                                   value = 7,width = "100%",step = 0.2)
                                  ),
                                  div(style="float:left;margin-left: 25px;",
                                      numericInput(inputId = "tstrend_ps_plot_height",label = 'Plot height (inch)',
                                                   value = 4.5,width = "100%",step = 0.2)
                                  ),
                                  div(style="float:left;margin-left: 25px;",
                                      numericInput(inputId = "tstrend_ps_plot_res",label = 'PNG plot resolution',value = '150',
                                                   min = 0,max = 600,step = 60,width = "100%")
                                  ),
                                  actionButton(inputId = 'tstrend_save_abundance_plot',label = 'Save',
                                               icon = icon('download',lib = 'font-awesome'),
                                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; 
                       margin-top: 25px;margin-left: 25px;"),
                                  bsTooltip(id = "tstrend_save_abundance_plot", 
                                            title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                                            placement = "bottom", options = list(container = "body")),
                                  actionButton(inputId = 'tstrend_save_abundance_plot_view',label = 'View saved plot',
                                               icon = icon('eye',lib = 'font-awesome'),
                                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                       float: left; margin-top: 25px; margin-left: 25px"),
                                  bsTooltip(id = "tstrend_save_abundance_plot_view",
                                            title = "Preview saved plots to correct width and height.",
                                            placement = "bottom", options = list(container = "body"))
                           )
                         )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_tstrend',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>3D analysis</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Functional plot</i>'),
              actionButton(inputId = 'page_after_tstrend',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )

