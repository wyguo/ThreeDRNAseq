#========================>> 3D analysis <<=======================
tabItem('ddd',
        tags$head(
          tags$style(type='text/css', 
                     ".nav-tabs-custom .nav-tabs li.active {font-weight: bold} ")),
        ##---------------Step 1: Select factors of interest------------
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 100%; display: block; margin-left: auto; margin-right: auto;"
                   src="DDD.png"/>')
              )
          ),
        fluidRow(
          box(title = 'Step 1: Experimental design',
              width=12,status = 'primary', solidHeader = T,
              HTML('If the information is not incorrect, please go to "Data generation" page to check the selections of experimental design in the sample table.'),
              br(),
              DT::dataTableOutput('factor_table')
              
          ),
          ##---------------Step 2: Set contrast groups------------
          box(title = 'Step 2: Set contrast groups',
              width=6,status = 'primary', solidHeader = T,
              # checkboxGroupInput(inputId = 'contrast_type',label = 'Select contrast group type (multi-selection allowed)',
              #                    choices = c('Difference of pair-wise group',
              #                                'Difference of group mean','Difference of group difference'),
              #                    selected = 'Difference of pair-wise group') %>%
              #   helper(icon = "question-circle",
              #          colour = NULL,
              #          type = 'markdown',
              #          title = 'Set contrast groups',
              #          size = 'l',
              #          content = 'set_contrast',
              #          style="font-size: 2.0rem;margin-right:50px;"),
              # hr(),
              radioButtons(inputId = 'contrast_type_pairwise',label = 'Difference of pair-wise group',
                                choices = c('Yes','No'),selected = 'Yes',inline = T) %>%
                helper(icon = "question-circle",
                       colour = NULL,
                       type = 'markdown',
                       title = 'Set contrast groups',
                       size = 'l',
                       content = 'set_contrast',
                       style="font-size: 2.0rem;margin-right:50px;"),
              selectorUI(1),
              tags$div(id = 'placeholder'),
              actionButton(inputId = 'insertParamBtn', label = "Add a line",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                           icon = icon('plus')),
              actionButton(inputId = 'removeParamBtn', label = "Remove a line",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4",
                           icon = icon('minus')),
              p(),
              HTML('<strong>Note: </strong> Select a single group label in each side of "VS" to make pair-wise 
                   contrast groups, e.g. "B VS A" for "B-A", "C VS A" for "C-A", etc.'),
              hr(),
              radioButtons(inputId = 'contrast_type_mean',label = 'Difference of multiple group mean',
                           choices = c('Yes','No'),selected = 'No',inline = T) %>%
                helper(icon = "question-circle",
                       colour = NULL,
                       type = 'markdown',
                       title = 'Set contrast groups',
                       size = 'l',
                       content = 'set_contrast',
                       style="font-size: 2.0rem;margin-right:50px;"),
              selectorUI_mean(1),
              tags$div(id = 'placeholder_mean'),
              bsButton(inputId = 'insertParamBtn_mean', label = "Add a line",
                           style="primary",
                           icon = icon('plus')),
              bsButton(inputId = 'removeParamBtn_mean', label = "Remove a line",
                           style="primary",
                           icon = icon('minus')),
              p(),
              HTML('<strong>Note: </strong> Select multiple group labels in each side of "VS" to make
                   contrast groups to compare group means, e.g. "A,B VS C,D" for "(A+B)/2-(C+D)/2", 
                   "A,B VS C,D,E" for "(A+B)/2-(C+D+E)/3", etc.'),
              hr(),
              radioButtons(inputId = 'contrast_type_diff',label = 'Difference of pair-wise group difference',
                           choices = c('Yes','No'),selected = 'No',inline = T) %>%
                helper(icon = "question-circle",
                       colour = NULL,
                       type = 'markdown',
                       title = 'Set contrast groups',
                       size = 'l',
                       content = 'set_contrast',
                       style="font-size: 2.0rem;margin-right:50px;"),
              selectorUI_pgdiff(1),
              tags$div(id = 'placeholder_pgdiff'),
              bsButton(inputId = 'insertParamBtn_pgdiff', label = "Add a line",
                           style="primary",
                           icon = icon('plus')),
              bsButton(inputId = 'removeParamBtn_pgdiff', label = "Remove a line",
                           style="primary",
                           icon = icon('minus')),
              p(),
              HTML('<strong>Note: </strong> Select two group labels in each side of "VS" to make
                   contrast groups to compare pair-wise group differences, 
                   e.g. "A,B VS C,D" for "(A-B)-(C-D)", "A,B VS C,E" for "(A-B)-(C-E)", etc.'),
              HTML('<p>&nbsp;</p>'),
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
          ),
          box(title='Factors of experimental design for contrast groups',
              width = 6,status = 'primary', solidHeader = F,
              HTML('<img style="width: 95%; display: block; margin-left: auto; margin-right: auto;" 
                   src="contrast_design.png"/>')
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
                  # selectInput(inputId = "DE_pipeline", label = NULL,
                  #             choices = c("limma","glmQL"),
                  #             selected = 'limma',width = 160
                  # )
                  selectInput(inputId = "DE_pipeline", label = NULL,
                              # choices = list(`limma-voom:`=list('limma'),`edgeR:` = c("glmQL", "glm")),
                              choices = 'limma',
                              selected = 'limma',width = 160
                  )
                  ),
              div(style="display:inline-block;vertical-align:middle;height:30px",
                  selectInput(inputId = "mean_variance_method", label = NULL,
                              # choices = list(`limma-voom:`=list('limma'),`edgeR:` = c("glmQL", "glm")),
                              choices = c('voom','voomWeights'),
                              selected = 'voom',width = 160
                  )
              ),
              br(),
              br(),
              HTML('<p align="justify">In addition to batch effects that occur for all the samples in a certain biological replicate, RNA-seq data
                   may have variations in sample quality due to, for example, degradation or contamination of specific samples. These problematic
                   samples are often shown as outliers in the PCA plot. In this case, the <strong>limma-voomWeights</strong> pipeline can be used to balance the
                   outliners. Note: <strong>limma-voom</strong> is more stringent than <strong>limma-voomWeights</strong>. User can select a proper pipeline based on the data
                   details. </p>'),
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
                   gene abundance. \\(\\Delta\\)PS is the PS differences of conditions based on the contrast groups.</p>'),
              div(style="display:inline-block;vertical-align:middle;float:right; padding:4px;",
                  actionButton('run.DE','Run',icon("send outline icon"),
                               style="color: #fff; background-color: #428bca; border-color: #2e6da4; font-size:150%; width:100%; height:150%; ")
              ),
              br(style="clear:both")
              ),
          box(title='DE/DAS/DTU testing statistics',
              width = 6,status = 'primary', solidHeader = F,
              HTML('<img style="width: 100%; display: block; margin-left: auto; margin-right: auto;" 
                   src="DDD_test.png"/>'),
              hr(),
              HTML('<p align="justify"><strong>Limma pipeline</strong> is recommended for 3D analysis. We compared different expression comparision pipelines, 
                    such as Sleuth (Pimentel et al., 2017), DESeq2 (Love et al., 2014) and EdgeR (Robinson et al., 2011), Limma is more capable
                   for complex experimental design and has more robust predictions without losing stringency (Smyth et al. 2013). At differential gene expression level, 
                   limma and edgeR have comparable performance. Sleuth and DESeq2 have no function of differential alternative splicing analysis, and 
                   using edgeR pipeline will give greatly reduced predictions of DAS genes/DTU transcripts compared with Limma. </p>')
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
          ##---------------Step 5: profile plot----------------
          # column(width = 12,
          box(title='Step 5: Profile plot',
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
                       selectInput(inputId = 'profile_slice_group',label = 'Slice profile plot on group',
                                   choices = NULL,
                                   selected = NULL),
                       div(style="float:left; ",
                           numericInput(inputId = "function.profile.x.hjust",label = 'Horizontal-just x-labs',value = '0.5',
                                        min = 0,max = 1,step = 0.1,width = "100%")
                       ),
                       div(style="float:left;margin-left: 25px;",
                           numericInput(inputId = "function.profile.x.vjust",label = 'Vertical-just x-lab',value = '0.5',
                                        min = 0,max = 1,step = 0.1,width = "100%")
                       ),
                       div(style="float:left;margin-left: 25px;",
                           numericInput(inputId = "function.profile.x.rotate",label = 'Rotate x-lab',value = '0',
                                        min = 0,max = 360,step = 15,width = "100%")
                       ),
                       br(style="clear:both"),
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
                       HTML('<strong>Note:</strong> First line of the csv file is treated as header.'),
                       p(),
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
              ),
              column(12,
                     verbatimTextOutput('ddd_profile_plot_folder_text'),
                     HTML('<strong>Note:</strong> If the App is not running locally, these plots can be downloaded in the final step "Generate report".')
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
          # actionButton(inputId = 'save.ps.plot',label = 'Save PS plot',
          #              icon = icon('download',lib = 'font-awesome'),
          #              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right; margin-top: 25px; margin-left: 25px"),
          div(style="float:left;margin-left: 25px;",
              numericInput(inputId = "ps.plot.width",label = 'Plot width (inch)',
                           value = 7,width = "100%",step = 0.2)
          ),
          div(style="float:left;margin-left: 25px;",
              numericInput(inputId = "ps.plot.height",label = 'Plot height (inch)',
                           value = 4.5,width = "100%",step = 0.2)
          ),
          div(style="float:left;margin-left: 25px;",
              numericInput(inputId = "ps.plot.res",label = 'PNG plot resolution',value = '150',
                           min = 0,max = 600,step = 60,width = "100%")
          ),
          actionButton(inputId = 'save_abundance_plot',label = 'Save',
                       icon = icon('download',lib = 'font-awesome'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left;
                       margin-top: 25px;margin-left: 25px;"),
          bsTooltip(id = "save_abundance_plot", 
                    title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                    placement = "bottom", options = list(container = "body")),
          actionButton(inputId = 'save_abundance_plot_view',label = 'View saved plot',
                       icon = icon('eye',lib = 'font-awesome'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                       float: left; margin-top: 25px; margin-left: 25px"),
          bsTooltip(id = "save_abundance_plot_view",
                    title = "Preview saved plots to correct width and height.",
                    placement = "bottom", options = list(container = "body"))
          ),
        HTML('<p>&nbsp;</p>'),
        fluidRow(
          ##---------------Step 6: Significant 3D numbers------------
          # uiOutput('show.top.stat.table'),
          box(title='Step 6: Significant 3D numbers',
              width = 12,status = 'primary', solidHeader = T,
              h4('Number of 3D genes/transcripts in each contrast group'),
              tableOutput('DDD.numbers')
          ),
          # ##---------- plot transcript per genes ------------
          # column(width = 3,
          #        box(title='Number of transcripts per gene',
          #            width = NULL,status = 'primary', solidHeader = F,
          #            spectrumInput(
          #              inputId = "tpg_color",
          #              label = "Pick a color for the bars:",
          #              choices = distinct.color(50),
          #              options = list(`toggle-palette-more-text` = "Show more"),
          #              width = "50%"
          #            ),
          #            numericInput(inputId = "tpg_number_x_hjust",label = 'Horizontal-just x-labs',value = '0.5',
          #                         min = 0,max = 1,step = 0.1,width = "100%"),
          #            numericInput(inputId = "tpg_number_x_vjust",label = 'Vertical-just x-labs',value = '0.5',
          #                         min = 0,max = 1,step = 0.1,width = "100%"),
          #            numericInput(inputId = "tpg_number_x_rotate",label = 'Rotate x-labs',value = '0',
          #                         min = 0,max = 360,step = 15,width = "100%"),
          #            bsTooltip(id = 'tpg_number_x_rotate',
          #                      title = "Adjust the labels on x-axis.",
          #                      placement = "bottom",
          #                      options = list(container = "body")),
          #            actionButton('tpg_plot_button','Plot',icon("send outline icon"),class="btn btn-primary",
          #                         style="color: #fff; background-color: #428bca; border-color: #2e6da4")
          #        )
          # ),
          #   column(width = 9,
          #          box(title= NULL,
          #              width = NULL,status = 'primary', solidHeader = T,
          #              plotOutput('tpg_plot'),
          #              hr(),
          #              column(12,
          #                     div(style="float:left;margin-left: 0px;",
          #                         numericInput(inputId = "tpg_number_width",label = 'Plot width (inch)',
          #                                      value = 8,width = "100%",step = 0.2)
          #                     ),
          #                     div(style="float:left;margin-left: 25px;",
          #                         numericInput(inputId = "tpg_number_height",label = 'Plot height (inch)',
          #                                      value = 5,width = "100%",step = 0.2)
          #                     ),
          #                     div(style="float:left;margin-left: 25px;",
          #                         numericInput(inputId = "tpg_number_res",label = 'PNG plot resolution',value = '150',
          #                                      min = 0,max = 600,step = 60,width = "100%")
          #                     ),
          #                     actionButton(inputId = 'tpg_save_button',label = 'Save',icon = icon('download',lib = 'font-awesome'),
          #                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
          #                                  float: left; margin-top: 25px; margin-left: 25px"),
          #                     actionButton(inputId = 'tpg_save_button_view',label = 'View saved plot',
          #                                  icon = icon('eye',lib = 'font-awesome'),
          #                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
          #                                  float: left; margin-top: 25px; margin-left: 25px"),
          #                     bsTooltip(id = "tpg_save_button_view",
          #                               title = "Preview saved plots to correct width and height.",
          #                               placement = "bottom", options = list(container = "body"))
          #                     )
          #              )
          #          ),
          column(width = 3,
                 box(title='Venn diagrams of 3D lists',
                     width = 13,status = 'primary', solidHeader = F,
                     selectInput(
                       inputId = 'across.contrast.group',
                       label = 'Contrast groups (multi-select allowed)',
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
                 tabBox(title = 'Venn diagrams',width = 13,
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
                 actionButton(inputId = 'save_3D_euler_plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: 
                              #2e6da4; float:left;margin-top:25px;margin-left:25px"),
                 bsTooltip(id = "save_3D_euler_plot", 
                           title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                           placement = "bottom", options = list(container = "body")),
                 actionButton(inputId = 'save_3D_euler_plot_view',label = 'View saved plot',
                              icon = icon('eye',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_3D_euler_plot_view",
                           title = "Preview saved plots to correct width and height.",
                           placement = "bottom", options = list(container = "body"))
          )
        ),
        HTML('<p>&nbsp;</p>'),
        fluidRow(
          column(width = 12,
                 tabBox(width = 13,
                        # Title can include an icon
                        # title = tagList(shiny::icon("gear"), "Loaded tables"),
                        title = "Number of up-down regulated targets",side = 'left',
                        tabPanel("DE genes",
                                 plotOutput("up.down.bar.plot.DE.genes")
                        ),
                        # tabPanel("DAS genes",
                        #          plotOutput("up.down.bar.plot.DAS.genes")
                        # ),
                        tabPanel("DE transcripts",
                                 plotOutput("up.down.bar.plot.DE.trans")
                        )
                        # tabPanel("DTU transcripts",
                        #          plotOutput("up.down.bar.plot.DTU.trans")
                        # )
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
                 actionButton(inputId = 'save_up_down_bar_plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_up_down_bar_plot", 
                           title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                           placement = "bottom", options = list(container = "body")),
                 actionButton(inputId = 'save_up_down_bar_plot_view',label = 'View saved plot',
                              icon = icon('eye',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_up_down_bar_plot_view",
                           title = "Preview saved plots to correct width and height.",
                           placement = "bottom", options = list(container = "body"))
          )
          ),
        HTML('<p>&nbsp;</p>'),
        fluidRow(
          column(width = 3,
                 box(title='Volcano plot',
                     width = 13,status = 'primary', solidHeader = F,
                     sliderInput(inputId = 'Volcano_top_label',label = 'Labels of top n points',
                                 min = 0,step = 1,max = 100,value = 10),
                     HTML('Top n labels of the points will be calculatd according to their distance to coordinate (0,0)'),
                     p(),
                     sliderInput(inputId = 'Volcano_point_size',label = 'Point size',
                                 min = 0,step = 0.3,max = 3,value = 1),
                     p(),
                     spectrumInput(
                       inputId = "color4Volcano0",
                       label = "Color of not significant points:",
                       choices = c('black',distinct.color(50)[1:50]),
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "200px"
                     ),
                     spectrumInput(
                       inputId = "color4Volcano1",
                       label = "Color of significant points:",
                       choices = c('red',distinct.color(50)[1:50]),
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "200px"
                     ),
                     p(),
                     actionButton('Volcano_plot_button','Plot',icon("send outline icon"),class="btn btn-primary",
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                     bsTooltip(id = "Volcano_plot_button",
                               title = "Depending on the data size, it may take quite a while to generate the plots.",
                               placement = "bottom", options = list(container = "body"))
                 )
          ),
          column(width = 9,
                 tabBox(width = 13,
                        # Title can include an icon
                        # title = tagList(shiny::icon("gear"), "Loaded tables"),
                        title = "Volcano plot",side = 'left',
                        tabPanel("DE genes",
                                 plotOutput("DE_genes_Volcano")
                        ),
                        tabPanel("DAS genes",
                                 plotOutput("DAS_genes_Volcano")
                        ),
                        tabPanel("DE transcripts",
                                 plotOutput("DE_trans_Volcano")
                        ),
                        tabPanel("DTU transcripts",
                                 plotOutput("DTU_trans_Volcano")
                        )
                 ),
                 div(style="float:left;margin-left: 0px;",
                     numericInput(inputId = "Volcano_width",label = 'Plot width (inch)',
                                  value = 8,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "Volcano_height",label = 'Plot height (inch)',
                                  value = 5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-left: 25px;",
                     numericInput(inputId = "Volcano_res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 actionButton(inputId = 'save_Volcano_plot',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_Volcano_plot", 
                           title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                           placement = "bottom", options = list(container = "body")),
                 actionButton(inputId = 'save_Volcano_plot_view',label = 'View saved plot',
                              icon = icon('eye',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4;
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_Volcano_plot_view",
                           title = "Preview saved plots to correct width and height.",
                           placement = "bottom", options = list(container = "body"))
          )
        ),
        HTML('<p>&nbsp;</p>'),
        ##---------------Step 7: Transcriptional vs alternative splicing regulation------------
        fluidRow(
          box(title='Step 7: Transcriptional vs alternative splicing regulation',
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
          column(width = 12,
                 tabBox(width = 13,title = 'Number of union set of all contrast groups',side = 'left',
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
                 actionButton(inputId = 'save_union_flow_chart',label = 'Save',
                              icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_union_flow_chart", 
                           title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                           placement = "bottom", options = list(container = "body")),
                 actionButton(inputId = 'save_union_flow_chart_view',label = 'View saved plot',
                              icon = icon('eye',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_union_flow_chart_view",
                           title = "Preview saved plots to correct width and height.",
                           placement = "bottom", options = list(container = "body"))
          ),
          HTML("<p>&nbsp;</p>"),
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
                        title = "",side = 'left',
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
                 actionButton(inputId = 'save_transvsAS_euler_plot',label = 'Save',
                              icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: 
                              #2e6da4; float:left; margin-top:25px;margin-left: 25px"),
                 bsTooltip(id = "save_transvsAS_euler_plot", 
                           title = 'Save plot to working directory. If App is not running locally, figures can be downloaded in the last page "Generate report"',
                           placement = "bottom", options = list(container = "body")),
                 actionButton(inputId = 'save_transvsAS_euler_plot_view',label = 'View saved plot',
                              icon = icon('eye',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                              float: left; margin-top: 25px; margin-left: 25px"),
                 bsTooltip(id = "save_transvsAS_euler_plot_view",
                           title = "Preview saved plots to correct width and height.",
                           placement = "bottom", options = list(container = "body"))
                 
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_ddd',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Data pre-processing</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Time-series trend</i>'),
              actionButton(inputId = 'page_after_ddd',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )


