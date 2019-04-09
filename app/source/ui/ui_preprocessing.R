#=======================>> Data pre-processing <<============================
tabItem('preprocessing',
        ##----------Workflow------------
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 100%; display: block; margin-left: auto; margin-right: auto;" 
                   src="data_processing.png"/>')
              )
          ),
        fluidRow(
          column(width = 7,offset = 0, 
                 style='padding:0px; padding-right:0px; padding-left:0px,margin-left: 0px; margin-right:0px',
                 ##----------Step 1: Filter low expression------------
                 box(title = 'Step 1: Filter low expressed transcripts and genes',
                     width = 12,status = 'primary',solidHeader = T,
                     sliderInput(inputId = 'cpm_cut',label = 'CPM cut-off',
                                 min = 0,max = 10,value = 1),
                     bsTooltip(id = "cpm_cut", 
                               title = "Set a count per million reads (CPM) cut-off",
                               placement = "bottom", options = list(container = "body")),
                     sliderInput(inputId = 'cpm_samples_n',label = 'Sample number cut-off',
                                 min = 0,max = 10,value = 1),
                     bsTooltip(id = "cpm_samples_n", 
                               title = "Set the minimum numbers of samples to apply the CPM cut-off",
                               placement = "bottom", options = list(container = "body")),
                     div(style="float:left; margin-right: 0%",
                         actionButton('run_filter','Filter & Mean-variance trend plot',icon("filter"),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                         # actionButton('mv_plot_button','Mean-variance plot',icon("pencil"),
                         #              style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                     ) %>%
                       helper(icon = "question-circle",
                              colour = NULL,
                              type = 'markdown',
                              title = 'Mean-variance trend plot to determine low expression cut-offs',
                              size = 'l',
                              content = 'filter_low',
                              style="font-size: 2.0rem;margin-right:50px;"),
                     bsTooltip(id = "run_filter", 
                               title = 'Click "Filter & Mean-variance trend plot" button to apply low exprssion filter and make mean-variance trend plot',
                               placement = "bottom", options = list(container = "body")),
                     HTML('&nbsp'),
                     br(),
                     hr(),
                     HTML('<strong>Note:</strong> (1) An expressed transcript must have &#x2265 <i>n</i> samples &#x2265 <i>m</i> CPM (Count Per Million reads) expression.
                          (2) An expressed gene must have at least one expressed transcript.') 
                     
                     )
                 ),
          column(width = 5,offset = 0, 
                 style='padding:0px; padding-right:0px; padding-left:0px,margin-left: 0px; margin-right:0px',
                 box(title = NULL,
                     width=12,status = 'primary', solidHeader = T,
                     HTML('RNA-seq data information'),
                     DT::dataTableOutput("RNAseq_datainfo")
                 )
          ),
          column(width = 12,
                 tabBox(title = 'Mean-variance trend plot',id = 'mv.plot', width = 13,
                        tabPanel(title = 'Transcript level',
                                 plotOutput(outputId = 'mv.trans.plot',click ='mv.trans.plot.click')
                        ),
                        tabPanel(title = 'Gene level',
                                 plotOutput(outputId = 'mv.genes.plot',click ='mv.genes.plot.click')
                        )
                 ),
                 div(style="float:left;margin-right: 25px;",
                     numericInput(inputId = "mv.plot.width",label = 'Plot width (inch)',
                                  value = 11.5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-right: 25px;",
                     numericInput(inputId = "mv.plot.height",label = 'Plot height (inch)',
                                  value = 4.5,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-right: 25px;",
                     numericInput(inputId = "mv.plot.res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 actionButton(inputId = 'save.mv.plot',label = 'Save',
                              icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; margin-top: 25px;")
          )
          ),
        br(),
        br(),
        ##----------Step 2: Principal component analysis------------
        fluidRow(
          column(width = 3,
                 box(title = 'Step 2: Principal component analysis (PCA)',width = 13,
                     status = 'primary',solidHeader = T,
                     div(style="display:inline-block;vertical-align:middle;width:120px",
                         selectInput(inputId = 'dim1',label = 'X-axis',width = '100%',
                                     choices = c('PC1','PC2','PC3','PC4'),selected = 'PC1')
                     ),
                     div(style="display:inline-block;vertical-align:middle;width:120px",
                         selectInput(inputId = 'dim2',label = 'Y-axis',width = '100%',
                                     choices = c('PC1','PC2','PC3','PC4'), selected = 'PC2')
                     ),
                     br(),
                     radioButtons(inputId = 'pca.plot.type',label = 'Plot type',
                                  choices = c('all samples','average expression'),
                                  selected = 'all samples',inline = T),
                     selectInput(inputId = 'pca_color_option',
                                 label = 'Colour samples based on',multiple = T,
                                 choices = NULL),
                     radioButtons(inputId = 'pca.add.circle',label = 'Add polygon to clusters',
                                  choices = c('none','ellipse','polygon'),
                                  selected = 'polygon',inline = T),
                     div(style="display:inline-block;vertical-align:middle;height:30px;",
                         actionButton('plot.pca.button','Plot',icon("pencil",lib = 'font-awesome'),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                     )%>%
                       helper(icon = "question-circle",
                              colour = NULL,
                              type = 'markdown',
                              title = 'PCA plot and batch effect estimation',
                              size = 'l',
                              content = 'PCA_and_batch_effects',
                              style="font-size: 2.0rem;margin-right:50px;"),
                     br(),
                     HTML('<ol>
                          <li style="margin-left: -20px;">Select PCs (X- and Y-axis) to visualise</li>
                          <li style="margin-left: -20px;">Select "Plot type", either to show all samples or average expression of conditions</li>
                          <li style="margin-left: -20px;">Select a factor/facotrs to colour samples into groups (multi-select allowed)</li>
                          <li style="margin-left: -20px;">Select "ellipse" or "polygon" to highlight the groups</li>
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
                 column(width = 12,
                        div(style="float:left;margin-right: 25px; margin-left: 0px",
                            numericInput(inputId = "pca.plot.x.limit.lower",label = 'x lower limit',
                                         value = -200,width = "100%",step = 0.2)
                        ),
                        div(style="float:left;margin-right: 25px;",
                            numericInput(inputId = "pca.plot.x.limit.upper",label = 'x upper limit',
                                         value = 200,width = "100%",step = 0.2)
                        )
                        ),
                 column(width = 12,
                        div(style="float:left;margin-right: 25px;",
                            numericInput(inputId = "pca.plot.width",label = 'Plot width (inch)',
                                         value = 9,width = "100%",step = 0.2)
                        ),
                        div(style="float:left;margin-right: 25px;",
                            numericInput(inputId = "pca.plot.height",label = 'Plot height (inch)',
                                         value = 8,width = "100%",step = 0.2)
                        ),
                        div(style="float:left;margin-right: 25px;",
                            numericInput(inputId = "pca.plot.res",label = 'PNG plot resolution',value = '150',
                                         min = 0,max = 600,step = 60,width = "100%")
                        ),
                        actionButton(inputId = 'save.pca.plot',label = 'Save',
                                     icon = icon('download',lib = 'font-awesome'),
                                     style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; margin-top: 25px")
                        )
          )
          ),
        br(),
        br(),
        ##----------Step 3: Batch effect estimation------------
        fluidRow(column(width = 3,
                        box(title = 'Step 3: Estimate batch effects',width = 13,
                            status = 'primary',solidHeader = T,
                            radioButtons(inputId = 'has_batcheffect',label = 'Data has batch effects?',
                                         choices = c('No','Yes'),selected = 'No',inline = T),
                            HTML('<h5><strong><font color="red">This step is only executed when there are distinct 
                                 bath effects in the data:
                                 <ul style="padding-left: 20px">
                                 <li>The biological replicates of the same conditions stay in separate clusters in the PCA plot.</li/>
                                 <li> There is experimental explanation of this separation (e.g. bio-reps in different time/labs).</li>
                                 </ul>
                                 Otherwise skip to Step 4 to avoid data over-modification.</font></strong></h4>'),
                            p() %>%
                              helper(icon = "question-circle",
                                     colour = NULL,
                                     type = 'markdown',
                                     title = 'PCA plot and batch effect estimation',
                                     size = 'l',
                                     content = 'PCA_and_batch_effects',
                                     style="font-size: 2.0rem;margin-right:50px;"),
                            radioButtons(inputId = 'RUVseq_method',label = 'RUVSeq method',
                                         choices = c('RUVr','RUVs','RUVg'),selected = 'RUVr',
                                         inline = T),
                            # actionButton('remove_batch','Remove batch effects',icon("cut"),
                            #              style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                            bsButton("remove_batch", label = "Remove batch effects", icon = icon("ban"),
                                     style="primary"),
                            bsTooltip(id = "remove_batch", 
                                      title = "Please make sure there are batch effects between biological replicates before click this button",
                                      placement = "bottom", options = list(container = "body")),
                            br(),
                            br(),
                            bsButton('update_pca_plot','Update PCA plot',icon("ban"),
                                         style="primary"),
                            bsTooltip(id = "update_pca_plot", 
                                      title = "The parameters to plot PCA are inherited from previous PCA panel.",
                                      placement = "bottom", options = list(container = "body"))
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
        
        ##----------Step 4: Normalization------------
        fluidRow(
          column(width = 3,
                 box(title = 'Step 4: Data normalization',width = 13,
                     status = 'primary',solidHeader = T,
                     radioButtons(inputId = 'norm_method',label = 'Normalization methods',
                                  choices = c('TMM','RLE','upperquartile'),
                                  selected = 'TMM',inline = T),
                     actionButton('run_norm','Run',icon("send outline icon",lib = 'font-awesome'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                     actionButton('norm.plot.button','Plot distribution',icon("pencil"),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4"),
                     br(),
                     br(),
                     HTML('<div align="justify">
                          Three normalisation methods, TMM (weighted trimmed mean of M-values), RLE (relative log expression) and 
                          upperquartile (upper-quartile), have comparable performance (Maza, 2016). TMM is 
                          more widely used in differential expression studies of RNA-seq data.</div>')
                     # br(),
                     # br(),
                     # verbatimTextOutput("norm.info")
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
                 div(style="float:left;margin-right: 25px;",
                     numericInput(inputId = "dist.plot.res",label = 'PNG plot resolution',value = '150',
                                  min = 0,max = 600,step = 60,width = "100%")
                 ),
                 div(style="float:left;margin-right: 25px;",
                     numericInput(inputId = "dist.plot.height",label = 'Plot height (inch)',
                                  value = 7,width = "100%",step = 0.2)
                 ),
                 div(style="float:left;margin-right: 25px;",
                     numericInput(inputId = "dist.plot.width",label = 'Plot width (inch)',
                                  value = 8,width = "100%",step = 0.2)
                 ),
                 actionButton(inputId = 'save.dist.plot',label = 'Save',
                              icon = icon('download',lib = 'font-awesome'),
                              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; margin-top: 25px"),
                 br(),
                 br(),
                 br(),
                 br()
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_preprocessing',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Data generation</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>3D analysis</i>'),
              actionButton(inputId = 'page_after_preprocessing',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )