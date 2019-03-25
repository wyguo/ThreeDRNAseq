#=======================>> Data generation panel <<============================
tabItem("generation",
        tags$head(tags$style("#TxtOut {white-space: normal;}")),
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 100%; display: block; margin-left: auto; margin-right: auto;" 
                   src="data_generation.png"/>')
              )
          ),
        fluidRow(
          box(title = 'App working directory',
                    width=12,status = 'primary', solidHeader = T,
          textOutput('data_folder_path_text')
          )
        ),
        ###upload intermediate data is disabled
        # fluidRow(
        #   actionButton(inputId = 'load_intermediate_data_btn',label = 'Load intermediate data')
        # ),
        fluidRow(
          ##----------Step 1: Use a different folder to save results?------------
          box(title = 'Step 1: Use a different folder to save results?',
              width=6,status = 'primary', solidHeader = T,
              column(12,
                     # wellPanel(
                     # style = "background-color: #ffffff;",
                     HTML('<strong>Folder to save results:</strong>'),
                     # textOutput('data_folder_path_text'),
                     tags$p(),
                     radioButtons(inputId = 'use_diff_folder',label = 'Use a different folder',
                                  choices = c('No','Yes'),selected = 'No',inline = T),
                     uiOutput('data_folder_button_ui')
                     # )
              )
          ),
          ##----------Step 2: How to generate data of analysis?------------
          box(title = 'Step 2: How to generate data of analysis?',
              width=6,status = 'primary', solidHeader = T,
              column(12,
                     # wellPanel(
                     # style = "background-color: #ffffff;",
                     radioButtons(inputId = 'generate_new_data',label = 'How to generate data?',
                                  choices = c('Generate new data','Upload intermediate data'),selected = 'Generate new data',inline = T),
                     uiOutput('upload_data_ui'),
                     HTML('<h5 align="justify">If you have already run the 3D RNA-seq analysis once and the intermediate data has been saved as
                          "intermediate_data.RData" object, which can be uploaded here. <font color="red"><strong>By doing this, you can skip the following
                          data genertion steps and visualise/redo the analysis directly from "Data pre-processing" until to "Generate report" step.</strong></font></h4>')
                     # )
                     )
              )
          ),
        fluidRow(
          ##----------Step 3: Inputs of 3D analysis------------
          box(title = 'Step 1: Input data of 3D analysis',
              width=6,status = 'primary', solidHeader = T,
              column(width = 12,
                     # wellPanel(
                     # style = "background-color: #ffffff;",
                     p() %>%
                       helper(icon = "question-circle",
                              colour = NULL,
                              type = 'markdown',
                              title = 'Input data of 3D RNA-seq App',
                              size = 'l',
                              content = 'input_data',
                              style="font-size: 2.0rem;margin-right:50px;"),
                     br(),
                     fileInput("sample_file_button", ">Select sample information csv file (comma delimited)",
                               accept = c(
                                 ".csv")
                     ),
                     div(style="display:inline-block;vertical-align:middle;margin-left:15px",
                         textOutput('sample_file_path_text')
                     )
                     # )
              ),
              column(width = 12,
                     # wellPanel(
                     # style = "background-color: #ffffff;",
                     radioButtons(inputId = 'mapping_file_type',label = 'Transcript-gene mapping in',
                                  choices = c('csv (comma delimited)','gtf'),selected = 'csv (comma delimited)',inline = T),
                     fileInput("mapping_file_button", ">Select transcript-mapping file",
                               accept = c(
                                 ".csv",".gtf")
                     ),
                     div(style="display:inline-block;vertical-align:middle;margin-left:15px",
                         textOutput('mapping_file_path_text')
                     )
                     # )
              ),
              column(width = 12,
                     HTML('><strong> Upload transcript quantification</strong> '),
                     br(),  
                     br(),
                     fileInput("quant_zip_button", " >>(1): Select quantification zip file to upload",
                               accept = c(
                                 ".zip")
                     ),
                     HTML('<strong>Note:</strong> if the App is run in a local PC and the transcript quantifcation folder is already in the App working
                          directory, please skip step (1) and select the quantification folder in step (2) directly.'),
                     br(),
                     br(),
                     HTML('<strong> >> (2) Select directory of transcript quantification folder</strong>'),
                     br(),
                     shinyDirButton(id = 'quant_folder_button',label = 'Select a folder',title = 'Select the folder of quantification',
                                    buttonType = 'primary',
                                    icon = icon('folder')),
               
                     br(),
                     br(),
                     verbatimTextOutput('qaunt_folder_path_text')
              )
          ),
          ##----------Step 4: Select factors of interest------------
          box(title = 'Step 2: Select factors of experimental design',
              width=6,status = 'primary', solidHeader = T,
              column(width = 12,
                     radioButtons(inputId = 'has_srep',label = '>Does the data have sequencing replicates?',
                                  choices = c('No','Yes'),selected = 'No',inline = T) %>%
                       helper(icon = "question-circle",
                              colour = NULL,
                              type = 'markdown',
                              title = 'Filter low expression based mean-variance trend plot',
                              size = 'l',
                              content = 'select_factor',
                              style="font-size: 2.0rem;margin-right:50px;"),
                     selectInput(inputId = 'factor_column',label = '>Select factor column/columns of interest (multi-select allowed)',
                                 choices = NULL,
                                 selected = NULL,
                                 multiple = T,
                                 width = NULL),
                     selectInput(inputId = 'brep_column',label = '>Select biological replicate column',
                                 choices = NULL,
                                 selected = NULL,
                                 multiple = F,
                                 width = NULL),
                     uiOutput('show_srep_column'),
                     selectInput(inputId = 'quant_column',
                                 label = '>Select quantification folder name column',
                                 choices = NULL,
                                 selected = NULL,
                                 multiple = F,
                                 width = NULL),
                     radioButtons(inputId = 'quant_method',label = '>Transcripts are quantified by:',
                                  choices = c('salmon','kallisto'),
                                  selected = 'salmon',inline = T),
                     actionButton(inputId = 'samples_update',label = 'Add sample information',icon = icon('update'),
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left"),
                     bsTooltip(id = "samples_update", 
                               title = "Click to add selected information to sample table",
                               placement = "bottom", options = list(container = "body"))
              )
          )
        ),
        fluidRow(
          ##----------Step 5: Generate read counts and TPMs------------
          box(title='Step 3: Generate read counts and TPMs',
              width = 6,status = 'primary', solidHeader = T,
              column(width = 12,
                     div(style="display:inline-block;vertical-align:middle",
                         selectInput(inputId = "tximport_method",
                                     label = "Choose a tximport method:",
                                     choices = c("lengthScaledTPM", "scaledTPM", "no"),
                                     selected = "lengthScaledTPM",
                                     width = 200
                         )
                     ),
                     # bsTooltip(id = "tximport_method", 
                     #           title = "Set a count per million reads (CPM) cut-off",
                     #           placement = "bottom", options = list(container = "body")),
                     div(style="display:inline-block;vertical-align:middle;height:30px",
                         actionButton('tximport_run','Run',icon("send outline icon"),
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                     ),
                     HTML('<ul style="padding-left: 20px">
                          <li>lengthScaledTPM: scaled using the average transcript length over samples and then the library size (<font color="red">recommended</font>). </li>
                          <li>scaledTPM: scaled up to library size.</li>
                          <li>no: no adjustment.</li>
                          </ul>'),
                     HTML('More details can be found in the <a href="http://www.bioconductor.org/packages/release/bioc/manuals/tximport/man/tximport.pdf" target="_blank">tximport user manual</a>.')
                     )
              )
              ),
        fluidRow(
          # box(title = 'Workflow',
          #     width=6,status = 'primary', solidHeader = T,
          #     HTML('<img style="width: 100%; display: block; margin-left: auto; margin-right: auto;" 
          #          src="contrast_design.png"/>')
          #     ),
          tabBox(title = 'Visualise sample information and transcript-gene mapping',id = 'uploaded_data', width = 12,side = 'right',selected = 'Sample information',
                 tabPanel(title = 'Part of mapping',
                          DT::dataTableOutput('mapping_table')
                 ),
                 tabPanel(title = 'Sample information',
                          DT::dataTableOutput('sample_table')
                 )
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
          actionButton(inputId = 'page_before_generation',label = '',icon = icon('arrow-left'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
          HTML('<i>Introduction</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
          HTML('<i>Data pre-processing</i>'),
          actionButton(inputId = 'page_after_generation',label = '',icon = icon('arrow-right'),
                       style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )