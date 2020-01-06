
#========================>> report <<=======================
tabItem('report',
        fluidRow(
          # box(title = '3D RNA-seq analysis parameter summary',
          #     width=6,status = 'primary', solidHeader = T,
          #     HTML('The 3D analysis parameters in this table will be used to generate a report. If any parameters are incorrect, please double click corresponding
          #          cells in the table and type in the correct values. Then click "Update parameters" button to refresh the parameters for report.'),
          #     br(),
          #     br(),
          #     DT::dataTableOutput('params_table_panel'),
          #     br(),
          #     actionButton(inputId = 'update_params_list',label = 'Update parameters',
          #                  icon("refresh",lib="font-awesome"),
          #                  style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: right")
          #     ),
          box(title = 'Save data and results',
              width=6,status = 'primary', solidHeader = T,
              # tags$head(tags$style("#TxtOut {white-space: normal;}")),
              column(width = 12,
                     wellPanel(
                       actionButton(inputId = 'save_ddd_data_button',label = 'Save all 3D analysis results',
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4;",
                                    icon = icon('download')),
                       hr(),
                       HTML('Results are saved to 3D App working directory: 
                            <ul><li>
                            Intermediate data are saved in "data" folder.
                            </li>
                            <li>
                            Significant 3D lists and statistics in .csv (comma delimited) are saved in "result" folder.
                            </li>
                            <li>
                            Figures are saved in "figure" folder.
                            </li>
                            </ul>')
                     )
              ),
              column(width = 12,
                     wellPanel(
                       actionButton(inputId = 'generate_report',label = 'Generate report',
                                    style="color: #fff; background-color: #428bca; border-color: #2e6da4;",
                                    icon = icon('download')),
                       hr(),
                       HTML('Reports in word, pdf and html formats are saved in "report" folder in 3D App working directory ')
                       # br(),
                       # br(),
                       # HTML('Reports in html, word and pdf format are saved in "report" folder. Reports are also availabe to download by following links:'),
                       # br(),
                       # tags$a(href='report/3D_report.html', '3D RNA-seq analysis report.html',target='_blank',download = '3D_report.html'),
                       # br(),
                       # tags$a(href='report/3D_report.pdf', '3D RNA-seq analysis report.pdf', target='_blank',download = '3D_report.pdf'),
                       # br(),
                       # tags$a(href='report/3D_report.docx', '3D RNA-seq analysis report.docx', target='_blank',download = '3D_report.docx')
                     )
              ),
              column(width = 12,
                     wellPanel(
                       # actionButton("act_download_zip_results", label = "Click to download results", 
                       #              icon = icon("download"),
                       #              style="color: #fff; background-color: #428bca; border-color: #2e6da4; float:left;"),
                       downloadButton('download_zip_results', 'Click to download results',class="btn btn-primary",
                                      style="color: #fff; background-color: #428bca; border-color: #2e6da4;float:left"),
                       br(),
                       br(),
                       hr(),
                       HTML('<h5 align="justify">If the 3D RNA-seq App is running on a local PC, all the results (figures, tables, reports and intermediate .RData objects) are already saved in the working directory.
                            If the App is running on our server, all the results can be zipped and downloaded by clicking above button.</h5>')
                       )
                     )   
              )
          ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_report',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>TSIS</i>')
          )
        )
        )