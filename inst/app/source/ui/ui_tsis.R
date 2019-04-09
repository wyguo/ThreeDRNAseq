#========================>> TSIS <<=======================
tabItem('TSIS',
        fluidRow(
          box(title = 'Workflow',
              width=12,status = 'primary', solidHeader = T,
              HTML('<img style="width: 65%; display: block; margin-left: auto; margin-right: auto;" 
                   src="ISs.png"/>')
              )
          ),
        fluidRow(
          column(width = 12,
                 box(title='Isoform switch analysis',
                     width = NULL,status = 'primary', solidHeader = T,
                     HTML('<div align="justify">
                          Transcript isoform switches (ISs) occur when a pair of alternatively spliced isoforms reverse the order of their 
                          relative expression levels as shown in Figure 1.
                          <b>IsokTSP</b> is a method to detect transcript ISs between pair-wise conditions (Sebestyen et al., 2015) (Figure 1A) while 
                          <b>TSIS</b> (time-series IS) is used to identify ISs in time-series data (Figure 1B). 
                          To have robust results, isokTSP and TSIS characterize the transcript switch by 1) defining 
                          the isoform switch conditions/time-points for any pair of transcript isoforms within a gene, 2) describing the switch 
                          using five different features or metrics, 3) filtering the results with user provided specifications and 4) 
                          visualizing the results using different plots for the user to examine further details of the switches (see details in our publication Guo et al., 2017). 
                          </div>'),
                     HTML('<div align="justify">
                          <br>
                          <b>IsokTSP: </b>
                          In 3D RNA-seq analysis, isokTSP is used to identify ISs between pair-wise conditions that make up contrast groups of expression comparisons. To find IS 
                          point, the average TPMs of isoforms are treated as y-axis coordinates and 1 and 2 are used as x-axis coordinates at pair-wise conditions (Figure 1A). Then two lines,
                          which represent two isoforms, can be fitted based on these coordinates. 
                          The interaction (if exists) of two lines is defined as the IS point.</div>'),
                     HTML('<div align="justify">
                          <br>
                          <b>TSIS:</b> To find all IS points in time-series expression, average TPMs or spline fitted values of TPMs are used as y-axis coordinates and integers, 
                          1, 2, 3, ..., are x-axis coordinates. Lines of isoforms are fitted between each pair of consecutive time-points. The IS points can be identified from 
                          interactions of lines.</div>'),
                     HTML('<div align="justify">
                          <br>
                          <b>Significant IS:</b> Five metrics are used to evaluate the significance of ISs. (1) the probability of switch (i.e., the frequency of samples reversing 
                          their relative abundance at the switches); (2) the sum of the average TPM differences of the two isoforms in both intervals before and after the switch point; 
                          (3) adjusted p-value (i.e. FDR), which is the significance of the differences between the switched isoform abundances before and after the switch; 
                          (4) the time-points at both intervals before and after switch (set to 1 in isokTSP); and (5) Pearson correlation of two isoforms. Users can apply filters based
                          on these metrics to generate significant ISs (Guo et al., 2017).</div>'),
                     HTML('<p><img style="width: 50%; display: block; margin-left: auto; margin-right: auto;" 
                          src="TSIS.png"/></p>'),
                     HTML('<b>Figure 1:</b> Isoform switch analysis methods. Expression data with 3 replicates for each condition/time-point 
                          is simulated for isoforms \\(iso_{i}\\) and \\(iso_j\\) (blue and red circles). The points in the plots represent 
                          samples and the black lines connect the average of samples. (A) is the isokTSP method for comparisons 
                          of two conditions \\(c_1\\) and \\(c_2\\). The Time-Series Isoform Switch (TSIS) tool is designed for detection and 
                          characterization of isoform switches for time series data shown in (B). The time-series with 6 time-points is 
                          divided into 4 intervals by the intersection points of average expression. '),
                     hr(),
                     HTML('<b>If you use TSIS in your work, please cite:</b>'),
                     br(),
                     HTML('Wenbin Guo, Cristiane P. G. Calixto, John W.S. Brown, Runxuan Zhang, "TSIS: an R package to infer alternative 
                          splicing isoform switches for time-series data", Bioinformatics, https://doi.org/10.1093/bioinformatics/btx411, 2017.'),
                     br(),
                     HTML('<b>The paper for TSIS application in RNA-seq data from study of Arabidopsis in response to cold:</b>'),
                     br(),
                     HTML('Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.')
                     )
                     )
                     ),
        fluidRow(
          column(width = 12,
                 box(title='Step 1: Set parameters of IS analysis',
                     width = NULL,status = 'primary', solidHeader = T,
                     h4('Select TSIS or isokTSP'),
                     fluidRow(
                       column(12,
                              radioButtons(inputId = 'TSISorisokTSP',label = 'Select a type',
                                           choices = c('TSIS','isokTSP'),selected = 'isokTSP',inline = T),
                              HTML('<b>TSIS</b>: Time-Series Isoform Switch across sequencial time-points; 
                            <b>isokTSP</b>: Pair-Wise Isoform Switch between conditions of contrast groups.')
                              )
                     ),
                     hr(),
                     h4('Scoring parameters'),
                     fluidRow(
                       column(3,
                              uiOutput('show.method.intersection')
                       ),
                       column(3,
                              uiOutput('show.spline.df')
                       ),
                       p(),
                       column(12,
                           actionButton('scoring','Scoring',icon("send outline icon"),class="btn btn-primary",
                                        style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                           )
                           
                     ),
                     HTML('<div align="justify">
                          Press Scoring button to implement the scoring of isoform switches. The details of parameters:
                          <ul><li><b>Method for intersections:</b> Using either mean values or natural spline fitted smooth curves of time-series expression
                          to determine the intersection points of isoforms.</li>
                          <li><b>Degree of spline:</b> If use spline method, the higher degree leads to more break points of the time-series.</li>
                          </ul></div>'),
                     hr(),
                     h4('Filtering parameters'),
                     fluidRow(
                       column(3,
                              numericInput('TSIS_prob_cut',label='Probability cutoff:',value=0.5)
                       ),
                       column(3,
                              numericInput('TSIS_diff_cut',label='Difference cutoff:',value=1)
                       ),
                       column(3,
                              numericInput('TSIS_adj_pval_cut',label='Adjusted p-value cutoff:',value=0.01)
                       ),
                       column(3,
                              #
                              selectInput('TSIS_time_point_cut',label='Min time in interval:',choices = 1:10,selected = 2,multiple = F)
                       ),
                       column(3,
                              numericInput('TSIS_cor_cut',label='Correlation cutoff:',value=0)
                       ),
                       column(3,
                              div(style="display:inline-block;vertical-align:middle; margin-top:25px",
                                  actionButton('filtering','Filtering',icon("send outline icon"),class="btn btn-primary",
                                               style="color: #fff; background-color: #428bca; border-color: #2e6da4"))
                       )
                     ),
                     column(12,
                            br()
                     ),
                     HTML('Press Filtering button to filter the scores. The details of parameters:
                          <ul><li><b>Probability cutoff:</b> The isoform switch probability/frequency cut-off for the column "prob" in the output table. </li>
                          <li><b>Difference cutoff:</b> The isoform switch difference cut-off for the column "diff" in the output table. </li>
                          <li><b>P-value cutoff:</b> The p-value cut-off of both columns "before.pval" and "after.pval" in the output table. </li>
                          <li><b>Min time in interval:</b> The minimum time points for both columns "before.t.points" and "after.t.points" in the output table.</li>
                          <li><b>Correlation cutoff:</b> The cut-off for Pearson correlation of isoform pairs.</li>
                          </ul>')
                     )
                     ),
          column(width = 12,
                 box(title='TSIS results',
                     width = NULL,status = 'primary', solidHeader = F,
                     DT::DTOutput("score.table"),
                     # actionButton(inputId = 'refresh_scores_table',label = 'Refresh the table',icon = icon('refresh'),style='primary'),
                     HTML('The following table shows the feature scores of isoform switches. The columns in the output table:
                          <ul><li><b>iso1, iso2:</b> the isoform pairs. </li>
                          <li><b>iso1.mean.ratio, iso2.mean.ratio:</b> The mean ratios of isoforms to their gene. </li>
                          <li><b>before.interval, after.interval:</b> The intervals before and after switch points. </li>
                          <li><b>x.value, y.value:</b> The values of x axis (time) and y axis (expression) coordinates of the switch points.</li>
                          <li><b>prob:</b> The probability/frequency of switch.</li>
                          <li><b>diff:</b> The sum of average sample differences before and after switch.</li>
                          <li><b>before.pval, after.pval:</b> The paired t-test p-values of the samples in the intervals before and after switch points.</li>
                          <li><b>before.t.points, after.t.points:</b> The number of time points in intervals before and after the switch points.</li>
                          <li><b>cor:</b> The Pearson correlation of iso1 and iso2.</li>
                          </ul>')
                     )
                     )
                     ),
        fluidRow(
          ##---------- plot switch number ------------
          column(width = 3,
                 box(title='Step 2: Number of significant switches',
                     width = NULL,status = 'primary', solidHeader = T,
                     # radioButtons("densityplot.type", label = HTML("<b>Select plot type:</b>"),
                     #              choices = list("Frequency" = 'frequency',"Density" = 'density'),
                     #              selected = 'frequency',inline = F),
                     # radioButtons("show.density.line", label = HTML("<b>Density/frequency line:</b>"),
                     #              choices = list("FALSE" = 'FALSE',"TRUE" = 'TRUE'),
                     #              selected = 'FALSE',inline = F),
                     spectrumInput(
                       inputId = "switch.density.color.pick",
                       label = "Pick a color for the bars:",
                       choices = distinct.color(50),
                       options = list(`toggle-palette-more-text` = "Show more"),
                       width = "50%"
                     ),
                     numericInput(inputId = "is_number_x_hjust",label = 'Horizontal-just x-labs',value = '0.5',
                                  min = 0,max = 1,step = 0.1,width = "100%"),
                     numericInput(inputId = "is_number_x_vjust",label = 'Vertical-just x-lab',value = '0.5',
                                  min = 0,max = 1,step = 0.1,width = "100%"),
                     numericInput(inputId = "is_number_x_rotate",label = 'Rotate x-lab',value = '0',
                                  min = 0,max = 360,step = 15,width = "100%"),
                     bsTooltip(id = 'is_number_x_rotate',
                               title = "Adjust the labels on x-axis.",
                               placement = "bottom", 
                               options = list(container = "body")),
                     actionButton('plotISnumber','Plot',icon("send outline icon"),class="btn btn-primary",
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4")
                 )
          ),
          column(width = 9,
                 box(title= NULL,
                     width = NULL,status = 'primary', solidHeader = T,
                     plotOutput('switch.number'),
                     hr(),
                     column(12,
                            div(style="float:left;margin-left: 0px;",
                                numericInput(inputId = "is.number.width",label = 'Plot width (inch)',
                                             value = 8,width = "100%",step = 0.2)
                            ),
                            div(style="float:left;margin-left: 25px;",
                                numericInput(inputId = "is.number.height",label = 'Plot height (inch)',
                                             value = 5,width = "100%",step = 0.2)
                            ),
                            div(style="float:left;margin-left: 25px;",
                                numericInput(inputId = "is.number.res",label = 'PNG plot resolution',value = '150',
                                             min = 0,max = 600,step = 60,width = "100%")
                            ),
                            actionButton(inputId = 'plotISnumber.save',label = 'Save',icon = icon('download',lib = 'font-awesome'),
                                         style="color: #fff; background-color: #428bca; border-color: #2e6da4; 
                                         float: left; margin-top: 25px; margin-left: 25px")
                            )
                     )
                 )
          ),
        ##---------- plot isoform switch ------------
        fluidRow(
          box(title= 'Step 3: Plot significant switches',
              width = 3,status = 'primary', solidHeader = T,
              h4('Input isoform names:'),
              textInput('iso1', label = 'Isoform 1: ', value = ''),
              textInput('iso2', label = 'Isoform 2: ', value = ''),
              uiOutput('show.contrst.button'),
              selectInput('ribbon.plot','Plot types',c('Error bar','Ribbon')),
              radioButtons("show.scores", label = h4("Show statistics:"),
                           choices = list("TRUE" = 'TRUE', "FALSE" = "FALSE"),
                           selected = 'TRUE',inline = T),
              div(style="display: inline-block; margin-right: 5px;",
                  spectrumInput(
                    inputId = "color.iso1",
                    label = "Color for iso1:",
                    choices = distinct.color(50),
                    options = list(`toggle-palette-more-text` = "Show more"),
                    width = "100%"
                  )
              ),
              div(style="display: inline-block; margin-right: 5px;",
                  spectrumInput(
                    inputId = "color.iso2",
                    label = "Color for iso2:",
                    choices = distinct.color(51)[-1],
                    options = list(`toggle-palette-more-text` = "Show more"),
                    width = "100%"
                  )
              ),
              br(),
              actionButton('plot.IS','Plot',icon("send outline icon"),class="btn btn-primary",
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4")
          ),
          box(title= NULL,
              width = 9,status = 'primary', solidHeader = F,
              plotOutput('IS.plot',width = 'auto',height = '550px'),
              br(),
              div(style="float:left;margin-right: 25px;",
                  numericInput(inputId = "tsis.plot.width",label = 'Plot width (inch)',
                               value = 8,width = "100%",step = 0.2)
              ),
              div(style="float:left;margin-right: 25px;",
                  numericInput(inputId = "tsis.plot.height",label = 'Plot height (inch)',
                               value = 5,width = "100%",step = 0.2)
              ),
              div(style="float:left;margin-right: 25px;",
                  numericInput(inputId = "tsis.plot.res",label = 'PNG plot resolution',value = '150',
                               min = 0,max = 600,step = 60,width = "100%")
              ),
              actionButton(inputId = 'save_tsis_plot',label = 'Save',
                           icon = icon('download',lib = 'font-awesome'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4; float: left; margin-top: 25px")
              # bsTooltip(id = "save_tsis_plot", 
              #           title = "The plot will be saved to working directory. Users can download it at final step.",
              #           placement = "bottom", options = list(container = "body")),
          )
        ),
        fluidRow(
          box(title= 'Plot all significant isoform switches',
              width = 6,status = 'primary', solidHeader = F,
              column(width = 12,
                     verbatimTextOutput('number_of_IS'),
                     radioButtons(inputId = 'tsis.multiple.plot.format',label = 'Format',
                                  choices = c('all','png','pdf'),selected = 'all',inline = T),
                     div(style="float:left;margin-right: 25px;",
                         numericInput(inputId = "tsis.all.plot.width",label = 'Plot width (inch)',
                                      value = 8,width = "100%",step = 0.2)
                     ),
                     div(style="float:left;margin-right: 25px;",
                         numericInput(inputId = "tsis.all.plot.height",label = 'Plot height (inch)',
                                      value = 5,width = "100%",step = 0.2)
                     ),
                     div(style="float:left;margin-right: 25px;",
                         numericInput(inputId = "tsis.all.plot.res",label = 'PNG plot resolution',value = '150',
                                      min = 0,max = 600,step = 60,width = "100%")
                     )
              ),
              column(width = 12,
                     actionButton(inputId = 'plot_all_IS_btn',label = 'Plot all',
                                  style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
                     br(),
                     br(),
                     HTML('Note: it may take quite a while to save all the plots. Please wait.')
              )
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-left:15px; float: left;",
              actionButton(inputId = 'page_before_tsis',label = '',icon = icon('arrow-left'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;"),
              HTML('<i>Functional plot</i>')
          ),
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Generate report</i>'),
              actionButton(inputId = 'page_after_tsis',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )