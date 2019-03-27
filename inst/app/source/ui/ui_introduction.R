tabItem("introduction",
        fluidRow(
          column(width = 12,
                 box(title='Introduction',
                     width = NULL,status = 'primary', solidHeader = T,
                     # shiny::uiOutput('page1'),
                     # DT::dataTableOutput("info")
                     HTML('<img style="width: 40%; display: block; margin-left: auto; margin-right: auto;" src="logo.png"/>'),
                     br(),
                     HTML('
                          
<div align="justify">
                          <strong>The ThreeDRNAseq (3D RNA-seq)</strong> R package provides an interactive graphical user interface (GUI) for RNA-seq data differential expression (DE), differential alternative splicing (DAS) and differential transcript usage (DTU) (3D) analyses based on two popular pipelines: <a href="https://bioconductor.org/packages/release/bioc/html/limma.html" target="_blank">limma</a> (Ritchie et al., 2015) and <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">edgeR</a> (Robinson et al., 2010). The 3D RNA-seq GUI is based on R <a href="https://shiny.rstudio.com/" target="_blank">shiny</a> App and enables a complete RNA-seq analysis to be done within only <strong>3 Days (3D)</strong>. The input required by 3D RNA-seq is a folder containing transcript quantifications. Therefore, to perform the 3D analysis,
                          <ul>
                          <li>The first step is to generate transcript quantifications from tools such as <a href="https://combine-lab.github.io/salmon/" target="_blank">Salmon</a> (Patro et al., 2017) and <a href="https://pachterlab.github.io/kallisto/" target="_blank">Kallisto</a> (Bray et al., 2016). Depending on the size of the dataset, the transcript quantification procedure often takes <strong>1-2 Days</strong>.</li>
                          <li>The second step is to upload transcript read counts and perform 3D analysis on the App. The 3D RNA-seq analysis pipeline starts with a number of steps to pre-process the data and reduce noise (e.g. low expression filters based on expression mean-variance trend, visualise data variation, batch effect estimation, data normalisation, etc.). The optimal parameters for each step are determined by quality control plots (e.g. mean-variance trend plots, PCA plots, data distribution plots, etc.). Stringent filters are applied to statistical testing of 3D genes/transcripts to control false positives. After the analysis, publication quality plots (e.g. expression profile plots, heatmap, GO annotation plots, etc.) and reports can be generated. The entire 3D RNA-seq analysis takes only <strong>1 Day or less</strong> and all actions are performed by simple mouse clicks on the App.</li>
                          </ul>
                          The 3D RNA-seq App is designed for use by biologists to analyse their own RNA-seq data or for bioinformaticians with little experience of gene/transcript expression analysis. It handles simple pair-wise analyses and complex experimental designs such as time-series, developmental series and multiple conditions. It employs the state-of-the-art methods/statistics and generate 3D results quickly and accurately.
                          <p></p>
                          The 3D RAN-seq App was developed by Dr. Wenbin Guo from research into the analysis of time-series RNA-seq data (Calixto et al., 2018) with help from Dr. Cristiane Calixto and Nikoleta Tzioutziou and guidance from Prof. John W.S. Brown and Dr. Runxuan Zhang from the University of Dundee - School of Life Sciences and the James Hutton Institute -  Information and Computational Sciences. We acknowledge Dr. Iain Milne and Gordon Stephen for technical support. To use our pipeline in your work, please cite:
                          
                          <p><a href="http://www.plantcell.org/content/30/7/1424" target="_blank">Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.</a></p>
                          
                          <h3 id="pipeline">3D analysis pipeline</h3>
                          <img style="width: 70%; display: block; margin-left: auto; margin-right: auto;" src="pipeline.png"/>
                          <h3 id="get-help">How to get help</h3>
                          <h4 id="user-manuals">User manuals</h4>
                          3D RNA-seq App "easy-to-use" manual:
                          
                          <a href="https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md" target="_blank">https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md</a>
                          
                          <h4 id="tooltips">Tooltips</h4>
                          In the GUI, users can click tooltips <i style="color: #08088A;" class="fa fa-question-circle fa-lg"></i> in specific steps for help information.
                          <h4 id="contact">Contact us</h4>
                          3D RNA-seq App is developed and maintained by Dr. Wenbin Guo from Plant Sciences Division, School of Life Sciences, University of Dundee.
                          If you have any questions and suggestions, please contact:
                          <ul>
                          <li> Dr. Wenbin Guo: wenbin.guo@hutton.ac.uk</li>
                          <li> Dr. Runxuan Zhang: runxuan.zhang@hutton.ac.uk</li>
                          </ul>
                          <h3 id="references">References</h3>
                          <ul>
                          <li>Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525-527.</li>
                          <li>Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H.G., Zhang,R., and Brown,J.W.S. (2018) Rapid and Dynamic Alternative Splicing Impacts the Arabidopsis Cold Response Transcriptome. Plant Cell, 30, 1424-1444.</li>
                          <li>Guo,W., Calixto,C.P.G., Brown,J.W.S., and Zhang,R. (2017) TSIS: An R package to infer alternative splicing isoform switches for time-series data. Bioinformatics, 33, 3308–3310.</li>
                          <li>Law,C.W., Chen,Y., Shi,W., and Smyth,G.K. (2014) voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol, 15, R29.</li>
                          <li>Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417-419.</li>
                          <li>Risso,D., Ngai,J., Speed,T.P., and Dudoit,S. (2014) Normalization of RNA-seq data using factor analysis of control genes or samples. Nat. Biotechnol., 32, 896-902.</li>
                          <li><p>Robinson,M.D., McCarthy,D.J., and Smyth,G.K. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139-40.</p></li>
                          <li>Ritchie,M.E., Phipson,B., Wu,D., Hu,Y., Law,C.W., Shi,W., and Smyth,G.K. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 43, e47.</li>
                          </ul>
                          </div>
                          ')
                 )
          )
        ),
        fluidRow(
          div(style="display:inline-block;vertical-align:middle;margin-right:15px; float: right;",
              HTML('<i>Data generation</i>'),
              actionButton(inputId = 'page_after_introduction',label = '',icon = icon('arrow-right'),
                           style="color: #fff; background-color: #428bca; border-color: #2e6da4;")
          )
        )
        )