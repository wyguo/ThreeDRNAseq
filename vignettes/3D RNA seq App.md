<style>
body {
text-align: justify}
</style>
<style>
body {
text-align: justify}
</style>
<!-- <script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script> -->
Introduction
------------

The DDD-GUI is used to perform differential expression (DE) gene and transcript, differental alternaive spliced (DAS) gene and differential transcript usage (DTU) transcript analysis. To use our pipeline in your work, please cite:

<a href="http://www.plantcell.org/content/30/7/1424" target="_blank">Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.</a>

Installation
------------

Required R packages: shiny, shinydashboard, shinyFiles, tximport, edgeR, plotly, RUVSeq, eulerr, gridExtra, grid, ComplexHeatmap, Gmisc, limma

Prepare input data
------------------

The following datasets are required in the pipeline (<a href='#Figure1'>Figure 1</a>):

-   Gene-transcript mapping information in a csv file "mapping.csv" (<a href='#Figure1'>Figure 1A</a>). In the spreadsheet, the first column is the transcript list with column name "TXNAME" and the second column is the gene list with column name "GENEID".
-   Sample information in a csv file "samples.csv" (<a href='#Figure1'>Figure 1B</a>). The spreadsheet must include columns (hilighted in read colour):
    -   brep: biological replicates are labelled as "brep1", "brep2", "brep3",...
    -   srep: if there are sequencing replicates in the dataset, they must labelled as "srep1", 'srep2', 'srep3',..., otherwise, labbed as "srep1","srep1", "srep1",...
    -   condition: conditions of the samples.
    -   path: the complete path of <a href="https://combine-lab.github.io/salmon/" target="_blank">salmon</a> quantification files "quat.sf".
-   The <a href="https://combine-lab.github.io/salmon/" target="_blank">salmon</a> \[@salmon\] quantification output files (<a href='#Figure1'>Figure 1C</a>).
-   Contrast groups for expression comparison in a csv file "contrast.csv" (<a href='#Figure1'>Figure 1D</a>). In the spreadsheet, the first column is the treatment conditions with column name "Treatment" and the second column is the control conditions with column name "Control".

Note: the analysis will read information according to the column names in the spreadsheets. Please use the same column names as in examples.

![](inst/figure/input_csv.png)

<p id="Figure1">
<strong>Figure 1:</strong> User provided datasets.
</p>
Output data
-----------

Four folders will be created in the working directory to save the analysis results:

-   data folder: the intermediate datasets in the analysis will be saved in the "data" folder in R data format (.RData). R users can open and process these datasets in R software.
-   result folder: the DE, DAS and DTU results will be saved as csv files in the "result" folder.
-   figure folder: the figures/plots will be saved in the "figure" folder.
-   report folder: the final report in html, pdf and word format generated from the GUI will be saved in the "report" folder.

DDD-GUI
-------

The GUI includes 7 pages in the sidebar menu: Introduction, Data generation, Data pre-processing, DE DAS and DTU analysis, Result summary, Advanced plot and Generate report.
<center>
<img src="inst/figure/menu.png" align="middle" />
</center>
### Data generation

1.  To begin with, users can choose a working directory to save the output data, results and figures.

![](inst/figure/wkd.png) **Figure 2**: Choose working directory
<hr>
1.  The whole pipeline has a number of steps and each step is related to several datasets. This table is the summary of datasets required for the analysis. To generate/upload the datasets, the following options are provided:

    -   If the datasets (in .RData format) do not exist in the "data" folder of working directory, users can generate these datasets step-by-step.
    -   If all the required datasets are generated and saved in the "data" folder, users can upload by clicking "Load all" button.
    -   To perform analysis or visualise results in a specific step, users also can unload the required datasets for individual steps.

![](inst/figure/load_data.png) **Figure 3**: Data information for the analysis. The columns:

-   Data -- the .RData file names in the "data" folder.
-   In data folder -- logical, whether the required datasets exist in the "data" folder.
-   In workspace -- logical, whether the required datasets have been uploaded for analysis.
-   Description -- description of the datasets.
-   Generate from -- the steps to generate the required datasets.
-   Use in -- the steps to use the datasets.

<hr>
1.  Generate gene and transcript expression. The csv files "mapping.csv" and "samples.csv" (see <a href='#Figure1'> Figure 1A</a> and <a href='#Figure1'>B</a>) are uploaded (<a href='#Figure4'>Figure 4A</a>) and visualised (<a href='#Figure4'>Figure 4C</a>) in this interface. The <a href="https://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank">tximport</a> R package \[@tximport\] is used to generate gene and transcript expression of read counts and TPMs (transcript per million reads) based on three options: "no" (default), "scaledTPM" and "lengthScaledTPM". The "lengthScaledTPM" is recommended since it scales the expression using the transcript length and then the library size over samples (see the <a href="https://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank">tximport</a> documentation for details). The generated gene and transcript expression will be save as R data object txi\_genes.RData and txi\_trans.RData in the "data" folder, respectively. Spreadsheets of read counts and TPMs will be saved as csv files in the "result" folder.

![](inst/figure/expression_generation.png)

<p id="Figure4">
<strong>Figure 4:</strong> Generate gene and transcript expression.
</p>
### Data pre-processing

Once the read counts and TPMs are generated, the data is pre-processed with steps: Merge sequencing replicates (if exist), Remove low expressed transcripts/genes, Estimate batch effects (if exist) and Data normalisation. In each step, quality control plots are generated to optimise the paramters for pre-processing (<a href='#Figure5'>Figure 5</a>).

![](inst/figure/data_preprocessing.png)

<p id="Figure5">
<strong>Figure 5:</strong> Data pre-processing. The figures were generated based on the RNA-seq data from study of Arabidopsis in response to cold \[@mainrnaseq\]
</p>
#### Merge sequencing replicates

The sequencing replicates are generated by sequencing the same tissue several times, which do not add too much information to the technical variation. Sequencing replicates from the same biological replicate are often added to increase the depth. In the panel of <a href='#Figure6'>Figure 6A</a>, the sequencing replicate information will be extracted from the "srep" column of the "samples.csv" file (<a href='#Figure1'> Figure 1B</a>) to indicate whether to merge the sequencing replicates.

#### Filter low expressed genes/transcripts

The low expressed transcripts are filtered based on count per million (CPM) reads:
$$ CPM=\\frac{X\_i}{N}\\times 10^6 $$
 where *X*<sub>*i*</sub> is read count of transcript *i* and *N* is the library size of the sample.

-   An expressed transcript must have ≥*n* samples with expression ≥*m* CPM.
-   An expressed gene must have at least one expressed transcript.

Note: The sample number cut-off *n* and CPM cut-off *m* are determined by the mean-variance trend plot (<a href='#Figure6'>Figure 6B</a> and <a href='#Figure6'>Figure 6C</a>). Read counts are assumed to be negative binomial distributed of which the mean and variance has a relation:
$$ Var(\\log\_2X)=\\frac{1}{\\mu}+\\phi $$
 where *X* is read count, *μ* is the mean and *ϕ* is the overdispersion. The expression variance decreasing monotonically with the increasing of the mean. In real RNA-seq data, the low expressed transcripts confound with noise and the expression distribution is different to the expressed transcripts. In the mean-variance trend plot, the low expressed transcripts cause a drop trend in low expression region (<a href='#Figure5'>Figure 5</a> and <a href='#Figure6'>Figure 6C</a>). Therefore, the cut-offs *n* and *m* to filter the low expression can be optimised until the drop trend disappeared in the meam-variance trend plot.

![](inst/figure/merge_filter.png)
<p id='Figure6'>
<strong>Figure 6:</strong> Merge sequencing replicates and filter low expressed genes/transcripts.
<p>
#### Principal component analysis (PCA) and batch effect estimation

PCA is a mathematical method to reduce expression data dimensionality while remaining most of the data variance. The reduction is performed by projecting the data to directions or principal components(PCs) of highest data varablity. For example. Therefore, the data main variance is accessible by investigate a few number of PCs rather than thousands of variables. In omic data analysis, the first two to four principal componets are often used to visualise the similarities and differences of samples, thereby determing the groups of samples. In the PCA scatter plot, the samples from the same condition often stay close, but far from the samples of other conditions. It can be used as evidence for data quality checking.

PCA plot also can be used to identify whether the RNA-seq data contain batch effects, which are caused by biological replications in different laboratory conditions. Compared with random noise, batch effects can be distinguished owning to the systematic biases of signals across replicates. For example, <a href='#Figure7'>Figure 7A</a> shows the PCA plot (PC1 vs PC2) of a RNA-seq data with 3 time-points T2, T10 and T19. Each time-point has 3 biological replicates. It can be seen that the first bioloigcal replicate (brep1) stay in a seperate cluster to other two replicates, which indicates a clear batch effect of the data. In this pipeline, the <a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" target="_blank">RUVSeq</a> R package \[@RUVSeq\] is used to estimate the batch effects. The RUVSeq algorithm generates batch effect terms which can be incorporated in the design matrix of linear regression model for DE and DAS analysis, i.e.

observed expression = baseline effects + batch effects + noise

It also generates a modified expression matrix in which the batch effects has been removed from the data. To avoid data over-modification, this dataset is only used to make PCA plot, but not used for downstream DE and DAS analysis (<a href='#Figure7'>Figure 7B</a>).

In this panel (<a href='#Figure7'>Figure 7</a>), users can select and visualise different PCs based on transcript level or gene level exprssion, or the average expression of biological replicates. The scatter points can be grouped and coloured according to biological replicates or conditions. Ellipses or polygons can be added to the plots to highlight the grouped clusters. Three options, RUVr, RUVs and RUVg, in the <a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" target="_blank">RUVSeq</a> R package are provided to estimate the batch effects (see the package [documentation](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html) for details). Users can "Update PCA plot" to visualise whether the batch effects are removed properly. If no clear batch effects in the PCA plot, please skip this step.

![](inst/figure/PCA.png)
<p id='Figure7'>
<strong>Figure 7:</strong> PCA plot before (A) and after (B) remove batch effects.
<p>
#### Data normalisation

To enable unbiased comparisons, read counts must be normalised on the basis of sequencing depth across samples. Normalisation methods Trimmed Mean of M-values (TMM), Relative Log Expression (RLE) and upper-quartile are provided to reduce the effects from the systematic technical biases across samples \[@TMM; @Norm\]. Boxplots are used to visualise the expression distribution across samples before and after normalisation.

![](inst/figure/normalisation.png) Figure 8: Data normalisation.

DE, DAS and DTU analysis
------------------------

### Definition of DE, DAS and DTU in this pipeline

![](inst/figure/DDD_definition.png)

**Statistics:**

-   Adjusted p-value: in the expression comparitive analysis, each gene/transcript is fitted with a statistical model, reporting a p-value. However the testing is performed over a substantial population of genes/transcripts and each testing has a probability of making errors. Hence, all p-values in the analsysis are adjusted by controlling the false discovery rate (FDR). Please see the R function <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html',target='_blank'>p.adjust</a> for details.
-   *L*<sub>2</sub>*F**C*: log<sub>2</sub> fold change of expression based on contrast groups.
-   *Δ**P**S*: Transcript level PS (percent of splice) is defined as the ratio of transcript average abundance of conditions divided by the average gene abundance. *Δ**P**S* is the PS differences of conditions based on the contrast groups. The abundance to calculate PS values can be TPMs or read counts.

**DE gene/transcript analysis:** A gene/transcript is identified as DE in a contrast group if *L*<sub>2</sub>*F**C* of expression ≥ cut-off and with adjusted p-value &lt; cut-off.

Two pipelines are provided for DE analysis: <a target="_blank" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29">"limma-voom"</a> and <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/">"edgeR"</a>. The <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/">"edgeR"</a> pipeline includes two methods: <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmQLFit">"glmQL"</a> (Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests) and <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmFit">"glm"</a> (Genewise Negative Binomial Generalized Linear Models). In limma pipeline, *L*<sub>2</sub>*F**C**s* are calculated based count per million reads (CPMs) while in edgeR, *L*<sub>2</sub>*F**C**s* are based on read counts. From several real RNA-seq data analyses, high consensus is achieved between <a target="_blank" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29">"limma-voom"</a> and <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmQLFit">"glmQL"</a> (&gt;90%) and they have more stringent controls of false discovery rate than the <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmFit">"glm"</a> method.

**DAS gene and DTU transcript analysis:** To test differential expression, the gene level expression is compared to transcript level expression in the contrast groups. A gene is DAS in a contrast group if adjusted p-value &lt; cut-off and at least one transcript of the gene with *Δ**P**S*≥ cut-off. A transcript is DTU if adjusted p-value &lt; cut-off and *Δ**P**S*≥ cut-off.

**Two pipelines for expression comparative analysis:**

-   *limma pipeline*: To identify DAS genes, the expression of each transcript is compared to the weighted average expression of all the transcripts for the same gene (weight on transcript expression variance; the average expression can be treated as gene level expression). P-value of each test is converted to genewise p-value by using F-test or Simes method. To identify DTU transcript, each transcript is compared to the weighted average of all the other transcripts for the same gene. See the <a href='http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/diffSplice.html' target='_blank'>diffSplice</a> function in <a href='https://bioconductor.org/packages/release/bioc/html/limma.html' target='_blank'>limma</a> package for details.
-   *edgeR pipeline*: To identify DTU transcripts, log-fold-change of each transcript is compared to log-fold-change of entire gene (sum of all transcripts). The DTU transcript test p-values are summarised to genewise p-value by using F-test or Simes method to report DAS genes. See <a href='https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/diffSpliceDGE' target='_blank'>diffSpliceDGE</a> function in <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/">"edgeR"</a> R package for details.

In the GUI, contrast groups in csv file (<a href='#Figure1'>Figure 1D</a>) can be uploaded in panel <a href='#Figure9'>Figure 9A</a>. DE, DAS and DTU analysis pipeline, *Δ**P**S* and cut-offs of test statistics can be selected/generated in panels <a href='#Figure9'>Figure 9B</a> and <a href='#Figure9'>C</a>. Subset of *Δ**P**S* can be visualised in panel <a href='#Figure9'>Figure 9D</a>.

![](inst/figure/ddd.png)
<p id='Figure9'>
<strong>Figure 9:</strong> DE and DAS analysis.
<p>
Result summary
--------------

This page includes panels to show the DE, DAS and DTU analysis results (<a href='#Figure10'>Figure 10</a>-<a href='#Figure10'>12</a>):

-   The test statistics in different contrast groups, e.g. adjusted p-value and *L*<sub>2</sub>*F**C*.
-   The number of genes and transcripts with significant expression changes in contrast groups. The numbers of DE only genes, DAS only genes, bot DE and DAS genes, DE only transcripts, DTU only transcripts and both DE and DTU transcripts. The numbers and their comparisons are present as Euler diagrams and flow charts.
-   Up- and down-regulation gene and transcript numbers are shown as bar plots.
-   All the tables and plots can be saved to local folder in the working directory.

![](inst/figure/ddd_summary.png)
<p id='Figure10'>
<strong>Figure 10:</strong> Number of DE genes, DAS genes, DE transcripts and DTU transcripts, and bar plots of up-down regulation numbers.
<p>
![](inst/figure/venn.png)
<p id='Figure11'>
<strong>Figure 11:</strong> Euler diagram of numbers between contrast groups (A) and numbers to compare DE vs DAS genes and DT vs DTU transcripts (B).
<p>
![](inst/figure/flowchat.png)
<p id='Figure12'>
<strong>Figure 12:</strong> Flowchart to show numbers of DE and/or DAS genes and DE and/or DTU transcripts.
<p>
<!-- ## Session information -->
<!-- sessionInfo() -->
Advanced plot
-------------

### Heatmap

Users can make heatmap for DE genes, DAS genes, DE transcripts and DTU transcripts identified in the analysis (<a href='#Figure13'>Figure 13A</a>). Users can also upload a gene or transcript list in csv file to make heatmap for specific list (<a href='#Figure13'>Figure 13B</a>). To make the heatmap, the average TPMs of conditions for the targets are standardized into z-scores. The <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/dist.html' target='_blank'>dist</a> and <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html' target='_blank'>hclust</a> R functions are used to perform the hierarchical clustering analysis. The heamaps and the target list in each cluster of the heatmaps can be saved to local folder.

![](inst/figure/heatmap.png)
<p id='Figure13'>
<strong>Figure 13:</strong> Heatmap panel (A) and example user provided gene list in csv file for heatmap (B).
<p>
### Profile plot

Gene and transcript expression profiles (TPMs or read counts) and PS (percent spliced) can be visualised by typing a gene name (<a href='#Figure14'>Figure 14A</a>). Users can also generate a number of plots by provides a gene list (<a href='#Figure13'>Figure 13B</a>) and all the plots will be saved to "figure" folder in the working directory (<a href='#Figure14'>Figure 14B</a>).

![](inst/figure/profiles.png)
<p id='Figure14'>
<strong>Figure 14:</strong> Plot of expression profiles and PS across conditions
<p>
### GO annotation plot

Users can generate DE and DAS gene list by click "Generate" button (<a href='#Figure15'>Figure 15B</a>). These gene lists can be uploaded to Gene Ontology (GO) analysis tools (e.g. <a href='https://david.ncifcrf.gov/' target='_blank'>DAVID</a> and <a href='http://bioinfo.cau.edu.cn/agriGO/' target='_blank'>agriGO</a>) to generate GO annotation. To visualise the annotation plot, a csv file need to be uploaded to the GUI (<a href='#Figure15'>Figure 15B</a>). The file includes a "Category" column of CC (cellular component), BP (biological process) and MF (molecular function), a column of "Term" of GO annotation and columns of statistics to report the annotation enrichment, e.g. count, FDR, -log10(FDR), etc. (<a href='#Figure15'>Figure 15A</a>). In the panel, users can select different statistics to visulise. The selected gene list type will present in the plot title to distinguish whether the provided gene list is DE or DAS genes.

![](inst/figure/go.png)
<p id='Figure15'>
<strong>Figure 15:</strong> GO annotation plot
<p>
Generate report
---------------

The table in this page lists the location of save files and the parameters for the analysis. If the information is wroing, users can double click the cells in the table to correct to values. Based on this information, report in three format, word, pdf and html, will be generated and saved in the "report" folder of working directory.

![](inst/figure/report.png)
<p id='Figure16'>
<strong>Figure 16:</strong> Generate report
<p>
Session information
-------------------

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 7 x64 (build 7601) Service Pack 1
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] Gmisc_1.6.4                 htmlTable_1.12             
    ##  [3] Rcpp_0.12.18                ComplexHeatmap_1.18.1      
    ##  [5] gridExtra_2.3               eulerr_4.1.0               
    ##  [7] RUVSeq_1.14.0               EDASeq_2.14.1              
    ##  [9] ShortRead_1.38.0            GenomicAlignments_1.16.0   
    ## [11] SummarizedExperiment_1.10.1 DelayedArray_0.6.6         
    ## [13] matrixStats_0.54.0          Rsamtools_1.32.3           
    ## [15] GenomicRanges_1.32.6        GenomeInfoDb_1.16.0        
    ## [17] Biostrings_2.48.0           XVector_0.20.0             
    ## [19] IRanges_2.14.11             S4Vectors_0.18.3           
    ## [21] BiocParallel_1.14.2         Biobase_2.40.0             
    ## [23] BiocGenerics_0.26.0         plotly_4.8.0.9000          
    ## [25] ggplot2_3.0.0               edgeR_3.22.3               
    ## [27] limma_3.36.3                tximport_1.11.0            
    ## [29] shinyFiles_0.7.1            shinydashboard_0.7.0       
    ## [31] shiny_1.1.0                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_1.3-2       rjson_0.2.20           hwriter_1.3.2         
    ##  [4] rprojroot_1.3-2        circlize_0.4.4         base64enc_0.1-3       
    ##  [7] GlobalOptions_0.1.0    fs_1.2.6               rstudioapi_0.7        
    ## [10] bit64_0.9-7            AnnotationDbi_1.42.1   splines_3.5.1         
    ## [13] R.methodsS3_1.7.1      DESeq_1.32.0           geneplotter_1.58.0    
    ## [16] knitr_1.20             Formula_1.2-3          jsonlite_1.5          
    ## [19] annotate_1.58.0        cluster_2.0.7-1        R.oo_1.22.0           
    ## [22] compiler_3.5.1         httr_1.3.1             backports_1.1.2       
    ## [25] assertthat_0.2.0       Matrix_1.2-14          lazyeval_0.2.1        
    ## [28] later_0.7.4            acepack_1.4.1          htmltools_0.3.6       
    ## [31] prettyunits_1.0.2      tools_3.5.1            bindrcpp_0.2.2        
    ## [34] gtable_0.2.0           glue_1.3.0             GenomeInfoDbData_1.1.0
    ## [37] dplyr_0.7.6            rtracklayer_1.40.6     stringr_1.3.1         
    ## [40] mime_0.5               XML_3.98-1.16          zlibbioc_1.26.0       
    ## [43] MASS_7.3-50            scales_1.0.0           aroma.light_3.10.0    
    ## [46] hms_0.4.2              promises_1.0.1         RColorBrewer_1.1-2    
    ## [49] yaml_2.2.0             memoise_1.1.0          rpart_4.1-13          
    ## [52] biomaRt_2.36.1         latticeExtra_0.6-28    stringi_1.1.7         
    ## [55] RSQLite_2.1.1          genefilter_1.62.0      checkmate_1.8.5       
    ## [58] GenomicFeatures_1.32.2 shape_1.4.4            rlang_0.2.2           
    ## [61] pkgconfig_2.0.2        bitops_1.0-6           evaluate_0.11         
    ## [64] lattice_0.20-35        purrr_0.2.5            bindr_0.1.1           
    ## [67] htmlwidgets_1.2        bit_1.1-14             tidyselect_0.2.4      
    ## [70] plyr_1.8.4             magrittr_1.5           R6_2.2.2              
    ## [73] Hmisc_4.1-1            DBI_1.0.0              foreign_0.8-71        
    ## [76] pillar_1.3.0           withr_2.1.2            nnet_7.3-12           
    ## [79] forestplot_1.7.2       abind_1.4-5            survival_2.42-6       
    ## [82] RCurl_1.95-4.11        tibble_1.4.2           crayon_1.3.4          
    ## [85] rmarkdown_1.10         GetoptLong_0.1.7       progress_1.2.0        
    ## [88] locfit_1.5-9.1         data.table_1.11.6      blob_1.1.1            
    ## [91] digest_0.6.17          xtable_1.8-3           tidyr_0.8.1           
    ## [94] httpuv_1.4.5           R.utils_2.7.0          munsell_0.5.0         
    ## [97] viridisLite_0.3.0

References
----------
