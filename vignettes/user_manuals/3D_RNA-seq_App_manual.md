3D RNA-seq App
==============
<i><h4>Interactive shiny App user manual</h4></i>
<i><h4>Wenbin Guo</h4></i>
<i><h4>26 November 2018</h4></i>

Table of contents
-----------------

-   [Table of contents](#table-of-contents)
-   [How to use this manual](#how-to-use-this-manual)
-   [Installation](#installation)
    -   [Install ThreeDRNAseq package](#install-threedrnaseq-package)
    -   [Install dependency packages](#install-dependency-packages)
-   [Example data](#example-data)
-   [Prepare input data](#prepare-input-data)
-   [Output data](#output-data)
-   [3D RNA-seq App](#d-rna-seq-app)
    -   [Data generation](#data-generation)
    -   [Data pre-processing](#data-pre-processing)
    -   [3D RNA-seq analysis](#d-rna-seq-analysis)
    -   [Result summary](#result-summary)
    -   [Advanced plot](#advanced-plot)
    -   [Generate report](#generate-report)
-   [Saved files in local directory](#saved-files-in-local-directory)
    -   [Saved files in the "data" folder](#saved-files-in-the-data-folder)
    -   [Saved files in the "result" folder](#saved-files-in-the-result-folder)
    -   [Saved files in the "figure" folder](#saved-files-in-the-figure-folder)
    -   [Saved files in the "report" folder](#saved-files-in-the-report-folder)
-   [References](#references)
-   [Session information](#session-information)

<div align="justify">

How to use this manual
----------------------

The 3D RNA-seq pipeline includes complex steps of expression data generation, data pre-processing, 3D analysis and result summary. This manual provides step-by-step instructions to use the 3D RNA-seq App. To use this App in your work, please cite:

<a href="http://www.plantcell.org/content/30/7/1424" target="_blank">Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.</a>

Command-line based 3D analysis manual: <a href="https://github.com/wyguo/ThreeDRNAseq/tree/master/vignettes/user_manuals" target="_blank">
<button class="button">
<font color="black">Click me!</font>
</button>
</a>

Installation and Run
--------------------

ThreeDRNAseq R package can be installed from Github by using <a href='https://cran.r-project.org/web/packages/devtools/index.html' target='_blank'>devtools</a> R package

### Install ThreeDRNAseq package

``` r
##################################################################################################
## use devtools R package to install ThreeDRNAseq from Github
###---> If devtools is not installed, please install
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages('devtools')

###---> Install ThreeDRNAseq
devtools::install_github('wyguo/ThreeDRNAseq')
```

### Install dependency packages

``` r
#######################################################################################################
## Install packages of dependency
###---> Install packages from Cran
cran.package.list <- c('shiny','shinydashboard','rhandsontable','shinyFiles','shinyjs','DT',
                       'plotly','ggplot2','eulerr','ggrepel',
                       'gridExtra','fastcluster','Gmisc','rmarkdown')
for(i in cran.package.list){
   if(!(i %in% rownames(installed.packages()))){
     message('Installing package: ',i)
     install.packages(i,dependencies = T)
   } else next
}

###---> Install packages from Bioconductor
bioconductor.package.list <- c('tximport','edgeR','limma','RUVSeq','ComplexHeatmap','rhdf5')
for(i in bioconductor.package.list){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!(i %in% rownames(installed.packages()))){
    message('Installing package: ',i)
    BiocManager::install(i, version = "3.8")
  } else next
}
```
### Run 3D RNA-seq App
Please type the command in your RStudio and the App will be launched in default web browser (Google Chrome is recommended for image visualisation issues).
```r
library(ThreeDRNAseq)
ThreeDRNAseq.app()
```

<table>
<caption>The main R packages used in the 3D RNA-seq analysis. If any other R packages are missing from your PC, please install accordingly.</caption>
<colgroup>
<col width="21%" />
<col width="14%" />
<col width="64%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Usage</th>
<th align="left">Package</th>
<th align="left">Link</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Make the shiny app</td>
<td align="left">shiny</td>
<td align="left"><a href="https://shiny.rstudio.com/" class="uri">https://shiny.rstudio.com/</a></td>
</tr>
<tr class="even">
<td align="left"></td>
<td align="left">shinydashboard</td>
<td align="left"><a href="https://rstudio.github.io/shinydashboard/" class="uri">https://rstudio.github.io/shinydashboard/</a></td>
</tr>
<tr class="odd">
<td align="left"></td>
<td align="left">rhandsontable</td>
<td align="left"><a href="https://jrowen.github.io/rhandsontable/" class="uri">https://jrowen.github.io/rhandsontable/</a></td>
</tr>
<tr class="even">
<td align="left"></td>
<td align="left">shinyFiles</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/shinyFiles/index.html" class="uri">https://cran.r-project.org/web/packages/shinyFiles/index.html</a></td>
</tr>
<tr class="odd">
<td align="left"></td>
<td align="left">shinyjs</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/shinyjs/index.html" class="uri">https://cran.r-project.org/web/packages/shinyjs/index.html</a></td>
</tr>
<tr class="even">
<td align="left"></td>
<td align="left">DT</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/DT/index.html" class="uri">https://cran.r-project.org/web/packages/DT/index.html</a></td>
</tr>
<tr class="odd">
<td align="left">Read csv file</td>
<td align="left">readr</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/readr/index.html" class="uri">https://cran.r-project.org/web/packages/readr/index.html</a></td>
</tr>
<tr class="even">
<td align="left">Expression generator</td>
<td align="left">tximport</td>
<td align="left"><a href="http://bioconductor.org/packages/release/bioc/html/tximport.html" class="uri">http://bioconductor.org/packages/release/bioc/html/tximport.html</a></td>
</tr>
<tr class="odd">
<td align="left">Batch effect estimation</td>
<td align="left">RUVSeq</td>
<td align="left"><a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" class="uri">https://bioconductor.org/packages/release/bioc/html/RUVSeq.html</a></td>
</tr>
<tr class="even">
<td align="left">3D analysis</td>
<td align="left">limma</td>
<td align="left"><a href="https://bioconductor.org/packages/release/bioc/html/limma.html" class="uri">https://bioconductor.org/packages/release/bioc/html/limma.html</a></td>
</tr>
<tr class="odd">
<td align="left"></td>
<td align="left">edgeR</td>
<td align="left"><a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" class="uri">https://bioconductor.org/packages/release/bioc/html/edgeR.html</a></td>
</tr>
<tr class="even">
<td align="left">Heatmap</td>
<td align="left">ComplexHeatmap</td>
<td align="left"><a href="https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html" class="uri">https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html</a></td>
</tr>
<tr class="odd">
<td align="left">Cluster of heatmap</td>
<td align="left">fastcluster</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/fastcluster/index.html" class="uri">https://cran.r-project.org/web/packages/fastcluster/index.html</a></td>
</tr>
<tr class="even">
<td align="left">Flow charts</td>
<td align="left">Gmisc</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/Gmisc/index.html" class="uri">https://cran.r-project.org/web/packages/Gmisc/index.html</a></td>
</tr>
<tr class="odd">
<td align="left">Euler plots</td>
<td align="left">eulerr</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/eulerr/index.html" class="uri">https://cran.r-project.org/web/packages/eulerr/index.html</a></td>
</tr>
<tr class="even">
<td align="left">Other plots</td>
<td align="left">ggplot2</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/ggplot2/index.html" class="uri">https://cran.r-project.org/web/packages/ggplot2/index.html</a></td>
</tr>
<tr class="odd">
<td align="left"></td>
<td align="left">plotly</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/plotly/index.html" class="uri">https://cran.r-project.org/web/packages/plotly/index.html</a></td>
</tr>
<tr class="even">
<td align="left"></td>
<td align="left">gridExtra</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/gridExtra/index.html" class="uri">https://cran.r-project.org/web/packages/gridExtra/index.html</a></td>
</tr>
<tr class="odd">
<td align="left"></td>
<td align="left">grid</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/grid/index.html" class="uri">https://cran.r-project.org/web/packages/grid/index.html</a></td>
</tr>
<tr class="even">
<td align="left">Generate report</td>
<td align="left">rmarkdown</td>
<td align="left"><a href="https://cran.r-project.org/web/packages/rmarkdown/index.html" class="uri">https://cran.r-project.org/web/packages/rmarkdown/index.html</a></td>
</tr>
</tbody>
</table>

[Go back to Table of contents](#table-of-contents)

Example data
------------

Example data can be downloaded from: xxx

The example RNA-seq data is a subset from study of Arabidopsis in response to cold (Calixto et al., 2018). Datasets of three time-points (red circles) were extracted from the whole study:

<center>
<img src="https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/fig/example_exp.png?raw=true">
</center>
<br>

| Time-points | Description                                                                    |
|:------------|:-------------------------------------------------------------------------------|
| T2          | Second time-point of Day 0, temperature 20<sup>*o*</sup>*C*                    |
| T10         | Second time-point of Day 1, temperature 4<sup>*o*</sup>*C* of cold transition  |
| T19         | Second time-point of Day 4, temperature 4<sup>*o*</sup>*C* of cold acclimation |

Each time-point has 3 biological replicates and each biological replicate has 3 sequencing replicates. Transcript quantification was generated by using Salmon (Patro et al., 2017) and AtRTD2-QUASI (Zhang et al., 2017) as transcriptome reference. Expression comparisons were designed as "T10 vs T2" and "T19-T2" to identify the cold response genes and transcripts at both transcription and AS level.

Data download links:
- Google drive: <a href="https://drive.google.com/open?id=1VFxWwi6q3Lv1o79RlKDtk1q-BJXa81wD" target="_blank">https://drive.google.com/open?id=1VFxWwi6q3Lv1o79RlKDtk1q-BJXa81wD</a>
<br>

[Go back to Table of contents](#table-of-contents)

Prepare input data
------------------

The following datasets are required in the pipeline (<a href='#Figure1'>Figure 1</a>):

-   Gene-transcript mapping information in a csv file **"mapping.csv"** (<a href='#Figure1'>Figure 1A</a>). In the spreadsheet, the first column is the transcript list with column name "TXNAME" and the second column is the gene list with column name "GENEID".
-   Sample information in a csv file **"samples.csv"** (<a href='#Figure1'>Figure 1B</a>). The spreadsheet must include columns (highlighted in read colour):
    -   brep: biological replicates are labelled as "brep1", "brep2", "brep3",... to distinguish the samples.
    -   srep: if there are sequencing replicates (the same biological replicate is sequenced multiple times) in the dataset, they must labelled as "srep1", 'srep2', 'srep3',..., otherwise, labelled as "srep1","srep1", "srep1",... for all the samples.
    -   condition: conditions of the samples.
    -   path: the complete path of transcript quantification, e.g. the "quant.sf" files of <a href="https://combine-lab.github.io/salmon/" target="_blank">Salmon</a> (Patro et al., 2017) and "abundance.h5" files of <a href="https://pachterlab.github.io/kallisto/about" target="_blank">Kallisto</a> (Bray et al., 2016)
-   The **transcript quantification outputs**. The folder names of the quantification outputs must match to the path information in the "samples.csv" files. For example, <a href='#Figure1'>Figure 1C</a> shows the <a href="https://combine-lab.github.io/salmon/" target="_blank">Salmon</a> (Patro et al., 2017) quantification output files.
-   Contrast groups for expression comparison in a csv file **"contrast.csv"** (<a href='#Figure1'>Figure 1D</a>). In the spreadsheet, the first column is the treatment conditions with column name "Treatment" and the second column is the control conditions with column name "Control".

**Note:** the analysis will read information according to the column names in the spreadsheets. Please use the same column names as in examples.

![](3D_App_figure/input_csv.png)

<p id="Figure1">
<strong>Figure 1:</strong> User provided datasets.
</p>
<p id="Section_output_data">
</p>
<br>

[Go back to Table of contents](#table-of-contents)

Output data
-----------

Four folders will be created in the working directory to save the analysis results:

-   data folder: the intermediate datasets in the analysis will be saved in the "data" folder in R data format (.RData). Advanced R users can open and process these objects in R software. For basic R users, please keep in mind these ".RData" objects are used to store the intermediate datasets.
-   result folder: the DE, DAS and DTU results will be saved as csv files in the "result" folder.
-   figure folder: the figures/plots will be saved in the "figure" folder.
-   report folder: the final report in html, pdf and word formats generated from the App will be saved in the "report" folder.

The detailed descriptions of files saved in the above folders can be found in Section <a href='#Section_saved_files'>"Saved files in local directory"</a>.

[Go back to Table of contents](#table-of-contents)

3D RNA-seq App
--------------

The App includes 7 pages in the sidebar menu: Introduction, Data generation, Data pre-processing, DE DAS and DTU analysis, Result summary, Advanced plot and Generate report.
<center>
<img src="3D_App_figure/menu.png" align="middle" />
</center>
<br>

### Data generation

#### Set directory and load datasets

To begin with, users can choose a working directory to save the output data, results and figures.

![](3D_App_figure/wkd.png) **Figure 2**: Choose working directory. Please click the "triangles" in the left-side "Directories" box to access and select sub-folders. The click and select in the right-side "Content" box is invalid.

The whole pipeline has a number of steps and each step is related to several intermediate datasets (see Section "Output data"). This table summarise the generation and usage steps of these datasets and whether (True or False) they have been loaded in workspace ready for use. To avoid running the pipeline from very beginning every time the App get refreshed, the following options are provided:

-   If all/part of the required datasets (in .RData format) do not exist in the "data" folder of working directory (see Section: <a href='#Section_saved_files'>"Saved files in local directory"</a>), users can generate these datasets from corresponding steps.
-   If all/part of the required datasets (in .RData format) already exist in the "data" folder, users can click "Load all" button to load them.
-   To perform analysis or visualise results in a specific step, users also can load the required datasets for individual steps.

![](3D_App_figure/load_data.png) **Figure 3**: Data information for the analysis. The columns:

-   Data -- the .RData file names in the "data" folder.
-   In data folder -- logical, whether the required datasets exist in the "data" folder.
-   In workspace -- logical, whether the required datasets have been uploaded for analysis.
-   Description -- description of the datasets.
-   Generate from -- the steps to generate the required datasets.
-   Use in -- the steps to use the datasets.

#### Generate gene and transcript expression

The csv files "mapping.csv" and "samples.csv" (see <a href='#Figure1'> Figure 1A</a> and <a href='#Figure1'>B</a>) are loaded (<a href='#Figure4'>Figure 4A</a>) and visualised (<a href='#Figure4'>Figure 4C</a>) in this step. The <a href="https://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank">tximport</a> R package (Soneson et al., 2016) is used to generate gene and transcript expression of read counts and TPMs (transcript per million reads) based on three options: "no" (default), "scaledTPM" and "lengthScaledTPM". The "lengthScaledTPM" is recommended since it scales the expression using the transcript length and then the library size over samples (see the <a href="https://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank">tximport</a> documentation for details). The generated gene and transcript expression will be save as R data object txi\_genes.RData and txi\_trans.RData in the "data" folder, respectively. Spreadsheets of read counts and TPMs will be saved as csv files in the "result" folder.

![](3D_App_figure/expression_generation.png)

<p id="Figure4">
<strong>Figure 4:</strong> Generate gene and transcript expression.
</p>
<br>

[Go back to Table of contents](#table-of-contents)

### Data pre-processing

Once the read counts and TPMs are generated, the data is pre-processed with steps: Merge sequencing replicates (if exist), Remove low expressed transcripts/genes, Estimate batch effects (if exist) and Data normalisation. In each step, quality control plots are generated to optimise the parameters for pre-processing (<a href='#Figure5'>Figure 5</a>).

![](3D_App_figure/data_preprocessing.png)

<p id="Figure5">
<strong>Figure 5:</strong> Data pre-processing. The figures were generated based on the RNA-seq data from study of Arabidopsis in response to cold (Calixto et al., 2018).
</p>
<br>

#### Merge sequencing replicates

The sequencing replicates are generated by sequencing the same biological replicate multiple times, which do not add too much information to the technical variation. Sequencing replicates from the same biological replicate are often added to increase the depth. In the panel of <a href='#Figure6'>Figure 6A</a>, the sequencing replicate information will be extracted from the "srep" column of the "samples.csv" file (<a href='#Figure1'> Figure 1B</a>) to indicate whether to merge the sequencing replicates.

#### Filter low expression

The low expressed transcripts are filtered based on count per million (CPM) reads:
$$ CPM=\\frac{X\_i}{N}\\times 10^6 $$
 where *X*<sub>*i*</sub> is read count of transcript *i* and *N* is the library size of the sample.

-   An expressed transcript must have ≥*n* samples with expression ≥*m* CPM.
-   An expressed gene must have at least one expressed transcript.

The sample number cut-off *n* and CPM cut-off *m* are determined by the mean-variance trend plot (<a href='#Figure6'>Figure 6B</a> and <a href='#Figure6'>Figure 6C</a>). Read counts are assumed to be negative binomial distributed of which the mean and variance has a relation:
$$ Var(\\log\_2X)=\\frac{1}{\\mu}+\\phi $$
 where *X* is read count, *μ* is the mean and *ϕ* is the overdispersion. The expression variance decreasing monotonically with the increasing of the mean. In real RNA-seq data, the low expressed transcripts confound with noise and the expression distribution is different to the expressed transcripts. In the mean-variance trend plot, the low expressed transcripts cause a drop trend in low expression region (<a href='#Figure5'>Figure 5</a> and <a href='#Figure6'>Figure 6C</a>). Therefore, the cut-offs *n* and *m* to filter the low expression can be optimised until the drop trend in low expression region disappeared in the mean-variance trend plot.

![](3D_App_figure/merge_filter.png)
<p id='Figure6'>
<strong>Figure 6:</strong> Merge sequencing replicates and filter low expressed genes/transcripts. (D) is the scatter plot of mean-variance trends before and after filtering the low expressed transcripts. The fitted lines of the trends are also presented (Yellow line: before filtering; Red line: after filtering).
<p>
<br>

#### Principal component analysis (PCA) and batch effect estimation

PCA is a mathematical method to reduce expression data dimensionality while remaining most of the data variance. The reduction is performed by projecting the data to directions or principal components (PCs) from highest to lowest data variability. Therefore, the data main variance is accessible by investigate top few number of PCs rather than thousands of variables. In Omic data analysis, the first two to four principal components are often used to visualise the similarities and differences of samples, thereby determining the grouping of samples. In the PCA scatter plot, the samples from the same condition often stay close, but far from the samples of other conditions. It can be used as evidence for data quality checking.

PCA plot also can be used to identify whether the RNA-seq data contain batch effects, which are caused by biological replications in different laboratory conditions. Compared with random noise, batch effects can be distinguished owning to the systematic biases of signals across replicates. For example, <a href='#Figure7'>Figure 7A</a> shows the PCA plot (PC1 vs PC2) of a RNA-seq data with 3 time-points T2, T10 and T19. Each time-point has 3 biological replicates. It can be seen that the first biological replicate (brep1) stays in a separate cluster to other two replicates, which indicates a clear batch effect of the data. In this pipeline, the <a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" target="_blank">RUVSeq</a> R package (Risso et al., 2014) is used to estimate the batch effects. The RUVSeq algorithm generates batch effect terms which can be incorporated in the design matrix of linear regression model for 3D analysis, i.e.

observed expression = baseline effects + batch effects + noise

It also generates a pseudo read count matrix in which the batch effects have been removed from the data. To avoid data over-modification, this matrix is only used to make PCA plot, but not used for downstream 3D analysis (<a href='#Figure7'>Figure 7B</a>).

In this panel (<a href='#Figure7'>Figure 7</a>), users can select and visualise different PCs based on transcript level or gene level expression, or the average expression of biological replicates. The scatter points can be grouped and coloured according to biological replicates or conditions. Ellipses or polygons can be added to the plots to highlight the grouped clusters. Three options, RUVr, RUVs and RUVg, in the <a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" target="_blank">RUVSeq</a> R package are provided to estimate the batch effects (see the package [documentation](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html) for details). Users can "Update PCA plot" to visualise whether the batch effects are removed properly. **If there are no clear batch effects in the PCA plot, please skip this step.**

![](3D_App_figure/PCA.png)
<p id='Figure7'>
<strong>Figure 7:</strong> PCA plot before (A) and after (B) remove batch effects.
<p>
<br>

#### Data normalisation

To enable unbiased comparisons across samples, read counts must be normalised on the basis of sequencing depths. Normalisation methods Trimmed Mean of M-values (TMM), Relative Log Expression (RLE) and upper-quartile can be used to reduce the effects from the systematic technical biases across samples (Bullard et al., 2010). Box plots are used to visualise the expression distribution of raw read counts and normalised expression across samples.

![](3D_App_figure/normalisation.png) Figure 8: Data normalisation. <br>

[Go back to Table of contents](#table-of-contents)

### 3D RNA-seq analysis

#### Definition of 3D in this pipeline

![](3D_App_figure/DDD_definition.png) <br>

#### Statistics

-   Contrast groups: refer to the conditions to compare, e.g. "B-A" and "C-A" compares the conditions "B" and "C" to condition "A".
-   Adjusted p-value: in the expression comparative analysis, each gene/transcript is fitted with a statistical model, reporting a p-value. However, the testing is performed over a substantial number of genes/transcripts and each testing has a probability of making errors. Hence, all p-values in the analysis are adjusted to control the false discovery rate (FDR). Please see the R function <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html',target='_blank'>p.adjust</a> for details.
-   *L*<sub>2</sub>*F**C*: log<sub>2</sub> fold change of expression based on contrast groups.
-   *Δ**P**S*: Transcript level PS (percent of splice) is defined as the ratio of transcript average abundance of conditions divided by the average gene abundance. *Δ**P**S* is the PS differences of conditions based on the contrast groups. We recommend to use TPMs to calculate PS values. The option to generate PS from read counts is also provided (<a href='#Figure9'>Figure 9C</a>).

#### DE gene/transcript analysis

A gene/transcript is identified as DE in a contrast group if *L*<sub>2</sub>*F**C* of expression ≥ cut-off and with adjusted p-value &lt; cut-off.

Two pipelines are provided for DE analysis: <a target="_blank" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29">"limma-voom"</a> and <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/">"edgeR"</a>. The <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/">"edgeR"</a> pipeline includes two methods: <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmQLFit">"glmQL"</a> (Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests) and <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmFit">"glm"</a> (Genewise Negative Binomial Generalized Linear Models). In limma pipeline, *L*<sub>2</sub>*F**C**s* are calculated based count per million reads (CPMs) while in edgeR, *L*<sub>2</sub>*F**C**s* are based on read counts. From several real RNA-seq data analyses, high consensus is achieved between <a target="_blank" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29">"limma-voom"</a> and <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmQLFit">"glmQL"</a> (&gt;90%) and they have more stringent controls of false discovery rate than the <a target="_blank" href="https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/glmFit">"glm"</a> method.

#### DAS gene and DTU transcript analysis

To test differential alternative splicing (AS), the gene level expression is compared to transcript level expression in the contrast groups. A gene is DAS in a contrast group if adjusted p-value &lt; cut-off and at least one transcript of the gene with *Δ**P**S*≥ cut-off. A transcript is DTU if adjusted p-value &lt; cut-off and *Δ**P**S*≥ cut-off.

**Two pipelines for AS analysis:**

-   *limma pipeline*: To identify DAS genes, the expression of each transcript is compared to the weighted average expression of all the transcripts for the same gene (weight on transcript expression variance; the average expression can be treated as gene level expression). P-value of each test is converted to genewise p-value by using F-test or Simes method. To identify DTU transcript, each transcript is compared to the weighted average of all the other transcripts for the same gene. See the <a href='http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/diffSplice.html' target='_blank'>diffSplice</a> function in <a href='https://bioconductor.org/packages/release/bioc/html/limma.html' target='_blank'>limma</a> package (Ritchie et al., 2015) for details.
-   *edgeR pipeline*: To identify DTU transcripts, log-fold-change of each transcript is compared to log-fold-change of entire gene (sum of all transcripts). The DTU transcript test p-values are summarised to genewise p-value by using F-test or Simes method to report DAS genes. See <a href='https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/diffSpliceDGE' target='_blank'>diffSpliceDGE</a> function in <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/">"edgeR"</a> R package (Robinson et al., 2011) for details.

To perform the analysis, the contrast groups in csv file (<a href='#Figure1'>Figure 1D</a>) can be loaded in panel <a href='#Figure9'>Figure 9A</a>. 3D analysis pipeline, *Δ**P**S* and cut-offs of test statistics can be selected/generated in panels <a href='#Figure9'>Figure 9B</a> and <a href='#Figure9'>C</a>. Subset of *Δ**P**S* can be visualised in panel <a href='#Figure9'>Figure 9D</a>.

![](3D_App_figure/ddd.png)
<p id='Figure9'>
<strong>Figure 9:</strong> 3D RNA-seq analysis panels.
<p>
<br>

[Go back to Table of contents](#table-of-contents)

### Result summary

This page includes panels to show the 3D analysis results (<a href='#Figure10'>Figure 10</a>-<a href='#Figure10'>12</a>):

-   The test statistics in different contrast groups, e.g. adjusted p-value and *L*<sub>2</sub>*F**C*.
-   The number of genes and transcripts with significant expression changes in contrast groups. The numbers of DE only genes, DAS only genes, both DE and DAS genes, DE only transcripts, DTU only transcripts and both DE and DTU transcripts. The numbers and their comparisons are present as Euler diagrams and flow charts.
-   Up- and down-regulation gene and transcript numbers are shown as bar plots.
-   All the tables and plots can be saved to local folder in the working directory (see Section: <a href='#Section_saved_files'>"Saved files in local directory"</a>).

![](3D_App_figure/ddd_summary.png)
<p id='Figure10'>
<strong>Figure 10:</strong> Number of DE genes, DAS genes, DE transcripts and DTU transcripts, and bar plots of up-down regulation numbers.
<p>
![](3D_App_figure/venn.png)
<p id='Figure11'>
<strong>Figure 11:</strong> Euler diagram of numbers between contrast groups (A) and numbers to compare DE vs DAS genes and DT vs DTU transcripts (B).
<p>
![](3D_App_figure/flowchat.png)
<p id='Figure12'>
<strong>Figure 12:</strong> Flow charts to show numbers of DE and/or DAS genes and DE and/or DTU transcripts.
<p>
<br>

[Go back to Table of contents](#table-of-contents)

### Advanced plot

##### Heatmap

Users can make heatmap for DE genes, DAS genes, DE transcripts and DTU transcripts identified from the analysis (<a href='#Figure13'>Figure 13A</a>). Users can also upload a gene or transcript list in csv file to make heatmap for specific list (<a href='#Figure13'>Figure 13B</a>). To make the heatmap, the average TPMs of conditions for the targets are standardized into z-scores. In corresponding plots, the DAS only genes or DTU only transcripts are excluded since they only have AS changes rather than significant abundance changes. The <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/dist.html' target='_blank'>dist</a> and <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html' target='_blank'>hclust</a> R functions are used to perform the hierarchical clustering analysis. The heatmaps and the target list in each cluster of the heatmaps can be saved to local folder.

![](3D_App_figure/heatmap.png)
<p id='Figure13'>
<strong>Figure 13:</strong> Heatmap panel (A) and example of user provided gene list in csv file to make heatmap (B). To avoid confliction, if the radio button "Select target" is selected, the "Or Choose target list csv file" function is enabled; If the radio button "Upload target list" is selected, the "Select targets" function is enabled.
<p>
<br>

##### Profile plot

Gene and transcript expression profiles (TPMs or read counts) and PS (percent spliced) can be visualised by typing a gene name (<a href='#Figure14'>Figure 14A</a>). Users can also generate a number of plots by provides a gene list (<a href='#Figure13'>Figure 13B</a>) and all the plots will be saved to "Multiple profile plot" folder in the "figure" folder of working directory (<a href='#Figure14'>Figure 14B</a>).

![](3D_App_figure/profiles.png)
<p id='Figure14'>
<strong>Figure 14:</strong> Plot of expression profiles and PS across conditions
<p>
<br>

##### GO annotation plot

Users can generate DE and DAS gene list by click "Generate" button (<a href='#Figure15'>Figure 15B</a>). These gene lists can be loaded to Gene Ontology (GO) analysis tools/databases (e.g. <a href='https://david.ncifcrf.gov/' target='_blank'>DAVID</a> and <a href='http://bioinfo.cau.edu.cn/agriGO/' target='_blank'>agriGO</a>) to generate GO annotation. A csv file with GO annotation information is required to generate the annotation plot (<a href='#Figure15'>Figure 15B</a>). The file includes a "Category" column of CC (cellular component), BP (biological process) and MF (molecular function), a column of "Term" of GO annotation and columns of statistics to report the annotation enrichment, e.g. count, FDR, -log10(FDR), etc. (<a href='#Figure15'>Figure 15A</a>). In the panel, users can select different statistics to visualise. The selected gene list type will present in the plot title to distinguish whether the provided gene list is DE or DAS genes.

![](3D_App_figure/go.png)
<p id='Figure15'>
<strong>Figure 15:</strong> GO annotation plot panel.
<p>
<br>

[Go back to Table of contents](#table-of-contents)

### Generate report

The table in this page lists the parameters for 3D analysis and the locations of saved files, which will be used to generate the report. If the information is wrong, users can double click the cells in the table to correct to values. Based on this information, reports in three formats, word, pdf and html, will be generated and saved in the "report" folder of working directory.

![](3D_App_figure/report.png)
<p id='Figure16'>
<strong>Figure 16:</strong> Generate report panel.
<p>
<p id="Section_saved_files">
</p>
<br>

Saved files in local directory
------------------------------

### Saved files in the "data" folder

<table>
<caption>The intermediate datasets in &quot;.RData&quot; format are saved in the &quot;data&quot; folder. Advanced R users can access these &quot;.RData&quot; objects by R command line. For the returned objects from R functions/packages, please go to corresponding documentations for details.</caption>
<colgroup>
<col width="10%" />
<col width="89%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">File names</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">contrast.RData</td>
<td align="left">Contrast group</td>
</tr>
<tr class="even">
<td align="left">DAS_genes.RData</td>
<td align="left">Testing statistics of DAS genes.</td>
</tr>
<tr class="odd">
<td align="left">DE_genes.RData</td>
<td align="left">Testing statistics of DE genes.</td>
</tr>
<tr class="even">
<td align="left">DE_trans.RData</td>
<td align="left">Testing statistics of DE transcripts.</td>
</tr>
<tr class="odd">
<td align="left">deltaPS.RData</td>
<td align="left">Matrix of deltaPS values.</td>
</tr>
<tr class="even">
<td align="left">DTU_trans.RData</td>
<td align="left">Testing statistics of DTU transcripts.</td>
</tr>
<tr class="odd">
<td align="left">genes_3D_stat.RData</td>
<td align="left">The returned object of DE gene analysis by using limma.pipeline/edgeR.pipeline functions in ThreeDRNAseq R package.</td>
</tr>
<tr class="even">
<td align="left">genes_batch.RData</td>
<td align="left">The returned object of gene level batch effect estimation by using RUVSeq R package.</td>
</tr>
<tr class="odd">
<td align="left">genes_dge.RData</td>
<td align="left">The returned object of gene expression normalization by using calcNormFactors function in edgeR.</td>
</tr>
<tr class="even">
<td align="left">mapping.RData</td>
<td align="left">Matrix of transcript-gene mapping.</td>
</tr>
<tr class="odd">
<td align="left">PS.RData</td>
<td align="left">Matrix of transcript percent spliced (PS), which is used to generate deltaPS based on contrast groups.</td>
</tr>
<tr class="even">
<td align="left">samples.RData</td>
<td align="left">The sample information before merging sequencing replicates.</td>
</tr>
<tr class="odd">
<td align="left">samples_new.RData</td>
<td align="left">The sample information after merging sequencing replicates. If there are no sequencing replicates, the &quot;sample_new&quot; object has the same number of rows with the &quot;samples&quot; object.</td>
</tr>
<tr class="even">
<td align="left">target_high.RData</td>
<td align="left">A list object with two element: expressed transcripts (trans_high) and expressed genes (genes_high).</td>
</tr>
<tr class="odd">
<td align="left">trans_3D_stat.RData</td>
<td align="left">The returned object of DAS gene, DE and DTU transcripts analysis by using limma.pipeline/edgeR.pipeline functions in ThreeDRNAseq R package.</td>
</tr>
<tr class="even">
<td align="left">trans_batch.RData</td>
<td align="left">The returned object of transcript level batch effect estimation by using RUVSeq R package..</td>
</tr>
<tr class="odd">
<td align="left">trans_dge.RData</td>
<td align="left">The returned object of transcript expression normalization by using calcNormFactors function in edgeR.</td>
</tr>
<tr class="even">
<td align="left">txi_genes.RData</td>
<td align="left">The returned object of gene expression by using tximport R package.</td>
</tr>
<tr class="odd">
<td align="left">txi_trans.RData</td>
<td align="left">The returned object of transcript expression by using tximport R package.</td>
</tr>
</tbody>
</table>

### Saved files in the "result" folder

<table>
<caption>The csv files in the result folder of working directory. Raw read counts/TPMs are the datasets before any data pre-processing, e.g. sequencing replicate merge and low expression filters.</caption>
<colgroup>
<col width="24%" />
<col width="75%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">File.names</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">counts_genes.csv</td>
<td align="left">Gene level raw read counts</td>
</tr>
<tr class="even">
<td align="left">counts_trans.csv</td>
<td align="left">Transcript level raw read counts</td>
</tr>
<tr class="odd">
<td align="left">DAS genes.csv</td>
<td align="left">Testing statistics of DAS genes.</td>
</tr>
<tr class="even">
<td align="left">data.info.csv</td>
<td align="left">Data information during data pre-processing.</td>
</tr>
<tr class="odd">
<td align="left">DE DAS DTU numbers.csv</td>
<td align="left">Numbers of 3D genes/transcripts in contrast groups.</td>
</tr>
<tr class="even">
<td align="left">DE genes.csv</td>
<td align="left">Testing statistics of DE genes.</td>
</tr>
<tr class="odd">
<td align="left">DE transcripts.csv</td>
<td align="left">Testing statistics of DE transcripts.</td>
</tr>
<tr class="even">
<td align="left">DE vs DAS gene number.csv</td>
<td align="left">Numbers of DE vs DAS genes in contrast groups.</td>
</tr>
<tr class="odd">
<td align="left">DE vs DTU transcript number.csv</td>
<td align="left">Numbers of DE vs DTU transcripts in contrast groups.</td>
</tr>
<tr class="even">
<td align="left">DTU transcripts.csv</td>
<td align="left">Testing statistics of DTU transcripts.</td>
</tr>
<tr class="odd">
<td align="left">Parameter summary.csv</td>
<td align="left">Methods/Parameters/Cut-offs used for 3D RNA-seq analysis.</td>
</tr>
<tr class="even">
<td align="left">Target in each cluster heatmap * DE&amp;DAS genes.csv</td>
<td align="left">DE&amp;DAS gene lists in clusters of DAS gene heatmap. The DASonly genes are excluded since they have no significant abundance changes across samples.</td>
</tr>
<tr class="odd">
<td align="left">Target in each cluster heatmap * DE genes.csv</td>
<td align="left">DE gene list in clusters of DE gene heatmap.</td>
</tr>
<tr class="even">
<td align="left">Target in each cluster heatmap * DE trans.csv</td>
<td align="left">DE&amp;DTU transcript lists in clusters of DTU transcript heatmap. The DTUonly transcripts are excluded since they have no significant abundance changes across samples.</td>
</tr>
<tr class="odd">
<td align="left">Target in each cluster heatmap * DE&amp;DTU trans.csv</td>
<td align="left">DE transcript list in clusters of DE transcript heatmap.</td>
</tr>
<tr class="even">
<td align="left">TPM_genes.csv</td>
<td align="left">Gene level raw TPMs</td>
</tr>
<tr class="odd">
<td align="left">TPM_trans.csv</td>
<td align="left">Transcript level raw TPMs</td>
</tr>
</tbody>
</table>

### Saved files in the "figure" folder

<table>
<caption>Figures are saved to &quot;figure&quot; folder in both &quot;png&quot; and &quot;pdf&quot; formats.</caption>
<colgroup>
<col width="54%" />
<col width="45%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">File.names</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">DAS genes GO annotation plot.png/.pdf</td>
<td align="left">DAS genes GO annotation plot</td>
</tr>
<tr class="even">
<td align="left">DAS genes updown regulation numbers.png/.pdf</td>
<td align="left">DAS genes updown regulation numbers</td>
</tr>
<tr class="odd">
<td align="left">DE genes GO annotation plot.png/.pdf</td>
<td align="left">DE genes GO annotation plot</td>
</tr>
<tr class="even">
<td align="left">DE genes updown regulation numbers.png/.pdf</td>
<td align="left">DE genes updown regulation numbers</td>
</tr>
<tr class="odd">
<td align="left">DE transcripts updown regulation numbers.png/.pdf</td>
<td align="left">DE transcripts updown regulation numbers</td>
</tr>
<tr class="even">
<td align="left">DTU transcripts updown regulation numbers.png/.pdf</td>
<td align="left">DTU transcripts updown regulation numbers</td>
</tr>
<tr class="odd">
<td align="left">Gene data distribution.png/.pdf</td>
<td align="left">Gene data distribution</td>
</tr>
<tr class="even">
<td align="left">Gene mean-variance trend.png/.pdf</td>
<td align="left">Gene mean-variance trend</td>
</tr>
<tr class="odd">
<td align="left">Gene PCA Average expression.png/.pdf</td>
<td align="left">Gene PCA Average expression</td>
</tr>
<tr class="even">
<td align="left">Gene PCA batch effect removed Bio-reps.png/.pdf</td>
<td align="left">Gene PCA batch effect removed Bio-reps</td>
</tr>
<tr class="odd">
<td align="left">Gene PCA Bio-reps.png/.pdf</td>
<td align="left">Gene PCA Bio-reps</td>
</tr>
<tr class="even">
<td align="left">Heatmap DAS genes.png/.pdf</td>
<td align="left">Heatmap DAS genes</td>
</tr>
<tr class="odd">
<td align="left">Heatmap DE genes.png/.pdf</td>
<td align="left">Heatmap DE genes</td>
</tr>
<tr class="even">
<td align="left">Heatmap DE transcripts.png/.pdf</td>
<td align="left">Heatmap DE transcripts</td>
</tr>
<tr class="odd">
<td align="left">Heatmap DTU transcripts.png/.pdf</td>
<td align="left">Heatmap DTU transcripts</td>
</tr>
<tr class="even">
<td align="left">Transcript data distribution.png/.pdf</td>
<td align="left">Transcript data distribution</td>
</tr>
<tr class="odd">
<td align="left">Transcript mean-variance trend.png/.pdf</td>
<td align="left">Transcript mean-variance trend</td>
</tr>
<tr class="even">
<td align="left">Transcript PCA Average expression.png/.pdf</td>
<td align="left">Transcript PCA Average expression</td>
</tr>
<tr class="odd">
<td align="left">Transcript PCA batch effect removed Bio-reps.png/.pdf</td>
<td align="left">Transcript PCA batch effect removed Bio-reps</td>
</tr>
<tr class="even">
<td align="left">Transcript PCA Bio-reps.png/.pdf</td>
<td align="left">Transcript PCA Bio-reps</td>
</tr>
<tr class="odd">
<td align="left">Union set DE genes vs DAS genes.png/.pdf</td>
<td align="left">Flow chart -Union set DE genes vs DAS genes</td>
</tr>
<tr class="even">
<td align="left">Union set DE transcripts vs DTU transcripts.png/.pdf</td>
<td align="left">Flow chart -Union set DE transcripts vs DTU transcripts</td>
</tr>
</tbody>
</table>

### Saved files in the "report" folder

| File.names  | Description |
|:------------|:------------|
| report.docx | word format |
| report.html | html format |
| report.pdf  | pdf format  |

[Go back to Table of contents](#table-of-contents)

References
----------

Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525–527.

Bullard,J.H., Purdom,E., Hansen,K.D., and Dudoit,S. (2010) Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics, 11, 94.

Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H.G., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. bioRxiv, 251876.

Law,C.W., Chen,Y., Shi,W., and Smyth,G.K. (2014) voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol, 15, R29.

Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.

Risso,D., Ngai,J., Speed,T.P., and Dudoit,S. (2014) Normalization of RNA-seq data using factor analysis of control genes or samples. Nat. Biotechnol., 32, 896–902.

Ritchie,M.E., Phipson,B., Wu,D., Hu,Y., Law,C.W., Shi,W., and Smyth,G.K. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 43, e47.

Robinson,M., Mccarthy,D., Chen,Y., and Smyth,G.K. (2011) edgeR : differential expression analysis of digital gene expression data User ’ s Guide. Most, 23, 1–77.

Soneson,C., Love,M.I., and Robinson,M.D. (2016) Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4, 1521.

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.5.1  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
    ##  [5] tools_3.5.1     htmltools_0.3.6 yaml_2.2.0      Rcpp_1.0.0     
    ##  [9] stringi_1.1.7   rmarkdown_1.10  knitr_1.20      stringr_1.3.1  
    ## [13] digest_0.6.18   evaluate_0.12

[Go back to Table of contents](#table-of-contents)

</div>
