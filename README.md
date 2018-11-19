<img src="vignettes/fig/header.png" align="right"/> <br> <br> <br>

Description
===========

This package provides an interactive graphical user interface (GUI) for RNA-seq data differential expression (DE), differential alternative splicing (DAS) and differential transcript usage (DTU) analyses based on two popular pipelines: limma (Smyth et al. 2013) and edgeR (Robinson et al., 2010). The 3D RNAseq GUI is based on R shiny App and enables command-line-free analysis. To perform analysis, the first step is to generate transcript quantification from quantification tools, such as Salmon (Patro et al., 2017) and Kallisto (Bray et al., 2016). Then users can do mouse click on the App to upload transcript read counts, perform DE and DAS analysis, and make beautiful plots, e.g. expression mean-variance trend plots, PCA plots, heatmap, GO annotation plots, etc. The GUI has steps of proper data pre-processing and false positive controls. These lead to robust DE DAS predictions. All the results and plots can be saved to a report in html, pdf or word format. The analysis pipeline has been successfully used in different RNA-seq studies from Arabidopsis (Calixto et al., 2018), barley (Bull et al., 2017) and potato.

Installation and loading
========================

Install ThreeDRNAseq package
----------------------------

The package can be installed from Github by using <a href='https://cran.r-project.org/web/packages/devtools/index.html' target='_blank'>devtools</a> R package

``` r
#######################################################################################################
## use devtools R package to install ThreeDRNAseq from Github
###---> If devtools is not installed, please install
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages('devtools')
devtools::instasll_github('wyguo/ThreeDRNAseq')
```

Install dependency packages
---------------------------

``` r
#######################################################################################################
## Install packages of dependency
###---> Install packages from Cran
cran.package.list <- c('shiny','shinydashboard','shinydashboard','shinyFiles','plotly','eulerr',
                       'gridExtra','Gmisc')
for(i in cran.package.list){
   if(!requireNamespace(i, quietly = TRUE)){
     message('Installing package: ',i)
     install.packages(i,dependencies = T)
   } else next
}

###---> Install packages from Bioconductor
bioconductor.package.list <- c('tximport','edgeR','limma','RUVSeq','ComplexHeatmap')
for(i in bioconductor.package.list){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!requireNamespace(i, quietly = TRUE)){
    message('Installing package: ',i)
    BiocManager::install("tximport", version = "3.8")
  } else next
}
```

Note: if any other R packages are missing from your PC, please install accordingly.

Launch the 3D RNA-seq App
=========================

Once all the dependency packages are installed, in the RStudio, typing

``` r
library(ThreeDRNAseq)
ThreeDRNAseq.app()
```

User manuals
============

Two versions of step-by-step user manuals are provided:

-   3D RNA-seq App <a href="xxx" target="_blank">user manual</a> (command-line free).
-   Command-line <a href="xxx" target="_blank">user manual</a>. Advanced R users also can use command line to perform the analysis.

References
==========

Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525–527.

Bull,H., Casao,M.C., Zwirek,M., Flavell,A.J., Thomas,W.T.B.B., Guo,W., Zhang,R., Rapazote-Flores,P., Kyriakidis,S., Russell,J., Druka,A., McKim,S.M., and Waugh,R. (2017) Barley SIX-ROWED SPIKE3 encodes a putative Jumonji C-type H3K9me2/me3 demethylase that represses lateral spikelet fertility. Nat. Commun., 8, 936.

Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell, tpc.00177.2018.

Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.

Robinson,M.D., McCarthy,D.J., and Smyth,G.K. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–40.

Smyth,G.K., Ritchie,M., Thorne,N., Wettenhall,J., and Shi,W. (2013) limma:Linear Models for Microarray Data User’s Guide(Now Including RNA-Seq Data Analysis). R Man., 1–123.
