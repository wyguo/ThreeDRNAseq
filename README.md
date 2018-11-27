
<img src="vignettes/fig/header.gif" align="right"/>

Description
===========

<div align="justify">
ThreeDRNAseq (3D RNA-seq) R package provides an interactive graphical user interface (GUI) for RNA-seq data differential expression (DE), differential alternative splicing (DAS) and differential transcript usage (DTU) (3D) analyses based on two popular pipelines: <a href="https://bioconductor.org/packages/release/bioc/html/limma.html" target="_blank">limma</a> (Smyth et al. 2013) and <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">edgeR</a> (Robinson et al., 2010). The 3D RNA-seq GUI is based on R <a href="https://shiny.rstudio.com/" target="_blank">shiny</a> App and enables the RNA-seq analysis to be done only with in 3 Days (3D). To perform the 3D analysis,

-   The first step is to generate transcript quantification from tools, such as <a href="https://combine-lab.github.io/salmon/" target="_blank">Salmon</a> (Patro et al., 2017) and <a href="https://pachterlab.github.io/kallisto/" target="_blank">Kallisto</a> (Bray et al., 2016). The transcript quantification procedure often takes **1-2 Days**.
-   Then users can do mouse click on the App to upload transcript read counts, perform 3D analysis, and make beautiful plots, e.g. expression mean-variance trend plots, PCA plots, heatmap, GO annotation plots, etc.. The results of significance and all the plots can be wrapped into html, pdf and/or word reports by one-button-click. The entire 3D mouse-click analysis only takes **1 Day or even less**.

The 3D RNA-seq analysis pipeline has steps of proper data pre-processing, e.g. low expression filters based on expression mean-variance trend, visualise data variation, batch effect estimation, data normalisation, etc.. The optimal parameters for each step are determined by quality control plots, e.g. mean-variance trend plots, PCA plots, data distribution plots, etc.. The pipeline also has stringent controls of false positive, leading to robust 3D predictions. The pipeline has been successfully applied in different RNA-seq studies from Arabidopsis (Calixto et al., 2018), barley (Bull et al., 2017) and potato to identify novel condition-responsive genes/transcripts, especially those with significant alternative splicing changes. To use our pipeline in your work, please cite:

[Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.](http://www.plantcell.org/content/30/7/1424){target="_blank"}


User manuals
============

Two versions of step-by-step user manuals are provided:

-   3D RNA-seq App manual (command-line free) for easy to use:

[https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md](https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md){target="_blank"}

-   3D RNA-seq command-line based manual for advanced R users:

[https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/Command_line_based_manul.md"](https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/Command_line_based_manul.md){target="_blank"}

References
==========

Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525–527.

Bull,H., Casao,M.C., Zwirek,M., Flavell,A.J., Thomas,W.T.B.B., Guo,W., Zhang,R., Rapazote-Flores,P., Kyriakidis,S., Russell,J., Druka,A., McKim,S.M., and Waugh,R. (2017) Barley SIX-ROWED SPIKE3 encodes a putative Jumonji C-type H3K9me2/me3 demethylase that represses lateral spikelet fertility. Nat. Commun., 8, 936.

Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell, tpc.00177.2018.

Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.

Robinson,M.D., McCarthy,D.J., and Smyth,G.K. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–40.

Smyth,G.K., Ritchie,M., Thorne,N., Wettenhall,J., and Shi,W. (2013) limma:Linear Models for Microarray Data User’s Guide(Now Including RNA-Seq Data Analysis). R Man., 1–123.
</div>
