3D RNA-seq App
==============
<i><h4>Command-line based user manual</h4></i>
<i><h4>Wenbin Guo</h4></i>
<i><h4>28 May 2019</h4></i>
<i><h4>Division of Plant Sciences, University of Dundee at the James Hutton Institute, Dundee DD2 5DA, UK</h4></i>

Table of contents
-----------------

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Load R package](#load-r-package)
-   [Example data](#example-data)
-   [Set pipeline parameters](#set-pipeline-parameters)
-   [Data generation](#data-generation)
    -   [Read files](#read-files)
    -   [Run tximport](#run-tximport)
-   [Data pre-processing](#data-pre-processing)
    -   [Step 1: Merge sequencing replicates](#step-1-merge-sequencing-replicates)
    -   [Step 2: Filter low expression transcripts/genes](#step-2-filter-low-expression-transcriptsgenes)
    -   [Step 3: Principal component analysis (PCA)](#step-3-principal-component-analysis-pca)
    -   [Step 4: Batch effect estimation](#step-4-batch-effect-estimation)
    -   [Step 4: Data normalization](#step-4-data-normalization)
    -   [Data information](#data-information)
-   [DE, DAS and DTU analysis](#de-das-and-dtu-analysis)
    -   [Step 1: Set contrast group](#step-1-set-contrast-group)
    -   [Step 2: DE genes](#step-2-de-genes)
    -   [Step 3: DAS genes, DE and DTU transcripts](#step-3-das-genes-de-and-dtu-transcripts)
    -   [Step 4 Result summary](#step-4-result-summary)
        -   [Significant 3D targets and testing statistics](#significant-3d-targets-and-testing-statistics)
        -   [Significant 3D target numbers](#significant-3d-target-numbers)
    -   [Step 5 Make plot](#step-5-make-plot)
        -   [Up- and down-regulation](#up--and-down-regulation)
        -   [Volcano plot](#volcano-plot)
        -   [Flow chart](#flow-chart)
        -   [Comparisons between contrast groups](#comparisons-between-contrast-groups)
        -   [Comparison of transcription and AS targets](#comparison-of-transcription-and-as-targets)
-   [Functional plot](#functional-plot)
    -   [Heatmap of 3D targets](#heatmap-of-3d-targets)
    -   [Profile plot](#profile-plot)
    -   [GO annotation plot of DE/DAS genes](#go-annotation-plot-of-dedas-genes)
-   [Isoform switch analysis](#isoform-switch-analysis)
    -   [Generate switch scores](#generate-switch-scores)
    -   [Isoform switch plot](#isoform-switch-plot)
        -   [Switch numbers](#switch-numbers)
        -   [Isoform switch plot](#isoform-switch-plot-1)
-   [Generate report](#generate-report)
    -   [Summary parameters](#summary-parameters)
    -   [Generate report](#generate-report-1)
    -   [Save significant results](#save-significant-results)
-   [Saved files in local directory](#saved-files-in-local-directory)
    -   [Saved files in the "data" folder](#saved-files-in-the-data-folder)
    -   [Saved files in the "result" folder](#saved-files-in-the-result-folder)
    -   [Saved files in the "figure" folder](#saved-files-in-the-figure-folder)
    -   [Saved files in the "report" folder](#saved-files-in-the-report-folder)
-   [References](#references)
-   [Session information](#session-information)

Introduction
------------

Instead of the 3D RNA-seq App, advanced R users can use command line to perform 3D analysis. This file includes the R codes, which corresponds to the steps in the 3D RNA-seq App. Users can refer to the details in the manual:

<a href="https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md" target="_blank">https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md</a>

To use our pipeline in your work, please cite:

- Guo,W. et al. (2020) 3D RNA-seq: a powerful and flexible tool for rapid and accurate differential expression and alternative splicing analysis of RNA-seq data for biologists. RNA Biol., DOI: 10.1080/15476286.2020.1858253.
- Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H.G., Zhang,R., and Brown,J.W.S. (2018) Rapid and Dynamic Alternative Splicing Impacts the Arabidopsis Cold Response Transcriptome. Plant Cell, 30, 1424–1444.

**3D analysis pipeline**

![](3D_App_figure/pipeline.png)

Installation
------------

***Install dependency packages***

``` r
#####################################################################################
## Install packages of dependency
##----->> Install packages from Cran
cran.package.list <- c("plotly","ggplot2","eulerr",
                       "gridExtra","grid","fastcluster","base64enc","ggrepel")
for(i in cran.package.list){
   if(!(i %in% rownames(installed.packages()))){
     message('Installing package: ',i)
     install.packages(i,dependencies = T)
   } else next
}

##----->> Install packages from Bioconductor
bioconductor.package.list <- c('tximport','edgeR','limma','RUVSeq',
                               'ComplexHeatmap','rhdf5')
for(i in bioconductor.package.list){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!(i %in% rownames(installed.packages()))){
    message('Installing package: ',i)
    BiocManager::install(i, version = "3.8",dependencies = T)
  } else next
}
```

***Install ThreeDRNAseq R package***

ThreeDRNAseq R package can be installed from Github by using <a href='https://cran.r-project.org/web/packages/devtools/index.html' target='_blank'>devtools</a> R package

``` r
##################################################################################################
## use devtools R package to install ThreeDRNAseq from Github
##----->> If devtools is not installed, please install
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages('devtools',dependencies = TRUE)

##----->> Install ThreeDRNAseq
if(!requireNamespace("ThreeDRNAseq", quietly = TRUE))
  devtools::install_github('wyguo/ThreeDRNAseq')
```

**NOTE**: To run the ThreeDRNAseq App on a local PC, if any other denpendency R packages are missing, please install them. 

Load R package
--------------

Assume all the required R packages are installed.

``` r
###---> ThreeDRNAseq R package
library(ThreeDRNAseq)

###---> Denpendency R package
library(tximport)
library(edgeR)
library(limma)
library(RUVSeq)
library(eulerr)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)

options(stringsAsFactors=F)
setwd("D:/PhD project/R projects/test round 2019/ThreeDRNAseq_original/example_data")
```

Example data
------------

**Download link**: <https://www.dropbox.com/s/8vwuz6u2yl7v9qx/3D%20RNA-seq%20App%20example%20data.zip?dl=0>

**Description**: This example is a sub-dataset from a time-series study of Arabidopsis plants exposed to cold (Calixto et al., 2018). RNA-seq data of 6 time-points were extracted from the whole dataset. The time-points are 3 and 6 hours after dusk at 20<sup>*o*</sup>*C*, the first day of transition to 4<sup>*o*</sup>*C* and the fourth day of acclimation to 4<sup>*o*</sup>*C* (red boxes in Figure 2). Each time-point has 3 biological replicates and each biological replicate has 3 sequencing replicates. Transcript quantifications were generated using Salmon (Patro et al., 2017) and AtRTD2-QUASI (Zhang et al., 2017) as the reference transcriptome. To further reduce the data size, the expression of 8,944 transcripts (from 2,000 genes ) were extracted from the Salmon quantification to identify the cold response genes and transcripts at both transcriptional and AS level.

![](3D_App_figure/example_exp.png)

<p id="fig_example_exp">
<strong>Figure 2:</strong> Experimental design of time-series RNA-seq data from study of Arabidopsis cold response. For this example, a subset of samples, genes and transcripts were extracted from the whole dataset to reduce data size.
</p>
Set pipeline parameters
-----------------------

``` r
##save to object
DDD.data <- list()
################################################################################
##----->> Set folders to read and save results
data.folder <- file.path(getwd(),'data') # for .RData format results
result.folder <- file.path(getwd(),'result') # for 3D analysis results in csv files
figure.folder <- file.path(getwd(),'figure')# for figures
report.folder <- file.path(getwd(),'report')

DDD.data$data.folder <- data.folder
DDD.data$result.folder <- result.folder
DDD.data$figure.folder <- figure.folder
DDD.data$report.folder <- report.folder

if(!file.exists(data.folder))
  dir.create(path = data.folder,recursive = T)
if(!file.exists(result.folder))
  dir.create(path = result.folder,recursive = T)
if(!file.exists(figure.folder))
  dir.create(path = figure.folder,recursive = T)
if(!file.exists(report.folder))
  dir.create(path = report.folder,recursive = T)

### Set the input data folder
##----->> folder of input files
input.folder <- '3D RNA-seq App example data'
quant.folder <- '3D RNA-seq App example data/quant'

################################################################################
##----->> parameters of tximport to generate read counts and TPMs
quant_method <- 'salmon' # abundance generator
tximport_method <- 'lengthScaledTPM' # method to generate expression in tximport

################################################################################
##----->> parameters for data pre-processing
### has sequencign replicates?
has_srep <- T

### parameter for low expression filters
cpm_cut <- 1
cpm_samples_n <- 3

### parameter for batch effect estimation
has_batcheffect <- T
RUVseq_method <- 'RUVr' # RUVseq_method is one of 'RUVr', 'RUVs' and 'RUVg'

### data normalisation parameter
norm_method <- 'TMM' ## norm_method is one of 'TMM','RLE' and 'upperquartile'

################################################################################
##----->> parameters for 3D analysis
pval_adj_method <- 'BH'
pval_cut <- 0.01
l2fc_cut <- 1
DE_pipeline <- 'limma'
deltaPS_cut <- 0.1
DAS_pval_method <- 'F-test'

################################################################################
##----->> heatmap
dist_method <- 'euclidean'
cluster_method <- 'ward.D'
cluster_number <- 10

################################################################################
##----->> TSIS
TSISorisokTSP <- 'isokTSP'
TSIS_method_intersection <- method_intersection <- 'mean'
TSIS_spline_df <- spline_df <- NULL
TSIS_prob_cut <- 0.5
TSIS_diff_cut <- 1
TSIS_adj_pval_cut <- 0.05
TSIS_time_point_cut <- 1
TSIS_cor_cut <- 0
```

Data generation
---------------

Use the <a href="http://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank">tximport</a> R package (Soneson et al., 2016) to generate the transcript and gene read counts and TPMs (transcript per million reads) from transcript quantification outputs, such as Salmon (Patro et al., 2017) and Kallisto (Bray et al., 2016).

### Read files

``` r
################################################################################
##----->> Meta table includes sample information, e.g. conditions, bio-reps, seq-reps, abundance paths, etc.
metatable <- read.csv(file.path(getwd(),'metatable.csv'))
##select the columns of experimental design
factor_col <- c('time','day')
brep_col <- 'bio_rep'
srep_col <- 'seq_rep'
quant_col <- 'quant_files'

##arrange information in the metatable
metatable$label <- as.vector(interaction(metatable[,factor_col]))
metatable$sample.name <- as.vector(interaction(metatable[,c(factor_col,brep_col,srep_col)]))
metatable$quant.folder <-  file.path(quant.folder,metatable$quant_files,
                                     ifelse(quant_method=='salmon','quant.sf','abundance.h5'))

##----->> Transcript-gene association mapping
mapping <-read.csv(file.path(getwd(),'mapping.csv'))
mapping <- data.frame(as.matrix(mapping),stringsAsFactors = F)
rownames(mapping) <- mapping$TXNAME
```

### Run tximport

``` r
################################################################################
##----->> Generate gene expression
##
txi_genes <- tximport(metatable$quant.folder,dropInfReps = T,
                      type = quant_method, tx2gene = mapping,
                      countsFromAbundance = tximport_method)

## give colunames to the datasets
colnames(txi_genes$counts) <-
  colnames(txi_genes$abundance) <-
  colnames(txi_genes$length) <-metatable$sample.name

## save the data
write.csv(txi_genes$counts,file=paste0(result.folder,'/counts_genes.csv'))
write.csv(txi_genes$abundance,file=paste0(result.folder,'/TPM_genes.csv'))
save(txi_genes,file=paste0(data.folder,'/txi_genes.RData'))

################################################################################
##----->> Generate transcripts expression
txi_trans<- tximport(metatable$quant.folder, 
                     type = quant_method, tx2gene = NULL,
                     countsFromAbundance = tximport_method,
                     txOut = T,dropInfReps = T)

## give colunames to the datasets
colnames(txi_trans$counts) <- 
  colnames(txi_trans$abundance) <-
  colnames(txi_trans$length) <-metatable$sample.name

## save the data
write.csv(txi_trans$counts,file=paste0(result.folder,'/counts_trans.csv'))
write.csv(txi_trans$abundance,file=paste0(result.folder,'/TPM_trans.csv'))
save(txi_trans,file=paste0(data.folder,'/txi_trans.RData'))

################################################################################
##extract gene and transcript read counts
genes_counts <- txi_genes$counts
trans_counts <- txi_trans$counts
trans_TPM <- txi_trans$abundance
```

Data pre-processing
-------------------

### Step 1: Merge sequencing replicates

If the RNA-seq data includes sequencing replicates (i.e. the same biological replicate/library is sequenced multiple times), the sequencing replicates are merged to increase the library depths, since they do not provided useful information for biological/technical variability.

``` r
##If no sequencing replicates, genes_counts and trans_counts remain the same by 
if(has_srep){
  idx <- paste0(metatable$label,'.',metatable[,brep_col])
  genes_counts <- sumarrays(genes_counts,group = idx)
  trans_counts <- sumarrays(trans_counts,group = idx)
  metatable_new <- metatable[metatable[,srep_col]==metatable[,srep_col][1],]
} else {
  metatable_new <- metatable
}
```

### Step 2: Filter low expression transcripts/genes

To filter the low expressed transcripts across the samples, the expression mean-variance trend is used to determine the CPM (count per million reads) cut-off for low expression and the number of samples to cut. A gene is expressed if any of its transcript is expressed. Please see the details in the 3D RNA-seq App user manual.

``` r
################################################################################
##----->> Do the filters
target_high <- low.expression.filter(abundance = trans_counts,
                                     mapping = mapping,
                                     abundance.cut = cpm_cut,
                                     sample.n = cpm_samples_n,
                                     unit = 'counts',
                                     Log=F)
##save expressed genes and transcripts
save(target_high,file=paste0(data.folder,'/target_high.RData'))

################################################################################
##----->> Mean-variance plot
## transcript level

counts.raw = trans_counts[rowSums(trans_counts>0)>0,]
counts.filtered = trans_counts[target_high$trans_high,]
mv.trans <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = metatable_new$label)
### make plot
fit.raw <- mv.trans$fit.raw
fit.filtered <- mv.trans$fit.filtered
mv.trans.plot <- function(){
  par(mfrow=c(1,2))
  plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
  title('\n\nRaw counts (transcript level)')
  plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     l = fit.filtered$l,lwd=2,col='black')
  title('\n\nFiltered counts (transcript level)')
  lines(fit.raw$l, col = "gold",lty=4,lwd=2)
  legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
         legend = c('low-exp removed','low-exp kept'))
}
mv.trans.plot()

### save to figure folder
png(filename = paste0(figure.folder,'/Transcript mean-variance trend.png'),
    width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
mv.trans.plot()
dev.off()

pdf(file = paste0(figure.folder,'/Transcript mean-variance trend.pdf'),
    width = 25/2.54,height = 12/2.54)
mv.trans.plot()
dev.off()

################################################################################
## gene level
counts.raw = genes_counts[rowSums(genes_counts>0)>0,]
counts.filtered = genes_counts[target_high$genes_high,]
mv.genes <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = metatable_new$label)
### make plot
fit.raw <- mv.genes$fit.raw
fit.filtered <- mv.genes$fit.filtered
mv.genes.plot <- function(){
  par(mfrow=c(1,2))
  plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                     l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
  title('\n\nRaw counts (gene level)')
  plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                     l = fit.filtered$l,lwd=2,col='black')
  title('\n\nFiltered counts (gene level)')
  lines(fit.raw$l, col = "gold",lty=4,lwd=2)
  legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
         legend = c('low-exp removed','low-exp kept'))
}
mv.genes.plot()

### save to figure folder
png(filename = paste0(figure.folder,'/Gene mean-variance trend.png'),
    width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
mv.genes.plot()
dev.off()

pdf(file = paste0(figure.folder,'/Gene mean-variance trend.pdf'),
    width = 25/2.54,height = 12/2.54)
mv.genes.plot()
dev.off()
```

### Step 3: Principal component analysis (PCA)

The PCA plot is used to visualize expression variation across conditions/samples and to determine whether the data has batch effects between different biological replicates.

``` r
################################################################################
##----->> trans level
data2pca <- trans_counts[target_high$trans_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'Transcript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()


################################################################################
##----->> genes level
data2pca <- genes_counts[target_high$genes_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots

groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'genescript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()
```

### Step 4: Batch effect estimation

The batch effect estimation step is only executed when there are distinct bath effects in the data: (1) The biological replicates of the same conditions stay in separate clusters in the PCA plot and There is experimental explanation of this separation (e.g. bio-reps in different time/labs). Otherwise skip this step to avoid data over-modification. The <a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" target="_blank">RUVSeq</a> (Risso et al., 2014) is used to estimate the batch effects. It will generate modified read counts, where the batch effects have been removed, for batch-free PCA plots and batch effect term which can be passed to the linear regression model for 3D analysis.

``` r
design <- condition2design(condition = metatable_new$label,
                           batch.effect = NULL)

################################################################################
##----->> trans level
trans_batch <- remove.batch(read.counts = trans_counts[target_high$trans_high,],
                            condition = metatable_new$label,
                            design = design,
                            contrast=NULL,
                            group = metatable_new$label,
                            method = RUVseq_method)
save(trans_batch,file=paste0(data.folder,'/trans_batch.RData')) 

################################################################################
##----->> genes level
genes_batch <- remove.batch(read.counts = genes_counts[target_high$genes_high,],
                            condition = metatable_new$label,
                            design = design,
                            contrast=NULL,
                            group = metatable_new$label,
                            method = RUVseq_method)
save(genes_batch,file=paste0(data.folder,'/genes_batch.RData')) 


################################################################################
## DO the PCA again
################################################################################

##----->> trans level
data2pca <- trans_batch$normalizedCounts[target_high$trans_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'Transcript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA batch effect removed Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA batch effect removed Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'Transcript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA batch effect removed average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA batch effect removed average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()


################################################################################
##----->> genes level
data2pca <- genes_batch$normalizedCounts[target_high$genes_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                  groups = groups,plot.title = 'genescript PCA: bio-reps',
                  ellipse.type = ellipse.type,
                  add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA batch effect removed Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA batch effect removed Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'genescript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA batch effect removed average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA batch effect removed average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()
```

### Step 4: Data normalization

The same transcript/gene in deeper sequenced samples have more mapped reads. To make genes/transcripts comparable across samples, the "TMM", "RLE" or "upperquantile" method can be used to normalize the expression.

``` r
################################################################################
##----->> trans level
dge <- DGEList(counts=trans_counts[target_high$trans_high,],
               group = metatable_new$label,
               genes = mapping[target_high$trans_high,])
trans_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(trans_dge,file=paste0(data.folder,'/trans_dge.RData'))

################################################################################
##----->> genes level
dge <- DGEList(counts=genes_counts[target_high$genes_high,],
               group = metatable_new$label)
genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(genes_dge,file=paste0(data.folder,'/genes_dge.RData'))

################################################################################
##----->> distribution plot
sample.name <- paste0(metatable_new$label,'.',metatable_new[,brep_col])
condition <- metatable_new$label

###--- trans level
data.before <- trans_counts[target_high$trans_high,]
data.after <- counts2CPM(obj = trans_dge,Log = T)
g <- boxplotNormalised(data.before = data.before,
                        data.after = data.after,
                        condition = condition,
                        sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Transcript expression distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript expression distribution.pdf'),
    width = 20/2.54,height = 20/2.54)
do.call(grid.arrange,g)
dev.off()

###--- genes level
data.before <- genes_counts[target_high$genes_high,]
data.after <- counts2CPM(obj = genes_dge,Log = T)
g <- boxplotNormalised(data.before = data.before,
                        data.after = data.after,
                        condition = condition,
                        sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Gene expression distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene expression distribution.pdf'),
    width = 20/2.54,height = 20/2.54)
do.call(grid.arrange,g)
dev.off()
```

### Data information

``` r
RNAseq_info <- data.frame(
  Description=c('Raw transcripts',
                'Raw genes',
                'Samples',
                'Samples after merging seq-reps',
                'Condition of interest',
                'CPM cut-off',
                'Min samples to CPM cut-off',
                'Expressed transcripts',
                'Expressed genes'),
  Number=c(length(mapping$TXNAME),
           length(unique(mapping$GENEID)),
           nrow(metatable),
           nrow(metatable_new),
           length(unique(metatable$label)),
           cpm_cut,
           cpm_samples_n,
           length(target_high$trans_high),
           length(target_high$genes_high))
)
DDD.data$RNAseq_info <- RNAseq_info

RNAseq_info
```

DE, DAS and DTU analysis
------------------------

Three pipelines, <a href="https://bioconductor.org/packages/release/bioc/html/limma.html" target="_blank">limma</a> (Ritchie et al., 2015), <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">glmQL (edgeR)</a> and <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">glm (edgeR)</a> (Robinson et al., 2010) are provided for DE gene and transcripts analysis. "limma" and "edgeR" have comparable performance and they have stringent controls of false positives than the "glm".

### Step 1: Set contrast group

A contrast group is a user-defined comparison between two samples or two groups of samples:

1.  Pair-wise comparison of samples: e.g. B-A compares group B to group A; C-B compares group C to group B (Figure A).
2.  Compare mean of multiple samples: e.g. (WT.A+WT.B)/2-(MU.A+MU.B)/2 compares the mean of group WT.A and WT.B to the mean of group MU.A and MU.B (Figure B).
3.  Compare difference of two differences (interactions): e.g. (WT.A-WT.B)-(MU.A-MU.B) compares the difference (*L*<sub>2</sub>*F**C*) of group WT.A and WT.B to the difference (*L*<sub>2</sub>*F**C*) of group MU.A and MU.B (Figure B).

![](3D_App_figure/select_factor.png)

``` r
################################################################################
##----->> pair-wise contrast groups
contrast_pw <- c('T2.Day1-T2.Day0','T2.Day4-T2.Day0','T3.Day1-T3.Day0','T3.Day4-T3.Day0')

##----->> group mean contrast groups
contrast_mean <- c('(T2.Day1+T3.Day1)/2-(T2.Day0+T3.Day0)/2','(T2.Day4+T3.Day4)/2-(T2.Day0+T3.Day0)/2')

##----->> group differences contrast groups
contrast_pgdiff <- c('(T2.Day1-T3.Day1)-(T2.Day0-T3.Day0)','(T2.Day4-T3.Day4)-(T2.Day0-T3.Day0)')

##----->> put together
contrast <- unique(c(contrast_pw,contrast_mean,contrast_pgdiff))

DDD.data$contrast_pw <- contrast_pw
DDD.data$contrast_mean <- contrast_mean
DDD.data$contrast_pgdiff <- contrast_pgdiff
DDD.data$contrast <- contrast
```

### Step 2: DE genes

Principle of 3D analysis: ![](3D_App_figure/DDD_definition.png)

``` r
batch.effect <- genes_batch$W
# batch.effect <- NULL ## if has no batch effects
design <- condition2design(condition = metatable_new$label,
                           batch.effect = batch.effect)

################################################################################
if(DE_pipeline == 'limma'){
##----->> limma pipeline
genes_3D_stat <- limma.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 adjust.method = pval_adj_method)
}

if(DE_pipeline == 'glmQL'){
##----->> edgeR glmQL pipeline
genes_3D_stat <- edgeR.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 method = 'glmQL',
                                 adjust.method = pval_adj_method)
}

if(DE_pipeline == 'glm'){
##----->> edgeR glm pipeline
genes_3D_stat <- edgeR.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 method = 'glm',
                                 adjust.method = pval_adj_method)
}
## save results
DDD.data$genes_3D_stat <- genes_3D_stat
```

### Step 3: DAS genes, DE and DTU transcripts

``` r
################################################################################
##----->> generate deltaPS
deltaPS <- transAbundance2PS(transAbundance =txi_trans$abundance[target_high$trans_high,],
                             PS = NULL,
                             contrast = contrast,
                             condition = metatable$label,
                             mapping = mapping[target_high$trans_high,])

DDD.data$PS <- PS <- deltaPS$PS
DDD.data$deltaPS <- deltaPS <- deltaPS$deltaPS


################################################################################
##----->> DAS genes,DE and DTU transcripts
batch.effect <- genes_batch$W
# batch.effect <- NULL ## if has no batch effects
design <- condition2design(condition = metatable_new$label,
                           batch.effect = batch.effect)

################################################################################
if(DE_pipeline == 'limma'){
##----->> limma pipeline
trans_3D_stat <- limma.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 adjust.method = pval_adj_method)
}

if(DE_pipeline == 'glmQL'){
##----->> edgeR glmQL pipeline
trans_3D_stat <- edgeR.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 method = 'glmQL',
                                 adjust.method = pval_adj_method)
}

if(DE_pipeline == 'glm'){
##----->> edgeR glm pipeline
trans_3D_stat <- edgeR.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 method = 'glm',
                                 adjust.method = pval_adj_method)
}
## save results
DDD.data$trans_3D_stat <- trans_3D_stat
```

### Step 4 Result summary

After the the 3D analysis, the following information is summarized and will appear at the bottom of the page:

1.  The test statistics in different contrast groups, e.g. adjusted p-value and *L*<sub>2</sub>*F**C*.
2.  The number of genes and transcripts with significant expression changes in contrast groups.
3.  The number of up- and down-regulated DE genes and transcripts.
4.  The numbers of genes/transcripts regulated only by transcription (DE), only by alternative splicing (DAS/DTU) and by both transcription and alternative splicing (DE+DAS/DE+DTU).

These summaries can be filtered, customized and the figures can be saved with specified formats and sizes.

#### Significant 3D targets and testing statistics

``` r
################################################################################
##----->> Summary DE genes
DE_genes <- summaryDEtarget(stat = genes_3D_stat$DE.stat,
                              cutoff = c(adj.pval=pval_cut,
                                         log2FC=l2fc_cut))
DDD.data$DE_genes <- DE_genes

################################################################################
## summary DAS genes, DE and DTU trans
##----->> DE trans
DE_trans <- summaryDEtarget(stat = trans_3D_stat$DE.stat,
                              cutoff = c(adj.pval=pval_cut,
                                         log2FC=l2fc_cut))
DDD.data$DE_trans <- DE_trans

##----->> DAS genes
if(DAS_pval_method=='F-test') {
  DAS.stat <- trans_3D_stat$DAS.F.stat
} else {
  DAS.stat <- trans_3D_stat$DAS.Simes.stat
}

lfc <- genes_3D_stat$DE.lfc
lfc <- reshape2::melt(as.matrix(lfc))
colnames(lfc) <- c('target','contrast','log2FC')
DAS_genes <- summaryDAStarget(stat = DAS.stat,
                                lfc = lfc,
                                cutoff=c(pval_cut,deltaPS_cut))
DDD.data$DAS_genes <- DAS_genes

##----->> DTU trans
lfc <- trans_3D_stat$DE.lfc
lfc <- reshape2::melt(as.matrix(lfc))
colnames(lfc) <- c('target','contrast','log2FC')
DTU_trans <- summaryDAStarget(stat = trans_3D_stat$DTU.stat,
                                lfc = lfc,cutoff = c(adj.pval=pval_cut,
                                                     deltaPS=deltaPS_cut))
DDD.data$DTU_trans <- DTU_trans

################################################################################
## save csv
write.csv(DE_genes,file=paste0(result.folder,'/DE genes.csv'),row.names = F)
write.csv(DAS_genes,file=paste0(result.folder,'/DAS genes.csv'),row.names = F)
write.csv(DE_trans,file=paste0(result.folder,'/DE transcripts.csv'),row.names = F)
write.csv(DTU_trans,file=paste0(result.folder,'/DTU transcripts.csv'),row.names = F)
```

#### Significant 3D target numbers

``` r
################################################################################
##----->> target numbers
DDD_numbers <- summary3Dnumber(DE_genes = DE_genes,
                                 DAS_genes = DAS_genes,
                                 DE_trans = DE_trans,
                                 DTU_trans=DTU_trans,
                                 contrast = contrast)
DDD_numbers
write.csv(DDD_numbers,file=paste0(result.folder,'/DE DAS DTU numbers.csv'),
          row.names = F)
DDD.data$DDD_numbers <- DDD_numbers
################################################################################
##----->> DE vs DAS
DEvsDAS_results <- DEvsDAS(DE_genes = DE_genes,
                           DAS_genes = DAS_genes,
                           contrast = contrast)
DEvsDAS_results
DDD.data$DEvsDAS_results <- DEvsDAS_results
write.csv(DEvsDAS_results,file=paste0(result.folder,'/DE vs DAS gene number.csv'),
          row.names = F)


################################################################################
##----->> DE vs DTU
DEvsDTU_results <- DEvsDTU(DE_trans = DE_trans,
                           DTU_trans = DTU_trans,
                           contrast = contrast)
DEvsDTU_results
DDD.data$DEvsDTU_results <- DEvsDTU_results
write.csv(DEvsDTU_results,file=paste0(result.folder,'/DE vs DTU transcript number.csv'),row.names = F)
```

### Step 5 Make plot

#### Up- and down-regulation

``` r
################################################################################
##----->> DE genes
idx <- factor(DE_genes$contrast,levels = contrast)
targets <-  split(DE_genes,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down-regulated','up-regulated'),number=0)
  } else {
    x <- data.frame(contrast=i,table(targets[[i]]$up.down))
    colnames(x) <- c('contrast','regulation','number')
  }
  x
})
data2plot <- do.call(rbind,data2plot)
g.updown <- plotUpdown(data2plot,plot.title = 'DE genes',contrast = contrast)
print(g.updown)

### save to figure
png(paste0(figure.folder,'/DE genes up and down regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

pdf(paste0(figure.folder,'/DE genes up and down regulation numbers.pdf'),
    width = length(contrast)*5/2.54,10/2.54)
print(g.updown)
dev.off()

################################################################################
##----->> DE trans
idx <- factor(DE_trans$contrast,levels = contrast)
targets <-  split(DE_trans,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down-regulated','up-regulated'),number=0)
  } else {
    x <- data.frame(contrast=i,table(targets[[i]]$up.down))
    colnames(x) <- c('contrast','regulation','number')
  }
  x
})
data2plot <- do.call(rbind,data2plot)
g.updown <- plotUpdown(data2plot,plot.title = 'DE trans',contrast = contrast)
print(g.updown)

### save to figure
png(paste0(figure.folder,'/DE transcripts up and down regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

pdf(paste0(figure.folder,'/DE transcripts up and down regulation numbers.pdf'),
    width = length(contrast)*5/2.54,10/2.54)
print(g.updown)
dev.off()
```

#### Volcano plot

``` r
top.n <- 10
size <- 1
col0 <- 'black'
col1 <- 'red'
idx <- c('DE genes','DAS genes','DE transcripts','DTU transcripts')
title.idx <- c(paste0("Volcano plot: DE genes (Low expression filtered; \nAdjusted p<",
                      pval_cut,'; |L2FC|>=',l2fc_cut,'; Labels: top ',
                      top.n,' distance to (0,0))'),
               paste0("Volcano plot: DAS genes (Low expression filtered; \nAdjusted p<",
                      pval_cut,'; |MaxdeltaPS|>=',deltaPS_cut,'; Labels: top ',
                      top.n,' distance to (0,0))'),
               paste0("Volcano plot: DE transcripts (Low expression filtered; \nAdjusted p<",
                      pval_cut,'; |L2FC|>=',l2fc_cut,'; Labels: top ',
                      top.n,' distance to (0,0))'),
               paste0("Volcano plot: DTU transcripts (Low expression filtered; \nAdjusted p<",
                      pval_cut,'; |deltaPS|>=',deltaPS_cut,'; Labels: top ',
                      top.n,' distance to (0,0))')
)
names(title.idx) <- idx
g <- lapply(idx,function(i){
  if(i=='DE genes'){
    DDD.stat <- genes_3D_stat$DE.stat
    data2plot <- data.frame(target=DDD.stat$target,
                            contrast=DDD.stat$contrast,
                            x=DDD.stat$log2FC,
                            y=-log10(DDD.stat$adj.pval))
    data2plot$significance <- 'Not significant'
    data2plot$significance[DDD.stat$adj.pval< pval_cut & 
                             abs(DDD.stat$log2FC)>=l2fc_cut] <- 'Significant'
    q <- plotVolcano(data2plot = data2plot,xlab = 'log2FC of genes',ylab='-log10(FDR)',
                     title = title.idx[i],
                     col0 = col0,col1 = col1,size = size,top.n = top.n)
  }
  
  if(i=='DAS genes'){
    if(DAS_pval_method=='F-test')
      DDD.stat <- trans_3D_stat$DAS.F.stat else DDD.stat <- trans_3D_stat$DAS.simes.stat
      data2plot <- data.frame(target=DDD.stat$target,
                              contrast=DDD.stat$contrast,
                              x=DDD.stat$maxdeltaPS,
                              y=-log10(DDD.stat$adj.pval))
      data2plot$significance <- 'Not significant'
      data2plot$significance[DDD.stat$adj.pval< pval_cut & 
                               abs(DDD.stat$maxdeltaPS)>=deltaPS_cut] <- 'Significant'
      q <- plotVolcano(data2plot = data2plot,
                       xlab = 'MaxdeltaPS of transcripts in each gene',ylab='-log10(FDR)',
                       title = title.idx[i],col0 = col0,col1 = col1,size = size,
                       top.n = top.n)
  }
  
  if(i=='DE transcripts'){
    DDD.stat <- trans_3D_stat$DE.stat
    data2plot <- data.frame(target=DDD.stat$target,
                            contrast=DDD.stat$contrast,
                            x=DDD.stat$log2FC,
                            y=-log10(DDD.stat$adj.pval))
    data2plot$significance <- 'Not significant'
    data2plot$significance[DDD.stat$adj.pval< pval_cut & 
                             abs(DDD.stat$log2FC)>=l2fc_cut] <- 'Significant'
    q <- plotVolcano(data2plot = data2plot,xlab = 'log2FC of transcripts',
                     ylab='-log10(FDR)',title = title.idx[i],
                     col0 = col0,col1 = col1,size = size,top.n = top.n)
  }
  
  if(i=='DTU transcripts'){
    DDD.stat <- trans_3D_stat$DTU.stat 
    data2plot <- data.frame(target=DDD.stat$target,
                            contrast=DDD.stat$contrast,
                            x=DDD.stat$deltaPS,
                            y=-log10(DDD.stat$adj.pval))
    data2plot$significance <- 'Not significant'
    data2plot$significance[DDD.stat$adj.pval< pval_cut & 
                             abs(DDD.stat$deltaPS)>=deltaPS_cut] <- 'Significant'
    q <- plotVolcano(data2plot = data2plot,
                     xlab = 'deltaPS of transcripts',ylab='-log10(FDR)',
                     title = title.idx[i],
                     col0 = col0,col1 = col1,size = size,top.n = top.n)
  }
  q
})
names(g) <- idx

######################################################################
##----->> save plot
lapply(names(g),function(i){
  message(i)
  png(paste0(figure.folder,'/',i,' volcano plot.png'),
      width = 10,height = 6,units = 'in',
      res = 150)
  print(g[[i]])
  dev.off()
  
  pdf(paste0(figure.folder,'/',i,' volcano plot.pdf'),
      width = 10,height = 6)
  print(g[[i]])
  dev.off()
})
```

#### Flow chart

``` r
################################################################################
##----->> DE vs DAS genes
DE.genes <- unique(DE_genes$target)
DAS.genes <- unique(DAS_genes$target)
genes.flow.chart <- function(){
  plotFlowChart(expressed =target_high$genes_high,
                x = DE.genes,
                y = DAS.genes,
                type = 'genes',
                pval.cutoff = pval_cut,lfc.cutoff = l2fc_cut,
                deltaPS.cutoff = deltaPS_cut)
}
genes.flow.chart()


png(filename = paste0(figure.folder,'/Union set DE genes vs DAS genes.png'),
    width = 22/2.54,height = 13/2.54,units = 'in',res = 300)
genes.flow.chart()
dev.off()

pdf(file = paste0(figure.folder,'/Union set DE genes vs DAS genes.pdf'),
    width = 22/2.54,height = 13/2.54)
genes.flow.chart()
dev.off()

################################################################################
##----->> DE vs DTU transcripts
DE.trans<- unique(DE_trans$target)
DTU.trans <- unique(DTU_trans$target)

trans.flow.chart <- function(){
  plotFlowChart(expressed = target_high$trans_high, 
                x = DE.trans, 
                y = DTU.trans, 
                type = 'transcripts', 
                pval.cutoff = pval_cut,
                lfc.cutoff = l2fc_cut,
                deltaPS.cutoff = deltaPS_cut)
}

trans.flow.chart()

png(filename = paste0(figure.folder,'/Union set DE transcripts vs DTU transcripts.png'),
    width = 22/2.54,height = 13/2.54,units = 'in',res = 300)
trans.flow.chart()
dev.off()

pdf(file = paste0(figure.folder,'/Union set DE transcripts vs DTU transcripts.pdf'),
    width = 22/2.54,height = 13/2.54)
trans.flow.chart()
dev.off()
```

#### Comparisons between contrast groups

``` r
################################################################################
##----->> DE genes
contrast2plot <- c("T2.Day1-T2.Day0","T2.Day4-T2.Day0")
targets <- lapply(contrast2plot,function(i){
  subset(DE_genes,contrast==i)$target
})
names(targets) <- contrast2plot
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DE genes', gp=gpar(cex=1.2)))

################################################################################
##----->> DAS genes
targets <- lapply(contrast2plot,function(i){
  subset(DAS_genes,contrast==i)$target
})
names(targets) <- contrast2plot
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DAS genes', gp=gpar(cex=1.2)))

################################################################################
##----->> DE transcripts
targets <- lapply(contrast2plot,function(i){
  subset(DE_trans,contrast==i)$target
})
names(targets) <- contrast2plot
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DE transcripts', gp=gpar(cex=1.2)))
    
################################################################################
##----->> DTU transcripts
targets <- lapply(contrast2plot,function(i){
  subset(DTU_trans,contrast==i)$target
})
names(targets) <- contrast2plot
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DTU transcripts', gp=gpar(cex=1.2)))
```

#### Comparison of transcription and AS targets

``` r
contrast.idx <- contrast[1]
################################################################################
##----->> DE vs DAS genes
x <- unlist(DEvsDAS_results[DEvsDAS_results$Contrast==contrast.idx,-1])
if(length(x)==0){
  message('No DE and/or DAS genes')
} else {
  names(x) <- c('DE','DE&DAS','DAS')
  g <- plotEulerDiagram(x = x,fill = gg.color.hue(2))
  g
  grid.arrange(g,top=textGrob('DE vs DAS genes', gp=gpar(cex=1.2)))
}

################################################################################
##----->> DE vs DTU transcripts
x <- unlist(DEvsDTU_results[DEvsDTU_results$Contrast==contrast.idx,-1])
if(length(x)==0){
  message('No DE and/or DTU transcripts')
} else {
  names(x) <- c('DE','DE&DTU','DTU')
  g <- plotEulerDiagram(x = x,fill = gg.color.hue(2))
  g
  grid.arrange(g,top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
}
```

Functional plot
---------------

### Heatmap of 3D targets

The targets are firstly clustered into groups by using the <a href="https://cran.r-project.org/web/packages/fastcluster/index.html" target="_blank">fastcluster</a> R package, and then made into heatmap by using the <a href="https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html" target="_blank">ComplexHeatmap</a> R package

``` r
################################################################################
##----->> DE genes
targets <- unique(DE_genes$target)
data2heatmap <- txi_genes$abundance[targets,]
column_title <- paste0(length(targets),' DE genes')
data2plot <- rowmean(x = t(data2heatmap),
                     group = metatable$label,
                     reorder = F)
data2plot <- t(scale(data2plot))
hc.dist <- dist(data2plot,method = dist_method)
hc <- fastcluster::hclust(hc.dist,method = cluster_method)
clusters <- cutree(hc, k = cluster_number)
clusters <- reorderClusters(clusters = clusters,dat = data2plot)

### save the target list in each cluster to result folder
x <- split(names(clusters),clusters)
x <- lapply(names(x),function(i){
  data.frame(Clusters=i,Targets=x[[i]])
})
x <- do.call(rbind,x)
colnames(x) <- c('Clusters','Targets')

###############################

g <- Heatmap(as.matrix(data2plot), name = 'Z-scores', 
             cluster_rows = TRUE,
             clustering_method_rows=cluster_method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DE genes.png'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54,
    units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DE genes.pdf'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()

################################################################################
##----->> DAS genes
targets <- intersect(unique(DE_genes$target),unique(DAS_genes$target))
data2heatmap <- txi_genes$abundance[targets,]
column_title <- paste0(length(targets),' DE&DAS genes')
data2plot <- rowmean(x = t(data2heatmap),group = metatable$label,reorder = F)
data2plot <- t(scale(data2plot))
hc.dist <- dist(data2plot,method = dist_method)
hc <- fastcluster::hclust(hc.dist,method = cluster_method)
clusters <- cutree(hc, k = cluster_number)
clusters <- reorderClusters(clusters = clusters,dat = data2plot)


### save the target list in each cluster to result folder
x <- split(names(clusters),clusters)
x <- lapply(names(x),function(i){
  data.frame(Clusters=i,Targets=x[[i]])
})
x <- do.call(rbind,x)
colnames(x) <- c('Clusters','Targets')
write.csv(x,file=paste0(result.folder,'/Target in each cluster heatmap ', column_title,'.csv'),
          row.names = F)
###############################

g <- Heatmap(as.matrix(data2plot), name = 'Z-scores', 
             cluster_rows = TRUE,
             clustering_method_rows=cluster_method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DAS genes.png'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DAS genes.pdf'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()


################################################################################
##----->> DE trans
targets <- unique(DE_trans$target)
data2heatmap <- txi_trans$abundance[targets,]
column_title <- paste0(length(targets),' DE trans')
data2plot <- rowmean(x = t(data2heatmap),
                     group = metatable$label,
                     reorder = F)
data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot,method = dist_method)
hc <- fastcluster::hclust(hc.dist,method = cluster_method)
clusters <- cutree(hc, k = cluster_number)
clusters <- reorderClusters(clusters = clusters,dat = data2plot)

### save the target list in each cluster to result folder
x <- split(names(clusters),clusters)
x <- lapply(names(x),function(i){
  data.frame(Clusters=i,Targets=x[[i]])
})
x <- do.call(rbind,x)
colnames(x) <- c('Clusters','Targets')
write.csv(x,file=paste0(result.folder,'/Target in each cluster heatmap ', column_title,'.csv'),
          row.names = F)
###############################

g <- Heatmap(as.matrix(data2plot), name = 'Z-scores', 
             cluster_rows = TRUE,
             clustering_method_rows=cluster_method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DE transcripts.png'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DE transcripts.pdf'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()


################################################################################
##----->> DTU trans
dist_method <- 'euclidean'
cluster_method <- 'ward.D'
cluster_number <- 10
targets <- intersect(unique(DE_trans$target),unique(DTU_trans$target))
data2heatmap <- txi_trans$abundance[targets,]
column_title <- paste0(length(targets),' DE&DTU trans')
data2plot <- rowmean(x = t(data2heatmap),
                     group = metatable$label,
                     reorder = F)
data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot,method = dist_method)
hc <- fastcluster::hclust(hc.dist,method = cluster_method)
clusters <- cutree(hc, k = cluster_number)
clusters <- reorderClusters(clusters = clusters,dat = data2plot)

### save the target list in each cluster to result folder
x <- split(names(clusters),clusters)
x <- lapply(names(x),function(i){
  data.frame(Clusters=i,Targets=x[[i]])
})
x <- do.call(rbind,x)
colnames(x) <- c('Clusters','Targets')
write.csv(x,file=paste0(result.folder,'/Target in each cluster heatmap ', column_title,'.csv'),
          row.names = F)
###############################

g <- Heatmap(as.matrix(data2plot), name = 'Z-scores', 
             cluster_rows = TRUE,
             clustering_method_rows=cluster_method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DTU transcripts.png'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DTU transcripts.pdf'),
    width = pmax(10,1*length(unique(metatable$label)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
```

### Profile plot

``` r
################################################################################
##----->> generate annotation of DE and/or DAS, DE and/or DTU
DE.genes <- unique(DE_genes$target)
DAS.genes <- unique(DAS_genes$target)
DE.trans <- unique(DE_trans$target)
DTU.trans <- unique(DTU_trans$target)

genes.ann <- set2(DE.genes,DAS.genes)
names(genes.ann) <- c('DEonly','DE&DAS','DASonly')
genes.ann <- plyr::ldply(genes.ann,cbind)[,c(2,1)]
colnames(genes.ann) <- c('target','annotation')


trans.ann <- set2(DE.trans,DTU.trans)
names(trans.ann) <- c('DEonly','DE&DTU','DTUonly')
trans.ann <- plyr::ldply(trans.ann,cbind)[,c(2,1)]
colnames(trans.ann) <- c('target','annotation')

################################################################################
##give a gene name
gene <- 'AT5G24470'
##----->> profile plot of TPM or read counts
g.pr <- plotAbundance(data.exp = txi_trans$abundance[target_high$trans_high,],
                       gene = gene,
                       mapping = mapping[target_high$trans_high,],
                       genes.ann = genes.ann,
                       trans.ann = trans.ann,
                       trans.expressed = NULL,
                       reps = metatable$label,
                       y.lab = 'TPM')
g.pr

################################################################################
##----->> PS plot
g.ps <- plotPS(data.exp = txi_trans$abundance[target_high$trans_high,],
                gene = gene,
                mapping = mapping[target_high$trans_high,],
                genes.ann = genes.ann,
                trans.ann = trans.ann,
                trans.expressed = NULL,
                reps = metatable$label,
                y.lab = 'PS')
g.ps
```

### GO annotation plot of DE/DAS genes

To make the plot, users must provide GO annotation table in csv file, with first column of "Category" (i.e. CC, BP and MF), second column of GO "Term" and the remaining columns of statistics of significance.

``` r
################################################################################
##----->> DE genes
go.table <- readr::read_csv(paste0(input.folder,
                                   '/DE genes GO annotation simplified.csv'))

g <- plotGO(go.table = go.table,col.idx = '(-)log10(FDR)',
            plot.title = 'GO annotation: DE genes')

### save to figure
png(paste0(figure.folder,'/DE genes GO annotation plot.png'),
    width = 20/2.54,height = 21/2.54,
    units = 'in',res = 300)
print(g)
dev.off()

pdf(paste0(figure.folder,'/DE genes GO annotation plot.pdf'),
    width = 20/2.54,height = 21/2.54)
print(g)
dev.off()
```

``` r
################################################################################
##----->> DAS genes
go.table <- readr::read_csv(paste0(input.folder,'/DAS genes GO annotation simplified.csv'))
g <- plotGO(go.table = go.table,col.idx = '(-)log10(FDR)',
         plot.title = 'GO annotation: DAS genes')

### save to figure
png(paste0(figure.folder,'/DAS genes GO annotation plot.png'),
    width = 20/2.54,height = 8/2.54,
    units = 'in',res = 300)
print(g)
dev.off()

pdf(paste0(figure.folder,'/DAS genes GO annotation plot.pdf'),
    width = 20/2.54,height = 8/2.54)
print(g)
dev.off()
```

Isoform switch analysis
-----------------------

### Generate switch scores

``` r
contrast <- contrast_pw
condition <- metatable$label
DAS_genes <- DAS_genes
TSIS.mapping.full <- mapping
rownames(TSIS.mapping.full) <- TSIS.mapping.full$TXNAME
TSIS.mapping.full <- TSIS.mapping.full[target_high$trans_high,]

###################################################################
##----->> isoform swtich analysis
if(TSISorisokTSP == 'isokTSP'){
  message('Peform isokTSP analysis ...')
  if(is.null(contrast) | is.null(condition) | is.null(metatable) | 
     is.null(DAS_genes) | is.null(TSIS.mapping.full) | 
     is.null(trans_TPM))
    return(NULL)
  
  scores.results <- list()
  for(i in seq_along(contrast)){
    contrast.idx <- contrast[i]
    metatable.idx <- unlist(strsplit(contrast.idx,'-'))
    if(any((metatable.idx %in% condition)==F)){
      message(paste0('IS analysis of contrast group: ',contrast.idx, ' is not applied'))
      next
    }
    
    message(paste0('IS analysis of contrast group: ',contrast.idx))
    col.idx <- which(condition %in% metatable.idx)
    DAS.genes <- unique(DAS_genes$target[DAS_genes$contrast==contrast.idx])
    if(length(DAS.genes)==0)
      next
    TSIS.mapping <- TSIS.mapping.full[TSIS.mapping.full$GENEID %in% DAS.genes,]
    TSIS.mapping <- TSIS.mapping[,c('GENEID','TXNAME')]
    TSIS.tpm <- trans_TPM[TSIS.mapping$TXNAME,col.idx]
    times <- metatable$label[col.idx]
    
    scores<-iso.switch.pairwise(data.exp=TSIS.tpm,
                                mapping = TSIS.mapping,
                                times=times,
                                rank=F,
                                min.difference = 1,
                                spline = F,
                                verbose = T,shiny = F)
    
    scores <- data.frame(contrast=contrast.idx,scores,row.names = NULL,check.names = F)
    scores.results <- c(scores.results,setNames(list(scores),contrast.idx))
  }
  scores.results <- do.call(rbind,scores.results)
  rownames(scores.results) <- NULL
  scores <- scores.results
} else {
  DAS.genes <- unique(DAS_genes$target)
  if(length(DAS.genes)==0)
    next
  TSIS.mapping <- TSIS.mapping.full[TSIS.mapping.full$GENEID %in% DAS.genes,]
  TSIS.mapping <- TSIS.mapping[,c('GENEID','TXNAME')]
  TSIS.tpm <- trans_TPM[TSIS.mapping$TXNAME,]
  times <- metatable$label
  scores.results <- iso.switch(data.exp = TSIS.tpm,
                                     mapping = TSIS.mapping,
                                     times=times,
                                     rank=F,
                                     min.t.points = 1,
                                     min.difference = 1,
                                     spline = (method_intersection == 'Spline'),
                                     spline.df = spline_df,
                                     verbose = T)
  # scores.results <- roundDF(scores,digits = 6)
  
  rownames(scores.results) <- NULL
  scores <- scores.results
}

###################################################################
##----->> filter
scores.filtered <-score.filter(
  scores = scores,
  prob.cutoff = TSIS_prob_cut,
  diff.cutoff = TSIS_diff_cut,
  t.points.cutoff = ifelse(TSISorisokTSP == 'isokTSP',1,TSIS_time_point_cut),
  pval.cutoff = 0.05,
  FDR.cutoff = TSIS_adj_pval_cut,
  cor.cutoff = TSIS_cor_cut,
  data.exp = NULL,
  mapping = NULL,
  sub.isoform.list = NULL,
  sub.isoform = F,
  max.ratio = F,
  x.value.limit = c(1,length(unique(metatable$label)))
)
DDD.data$scores_filtered <- scores.filtered
```

### Isoform switch plot

#### Switch numbers

``` r
x <- data.frame(table(factor(scores.filtered$contrast,levels = contrast_pw)))
g <- plotPWISnumber(x)
g

png(paste0(figure.folder,'/Isoform switch number.png'),
    width = 6,height = 4,
    units = 'in',res = 300)
print(g)
dev.off()
```

#### Isoform switch plot

``` r
iso1 <- 'AT5G65060_ID6'
iso2 <- 'AT5G65060_s1'
scores2plot <- scores.filtered
times0 <- metatable$label
data2plot <- trans_TPM[c(iso1,iso2),]
select.contrast <- contrast[1]

if(TSISorisokTSP == 'isokTSP'){
  if(!('contrast' %in% colnames(scores2plot))){
    message('Please select correct type of isoform switch.')
    return(NULL)
  }
  contrast.idx <- select.contrast
  contrast.idx <- unlist(strsplit(contrast.idx,'-'))[2:1]
  times.idx <- which(times0 %in% contrast.idx)
  times0 <- times0[times.idx]
  data2plot <- data2plot[,times.idx]
  scores2plot <- scores2plot[scores2plot$contrast==select.contrast,]
}

if(!is.numeric(times0)){
  times <- as.numeric(factor(times0,levels = unique(times0)))
  xlabs <- paste0(unique(times0),' (x=',unique(times),')')
} else {
  times <- times0
  xlabs <- unique(times0)
}

g <- plotTSIS(data2plot = data2plot,
              scores = scores2plot,
              iso1 = iso1,
              ribbon.plot = F,
              iso2 = iso2,
              y.lab = 'TPM',
              spline = F,
              spline.df = 9,
              times = times,
              errorbar.width = 0.03,
              x.lower.boundary = ifelse(is.numeric(times),min(times),1),
              x.upper.boundary = ifelse(is.numeric(times),max(times),
                                        length(unique(times))),
              show.scores = T,
              show.region = F)+
  scale_x_continuous(breaks = unique(times),labels = xlabs)
g
```

Generate report
---------------

### Summary parameters

``` r
params_list <- list()
params_list$condition_n = length(unique((metatable_new$label)))
params_list$brep_n = length(unique(metatable[,brep_col]))
params_list$srep_n = length(unique(metatable[,srep_col]))
params_list$samples_n = nrow(metatable_new)
params_list$has_srep = has_srep
params_list$quant_method = quant_method
params_list$tximport_method = tximport_method
params_list$cpm_cut = cpm_cut
params_list$cpm_samples_n = cpm_samples_n
params_list$norm_method = norm_method
params_list$has_batcheffect = has_batcheffect
params_list$RUVseq_method = RUVseq_method
params_list$contrast = contrast
params_list$DE_pipeline = DE_pipeline
params_list$pval_adj_method = pval_adj_method
params_list$pval_cut = pval_cut
params_list$l2fc_cut = l2fc_cut
params_list$deltaPS_cut = deltaPS_cut
params_list$DAS_pval_method = DAS_pval_method

##heatmap
params_list$dist_method <- dist_method
params_list$cluster_method <- cluster_method
params_list$cluster_number <- cluster_number

##TSIS
params_list$TSISorisokTSP <- TSISorisokTSP
params_list$TSIS_method_intersection <- method_intersection
params_list$TSIS_spline_df <- spline_df
params_list$TSIS_prob_cut <- TSIS_prob_cut
params_list$TSIS_diff_cut <- TSIS_diff_cut
params_list$TSIS_adj_pval_cut <- TSIS_adj_pval_cut
params_list$TSIS_time_point_cut <- ifelse(TSISorisokTSP == 'isokTSP',1,TSIS_time_point_cut)
params_list$TSIS_cor_cut <- TSIS_cor_cut

DDD.data$conditions <- metatable$label
DDD.data$params_list <- params_list
save(DDD.data,file=paste0(data.folder,'/DDD.data.RData'))
```

### Generate report

``` r
if(!file.exists('3D_report.Rmd'))
  file.copy(from=file.path(system.file("app", package = "ThreeDRNAseq"),
                           '3D_report.Rmd'), 
            to=getwd(), 
            overwrite = T, recursive = T, 
            copy.mode = T)

for (i in c('html_document','word_document','pdf_document')) {
  rmarkdown::render(input = '3D_report.Rmd',              
                    output_format = 'html_document',
                    output_dir = report.folder,
                    intermediates_dir = report.folder,
                    knit_root_dir = report.folder,
                    params = c(DDD.data=list(DDD.data))
  )
}
```

### Save significant results

``` r
####save results 
idx <- c('DE_genes','DAS_genes','DE_trans','DTU_trans','samples','contrast',
         'DDD_numbers','DEvsDAS_results','DEvsDTU_results','RNAseq_info')
idx.names <-gsub('_',' ',idx)
idx.names <- gsub('trans','transcripts',idx.names)
idx.names[1:4] <- paste0('Significant ',idx.names[1:4],' list and statistics')

idx <- c(idx,'scores','scores_filtered')
idx.names <- c(idx.names,'Raw isoform switch scores',
               'Significant isoform switch scores')
for(i in seq_along(idx)){
  if(is.null(DDD.data[[idx[i]]]))
    next
  write.csv(x = DDD.data[[idx[i]]],file = paste0(DDD.data$result.folder,'/',idx.names[i],'.csv'),row.names = F)
}
### save transcript-gene mapping
write.csv(x = DDD.data$mapping,
          file = paste0(DDD.data$result.folder,'/Transcript and gene mapping.csv'),
          row.names = F,na = '')

### save 3d list
threeD.list <-lapply(idx[1:4],function(i){
  unique(DDD.data[[i]]$target)
})

if(!any(sapply(threeD.list,is.null))){
  n <- max(sapply(threeD.list, length))
  threeD.list <- lapply(threeD.list, function(x){
    y <- rep(NA,n)
    y[1:length(x)] <- x
    y
  })
  names(threeD.list) <- idx[1:4]
  threeD <- do.call(cbind,threeD.list)
  write.csv(x = threeD,
            file = paste0(DDD.data$result.folder,'/DDD genes and transcript lists 
                          across all contrast groups.csv'),
            row.names = F,na = '')
}

##save all gene/transcript statistics
write.csv(x = DDD.data$genes_3D_stat$DE.stat,
          file = paste0(DDD.data$result.folder,'/DE gene testing statistics.csv'),
          row.names = F,na = '')

write.csv(x = DDD.data$trans_3D_stat$DE.stat,
          file = paste0(DDD.data$result.folder,
                        '/DE transcripts testing statistics.csv'),
          row.names = F,na = '')

if(DDD.data$params_list$DAS_pval_method=='F-test'){
  write.csv(x = DDD.data$trans_3D_stat$DAS.F.stat,
            file = paste0(DDD.data$result.folder,
                          '/DAS genes testing statistics.csv'),
            row.names = F,na = '')
}

if(DDD.data$params_list$DAS_pval_method=='Simes'){
  write.csv(x = DDD.data$trans_3D_stat$DAS.simes.stat,
            file = paste0(DDD.data$result.folder,
                          '/DAS genes testing statistics.csv'),
            row.names = F,na = '')
}

write.csv(x = DDD.data$trans_3D_stat$DTU.stat,
          file = paste0(DDD.data$result.folder,
                        '/DTU transcripts testing statistics.csv'),
          row.names = F,na = '')
```

Saved files in local directory
------------------------------

### Saved files in the "data" folder

The intermediate datasets in ".RData" format are saved in the "data" folder. Advanced R users can access these ".RData" objects by R command line. For the returned objects from R functions/packages, please go to corresponding documentations for details.

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

References
----------

Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525–527.

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
    ##  [1] compiler_3.5.1  magrittr_1.5    tools_3.5.1     htmltools_0.3.6
    ##  [5] yaml_2.2.0      Rcpp_1.0.0      stringi_1.2.4   rmarkdown_1.13 
    ##  [9] knitr_1.23      stringr_1.3.1   xfun_0.7        digest_0.6.18  
    ## [13] evaluate_0.13
