-   [Table of contents](#table-of-contents)
-   [Introduction](#introduction)
-   [Install and load R packages](#install-and-load-r-packages)
    -   [Install ThreeDRNAseq package](#install-threedrnaseq-package)
    -   [Install dependency packages](#install-dependency-packages)
    -   [Load R package](#load-r-package)
-   [Example data](#example-data)
-   [Set pipeline parameters](#set-pipeline-parameters)
-   [Data generation](#data-generation)
-   [Data pre-processing](#data-pre-processing)
    -   [Step 1: Merge sequencing replicates](#step-1-merge-sequencing-replicates)
    -   [Step 2: Filter low expression transcripts/genes](#step-2-filter-low-expression-transcriptsgenes)
    -   [Step 3: Principal component analysis (PCA)](#step-3-principal-component-analysis-pca)
    -   [Step 3-continued: Batch effect estimation](#step-3-continued-batch-effect-estimation)
    -   [Step 4: Data normalization](#step-4-data-normalization)
    -   [Data information](#data-information)
-   [DE, DAS and DTU analysis](#de-das-and-dtu-analysis)
    -   [Step 1: Load contrast group](#step-1-load-contrast-group)
    -   [Step 2: DE genes](#step-2-de-genes)
    -   [Step 3: DAS genes, DE and DTU transcripts](#step-3-das-genes-de-and-dtu-transcripts)
-   [Result summary](#result-summary)
    -   [Significant 3D targets and testing statistics](#significant-3d-targets-and-testing-statistics)
    -   [Significant 3D target numbers](#significant-3d-target-numbers)
    -   [Make plot](#make-plot)
-   [Advanced plots](#advanced-plots)
    -   [Heatmap of 3D targets](#heatmap-of-3d-targets)
    -   [Profile plot](#profile-plot)
    -   [GO annotation plot of DE/DAS genes](#go-annotation-plot-of-dedas-genes)
-   [Generate report](#generate-report)
-   [Saved files in local directory](#saved-files-in-local-directory)
    -   [Saved files in the "data" folder](#saved-files-in-the-data-folder)
    -   [Saved files in the "result" folder](#saved-files-in-the-result-folder)
    -   [Saved files in the "figure" folder](#saved-files-in-the-figure-folder)
    -   [Saved files in the "report" folder](#saved-files-in-the-report-folder)
-   [References](#references)
-   [Session information](#session-information)

Table of contents
-----------------

Introduction
------------

Instead of mouse click analysis in the <a href="xxx" target="_blank">3D RNA-seq App</a>, advanced R users also can use command line to do the same analysis. This file includes step-by-step coding, which corresponds to the steps in the 3D RNA-seq App. Users can refer to the details in the manual:xxx.

To use our pipeline in your work, please cite:

<a href="http://www.plantcell.org/content/30/7/1424" target="_blank">Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.</a>

Install and load R packages
---------------------------

### Install ThreeDRNAseq package

``` r
#######################################################################################################
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
cran.package.list <- c('shiny','shinydashboard','shinyFiles','plotly','eulerr',
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
</tbody>
</table>

### Load R package

Assume all the required R packages are installed.

``` r
library(ThreeDRNAseq)

library(shiny)
library(shinydashboard)
library(rhandsontable)
library(shinyFiles)
library(tximport)
library(edgeR)
library(limma)
library(plotly)
library(RUVSeq)
library(eulerr)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(Gmisc)

options(stringsAsFactors=F)
```

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

Set pipeline parameters
-----------------------

``` r
################################################################################
### Set folder names to save results
##----->> save results to local folders
data.folder <- 'data' # for .RData format results
result.folder <- 'result' # for 3D analysis results in csv files
figure.folder <- 'figure' # for figures
report.folder <- 'report'
### Set the input data folder
##----->> input data folder of samples.csv, mapping.csv, contrast.csv, etc.
input.folder <- 'examples/Arabidopsis_cold'

################################################################################
### parameters of tximport to generate read counts and TPMs
abundance.from <- 'salmon' # abundance generator
tximport.method <- 'lengthScaledTPM' # method to generate expression in tximport

################################################################################
### parameter for low expression filters
cpm.cut <- 1
sample.n.cut <- 1

################################################################################
### data normalisation parameter
norm.method <- 'TMM' ## norm.method is one of 'TMM','RLE' and 'upperquartile'

################################################################################
### parameters for 3D analysis
p.adjust.method <- 'BH'
pval.cutoff <- 0.01
lfc.cutoff <- 1
DE.pipeline <- 'limma'

################################################################################
### parameters for clusters in heatmap
dist.method <- 'euclidean'
cluster.method <- 'ward.D'
cluster.number <- 10
```

Data generation
---------------

Use the <a href="http://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank">tximport</a> R package (Soneson et al., 2016) to generate the transcript and gene level read counts and TPMs (transcript per million reads) from outputs of quantification tools, such as Salmon (Patro et al., 2017) and Kallisto (Bray et al., 2016).

``` r
################################################################################
##----->> Sample information, e.g. conditions, bio-reps, seq-reps, abundance paths, etc.
samples <- read.csv(paste0(input.folder,'/samples.csv'))
samples$condition <- gsub(pattern = '_','.',samples$condition)
samples$sample.name <- paste0(samples$condition,'_',samples$brep,'_',samples$srep)
save(samples,file=paste0(data.folder,'/samples.RData'))
data.files <-  as.character(samples$path) #Complete path of abundance, e.g. quant.sf files of salmon.

##----->> Transcript-gene mapping
mapping <- read.csv(paste0(input.folder,'/mapping.csv'))
rownames(mapping) <- mapping$TXNAME
save(mapping,file=paste0(data.folder,'/mapping.RData'))

################################################################################
##----->> Generate gene expression
##
txi_genes <- tximport(data.files,dropInfReps = T,
                      type = abundance.from, tx2gene = mapping,
                      countsFromAbundance = tximport.method)

## give colunames to the datasets
colnames(txi_genes$counts)<- samples$sample.name
colnames(txi_genes$abundance)<-samples$sample.name
colnames(txi_genes$length)<-samples$sample.name

## save the data
write.csv(txi_genes$counts,file=paste0(result.folder,'/counts_genes.csv'))
write.csv(txi_genes$abundance,file=paste0(result.folder,'/TPM_genes.csv'))
save(txi_genes,file=paste0(data.folder,'/txi_genes.RData'))

################################################################################
##----->> Generate transcripts expression
txi_trans<- tximport(data.files, type = abundance.from, tx2gene = NULL,
                     countsFromAbundance = tximport.method,
                     txOut = T,dropInfReps = T)

## give colunames to the datasets
colnames(txi_trans$counts)<- samples$sample.name
colnames(txi_trans$abundance)<-samples$sample.name
colnames(txi_trans$length)<-samples$sample.name

## save the data
write.csv(txi_trans$counts,file=paste0(result.folder,'/counts_trans.csv'))
write.csv(txi_trans$abundance,file=paste0(result.folder,'/TPM_trans.csv'))
save(txi_trans,file=paste0(data.folder,'/txi_trans.RData'))

################################################################################
genes_counts=txi_genes$counts
trans_counts=txi_trans$counts
```

Data pre-processing
-------------------

### Step 1: Merge sequencing replicates

If the RNA-seq data includes sequencing replicates (i.e. the same biological replicate/library is sequenced multiple times), the sequencing replicates are merged to increase the library depths, since they do not provided useful information for biological/technical variability.

``` r
##If no sequencing replicates, genes_counts and trans_counts remain the same by 
##runing the command
idx <- paste0(samples$condition,'_',samples$brep)
genes_counts <- sumarrays(genes_counts,group = idx)
trans_counts <- sumarrays(trans_counts,group = idx)

## update the samples
samples_new <- samples[samples$srep==samples$srep[1],]
save(samples_new,file=paste0(data.folder,'/samples_new.RData'))
```

### Step 2: Filter low expression transcripts/genes

To filter the low expressed transcripts across the samples, the expression mean-variance trend is used to determine the CPM (count per million reads) cut-off for low expression and the number of samples to cut. A gene is expressed if any of its transcript is expressed. Please see the details in the 3D RNA-seq App user manual.

``` r
################################################################################
##----->> Do the filters
target_high <- low.expression.filter(abundance = trans_counts,
                                     mapping = mapping,
                                     abundance.cut = cpm.cut,
                                     sample.n = sample.n.cut,
                                     unit = 'counts',
                                     Log=F)
save(target_high,file=paste0(data.folder,'/target_high.RData'))

################################################################################
##----->> Mean-variance plot
## transcript level
condition <- paste0(samples_new$condition)
counts.raw = trans_counts[rowSums(trans_counts>0)>0,]
counts.filtered = trans_counts[target_high$trans_high,]
mv.trans <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = condition)
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
condition <- paste0(samples_new$condition)
counts.raw = genes_counts[rowSums(genes_counts>0)>0,]
counts.filtered = genes_counts[target_high$genes_high,]
mv.genes <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = condition)
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

The PCA plot is a rough visualization of expression variation across conditions/samples and to determine whether the data has batch effects between different biological replicates.

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
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- samples_new$brep ## colour on biological replicates
# groups <- samples_new$condition ## colour on condtions
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
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- samples_new$brep
data2pca.ave <- rowmean(data2pca,samples_new$condition,reorder = F)
groups <- unique(samples_new$condition)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'Transcript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA Average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA Average expression.pdf'),
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
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- samples_new$brep ## colour on biological replicates
# groups <- samples_new$condition ## colour on condtions
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
groups <- samples_new$brep
data2pca.ave <- rowmean(data2pca,samples_new$condition,reorder = F)
groups <- unique(samples_new$condition)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'genescript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA Average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA Average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()
```

### Step 3-continued: Batch effect estimation

Only if see clear batch effects on the PCA plots, otherwise skip this step. The <a href="https://bioconductor.org/packages/release/bioc/html/RUVSeq.html" target="_blank">RUVSeq</a> (Risso et al., 2014) is used to estimate the batch effects. It will generate modified read counts, where the batch effects have been removed, for batch-free PCA plots and batch effect term which can be passed to the linear regression model for 3D analysis.

``` r
ruvseq.method <- 'RUVr' # ruvseq.method is one of 'RUVr', 'RUVs' and 'RUVg'
design <- condition2design(condition = samples_new$condition,
                           batch.effect = NULL)


################################################################################
##----->> trans level
trans_batch <- remove.batch(read.counts = trans_counts[target_high$trans_high,],
                            condition = samples_new$condition,
                            design = design,
                            contrast=NULL,
                            group = samples_new$condition,
                            method = ruvseq.method)
save(trans_batch,file=paste0(data.folder,'/trans_batch.RData')) 

################################################################################
##----->> genes level
genes_batch <- remove.batch(read.counts = genes_counts[target_high$genes_high,],
                            condition = samples_new$condition,
                            design = design,
                            contrast=NULL,
                            group = samples_new$condition,
                            method = ruvseq.method)
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
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- samples_new$brep ## colour on biological replicates
# groups <- samples_new$condition ## colour on condtions
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
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- samples_new$brep
data2pca.ave <- rowmean(data2pca,samples_new$condition,reorder = F)
groups <- unique(samples_new$condition)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'Transcript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA Average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA Average expression.pdf'),
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
groups <- samples_new$brep ## colour on biological replicates
# groups <- samples_new$condition ## colour on condtions
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
groups <- samples_new$brep
data2pca.ave <- rowmean(data2pca,samples_new$condition,reorder = F)
groups <- unique(samples_new$condition)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                  groups = groups,plot.title = 'genescript PCA: average expression',
                  ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA Average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA Average expression.pdf'),
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
               group = samples_new$condition,
               genes = mapping[target_high$trans_high,])
trans_dge <- suppressWarnings(calcNormFactors(dge,method = norm.method))
save(trans_dge,file=paste0(data.folder,'/trans_dge.RData'))

################################################################################
##----->> genes level
dge <- DGEList(counts=genes_counts[target_high$genes_high,],
               group = samples_new$condition)
genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm.method))
save(genes_dge,file=paste0(data.folder,'/genes_dge.RData'))

################################################################################
##----->> distribution plot
sample.name <- paste0(samples_new$condition,'.',samples_new$brep)
condition <- samples_new$condition

###--- trans level
data.before <- trans_counts[target_high$trans_high,]
data.after <- counts2CPM(obj = trans_dge,Log = T)
g <- boxplot.normalised(data.before = data.before,
                        data.after = data.after,
                        condition = condition,
                        sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Transcript data distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript data distribution.pdf'),
    width = 20/2.54,height = 20/2.54)
do.call(grid.arrange,g)
dev.off()

###--- genes level
data.before <- genes_counts[target_high$genes_high,]
data.after <- counts2CPM(obj = genes_dge,Log = T)
g <- boxplot.normalised(data.before = data.before,
                        data.after = data.after,
                        condition = condition,
                        sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Gene data distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene data distribution.pdf'),
    width = 20/2.54,height = 20/2.54)
do.call(grid.arrange,g)
dev.off()
```

### Data information

``` r
data.info <- data.frame(
  raw.number = c(length(mapping$TXNAME),length(unique(mapping$GENEID))),
  sample = rep(nrow(samples),2),
  condition = rep(length(unique(samples$condition)),2),
  brep = rep(length(unique(samples$brep)),2),
  srep = rep(length(unique(samples$srep)),2),
  merged.sample = c(ncol(trans_counts),ncol(genes_counts)),
  cpm.cut = c(cpm.cut,'--'),
  sample.n.cut = c(sample.n.cut,'--'),
  after.filter = c(length(target_high$trans_high),length(target_high$genes_high))
)
rownames(data.info) <- c('Transcript','Gene')
data.info
write.csv(data.info,file=paste0(result.folder,'/data.info.csv'),row.names = T)
```

DE, DAS and DTU analysis
------------------------

Three pipelines, <a href="https://bioconductor.org/packages/release/bioc/html/limma.html" target="_blank">limma</a> (Ritchie et al., 2015), <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">glmQL (edgeR)</a> and <a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">glm (edgeR)</a> (Robinson et al., 2010) are provided for DE gene and transcripts analysis. "limma" and "edgeR" have comparable performance and they have stringent controls of false positives than the "glm".

### Step 1: Load contrast group

``` r
# The contrast groups passed to the analysis are in a vector object, e.g. \code{contrast <- c('B-A','C-A')} means compare conditions "B" and "C" to condition "A".
contrast <- read.csv(paste0(input.folder,'/contrast.csv'))
contrast <- paste0(contrast[,1],'-',contrast[,2])
###---generate proper contrast format
contrast <- gsub('_','.',contrast)
save(contrast,file=paste0(data.folder,'/contrast.RData'))
```

### Step 2: DE genes

``` r
batch.effect <- genes_batch$W
# batch.effect <- NULL ## if has no batch effects
design <- condition2design(condition = samples_new$condition,
                           batch.effect = batch.effect)

################################################################################
if(DE.pipeline == 'limma'){
##----->> limma pipeline
genes_3D_stat <- limma.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 adjust.method = p.adjust.method)
}

if(DE.pipeline == 'glmQL'){
##----->> edgeR glmQL pipeline
genes_3D_stat <- edgeR.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 method = 'glmQL',
                                 adjust.method = p.adjust.method)
}

if(DE.pipeline == 'glm'){
##----->> edgeR glm pipeline
genes_3D_stat <- edgeR.pipeline(dge = genes_dge,
                                 design = design,
                                 deltaPS = NULL,
                                 contrast = contrast,
                                 diffAS = F,
                                 method = 'glm',
                                 adjust.method = p.adjust.method)
}
## save results
save(genes_3D_stat,file=paste0(data.folder,'/genes_3D_stat.RData'))
```

### Step 3: DAS genes, DE and DTU transcripts

``` r
################################################################################
##----->> generate deltaPS
transAbundance <- txi_trans$abundance[target_high$trans_high,]
deltaPS <- transAbundance2PS(transAbundance = transAbundance,
                             PS = NULL,
                             contrast = contrast,
                             condition = samples$condition,
                             mapping = mapping[target_high$trans_high,])
PS <- deltaPS$PS
deltaPS <- deltaPS$deltaPS
## save 
save(deltaPS,file=paste0(data.folder,'/deltaPS.RData')) 
save(PS,file=paste0(data.folder,'/PS.RData'))

################################################################################
##----->> DAS genes,DE and DTU transcripts
batch.effect <- genes_batch$W
# batch.effect <- NULL ## if has no batch effects
design <- condition2design(condition = samples_new$condition,
                           batch.effect = batch.effect)

################################################################################
p.adjust.method <- 'BH'
pval.cutoff <- 0.01
lfc.cutoff <- 1
deltaPS.cutoff <- 0.1

if(DE.pipeline == 'limma'){
##----->> limma pipeline
trans_3D_stat <- limma.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 adjust.method = p.adjust.method)
}

if(DE.pipeline == 'glmQL'){
##----->> edgeR glmQL pipeline
trans_3D_stat <- edgeR.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 method = 'glmQL',
                                 adjust.method = p.adjust.method)
}

if(DE.pipeline == 'glm'){
##----->> edgeR glm pipeline
trans_3D_stat <- edgeR.pipeline(dge = trans_dge,
                                 design = design,
                                 deltaPS = deltaPS,
                                 contrast = contrast,
                                 diffAS = T,
                                 method = 'glm',
                                 adjust.method = p.adjust.method)
}
## save results
save(trans_3D_stat,file=paste0(data.folder,'/trans_3D_stat.RData'))
```

Result summary
--------------

### Significant 3D targets and testing statistics

``` r
################################################################################
##----->> Summary DE genes
DE_genes <- summaryDEtarget(stat = genes_3D_stat$DE.stat,
                              cutoff = c(adj.pval=pval.cutoff,
                                         log2FC=lfc.cutoff))

################################################################################
## summary DAS genes, DE and DTU trans
##----->> DE trans
DE_trans <- summaryDEtarget(stat = trans_3D_stat$DE.stat,
                              cutoff = c(adj.pval=pval.cutoff,
                                         log2FC=lfc.cutoff))

##----->> DAS genes
DAS.p.method <- 'F-test' ## DAS.p.method is one of 'F-test' and 
if(DAS.p.method=='F-test') {
  DAS.stat <- trans_3D_stat$DAS.F.stat
} else {
  DAS.stat <- trans_3D_stat$DAS.Simes.stat
}

lfc <- genes_3D_stat$DE.lfc
lfc <- reshape2::melt(as.matrix(lfc))
colnames(lfc) <- c('target','contrast','log2FC')
DAS_genes <- summaryDAStarget(stat = DAS.stat,
                                lfc = lfc,
                                cutoff=c(pval.cutoff,deltaPS.cutoff))

##----->> DTU trans
lfc <- trans_3D_stat$DE.lfc
lfc <- reshape2::melt(as.matrix(lfc))
colnames(lfc) <- c('target','contrast','log2FC')
DTU_trans <- summaryDAStarget(stat = trans_3D_stat$DTU.stat,
                                lfc = lfc,cutoff = c(adj.pval=pval.cutoff,
                                                     deltaPS=deltaPS.cutoff))

################################################################################
## save results
save(DE_genes,file=paste0(data.folder,'/DE_genes.RData'))
save(DE_trans,file=paste0(data.folder,'/DE_trans.RData')) 
save(DAS_genes,file=paste0(data.folder,'/DAS_genes.RData')) 
save(DTU_trans,file=paste0(data.folder,'/DTU_trans.RData')) 

## csv
write.csv(DE_genes,file=paste0(result.folder,'/DE genes.csv'),row.names = F)
write.csv(DAS_genes,file=paste0(result.folder,'/DAS genes.csv'),row.names = F)
write.csv(DE_trans,file=paste0(result.folder,'/DE transcripts.csv'),row.names = F)
write.csv(DTU_trans,file=paste0(result.folder,'/DTU transcripts.csv'),row.names = F)
```

### Significant 3D target numbers

``` r
################################################################################
##----->> target numbers
DDD.numbers <- summary3Dnumber(DE_genes = DE_genes,
                                 DAS_genes = DAS_genes,
                                 DE_trans = DE_trans,
                                 DTU_trans=DTU_trans,
                                 contrast = contrast)
DDD.numbers
write.csv(DDD.numbers,file=paste0(result.folder,'/DE DAS DTU numbers.csv'),
          row.names = F)

################################################################################
##----->> DE vs DAS
DEvsDAS.results <- DEvsDAS(DE_genes = DE_genes,
                           DAS_genes = DAS_genes,
                           contrast = contrast)
DEvsDAS.results
write.csv(DEvsDAS.results,file=paste0(result.folder,'/DE vs DAS gene number.csv'),
          row.names = F)


################################################################################
##----->> DE vs DTU
DEvsDTU.results <- DEvsDTU(DE_trans = DE_trans,
                           DTU_trans = DTU_trans,
                           contrast = contrast)
DEvsDTU.results
write.csv(DEvsDTU.results,file=paste0(result.folder,'/DE vs DTU transcript number.csv'),row.names = F)
```

### Make plot

#### Up- and down-regulation

``` r
################################################################################
##----->> DE genes
idx <- factor(DE_genes$contrast,levels = contrast)
targets <-  split(DE_genes,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
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
png(paste0(figure.folder,'/DE genes updown regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

pdf(paste0(figure.folder,'/DE genes updown regulation numbers.pdf'),
    width = length(contrast)*5/2.54,10/2.54)
print(g.updown)
dev.off()


################################################################################
##----->> DAS genes
idx <- factor(DAS_genes$contrast,levels = contrast)
targets <-  split(DAS_genes,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
  } else {
    x <- data.frame(contrast=i,table(targets[[i]]$up.down))
    colnames(x) <- c('contrast','regulation','number')
  }
  x
})
data2plot <- do.call(rbind,data2plot)
g.updown <- plotUpdown(data2plot,plot.title = 'DAS genes',contrast = contrast)
print(g.updown)

### save to figure
png(paste0(figure.folder,'/DAS genes updown regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

pdf(paste0(figure.folder,'/DAS genes updown regulation numbers.pdf'),
    width = length(contrast)*5/2.54,10/2.54)
print(g.updown)
dev.off()


################################################################################
##----->> DE trans
idx <- factor(DE_trans$contrast,levels = contrast)
targets <-  split(DE_trans,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
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
png(paste0(figure.folder,'/DE transcripts updown regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

pdf(paste0(figure.folder,'/DE transcripts updown regulation numbers.pdf'),
    width = length(contrast)*5/2.54,10/2.54)
print(g.updown)
dev.off()

################################################################################
##----->> DTU trans
idx <- factor(DTU_trans$contrast,levels = contrast)
targets <-  split(DTU_trans,idx)
data2plot <- lapply(contrast,function(i){
  if(nrow(targets[[i]])==0){
    x <- data.frame(contrast=i,regulation=c('down_regulate','up_regulate'),number=0)
  } else {
    x <- data.frame(contrast=i,table(targets[[i]]$up.down))
    colnames(x) <- c('contrast','regulation','number')
  }
  x
})
data2plot <- do.call(rbind,data2plot)
g.updown <- plotUpdown(data2plot,plot.title = 'DTU trans',contrast = contrast)
print(g.updown)

### save to figure
png(paste0(figure.folder,'/DTU transcripts updown regulation numbers.png'),
    width = length(contrast)*5/2.54,10/2.54,units = 'in',res = 300)
print(g.updown)
dev.off()

pdf(paste0(figure.folder,'/DTU transcripts updown regulation numbers.pdf'),
    width = length(contrast)*5/2.54,10/2.54)
print(g.updown)
dev.off()
```

#### Flow chart

``` r
################################################################################
##----->> DE vs DAS genes
DE.genes <- unique(DE_genes$target) 
DAS.genes <- unique(DAS_genes$target) 

genes.flow.chart <- function(){
  plotFlowChart(expressed = target_high$genes_high, 
                x = DE.genes, 
                y = DAS.genes, 
                type = 'genes', 
                pval.cutoff = pval.cutoff,
                lfc.cutoff = lfc.cutoff,
                deltaPS.cutoff = deltaPS.cutoff)
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
                pval.cutoff = pval.cutoff,
                lfc.cutoff = lfc.cutoff,
                deltaPS.cutoff = deltaPS.cutoff)
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
targets <- lapply(contrast,function(i){
  subset(DE_genes,contrast==i)$target
})
names(targets) <- contrast
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DE genes', gp=gpar(cex=1.2)))

################################################################################
##----->> DAS genes
targets <- lapply(contrast,function(i){
  subset(DAS_genes,contrast==i)$target
})
names(targets) <- contrast
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DAS genes', gp=gpar(cex=1.2)))

################################################################################
##----->> DE transcripts
targets <- lapply(contrast,function(i){
  subset(DE_trans,contrast==i)$target
})
names(targets) <- contrast
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DE transcripts', gp=gpar(cex=1.2)))
    
################################################################################
##----->> DTU transcripts
targets <- lapply(contrast,function(i){
  subset(DTU_trans,contrast==i)$target
})
names(targets) <- contrast
g <- plotEulerDiagram(x = targets)
grid.arrange(g,top=textGrob('DTU transcripts', gp=gpar(cex=1.2)))
```

#### Comparison of transcription and AS targets

``` r
contrast.idx <- contrast[1]
################################################################################
##----->> DE vs DAS genes
x <- unlist(DEvsDAS.results[DEvsDAS.results$Contrast==contrast.idx,-1])
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
x <- unlist(DEvsDTU.results[DEvsDTU.results$Contrast==contrast.idx,-1])
if(length(x)==0){
  message('No DE and/or DTU transcripts')
} else {
  names(x) <- c('DE','DE&DTU','DTU')
  g <- plotEulerDiagram(x = x,fill = gg.color.hue(2))
  g
  grid.arrange(g,top=textGrob('DE vs DTU transcripts', gp=gpar(cex=1.2)))
}
```

Advanced plots
--------------

### Heatmap of 3D targets

The targets are firstly clustered into groups by using the <a href="https://cran.r-project.org/web/packages/fastcluster/index.html" target="_blank">fastcluster</a> R package, and then made into heatmap by using the <a href="https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html" target="_blank">ComplexHeatmap</a> R package

``` r
################################################################################
##----->> DE genes
targets <- unique(DE_genes$target)
data2heatmap <- txi_genes$abundance[targets,]
column_title <- paste0(length(targets),' DE genes')
data2plot <- rowmean(x = t(data2heatmap),
                     group = samples$condition,
                     reorder = F)
data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot,method = dist.method)
hc <- fastcluster::hclust(hc.dist,method = cluster.method)
clusters <- cutree(hc, k = cluster.number)
clusters <- reorder.clusters(clusters = clusters,dat = data2plot)

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
             clustering_method_rows=cluster.method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DE genes.png'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DE genes.pdf'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()

################################################################################
##----->> DAS genes
dist.method <- 'euclidean'
cluster.method <- 'ward.D'
cluster.number <- 10
targets <- intersect(unique(DE_genes$target),unique(DAS_genes$target))
data2heatmap <- txi_genes$abundance[targets,]
column_title <- paste0(length(targets),' DE&DAS genes')
data2plot <- rowmean(x = t(data2heatmap),
                     group = samples$condition,
                     reorder = F)
data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot,method = dist.method)
hc <- fastcluster::hclust(hc.dist,method = cluster.method)
clusters <- cutree(hc, k = cluster.number)
clusters <- reorder.clusters(clusters = clusters,dat = data2plot)

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
             clustering_method_rows=cluster.method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DAS genes.png'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DAS genes.pdf'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()


################################################################################
##----->> DE trans
dist.method <- 'euclidean'
cluster.method <- 'ward.D'
cluster.number <- 10
targets <- unique(DE_trans$target)
data2heatmap <- txi_trans$abundance[targets,]
column_title <- paste0(length(targets),' DE trans')
data2plot <- rowmean(x = t(data2heatmap),
                     group = samples$condition,
                     reorder = F)
data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot,method = dist.method)
hc <- fastcluster::hclust(hc.dist,method = cluster.method)
clusters <- cutree(hc, k = cluster.number)
clusters <- reorder.clusters(clusters = clusters,dat = data2plot)

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
             clustering_method_rows=cluster.method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DE transcripts.png'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DE transcripts.pdf'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()


################################################################################
##----->> DTU trans
dist.method <- 'euclidean'
cluster.method <- 'ward.D'
cluster.number <- 10
targets <- intersect(unique(DE_trans$target),unique(DTU_trans$target))
data2heatmap <- txi_trans$abundance[targets,]
column_title <- paste0(length(targets),' DE&DTU trans')
data2plot <- rowmean(x = t(data2heatmap),
                     group = samples$condition,
                     reorder = F)
data2plot <- t(scale(data2plot))

hc.dist <- dist(data2plot,method = dist.method)
hc <- fastcluster::hclust(hc.dist,method = cluster.method)
clusters <- cutree(hc, k = cluster.number)
clusters <- reorder.clusters(clusters = clusters,dat = data2plot)

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
             clustering_method_rows=cluster.method,
             row_dend_reorder = T,
             show_row_names = FALSE, 
             show_column_names = ifelse(ncol(data2plot)>10,F,T),
             cluster_columns = FALSE,
             split=clusters,
             column_title= column_title)

draw(g,column_title='Conditions',column_title_side = "bottom")

### save to figure
png(paste0(figure.folder,'/Heatmap DTU transcripts.png'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54,units = 'in',res = 300)
draw(g,column_title='Conditions',column_title_side = "bottom")
dev.off()
pdf(paste0(figure.folder,'/Heatmap DTU transcripts.pdf'),
    width = pmax(10,1*length(unique(samples$condition)))/2.54,height = 20/2.54)
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
                       reps = samples$condition,
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
                reps = samples$condition,
                y.lab = 'PS')
g.ps
```

### GO annotation plot of DE/DAS genes

To make the plot, users must provide GO annotation table in csv file, with first column of "Category" (i.e. CC, BP and MF), second column of GO "Term" and the remaining columns of statistics of significance.

``` r
################################################################################
##----->> DE genes
go.table <- suppressMessages(readr::read_csv(paste0(input.folder,'/DE genes GO terms.csv')))
go.table <- go.table[,c(1,2,5)]
g <- plotGO(go.table = go.table,col.idx = '-log10(FDR)',
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
go.table <- suppressMessages(readr::read_csv(paste0(input.folder,'/DAS genes GO terms.csv')))
go.table <- go.table[,c(1,2,5)]
g <- plotGO(go.table = go.table,col.idx = '-log10(FDR)',
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

Generate report
---------------

``` r
################################################################################
##----->> Summary parameters 
para <- rbind(
  c('Folders','Directory',getwd()),
  c('','Data',paste0(getwd(),'/',data.folder)),
  c('','Results',paste0(getwd(),'/',result.folder)),
  c('','Figures',paste0(getwd(),'/',figure.folder)),
   c('','Reports',paste0(getwd(),'/',report.folder)),
  c('Data generation','tximport method',tximport.method),
  c('Data pre-processing','Sequencing replicates',length(unique(samples$srep))),
  c('','Sequencing replicate merged',ifelse(nrow(samples)>nrow(samples_new), 'Yes','No')),
  c('','Low expression CPM cut-off',cpm.cut),
  c('','Sample number for CPM cut-off', sample.n.cut),
  c('','Batch effect estimation',
    ifelse(length(list.files(data.folder,'*batch.RData'))>0,'Yes','No')),
  c('','Batch effect estimation method',ruvseq.method),
  c('','Normalization method',norm.method),
  c('DE DAS and DTU','Pipeline',DE.pipeline),
  c('','AS function',ifelse(DE.pipeline=='limma',
                            'limma::diffSplice','edgeR::diffSpliceDGE')),
  c('','P-value adjust method',p.adjust.method),
  c('','Adjusted p-value cut-off',pval.cutoff),
  c('','Log2 fold change cut-off',lfc.cutoff),
  c('','delta PS cut-off',deltaPS.cutoff),
  c('Heatmap','Distance method',dist.method),
  c('','Cluster method',cluster.method),
  c('','Number of clusters',cluster.number)
)
para <- data.frame(para)
colnames(para) <- c('Step','Description','Parameter (Double click to edit if not correct)')
para
write.csv(para,file=paste0(result.folder,'/Parameter summary.csv'),row.names = F)

################################################################################
#----->> download report.Rmd file from Github
report.url <- 'https://raw.githubusercontent.com/wyguo/ThreeDRNAseq/master/vignettes/report.Rmd'
if(!file.exists('report.Rmd'))
  download.file(url = report.url,destfile = 'report.Rmd',method = 'auto')
generate.report(input.Rmd = 'report.Rmd',report.folder = report.folder,type = 'all')
```

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
    ##  [1] compiler_3.5.1  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
    ##  [5] tools_3.5.1     htmltools_0.3.6 yaml_2.2.0      Rcpp_0.12.18   
    ##  [9] stringi_1.1.7   rmarkdown_1.10  knitr_1.20      stringr_1.3.1  
    ## [13] digest_0.6.18   evaluate_0.11
