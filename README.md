---
title: "ThreeDRNAseq: a shiny App for differential expression and differential alternative splicing analysis"
subtitle: 'User manual'
author: "Wenbin Guo"
date: "2018-11-17"
output: 
  BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{Introduction to ThreeDRNAseq}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Description
This package provides an interactive graphical user interface (GUI) for 
RNA-seq data differential expression (DE), differential alternative splicing (DAS) 
and differential transcript usage (DTU) analyses based on two popular pipelines: 
limma and edgeR. The 3DRNAseq GUI is based on R shiny App and enables command-line-free 
analysis. To perform analysis, the first step is to generate transcript quantification 
from quantification tools, such as Salmon and Kallisto. Then users can do mouse click on 
the App to upload transcript read counts, perform DE and DAS analysis, and make beautiful plots, 
e.g. expression mean-variance trend plots, PCA plots, heatmap, GO annotation plots, etc. 
The GUI has steps of proper data pre-processing and false positive controls. These lead to 
robust DE DAS predictions. All the results and plots can be saved to a report in html, 
pdf or word format. The analysis pipeline has been successfully used in different RNA-seq studies 
from Arabidopsis, barley and potato. 

# Installation and loading
