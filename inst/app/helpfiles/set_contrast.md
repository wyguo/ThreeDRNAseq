---
output:
  html_document: default
  word_document: default
---
<div align="justify">


3D RNA-seq App provides a flexible way where users can generate contrast groups of simple pair-wise analyses and complex experimental design such as time-series, development series and multiple conditions from grouping labels. For example, in Supplemental Figure 4B, contrast groups can be set as WT.B-WT.A (WT.B vs WT.A) and MU.B-MU.A (MU.B vs MU.A) to compare expression of condition WT.B to WT.A and MU.B to MU.A, respectively. Contrast group settings of pair-wise conditions suit to most experimental design in RNA-seq studies. In addition, users can set contrast group (MU.B+WT.B)/2-(MU.A+WT.A)/2 to compare the mean of multiple conditions or contrast group (MU.B-WT.B)-(MU.A-WT.A) to compare the difference of two $L_2Fs$ of pair-wise conditions (Figure B).

<!--![](select_factor.png)-->
<img src="select_factor.png" style="width: 80%; display: block; margin-left: auto; margin-right: auto;">

**Figure**: Generate group labels of factor levels of experimental design and set contrast groups for expression comparisons.  Meta-table of single factor (A) and multiple factors (B).
