#' Example data of the xxx package
#' The example dataset includes 200 transcripts from 84 genes in Arabidopsis. TPM and read count expression 
#' were extracted from three time-points T2 (\eqn{20^oC}), T10 (Day 1 of \eqn{4^oC} transition) and T19 (Day 4 of \eqn{4^oC} acclimation)
#' from study of Arabidopsis in response to cold (Calixto et al., 2018).
#' 
#' @docType data
#' 
#' @usage  data(exp.data)
#' 
#' @format a list object with elements: trans_TPM (transcript level TPM data.frame), trans_counts (transcript level read count data.frame),
#' mapping (transcript-gene mapping information), samples (condition, biological and sequencing replicate information), 
#' genes.ann (DE and/or DAS gene information) and trans.ann (DE and/or DTU transcript information).
#' 
#' @examples 
#' data(exp.data)
#' lapply(exp.data, head)
#' 
#' @references Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. 
#' (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell, tpc.00177.2018.

# load('data/txi_genes.RData')
# load('data/txi_trans.RData')
# load('data/samples.RData')
# load('data/mapping.RData')
# load('data/target_high.RData')
# 
# genes <- c(intersect(DE.genes,DAS.genes)[1:40],setdiff(DE.genes,DAS.genes)[1:30],setdiff(DAS.genes,DE.genes)[1:30])
# mapping <- mapping[target_high$trans_high,]
# trans <- mapping$TXNAME[mapping$GENEID %in% genes]
# length(trans)
# mapping <- mapping[trans,]
# trans_counts <- txi_trans$counts[trans,]
# trans_TPM <- txi_trans$abundance[trans,]
# samples <- samples[,c('condition','brep','srep')]
# genes.ann <- genes.ann[genes,]
# trans.ann <- na.omit(trans.ann[genes,])
# genes.ann0 <- genes.ann
# trans.ann0 <- trans.ann
# exp.data <- list(trans_TPM=trans_TPM,trans_counts=trans_counts,mapping=mapping,samples=samples,genes.ann=genes.ann,trans.ann=trans.ann)
# save(exp.data,file='data/exp.data.RData')

# 
# trans <- target_high$trans_high[1:200]
# genes <- substr(trans,1,9)
# mapping <- data.frame(TXNAME=trans,GENEID=genes)
# genes <- unique(genes)
# trans_counts <- txi_trans$counts[trans,]
# genes_counts <- txi_genes$counts[genes,]
# trans_TPM <- txi_trans$abundance[trans,]
# genes_TPM <- txi_genes$abundance[genes,]
# samples <- samples[,c('condition','brep','srep')]

# exp.data <- list(trans_TPM=trans_TPM,trans_counts=trans_counts,mapping=mapping,samples=samples)
# save(exp.data,file='data/exp.data.RData')

