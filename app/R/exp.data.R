#' Example data of the ThreeDRNAseq package
#' The example dataset includes 200 transcripts and 84 genes, which a sub-dataset from RNA-seq study of Arabidopsis in response to 
#' cold (Calixto et al., 2018). Transcript level TPM and read count expression were extracted from three time-points T2 (\eqn{20^oC}), T10 (Day 1 of \eqn{4^oC} transition) 
#' and T19 (Day 4 of \eqn{4^oC} acclimation). The \code{exp.data} is a list object with elements:
#' \itemize{
#'   \item{\code{trans_TPM:} transcript level TPM.}
#'   \item{\code{trans_counts:} transcript level read counts.} 
#'   \item{\code{mapping:} transcript-gene mapping.}
#'   \item{\code{samples:} the sample information.}
#'   \item{\code{go.table:} GO annotation of DE genes.}
#' }
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
#' names(exp.data)
#' lapply(exp.data, head)
#' 
#' @references Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. 
#' (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell, tpc.00177.2018.
