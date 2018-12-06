#' Reorder the clusters according to their average exprssion 
#' @param clusters a vector of clusters of targets
#' @param data expression data of these target
#' @return a vector of ordered clusters
#' @export
reorderClusters <- function(clusters,dat){
  targets <- rownames(dat)
  names(clusters) <- targets
  idx0 <- rowmean(dat,group = clusters,reorder = T)
  para <- lapply(1:ncol(idx0),function(i) idx0[,i])
  idx0 <- idx0[do.call(order,para),]
  idx0 <- idx0[rev(rownames(idx0)),]
  order.idx <- rownames(idx0)
  idx <- 1:nrow(idx0)
  names(idx) <- order.idx
  idx <- idx[as.character(clusters)]
  names(idx) <- names(clusters)
  idx
}