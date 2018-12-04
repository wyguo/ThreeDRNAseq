#' Sum Over Replicate Arrays
#' 
#' @param x a matrix-like object.
#' @param group grouping of sample identifier.
#' @return A data object of the same class as x with a column for sums according grouping.
#' @examples 
#' set.seed(100)
#' x<- matrix(rnorm(8*4),8,4)
#' colnames(x) <- c("a","a","b","b")
#' sumarrays(x)
#' data.frame(a=x[,1]+x[,2],b=x[,3]+x[,4])
#' @export

sumarrays<-function(x,group=NULL){
  if(is.null(group))
    group<-colnames(x)
  colnames(x)<-group
  
  x<-rowsum(t(x),group = group)
  #order the columns as the input group odering
  x<-data.frame(t(x)[,unique(group)])
  colnames(x) <- unique(group)
  return(x)
}

