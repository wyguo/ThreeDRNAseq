#' Fit a time-series with spline curve
#'
#' @param x a vector of time-series expression
#' @param times a numeric vector of time labellings of all the relicated samples, e.g. 1,1,1,2,2,2,3,3,3,...
#' @param nrep number of replicates.
#' @param df degree of freedom used in \code{\link{ns}} in \code{\link{splines}} package.
#' @param ... additional arguments passed to \code{\link{predict}}.
#'
#' @return predictions results returned from \code{\link{predict}}.
#'
#' @seealso \code{\link{lm}}, \code{\link{ns}}, \code{\link{predict}}
#'
#' @export
#'

ts.spline<-function(x,times,df=5,...){
  df<-min(df,length(times)-2)
  times<-as.numeric(times)
  x<-data.frame(times=times,value=as.numeric(x))
  fit<-lm(value~splines::ns(times,df=df),data=x)
  predict(fit,data.frame(times=unique(times)),...)
}

