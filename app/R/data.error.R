#' Calculate standard deviation or standard errors
#'
#' Calculate standard deviation or standard errors from a numerical vector.
#'
#' @param x a numerical vector
#' @param error.type method used to calculate data error. Options include "sd" for standard deviation and "stderr" for standard error.
#'
#' @return data error
#' @seealso \code{\link{sd}}
#'
#' @examples
#'
#' set.seed(2000)
#' x=rnorm(1000,mean=0,sd=1)
#' ##standard error
#' data.error(x,error.type = 'stderr')
#' ##standard deviation
#' data.error(x,error.type = 'sd')
#'
#' @export
#'

data.error <- function (x, error.type = c("stderr","sd")){
  error.type <- match.arg(error.type,c("stderr","sd"))
  if (error.type == "stderr") 
    error <- sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))
  if (error.type == "sd") 
    error <- sd(na.omit(x))
  return(error)
}
