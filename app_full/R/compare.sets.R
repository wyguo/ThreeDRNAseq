#' Generate all possible logical relations between set x and y.
#' @param x,y vectors of set x and set y.
#' @return a list of three elements: "x.only", "xy" and "y.only".
#' @examples 
#' x <- letters[1:10]
#' y <- letters[5:15]
#' z <- set2(x,y)
#' combo <- c(x=length(z$x.only),y=length(z$y.only),"x&y"=length(z$xy))
#' plotEulerDiagram(combo)
#' 
#' @export
#' @seealso \code{\link{plotEulerDiagram}} and \code{\link{eulerr::euler}}
#' 
set2 <- function(x,y){
  x.only <- setdiff(x,y)
  xy <- intersect(x,y)
  y.only <- setdiff(y,x)
  results <- list(x.only=x.only,xy=xy,y.only=y.only)
  attributes(results) <- list(
    x.only='x-x&y',
    xy='x&y',
    y.only='y-x&y'
  )
  names(results) <- c('x.only','xy','y.only')
  return(results)
}


