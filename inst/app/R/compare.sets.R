#' Generate all possible logical relations between set x and y.
#' @param x,y,z vectors of set x, set y and set z.
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

#############################################################################################
##-----------------------------------------------------------------------------------------
#' @rdname set2
#' @export

set3 <- function(x,y,z){
  a5 <- Reduce(intersect, list(x,y,z))
  a2 <- setdiff(intersect(x,y),a5)
  a4 <- setdiff(intersect(x,z),a5)
  a6 <- setdiff(intersect(y,z),a5)
  a1 <- setdiff(x,c(a2,a4,a5))
  a3 <- setdiff(y,c(a2,a5,a6))
  a7 <- setdiff(z,c(a4,a5,a6))
  results <- list(a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6,a7=a7)
  # attributes(results) <- list(
  #   a1='x-(y&z)',
  #   a2='(x&y)-z',
  #   a3='y-(x&z)',
  #   a4='(x&z)-y',
  #   a5='x&y&z',
  #   a6='(y&z)-x',
  #   a7='z-(x&y)'
  # )
  name.idx <- c('x-(y&z)','(x&y)-z','y-(x&z)','(x&z)-y','x&y&z','(y&z)-x','z-(x&y)')
  name.idx <- paste0(paste0('a',1:7,':'),name.idx)
  names(results) <- name.idx
  return(results)
}

