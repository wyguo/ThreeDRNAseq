#' Generate all possible logical relations between set x and y.
#' @param x,y vectors of set x and set y.
#' @return a list of three elements: "x.only", "xy" and "y.only".
#' @examples 
#' x <- letters[1:10]
#' y <- letters[5:15]
#' z <- set2(x,y)
#' combo <- c(x=length(z$x.only),y=length(z$y.only),"x&y"=length(z$xy))
#' plot.euler.diagram(combo)
#' 
#' @export
#' @seealso \code{\link{plot.euler.diagram}} and \code{\link{eulerr::euler}}
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


#' Plot Euler diagram
#' @param x 1) a list of individuals in different sets for comparisons or 2) a vecotr of numbers of set relations.
#' @param fill a vector of colours to fill the diagram. 
#' @param shape geometric shape used in the diagram.
#' @param ... arguments passed down to the \code{\link{eulerr::euler}} function.
#' @return a Euler diagram
#' @examples 
#' x <- letters[1:10]
#' y <- letters[5:15]
#' z <- set2(x,y)
#' combo <- c(x=length(z$x.only),y=length(z$y.only),"x&y"=length(z$xy))
#' plot.euler.diagram(list(x=x,y=y))
#' plot.euler.diagram(combo)
#' @export
#' @seealso \code{\link{eulerr::euler}}.
#' 
plot.euler.diagram <- function(x,
                               fill = gg.color.hue(length(x)),
                               shape = c("ellipse", "circle"),...){
  shape <- match.arg(shape,c("ellipse", "circle"))
  fit <- euler(x,...)
  # grid.newpage()
  g <- plot(fit,quantities = TRUE,
            labels = list(font =1),
            fill=fill,shape = shape)
  g
}