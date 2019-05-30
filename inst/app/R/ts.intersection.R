#' Searching for intersection points of two time-series
#'
#' @param x1 a vector of values for first time-series
#' @param x2 a vector of values for second time-series
#'
#' @return a data frame of intersection points. The first and second columns are x and y coordinates values of intersections, respectively.
#'
#' @export
#'
ts.intersection<-function(x1,x2,times){
  y <- sort(unique(times),decreasing = F)
  d1 <- x1
  d2 <- x2
  poi <- which(diff(d1 > d2) != 0) 
  if(length(poi)==0)
    return(NULL)
  intersect.points <- y[poi]
  slope1 <- (d1[poi]-d1[poi+1])/(y[poi]-y[poi+1])
  slope2 <- (d2[poi]-d2[poi+1])/(y[poi]-y[poi+1])
  x.points <- intersect.points + ((d2[poi] - d1[poi])/(slope1 - slope2))
  y.points <- d1[poi] + (slope1 * (x.points - intersect.points))
  intersection.points <- data.frame(x.points=x.points,y.points=y.points,row.names = NULL)
  # dup.idx <- which(duplicated(intersection.points) | duplicated(intersection.points,fromLast = T))
  # if(length(dup.idx)>0)
  #   intersection.points<-intersection.points[-dup.idx,]

  intersection.points<-intersection.points[!duplicated(round(intersection.points,5)),]
  return(intersection.points)
}
