#' Round the numeric columns of a data frame
#' @param df a data fram object
#' @param digits an integer indicating the number of decimal places
#' @param ... arguments to be passed to \code{\link{round}}.
#' @seealso  \code{\link{round}}.
#' @export
#'
roundDF <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}