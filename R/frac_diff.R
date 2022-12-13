#' @title frac_diff
#' @description Takes fractional differences of order d
#' @param x A (vector) time series
#' @param d The (scalar) integration order

#' @return The fractionally differenced x
#' @examples frac_diff(rnorm(100), -0.75)
#'
#' @export
frac_diff <- function(x, d)
{
    N  <-  length(x)
    if (any(is.na(x))) {warning("NA values in frac_diff: treated as zero!"); x[is.na(x)] <- 0}
    y  <- stats::filter(c(rep(0,N),x),
                        filter=c(1,cumprod(((1:N)-1-d)/(1:N))),
                        sides = 1)[-(1:N)]
    if (is.ts(x)) y <- ts(y,start=start(y),frequency=frequency(y))
    return(y)
}

