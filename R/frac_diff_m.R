#' @title frac_diff_m
#' @description Multivariate fractional differencing function
#' @param x A (n x p) matrix, where each column corresponds to a time series
#' @param d Either a scalar or a vector of length p defining the amount of differencing

#' @return The fractionally differenced x
#' @examples frac_diff_m(matrix(rnorm(200), 100, 2), c(-0.75, -1.25))
#'
#' @export
frac_diff_m <- function(x,d)
{
    y <- as.data.frame(x)
    y <- mapply(frac_diff,x=y,d=d)
    if (is.ts(x)) y <- ts(y,start=start(x),frequency=frequency(x))
    if (is.matrix(x)) y <- as.matrix(y)
    return(y)
}
