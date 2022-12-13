#' @title ma_inf
#' @description Returns the Wald coefficients of the fractional differencing polynomial
#' @param d A (scalar) integration order
#' @param N An integer giving the length of the expansion

#' @return The power series expansion (Wald coefficients) of the fractional differencing polynomial
#' @examples ma_inf(0.75, 10)
#'
#' @export

# Coefficients of fractional difference operator
ma_inf <- function(d, N) c(1,cumprod(((1:N)-1-d)/(1:N)))
