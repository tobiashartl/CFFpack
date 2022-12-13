XToZ <- function(X) 0.5*(log(1+X)-log(1-X))

ZToX <- function(Z) (exp(2*Z)-1)/(exp(2*Z)+1)

ARToZ <- function(AR)
{
    X <- arfima::ARToPacf(AR)
    0.5*(log(1+X)-log(1-X))
}

ZToAR <- function(Z)
{
    X <- (exp(2*Z)-1)/(exp(2*Z)+1)
    arfima::PacfToAR(X)
}

MSE <- function(Z, p, d, n)
{
    ar <- ZToAR(Z[1:p])
    ma <- -ZToAR(Z[-(1:p)])
    log(sum((n:1)*(ARMAtoMA(ar,ma,n)-ma_inf(-d,n)[-1])^2)+0.0001)
}

ARToZi1 <- function(AR)
{
    ARs <- -c(polynom::polynomial(c(1,-AR))/polynom::polynomial(c(1,-1)))[-1]
    X <- arfima::ARToPacf(AR)
    0.5*(log(1+X)-log(1-X))
}

ZToARi1 <- function(Z)
{
    if (any(is.na(Z))) return(NA)
    X <- (exp(2*Z)-1)/(exp(2*Z)+1)
    ARs <- arfima::PacfToAR(X)
    AR <- -c(polynom::polynomial(c(1,-ARs))*polynom::polynomial(c(1,-1)))[-1]
    if (length(AR)<length(Z)+1) AR <- c(AR,rep(0,length(Z)+1-length(AR)))
    return(AR)
}

MSEi1 <- function(Z, p, d, n)
{
    ar <- ZToARi1(Z[seq_len(p-1)])
    ma <- -ZToAR(Z[p:(2*p-1)])
    log(sum((n:1)*(ARMAtoMA(ar,ma,n)-ma_inf(-d,n)[-1])^2)+0.0001)
}

sinlink <- function(x,x0,x1)
{
    y <- sin(((x-x0)/(x1-x0)-0.5)*pi)/2+0.5
    y[x<x0] <- 0
    y[x>x1] <- 1
    y
}
