vech <- function(x) x[lower.tri(x,diag=TRUE)]

vecl <- function(x) x[lower.tri(x)]

vech2mat <- function(x)
{
    k <- which(length(x)==((1:500)*(2:501)/2))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=TRUE)] <- x
    X <- X+t(X)-diag(diag(X))
    X
}

vecl2mat <- function(x)
{
    k <- which(length(x)==((1:500)*(0:499)/2))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=FALSE)] <- x
    X <- X+t(X)
    diag(X) <- 1
    X
}


cholvech <- function(x) vech(t(chol(x)))

cholvech2mat <- function(x)
{
    k <- which(length(x)==((1:500)*(2:501)/2))
    if (any(is.na(x))) return(matrix(NA,k,k))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=TRUE)] <- x
    X%*%t(X)
}

mlog <- function(A)
{
    if (length(A)==1) return(log(A))
    L <- diag(eigen(A)$values)
    V <- eigen(A)$vectors
    V%*%diag(log(diag(L)))%*%t(V)
}

mexp <- function(A)
{
    if (length(A)==1) return(exp(A))
    L <- diag(eigen(A)$values)
    V <- eigen(A)$vectors
    V%*%diag(exp(diag(L)))%*%t(V)
}

mlogvech <- function(x) vech(mlog(x))

mlogvech2mat <- function(x)
{
    if (length(x)==1) return(matrix(exp(x),1,1))
    k <- which(length(x)==((1:500)*(2:501)/2))
    if (any(is.na(x))) return(matrix(NA,k,k))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=TRUE)] <- x
    X <- X+t(X)-diag(diag(X))
    mexp(X)
}

ldlvec <- function(x)
{
    ch <- chol(x)
    dd <- diag(ch)
    L <- t(ch/dd)
    DD <- dd^2
    c(log(DD),L[lower.tri(L)])
}

ldlvec2mat <- function(x)
{
    if (length(x) == 1)
        return(matrix(exp(x), 1, 1))
    k <- which(length(x) == ((1:500) * (2:501)/2))
    L <- diag(k)
    L[lower.tri(L)] <- x[-(1:k)]
    L%*%diag(exp(x[1:k]))%*%t(L)
}

bdiag <- function (...)
{
    if (nargs() == 1)  x <- as.list(...) else x <- list(...)
    n <- length(x)
    if (n == 0)	return(NULL)
    x <- lapply(x, function(y) if (length(y))
        as.matrix(y)
        else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2],
                                                 (y[3] + 1):y[4]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}

mroot <- function(x)
{
    if (is.null(dim(x))|all(dim(x)==c(1,1))) return(sqrt(x))
    e <- eigen(x)
    e$vectors%*%diag(sqrt(e$values))%*%t(e$vectors)
}

# allow number (matrix dimension in lower.tri)
lower.tri <- function (x, diag = FALSE)
{
    if (is.integer(x)&&length(x)==1)  x <- diag(x)
    if (diag) row(x) >= col(x)
    else row(x) > col(x)
}


# Generalization of "lower.tri" for non-quadratic matrices
Lower.tri <- function(d1,d2=d1,diag=FALSE)
{
    if (is.matrix(d1)) {d2 <- ncol(d1); d1 <- nrow(d1)}
    rbind(lower.tri(diag(d2),diag=diag),matrix(TRUE,d1-d2,d2))
}


# function to calculate companion form of AR(p) or VAR(p) matrix (k x kp)
# includes eigenvalue calculation and checking of stability
toComp <- function(x)
{
    if (is.null(dim(x))) x <- matrix(x,1)
    k <- nrow(x)
    p <- ncol(x)/k
    if (p>1) x <- rbind(x, cbind(diag((p-1)*k),matrix(0,(p-1)*k,k)))
    eigv <- eigen(x,only.values=TRUE)$values
    return(list(CompMat=x,eigv=eigv,stable=all(round(abs(eigv), digits = 2)<1)))
}

# rotation to make a matrix lower triangular
rotate_triang <- function(A)
{
    k <- NCOL(A)
    W <- qr.Q(qr(t(A[1:k,])))
    AW <- A%*%W
    # AW[1:k,1:k][upper.tri(diag(k))] <- 0
    list(W=W,A=AW)
}
