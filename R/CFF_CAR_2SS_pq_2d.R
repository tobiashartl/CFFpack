#' @title CFF_CAR_r_trans_pq
#' @description Function places restrictions on transition matrix to ensure that
#' factors are fractionally integrated noise or AR(l) processes
#' @param theta Contains the parameters. These are the the integration orders d,
#' and the AR coefficients Phi
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.

#' @return Returns the restricted transition matrix of the fractionally integrated
#' factor model
CFF_CAR_r_trans_pq <- function(theta,f1,f2,l,pq) {
    lf1 <- length(f1)
    sf1 <- sum(f1)
    d <- unlist(mapply(rep,theta[1:lf1],f1))
    if (sf1>1) Pi <- matrix(apply(t(sapply(d,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ar)),2,diag),sum(f1))
    if (sf1==1) Pi <- matrix(get(paste("d2arma",pq,pq,sep=""))(d)$ar,sf1)

    if (f2==1) Phi <- matrix(theta[-(1:lf1)],1)
    if (f2>1) Phi <- diag(theta[(lf1+1):(lf1+f2)])
    if (f2>1 & l>1) for (j in 2:l) Phi <- cbind(Phi,diag(theta[(lf1+1+f2*(j-1)):(f2*j+lf1)]))

    Psi <- bdiag(Pi,Phi)
    return(c(Psi))
}



#' @title CFF_CAR_r_obs_pq
#' @description Function places triangular restrictions on loadings to ensure
#' identification of the factors
#' @param gamma Contains the parameters. These are the integration orders,
#' the factor loadings for the long memory factors, and the loadings for the
#' short memory factors
#' @param p Number of variables
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.

#' @return Returns the restricted factor loadings
CFF_CAR_r_obs_pq <- function(gamma,p,f1,f2,pq) {
    lf1 <- length(f1)
    d <- unlist(mapply(rep,gamma[1:lf1],f1))
    gamma <- gamma[-(1:lf1)]
    a0 <- gamma[1:p]
    gamma <- gamma[-(1:p)]

    if (f1[1]>1) {
        a1 <- matrix(0,p,f1[1])
        a1[lower.tri(matrix(0,p,f1[1]),diag=TRUE)] <- gamma[1:(f1[1]*p-f1[1]*(f1[1]-1)/2)]
        gamma <- gamma[-(1:(f1[1]*p-f1[1]*(f1[1]-1)/2))]
    } else {
        a1 <- gamma[1:p]
        gamma <- gamma[-(1:p)]
    }
    if (length(f1)>1) {
        for (i in 2:length(f1)) {
            if (f1[i]>1) {
                ai <- matrix(0,p,f1[i])
                ai[lower.tri(matrix(0,p,f1[i]),diag=TRUE)] <- gamma[1:(f1[i]*p-f1[i]*(f1[i]-1)/2)]
                gamma <- gamma[-(1:(f1[i]*p-f1[i]*(f1[i]-1)/2))]
            } else {
                ai <- gamma[1:p]
                gamma <- gamma[-(1:p)]
            }
            a1 <- cbind(a1,ai)
        }
    }

    A1 <- matrix(rep(c(a1),pq+1)*rep(c(t(sapply(d,function(x) c(1,get(paste("d2arma",pq,pq,sep=""))(x)$ma)))),each=p),p)

    if (f2>1)
    {
        a2 <- matrix(0,p,f2)
        a2[lower.tri(a2,diag=TRUE)] <- gamma[1:(f2*p-f2*(f2-1)/2)]
    } else a2 <- gamma[1:p]

    return(c(a0,A1,a2))
}



#' @title CFF_CAR_2_SS_pq
#' @description Function builds the CFF-CAR model as a state space model
#' @param y Observable data
#' @param d Vector of integration orders
#' @param Phi Matrix of AR coefficients
#' @param A Matrix of factor loadings
#' @param H Observation noise variance matrix
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.

#' @return Returns the state space form of the CFF-CAR model
CFF_CAR_2_SS_pq <- function(y, d,Phi,A,H,f1,f2,l,pq) {
    p <- ncol(y)
    lf1 <- length(f1)
    sf1 <- sum(f1)
    d <- unlist(mapply(rep,d,f1))

    # Construct transition matrix
    if (sf1>1) Tt <- matrix(apply(t(sapply(d,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ar)),2,diag),sf1)
    if (sf1==1) Tt <- matrix(get(paste("d2arma",pq,pq,sep=""))(d)$ar,sf1)
    Tt <- toComp(cbind(Tt,matrix(0,sf1,sf1)))$CompMat
    PHI <- toComp(cbind(Phi,matrix(0,f2,f2)))$CompMat
    Tt <- bdiag(1,Tt,PHI)

    # Observation Matrix
    A1 <- A[,2:(sf1+1)]
    A1 <- matrix(rep(c(A1),pq+1)*rep(c(t(sapply(d,function(x) c(1,get(paste("d2arma",pq,pq,sep=""))(x)$ma)))),each=p),p)
    Zt <- cbind(A[,1],A1,A[,-(1:(sf1+1))],matrix(0,p,f2*l))

    # Transition noise coefficients and variance
    Qt <- diag(sf1+f2)
    Rt <- matrix(0,1+sf1*(pq+1)+f2*(l+1),sf1+f2)
    Rt[matrix(c(2:(sf1+1),(sf1*(pq+1)+2):(sf1*(pq+1)+f2+1),1:(sf1+f2)),ncol=2)] <- 1

    # Observation noise variance
    Ht <- H

    # Starting values
    a1 <- c(1,rep(0,sf1*(pq+1)+f2*(l+1)))
    P1 <- bdiag(diag(1+sf1*(pq+1))*0,
                matrix(solve(diag((f2*(l+1))^2)-PHI%x%PHI)%*%c(bdiag(diag(f2),diag(f2*l)*0)),(l+1)*f2,(l+1)*f2))

    model <- SSModel(y~-1+SSMcustom(Z = Zt, T = Tt, R = Rt,
                                    Q = Qt, a1 = a1, P1 = P1), H = Ht)

    return(model)
}


#' @title CFF_KS_2_FACTORS_pq
#' @description Function smoothes out the factors
#' @param ksmo A matrix holding the smoothed states from the KFAS package
#' @param d Vector of integration orders
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.

#' @return Returns the smoothed out factors
CFF_KS_2_FACTORS_pq <- function(ksmo, d, f1, f2, pq){
    sf1 <- sum(f1)
    lf1 <- length(f1)
    d <- unlist(mapply(rep,d,f1))
    if (sf1>1) Rho <- t(matrix(apply(t(sapply(d,function(x) c(1,get(paste("d2arma",pq,pq,sep=""))(x)$ma))),2,diag),sum(f1)))
    if (sf1==1) Rho <- sapply(d,function(x) c(1,get(paste("d2arma",pq,pq,sep=""))(x)$ma))

    if (is.null(ksmo$alphahat)) {
        X <-  ksmo$a[,c(2:(1+(1+pq)*sf1),(2+(1+pq)*sf1):(1+(1+pq)*sf1+f2))]
    } else {
        X <-  ksmo$alphahat[,c(2:(1+(1+pq)*sf1),(2+(1+pq)*sf1):(1+(1+pq)*sf1+f2))]
    }
    X[,1:((1+pq)*sf1)] <- X[,1:((1+pq)*sf1)]%*%Rho
    return(X[,c(1:sf1,(sf1*(1+pq)+1):(sf1*(1+pq)+f2))])
}


# Alternative State-Space form of CFF-CAR model:
# Fractional factors are directly included as states
CFF_CAR_2_SSX_pq <- function(y, d,Phi,A,H,f1,f2,l,pq){
    p <- ncol(y)
    lf1 <- length(f1)
    sf1 <- sum(f1)
    dlong <- c(unlist(mapply(rep,d,f1)))

    # Construct transition matrix
    if (sf1>1) Tt <- matrix(apply(t(sapply(dlong,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ar)),2,diag),sf1)
    if (sf1==1) Tt <- matrix(get(paste("d2arma",pq,pq,sep=""))(dlong)$ar,sf1)
    Tt <- t(toComp(cbind(Tt,matrix(0,sf1,sf1)))$CompMat)
    PHI <- toComp(Phi)$CompMat
    Tt <- bdiag(1,Tt,PHI)

    # Observation Matrix

    Zt <- cbind(A[,1],A[,2:(sf1+1)],matrix(0,p,sf1*pq),A[,-(1:(sf1+1))])

    # Transition noise coefficients and variance
    Qt <- diag(sf1+f2)
    Rt <- matrix(0,1+sf1*(pq+1)+f2*l,sf1+f2)
    Rt[matrix(c(2:(sf1+1),(sf1*(pq+1)+2):(sf1*(pq+1)+f2+1),1:(sf1+f2)),ncol=2)] <- 1
    Rt[2:(sf1*(pq+1)+1),1:sf1] <- rbind(diag(sf1),t(matrix(apply(t(sapply(dlong,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ma)),2,diag),sf1)))

    # Observation noise variance
    Ht <- H

    # Starting values
    if (l>1) vareta2 <- bdiag(diag(f2),diag(f2*(l-1))*0)
    if (l==1) vareta2 <- diag(f2)
    a1 <- c(1,rep(0,sf1*(pq+1)+f2*l))
    P1 <- bdiag(diag(1+sf1*(pq+1))*0,
                matrix(solve(diag((f2*l)^2)-PHI%x%PHI)%*%c(vareta2),l*f2,l*f2))

    model <- SSModel(y~ -1 + SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1,
                                       P1 = P1), H = Ht)

    return(model)
}


