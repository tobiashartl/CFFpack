#' @title CFF_CAR_LL_pq
#' @description Negative log likelihood for the CFF-CAR model together with the score
#' @param parm Vector of parameters. It contains the memory parameters, the
#' autoregressive parameters for the short memory factors, the factor loadings,
#' and the variances of the noise terms in the observations equation
#' @param Y Observable data
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.
#' @return Negative log likelihood for the CFF-CAR model together with the
#' gradient of the LL
#' @export
CFF_CAR_LL_pq <- function(parm,Y, f1,f2,l,pq)
{
    n <- nrow(Y)
    p <- ncol(Y)
    lf1 <- length(f1)
    sf1 <- sum(f1)
    k <- sf1+f2 # no. factors
    kss <- 1+sf1*(pq+1)+f2*(l+1) # state dimension
    act <- c(2:(sf1+1),(sf1*(pq+1)+2):(sf1*(pq+1)+f2+1))
    lags <- c((1+sf1+1):(1+sf1+pq*sf1),(1+sf1*(pq+1)+f2+1):(1+sf1*(pq+1)+f2+f2*l))
    obs <- c(1:(1+sf1+pq*sf1),(1+sf1*(pq+1)+1):(1+sf1*(pq+1)+f2))
    trans2 <- c(2:(1+sf1*pq),(1+sf1*pq+sf1+1):(1+sf1*pq+sf1+f2*l))

    d <- parm[1:lf1]
    phi <- parm[(lf1+1):(lf1+f2*l)]
    a <- parm[(lf1+f2*l+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2)]
    dlong <- unlist(mapply(rep,d,f1))
    H <- diag(parm[(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+p)])
    #if (any(d>1 | d<0)) return(Inf)
    if (any(diag(H)<0)) return(Inf)

    # Construct transition matrix
    if (sf1>1) Tt <- matrix(apply(t(sapply(dlong,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ar)),2,diag),sf1)
    if (sf1==1) Tt <- matrix(get(paste("d2arma",pq,pq,sep=""))(d)$ar,sf1)
    Tt <- toComp(cbind(Tt,matrix(0,sf1,sf1)))$CompMat
    if (f2==1) Phi <- matrix(phi,1)
    if (f2>1) Phi <- diag(phi[1:f2])
    if (f2>1&l>1) for (j in 2:l) Phi <- cbind(Phi,diag(phi[(j-1)*f2+(1:f2)]))
    tocomp <- toComp(cbind(Phi,matrix(0,f2,f2)))
    if (any(!tocomp$stable)) return(Inf)
    PHI <- tocomp$CompMat
    Tt <- bdiag(1,Tt,PHI)

    # Observation Matrix
    A <- matrix(CFF_CAR_r_obs_pq(c(d,a),p,f1,f2,pq),p)
    Zt <- cbind(A,matrix(0,p,f2*l))

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

    model <- SSModel(Y~-1+SSMcustom(Z = Zt, T = Tt, R = Rt,
                                    Q = Qt, a1 = a1, P1 = P1), H = Ht)
    ksmo <- KFS(model)
    if (any(is.na(ksmo$alphahat))) return(Inf)

    Deriv <- function()
    {
        P <- ksmo$V
        ahat <- ksmo$alphahat

        S_1T_1T <- apply(P[obs,obs,],1:2,sum)+crossprod(ahat[,obs])    #MZ
        S_2T_1Tm1 <- apply(P[act,lags,],1:2,sum) + crossprod(ahat[,act],ahat[,lags])   #MT
        S_1Tm1_1Tm1 <- apply(P[lags,lags,],1:2,sum)+ crossprod(ahat[,lags])

        ehat <- (Y)- ahat[,obs] %*% t(Zt[,obs])
        MH <-  crossprod(ehat)+Zt[,obs]%*%apply(P[obs,obs,],1:2,sum)%*%t(Zt[,obs])-n*H
        Qinv <- diag(k)
        Hinv <- solve(H)

        theta <- c(d,phi)
        restr_trans <- numericDeriv(quote(CFF_CAR_r_trans_pq(theta,f1,f2,l,pq)),"theta")
        S1 <-  attr(restr_trans,"gradient")
        d_gamma <- c(d,a)
        restr_obs <- numericDeriv(quote(CFF_CAR_r_obs_pq(d_gamma,p,f1,f2,pq)),"d_gamma")
        S2 <- attr(restr_obs,"gradient")[,-(1:lf1)]
        S <- bdiag(S1,S2)
        S[(nrow(S1)+1):nrow(S),1:lf1] <- attr(restr_obs,"gradient")[,1:lf1]

        del_l_T <- c(Qinv%*%(S_2T_1Tm1-Tt[act,trans2]%*%S_1Tm1_1Tm1))
        del_l_Z <- Hinv%*%(t(Y)%*%(ahat[,obs])-Zt[,obs]%*%S_1T_1T)

        del_l_d_phi_a <- c(c(del_l_T,del_l_Z)%*%S)
        del_l_h <- diag(Hinv%*%MH%*%Hinv-0.5*diag(diag(Hinv%*%MH%*%Hinv)))

        del_l_PAR <- c(del_l_d_phi_a,del_l_h)
        return(-del_l_PAR)
    }
    LL <- -ksmo$logLik
    attr(LL,"gradient") <- matrix(Deriv(),1)
    return(LL)
}


#' @title CFF_CAR_LLval_pq
#' @description Negative log likelihood for the CFF-CAR model
#' @param parm Vector of parameters. It contains the memory parameters, the
#' autoregressive parameters for the short memory factors, the factor loadings,
#' and the variances of the noise terms in the observations equation
#' @param Y Observable data
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.
#' @return Negative log likelihood for the CFF-CAR model
CFF_CAR_LLval_pq <- function(parm,Y, f1,f2,l,pq)
{
    n <- nrow(Y)
    p <- ncol(Y)
    lf1 <- length(f1)
    sf1 <- sum(f1)
    k <- sf1+f2 # no. factors
    kss <- 1+sf1*(pq+1)+f2*(l+1) # state dimension
    act <- c(2:(sf1+1),(sf1*(pq+1)+2):(sf1*(pq+1)+f2+1))
    lags <- c((1+sf1+1):(1+sf1+pq*sf1),(1+sf1*(pq+1)+f2+1):(1+sf1*(pq+1)+f2+f2*l))
    obs <- c(1:(1+sf1+pq*sf1),(1+sf1*(pq+1)+1):(1+sf1*(pq+1)+f2))
    trans2 <- c(2:(1+sf1*pq),(1+sf1*pq+sf1+1):(1+sf1*pq+sf1+f2*l))

    d <- parm[1:lf1]
    phi <- parm[(lf1+1):(lf1+f2*l)]
    a <- parm[(lf1+f2*l+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2)]
    dlong <- unlist(mapply(rep,d,f1))
    H <- diag(parm[(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+p)])
    #if (any(d>1 | d<0)) return(Inf)
    if (any(diag(H)<0)) return(Inf)

    # Construct transition matrix
    if (sf1>1) Tt <- matrix(apply(t(sapply(dlong,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ar)),2,diag),sf1)
    if (sf1==1) Tt <- matrix(get(paste("d2arma",pq,pq,sep=""))(d)$ar,sf1)
    Tt <- toComp(cbind(Tt,matrix(0,sf1,sf1)))$CompMat
    if (f2==1) Phi <- matrix(phi,1)
    if (f2>1) Phi <- diag(phi[1:f2])
    if (f2>1&l>1) for (j in 2:l) Phi <- cbind(Phi,diag(phi[(j-1)*f2+(1:f2)]))
    tocomp <- toComp(cbind(Phi,matrix(0,f2,f2)))
    if (any(!tocomp$stable)) return(Inf)
    PHI <- tocomp$CompMat
    Tt <- bdiag(1,Tt,PHI)

    # Observation Matrix
    A <- matrix(CFF_CAR_r_obs_pq(c(d,a),p,f1,f2,pq),p)
    Zt <- cbind(A,matrix(0,p,f2*l))

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

    model <- KFAS::SSModel(Y ~ -1 +
                         SSMcustom(Z = Zt, T = Tt, R =Rt, Q = Qt, a1 = a1,
                                          P1 = P1),
                     H = Ht)

    return(-logLik(model))
}





#' @title CFF_CAR_LLderiv_pq
#' @description Derivative of the log likelihood for the CFF-CAR model
#' @param parm Vector of parameters. It contains the memory parameters, the
#' autoregressive parameters for the short memory factors, the factor loadings,
#' and the variances of the noise terms in the observations equation
#' @param Y Observable data
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.
#' @return Derivative of the log likelihood for the CFF-CAR model
#' @export
CFF_CAR_LLderiv_pq <- function(parm,Y, f1,f2,l,pq)
{
    n <- nrow(Y)
    p <- ncol(Y)
    lf1 <- length(f1)
    sf1 <- sum(f1)
    k <- sf1+f2 # no. factors
    kss <- 1+sf1*(pq+1)+f2*(l+1) # state dimension
    act <- c(2:(sf1+1),(sf1*(pq+1)+2):(sf1*(pq+1)+f2+1))
    lags <- c((1+sf1+1):(1+sf1+pq*sf1),(1+sf1*(pq+1)+f2+1):(1+sf1*(pq+1)+f2+f2*l))
    obs <- c(1:(1+sf1+pq*sf1),(1+sf1*(pq+1)+1):(1+sf1*(pq+1)+f2))
    trans2 <- c(2:(1+sf1*pq),(1+sf1*pq+sf1+1):(1+sf1*pq+sf1+f2*l))

    d <- parm[1:lf1]
    phi <- parm[(lf1+1):(lf1+f2*l)]
    a <- parm[(lf1+f2*l+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2)]
    dlong <- unlist(mapply(rep,d,f1))
    H <- diag(parm[(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+p)])

    # Construct transition matrix
    if (sf1>1) Tt <- matrix(apply(t(sapply(dlong,function(x) get(paste("d2arma",pq,pq,sep=""))(x)$ar)),2,diag),sf1)
    if (sf1==1) Tt <- matrix(get(paste("d2arma",pq,pq,sep=""))(d)$ar,sf1)
    Tt <- toComp(cbind(Tt,matrix(0,sf1,sf1)))$CompMat
    if (f2==1) Phi <- matrix(phi,1)
    if (f2>1) Phi <- diag(phi[1:f2])
    if (f2>1&l>1) for (j in 2:l) Phi <- cbind(Phi,diag(phi[(j-1)*f2+(1:f2)]))
    tocomp <- toComp(cbind(Phi,matrix(0,f2,f2)))
    if (any(!tocomp$stable)) return(Inf)
    PHI <- tocomp$CompMat
    Tt <- bdiag(1,Tt,PHI)

    # Observation Matrix
    A <- matrix(CFF_CAR_r_obs_pq(c(d,a),p,f1,f2,pq),p)
    Zt <- cbind(A,matrix(0,p,f2*l))

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

    model <- KFAS::SSModel(Y ~ -1 +
                               SSMcustom(Z = Zt, T = Tt, R =Rt, Q = Qt, a1 = a1,
                                         P1 = P1),
                           H = Ht)
    ksmo <- KFS(model)
    if (any(is.na(ksmo$alphahat))) return(Inf)
    #if (any(d>1 | d<0)) return(Inf)

    P <- ksmo$V
    ahat <- ksmo$alphahat

    S_1T_1T <- apply(P[obs,obs,],1:2,sum)+crossprod(ahat[,obs])    #MZ
    S_2T_1Tm1 <- apply(P[act,lags,],1:2,sum) + crossprod(ahat[,act],ahat[,lags])   #MT
    S_1Tm1_1Tm1 <- apply(P[lags,lags,],1:2,sum)+crossprod(ahat[,lags])

    ehat <- (Y)- ahat[,obs] %*% t(Zt[,obs])
    MH <-  crossprod(ehat)+Zt[,obs]%*%apply(P[obs,obs,],1:2,sum)%*%t(Zt[,obs])-n*H
    Qinv <- diag(k)
    Hinv <- solve(H)

    theta <- c(d,phi)
    restr_trans <- numericDeriv(quote(CFF_CAR_r_trans_pq(theta,f1,f2,l,pq)),"theta")
    S1 <-  attr(restr_trans,"gradient")
    d_gamma <- c(d,a)
    restr_obs <- numericDeriv(quote(CFF_CAR_r_obs_pq(d_gamma,p,f1,f2,pq)),"d_gamma")
    S2 <- attr(restr_obs,"gradient")[,-(1:lf1)]
    S <- bdiag(S1,S2)
    S[(nrow(S1)+1):nrow(S),1:lf1] <- attr(restr_obs,"gradient")[,1:lf1]

    del_l_T <- c(Qinv%*%(S_2T_1Tm1-Tt[act,trans2]%*%S_1Tm1_1Tm1))
    del_l_Z <- Hinv%*%(t(Y)%*%(ahat[,obs])-Zt[,obs]%*%S_1T_1T)

    del_l_d_phi_a <- c(c(del_l_T,del_l_Z)%*%S)
    del_l_h <- diag(Hinv%*%MH%*%Hinv-0.5*diag(diag(Hinv%*%MH%*%Hinv)))

    del_l_PAR <- c(del_l_d_phi_a,del_l_h)
    return(-del_l_PAR)
}

### Funktionen, die die Innovationen, die Innovationsvarianzen und die Likelihood als vektoren ausgeben -> um sie dann abzuleiten...

vec_resid <- function(parm,Y, f1,f2,l,pq)
{
    p <- ncol(Y)
    lf1 <- length(f1)

    d <- parm[1:lf1]
    phi <- parm[(lf1+1):(lf1+f2*l)]
    a <- parm[(lf1+f2*l+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2)]
    H <- diag(parm[(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+p)])
    if (f2==1) Phi <- matrix(phi,1)
    if (f2>1) Phi <- diag(phi[1:f2])
    if (f2>1&l>1) for (j in 2:l) Phi <- cbind(Phi,diag(phi[(j-1)*f2+(1:f2)]))
    A <- matrix(CFF_CAR_r_obs_pq(c(d,a),p,f1,f2,pq),p)[,c(1,2:(sum(f1)+1),(sum(f1)*(pq+1)+2):(sum(f1)*(pq+1)+1+f2))]


    model <- CFF_CAR_2_SS_pq(Y, d,Phi,A,H,f1,f2,l,pq)
    kfil <- KFS(model, smoothing="none",simplify=FALSE)
    return(c(kfil$v))
}

vec_sigma <- function(parm,Y, f1,f2,l,pq)
{
    p <- ncol(Y)
    lf1 <- length(f1)

    d <- parm[1:lf1]
    phi <- parm[(lf1+1):(lf1+f2*l)]
    a <- parm[(lf1+f2*l+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2)]
    H <- diag(parm[(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+p)])
    if (f2==1) Phi <- matrix(phi,1)
    if (f2>1) Phi <- diag(phi[1:f2])
    if (f2>1&l>1) for (j in 2:l) Phi <- cbind(Phi,diag(phi[(j-1)*f2+(1:f2)]))
    A <- matrix(CFF_CAR_r_obs_pq(c(d,a),p,f1,f2,pq),p)[,c(1,2:(sum(f1)+1),(sum(f1)*(pq+1)+2):(sum(f1)*(pq+1)+1+f2))]

    model <- CFF_CAR_2_SS_pq(Y, d,Phi,A,H,f1,f2,l,pq)
    kfil <- KFS(model, smoothing="none",simplify=FALSE)

    F <- sapply(1:n,function(t) model$Z[,,1]%*%kfil$P[,,t]%*%t(model$Z[,,1])+model$H[,,1])
    return(c(F))
}

vec_ll <- function(parm,Y, f1,f2,l,pq)
{
    p <- NCOL(Y)
    n <- NROW(Y)
    lf1 <- length(f1)

    d <- parm[1:lf1]
    phi <- parm[(lf1+1):(lf1+f2*l)]
    a <- parm[(lf1+f2*l+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2)]
    H <- diag(parm[(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+1):(lf1+f2*l+p+sum(f1*(p-f1)+f1*(f1+1)/2)+f2*(p-f2)+f2*(f2+1)/2+p)])
    if (f2==1) Phi <- matrix(phi,1)
    if (f2>1) Phi <- diag(phi[1:f2])
    if (f2>1&l>1) for (j in 2:l) Phi <- cbind(Phi,diag(phi[(j-1)*f2+(1:f2)]))
    A <- matrix(CFF_CAR_r_obs_pq(c(d,a),p,f1,f2,pq),p)[,c(1,2:(sum(f1)+1),(sum(f1)*(pq+1)+2):(sum(f1)*(pq+1)+1+f2))]

    model <- CFF_CAR_2_SS_pq(Y, d,Phi,A,H,f1,f2,l,pq)
    kfil <- KFS(model, smoothing="none",simplify=FALSE)

    F <- sapply(1:n,function(t) model$Z[,,1]%*%kfil$P[,,t]%*%t(model$Z[,,1])+model$H[,,1])
    F <- array(F,c(p,p,n))
    ldetf <- apply(F,3,function(x) log(det(x)))
    vFv <- sapply(1:n,function(t) c(kfil$v[,t]%*%solve(F[,,t])%*%kfil$v[,t]))
    l <- -0.5*(ldetf+vFv)
    return(l)
}
