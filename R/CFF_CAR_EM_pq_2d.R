#' @title EMstep_CFF_CAR_pq
#' @description Single Step of the EM algorithm
#' @param Y Observable data
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.
#' @param d Vector of integration orders
#' @param Phi Matrix of AR coefficients
#' @param A Matrix of factor loadings
#' @param H Observation noise variance matrix
#' @return EM-updated parameters and likelihood
EMstep_CFF_CAR_pq <- function(Y,f1,f2,l,pq, d,Phi,A,H)
{
    ## Dimensions
    n <- nrow(Y)
    p <- ncol(Y)
    lf1 <- length(f1)
    sf1 <- sum(f1)
    k <- sf1+f2 # no. factors
    kss <- 1+sf1*(pq+1)+f2*(l+1) # state dimension
    act <- c(2:(sf1+1),(sf1*(pq+1)+2):(sf1*(pq+1)+f2+1))

    ## State Space Form, Kalman Filter and Smoother
    model <- CFF_CAR_2_SS_pq(Y, d,Phi,A,H, f1,f2,l,pq)
    ksmo <- KFS(model)

    ## Moment matrices
    Xa <- rbind(ksmo$alphahat,ksmo$a[n+1,])[,act]   # Smoothed actual state E_n(X[1]),...,E_n(X[n+1])
    Xl <- rbind(ksmo$alphahat,ksmo$a[n+1,])[,-c(1,act)] # Smoothed lagged states
    Xr <- ksmo$alphahat[,1:(sf1*(pq+1)+f2+1)]   # Smoothed state E_n(X[1]),...,E_n(X[n]) used for regression y~X

    PP <- apply(ksmo$V,1:2,sum)
    Pal <- (PP+ksmo$P[,,n+1])[act,-c(1,act)]
    Pll <- (PP+ksmo$P[,,n+1])[-c(1,act),-c(1,act)]
    Prr <- PP[1:(sf1*(pq+1)+f2+1),1:(sf1*(pq+1)+f2+1)]   # Sum of smoothed state variances  Var_n(X[1])+...+Var_n(X[n+1])

    Qrr <- t(Xr)%*%Xr + Prr
    QYr <- t(Y)%*%Xr
    Qll <- t(Xl)%*%Xl+Pll
    Qal <- t(Xa)%*%Xl+ Pal

    # Linearization of restrictions: S (d',phi',a')' +s
    if (l==1) phi <- diag(Phi) else phi <- apply(array(Phi,c(f2,f2,l)),3,diag)
    theta <- c(d,phi)
    gamma <- A[,1]
    A <- as.matrix(A[,-1])
    for (i in 1:lf1)
    {
        gamma <- c(gamma,A[,1:f1[i]][Lower.tri(p,f1[i],diag=TRUE)])
        A <- as.matrix(A[,-(1:f1[i])])
    }
    gamma <- c(gamma,A[,1:f2][Lower.tri(p,f2,diag=TRUE)])
    d_gamma <- c(d,gamma)

    restr_trans <- numericDeriv(quote(CFF_CAR_r_trans_pq(theta,f1,f2,l,pq)),"theta")
    S1 <-  attr(restr_trans,"gradient")

    restr_obs <- numericDeriv(quote(CFF_CAR_r_obs_pq(d_gamma,p,f1,f2,pq)),"d_gamma")
    S2 <- attr(restr_obs,"gradient")[,-(1:lf1)]

    S <- bdiag(S1,S2)
    S[(nrow(S1)+1):nrow(S),1:lf1] <- attr(restr_obs,"gradient")[,1:lf1]
    s <- c(restr_trans,restr_obs)-S%*%c(theta,gamma)

    # Transition matrix estimate
    Q <- diag(k)

    thetagamma <- solve(t(S)%*%bdiag(Qll%x%solve(Q),Qrr%x%solve(H))%*%S,
                        c(t(S)%*%c(solve(Q)%*%(Qal-matrix(s[1:(k*(kss-k-1))],k,kss-k-1)%*%Qll),solve(H)%*%(QYr-matrix(s[(k*(kss-k-1)+1):(length(s))],nrow(QYr),nrow(Qrr))%*%Qrr))))

    d <- thetagamma[1:lf1]
    if (f2==1) Phi <- matrix(thetagamma[(lf1+1):(lf1+l)],1)
    if (f2>1) Phi <- diag(thetagamma[(lf1+1):(lf1+f2)])
    if (f2>1 & l>1) for (j in 2:l) Phi <- cbind(Phi,diag(thetagamma[(lf1+1+f2*(j-1)):(f2*j+lf1)]))

    AA <- matrix(CFF_CAR_r_obs_pq(c(d,thetagamma[-(1:(lf1+f2*l))]),p,f1,f2,pq),p)

    # Observation noise covariance matrix estimate
    E <- Y-Xr%*%t(AA)
    H <- diag(diag(t(E)%*%E+AA%*%Prr%*%t(AA)))/n

    A <- AA[,c(1:(sf1+1),(2+(pq+1)*sf1):(1+(pq+1)*sf1+f2))]

    cat(thetagamma,"\n")

    return(list(d=d,Phi=Phi,A=A,H=H,LL=ksmo$logLik))
}


###
#' @title EM_CFF_CAR_pq
#' @description EM algorithm for parameter estimation of the CFF-CAR model
#' @param Y Observable data
#' @param f1 A vector of length equal the number of fractionally integrated factor
#' groups. Each entry corresponds to the number of factors for each group, in
#' descending memory.
#' @param f2 A scalar that equals the number of short-memory factors
#' @param l Scalar order of the AR polynomials for the stationary factors
#' @param pq Scalar defining the order of the ARMA approximation for the fractionally
#' integrated factors. Typically equal to three.
#' @param d Vector of integration orders
#' @param Phi Matrix of AR coefficients
#' @param A Matrix of factor loadings
#' @param H Observation noise variance matrix
#' @param itermax Maximum number of iterations of the EM. Default is 1000
#' @param stoptol Tolerance for stopping the EM algorithm (in terms of the log likelihood).
#' Default is .1
#' @return EM-updated parameters and likelihood
#'
#'
#' @export
EM_CFF_CAR_pq <- function(Y,f1,f2,l,pq,d,Phi,A,H, itermax=1000, stoptol=0.1)
{
    ### Vektor of Likelihoods
    ll <- rep(NA,itermax)

    ### Iterate until convergence
    for (i in 1:itermax)
    {
        em <- EMstep_CFF_CAR_pq(Y=Y,f1=f1,f2=f2,l=l,pq=pq,d=d,Phi=Phi,A=A,H=H)
        A <- em$A
        Phi <- em$Phi
        H <- em$H
        d <- em$d

        cat("Iteration",i,"Log Likelihood", em$LL,"\n")
        cat("d:",d," Phi: ",c(apply(array(Phi,c(f2,f2,l)),3,diag))," A: ",c(A),"H:",diag(H),"\n")
        ll[i] <- em$LL

        if (i>1) if(ll[i]-ll[i-1]<stoptol) break
    }

    ### Return estimates
    return(list(
        ll = na.omit(ll),
        f1 = f1, f2 = f2, pq = pq, l = l, n = nrow(Y), p = ncol(Y),
        d = d, Phi = Phi, A = A, H = H))
}
