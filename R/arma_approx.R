#' @title ARMA_approx
#' @description For a given order (p, q) of an ARMA polynomial, the function finds the best ARMA coefficients to approximate the fractional differencing polynomial.
#' This is done by first by choosing the ARMA coefficients that minimize the MSE between the Wald representation of the ARMA process and the Wald representation of the fractional differencing polynomial for a fixed d.
#' The procedure repeats for a fine grid of integration orders d, such that ARMA approximations are estimated for each grid point.
#' To obtain a smooth function in d, the function uses splines to connect neighboring ARMA coefficients. For details, see Hartl & Jucknewitz (2022, ER):
#' Automatically stores the resulting functions e.g. for (p, q)= (3, 3), n=250 as 'd2arma33_n250.RData'
#' Approximate state space modelling of unobserved fractional components.
#' @param p AR-Order of ARMA polynomial. Three is suggested
#' @param q MA-Order of ARMA polynomial. Three is suggested
#' @param n Length of the underlying Wold representation to be approximated. Should be equal to the length of the time series of interest
#' @param R Number of grid points
#' @param ncl Number of cores to be used for parallel computation

#' @return
#' @examples ARMA_approx(n=250)
#'
#' @export
ARMA_approx <- function(p=3, q=3, n, R=10000, ncl=4){

    ds <- c(-0.5,seq(-0.45,1.2,0.1))

    ### random starting values
    Zstart <- cbind(runif(R,0,12),matrix(runif(R*(p+q-1),-5,5),R))
    Zstart <- array(Zstart,c(R,p+q,length(ds)))
    Zest1 <- array(NA, c(R,p+q,length(ds)))
    MSEest1 <- array(NA,c(R,length(ds)))

    # Initialize cluster and export functions
    cl <- makeCluster(ncl)
    clusterEvalQ(cl, library(parallel))
    clusterExport(cl=cl, varlist=c("ds", "Zstart", "Zest1", "MSEest1", "p", "q", "R", "n"), envir=environment())
    clusterExport(cl, ls("package:CFFpack"))


    Zglobal <- parSapplyLB(cl, 1:length(ds),
                           function(j)
                           {
                               for (i in 1:R)
                               {try({
                                   opt <- optim(fn=MSE, par=Zstart[i,,j],
                                                p=p, d=ds[j], n=n, control=list(maxit=10000))
                                   Zest1[i,,j] <- opt$par
                                   MSEest1[i,j] <- opt$value
                                   cat(ds[j],"mse",opt$value,"\n")
                               })}
                               return(Zest1[which.min(MSEest1[,j]),,j])
                           })

    Zstart <- t(Zglobal)

    Zest2 <- array(NA, c(length(ds),p+q,length(ds)))
    MSEest2 <- array(NA,c(length(ds),length(ds)))
    clusterExport(cl, c("Zstart", "Zest2", "Zglobal", "MSEest2"), envir=environment())

    Zglobal <- parSapplyLB(cl, 1:length(ds),
                           function(j){
                               for (i in 1:length(ds))
                               {try({
                                   opt <- optim(fn=MSE, par=Zstart[i,],
                                                p=p, d=ds[j], n=n, control=list(maxit=10000))
                                   Zest2[i,,j] <- opt$par
                                   MSEest2[i,j] <- opt$value
                                   cat(ds[j],"mse",opt$value,"\n")
                               })}
                               return(Zest2[which.min(MSEest2[,j]),,j])
                           })

    stopCluster(cl)



    ####################################################
    ### Fine grid
    ####################################################

    ds_fine <- seq(-0.5,1.2,0.01)
    Zglobal_fine1 <- matrix(NA,p+q,length(ds_fine))
    Zglobal_fine1[,1] <- Zglobal[,1]

    Zest_fine1 <- array(NA,c(4+length(ds),p+q,length(ds_fine)))
    MSEest_fine1 <- array(NA,c(4+length(ds),length(ds_fine)))

    for (s in 2:length(ds_fine))
    {
        jj0 <- max(which(ds<=ds_fine[s]))
        js_fine <- seq(s-4,s-1,1)
        js_fine[js_fine<1] <- 1
        Zstart_fine <- rbind(t(Zglobal),t(Zglobal_fine1[,js_fine]))

        for (i in 1:nrow(Zstart_fine))
        {try({
            opt <- optim(fn=MSE, par=Zstart_fine[i,],
                         p=p, d=ds_fine[s], n=n, control=list(maxit=10000))
            Zest_fine1[i,,s] <-  opt$par
            MSEest_fine1[i,s] <- opt$value
            cat(ds_fine[s],"mse",opt$value,"\n")
        })}
        Zglobal_fine1[,s] <- Zest_fine1[which.min(MSEest_fine1[,s]),,s]
    }

    which(is.na(Zest_fine1),arr.ind=TRUE)
    Zglobal_fine2 <- Zglobal_fine1

    Zest_fine2 <- array(NA,c(11,p+q,length(ds_fine)))
    MSEest_fine2 <- array(NA,c(11,length(ds_fine)))

    for (s in length(ds_fine):1)
    {
        js_fine <- seq(s-5,s+5,1)
        js_fine[js_fine>s] <- s
        js_fine[js_fine<1] <- 1
        Zstart_fine <- t(Zglobal_fine2[,js_fine])

        for (i in 1:nrow(Zstart_fine))
        {try({
            opt <- optim(fn=MSE, par=Zstart_fine[i,],
                         p=p, d=ds_fine[s], n=n, control=list(maxit=10000))
            Zest_fine2[i,,s] <-  opt$par
            MSEest_fine2[i,s] <- opt$value
            cat(ds_fine[s],"mse",opt$value,"\n")
        })}
        Zglobal_fine2[,s] <- Zest_fine2[which.min(MSEest_fine2[,s]),,s]
    }

    Zglobal_fine <- Zglobal_fine2
    j0 <- which(ds_fine==0)
    Zglobal_fine[,j0] <- rowMeans(Zglobal_fine[,c(j0-1,j0+1)])
    j1 <- which(ds_fine==1)
    Zglobal_fine[,j1] <- rowMeans(Zglobal_fine[,c(j1-1,j1+1)])

    results_i0 <- list("ds","ds_fine", "MSEest_fine1", "MSEest_fine2", "Zest_fine1", "Zest_fine2", "Zglobal", "Zglobal_fine", "Zglobal_fine1", "Zglobal_fine2")

    #save(list=c("ds","ds_fine", "MSEest_fine1", "MSEest_fine2", "Zest_fine1", "Zest_fine2", "Zglobal", "Zglobal_fine", "Zglobal_fine1", "Zglobal_fine2"),
    #     file=paste("grid_d2arma",p,q,"i0_n",n,".RData",sep=""))



    # ==========================================================================
    # Smoothing part for I(0)
    t1 <- 150


    # smooth
    x <- ds_fine[1:t1]
    Y <- t(Zglobal_fine)[1:t1,]

    m <- trunc(length(x)/8)
    k <- 4
    knots <- seq(min(x), max(x), length.out=m+2)
    outer_min <- min(x)-((k-1)):1*(max(x)-min(x))/(m+1)
    outer_max <- max(x)+(1:(k-1))*(max(x)-min(x))/(m+1)
    knots  <- c(outer_min,knots ,outer_max)
    est <- lm(Y ~ -1 + splines::splineDesign(x=x, knots=knots, ord=k))


    # create resulting function
    knots33i0 <- knots
    estpar33i0 <- est$coef
    d2X33i0 <- function(d)
    {
        k <- 4
        Z <- c(splines::splineDesign(x=d, knots=knots33i0, ord=k, deriv=0, outer.ok=TRUE)%*%estpar33i0)
        return(Z)
    }

    d2arma33i0 <- function(d)
    {
        k <- 4
        Z <- d2X33i0(d)
        ar <- ZToAR(Z[1:3])
        ma <- -ZToAR(Z[4:6])
        return(list(ar=ar,ma=ma))
    }


    # save
    #save(
    #    d2arma33i0, d2X33i0, estpar33i0, knots33i0, file=paste("d2arma33i0_n",n,".RData",sep=""))




    # =========================================================================
    # I(1) part
    ds <- c(0.8,seq(0.85,2.2,0.1))

    ### random starting values
    Zstart <- matrix(runif(R*(p+q-1),-6,6),R)
    Zstart <- array(Zstart,c(R,p+q-1,length(ds)))


    Zest1 <- array(NA, c(R,p+q-1,length(ds)))
    MSEest1 <- array(NA,c(R,length(ds)))

    cl <- makeCluster(ncl)
    clusterExport(cl, c("ds", "MSEi1", "Zstart", "p", "q", "n", "Zest1",
                        "MSEest1") ,envir=environment())
    clusterExport(cl, ls("package:CFFpack"))

    Zglobal <- parSapplyLB(cl, 1:length(ds),
                           function(j)
                           {
                               for (i in 1:R)
                               {try({
                                   opt <- optim(fn=MSEi1, par=Zstart[i,,j],
                                                p=p, d=ds[j], n=n, control=list(maxit=10000))
                                   Zest1[i,,j] <- opt$par
                                   MSEest1[i,j] <- opt$value
                                   cat(ds[j],"mse",opt$value,"\n")
                               })}
                               return(Zest1[which.min(MSEest1[,j]),,j])
                           })

    Zstart <- t(Zglobal)

    Zest2 <- array(NA, c(length(ds),p+q-1,length(ds)))
    MSEest2 <- array(NA,c(length(ds),length(ds)))
    clusterExport(cl, c("Zstart", "Zest2", "MSEest2"), envir=environment())

    Zglobal <- parSapplyLB(cl, 1:length(ds),
                           function(j){
                               for (i in 1:length(ds))
                               {try({
                                   opt <- optim(fn=MSEi1, par=Zstart[i,],
                                                p=p, d=ds[j], n=n, control=list(maxit=10000))
                                   Zest2[i,,j] <- opt$par
                                   MSEest2[i,j] <- opt$value
                                   cat(ds[j],"mse",opt$value,"\n")
                               })}
                               return(Zest2[which.min(MSEest2[,j]),,j])
                           })

    stopCluster(cl)


    ####################################################
    ### fine grid
    ####################################################

    ds_fine <- seq(0.8,2.2,0.01)

    Zglobal_fine1 <- matrix(NA,p+q-1,length(ds_fine))
    Zglobal_fine1[,1] <- Zglobal[,1]

    Zest_fine1 <- array(NA,c(4+length(ds),p+q-1,length(ds_fine)))
    MSEest_fine1 <- array(NA,c(4+length(ds),length(ds_fine)))

    for (s in 2:length(ds_fine))
    {
        jj0 <- max(which(ds<=ds_fine[s]))
        js_fine <- seq(s-4,s-1,1)
        js_fine[js_fine<1] <- 1
        Zstart_fine <- rbind(t(Zglobal),t(Zglobal_fine1[,js_fine]))

        for (i in 1:nrow(Zstart_fine))
        {try({
            opt <- optim(fn=MSEi1, par=Zstart_fine[i,],
                         p=p, d=ds_fine[s], n=n, control=list(maxit=10000))
            Zest_fine1[i,,s] <-  opt$par
            MSEest_fine1[i,s] <- opt$value
            cat(ds_fine[s],"mse",opt$value,"\n")
        })}
        Zglobal_fine1[,s] <- Zest_fine1[which.min(MSEest_fine1[,s]),,s]
    }

    which(is.na(Zest_fine1),arr.ind=TRUE)
    Zglobal_fine2 <- Zglobal_fine1

    Zest_fine2 <- array(NA,c(11,p+q-1,length(ds_fine)))
    MSEest_fine2 <- array(NA,c(11,length(ds_fine)))

    for (s in length(ds_fine):1)
    {
        js_fine <- seq(s-5,s+5,1)
        js_fine[js_fine>s] <- s
        js_fine[js_fine<1] <- 1
        Zstart_fine <- t(Zglobal_fine2[,js_fine])

        for (i in 1:nrow(Zstart_fine))
        {try({
            opt <- optim(fn=MSEi1, par=Zstart_fine[i,],
                         p=p, d=ds_fine[s], n=n, control=list(maxit=10000))
            Zest_fine2[i,,s] <-  opt$par
            MSEest_fine2[i,s] <- opt$value
            cat(ds_fine[s],"mse",opt$value,"\n")
        })}
        Zglobal_fine2[,s] <- Zest_fine2[which.min(MSEest_fine2[,s]),,s]
    }

    Zglobal_fine <- Zglobal_fine2
    j1 <- which(ds_fine==1)
    Zglobal_fine[,j1] <- rowMeans(Zglobal_fine[,c(j1-1,j1+1)])
    j2 <- which(ds_fine==2)
    Zglobal_fine[,j2] <- rowMeans(Zglobal_fine[,c(j2-1,j2+1)])


    results_i1 <- list("ds","ds_fine", "MSEest_fine1", "MSEest_fine2", "Zest_fine1", "Zest_fine2", "Zglobal", "Zglobal_fine", "Zglobal_fine1", "Zglobal_fine2")



    # ==========================================================================
    # Smoothing for I(1)

    j1 <- which(ds_fine==1)
    Zglobal_fine[,j1] <- rowMeans(Zglobal_fine[,c(j1-1,j1+1)])
    j2 <- which(ds_fine==2)
    Zglobal_fine[,j2] <- rowMeans(Zglobal_fine[,c(j2-1,j2+1)])

    cbind(ds_fine,t(Zglobal_fine))
    t1 <- 120


    # smooth

    x <- ds_fine[1:t1]
    Y <- t(Zglobal_fine)[1:t1,]

    m <- trunc(length(x)/8)
    k <- 4
    knots <- seq(min(x), max(x), length.out=m+2)
    outer_min <- min(x)-((k-1)):1*(max(x)-min(x))/(m+1)
    outer_max <- max(x)+(1:(k-1))*(max(x)-min(x))/(m+1)
    knots  <- c(outer_min,knots ,outer_max)
    est <- lm(Y ~ -1 + splines::splineDesign(x=x, knots=knots, ord=k))



    knots33i1 <- knots
    estpar33i1 <- est$coef

    d2X33i1 <- function(d)
    {
        k <- 4
        Z <- c(splines::splineDesign(x=d, knots=knots33i1, ord=k, deriv=0, outer.ok=TRUE)%*%estpar33i1)
        return(Z)
    }

    d2arma33i1 <- function(d)
    {
        k <- 4
        Z <- d2X33i1(d)
        ar <- ZToARi1(Z[1:2])
        ma <- -ZToAR(Z[-(1:2)])
        return(list(ar=ar,ma=ma))
    }




    # TOGETHER

    d2arma33 <- function(d)
    {
        p1 <- 0.99
        dd <- 0.05
        f <- d2arma33i0(d)
        fi1 <- d2arma33i1(d)
        s <- 1-sinlink(d,p1-dd,p1+dd)
        si1 <- sinlink(d,p1-dd,p1+dd)
        ar=s*f$ar+si1*fi1$ar
        ma=s*f$ma+si1*fi1$ma
        return(list(ar=ar,ma=ma))
    }


    save(d2arma33i0, d2X33i0, estpar33i0, knots33i0,
        d2arma33i1, d2X33i1, estpar33i1, knots33i1,
        d2arma33,
        file=paste("d2arma", p, q, "_n", n,".RData",sep=""))

}

