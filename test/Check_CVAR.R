devtools::load_all()
library(KFAS)

### Simulated model for ML
# config
f1 <- c(2, 2)
f2 <- c(2)
p <- 10
n <- 100
d <- c(1.3, 0.8)
ar <- diag(c(0.8, 0.8))


set.seed(42)
# Long memory factors
F1 <- frac_diff_m(matrix(rnorm(n*f1[1]), nrow = n), -d[1])
F2 <- frac_diff_m(matrix(rnorm(n*f1[2]), nrow = n), -d[2])
# Short memory factors
F3 <- sapply(1:f2, function(x) filter(c(rnorm(n)), ar[x,x], "recursive"))
# Measurment errors
U <- matrix(rnorm(n*p), n, p)
# Loadings
A <- matrix(runif((sum(f1, f2)*p), 1, 2), p)
A[1,2] <- A[1, 4] <- A[1, 6] <- 0

# Account for mean
A <- cbind(0, A)

Y <-cbind(F1, F2, F3)  %*% t(A[, -1]) + A[,1]  + U

### common dynamics
plot(Y[,1], type="l")
sapply(2:p, function(x) lines(Y[,x]))

# Task 1: Estimate
l=1
pq=3
em_pq <- EM_CFF_CAR_pq(Y,f1,f2,l=1,pq=3,d=d,Phi=ar,A=A,H=diag(1, p), itermax=100, stoptol=0.1)
parm_em <- c(em_pq$d,diag(em_pq$Phi),em_pq$A[em_pq$A!=0],diag(em_pq$H))

opt <- optim(fn=CFF_CAR_LLval_pq,gr=CFF_CAR_LLderiv_pq, par=parm_em,Y=Y, f1=f1,f2=f2,l=l,pq=pq, method="BFGS", control=list(trace=6,maxit=1000))
opt

parm <- opt$par

cbind(parm, c(d, diag(ar), A[em_pq$A != 0], rep(1, p)))

parm_true <- c(d, diag(ar), A[em_pq$A != 0], rep(1, p))
CFF_CAR_LLderiv_pq(parm_true, Y, f1, f2, l, pq)
cbind(c(attr(numericDeriv(quote(CFF_CAR_LLval_pq(parm_true, Y, f1, f2, l, pq)), "parm_true"), "gradient")),
      CFF_CAR_LLderiv_pq(parm_true, Y, f1, f2, l, pq)
)
