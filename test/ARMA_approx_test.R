devtools::load_all()
n <- 1000
p <- q <- 3
setwd("/Users/tobias/testpakete/CFFpack/data")
library(parallel)
cl <- 6
ARMA_approx(p=3, q=3, n=1000, R=10000,ncl = 6)


arma_appr <- d2arma33(0.7)
ma_d <- frac_diff(c(1, rep(0, 20)), -0.7)
ma <- ARMAtoMA(ar = arma_appr$ar, ma = c(arma_appr$ma), lag.max=20)
ma_d[-1] - ma



