#--------- SIMULATE NORMAL TIME SERIES DATA WITH EXP AUTOCORRELATION

sim.ts.expcorr <- function(n, tsill = 1, trange = 1)
{
  t.inc <- runif(n-1)
  tim <- cumsum(c(0,t.inc))
  tim.mat <- matrix(tim, ncol = 1) %x% matrix(rep(1, times = n), nrow = 1)
  cov.mat <- tsill*exp(-abs(tim.mat - t(tim.mat))/trange)
  data.frame(tim = tim, eps = t(chol(cov.mat)) %*% rnorm(n))
}

