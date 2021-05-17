# ----------- MINUS 2 PSEUDO-LOGLIKELIHOOD

m2PLL <- function(parms, yt, X, Z, n.per.term, grp.vec,
  Del.i, A.5, timevec, dX, n, EstMeth = "REML")
{
  if(any(abs(parms) > 10)) return(1e35)
  alpha <- exp(parms[1])
  g.parms <- exp(parms[2:length(parms)])
  fp <- length(X[1,])
  rp <- length(Z[1,])
  LG.out <- LoopGroup(alpha, yt, X, Z, grp.vec, Del.i, A.5, timevec)
  LX <- LG.out$LX
  LZ <- LG.out$LZ
  Ly <- LG.out$Ly
  logdet <- LG.out$logdet
  GZSZ <- diag(rep(1/g.parms, times = n.per.term)) + t(LZ) %*% LZ
  solve1 <- solve(GZSZ, t(LZ) %*% LX)
  vec1 <- t(Ly) %*% LX - t(Ly) %*% LZ %*% solve1
  XQX <- t(LX) %*% LX - t(LX) %*% LZ %*% solve1
  qQq.sum <- t(Ly) %*% Ly - t(Ly) %*% LZ %*%
    solve(GZSZ, t(LZ) %*% Ly) - vec1 %*% mginv(XQX) %*% t(vec1)
  XQX.det <- sum(log(Re(svd(XQX)$d)))
  Q.det <- logdet + sum(log(Re(svd(GZSZ)$d))) +
    sum(n.per.term*log(g.parms))
  if(EstMeth == "REML")
    outv <- (n - fp)*log(qQq.sum) + XQX.det + Q.det
  if(EstMeth == "ML")
     outv <- n*log(qQq.sum) + Q.det
  outv

}

