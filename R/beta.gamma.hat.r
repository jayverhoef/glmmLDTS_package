# ----------- Estimate Beta and Gamma

beta.gamma.hat <- function(parms, grp.vec, timevec,
  X, Z, yt, A.5, Del.i, n.per.term, lambda = 0,
  betahat.current = NULL)
{
  X <- as.matrix(X)
  alpha <- exp(parms[1])
  S.parms <- exp(parms[2:length(parms)])
  fp <- length(X[1,])
  rp <- length(Z[1,])
  LG.out <- LoopGroup(alpha, yt, X, Z, grp.vec, Del.i, A.5, timevec)
  LX <- LG.out$LX
  LZ <- LG.out$LZ
  Ly <- LG.out$Ly
  SZSigZ <- diag(rep(1/S.parms, times = n.per.term)) + t(LZ) %*% LZ
  solve1 <- solve(SZSigZ, t(LZ) %*% LX)
  solve2 <- solve(SZSigZ, t(LZ) %*% Ly)
  XQX <- t(LX) %*% LX - t(LX) %*% LZ %*% solve1
  XQy <- t(LX) %*% Ly - t(LX) %*% LZ %*% solve2
  ZQX <- t(LZ) %*% LX - t(LZ) %*% LZ %*% solve1
  ZQy <- t(LZ) %*% Ly - t(LZ) %*% LZ %*% solve2
  if(!is.null(betahat.current)) {
    v <- abs(betahat.current)
    XQXlam <- XQX + diag(lambda*v)
  }
  else
    XQXlam <- XQX + diag(as.matrix(rep(lambda, times = length(XQX[,1]))))
  covb <- mginv(XQXlam)
  covb1 <- mginv(XQX)
  beta.hat <- covb %*% XQy
#  beta.hat <- it.solve(XQXlam,XQy)
  gamma.hat <- diag(rep(S.parms, times = n.per.term)) %*%
    (ZQy - ZQX %*% beta.hat)
  qQq <- t(Ly) %*% Ly - t(Ly) %*% LZ %*% solve2 -
    t(XQy) %*% beta.hat
  list(beta.hat = beta.hat, gamma.hat = gamma.hat,
    covb = covb1, qQq = qQq)
}

