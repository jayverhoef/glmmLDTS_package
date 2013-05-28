# ----------- Loop through groups and store quantities for likelihood evaluation

LoopGroup <- function(alpha, yt, X, Z, grp.vec, Del.i, A.5, timevec)
{
  grp.int <- as.integer(grp.vec)
  LX <- NULL
  LZ <- NULL
  Ly <- NULL
  logdet <- 0
  for (i in 1:max(grp.int)) {
    ind <- grp.int == i
    tvi <- timevec[ind]
    A5i <- A.5[ind]
    Dti <- Del.i[ind]
    Xi <- matrix(X[ind,], nrow = sum(ind))
    Zi <- matrix(Z[ind,], nrow = sum(ind))
    ni <- sum(ind)
    yti <- yt[ind]
    dij <- abs(tvi[1:(ni - 1)] - tvi[2:ni])/alpha
    LM <- Lmats(ni, dij, Xi, Zi, yti, A5i, Dti)
    LX <- rbind(LX,LM$LX)
    LZ <- rbind(LZ,LM$LZ)
    Ly <- rbind(Ly,LM$Ly)
    logdet <- logdet + LM$logdet + sum(log(A5i)) + sum(log(Dti))
  }
    list(LX = LX, LZ = LZ, Ly = Ly, logdet = logdet)
}

