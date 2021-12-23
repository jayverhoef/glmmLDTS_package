# ------- COMPUTE LX, LZ, AND Ly AS BELOW EQ. (16) IN MY PAPER
# ------- ALSO COMPUTE LOG DETERMINANT AS IN EQ. (17) OF MY PAPER

Lmats <- function(n, dij, X, Z, y, A.5, Del.i) {
  nf <- length(X[1,])
  nr <- length(Z[1,])
  LX <- matrix(0, nrow = n, ncol = nf)
  LZ <- matrix(0, nrow = n, ncol = nr)
  Ly <- matrix(0, nrow = n, ncol = 1)
  X <- X/(A.5*Del.i)
  Z <- Z/(A.5*Del.i)
  y <- y/(A.5*Del.i)
  logdet <- 0
  if(n > 1) {
    for (i in 1:(n - 1)) {
      den <- sqrt(1 - exp(-2*dij[i]))
      LX[i,] <- (X[i,] - exp(-dij[i])*X[i + 1,])/den
      LZ[i,] <- (Z[i,] - exp(-dij[i])*Z[i + 1,])/den
      Ly[i] <- (y[i] - exp(-dij[i])*y[i + 1])/den
      logdet <- logdet + log(1 - exp(-2*dij[i]))
    }
  }
  LX[n,] <- X[n,]
  LZ[n,] <- Z[n,]
  Ly[n] <- y[n]

  list(LX = as.matrix(LX), LZ = as.matrix(LZ), Ly = Ly, logdet = logdet)
}

