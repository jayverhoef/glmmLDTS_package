# ----------- CREATE 1D PARAMETER FUNCTION FROM LIKELIHOOD
# ----------- MINUS 2 LOG-PSEUDO-LIKELIHOOD

m2LL.1D <- function(parm, pos, parms, yt, X, Z, n.per.term,
  grp.vec, Del.i, A.5, timevec, dX, n, EstMeth = "REML")
{
  parms[pos] <- parm
  m2PLL(parms, yt, X, Z, n.per.term, grp.vec,
    Del.i, A.5, timevec, dX, n, EstMeth = EstMeth)
}

