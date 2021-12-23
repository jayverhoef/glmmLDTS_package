# ------- A FUNCTION TO FIT A LOGIT MODEL BACK ON ORIGINAL SCALE

logit.poly.fit <- function(x.start, x.end, poly.parm.vector)
{
  rang <- x.end - x.start
  x <- x.start + rang*(0:1000)/1000
  fit <- rep(0, times = 1001)
  for (i in 1:length(poly.parm.vector))
    fit <-  fit + poly.parm.vector[i]*x^(i - 1)
  cbind(x,exp(fit)/(1 + exp(fit)))
}

