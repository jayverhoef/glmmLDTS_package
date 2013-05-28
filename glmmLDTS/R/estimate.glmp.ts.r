#--------- CREATE ESTIMATES FROM A QLM.TS OBJECT

estimate.glmp.ts <- function(glmp.ts.object,
  fixed.effect.list = NULL)
{
	work <- glmp.ts.object$fixed.effects[,1:4]
	lin.comb.wts <- matrix(rep(0, times =
    length(work[,1])), ncol = 1)
	if(work[1,"effect"] == "intercept")
    lin.comb.wts[1] <- 1
  for(i in 1:length(fixed.effect.list)) {
    # if numeric
    if(work[work[,"effect"] ==
         fixed.effect.list[[i]]$effect,"levels"][1] == "") {
      lin.comb.wts[work[,"effect"] == fixed.effect.list[[i]]$effect] <-
       as.numeric(fixed.effect.list[[i]]$weight)
    }
    # if categorical (factor)
    else {
      for(j in 1:length(fixed.effect.list[[i]]$level)) {
        lin.comb.wts[ work[,"effect"] == fixed.effect.list[[i]]$effect &
          work[,"levels"] == as.character(fixed.effect.list[[i]]$level[j])] <-
          fixed.effect.list[[i]]$weight[j]
      }
    }
  }
	ind <- !is.na(glmp.ts.object$fixed.effects[,"std.err"])
	lin.comb.wts <- matrix(lin.comb.wts[ind], ncol = 1)
  lin.comb.est <- sum(lin.comb.wts*matrix(work[ind,"estimate"], ncol = 1))
  lin.comb.var <- t(lin.comb.wts) %*% as.matrix(glmp.ts.object$covb) %*%
    lin.comb.wts
  lin.comb.025 <- lin.comb.est - 1.96*sqrt(lin.comb.var)
  lin.comb.975 <- lin.comb.est + 1.96*sqrt(lin.comb.var)
  antilogit.est <- exp(lin.comb.est)/(1 + exp(lin.comb.est))
  antilogit.var <- lin.comb.var*exp(2*lin.comb.est)/(1 + exp(lin.comb.est))^4
  antilogit.025 <- exp(lin.comb.025)/(1 + exp(lin.comb.025))
  antilogit.975 <- exp(lin.comb.975)/(1 + exp(lin.comb.975))
  corr.fact.est <- (1 + exp(lin.comb.est))/exp(lin.comb.est)
  corr.fact.var <- lin.comb.var/exp(2*lin.comb.est)
  corr.fact.025 <- (1 + exp(lin.comb.975))/exp(lin.comb.975)
  corr.fact.975 <- (1 + exp(lin.comb.025))/exp(lin.comb.025)
  data.frame(
    scale = c("logit", "probability", "reciprocal probability"),
    estimate = c(lin.comb.est, antilogit.est, corr.fact.est),
    variance = c(lin.comb.var, antilogit.var, corr.fact.var),
    lower.95CI = c(lin.comb.025, antilogit.025, corr.fact.025),
    upper.95CI = c(lin.comb.975, antilogit.975, corr.fact.975)
    )

}

