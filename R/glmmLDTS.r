glmmLDTS <- function(fixed.formula,
  random.formula,
  data,
  timecol,
  EstMeth = "REML",
  distribution = "binomial",
  trialscol = NULL,
  group.vec = NULL,
  ridge.reg = NULL,
  lambda = 0)
{
  if(!is.null(ridge.reg)) {
    if(!(ridge.reg == "local" | ridge.reg == "global"))
    return("illegal specification of ridge.reg argument")
  }
  if(lambda < 0)
    return("lambda must be greater than or equal to 0")
  #store the start time
  start.time = Sys.time()
  # create a list of variables in the formula
	varlist <- as.character(attr(terms(fixed.formula, data = data),"variables"))
	# get the name of the response variable
	response.col <- varlist[2]

	# set the WARNINGS to NULL
	WARNINGS <- NULL

	# check that data are binomial if distribution = "binomial"
	if(distribution == "binomial") {
		response_data <- data[,response.col]
		if(all(response_data%%1 > 0)) {
			return("Distribution specified is binomial but non-integer data provided for response variable")
		}
	}

	# check for missing values in response.col, remove with warning
	if(any(is.na(data[,response.col]))){
    data <- data[!is.na(data[,response.col]),]
    WARNINGS <- rbind(WARNINGS,
      "Records with missing data for Response variable were removed")
  }

  # fixed effects design matrix
	X <- model.matrix(fixed.formula, data = data)
	# number of columns in the fixed effects design matrix
	nXc <- length(X[1,])
	# check and remove any columns that are all zeros
	ind.not0s <- rep(FALSE, times = nXc)
	for(i in 1:nXc) ind.not0s[i] <- (sum(abs(X[,i])) != 0)
	X <- X[,ind.not0s]
	if(sum(ind.not0s) < nXc)
    WARNINGS <- rbind(WARNINGS,
      "Some Design matrix columns were all 0")
  	# check for estimability
  	Xre <- rref(X)
	# create X with estimable functions
	keepXcol <- apply(abs(Xre),2,sum) == 1
	if(any(keepXcol == FALSE))
    WARNINGS <- rbind(WARNINGS,
      "Some Design matrix columns eliminated because they were nonestimable")
	if(any(keepXcol == FALSE)) X <- X[,keepXcol]
	# get the new number of columns in the fixed effects design matrix
	nXc <- length(X[1,])


  # random effects design matrix
	terms.list <- attr(terms(random.formula),"term.labels")
	Z0 <- NULL
	# vector to keep track of columns when terms change in Z0
	terms.ind <- NULL
	for (i in 1:length(terms.list)) {
		form1 <- formula(paste("~ ", terms.list[i], " - 1"))
		Z0.temp <- model.matrix(form1, data = data)
		terms.ind <- c(terms.ind, rep(i, times = length(Z0.temp[1,])))
		Z0 <- cbind(Z0,Z0.temp)
	}
	nZc <- length(Z0[1,])
	# check and remove any columns that are all zeros
	ind.not0s <- rep(FALSE, times = nZc)
	for(i in 1:nZc) ind.not0s[i] <- (sum(abs(Z0[,i])) != 0)
	Z0 <- Z0[,ind.not0s]
	if(sum(ind.not0s) < nZc)
    WARNINGS <- rbind(WARNINGS,
      "Some Random Effects Design matrix columns were all 0")
	terms.ind <- terms.ind[ind.not0s]
	# get the new number of columns in the random effects design matrix
	nZc <- sum(ind.not0s)
	# number of random terms
	nrt <- max(terms.ind)
	n.per.term <- NULL
	for(i in 1:nrt) n.per.term <- c(n.per.term, sum(terms.ind == i))

	# vector of response variable
	y <- data[,response.col]

	# sample size
	n <- length(y)

	# if response variable is binomial with n trials change y to proportion
	if(distribution == "binomial") {
	  if(!is.null(trialscol)){
		  trialsvec <- data[,trialscol]
		  y <- y/trialsvec
	  }
	  # if response is bernoulli, set trialsvec to all ones
	  else trialsvec <- rep(1, times = n)
	}

	#create a single vector of time column
	timevec <- as.numeric(data[,timecol])

	#create a single vector of groupings if necessary
  if (!is.null(group.vec)) grp.vec <- data[,group.vec]
  else grp.vec <- rep(1, times = n)
  grp.int <- as.integer(grp.vec)

	# rank of X
	dX <- sum(svd(X)$d>1e-10)

	# Initial parameter estimates, fixed effects

	beta.hat <- rep(0, times = nXc)
	if(distribution == "binomial") {
	  beta.hat[1] <- log((sum(y)/sum(trialsvec))/
		  (1 - sum(y)/sum(trialsvec)))
  }
  if(distribution == "Poisson") {
    beta.hat[1] <- log(mean(y))
  }
	# Initial parameter estimates, random effects
	gamma.hat <- rep(0, times = length(Z0[1,]))
  # Initial covariance parameters after profiling
  psi.hat <- c(0,rep(0,times = nrt))
  psi.new <- c(psi.hat)


		# ----------- START LOOPING HERE ---------------------------

    # create an indicator to stop looping
		stoploop <- 0
		# keep track of the number of iterations
		iter <- 0
		# keep track of number of inner iterations for beta
		inner.iter2 <- NULL
		# begin looping
		if(distribution == "binomial") y <- y*0.99999 + 0.000001
		while(stoploop == 0) {
      psi.current <- psi.new
      beta.current <- beta.hat
      gamma.current <- gamma.hat
			eta.hat <- X %*% beta.hat + Z0 %*% gamma.hat
      #diagonal elements of Delta~^{-1} of my manuscript
      if(distribution == "binomial")
			  Del.i <- as.vector((1 + exp(eta.hat))^2/exp(eta.hat))
      if(distribution == "Gaussian" | distribution == "normal")
			  Del.i <- as.vector(rep(1, times = length(eta.hat)))
      if(distribution == "Poisson")
			  Del.i <- as.vector(1/exp(eta.hat))
      #diagonal elements of A^(1/2) of my manuscript
      if(distribution == "binomial")
			  A.5 <- as.vector(sqrt(exp(eta.hat)/(1 +
          exp(eta.hat))^2/trialsvec))
      if(distribution == "Gaussian" | distribution == "normal")
			  A.5 <- as.vector(rep(1, times = length(eta.hat)))
      if(distribution == "Poisson")
			  A.5 <- as.vector(sqrt(exp(eta.hat)))
			#pseudo data y~ equation (3) in my manuscript
      if(distribution == "binomial")
			  yt <- Del.i*(y - exp(eta.hat)/(1 + exp(eta.hat))) + eta.hat
      if(distribution == "Gaussian" | distribution == "normal")
			  yt <- y
      if(distribution == "Poisson")
			  yt <- Del.i*(y - exp(eta.hat)) + eta.hat

      #loop through covariance parameters once each, and do a one-dimensional
      # optimization for each parameter, holding others fixed
      for(i in 1:length(psi.hat)) {
		psi.est <- optimize(m2LL.1D, pos = i, interval = c(-19, 19),
			parms = psi.hat, yt = yt, X = X, Z = Z0, n.per.term = n.per.term,
			grp.vec = grp.vec, Del.i = Del.i, A.5 = A.5, tol = .01,
			timevec = timevec, dX = dX, n = n, EstMeth = EstMeth)
      psi.hat[i] <- psi.est$minimum
      }

      stplp1 <- 0
      k.it <- 0
      while(stplp1 == 0) {
        # Estimate beta and gamma
        beta.curr <- beta.hat
        gamma.curr <- gamma.hat
        beta.hat <- 1.0*beta.hat
        gamma.hat <- 1.0*gamma.hat
			  eta.hat <- X %*% beta.hat + Z0 %*% gamma.hat
        #diagonal elements of Delta~^{-1} of my manuscript
        if(distribution == "binomial")
			    Del.i <- as.vector((1 + exp(eta.hat))^2/exp(eta.hat))
        if(distribution == "Gaussian" | distribution == "normal")
			    Del.i <- as.vector(rep(1, times = length(eta.hat)))
        if(distribution == "Poisson")
			    Del.i <- as.vector(1/exp(eta.hat))
        #diagonal elements of A^(1/2) of my manuscript
        if(distribution == "binomial")
			    A.5 <- as.vector(sqrt(exp(eta.hat)/(1 +
            exp(eta.hat))^2/trialsvec))
        if(distribution == "Gaussian" | distribution == "normal")
			    A.5 <- as.vector(rep(1, times = length(eta.hat)))
        if(distribution == "Poisson")
			    A.5 <- as.vector(sqrt(exp(eta.hat)))
			  #pseudo data y~ equation (3) in my manuscript
        if(distribution == "binomial")
			    yt <- Del.i*(y - exp(eta.hat)/(1 + exp(eta.hat))) + eta.hat
        if(distribution == "Gaussian" | distribution == "normal")
			    yt <- y
        if(distribution == "Poisson")
			    yt <- Del.i*(y - exp(eta.hat)) + eta.hat
			  if(is.null(ridge.reg)) beta.vec <- NULL
			  else {
			    if(ridge.reg == "global") beta.vec <- 0*beta.hat + 1
			    if(ridge.reg == "local") beta.vec <- beta.hat
        }
        bgh <- beta.gamma.hat(psi.hat, grp.int, timevec,
          X, Z0, yt, A.5, Del.i, n.per.term, lambda, beta.vec)
        beta.hat <- bgh$beta.hat
        gamma.hat <- bgh$gamma.hat
			  if(all(abs(beta.hat - beta.curr)/beta.curr < 1e-4))
          stplp1 <- 1
        if(k.it >= 100) stplp1 <- 1
        k.it <- k.it + 1
      }
      inner.iter2 <- cbind(inner.iter2, k.it)

      #convergence criteria on the covariance parameters
			psi.new <- c(exp(psi.hat))
#			if( any(max(abs(psi.new - psi.current)/psi.current) < 1e-5 |
#				(iter >= 20)) ) stoploop <- 1
      #convergence criteria on the fixed effect parameters
			if(all(abs(psi.new - psi.current)/psi.current < 1e-5))
          stoploop <- 1
			if (iter >= 100) stoploop <- 1

			iter <- iter + 1
		}

		# ----------- DONE LOOPING HERE ---------------------------

		# autocorrelation parameter, need to exponentiate as per
		# estimation function
    alpha <- exp(psi.hat[1])
		# random effects variances, divided by sigma2,
		# need to exponentiate as per estimation function
    S.parms <- exp(psi.hat[2:length(psi.hat)])
    # estimate sigma^2
    if(EstMeth == "REML")
		  sigma2 <- bgh$qQq/(n - nXc)
    if(EstMeth == "ML")
		  sigma2 <- bgh$qQq/n
		# get covariance matrix
    covb <- bgh$covb
    # scale covariance matrix for fact that sigma^2 was factored out
    covb <- as.numeric(sigma2)*covb
	colnames(covb) <- colnames(X)
	rownames(covb) <- colnames(X)
    # adjust random effects variance for fact that sigma^2 was factored
    nu <- sigma2*S.parms

    # fit for each observation based on fixed effects
	mu.t <- as.matrix(X) %*% matrix(beta.hat, nrow = nXc)
	temp1 <- covb %*% t(as.matrix(X))
	# standard errors of fits
	mu.se <- sqrt(rowSums(as.matrix(X) * t(temp1)))
	# go back to original scale and transform conf. intervals
	mu <- exp(mu.t)/(1 + exp(mu.t))
	mu.lci <- exp(mu.t - 1.96*mu.se)/(1 + exp(mu.t - 1.96*mu.se))
	mu.uci <- exp(mu.t + 1.96*mu.se)/(1 + exp(mu.t + 1.96*mu.se))

    # unstandardize beta hat and return new beta.hat and covb
	bhat.se <- sqrt(diag(covb))

    # create fixed effects table
	fixed.eff.est <- data.frame(parameters = beta.hat,
		std.err = bhat.se,
		df = n - dX,
		t.value = beta.hat/bhat.se,
		prob.t = round(100000*(1- pt(abs(beta.hat/bhat.se),
			df = n - dX))*2)/100000
		)
	rownames(fixed.eff.est) <- colnames(X)

	  fixed.effects.estimates <- fixed.effect.table(formula = fixed.formula,
      data = data,
			fixed.effects.estimates = fixed.eff.est)
    if(!is.null(ridge.reg) & fixed.effects.estimates[1,1] == "intercept")
      fixed.effects.estimates[1,3] <- fixed.effects.estimates[1,3] +
         mean(gamma.hat)
    # create a table of fits
	  fit.table <- data.frame(time = data[,timecol], Y = y,
		  mu = mu, mu.lci = mu.lci, mu.uci = mu.uci)

	  typeIII.hypoth <- hypoth(fixed.formula, data, log(c(alpha, S.parms)),
      grp.int, timevec, Z0, yt, A.5, Del.i, n.per.term)

    xout <- list(
      dataset = data,
      WARNINGS = WARNINGS,
      fixed.formula = fixed.formula,
      random.formula = random.formula,
      sample.size = n,
      timecol = timecol,
      trialscol = trialscol,
      group.vec = group.vec,
      ridge.reg = ridge.reg,
      lambda = lambda,
      start.time = start.time,
      end.time = Sys.time(),
      R.cov.parameters = data.frame(variance = sigma2, trange = alpha),
      G.cov.parameters = data.frame(nu = nu),
      fixed.effects = fixed.effects.estimates,
      random.effects = gamma.hat,
      typeIII.hypoth = typeIII.hypoth,
      fit.table = fit.table,
      covb = covb,
      outer.iterations = iter,
      inner.iterations2 = inner.iter2
    )

	class(xout) <- "glmmLDTS"
	xout

}

