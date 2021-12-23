#------------ FUNCTION FOR COMPUTING TYPE III GENERAL LINEAR HYPOTHESES

hypoth <- function(formula, data, parms, grp.vec, timevec,
  Z, yt, A.5, Del.i, n.per.term)
{
	termlabs <- attr(terms(formula),"term.labels")
	if(length(termlabs) == 0) return(NULL)
	z <- model.frame(formula, data)
	z <- data[,names(z)[1]]
	z.col <- as.character(attr(terms(formula, data = data),"variables"))[2]
	form1 <- termlabs[1]
	if(length(termlabs) > 1) for(i in 2:length(termlabs)) form1 <- paste(form1, "+", termlabs[i])
	formula <- as.formula(paste(z.col, "~", form1))
	X <- X.design(formula = formula, data = data, Xmethod = "sum.to.0")
	terms.table <- attr(terms(formula),"factors")
	n <- length(z)
	p <- length(X[1,])
  bgh <-  beta.gamma.hat(parms, grp.vec, timevec,
    X, Z, yt, A.5, Del.i, n.per.term)
	sigma2 <- bgh$qQq/(n - p)
  covb.sum0 <- bgh$covb
  covb.sum0 <- as.numeric(sigma2)*covb.sum0
	b.sum0 <- bgh$beta.hat
	nterms <- length(terms.table[1,])
	cumcol <- 1
	storage <- matrix(0, ncol = 4, nrow = nterms)
	nfactor <- 0
	for(i in 1:nterms) {
		facts <- rownames(terms.table)[terms.table[,i] > 0]
		nfacts <- length(facts)
		nxcol <- 1
		for(j in 1:nfacts)
			if(is.factor(data[,facts[j]]))
				nxcol <- nxcol*(length(levels(data[,facts[j]])) - 1)
		K <- matrix(0, nrow = p, ncol = nxcol)
		K[(cumcol + 1):(cumcol + nxcol), 1:nxcol] <- diag(nxcol)
		cumcol <- cumcol + nxcol
		storage[i,1] <- nxcol
		storage[i,2] <- n - p
		storage[i,3] <- t(t(K) %*% b.sum0) %*% mginv(t(K) %*% covb.sum0 %*% K) %*%
			(t(K) %*% b.sum0)/sum(svd(K)$d>1e-10)
		storage[i,4] <- round(100000*(1 - pf(storage[i,3],storage[i,1], storage[i,2])))/100000
	}
	storage <- as.data.frame(storage)
	colnames(storage) <- c("Num.df", "Den.df", "F.value", "Prob.F")
	outpt <- data.frame(Effect = termlabs, storage)
	outpt
}

