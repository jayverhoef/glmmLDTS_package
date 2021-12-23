# -------PUT BETA HAT ESTIMATES BACK ON ORIGINAL UNSTANDARDIZED SCALE

unstdX.beta.hat <- function(X, beta.hat, covb) {
	bhat.se <- beta.hat*0
	B <- matrix(0, ncol = length(X[1,]), nrow =length(X[1,]))
	B[1,1] <- 1
	for(i in 2:length(X[1,])) {
		X.std.means <- mean(X[,i])
		X.std.stdev <- sqrt(var(X[,i]))
		B[i,i] <- 1/X.std.stdev
		B[1,i] <- -X.std.means/X.std.stdev
	}
	beta.hat1 <- B %*% beta.hat
	covb1 <- B %*% covb %*% t(B)
	list(beta.hat = beta.hat1, covb = covb1)
}

