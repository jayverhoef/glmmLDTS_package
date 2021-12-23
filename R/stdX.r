# ----------- STANDARDIZE THE COLUMNS OF THE DESIGN MATRIX

stdX <- function(X) {
	X.std <- X
	for(i in 2:length(X[1,])) {
		X.std.means <- mean(X[,i])
		X.std.stdev <- sqrt(var(X[,i]))
		X.std[,i] <- (X[,i] - X.std.means)/X.std.stdev
	}
	X.std
}

