# ------- A FUNCTION TO CREATE A SAS-LIKE TABLE OF FIXED EFFECTS

fixed.effect.table <-
function(formula, data, fixed.effects.estimates)
{
  fee <- fixed.effects.estimates
  form1 <- formula
  fee.i <- 0
  fee.xtra <- NULL
  levels.lab <- NULL
  if(attr(terms(form1),"intercept") == 1) {
    levels.lab <- data.frame(
      factor = "intercept",
      levels = "",
      fee[1,]
     )
     fee.i <- 1
  }
  terms.table <- attr(terms(form1),"factors")
  if(length(attr(terms(form1),"term.labels")) > 0) {
    if(sum(terms.table[,1]) == 1 & attr(terms(form1),"intercept") == 0)
      terms.table[terms.table[,1] == 1,1] <- 2
    for (j in 1:length(terms.table[1,])) {
      ind.j <- terms.table[,j] > 0
      n.interact <- sum(ind.j)
      lab.i <- NULL
      for (i in 1:n.interact) {
        if(!any(rownames(terms.table)[ind.j][i] == colnames(data))) {
          levels.i <- matrix("")
          if(i == 1) lab.i <- matrix(levels.i, ncol = 1)
          if(i > 1) {
            lab.t <- NULL
            for (k in 2:i) lab.t <- cbind(lab.t,
              matrix(rep(lab.i[,k-1], times = length(levels.i)), ncol = 1))
              lab.tmp <- matrix(rep(levels.i, each = length(lab.i[,1])), ncol = 1)
              lab.i <- cbind(lab.t,lab.tmp)
          }
        }
        else if(is.factor(data[,rownames(terms.table)[ind.j][i]])) {
          levels.i <- levels(data[,rownames(terms.table)[ind.j][i]])
          if(i == 1) lab.i <- matrix(levels.i, ncol = 1)
          if(i > 1) {
            lab.t <- NULL
            for (k in 2:i) lab.t <- cbind(lab.t,
              matrix(rep(lab.i[,k-1], times = length(levels.i)), ncol = 1))
              lab.tmp <- matrix(rep(levels.i, each = length(lab.i[,1])), ncol = 1)
              lab.i <- cbind(lab.t,lab.tmp)
          }
        }
        else {
          levels.i <- matrix("")
          if(i == 1) lab.i <- matrix(levels.i, ncol = 1)
          if(i > 1) {
            lab.t <- NULL
            for (k in 2:i) lab.t <- cbind(lab.t,
              matrix(rep(lab.i[,k-1], times = length(levels.i)), ncol = 1))
              lab.tmp <- matrix(rep(levels.i, each = length(lab.i[,1])), ncol = 1)
              lab.i <- cbind(lab.t,lab.tmp)
          }
        }
      }
      cum.ind <- rep(TRUE, times = length(lab.i[,1]))
      levels.lab.i <- NULL
      for (i in 1:n.interact) {
        if(!any(rownames(terms.table)[ind.j][i] == colnames(data)))
          cum.ind <- cum.ind & TRUE
        else if(is.factor(data[,rownames(terms.table)[ind.j][i]])){
          if (terms.table[ind.j,j][i] == 1) {
            cum.ind <- cum.ind & lab.i[,i] !=
            levels(data[,rownames(terms.table)[ind.j][i]])[1]
          }
        }
        else cum.ind <- cum.ind & TRUE
        if(i == 1) levels.lab.i <- lab.i[,1]
        if(i > 1) levels.lab.i <- paste(levels.lab.i, ",", lab.i[,i], sep = "")
      }
      n.fee <- sum(cum.ind)
      n.j <- length(cum.ind)
      fee.xtra <- matrix(c(0,NA,NA,NA,NA), nrow = n.j, ncol = 5, byrow = TRUE)
      fee.xtra[cum.ind,1:5] <- matrix(unlist(fee[(fee.i + 1):(fee.i + n.fee), 1:5]), ncol = 5)
      fee.xtra1 <- data.frame(factor = colnames(terms.table)[j], levels = levels.lab.i,
        parameters = fee.xtra[,1], std.err = fee.xtra[,2],
        df = fee.xtra[,3], t.value = fee.xtra[,4],
        prob.t = fee.xtra[,5])
      if(is.null(levels.lab)) levels.lab <- fee.xtra1
      else levels.lab <- rbind(levels.lab, fee.xtra1)
      fee.i <- fee.i + n.fee
    }
  }
  rownames(levels.lab) <- 1:length(levels.lab[,1])
  colnames(levels.lab) <- c("effect", "levels", "estimate", "std.err", "df",
    "t.value", "prob.t")
  levels.lab
}


