#------CREATE EVENLY SPACED CLASSES FROM NUMERIC VARIABLE
#------LABEL CLASSES WITH NUMERIC MIDPOINT OF CLASS

numeric.2.classes <- function(data.in, column.name, new.column.name,
  nclasses, n.digits)
{
  cat.range <- (max(data.in[,column.name]) -
    min(data.in[,column.name]))/nclasses + 1e-10
  n2c <- min(data.in[,column.name]) + cat.range *
    floor((data.in[,column.name] -
      min(data.in[,column.name]) + 1e-11)/cat.range) +
    cat.range/2
  n2c <- round(10^n.digits*n2c)/10^n.digits
  data.in[,new.column.name] <- as.factor(n2c)
  data.in
}

