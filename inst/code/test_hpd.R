hpd <- function (x, conf = 0.95) {
  n <- length(x)
  m <- round(n*(1 - conf))
  x <- sort(x)
  y <- x[seq(n-m+1,n)] - x[seq(1,m)]
  i <- which.min(y)
  return(c(x[i],x[n-m+i]))
}

set.seed(1)
n <- 1e4
x <- rnorm(n)
print(hpd(x))
