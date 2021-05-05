# This should give the same, or nearly the same, result as
# sum(dpois(x,y,log = TRUE)), except that terms that do not depend on
# the Poisson rates u are disregarded.
loglik_poisson <- function (x, y, e)
  return(sum(x*log(y + e) - y))

# TO DO: Explain here what this function does, and how to use it.
hpd <- function (x, conf.level = 0.95) {
  n <- length(x)
  m <- round(n*(1 - conf.level))
  x <- sort(x)
  y <- x[seq(n-m+1,n)] - x[seq(1,m)]
  i <- which.min(y)
  return(c(x[i],x[n-m+i]))
}

# TO DO: Explain here what this function does, and how to use it.
simulate_posterior_poisson <- function (x, L, f, ns = 1000, s = 1, e = 1e-15) {
  k <- length(f)
  r <- 0
  samples <- matrix(0,ns,k)
  for (i in 1:ns) {
    for (j in 1:k) {
      fnew    <- f
      fnew[j] <- exp(log(f[j]) + rnorm(1,sd = s))
      u       <- drop(L %*% f)
      unew    <- drop(L %*% fnew)
      # NOTE: gI think I need to fix this acceptance probability via
      # the change of variables formula to account for the fact that I
      # am sampling log(f), not f.
      a <- min(1,exp(loglik_poisson(x,unew,e) - loglik_poisson(x,u,e)))
      if (runif(1) < a) {
        f <- fnew
        r <- r + 1
      }
    }
    samples[i,] <- f
  }
  cat(sprintf("acceptance rate: %0.3f\n",r/(k*ns)))
  return(samples)
}
