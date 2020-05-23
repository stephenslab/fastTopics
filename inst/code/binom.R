# Fit model parameters (p0, p1) of binomial model x ~ Binom(n,p), where
# the binomial success rates are p = q*p1 + (1-q)*p0, and n = x + y.
# Input argument "e" is a small positive number added to the
# likelihood and gradient to avoid NaNs; specifically, logarithms of
# zero and division by zero.
fit_binom_optim <- function (x, y, q, e = 1e-15) {

  # Make sure none of the "weights" are exactly 0 or 1.
  q <- pmax(q,e)
  q <- pmin(q,1 - e)

  # Compute the binomial success rates given the model parameters.
  pbinom <- function (p0, p1, q)
    q*p1 + (1-q)*p0
    
  # Define function for computing the negative log-likelihood at x,
  # where par = c(p0,p1).
  f <- function (par) {
    p <- pbinom(par[1],par[2],q)
    return(-(sum(x*log(p+e) + y*log(1-p+e))))
  }

  # Define function for computing the gradient of the negative
  # log-likelihood at par, where par = c(p0,p1).
  g <- function (par) {
    p <- pbinom(par[1],par[2],q)
    u <- x/(p+e) - y/(1-p+e)
    return(-c(sum(u*(1-q)),
              sum(u*q)))
  }
  
  # Fit the binomial probabilities using the "limited-memory"
  # quasi-Newton method implemented in "optim".
  out <- optim(c(0.5,0.5),f,g,method = "L-BFGS-B",lower = c(0,0),
               upper = c(1,1),control = list(factr = 1e-10,maxit = 40))

  # Output MLEs of p0 and p1, and the other "optim" outputs.
  names(out$par) <- c("p0","p1")
  return(out$par)
}

# TO DO: Explain here what this function does, and how to use it.
fit_binom_em <- function (x, y, q) {

}

  # Compute the log-ratio of the likelihoods.
  #
  # TO DO: Explain this code in greater detail.
  #
  # This is the marginal probability density of x ~ Binom(n,p), with
  # success rate p drawn uniformly over the [0,1] interval. Note that
  # the density function is the same regardless of the number of
  # successes, n1.
  ## n     <- n0 + n1
  ## ll0   <- (lgamma(sum(n) + 1) - lgamma(sum(n) + 2))
  ## llmle <- -out$value
  ## loglik <- function (x, y)
  ##   -f(c(x,y))
  ## quadf <- function (X, Y) {
  ##   Z <- X
  ##   for (i in 1:length(X))
  ##     Z[i] <- loglik(X[i],Y[i])
  ##   return(exp(Z - llmle))
  ## }
  ## out$loglikratio <- sum(lchoose(n,n1)) + log(quad2d(quadf,0,1,0,1)) +
  ##                    llmle - ll0
  
