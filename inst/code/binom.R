# Compute the binomial success rates given the model parameters.
pbinom <- function (p0, p1, q)
  q*p1 + (1-q)*p0
    
# Compute maximum-likelihood estimates (MLEs) of parameters (p0, p1)
# in the binomial topic model. The binomial topic model is x ~
# Binom(n,p), where the binomial success rates are p = q*p1 +
# (1-q)*p0, and n = x + y. Input argument "e" is a small positive
# number added to the likelihood and gradient to avoid NaNs;
# specifically, logarithms of zero and division by zero.
fit_binom_optim <- function (x, y, q, e = 1e-15) {

  # Make sure none of the "weights" are exactly 0 or 1.
  q <- pmax(q,e)
  q <- pmin(q,1 - e)

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
               upper = c(1,1),control = list(factr = 1e-12,maxit = 100))

  # Output MLEs of p0 and p1, and the other "optim" outputs.
  names(out$par) <- c("p0","p1")
  return(out)
}

# Perform EM updates to compute maximum-likelihood estimates (MLEs) of
# parameters (p0, p1) in the binomial topic model. Specifically, the
# binomial topic model is x ~ Binom(n,p), where the binomial success
# rates are p = q*p1 + (1-q)*p0, and n = x + y. Input argument "e" is
# a small positive number added to the likelihood and gradient to
# avoid NaNs; specifically, logarithms of zero and division by zero.
fit_binom_em <- function (x, y, q, p0 = 0.5, p1 = 0.5, numiter = 40,
                          e = 1e-15) {
  loglik <- rep(0,numiter)
  for (iter in 1:numiter) {

    # E-step
    # ------
    p00 <- (1-p0)*(1-q)
    p01 <- (1-p1)*q
    p10 <- p0*(1-q)
    p11 <- p1*q
    p00 <- p00/(p00 + p01)
    p11 <- p11/(p10 + p11)
    
    # M-step
    # ------
    p0 <- sum(x*(1 - p11))/sum(y*p00)
    p1 <- sum(x*p11)/sum(y*(1 - p00))
    p0 <- p0/(1 + p0)
    p1 <- p1/(1 + p1)
    
    # Compute the log-likelihood at the current estimates of the model
    # parameters.
    p            <- pbinom(p0,p1,q)
    loglik[iter] <- sum(x*log(p+e) + y*log(1-p+e))
  }

  # Output the MLEs of p0 and p1, and the log-likelihood at each EM
  # iteration.
  p <- c(p0,p1)
  names(p) <- c("p0","p1")
  return(list(p = p,loglik = loglik))
}

compute_lfoldchange_helper <- function (X, F, L, k) {

  # Get the number of rows (n) and columns (m) of the counts matrix,
  # and the number of topics (k).
  n <- nrow(X)  
  m <- ncol(X)

  # Initialize storage for the expected counts.
  n0 <- rep(0,m)
  n1 <- rep(0,m)
  
  # Repeat for row (i) and column (j) of the counts matrix.
  for (i in 1:n)
    for (j in 1:m) {
      x <- X[i,j]
        
      # Compute the posterior topic assignment probabilities.
      p <- F[j,] * L[i,]
      p <- p / sum(p)

      # Update the expectations.
      n0[j] <- n0[j] + x*(1-p[k])
      n1[j] <- n1[j] + x*p[k]
    }

  # Output the expectations.
  return(list(n0 = n0,n1 = n1))
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
  
