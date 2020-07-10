# Compute the binomial success rates given the model parameters.
get_binom_probs <- function (p0, p1, q)
  (1-q)*p0 + q*p1
    
# This should give the same result as dbinom(x,s,p), except that terms
# that do not depend on the success probabilities p are ignored.
loglik_binom <- function (x, s, p, e = 1e-15)
  return(sum(x*log(p+e) + (s-x)*log(1-p+e)))

# Compute maximum-likelihood estimates (MLEs) of the parameters in the
# single-count binomial model: x ~ Binom(s,p), with binomial success
# rates p = p0*(1-q) + p1*q. Probabilities p0, p1 are estimated, and
# vectors s and q are provided. The likelihood is maximized using
# optim with method = "L-BFGS-B". # Input arguments p0 and p1 are
# initial estimates of the parameters.  Input "control" is a list of
# control parameters passed to optim.  Input argument "e" is a small
# positive scalar added to the likelihood and gradient to avoid NaNs;
# specifically, logarithms of zero and division by zero.
#
#' @importFrom stats optim
#' 
fit_binom_optim <- function (x, s, q, p0 = 0.5, p1 = 0.5, e = 1e-15,
                             control = list(factr = 1e-14, maxit = 100)) {

  # Make sure none of the "weights" are exactly zero or exactly one.
  q <- pmax(q,e)
  q <- pmin(q,1 - e)

  # Returns the negative log-likelihood, where par = c(p0,p1).
  f <- function (par) {
    p <- get_binom_probs(par[1],par[2],q)
    return(-loglik_binom(x,s,p,e))
  }
  
  # Returns the gradient of the negative log-likelihood, where par =
  # c(p0,p1).
  g <- function (par) {
    p <- get_binom_probs(par[1],par[2],q)
    u <- (s-x)/(1-p+e) - x/(p+e)
    return(c(sum(u*(1-q)),sum(u*q)))
  }
  
  # Fit the binomial probabilities using the limited-memory
  # quasi-Newton method implemented in optim.
  out <- optim(c(p0,p1),f,g,method = "L-BFGS-B",lower = c(0,0),
               upper = c(1,1),control = control)

  # Output MLEs of p0 and p1, and the other optim outputs.
  names(out$par) <- c("p0","p1")
  return(out)
}

# Perform EM updates to compute maximum-likelihood estimates (MLEs) of
# parameters (p0, p1) in the binomial topic model. Specifically, the
# binomial topic model is x ~ Binom(s,p), where the binomial success
# rates are p = (1-q)*p0 + q*p1. Input argument "e" is a small
# positive number added to the likelihood and gradient to avoid NaNs;
# specifically, logarithms of zero and division by zero.
fit_binom_em <- function (x, s, q, p0 = 0.5, p1 = 0.5, numiter = 40,
                          e = 1e-15) {

  # Monitor progress by computing the log-likelihood at each iteration.
  loglik <- rep(0,numiter)
  for (iter in 1:numiter) {

    # E-step
    # ------
    # Compute the posterior probabilities,
    #
    #   p00 = Pr(z = 0 | x = 0)
    #   p10 = Pr(z = 0 | x = 1)
    #   p01 = Pr(z = 1 | x = 0)
    #   p11 = Pr(z = 1 | x = 1)
    # 
    p00 <- (1-p0)*(1-q)
    p01 <- (1-p1)*q
    p10 <- p0*(1-q)
    p11 <- p1*q
    p00 <- p00/(p00 + p01)
    p11 <- p11/(p10 + p11)
    p10 <- 1 - p11
    p01 <- 1 - p00
    
    # M-step
    # ------
    # Update the binomial topic model parameters,
    # p0 = Pr(x = 1 | z = 0) and p1 = Pr(x = 1 | z = 1).
    p0 <- sum(x*p10)/sum((s-x)*p00)
    p1 <- sum(x*p11)/sum((s-x)*p01)
    p0 <- p0/(1 + p0)
    p1 <- p1/(1 + p1)
    
    # Compute the log-likelihood at the current estimates of the model
    # parameters (ignoring terms that don't depend on p0 or p1).
    p            <- get_binom_probs(p0,p1,q)
    loglik[iter] <- loglik_binom(x,s-x,p,e)
  }

  # Output the estimates of p0 and p1, and the log-likelihood at each EM
  # iteration.
  p        <- c(p0,p1)
  names(p) <- c("p0","p1")
  return(list(p = p,loglik = loglik))
}
