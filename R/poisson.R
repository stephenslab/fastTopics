# Compute the Poisson rates given the size factors (s), topic
# proportions (q), and parameters (f0, f1) of the Poisson model.
get_poisson_rates <- function (s, q, f0, f1)
  s*((1-q)*f0 + q*f1)

# This should give the same, or nearly the same, result as
# sum(dpois(x,u,log = TRUE)), except that terms that do not depend on
# the Poisson rates u are discarded.
loglik_poisson <- function (x, u, e = 1e-15)
  return(sum(x*log(u+e) - u))

# Compute maximum-likelihood estimates (MLEs) of the parameters in the
# single-count Poisson model: x ~ Poisson(s*u), with u given by u =
# f0*(1-q) + f1*q. Parameters f0, f1 are estimated, and vectors s, q
# are provided. The likelihood is maximized by optim with method =
# "L-BFGS-B".
#
# Input arguments f0 and f1 are initial estimates of the parameters.
# Input "control" is a list of control parameters passed to optim.
# Input argument "e" is a small positive scalar added to the
# likelihood and gradient to avoid logarithms of zero and division by
# zero.
#
#' @importFrom stats optim
#' 
fit_poisson_optim <- function (x, s, q, f0 = 1, f1 = 1, e = 1e-15,
                               control = list(factr = 1e-14, maxit = 100)) {

  # Make sure none of the topic proportions are exactly zero or
  # exactly one.
  q <- pmax(q,e)
  q <- pmin(q,1-e)

  # This is used to computes the negative log-likelihood, in which par
  # = c(f0,f1).
  f <- function (par) {
    u <- get_poisson_rates(s,q,par[1],par[2])
    return(-loglik_poisson(x,u,e))
  }
  
  # Returns the gradient of the negative log-likelihood, in which par
  # = c(f0,f1).
  g <- function (par) {
    u <- get_poisson_rates(s,q,par[1],par[2])
    y <- (u - x)/(u + e)
    return(c(sum(y*(1-q)),sum(y*q)))
  }
  
  # Fit the Poisson rates f0, f1 using the limited-memory quasi-Newton
  # method implemented in optim.
  out <- optim(c(f0,f1),f,g,method = "L-BFGS-B",lower = c(0,0),
               upper = c(Inf,Inf),control = control)

  # Output the MLEs of f0 and f1, and the other "optim" outputs.
  names(out$par) <- c("f0","f1")
  return(out)
}

# TO DO: Revise this description.
#
# Perform EM updates to compute maximum-likelihood estimates (MLEs) of
# parameters (p0, p1) in the binomial topic model. Specifically, the
# binomial topic model is x ~ Binom(n,p), where the binomial success
# rates are p = q*p1 + (1-q)*p0, and n = x + y. Input argument "e" is
# a small positive number added to the likelihood and gradient to
# avoid NaNs; specifically, logarithms of zero and division by zero.
fit_poisson_em <- function (x, s, q, f0 = 1, f1 = 1, numiter = 40, e = 1e-15) {

  # Make sure none of the topic proportions are exactly zero or
  # exactly one.
  q <- pmax(q,e)
  q <- pmin(q,1-e)

  # Precalculations.
  w0 <- s*(1-q)
  w1 <- s*q
  
  # Monitor progress by computing the log-likelihood at each iteration.
  loglik <- rep(0,numiter)
  for (iter in 1:numiter) {

    # E-STEP
    # ------
    u  <- w0*f0 + w1*f1 + e 
    z0 <- x*w0*f0 + e
    z1 <- x*w1*f1 + e
    z0 <- z0/u
    z1 <- z1/u
    
    # M-STEP
    # ------
    f0 <- sum(z0)/sum(w0)
    f1 <- sum(z1)/sum(w1)
    
    # Compute the log-likelihood at the current estimates of the model
    # parameters (ignoring terms that don't depend on f0 or f1).
    u            <- get_poisson_rates(s,q,f0,f1)
    loglik[iter] <- loglik_poisson(x,u,e)
  }

  # Output the estimates of f0 and f1, and the log-likelihood at each EM
  # iteration.
  f        <- c(f0,f1)
  names(f) <- c("f0","f1")
  return(list(f = f,loglik = loglik))
}
