# Compute maximum-likelihood estimates (MLEs) of the parameters in the
# "univariate" ("single-gene") Poisson model: x ~ Poisson(u), with
# Poisson rates u given by u = f0*(1-q) + f1*q. Parameters f0, f1 are
# estimated, and vector q is provided. The likelihood is maximized
# using optim; specifically, the low-memory BFGS quasi-Newton method
# is used.
#
# Input arguments p0 and p1 are initial estimates of the parameters.
# Input "control" is a list of control parameters passed to optim.
# Input argument "e" is a small positive scalar added to the
# likelihood and gradient to avoid NaNs; specifically, logarithms of
# zero and division by zero.
#
#' @importFrom stats optim
#' 
fit_poisson_optim <- function (x, q, f0 = 1, f1 = 1, e = 1e-15,
                               control = list(factr = 1e-14,maxit = 100)) {

  # Make sure none of the "weights" are exactly zero or exactly one.
  q <- pmax(q,e)
  q <- pmin(q,1 - e)

  # Returns the negative log-likelihood, where par = c(p0,p1).
  f <- function (par) {
    p <- pbinom(par[1],par[2],q)
    return(-loglik_binom(x,y,p,e))
  }
  
  # Returns the gradient of the negative log-likelihood, where par =
  # c(p0,p1).
  g <- function (par) {
    p <- pbinom(par[1],par[2],q)
    u <- x/(p+e) - y/(1-p+e)
    return(-c(sum(u*(1-q)),
              sum(u*q)))
  }
  
  # Fit the binomial probabilities using the "limited-memory"
  # quasi-Newton method implemented in "optim".
  out <- optim(c(p0,p1),f,g,method = "L-BFGS-B",lower = c(0,0),
               upper = c(1,1),control = control)

  # Output MLEs of p0 and p1, and the other "optim" outputs.
  names(out$par) <- c("p0","p1")
  return(out)
}

