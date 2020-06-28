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

  # Make sure all the topic proportions are positive.
  q <- pmax(q,e)

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

