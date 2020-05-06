# Fit model parameters p0,p1 of binomial model n1 ~ Binom(n,p), where
# the binomial success rates are p = q*p1 + (1-q)*p0, and n = n0 + n1.
# Input argument "e" is a small positive number added to the
# likelihood and gradient to avoid NaNs; specifically, logarithms of
# zero and division by zero.
#'
#' @importFrom stats optim
#' 
fit_binom_topic_model <- function (n0, n1, q, e = 1e-5) {

  # Make sure none of the "weights" are exactly 0 or 1.
  q <- pmax(q,e)
  q <- pmin(q,1 - e)

  # Compute the binomial success rates given the model parameters.
  pbinom <- function (p0, p1, q)
    q*p1 + (1-q)*p0
    
  # Define function for computing the negative log-likelihood at x,
  # where x = c(p0,p1).
  f <- function (x) {
    p <- pbinom(x[1],x[2],q)
    return(-(sum(n1*log(p+e) + n0*log(1-p+e))))
  }

  # Define function for computing the gradient of the negative
  # log-likelihood at x, where x = c(p0,p1).
  g <- function (x) {
    p <- pbinom(x[1],x[2],q)
    g <- n1/(p+e) - n0/(1-p+e)
    return(c(-sum(g*(1-q)),
             -sum(g*q)))
  }
  
  # Fit the binomial probabilities using the "limited-memory"
  # quasi-Newton method implemented in "optim".
  out <- optim(c(0.5,0.5),f,g,method = "L-BFGS-B",lower = c(0,0),
               upper = c(1,1),control = list(factr = 1e-8,maxit = 100))

  # Output MLEs of p0 and p1, and the other "optim" outputs.
  names(out$par) <- c("p0","p1")
  return(out)
}
