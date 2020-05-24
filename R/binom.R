#' @title Function Title Goes Here
#'
#' @description Description of function goes here.
#'
#' @export
#' 
binom_topic_analysis <- function (X, fit, verbose = TRUE) {
  # TO DO.
}

# TO DO: Explain here what this function does, and how to use it.
multinom2binom <- function (X, fit, e = 0.01, version = c("Rcpp","R")) {
    
  # Compute the expected counts.
  if (version == "Rcpp")
    out <- binom_stats_rcpp(X,fit$F,fit$L)
  else if (version == "R")
    out <- binom_stats(X,fit$F,fit$L)
  a0 <- out$a0
  a1 <- out$a1
  b0 <- drop(out$b0)
  b1 <- drop(out$b1)

  # Compute parameters p0 and p1 in the binomial topic model, for each
  # gene j (row of X), and for each topic k.
  P0 <- (a0 + e)/(b0 + e)
  P1 <- (a1 + e)/(b1 + e)
  rownames(P0) <- rownames(fit$F)
  rownames(P1) <- rownames(fit$F)
  colnames(P0) <- colnames(fit$F)
  colnames(P1) <- colnames(fit$F)
  return(list(P0 = P0,P1 = P1))
}

#' @title Function Title Goes Here
#'
#' @description Description of function goes here.
#'
#' @details Using this function requires some care; only minimal
#' argument checking is performed, and error messages may not be
#' helpful.
#'
#' @param X Describe input argument X.
#'
#' @param F Describe input argument F.
#'
#' @param method Describe input argument "method".
#'
#' @param e Describe input argument "e".
#'
#' @param control Describe input argument "control".
#'
#' @param verbose Describe input argument "verbose".
#'
#' @return Describe return value here.
#'
#' @importFrom progress progress_bar
#' 
#' @export
#' 
fit_binom_topic_model <- function (X, fit, P0, P1, method = c("em", "optim"),
                                   e = 1e-15,
                                   control = list(factr = 1e-14,maxit = 100),
                                   verbose = TRUE) {

  # Verify and process the input arguments.
  #
  # TO DO.
  #
  method <- match.arg(method)
  
  # Get the number of rows (n) and columns (m) of X, and get the
  # number of topics (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(fit$F)

  # Initialize the binomial topic model parameter estimates.
  if (missing(P0))
    P0 <- matrix(0.5,m,k)
  if (missing(P1))
    P1 <- matrix(0.5,m,k)

  # Fit a binomial topic model for each column (gene) j and topic k.
  s <- rowSums(X)
  if (verbose)
  pb <- progress_bar$new(total = m)
  for (i in 1:m) {
    if (verbose)
      pb$tick()
    x <- X[,i]
    for (j in 1:k) {
      q <- fit$L[,j]
      if (method == "em") {
        # TO DO.
      } else if (method == "optim") {
        out     <- fit_binom_optim(x,s - x,q,P0[i,j],P1[i,j],e,control)
        P0[i,j] <- out$par["p0"]
        P1[i,j] <- out$par["p1"]
      } 
    }
  }

  # Output the maximum-likelihood estimates (MLEs) of the binomial
  # topic model for each column (gene) j, and for each topic k.
  return(list(P0 = P0,P1 = P1))
}

# TO DO: Explain here what this function does, and how to use it.
binom_stats <- function (X, F, L) {

  # Get the number of rows (n) and columns (m) of the counts matrix,
  # and the number of topics (k).
  n <- nrow(X)  
  m <- ncol(X)
  k <- ncol(F)

  # Initialize storage for expectations and marginal probabilities.
  a0 <- matrix(0,m,k)
  a1 <- matrix(0,m,k)
  b0 <- rep(0,k)
  b1 <- rep(0,k)
  
  # Repeat for row (i) and column (j) of the counts matrix.
  for (i in 1:n)
    for (j in 1:m) {
      x <- X[i,j]
        
      # Compute the posterior topic assignment probabilities.
      p <- F[j,] * L[i,]
      p <- p / sum(p)

      # Update the expectations.
      a0[j,] <- a0[j,] + x*(1-p)
      a1[j,] <- a1[j,] + x*p
      b0     <- b0     + x*(1-p)
      b1     <- b1     + x*p
    }

  # Output the expectations (a0, a1) and marginal probabilities (b0,
  # b1).
  return(list(a0 = a0,
              a1 = a1,
              b0 = b0,
              b1 = b1))
}

# Compute maximum-likelihood estimates (MLEs) of parameters p0, p1 in
# the binomial topic model using optim; specifically, the low-memory
# BFGS quasi-Newton method. The binomial topic model is x ~
# Binom(n,p), where the binomial success rates are p = q*p1 +
# (1-q)*p0, and n = x + y. Inputs p0 and p1 are initial parameter
# estimates. Input "control" is a list of control parameters passed to
# optim. Input argument "e" is a small positive scalar added to the
# likelihood and gradient to avoid NaNs; specifically, logarithms of
# zero and division by zero.
#
#' @importFrom stats optim
#' 
fit_binom_optim <- function (x, y, q, p0 = 0.5, p1 = 0.5, e = 1e-15,
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

# This return the same value as dbinom(x,x+y,p), except that terms
# that do not depend on the success probabilities p are ignored.
loglik_binom <- function (x, y, p, e = 1e-15)
  return(sum(x*log(p+e) + y*log(1-p+e)))

