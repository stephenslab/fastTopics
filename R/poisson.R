# Compute the Poisson rates given the size factors (s), topic
# proportions (q), and parameters (f0, f1) of the Poisson model.
get_poisson_rates <- function (s, q, f0, f1)
  s*((1-q)*f0 + q*f1)

# This should give the same, or nearly the same, result as
# sum(dpois(x,y,log = TRUE)), except that terms that do not depend on
# the Poisson rates u are disregarded.
loglik_poisson <- function (x, y, e = 1e-15)
  return(sum(x*log(y+e) - y))

# For each column of the counts matrix, and for each topic, compute
# maximum-likelihood estimates (MLEs) of the parameters in the
# single-count (univariate) Poisson model.
fit_univar_poisson_models <-
  function (X, L, s = rep(1,nrow(X)), method = c("em-rcpp","em","optim"),
            e = 1e-15, numiter.em = 40,
            control.optim = list(factr = 1e-14,maxit = 100), verbose = TRUE) {
  method <- match.arg(method)

  # Get the number of topics (k) and the number of columns of X (m).
  m <- ncol(X)
  k <- ncol(L)
  
  # Make sure none of the topic proportions are exactly zero or
  # exactly one.
  L <- pmax(L,e)
  L <- pmin(L,1-e)

  # Compute maximum-likelihood estimates (MLEs) of the Poisson model
  # parameters (f0,f1) for each topic, and for each column of of the
  # counts matrix, X.
  if (method == "optim")
    out <- fit_univar_poisson_models_optim(X,L,s,e,control.optim,verbose)
  else if (method == "em") {
    out <- fit_univar_poisson_models_em(X,L,s,e,numiter.em,verbose)
  } else if (method == "em-rcpp") {
    if (is.sparse.matrix(X)) {
      # TO DO.
    } else {
      # TO DO.
    }
  }
  return(out)
}

# Implements fit_univar_poisson_models for method = "optim".
#
#' @importFrom progress progress_bar
fit_univar_poisson_models_optim <- function (X, L, s, e, control, verbose)  {
  m      <- ncol(X)
  k      <- ncol(L)
  F0     <- matrix(0,m,k)
  F1     <- matrix(0,m,k)
  loglik <- matrix(0,m,k)
  if (verbose)
    pb <- progress_bar$new(total = m)
  for (i in 1:m) {
    if (verbose)
      pb$tick()
    for (j in 1:k) {
      out         <- fit_poisson_optim(X[,i],s,L[,j],e = e,control = control)
      F0[i,j]     <- out$par["f0"]
      F1[i,j]     <- out$par["f1"]
      loglik[i,j] <- -out$value
    }
  }
  return(list(F0 = F0,F1 = F1,loglik = loglik))
}

# Implements fit_univar_poisson_models for method = "em".
#
#' @importFrom progress progress_bar
fit_univar_poisson_models_em <- function (X, L, s, e, numiter, verbose) {
  m      <- ncol(X)
  k      <- ncol(L)
  F0     <- matrix(0,m,k)
  F1     <- matrix(0,m,k)
  loglik <- matrix(0,m,k)
  if (verbose)
    pb <- progress_bar$new(total = m)
  for (i in 1:m) {
    if (verbose)
      pb$tick()
    for (j in 1:k) {
      out         <- fit_poisson_em(X[,i],s,L[,j],e = e,numiter = numiter)
      F0[i,j]     <- out$f["f0"]
      F1[i,j]     <- out$f["f1"]
      loglik[i,j] <- max(out$loglik)
    }
  }
  return(list(F0 = F0,F1 = F1,loglik = loglik))
}

# Compute maximum-likelihood estimates (MLEs) of the parameters in the
# single-count Poisson model: x ~ Poisson(s*u), with u given by u =
# f0*(1-q) + f1*q. Parameters f0, f1 are estimated, and vectors s, q
# are provided. The likelihood is maximized by optim with method =
# "L-BFGS-B". Input arguments f0 and f1 are initial estimates of the
# parameters. Input "control" is a list of control parameters passed
# to optim. Input argument "e" is a small positive scalar added to
# the likelihood and gradient to avoid logarithms of zero and division
# by zero.
#
#' @importFrom stats optim
#' 
fit_poisson_optim <- function (x, s, q, f0 = 1, f1 = 1, e = 1e-15,
                               control = list(factr = 1e-14, maxit = 100)) {

  # This is used to computes the negative log-likelihood, in which
  # par = c(f0,f1).
  f <- function (par) {
    su <- get_poisson_rates(s,q,par[1],par[2])
    return(-loglik_poisson(x,su,e))
  }
  
  # Returns the gradient of the negative log-likelihood, in which
  # par = c(f0,f1).
  g <- function (par) {
    u <- get_poisson_rates(s,q,par[1],par[2])/s
    y <- (s*u - x)/(u + e)
    return(c(sum(y*(1-q)),sum(y*q)))
  }
  
  # Fit the Poisson rates f0, f1 using the L-BFGS-B quasi-Newton method
  # implemented in optim.
  out <- optim(c(f0,f1),f,g,method = "L-BFGS-B",lower = c(0,0),
               control = control)

  # Output the MLEs of f0 and f1, and the other "optim" outputs.
  out$par        <- out$par
  names(out$par) <- c("f0","f1")
  return(out)
}

# Perform EM updates to compute maximum-likelihood estimates (MLEs) of
# the parameters in the single-count Poisson model: x ~ Poisson(s*u),
# with u given by u = f0*(1-q) + f1*q. Parameters f0, f1 are
# estimated, and vectors s, q are provided. Input arguments f0 and f1
# are initial estimates of the parameters. Input argument "e" is a
# small positive scalar added to the denominator in the E-step to
# avoid division by zero.
fit_poisson_em <- function (x, s, q, f0 = 1, f1 = 1, e = 1e-15, numiter = 40) {

  # Store a couple pre-calculations to simplify the calculations
  # below.
  a <- sum(s*(1-q))
  b <- sum(s*q)
  
  # Perform the EM updates. Progress is monitored by computing the
  # log-likelihood at each iteration.
  loglik <- rep(0,numiter)
  for (iter in 1:numiter) {

    # E-STEP
    # ------
    z0 <- f0*(1-q)
    z1 <- f1*q
    u  <- z0 + z1 + e
    z0 <- x*z0/u
    z1 <- x*z1/u
    
    # M-STEP
    # ------
    f0 <- sum(z0)/a
    f1 <- sum(z1)/b
    
    # Compute the log-likelihood at the current estimates of the model
    # parameters (ignoring terms that don't depend on f0 or f1).
    su           <- get_poisson_rates(s,q,f0,f1)
    loglik[iter] <- loglik_poisson(x,su,e)
  }

  # Output the estimates of f0 and f1, and the log-likelihood at each EM
  # iteration.
  f        <- c(f0,f1)
  names(f) <- c("f0","f1")
  return(list(f = f,loglik = loglik))
}

