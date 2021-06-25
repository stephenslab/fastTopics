# For each column of the counts matrix, compute MLEs of the parameters
# in the Poisson glm. Input argument "method" may be set to "glm", or
# to any valid setting for "method" in update_factors_poisson_nmf.
#
#' @importFrom stats glm.control
fit_poisson_models <- function (X, L, method, eps = 1e-15, numiter = 100,
                                tol = 1e-8, nc = 1) {

  # Get the number of topics (k) and the number of data columns (m).
  m <- ncol(X)
  k <- ncol(L)
  
  # For each column of X, compute MLEs of the Poisson model parameters.
  if (method == "glm") {
    control <- glm.control(epsilon = tol,maxit = numiter)
    F <- fit_poisson_models_glm(X,L,control)
  } else {
    control <- list(numiter = numiter,nc = nc,eps = eps)
    F <- matrix(0.5,m,k)
    F <- update_factors_poisson_nmf(X,F,L,1:m,method,control)
  }
  return(F)
}

# Implements fit_poisson_models for method = "glm".
fit_poisson_models_glm <- function (X, L, control)  {
  m <- ncol(X)
  k <- ncol(L)
  F <- matrix(0,m,k)
  for (i in 1:m)
    F[i,] <- fit_poisson_glm(X[,i],L,control = control)$coef
  return(F)
}

# Fit the Poisson glm, x[i] ~ Poisson(u[i]), using the "glm" function,
# in which the Poisson rates are u[i] = sum(L[i,]*f). Inputs x and L are
# the data: x is a vector of length n; and L is an n x k matrix. Input
# argument f is a vector of length k giving initial estimates of the
# coefficients.
#
#' @importFrom stats formula
#' @importFrom stats glm
#' @importFrom stats glm.control
#' @importFrom stats summary.glm
#' @importFrom stats poisson
#'
fit_poisson_glm <-
  function (x, L, f = rep(0.5,ncol(L)),
            control = glm.control(epsilon = 1e-10, maxit = 100)) {
  k <- length(f)
  dat <- as.data.frame(cbind(x,L))
  names(dat) <- c("x",paste0("f",1:k))
  model <- formula(paste("x ~",paste(names(dat)[-1],collapse = " + "),"- 1"))
  fit <- suppressWarnings(glm(model,family = poisson(link = "identity"),
                              data = dat,start = f,control = control))
  return(list(fit = fit,coef = summary.glm(fit)$coefficients[,"Estimate"]))
}

# Simulate draws from the posterior distribution of f via random-walk
# Metropolis on t = log(f). The posterior distribution is based on a
# uniform prior and the Poisson glm likelihood with identity link
# function. Input ns specifies the number of Monte Carlo samples to
# simulate. Input s determines the width (standard deviation) of the
# random walk, and input f is the initial state of the Markov
# chain. In the code below, t = log(f) is the current state of the
# Markov chain.
#
# The outputs are (1) "samples", an ns x k matrix of Monte Carlo
# samples of t = log(f), where k = length(f), and (2) "ar", the
# Metropolis acceptance rate.
#
# This is mainly used to test simulate_posterior_poisson_rcpp and is
# not actually used in practice because it is too slow.
#
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats dpois
#' 
simulate_posterior_poisson <- function (x, L, f, ns = 1000, s = 0.3,
                                        e = 1e-15) {
  k <- length(f)
  t <- log(f)
  ar <- 0
  samples <- matrix(0,ns,k)
  D <- matrix(rnorm(ns*k),ns,k)
  U <- matrix(runif(ns*k),ns,k)
  for (i in 1:ns) {
    for (j in 1:k) {

      # Randomly suggest moving to tj(new) = tj + d, where d ~ N(0,s).
      tnew    <- t
      d       <- s*D[i,j]
      tnew[j] <- t[j] + d

      # Compute the Metropolis acceptance probability, and move to the
      # new state according to this acceptance probability. Note that
      # the additional "d" in the acceptance probability is needed to
      # account for the fact that we are simulating log(f), not f; see
      # p. 11 of Devroye (1986) "Non-uniform random variate generation".
      u     <- drop(L %*% exp(t))
      unew  <- drop(L %*% exp(tnew))
      ll    <- sum(dpois(x,u + e,log = TRUE))
      llnew <- sum(dpois(x,unew + e,log = TRUE))
      a     <- exp((llnew - ll) + d)
      a     <- min(1,a)
      if (U[i,j] < a) {
        t  <- tnew
        ar <- ar + 1
      }
    }

    # Store the current state of the Markov chain.
    samples[i,] <- t
  }

  # Output the states of the Markov chain and the acceptance rate.
  return(list(samples = samples,ar = ar/(k*ns)))
}

# compute_lfc_le_old = function (lf) {
#   k <- length(lf)
#   b <- rep(0,k)
#   for (i in 1:k) {
#     x    <- lf[i] - lf
#     x[i] <- Inf
#     b[i] <- x[which.min(abs(x))]
#   }    
#   return(b)
# }

# TO DO: Explain here what this function does, and how (and when) to
# use it.
compute_lfc0 <- function (samples, f, f0, conf.level = 0.9) {
  k    <- length(f)
  low  <- rep(0,k)
  high <- rep(0,k)
  est  <- log(f) - log(f0)
  for (i in 1:k) {
    out     <- hpd(samples[,i] - log(f0),conf.level)
    low[i]  <- out[1]
    high[i] <- out[2]
  }
  return(rbind(est = est,low = low,high = high))
}

# TO DO: Explain here what this function does, and how (and when) to
# use it.
wcompute_lfc_pairwise <- function (samples, t, conf.level = 0.9) {
  # TO DO.
}

# TO DO: Explain here what this function does, and how (and when) to
# use it.
compute_lfc_le <- function (samples, t, conf.level = 0.9) {
  # TO DO.
}

# Return the covariance of t = log(f) in the Poisson glm, x[i] ~
# Poisson(u[i]), in which the Poisson rates are u[i] = L[i,]*f. The
# covariance calculations are based on a Laplace approximation to the
# likelihood. The input f should contain MLEs of the coefficients.
compute_poisson_covariance <- function (x, L, f) {
  u <- drop(L %*% f)
  return(chol2inv(chol(tcrossprod(f) * crossprod(sqrt(x)/u*L))))
}

# Add pseudocounts to the data for inference with the Poisson glm;
# specifically, add k rows of counts e to the data matrix, and replace
# the loadings matrix L with rbind(L,I), where I is the k x k identity
# matrix. This is used to extend ML estimation to MAP estimation for
# the Poisson glm, in which f[j] is assigned prior Gamma(a,b), with
# a=1+e, b=1, for j = 1,...,k.
add_pseudocounts <- function (X, L, e = 0.01) {
  m <- ncol(X)
  k <- ncol(L)
  return(list(X = rbind(X,matrix(e,k,m)),
              L = rbind(L,diag(k))))
}
