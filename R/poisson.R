# For each column of the counts matrix, compute maximum-likelihood
# estimates (MLEs) of the parameters in the Poisson glm. Input
# argument "method" may be set to "glm", or to any valid setting for
# "method" in update_factors_poisson_nmf.
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

# Implements fit_poisson_models for method = "glm".
fit_poisson_models_glm <- function (X, L, control)  {
  m <- ncol(X)
  k <- ncol(L)
  F <- matrix(0,m,k)
  for (i in 1:m)
    F[i,] <- fit_poisson_glm(X[,i],L,control = control)$coef
  return(F)
}

# Return log-fold change statistics lfc[j,k], and accompanying
# standard errors se[j,k] and z-scores z[j,k], under the Poisson glm
# x[i,j] ~ Pois(u[i,j]), in which the Poisson rates are u[i,j] =
# sum(L[i,]*F[j,]). The LFC calculations are determined by the "stat"
# input argument: when stat = "vsmean", lfc[j,k] = log2(F[j,k]/mu[j]);
# when stat = "le", lfc[j,k] is defined as the "least extreme"
# pairwise log-fold change; and when stat = k', lfc[j,k] is the
# pairwise log-fold change statistic lfc[j,k] = log2(F[j,k]/F[j,k']).
#
# TO DO: Revise this function.
# 
compute_lfc_stats <- function (X, F, L, mu, stat = "vsmean",
                               version = c("Rcpp","R")) {
  if (version == "Rcpp") {
    # TO DO.
  } else {
    m   <- ncol(X)
    k   <- ncol(L)
    lfc <- matrix(0,m,k)
    se  <- matrix(0,m,k)
    z   <- matrix(0,m,k)
    for (i in 1:m) {
      out     <- compute_lfc_stats_helper(X[,i],L,F[i,],mu[i],stat)
      lfc[i,] <- out$lfc
      se[i,]  <- out$se
      z[i,]   <- out$z
    }
  }
  return(list(lfc = lfc,se = se,z = z))
}

# Return log-fold change statistics lfc[j], and their standard errors
# se[j] and z-scores z[j], under the Poisson glm x[i] ~ Pois(u[i]), in
# which the Poisson rates are u[i] = sum(L[i,]*f). The LFC
# calculations are determined by the "stat" input argument: when stat
# = "vsmean", lfc[j] = log2(f[j]/mu); when stat = "le", lfc[j] is
# defined as the "least extreme" pairwise log-fold change; and when
# stat = k, lfc[j] is the pairwise log-fold change statistic lfc[j] =
# log2(f[j]/f[k]).
#
# TO DO: Revise or remove this function.
#
compute_lfc_stats_helper <- function (x, L, f, mu, stat) {
  S <- compute_poisson_covariance(x,L,f)
  if (stat == "vsmean")
    out <- compute_lfc_vs_mean(f,mu,S)
  else if (stat == "le")
    out <- compute_lfc_le(f,S)
  else
    out <- compute_lfc_pairwise(f,S,stat)
  lfc <- out$lfc
  se  <- out$se
  z   <- lfc/se
  z[lfc == 0] <- 0
  return(list(lfc = lfc/log(2),se = se/log(2),z = z))
}

# Return the log-fold change statistics lfc[i] = log(f[i]/mu), and
# their standard errors se[i], under the Poisson glm with Poisson rates
# L*f, given the MLE estimates (f), the mean or a comparable statistic
# such as the median (mu), and S, the covariance of log(f). Here, the
# log-fold change is defined using the natural logarithm (not the
# base-2 logarithm, which is the convention).
compute_lfc_vs_mean <- function (f, mu, S)
  return(list(lfc = log(f/mu),se = sqrt(diag(S))))

# Return the pairwise log-fold change statistics lfc[j] =
# log(f[j]/f[i]), and their standard errors se[j], under the Poisson
# glm with Poisson rates L*f, given the MLE estimates (f), and S, the
# covariance of log(f). Here, the log-fold change is defined using the
# natural logarithm (not the base-2 logarithm, which is the
# convention).
#
# TO DO: Revise or remove this function.
#
compute_lfc_pairwise <- function (f, S, i) {
  lfc    <- log(f/f[i])
  se     <- sqrt(diag(S) + S[i,i] - 2*S[i,])
  lfc[i] <- 0
  se[i]  <- 0
  return(list(lfc = lfc,se = se))
}

# Return the "least extreme" pairwise log-fold change statistics
# lfc[j] = log(f[j]/f[i]), and their standard errors, se[j], under the
# Poisson glm with Poisson rates L*f. By "least extreme", we mean that
# i != j is chosen so that the log-fold change is closest to
# zero. Here, the log-fold change is defined using the natural
# logarithm (not the base-2 logarithm, which is the convention). The
# inputs are: f, the MLE estimates; and S, the covariance of log(f).
#
# TO DO: Revise or remove this function.
#
compute_lfc_le <- function (f, S) {
  n   <- length(f)
  lfc <- rep(0,n)
  se  <- rep(0,n)
  for (i in 1:n) {
    out    <- compute_lfc_pairwise(f,S,i)
    j      <- order(abs(out$lfc))[2]
    lfc[i] <- -out$lfc[j]
    se[i]  <- out$se[j]
  }
  return(list(lfc = lfc,se = se))
}

# Return the covariance of log(f) in the Poisson glm, x[i] ~
# Poisson(u[i]), in which the Poisson rates are u[i] = L[i,]*f. The
# covariance calculations are based on a Laplace approximation to the
# likelihood. We assume that the input f contains MLEs of the
# coefficients.
compute_poisson_covariance <- function (x, L, f) {
  u <- drop(L %*% f)
  return(chol2inv(chol(tcrossprod(f) * crossprod(sqrt(x)/u*L))))
}

# Simulate draws from the posterior distribution of f via random-walk
# Metropolis on log(f). The posterior distribution is based on a
# uniform prior and the Poisson glm likelihood with. Input ns
# specifies the number of Monte Carlo samples to simulate. Input s
# determines the width (standard deviation) of the random walk, and
# input f is the initial state of the Markov chain. Here, t = log(f)
# is the current state of the Markov chain.
#
# The outputs are (1) "samples", an ns x k matrix of Monte Carlo
# samples of f, where k = length(f), and (2) "ar", the Metropolis
# acceptance rate.
#
# This is mainly used to test simulate_posterior_poisson_rcpp and is
# not actually used in practice because it is too slow for most
# practical purposes.
#
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats dpois
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
      # the additional d in the acceptance probability is needed to
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
    samples[i,] <- exp(t)
  }

  # Output the states of the Markov chain and the acceptance rate.
  return(list(samples = samples,ar = ar/(k*ns)))
}

# Compute the highest posterior density (HPD) interval from a vector
# of random draws from the distribution. See Chen & Shao (1999) for
# background on HPD intervals.
hpd <- function (x, conf.level = 0.95) {
  n <- length(x)
  m <- round(n*(1 - conf.level))
  x <- sort(x)
  y <- x[seq(n-m+1,n)] - x[seq(1,m)]
  i <- which.min(y)
  return(c(x[i],x[n-m+i]))
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

