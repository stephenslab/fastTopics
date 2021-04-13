# For each column of the counts matrix, compute maximum-likelihood
# estimates (MLEs) of the parameters in the single-count Poisson
# model. Input argument "method" may be set to "glm" or any valid
# setting for "method" in update_factors_poisson_nmf.
#
#' @importFrom stats glm.control
fit_poisson_models <- function (X, L, s = rep(1,nrow(X)), method, eps = 1e-15,
                                numiter = 100, tol = 1e-8, nc = 1) {

  # Get the number of topics (k) and the number of data columns (m).
  m <- ncol(X)
  k <- ncol(L)
  
  # Enure that none of the topic proportions are exactly zero or
  # exactly one.
  L <- pmax(L,eps)
  L <- pmin(L,1-eps)

  # For each column of X, compute MLEs of the Poisson model parameters.
  if (method == "glm") {
    control <- glm.control(epsilon = tol,maxit = numiter)
    F <- fit_poisson_models_glm(X,L,s,control)
  } else {
    control <- list(numiter = numiter,nc = nc,eps = eps)
    F <- matrix(0.5,m,k)
    F <- update_factors_poisson_nmf(X,F,s*L,1:m,method,control)
  }
  return(F)
}

# Implements fit_poisson_models for method = "glm".
fit_poisson_models_glm <- function (X, L, s, control)  {
  m <- ncol(X)
  k <- ncol(L)
  F <- matrix(0,m,k)
  for (i in 1:m)
    F[i,] <- fit_poisson_glm(X[,i],s*L,control = control)$coef
  return(F)
}

# Fit the Poisson glm in which the Poisson rates are l1*f1 + ... +
# lk*fk. Inputs x and L are the data; x is a vector of length n, and L
# is an n x k matrix. Input argument f is a vector of length k giving
# initial estimates of the coefficients f1,...,fk.
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

# Return the covariance of log(f) in the Poisson glm model. See the
# comments accompanying fit_poisson_glm for an explanation of the
# model and inputs. The covariance calculations are based on a Laplace
# approximation to the likelihood, and assume that input f contains
# MLEs of the coefficients.
compute_poisson_covariance <- function (x, L, f) {
  u <- drop(L %*% f)
  return(chol2inv(chol(tcrossprod(f) * crossprod(sqrt(x)/u*L))))
}

# TO DO: Explain here what this function does, and how to use it.
#
# Outputs:
#  - LFC estimate
#  - se of LFC estimate
#
compute_lfc_stats <- function (x, L, f, mu, stat) {

}

