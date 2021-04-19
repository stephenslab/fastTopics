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

# TO DO: Explain here what this function does, and how to use it.
# stat: "vsmean", "le", or an index.
compute_lfc_stats <- function (X, F, L, s, mu = colSums(X)/sum(s),
                               stat = "vsmean", version = c("Rcpp","R")) {
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
#  - z-score.
#
# These outputs use the base-2 logarithm.
#
compute_lfc_stats_helper <- function (x, L, f, mu, stat) {
  S <- compute_poisson_covariance(x,L,f)
  if (stat == "vsmean")
    out <- compute_lfc_vs_mean(f,mu,S)
  else if (stat == "le")
    out <- compute_lfc_le(f,S)
  else
    out <- compute_lfc_pairwise(f,S,stat)
  z <- with(out,lfc/se)
  z[out$lfc == 0] <- 0
  return(with(out,list(lfc = lfc/log(2),
                       se  = se/log(2),
                       z   = z)))
}

# Return the log-fold change statistics log(f[i]/mu) and their
# standard errors given the MLE estimates (f), the mean or a
# comparable statistic such as the median (mu), and S, the covariance
# of log(f). Here, the LFC is defined using the natural logarithm (not
# the base-2 logarithm, which is the convention).
compute_lfc_vs_mean <- function (f, mu, S)
  return(list(lfc = log(f/mu),se = sqrt(diag(S))))

# Return the pairwise log-fold change statistics log(f[j]/f[i]) and
# their standard errors given the MLE estimates (f) and S, the
# covariance of log(f). Here, the LFC is defined using the natural
# logarithm (not the base-2 logarithm, which is the convention).
compute_lfc_pairwise <- function (f, S, i) {
  lfc    <- log(f/f[i])
  se     <- sqrt(diag(S) + S[i,i] - 2*S[i,])
  lfc[i] <- 0
  se[i]  <- 0
  return(list(lfc = lfc,se = se))
}

# Return the "least extreme" pairwise log-fold change statistics
# log(f[j]/f[i]) and their standard errors; that is, i != j is chosen
# so that the log-fold change is closest to zero. Here, the LFC is
# defined using the natural logarithm (not the base-2 logarithm, which
# is the convention). The inputs are: f, the MLE estimates; and S, the
# covariance of log(f).
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

# TO DO: Explain here what this function does, and how to use it.
add_pseudocounts <- function (X, L, s, e = 0.01) {
  m <- ncol(X)
  k <- ncol(L)
  X <- rbind(X,matrix(e,k,m))
  L <- rbind(L,diag(k))
  s <- c(s,rep(1,k))
  return(list(X = X,L = L,s = s))
}
