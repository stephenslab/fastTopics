# Compute the Poisson rates given the topic proportions (q) and
# parameters (f0, f1) of the Poisson model.
get_poisson_rates <- function (q, f0, f1)
  (1-q)*f0 + q*f1

# This should give the same, or nearly the same, result as
# sum(dpois(x,y,log = TRUE)), except that terms that do not depend on
# the Poisson rates u are disregarded.
loglik_poisson <- function (x, y, e = 1e-15)
  return(sum(x*log(y+e) - y))

# For each column of the counts matrix, and for each topic, compute
# maximum-likelihood estimates (MLEs) of the parameters in the
# single-count (univariate) Poisson model.
fit_univar_poisson_models <-
  function (X, L, s = rep(1,nrow(X)),
            method = c("em-rcpp","em","optim"),
            e = 1e-15, numiter.em = 40,
            control.optim = list(factr = 1e-14,maxit = 100),
            verbose = TRUE) {
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
    if (is.sparse.matrix(X))
      out <- fit_univar_poisson_models_em_sparse_rcpp(X,L,s,e,numiter.em,
                                                      verbose)
    else
      out <- fit_univar_poisson_models_em_rcpp(X,L,s,e,numiter.em,verbose)
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
    u <- get_poisson_rates(q,par[1],par[2])
    return(-loglik_poisson(x,s*u,e))
  }
  
  # Returns the gradient of the negative log-likelihood, in which
  # par = c(f0,f1).
  g <- function (par) {
    u <- get_poisson_rates(q,par[1],par[2])
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
    u            <- get_poisson_rates(q,f0,f1)
    loglik[iter] <- loglik_poisson(x,s*u,e)
  }

  # Output the estimates of f0 and f1, and the log-likelihood at each EM
  # iteration.
  f        <- c(f0,f1)
  names(f) <- c("f0","f1")
  return(list(f = f,loglik = loglik))
}

# Compute z-scores and other statistics given maximum-likelihood
# estimates (MLEs) of the Poisson model parameters, f0 and f1, for
# each topic, and each count (column of the counts matrix). The inputs
# are the n x m counts matrix (X), the m x k matrix of topic
# proportions (L), and two m x k matrices containing MLEs of the
# Poisson model parameters (F0, F1). There are also two optional
# arguments: s, a vector of length n containing the "size factors" (by
# default, they are all 1); and s, a small, positive constant used to
# prevent any of the calculations giving Inf or NaN. The outputs are
# all m x k matrices, where m is the number of columns in X, and k is
# the number of topics: (base-2) log-fold change statistics (beta);
# their standard errors (se); z-scores (z); and two-tailed p-values
# computed from the z-score (pval).
compute_univar_poisson_zscores <- function (X, L, F0, F1, s = rep(1,nrow(X)),
                                            e = 1e-15) {

  # Get the number of columns in the counts matrix (m) and the number
  # of topics (k).
  m <- nrow(F0)
  k <- ncol(F0)

  # Initialize storage for the outputs.
  beta <- matrix(0,m,k)
  se   <- matrix(0,m,k)
  Z    <- matrix(0,m,k)
  pval <- matrix(0,m,k)

  # For each column of the counts matrix and for for each topic,
  # compute the z-score for the log-fold change statistic.
  for (i in 1:m)
    for (j in 1:k) {
      out       <- compute_poisson_zscore(X[,i],L[,j],s,F0[i,j],F1[i,j])
      beta[i,j] <- out["beta"]
      se[i,j]   <- out["se"]
      Z[i,j]    <- out["Z"]
      pval[i,j] <- out["pval"]
    }

  # Return the (base-2) log-fold change statistics (beta), their
  # standard errors (se), the z-scores (z), and the two-sided
  # p-values (pval).
  return(list(beta = beta,se = se,Z = Z,pval = pval))
}

# This does the same thing as compute_univar_poisson_zscores, but the
# computation is implemented much more efficiently.
#
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
compute_univar_poisson_zscores_fast <- function (X, L, F0, F1,
                                                 s = rep(1,nrow(X)),
                                                 e = 1e-15) {

  # Get the number of columns in the counts matrix (m) and the number
  # of topics (k).
  m <- nrow(F0)
  k <- ncol(F0)
    
  # Ensure that the Poisson model parameters (f0,f1) are positive.
  F0 <- pmax(F0,e)
  F1 <- pmax(F1,e)

  # Compute the standard errors.
  a <- matrix(colSums(X),m,k)
  b <- matrix(colSums(s*L),m,k,byrow = TRUE)
  if (is.sparse.matrix(X))
    c <- compute_poisson_beta_stat_sparse(X,L,F0,F1)
  else
    c <- compute_poisson_beta_stat(X,L,F0,F1)
  se <- compute_poisson_beta_se(F1,a,b,c)

  # Return the (base-2) log-fold change statistics (beta), the
  # standard errors (se), the z-scores (Z), and the two-tailed
  # p-values (pval).
  return(compute_poisson_zscore_helper(F0,F1,se))
}

# Given the counts (x), topic proportions (q), "size factors" (s), and
# estimates of the Poisson model parameters, f0 and f1, return: (1)
# the (base-2) log-fold change statistic beta = log2(f1/f0), the
# standard error of the log-fold change (se), the z-score (z), and the
# two-tailed p-value computed from the z-score (pval).
compute_poisson_zscore <- function (x, q, s, f0, f1, e = 1e-15) {

  # Ensure that the Poisson model parameters (f0,f1) are positive.
  f0 <- max(f0,e)
  f1 <- max(f1,e)
    
  # Compute the standard error. The last line should give the same
  # result as this line of code, but the computation is done in a more
  # numerically stable way:
  #
  #   se <- sqrt(diag(solve(rbind(c(a/f0^2,f1/f0*b),
  #                               c(f1/f0*b,f1^2*c)))))[2]
  #  
  u  <- get_poisson_rates(q,f0,f1)
  a  <- sum(x)
  b  <- sum(s*q)
  c  <- sum(x*(q/u)^2)
  se <- compute_poisson_beta_se(f1,a,b,c)

  # Return the (base-2) log-fold change statistic (beta), the standard
  # error of the log-fold change (se), the z-score (z), and the
  # two-tailed p-value (pval).
  return(unlist(compute_poisson_zscore_helper(f0,f1,se)))
}

# This is used by compute_univar_poisson_zscores_fast to compute the
# precisions for the log-fold change statistics (beta) when X is a
# dense matrix.
compute_poisson_beta_stat <- function (X, L, F0, F1) {
  m <- nrow(F0)
  k <- ncol(F0)
  c <- matrix(0,m,k)
  for (i in 1:k) {
    u     <- outer(1 - L[,i],F0[,i]) + outer(L[,i],F1[,i])
    c[,i] <- colSums(X*(L[,i]/u)^2)
  }
  return(c)
}

# This is used by compute_univar_poisson_zscores_fast to compute the
# precisions for the log-fold change statistics (beta) when X is a
# sparse matrix.
#
#' @importFrom Matrix colSums
compute_poisson_beta_stat_sparse <- function (X, L, F0, F1) {
  m <- nrow(F0)
  K <- ncol(F0)
  c <- matrix(0,m,K)
  for (k in 1:K) {
    f0    <- F0[,k]
    f1    <- F1[,k]
    out   <- summary(X)
    i     <- out$i
    j     <- out$j
    x     <- out$x
    u     <- (1 - L[i,k])*F0[j,k] + L[i,k]*F1[j,k]
    c[,k] <- colSums(sparseMatrix(i = i,j = j,x = x*(L[i,k]/u)^2,
                                  dims = dim(X)))
  }
  return(c)
}
  
# This is used by compute_poisson_zscore and
# compute_univar_poisson_zscores_fast to compute the standard error of
# the log-fold change statistics (beta) given the summary statstics
# (a, b, c). The inputs may be scalars, vectors or matrices.
compute_poisson_beta_se <- function (f1, a, b, c) {
  se <- suppressWarnings(sqrt(a)/(f1*sqrt(a*c - b^2)))
  se[a*c < b^2] <- NA
  return(se)
}

# This is used by compute_poisson_zscore and
# compute_univar_poisson_zscores_fast to prepare the final outputs for
# those two functions.
compute_poisson_zscore_helper <- function (f0, f1, se) {
  b <- log(f1/f0)
  Z <- b/se
  Z[is.na(se)] <- 0
  return(list(beta = b/log(2),
              se   = se/log(2),
              Z    = Z,
              pval = pfromz(Z)))
}
     

