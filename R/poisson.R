# Compute the Poisson rates given the topic proportions (q) and
# parameters (f0, f1) of the Poisson model.
get_poisson_rates <- function (q, f0, f1)
  (1-q)*f0 + q*f1

# This should give the same, or nearly the same, result as
# sum(dpois(x,y,log = TRUE)), except that terms that do not depend on
# the Poisson rates u are disregarded.
loglik_poisson <- function (x, y, e)
  return(sum(x*log(y + e) - y))

# For each column of the counts matrix, and for each topic, compute
# maximum-likelihood estimates (MLEs) of the parameters in the
# single-count (univariate) Poisson model. For method = "glm",
# accompanying statistics (standard errors, z-scores and -log10
# p-values) are also outputted.
fit_univar_poisson_models <-
  function (X, L, s = rep(1,nrow(X)), method = c("em-rcpp","em","optim","glm"),
            e = 1e-15, numiter = 100, tol = 1e-8, verbose = TRUE) {
  method <- match.arg(method)

  # Get the number of topics (k) and the number of columns of X (m).
  m <- ncol(X)
  k <- ncol(L)
  
  # Make sure none of the topic proportions are exactly zero or
  # exactly one.
  L <- pmax(L,e)
  L <- pmin(L,1-e)

  # Compute MLEs of the Poisson model parameters for each topic and
  # for each column of X.
  if (method == "glm") {
    control <- list(epsilon = tol,maxit = numiter)
    out <- fit_univar_poisson_models_glm(X,L,s,e,control,verbose)
  } else if (method == "optim") {
    control <- list(maxit = numiter,factr = tol*.Machine$double.eps)
    out <- fit_univar_poisson_models_optim(X,L,s,e,control,verbose)
  } else if (method == "em")
    out <- fit_univar_poisson_models_em(X,L,s,e,numiter,verbose)
  else if (method == "em-rcpp") {
    if (is.sparse.matrix(X))
      out <- fit_univar_poisson_models_em_sparse_rcpp(X,L,s,e,numiter,tol,
                                                      verbose)
    else
      out <- fit_univar_poisson_models_em_rcpp(X,L,s,e,numiter,tol,verbose)
  }
  return(out)
}

# Produces the same result as fit_univar_poisson_models for the
# special case of "hard" topic assignments; that is, when the topic
# proportions matrix, L, is entirely zeros and ones.
#
#' @importFrom Matrix colSums
fit_univar_poisson_models_hard <- function (X, L, s = rep(1,nrow(X)),
                                            e = 1e-15) {
  m      <- ncol(X)
  k      <- ncol(L)
  L      <- round(L)
  F0     <- matrix(0,m,k)
  F1     <- matrix(0,m,k)
  loglik <- matrix(0,m,k)
  for (j in 1:k) {
    i0         <- which(L[,j] == 0)
    i1         <- which(L[,j] == 1)
    F0[,j]     <- colSums(X[i0,])/sum(s[i0])
    F1[,j]     <- colSums(X[i1,])/sum(s[i1])
    loglik[,j] <- -(cost(t(X[i0,]),F0[,j],s[i0],e,"poisson") +
                    cost(t(X[i1,]),F1[,j],s[i1],e,"poisson"))
  }
  return(list(F0 = F0,F1 = F1,loglik = loglik))
}

# Implements fit_univar_poisson_models for method = "glm".
#
#' @importFrom progress progress_bar
fit_univar_poisson_models_glm <- function (X, L, s, e, control, verbose)  {
  m    <- ncol(X)
  k    <- ncol(L)
  F0   <- matrix(0,m,k)
  F1   <- matrix(0,m,k)
  beta <- matrix(0,m,k)
  Z    <- matrix(0,m,k)
  se   <- matrix(0,m,k)
  pval <- matrix(0,m,k)
  if (verbose)
    pb <- progress_bar$new(total = m)
  for (i in 1:m) {
    if (verbose)
      pb$tick()
    for (j in 1:k) {
      out       <- fit_poisson_glm(X[,i],s,L[,j],control = control)
      F0[i,j]   <- out["f0"]
      F1[i,j]   <- out["f1"]
      beta[i,j] <- out["beta"]
      se[i,j]   <- out["se"]
      Z[i,j]    <- out["z"]
      pval[i,j] <- out["pval"]
    }
  }
  return(list(F0 = F0,F1 = F1,beta = beta,se = se,Z = Z,pval = pval))
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

# Fit the Poisson generalized linear model (glm) in which the Poisson
# rates are s*(b0 + q*b). This is a reparameterization of the model x
# x ~ Poisson(s*u), with u = f0*(1-q) + f1*q; the glm "identity"
# reparameterization is recovered by b0 = f0, b = f1 - f0. Input
# arguments f0 and f1 are initial parameter estimates. The outputs are
# the estimates of f0 and f1, the log2-fold change statistic (beta),
# the standard error (se) of b = f1 - f0, the z-score (z), and the
# -log10 p-value (pval).
#
#' @importFrom stats glm
#' @importFrom stats glm.control
#' @importFrom stats summary.glm
#'
fit_poisson_glm <-
  function (x, s, q, f0 = 0.5, f1 = 0.5,
            control = glm.control(epsilon = 1e-10, maxit = 100)) {

  # Fit the generalized linear model.
  dat <- data.frame(x = x,b0 = s,b = s*q)
  fit <- suppressWarnings(glm(x ~ b0 + b - 1,
                              family = poisson(link = "identity"),
                              data = dat,start = c(f0,f1 - f0),
                              control = control))

  # Output the parameter estimates, log-fold changes and test statistics.
  ans <- summary.glm(fit)$coefficients
  b0  <- ans["b0","Estimate"]
  b   <- ans["b","Estimate"]
  f0  <- b0
  f1  <- b + b0
  return(c(f0   = f0,
           f1   = f1,
           beta = log2(f1/f0),
           se   = ans["b","Std. Error"],
           z    = ans["b","z value"],
           pval = -log10(ans["b","Pr(>|z|)"])))
}
    
# Compute MLEs of the parameters in the single-count Poisson model: x
# ~ Poisson(s*u), with u = f0*(1-q) + f1*q. Parameters f0, f1 are
# estimated, and vectors s, q are provided. The likelihood is
# maximized by optim with method = "L-BFGS-B". Input arguments f0 and
# f1 are initial estimates of the parameters. Input "control" is a
# list of control parameters passed to optim. Input argument "e" is a
# small positive scalar added to the likelihood and gradient to avoid
# logarithms of zero and division by zero, and it is the lower bound
# on f0 and f1.
#
#' @importFrom stats optim
#' 
fit_poisson_optim <- function (x, s, q, f0 = 0.5, f1 = 0.5, e = 1e-15,
                               control = list(factr = 1e5, maxit = 100)) {

  # This is used to compute the negative log-likelihood, in which
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
    return(c(sum(y*(1 - q)),sum(y*q)))
  }
  
  # Fit the Poisson rates f0, f1 using the L-BFGS-B quasi-Newton method
  # implemented in optim.
  out <- optim(c(f0,f1),f,g,method = "L-BFGS-B",lower = c(e,e),
               control = control)

  # Output the MLEs of f0 and f1, and the other "optim" outputs.
  names(out$par) <- c("f0","f1")
  return(out)
}

# Run EM for maximum-likelihood estimation of the parameters in the
# single-count Poisson model: x ~ Poisson(s*u), with u given by u =
# f0*(1-q) + f1*q. Parameters f0, f1 are estimated, and vectors s, q
# are provided. Input arguments f0 and f1 are initial estimates of the
# parameters. Input argument "e" is a small positive scalar added to
# the denominator in the E-step to avoid division by zero.
fit_poisson_em <- function (x, s, q, f0 = 0.5, f1 = 0.5, e = 1e-15,
                            numiter = 100) {

  # Store a couple pre-calculations to simplify the calculations
  # below.
  a <- sum(s*(1-q)) + e
  b <- sum(s*q) + e
  
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
  return(list(f = c(f0 = f0,f1 = f1),loglik = loglik))
}

# Compute z-scores and other statistics given MLEs of the Poisson
# model parameters for each topic and for each count (column of the
# counts matrix). The inputs are the n x m counts matrix (X), the m x
# k matrix of topic proportions (L), and two m x k matrices containing
# estimmates of the Poisson model parameters (F0, F1). There are also
# two optional arguments: s, a vector of length n containing the "size
# factors" (by default, they are all 1); and e, a small, positive
# constant used to prevent any of the calculations giving Inf or
# NaN. The outputs are all m x k matrices, where m is the number of
# columns in X, and k is the number of topics: (base-2) log-fold
# change statistics (beta); standard errors (se) and z-scores (z) for
# the unknowns b = f1 - f0; and -log10 two-tailed p-values (pval)
# computed from the z-scores.
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
  # compute the log-fold change (beta), standard error (se), z-score
  # and -log10 p-value.
  for (i in 1:m)
    for (j in 1:k) {
      out       <- compute_poisson_zscore(X[,i],L[,j],s,F0[i,j],F1[i,j],e)
      beta[i,j] <- out["beta"]
      se[i,j]   <- out["se"]
      Z[i,j]    <- out["z"]
      pval[i,j] <- out["pval"]
    }

  # Return the model parameters (F0, F1), the (base-2) log-fold change
  # statistics (beta), the standard errors (se), the z-scores (z), and
  # the -log10 two-sided p-values (pval).
  return(list(F0 = F0,F1 = F1,beta = beta,se = se,Z = Z,pval = pval))
}

# Given the counts (x), topic proportions (q), "size factors" (s), and
# estimates of the Poisson model parameters, f0 and f1, return: (1)
# the (base-2) log-fold change statistic beta = log2(f1/f0), the
# standard error (se) and z-score (z) for b = f1 - f0, and the -log10
# two-tailed p-value (pval) computed from the z-score.
compute_poisson_zscore <- function (x, q, s, f0, f1, e) {

  # Ensure that the Poisson model parameters are positive.
  f0 <- max(f0,e)
  f1 <- max(f1,e)
    
  # Compute the standard error.
  u  <- get_poisson_rates(q,f0,f1)
  a  <- sum(x/u^2)
  b  <- sum(x*q/u^2)
  c  <- sum(x*(q/u)^2)
  H  <- rbind(c(a,b),
              c(b,c))
  se <- tryCatch(sqrt(diag(solve(H)))[2],error = function (e) NA)
  
  # Return the (base-2) log-fold change statistic (beta), the standard
  # error (se) and z-score (z) for b, and the -log10 two-tailed
  # p-value (pval) computed from z.
  return(unlist(compute_poisson_zscore_helper(f0,f1,se)))
}

# The same as compute_univar_poisson_zscores, but the computation is
# implemented much more efficiently using matrix operations.
#
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
compute_univar_poisson_zscores_fast <-
  function (X, L, F0, F1, s = rep(1,nrow(X)), e = 1e-15) {

  # Get the number of rows (n) and columns (m) in the counts matrix,
  # and the number of topics (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(F0)
    
  # Ensure that the Poisson model parameters are positive.
  F0 <- pmax(F0,e)
  F1 <- pmax(F1,e)

  # Compute the standard errors.
  if (is.sparse.matrix(X)) {
    a <- compute_poisson_glm_stat_sparse(X,matrix(1,n,k),L,F0,F1)
    b <- compute_poisson_glm_stat_sparse(X,L,L,F0,F1)
    c <- compute_poisson_glm_stat_sparse(X,L^2,L,F0,F1)
  } else {
    a <- compute_poisson_glm_stat(X,matrix(1,n,k),L,F0,F1)
    b <- compute_poisson_glm_stat(X,L,L,F0,F1)
    c <- compute_poisson_glm_stat(X,L^2,L,F0,F1)
  }
  se <- compute_laplace_se(a,b,c)

  # Return the model parameters (F0, F1), the (base-2) log-fold change
  # statistics (beta), the standard errors (se), the z-scores (Z), and
  # the -log10 two-tailed p-values (pval).
  out <- compute_poisson_zscore_helper(F0,F1,se)
  names(out) <- c("F0","F1","beta","se","Z","pval")
  return(out)
}

# This is used by compute_univar_poisson_zscores_fast to compute the
# precisions for the Poisson glm model parameters when X is a dense
# matrix. Here we assume the matrices F0 and F1 contain only positive
# values. Input argument should be a matrix of the same dimension as
# L; the other input arguments are described in the comments
# accompanying the function compute_univar_poisson_zscores. The output
# is a matrix H in which H[i,j] = sum(X[,i]*Y[,j]/u^2), where u is the
# vector of Poisson rates u = (1 - L[,j])*F0[,j] + L[,j]*F1[,j].
compute_poisson_glm_stat <- function (X, Y, L, F0, F1) {
  m <- nrow(F0)
  k <- ncol(F0)
  H <- matrix(0,m,k)
  for (i in 1:k) {
    u     <- outer(1 - L[,i],F0[,i]) + outer(L[,i],F1[,i])
    H[,i] <- colSums(X*(Y[,i]/u^2))
  }
  return(H)
}

# This is used by compute_univar_poisson_zscores_fast to compute the
# precisions for the Poissonn glm model parameters when X is a sparse
# matrix. Here we assume the matrices F0 and F1 contain only positive
# values. See the comments accompanying the function
# compute_poisson_glm_stat for additional details.
#
#' @importFrom Matrix colSums
compute_poisson_glm_stat_sparse <- function (X, Y, L, F0, F1) {
  m <- nrow(F0)
  K <- ncol(F0)
  H <- matrix(0,m,K)
  for (k in 1:K) {
    f0    <- F0[,k]
    f1    <- F1[,k]
    out   <- summary(X)
    i     <- out$i
    j     <- out$j
    x     <- out$x
    u     <- (1 - L[i,k])*F0[j,k] + L[i,k]*F1[j,k]
    H[,k] <- colSums(sparseMatrix(i,j,x = x*(Y[i,k]/u^2),dims = dim(X)))
  }
  return(H)
}
  
# This is used by compute_univar_poisson_zscores_fast to compute the
# Laplace approximation to the standard error given the precisions
# (the entries of the 2 x 2 Hessian matrix).
compute_laplace_se <- function (a, b, c) {
  se <- suppressWarnings(sqrt(a)/(sqrt(a*c - b^2)))
  se[a*c <= b^2] <- NA
  return(se)
}

# This is used by the functions compute_poisson_zscore and
# compute_univar_poisson_zscores_fast to prepare the final outputs.
compute_poisson_zscore_helper <- function (f0, f1, se) {
  b <- f1 - f0
  z <- b/se
  z[is.na(se)] <- 0
  return(list(f0   = f0,
              f1   = f1,
              beta = log2(f1/f0),
              se   = se,
              z    = z,
              pval = -lpfromz(z)))
}
