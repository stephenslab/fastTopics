# This is the workhorse function used by de_analysis for computing the
# log-fold change (LFC) statistics using a topic model. It runs a
# simple Markov chain Monte Carlo algorithm to simulate the posterior
# distribution of the LFC statistics, then computes key posterior
# quantities from the simulated Monte Carlo samples.
#
# Most of the input arguments are explained in the documentation
# accompanying de_analysis. Please also bear in mind the following the
# additional points: (1) f0 should be a estimate of the paramter f0 in
# the the "null" model x ~ Poisson(u), with u = s*f0; (2) F and L
# should specify the parameters of the Poisson glm models; that is, F
# should be returned from fit_poisson_models(X,L,...) in which L =
# s*fit$L, and "fit" is a multinomial topic model fit; (3) inputs D
# and U should be ns x k matrices, where k = ncol(F) is the number of
# topics and ns is the number of Monte Carlo samples to simulate, and
# D should contain normally distributed random numbers simulated using
# rnorm(mean = 0,sd = 1, and U should contain uniformly distributed
# random numbers simulated using runif(min = 0,max = 1).
#
# The return value is a list containing five matrices of the same
# dimension as F. Several of these matrices contain posterior
# quantities estimated via MCMC. The matrices are: (1) "est", the
# point estimates of the LFC statistics; (2) "postmean", the posterior
# mean estimates (3) "lower", the estimated lower limits of the HPD
# intervals; (4) "upper", the estimated upper limits of the HPD
# intervals; (5) and "z", the z-scores determined from the LFC
# estimates and the HPD intervals.
#
# Note that all outputted LFC statistics are defined with the base-2
# logarithm.
#
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom Matrix colSums
#' @importFrom progress progress_bar
compute_lfc_stats <- function (X, F, L, f0,
                               D = matrix(rnorm(ns*k),1000,ncol(F)),
                               U = matrix(runif(ns*k),1000,ncol(F)),
                               M = matrix(sample(k,ns*k,replace=TRUE),ns,k)-1,
                               lfc.stat = "le", conf.level = 0.68,
                               rw = 0.3, e = 1e-15, verbose = TRUE) {

  # Get the number of counts matrix columns (m) and the number of
  # topics (k).
  m <- nrow(F)
  k <- ncol(F)

  # Allocate storage for the outputs.
  est      <- matrix(0,m,k)
  postmean <- matrix(0,m,k)
  lower    <- matrix(0,m,k)
  upper    <- matrix(0,m,k)
  ar       <- matrix(0,m,k)
  dimnames(est)      <- dimnames(F)
  dimnames(postmean) <- dimnames(F)
  dimnames(lower)    <- dimnames(F)
  dimnames(upper)    <- dimnames(F)
  dimnames(ar)       <- dimnames(F)

  # Fill in the outputs, row by row. The core computation is performed
  # by compute_lfc_helper.
  ls <- colSums(L)
  if (verbose)
    pb <- progress_bar$new(total = m)
  for (j in 1:m) {
    if (verbose)
      pb$tick()
    out <- compute_lfc_stats_helper(j,X,F,L,D,U,M,ls,f0,lfc.stat,
                                    conf.level,rw,e)
    est[j,]      <- out$dat["est",]
    postmean[j,] <- out$dat["postmean",]
    lower[j,]    <- out$dat["lower",]
    upper[j,]    <- out$dat["upper",]
    ar[j,]       <- out$ar
  }

  # Compute the z-scores, then output the LFC point estimates (est),
  # the posterior means (postmean), lower and upper limits of the HPD
  # intervals, the z-scores (z), and MCMC acceptance rates (ar).
  return(list(ar       = ar,
              est      = est/log(2),
              postmean = postmean/log(2),
              lower    = lower/log(2),
              upper    = upper/log(2),
              z        = compute_zscores(postmean,lower,upper)))
}

# This is the multithreaded variant of compute_lfc_stats. See the
# comments accompanying compute_lfc_stats for details.
#
#' @importFrom parallel splitIndices
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapply
compute_lfc_stats_multicore <- function (X, F, L, f0, D, U, M, lfc.stat,
                                         conf.level, rw, e, nc) {
    
  # Get the number of counts matrix columns (m) and the number of
  # topics (k).
  m <- nrow(F)
  k <- ncol(F)
  
  # Split the data among the nc requesed threads.
  cols <- splitIndices(m,nc)
  dat <- vector("list",nc)
  for (i in 1:nc) {
    j <- cols[[i]]
    dat[[i]] <- list(X = X[,j],F = F[j,],f0 = f0[j])
  }

  # Distribute the calculations using parLapply.
  parlapplyf <- function (dat, L, D, U, M, lfc.stat, conf.level, rw, e)
    compute_lfc_stats(dat$X,dat$F,L,dat$f0,D,U,M,lfc.stat,conf.level,rw,e,
                      verbose = FALSE)
  cl <- makeCluster(nc)
  ans <- parLapply(cl = cl,dat,parlapplyf,L,D,U,M,lfc.stat,conf.level,rw,e)
  stopCluster(cl)

  # Combine the individual compute_lfc_stats outputs, and output the
  # combined result.
  out <- list(ar       = matrix(0,m,k),
              est      = matrix(0,m,k),
              postmean = matrix(0,m,k),
              lower    = matrix(0,m,k),
              upper    = matrix(0,m,k),
              z        = matrix(0,m,k))              
  dimnames(out$ar)       <- dimnames(F)
  dimnames(out$est)      <- dimnames(F)
  dimnames(out$postmean) <- dimnames(F)
  dimnames(out$lower)    <- dimnames(F)
  dimnames(out$upper)    <- dimnames(F)
  dimnames(out$z)        <- dimnames(F)
  for (i in 1:nc) {
    j <- cols[[i]]
    out$ar[j,]       <- ans[[i]]$ar
    out$est[j,]      <- ans[[i]]$est
    out$postmean[j,] <- ans[[i]]$postmean
    out$lower[j,]    <- ans[[i]]$lower
    out$upper[j,]    <- ans[[i]]$upper
    out$z[j,]        <- ans[[i]]$z
  }
  return(out)
}

# This implements the core computation for compute_lfc_stats.
compute_lfc_stats_helper <- function (j, X, F, L, D, U, M, ls, f0,
                                      lfc.stat, conf.level, rw, e) {
  k <- ncol(F)
  if (is.sparse.matrix(X)) {
    dat <- get.nonzeros(X,j)
    out <- simulate_posterior_poisson_sparse_rcpp(dat$x,L[dat$i,],ls,
                                                  F[j,],D,U,M,rw,e)
  } else
    out <- simulate_posterior_poisson_rcpp(X[,j],L,F[j,],D,U,M,rw,e)
  if (lfc.stat == "vsnull")
    dat <- compute_lfc_vsnull(F[j,],f0[j],out$samples,conf.level)
  else if (lfc.stat == "le")
    dat <- compute_lfc_le(F[j,],out$samples,conf.level)
  else
    dat <- compute_lfc_pairwise(F[j,],lfc.stat,out$samples,conf.level)
  return(list(ar = out$ar,dat = dat))
}
  
# Compute posterior statistics for LFC comparing against the null
# parameter, LFC(j) = log(fj/f0). The inputs are: f, the estimates of
# the Poisson model parameters; f0, the estimate of the null model
# parameter; samples, an ns x k matrix of Monte Carlo samples used to
# approximate the posterior distribution of g = log(f), where k is the
# number of topics and ns is the number of samples; and conf.level, a
# number between 0 and 1 giving the desired size of the HPD interval.
#
# The return value is a 4 x k matrix containing LFC statistics: (1)
# point estimate (est), posterior mean (postmean), lower limit of the
# HPD interval (lower) and upper limit (upper).
# 
#' @importFrom Matrix colMeans
compute_lfc_vsnull <- function (f, f0, samples, conf.level) {
  est      <- log(f) - log(f0)
  samples  <- samples - log(f0)
  postmean <- colMeans(samples)
  ans      <- compute_hpd_intervals(samples,conf.level)
  return(rbind(est      = est,
               postmean = postmean,
               lower    = ans$lower,
               upper    = ans$upper))
}

# Compute "pairwise" LFC posterior statistics LFC(k) = log(fk,fj). The
# inputs are: f, the estimates of the Poisson model parameters; j, the
# topic to compare with; samples, an ns x k matrix of Monte Carlo
# samples used to approximate the posterior distribution of g =
# log(f), where k is the number of topics and ns is the number of
# samples; and conf.level, a number between 0 and 1 giving the desired
# size of the HPD interval.
#
# The return value is a 4 x k matrix containing LFC statistics: (1)
# point estimate (est), posterior mean (postmean), lower limit of the
# HPD interval (lower) and upper limit (upper). By definition column j
# of this matrix is all zeros.
#
#' @importFrom Matrix colMeans
compute_lfc_pairwise <- function (f, j, samples, conf.level) {
  k        <- length(f)
  est      <- log(f/f[j])
  samples  <- samples - samples[,j]
  postmean <- colMeans(samples)
  ans      <- compute_hpd_intervals(samples,conf.level,setdiff(1:k,j))
  out      <- rbind(est      = est,
                    postmean = postmean,
                    lower    = ans$lower,
                    upper    = ans$upper)
  out[,j]  <- 0
  return(out)
}

# Compute posterior estimates of the "least extreme" LFC statistics
# LFC(j) = log(fj/fk), in which k is the topic other than j that
# yields the LFC statistic closest to zero. The inputs are: f, the
# estimates of the Poisson model parameters; samples, an ns x k matrix
# of Monte Carlo samples used to approximate the posterior
# distribution of g = log(f), where k is the number of topics and ns
# is the number of samples; and conf.level, a number between 0 and 1
# giving the desired size of the HPD interval.
#
# The return value is a 4 x k matrix containing LFC statistics: (1)
# point estimate (est), posterior mean (postmean), lower limit of the
# HPD interval (lower) and upper limit (upper).
# 
#' @importFrom Matrix colMeans
compute_lfc_le <- function (f, samples, conf.level) {
  est      <- drop(le_diff_rcpp(matrix(log(f),1,length(f))))
  samples  <- le_diff_rcpp(samples)
  postmean <- colMeans(samples)
  ans      <- compute_hpd_intervals(samples,conf.level)
  return(rbind(est      = est,
               postmean = postmean,
               lower    = ans$lower,
               upper    = ans$upper))
}

# Output an HPD interval for each column of the samples matrix.
compute_hpd_intervals <- function (samples, conf.level,
                                   cols = seq(1,ncol(samples))) {
  k     <- ncol(samples)
  lower <- rep(0,k)
  upper <- rep(0,k)
  for (i in cols) {
    out      <- hpd(samples[,i],conf.level)
    lower[i] <- out[1]
    upper[i] <- out[2]
  }
  return(list(lower = lower,upper = upper))
}

# Compute z-scores given posterior mean estimates (postmean), and
# lower and upper limits of the HPD intervals (lower, upper).
compute_zscores <- function (postmean, lower, upper) {
  z <- postmean

  # First fix all the negative z-scores.
  i         <- which(postmean < 0)
  postmean1 <- postmean[i]
  upper1    <- upper[i]
  z1        <- postmean1/(upper1 - postmean1)
  z1[upper1 <= postmean1] <- as.numeric(NA)
  z[i]      <- z1

  # Next, fix all the positive z-scores.
  i         <- which(postmean > 0)
  postmean1 <- postmean[i]
  lower1    <- lower[i]
  z1        <- postmean1/(postmean1 - lower1)
  z1[lower1 >= postmean1] <- as.numeric(NA)
  z[i]      <- z1

  return(z)
}
