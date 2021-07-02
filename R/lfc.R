# This is the workhorse function used by de_analysis for computing the
# log-fold change (LFC) statistics using a topic model. Most of the
# input arguments are explained in the documentation accompanying
# de_analysis. Please also bear in mind the following the additional
# points: (1) f0 should be a estimate of the paramter f0 in the the
# "null" model x ~ Poisson(u), with u = s*f0; (2) F and L should
# specify the parameters of the Poisson glm models; that is, F should
# be returned from fit_poisson_models(X,L,...) in which L = s*fit$L,
# and "fit" is a multinomial topic model fit; (3) inputs D and U
# should be ns x k matrices, where k = ncol(F) is the number of topics
# and ns is the number of Monte Carlo samples to simulate, and D
# should contain normally distributed random numbers simulated using
# rnorm(mean = 0,sd = 1, and U should contain uniformly distributed
# random numbers simulated using runif(min = 0,max = 1).
#
# The return value is a list containing five matrices of the same
# dimension as F. The matrices are: (1) "est", the estimated LFC
# statistics; (2) "low", the estimated lower limits of the HPD
# intervals; (3) "high", the estimated upper limits of the HPD
# intervals; (4) "z", the z-scores determined from the LFC estimates
# and the HPD intervals; and (5) "lpval", the -log10 two-tailed
# p-values computed from the z-scores. Note that all outputted LFC
# statistics are defined with the base-2 logarithm.
#
# TO DO: Allow for calculation of different LFC statistics.
#
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom Matrix colSums
compute_lfc_stats <- function (X, F, L, f0,
                               D = matrix(rnorm(ns*k),1000,ncol(F)),
                               U = matrix(runif(ns*k),1000,ncol(F)),
                               lfc.stat = "vsnull", conf.level = 0.9,
                               rw = 0.3, e = 1e-15) {

  # Get the number of counts matrix columns (m) and the number of
  # topics (k).
  m <- nrow(F)
  k <- ncol(F)

  # Allocate storage for the outputs.
  est  <- matrix(0,m,k)
  mean <- matrix(0,m,k)
  low  <- matrix(0,m,k)
  high <- matrix(0,m,k)
  ar   <- rep(0,m)
  dimnames(est)  <- dimnames(F)
  dimnames(mean) <- dimnames(F)
  dimnames(low)  <- dimnames(F)
  dimnames(high) <- dimnames(F)
  names(ar)      <- rownames(F)

  # Fill in the outputs, row by row. The core computation is performed
  # by compute_lfc_helper.
  ls <- colSums(L)
  for (j in 1:m) {
    out <- compute_lfc_stats_helper(j,X,F,L,D,U,ls,f0,lfc.stat,conf.level,rw,e)
    est[j,]  <- out$dat["est",]
    mean[j,] <- out$dat["mean",]
    low[j,]  <- out$dat["low",]
    high[j,] <- out$dat["high",]
    ar[j]    <- out$ar
  }

  # Compute the z-scores and -log10 p-values.
  z <- est/(2*(mean - low))
  return(list(est   = est/log(2),
              low   = low/log(2),
              high  = high/log(2),
              z     = z,
              lpval = -lpfromz(z),
              ar    = ar))
}

# This is the multithreaded variant of compute_lfc_stats. See the
# comments accompanying compute_lfc_stats for details.
#
#' @importFrom parallel splitIndices
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
compute_lfc_stats_multicore <- function (X, F, L, f0, D, U, lfc.stat,
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

  # Distribution the calculations using pblapply.
  pblapplyf <- function (dat, L, D, U, lfc.stat, conf.level, rw, e)
    compute_lfc_stats(dat$X,dat$F,L,dat$f0,D,U,lfc.stat,conf.level,rw,e)
  pbo <- pboptions(type = "txt",style = 3,char = "=",txt.width = 70)
  cl <- makeCluster(nc)
  ans <- pblapply(cl = cl,dat,pblapplyf,L,D,U,lfc.stat,conf.level,rw,e)
  pboptions(pbo)
  stopCluster(cl)

  # Combine the individual compute_lfc_stats outputs, and output the
  # combined result.
  out <- list(est   = matrix(0,m,k),
              low   = matrix(0,m,k),
              high  = matrix(0,m,k),
              z     = matrix(0,m,k),
              lpval = matrix(0,m,k),
              ar    = rep(0,m))
  dimnames(out$est)   <- dimnames(F)
  dimnames(out$low)   <- dimnames(F)
  dimnames(out$high)  <- dimnames(F)
  dimnames(out$z)     <- dimnames(F)
  dimnames(out$lpval) <- dimnames(F)
  names(out$ar)       <- rownames(F)
  for (i in 1:nc) {
    j <- cols[[i]]
    out$est[j,]   <- ans[[i]]$est
    out$low[j,]   <- ans[[i]]$low
    out$high[j,]  <- ans[[i]]$high
    out$z[j,]     <- ans[[i]]$z
    out$lpval[j,] <- ans[[i]]$lpval
    out$ar[j]     <- ans[[i]]$ar
  }
  return(out)
}

# This implements the core computation for compute_lfc_stats.
compute_lfc_stats_helper <- function (j, X, F, L, D, U, ls, f0, lfc.stat,
                                      conf.level, rw, e) {
  k <- ncol(F)
  if (is.sparse.matrix(X)) {
    dat <- get.nonzeros(X,j)
    out <- simulate_posterior_poisson_sparse_rcpp(dat$x,L[dat$i,],ls,
                                                  F[j,],D,U,rw,e)
  } else
    out <- simulate_posterior_poisson_rcpp(X[,j],L,F[j,],D,U,rw,e)
  dat <- compute_lfc_vsnull(out$samples,F[j,],f0[j],conf.level)
  return(list(dat = dat,ar = out$ar))
}
  
# TO DO: Explain here what this function does, and how (and when) to
# use it.
#
#' @importFrom Matrix colMeans
compute_lfc_vsnull <- function (samples, f, f0, conf.level) {

  # Set up storage for some of the the outputs.
  k    <- length(f)
  mean <- rep(0,k)
  low  <- rep(0,k)
  high <- rep(0,k)
  
  est  <- log(f) - log(f0)
  mean <- colMeans(samples - log(f0))

  # Repeat for each topic.
  for (i in 1:k) {
    out     <- hpd(samples[,i] - log(f0),conf.level)
    #  if (any(is.na(out)))
    #   browser()
    low[i]  <- out[1]
    high[i] <- out[2]
  }

  # TO DO: Summarize the return values.
  return(rbind(est  = est,
               mean = mean,
               low  = low,
               high = high))
}

# TO DO: Explain here what this function does, and how (and when) to
# use it.
compute_lfc_pairwise <- function (samples, t, conf.level) {
  # TO DO.
}

# TO DO: Explain here what this function does, and how (and when) to
# use it.
compute_lfc_le <- function (samples, t, conf.level) {
  # TO DO.
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

