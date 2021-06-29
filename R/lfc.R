# This is the workhorse function used by de_analysis for computing the
# log-fold change (LFC) statistics. Most of the input arguments are
# explained in the documentation for de_analysis. A couple notes about
# the input arguments: (1) f0 should be a estimate of the paramter f0
# in the the "null" model x ~ Poisson(u), with u = s*f0; (2) F and L
# should specify the parameters of the Poisson glm models; that is, F
# should be returned from fit_poisson_models(X,L,...) in which L =
# s*fit$L, and "fit" is a multinomial topic model fit.
#
# TO DO: Describe the outputs.
#
# TO DO: Allow for calculation of different LFC statistics.
#
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
compute_lfc_stats <- function (X, F, L, f0, stat = "vsnull", ns = 1000,
                               conf.level = 0.9, rw = 0.3, e = 1e-15,
                               nc = 1, verbose = TRUE) {

  # Get the number of columns in the counts matrix (m) and the number
  # of topics (k).
  m <- ncol(X)
  k <- ncol(F)

  # Compute the log-fold change statistics.
  cl  <- makeCluster(nc)
  opb <- pboptions(type = "txt",style = 3,char = "=",txt.width = 70)
  out <- pblapply(1:m,compute_lfc_stats_helper,X,F,L,f0,ns,conf.level,rw,e)
  pboptions(opb)
  stopCluster(cl)

  # Allocate storage for the outputs.
  est  <- matrix(0,m,k)
  mean <- matrix(0,m,k)
  low  <- matrix(0,m,k)
  high <- matrix(0,m,k)
  dimnames(est)  <- dimnames(F)
  dimnames(mean) <- dimnames(F)
  dimnames(low)  <- dimnames(F)
  dimnames(high) <- dimnames(F)

  # Fill in the outputs.
  for (j in 1:m) {
    est[j,]  <- out[[j]]["est",]
    mean[j,] <- out[[j]]["mean",]
    low[j,]  <- out[[j]]["low",]
    high[j,] <- out[[j]]["high",]
  }

  # Compute the z-scores and -log10 p-values.
  z <- est/(2*(mean - low))
  return(list(est   = est/log(2),
              low   = low/log(2),
              high  = high/log(2),
              z     = z,
              lpval = -lpfromz(z)))
}

# This implements the core computation for compute_lfc_stats. See the
# comments accompanying function compute_lfc_stats for details.
#
#' @importFrom stats runif
#' @importFrom stats rnorm
compute_lfc_stats_helper <- function (j, X, F, L, f0, ns, conf.level, rw, e) {
  k <- ncol(F)
  D <- matrix(rnorm(ns*k),ns,k)
  U <- matrix(runif(ns*k),ns,k)
  samples <- simulate_posterior_poisson_rcpp(X[,j],L,F[j,],D,U,rw,e)$samples
  return(compute_lfc_vsnull(samples,F[j,],f0[j],conf.level))
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

