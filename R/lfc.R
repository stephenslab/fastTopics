# TO DO: Explain in detail what this function does, and how (and when)
# to use it.
#
# TO DO:
# - Allow for calculation of different LFC statistics.
#
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
compute_lfc_stats <- function (X, F, L, f0, stat = "vsnull", ns = 1000,
                               conf.level = 0.9, rw = 0.3, e = 1e-15, nc = 1,
                               verbose = TRUE) {

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
  cat("\n")

  # Compute the z-scores and -log10 p-values.
  z <- est/(2*(mean - low))
  return(list(est   = est/log(2),
              low   = low/log(2),
              high  = high/log(2),
              z     = z,
              lpval = -lpfromz(z)))
}

# TO DO: Explain here what this function does, and how to use it.
#' @importFrom stats runif
#' @importFrom stats rnorm
compute_lfc_stats_helper <- function (j, X, F, L, f0, ns, conf.level, rw, e) {
  k <- ncol(F)
  D <- matrix(rnorm(ns*k),ns,k)
  U <- matrix(runif(ns*k),ns,k)
  samples <- simulate_posterior_poisson_rcpp(X[,j],L,F[j,],D,U,rw,e)$samples
  return(compute_lfc_vsf0(samples,F[j,],f0[j],conf.level))
}
  
# TO DO: Explain here what this function does, and how (and when) to
# use it.
#
#' @importFrom Matrix colMeans
compute_lfc_vsf0 <- function (samples, f, f0, conf.level) {
  k    <- length(f)
  mean <- rep(0,k)
  low  <- rep(0,k)
  high <- rep(0,k)
  est  <- log(f) - log(f0)
  mean <- colMeans(samples - log(f0))
  for (i in 1:k) {
    out     <- hpd(samples[,i] - log(f0),conf.level)
    low[i]  <- out[1]
    high[i] <- out[2]
  }
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

