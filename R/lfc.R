# TO DO: Explain in detail what this function does, and how (and when)
# to use it.
#
# TO DO:
# - Implement parlapply version of this.
# - Allow for calculation of different LFC statistics.
# - Add p-value calculation.
#
#' @importFrom stats runif
#' @importFrom stats rnorm
compute_lfc_stats <- function (X, F, L, f0, ns = 1000, conf.level = 0.9,
                               s = 0.3, e = 1e-15, nc = 1) {
  n <- nrow(F)
  k <- ncol(F)

  # Allocate storage for the outputs.
  est  <- matrix(0,n,k)
  low  <- matrix(0,n,k)
  high <- matrix(0,n,k)
  dimnames(est)  <- dimnames(F)
  dimnames(low)  <- dimnames(F)
  dimnames(high) <- dimnames(F)
  
  # Repeat for each column of X.
  for (i in 1:n) {
    cat(i," ",sep="")
    D        <- matrix(rnorm(ns*k),ns,k)
    U        <- matrix(runif(ns*k),ns,k)
    samples  <- simulate_posterior_poisson_rcpp(X[,i],L,F[i,],D,U,s,e)$samples
    out      <- compute_lfc_vsf0(samples,F[i,],f0[i],conf.level)
    est[i,]  <- out["est",]
    low[i,]  <- out["low",]
    high[i,] <- out["high",]
  }
  cat("\n")

  # Compute the z-scores and -log10 p-values.
  z <- est/(high - low)
  return(list(est = est/log(2),low = low/log(2),high = high/log(2),
              z = z,lpval = -lpfromz(z)))
}

# TO DO: Explain here what this function does, and how (and when) to
# use it.
compute_lfc_vsf0 <- function (samples, f, f0, conf.level) {
  k    <- length(f)
  low  <- rep(0,k)
  high <- rep(0,k)
  est  <- log(f) - log(f0)
  for (i in 1:k) {
    out     <- hpd(samples[,i] - log(f0),conf.level)
    low[i]  <- out[1]
    high[i] <- out[2]
  }
  return(rbind(est  = est,
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

