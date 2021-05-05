# Used to check that x is a vector in which x[i+1] >= x[i] for all i.
expect_nondecreasing <- function (x)
  expect_equal(diff(x) >= 0,rep(TRUE,length(x) - 1))

# Used to check that x is a vector in which x[i+1] <= x[i] for all i.
expect_nonincreasing <- function (x)
  expect_equal(diff(x) <= 0,rep(TRUE,length(x) - 1))

# The Rmosek package on CRAN will not work with REBayes. This function
# is used for some of the tests to check whether the correct Rmosek
# package (the one downloaded from mosek.com) is installed.
skip_if_mixkwdual_doesnt_work <- function() {
  skip_if_not_installed("REBayes")
  skip_if_not_installed("Rmosek")
  skip_if(!is.element("mosek_lptoprob",getNamespaceExports("Rmosek")))
}

# A short convenience function to quickly generate a data set used for
# testing the Poisson NMF optimization methods.
generate_test_data <- function (n, m, k, lmax = 0.5, fmax = 0.5) {
    
  # Generate count data to factorize.
  out <- simulate_count_data(n,m,k,lmax,fmax)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Remove any rows that are entirely zeros.
  rows <- which(rowSums(X > 0) > 0)
  X    <- X[rows,]
  L    <- L[rows,]

  # Choose a row with the least number of nonzero counts, and set all
  # the counts to zero except the largest one.
  i <- which.min(rowSums(X > 0))
  j <- which.max(X[i,])
  X[i,-j] <- 0

  # Remove any columns that are entirely zeros.
  cols <- which(colSums(X > 0) > 0)
  X    <- X[,cols]
  F    <- F[cols,]
  
  # Choose a column with the least number of nonzero counts, and set
  # all the counts to zero except the largest one.
  j <- which.min(colSums(X > 0))
  i <- which.max(X[,j])
  X[-i,j] <- 0
  
  # Remove any rows that are entirely zeros.
  rows <- which(rowSums(X > 0) > 0)
  X    <- X[rows,]
  L    <- L[rows,]

  return(list(X = X,F = F,L = L))
}

# Generate data for testing the fitting of the Poisson mixture model.
generate_poismix_data <- function (n, x) {
  m <- length(x)
  L <- rand(n,m)
  w <- as.double(rpois(n,L %*% x))
  return(list(L = L,w = w))
}

# Given an n x k matrix of topic proportions, where n is the number of
# samples and k is the number of topics, force "hard" topic
# assignments; that is, return a topic proportions matrix with all
# ones and zeros, such that the largest topic proportion in each row
# is replaced with a 1.
force_hard_topic_assignments <- function (L)
  as.matrix(Matrix::sparseMatrix(i = 1:nrow(L),j = apply(L,1,which.max),
                                 x = 1,dims = dim(L)))

# Compute the log-likelihoods for a Poisson non-negative matrix
# factorization. This should give the same result as
# loglik_poisson_nmf, but the computation is done less efficiently
# using dpois.
loglik_poisson_nmf_with_dpois <- function (X, fit) {
  Y <- with(fit,tcrossprod(L,F))
  return(rowSums(stats::dpois(X,Y,log = TRUE)))
}

# Compute the deviance for a Poisson non-negative matrix
# factorization. This should give the same result as
# deviance_poisson_nmf, but the computation is done less efficiently
# using "poisson" from the stats package.
deviance_poisson_nmf_with_poisson <- function (X, fit) {
  Y <- with(fit,tcrossprod(L,F))
  return(rowSums(stats::poisson()$dev.resids(X,Y,1)))
}

# Compute the log-likelihoods for a multinomial topic model. This
# should give the same result as loglik_multinom_topic_model, but the
# computation is done less efficiently using dmultinom.
loglik_multinom_topic_model_with_dmultinom <- function (X, fit) {
  n             <- nrow(fit$L)
  loglik        <- rep(0,n)
  names(loglik) <- rownames(X)
  Y             <- with(fit,tcrossprod(L,F))
  for (i in 1:n)
    loglik[i] <- stats::dmultinom(X[i,],prob = Y[i,],log = TRUE)
  return(loglik)
}

# Iterate the given updates for the factors (F) and loadings (L) matrices.
iterate_updates <- function (X, F, L, numiter, update_factors = NULL,
                             update_loadings = NULL, factors_first = TRUE) {
  loglik <- rep(0,numiter)
  dev    <- rep(0,numiter)
  res    <- rep(0,numiter)
  for (i in 1:numiter) {
    if (factors_first) {
      if (!is.null(update_factors))
        F <- update_factors(X,F,L)
      if (!is.null(update_loadings))
        L <- update_loadings(X,F,L)
    } else {
      if (!is.null(update_loadings))
        L <- update_loadings(X,F,L)
      if (!is.null(update_factors))
        F <- update_factors(X,F,L)
    }
    fit        <- list(F = F,L = L)
    class(fit) <- c("poisson_nmf_fit","list")
    loglik[i]  <- sum(loglik_poisson_nmf(X,fit))
    dev[i]     <- sum(deviance_poisson_nmf(X,fit))
    res[i]     <- with(poisson_nmf_kkt(X,F,L),max(abs(rbind(F,L))))
  }
  return(list(F = F,L = L,loglik = loglik,dev = dev,res = res))
}
