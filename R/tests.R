# Used to check that x is a vector in which x[i+1] >= x[i] for all i.
expect_nondecreasing <- function (x)
  expect_equal(diff(x) >= 0,rep(TRUE,length(x) - 1))

# Used to check that x is a vector in which x[i+1] <= x[i] for all i.
expect_nonincreasing <- function (x)
  expect_equal(diff(x) <= 0,rep(TRUE,length(x) - 1))

# A short convenience function to quickly generate a data set used for
# testing.
generate_test_data <- function (n, m, k, lmax = 1, fmax = 1) {
    
  # Generate count data to factorize.
  out <- simulate_count_data(n,m,k,lmax,fmax)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Remove any rows and columns that are entirely zeros.
  rows <- which(rowSums(X > 0) > 0)
  cols <- which(colSums(X > 0) > 0)
  X    <- X[rows,cols]
  L    <- L[rows,]
  F    <- F[cols,]
  return(list(X = X,F = F,L = L))
}

# Generate data for testing the fitting of the Poisson mixture model.
generate_poismix_data <- function (n, x) {
  m <- length(x)
  L <- matrix(runif(n*m),n,m)
  w <- rpois(n,L %*% x)
  return(list(L = L,w = w))
}

# Iterate the given updates for the factors (F) and loadings (L) matrices.
iterate_updates <- function (X, F, L, numiter, update_factors,
                             update_loadings) {
  loglik <- rep(0,numiter)
  dev    <- rep(0,numiter)
  for (i in 1:numiter) {
    F         <- update_factors(X,F,L)
    L         <- update_loadings(X,F,L)
    loglik[i] <- sum(loglik_poisson_nmf(X,list(F = F,L = L)))
    dev[i]    <- sum(deviance_poisson_nmf(X,list(F = F,L = L)))
  }
  return(list(F = F,L = L,loglik = loglik,dev = dev))
}
