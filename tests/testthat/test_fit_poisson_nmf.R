context("fit_poisson_nmf")

library(Matrix)
library(RcppParallel)

test_that(paste("betanmf and pnmfem updates produce same result, and",
                "increase the likelihood"),{
                    
  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  n   <- 80
  m   <- 100
  out <- simulate_count_data(n,m,3,lmax = 0.3,fmax = 0.3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Remove rows and columns containing all zeros.
  rows <- which(rowSums(X > 0) > 0)
  cols <- which(colSums(X > 0) > 0)
  X    <- X[rows,cols]
  L    <- L[rows,]
  F    <- F[cols,]

  # There should be at least one row and one column in the counts
  # matrix containing exactly one nonzero count.
  expect_gt(sum(rowSums(X > 0) == 1),0)
  expect_gt(sum(colSums(X > 0) == 1),0)

  # Run 20 multiplicative updates.
  n       <- 20
  F1      <- F
  L1      <- L
  loglik1 <- rep(0,n)
  dev1    <- rep(0,n)
  for (i in 1:n) {
    F1         <- t(betanmf_update_factors(X,L1,t(F1)))
    L1         <- betanmf_update_loadings(X,L1,t(F1))
    loglik1[i] <- sum(loglik_poisson_nmf(X,list(F = F1,L = L1)))
    dev1[i]    <- sum(deviance_poisson_nmf(X,list(F = F1,L = L1)))
  }

  # Run 20 EM updates.
  F2      <- F
  L2      <- L
  loglik2 <- rep(0,n)
  dev2    <- rep(0,n)
  for (i in 1:n) {
    F2         <- pnmfem_update_factors(X,F2,L2)
    L2         <- pnmfem_update_loadings(X,F2,L2)
    loglik2[i] <- sum(loglik_poisson_nmf(X,list(F = F2,L = L2)))
    dev2[i]    <- sum(deviance_poisson_nmf(X,list(F = F2,L = L2)))
  }

  # Run 20 EM updates again, this time using multithreaded computations.
  nc      <- 4
  F3      <- F
  L3      <- L
  loglik3 <- rep(0,n)
  dev3    <- rep(0,n)
  setThreadOptions(numThreads = nc)
  for (i in 1:n) {
    F3         <- pnmfem_update_factors(X,F3,L3,nc = nc)
    L3         <- pnmfem_update_loadings(X,F3,L3,nc = nc)
    loglik3[i] <- sum(loglik_poisson_nmf(X,list(F = F3,L = L3)))
    dev3[i]    <- sum(deviance_poisson_nmf(X,list(F = F3,L = L3)))
  }
  
  # Store the counts as a sparse matrix.
  X <- as(X,"dgCMatrix")
  
  # Run 20 EM updates a third time, this time using the sparse counts
  # matrix.
  F4      <- F
  L4      <- L
  loglik4 <- rep(0,n)
  dev4    <- rep(0,n)
  for (i in 1:n) {
    F4         <- pnmfem_update_factors(X,F4,L4)
    L4         <- pnmfem_update_loadings(X,F4,L4)
    loglik4[i] <- sum(loglik_poisson_nmf(X,list(F = F4,L = L4)))
    dev4[i]    <- sum(deviance_poisson_nmf(X,list(F = F4,L = L4)))
  }

  # Run 20 EM updates one last time, using the sparse counts matrix,
  # and using multithreaded computations.
  F5      <- F
  L5      <- L
  loglik5 <- rep(0,n)
  dev5    <- rep(0,n)
  for (i in 1:n) {
    F5         <- pnmfem_update_factors(X,F5,L5)
    L5         <- pnmfem_update_loadings(X,F5,L5)
    loglik5[i] <- sum(loglik_poisson_nmf(X,list(F = F5,L = L5)))
    dev5[i]    <- sum(deviance_poisson_nmf(X,list(F = F5,L = L5)))
  }
  
  # All the updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(loglik1)
  expect_nondecreasing(loglik2)
  expect_nondecreasing(loglik3)
  expect_nondecreasing(loglik4)
  expect_nondecreasing(loglik5)
  expect_nonincreasing(dev1)
  expect_nonincreasing(dev2)
  expect_nonincreasing(dev3)
  expect_nonincreasing(dev4)
  expect_nonincreasing(dev5)
  
  # The updated factors and loadings should all be the same.
  expect_equivalent(F1,F2,tolerance = 1e-12)
  expect_equivalent(F1,F3,tolerance = 1e-12)
  expect_equivalent(F1,F4,tolerance = 1e-12)
  expect_equivalent(F1,F5,tolerance = 1e-12)
  expect_equivalent(L1,L2,tolerance = 1e-12)
  expect_equivalent(L1,L3,tolerance = 1e-12)
  expect_equivalent(L1,L4,tolerance = 1e-12)
  expect_equivalent(L1,L5,tolerance = 1e-12)

  # The likelihoods and deviances should all be the same.
  expect_equal(loglik1,loglik2,tolerance = 1e-10)
  expect_equal(loglik1,loglik3,tolerance = 1e-10)
  expect_equal(loglik1,loglik4,tolerance = 1e-10)
  expect_equal(loglik1,loglik5,tolerance = 1e-10)
  expect_equal(dev1,dev2,tolerance = 1e-10)
  expect_equal(dev1,dev3,tolerance = 1e-10)
  expect_equal(dev1,dev4,tolerance = 1e-10)
  expect_equal(dev1,dev5,tolerance = 1e-10)
})

