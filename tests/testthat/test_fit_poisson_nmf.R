context("fit_poisson_nmf")

library(Matrix)

test_that(paste("betanmf and pnmfem updates produce same result, and",
                "increase the likelihood"),{
                    
  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- simulate_count_data(80,100,3,lmax = 0.3,fmax = 0.3)
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

  # Run 20 EM updates a second time, this time with a sparse counts
  # matrix.
  X <- as(X,"dgCMatrix")
  F3      <- F
  L3      <- L
  loglik3 <- rep(0,n)
  dev3    <- rep(0,n)
  for (i in 1:n) {
    F3         <- pnmfem_update_factors(X,F3,L3)
    L3         <- pnmfem_update_loadings(X,F3,L3)
    loglik3[i] <- sum(loglik_poisson_nmf(X,list(F = F3,L = L3)))
    dev3[i]    <- sum(deviance_poisson_nmf(X,list(F = F3,L = L3)))
  }
  
  # The updated factors and loadings should be the same.
  expect_equivalent(F1,F2,tolerance = 1e-12)
  expect_equivalent(F1,F3,tolerance = 1e-12)
  expect_equivalent(L1,L2,tolerance = 1e-12)
  expect_equivalent(L1,L3,tolerance = 1e-12)
  expect_equal(loglik1,loglik2,tolerance = 1e-10)
  expect_equal(loglik1,loglik3,tolerance = 1e-10)
  expect_equal(dev1,dev2,tolerance = 1e-10)
  expect_equal(dev1,dev3,tolerance = 1e-10)

  # The updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(loglik1)
  expect_nondecreasing(loglik2)
  expect_nonincreasing(dev1)
  expect_nonincreasing(dev2)
})

