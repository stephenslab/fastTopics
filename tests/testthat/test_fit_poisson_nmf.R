context("fit_poisson_nmf")

library(Matrix)
library(RcppParallel)
library(NNLM)

test_that(paste("betanmf and pnmfem updates produce same result, and",
                "monotonically increase the likelihood"),{
                    
  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3,0.3,0.3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # There should be at least one row and one column in the counts
  # matrix containing exactly one nonzero count.
  expect_gt(sum(rowSums(X > 0) == 1),0)
  expect_gt(sum(colSums(X > 0) == 1),0)

  # Run 20 multiplicative updates.
  numiter <- 20
  fit1 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) t(betanmf_update_factors(X,L,t(F))),
                          function (X,F,L) betanmf_update_loadings(X,L,t(F)))

  # Run 20 EM updates.
  fit2 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) pnmfem_update_factors(X,F,L),
                          function (X,F,L) pnmfem_update_loadings(X,F,L))

  # Run 20 EM updates again, this time using multithreaded computations.
  nc <- 4
  setThreadOptions(numThreads = nc)
  fit3 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) pnmfem_update_factors(X,F,L,nc=nc),
                          function (X,F,L) pnmfem_update_loadings(X,F,L,nc=nc))
  
  # Store the counts as a sparse matrix.
  X <- as(X,"dgCMatrix")
  
  # Run 20 EM updates a third time, this time using the sparse counts
  # matrix.
  fit4 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) pnmfem_update_factors(X,F,L),
                          function (X,F,L) pnmfem_update_loadings(X,F,L))

  # Run 20 EM updates one last time, using the sparse counts matrix,
  # and using multithreaded computations.
  fit5 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) pnmfem_update_factors(X,F,L,nc=nc),
                          function (X,F,L) pnmfem_update_loadings(X,F,L,nc=nc))
  
  # All the updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(fit1$loglik)
  expect_nondecreasing(fit2$loglik)
  expect_nondecreasing(fit3$loglik)
  expect_nondecreasing(fit4$loglik)
  expect_nondecreasing(fit5$loglik)
  expect_nonincreasing(fit1$dev)
  expect_nonincreasing(fit2$dev)
  expect_nonincreasing(fit3$dev)
  expect_nonincreasing(fit4$dev)
  expect_nonincreasing(fit5$dev)
  
  # The updated factors and loadings should all be the same.
  expect_equivalent(fit1$F,fit2$F,tolerance = 1e-12)
  expect_equivalent(fit1$F,fit3$F,tolerance = 1e-12)
  expect_equivalent(fit1$F,fit4$F,tolerance = 1e-12)
  expect_equivalent(fit1$F,fit5$F,tolerance = 1e-12)
  expect_equivalent(fit1$L,fit2$L,tolerance = 1e-12)
  expect_equivalent(fit1$L,fit3$L,tolerance = 1e-12)
  expect_equivalent(fit1$L,fit4$L,tolerance = 1e-12)
  expect_equivalent(fit1$L,fit5$L,tolerance = 1e-12)

  # The likelihoods and deviances should all be the same.
  expect_equal(fit1$loglik,fit2$loglik,tolerance = 1e-10)
  expect_equal(fit1$loglik,fit3$loglik,tolerance = 1e-10)
  expect_equal(fit1$loglik,fit4$loglik,tolerance = 1e-10)
  expect_equal(fit1$loglik,fit5$loglik,tolerance = 1e-10)
  expect_equal(fit1$dev,fit2$dev,tolerance = 1e-10)
  expect_equal(fit1$dev,fit3$dev,tolerance = 1e-10)
  expect_equal(fit1$dev,fit4$dev,tolerance = 1e-10)
  expect_equal(fit1$dev,fit5$dev,tolerance = 1e-10)
})

test_that("ccd updates monotonically increase the likelihood",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3,0.3,0.3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Run 20 cyclic co-ordinate descent (CCD) updates.
  fit <- iterate_updates(X,F,L,20,
                         function (X,F,L) t(ccd_update_factors(X,L,t(F))),
                         function (X,F,L) ccd_update_loadings(X,L,t(F)))

  # The CCD updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(fit$loglik)
  expect_nonincreasing(fit$dev)
})

test_that("scd updates monotonically increase the likelihood",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  k   <- 3
  out <- generate_test_data(80,100,k,0.3,0.3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Run 20 sequential coordinate descent (SCD) updates.
  fit <- iterate_updates(X,F,L,20,
                         function (X,F,L) t(scd_update_factors(X,L,t(F),4)),
                         function (X,F,L) scd_update_loadings(X,L,t(F),4),
                         factors_first = FALSE)

  # The SCD updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(fit$loglik)
  expect_nonincreasing(fit$dev)
})

test_that("scd updates and nnmf from NNLM package produce same result",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  k   <- 3
  out <- generate_test_data(80,100,k,0.3,0.3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Run the SCD algorithm implemented in the NNLM package.
  numiter <- 20
  fit1 <- suppressWarnings(nnmf(X,k,init = list(W = L,H = t(F)),
                                method = "scd",loss = "mkl",rel.tol = 0,
                                max.iter = numiter,inner.max.iter = 4,
                                trace = 1,verbose = 2))

  # Run 20 sequential coordinate descent (SCD) updates.
  fit2 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) t(scd_update_factors(X,L,t(F),4)),
                          function (X,F,L) scd_update_loadings(X,L,t(F),4),
                          factors_first = FALSE)

  # Run 20 sequential coordinate descent (SCD) updates, this time
  # using the multithreaded computations.
  # TO DO.
  
  expect_equivalent(fit1$W,fit2$L,tolerance = 1e-12)
  expect_equivalent(t(fit1$H),fit2$F,tolerance = 1e-12)
})

