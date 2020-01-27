context("fit_poisson_nmf")

library(Matrix)
library(RcppParallel)
library(NNLM)

test_that(paste("betanmf and pnmfem updates produce same result, and",
                "monotonically increase the likelihood"),{
                    
  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

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
  expect_equivalent(fit1$F,fit2$F,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$F,fit3$F,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$F,fit4$F,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$F,fit5$F,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$L,fit2$L,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$L,fit3$L,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$L,fit4$L,tolerance = 1e-12,scale = 1)
  expect_equivalent(fit1$L,fit5$L,tolerance = 1e-12,scale = 1)

  # The likelihoods and deviances should all be the same.
  expect_equal(fit1$loglik,fit2$loglik,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$loglik,fit3$loglik,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$loglik,fit4$loglik,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$loglik,fit5$loglik,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$dev,fit2$dev,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$dev,fit3$dev,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$dev,fit4$dev,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$dev,fit5$dev,tolerance = 1e-10,scale = 1)
})

test_that(paste("ccd and scd updates produce the same result, and",
                "monotonically increase the likelihood, when initialized",
                "\"close enough\" to a stationary point"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # First find a reasonably good initialization using the SCD updates.
  fit0 <- iterate_updates(X,F,L,10,
                          function (X,F,L) t(scd_update_factors(X,L,t(F))),
                          function (X,F,L) scd_update_loadings(X,L,t(F)))
  F0   <- fit0$F
  L0   <- fit0$L
  
  # Run 20 cyclic co-ordinate descent (CCD) updates.
  numiter <- 20
  fit1 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) t(ccd_update_factors(X,L,t(F))),
                          function (X,F,L) ccd_update_loadings(X,L,t(F)))

  # Run 20 sequential coordinate descent (SCD) updates.
  fit2 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) t(scd_update_factors(X,L,t(F))),
                          function (X,F,L) scd_update_loadings(X,L,t(F)))

  # Redo the SCD updates with a sparse matrix.
  fit3 <- iterate_updates(as(X,"dgCMatrix"),F0,L0,numiter,
                          function (X,F,L) t(scd_update_factors(X,L,t(F))),
                          function (X,F,L) scd_update_loadings(X,L,t(F)))
  
  # Run 20 sequential coordinate descent (SCD) updates using the
  # multithreaded computations.
  nc <- 4
  setThreadOptions(numThreads = nc)
  fit4 <- iterate_updates(X,F0,L0,numiter,
            function (X,F,L) t(scd_update_factors(X,L,t(F),nc = nc)),
            function (X,F,L) scd_update_loadings(X,L,t(F),nc = nc))

  # Redo the SCD updates one more time with a sparse matrix and
  # multithreaded computations.
  # fit5 <- iterate_updates(as(X,"dgCMatrix"),F0,L0,numiter,
  #           function (X,F,L) t(scd_update_factors(X,L,t(F),nc = nc)),
  #           function (X,F,L) scd_update_loadings(X,L,t(F),nc = nc))
  
  # All the updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(fit1$loglik)
  expect_nondecreasing(fit2$loglik)
  expect_nondecreasing(fit3$loglik)
  expect_nondecreasing(fit4$loglik)
  # expect_nondecreasing(fit5$loglik)
  expect_nonincreasing(fit1$dev)
  expect_nonincreasing(fit2$dev)
  expect_nonincreasing(fit3$dev)
  expect_nonincreasing(fit4$dev)
  # expect_nonincreasing(fit5$dev)

  # The updated factors and loadings should be nearly the same.
  expect_equal(fit1$F,fit2$F,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$L,fit2$L,tolerance = 1e-8,scale = 1)
  expect_equal(fit2$F,fit3$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit3$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit4$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit4$L,tolerance = 1e-15,scale = 1)
  # expect_equal(fit2$F,fit5$F,tolerance = 1e-15,scale = 1)
  # expect_equal(fit2$L,fit5$L,tolerance = 1e-15,scale = 1)
  
  # The likelihoods and deviances should be nearly the same.
  expect_equal(fit1$loglik,fit2$loglik,tolerance = 1e-6,scale = 1)
  expect_equal(fit2$loglik,fit3$loglik,tolerance = 1e-12,scale = 1)
  expect_equal(fit2$loglik,fit4$loglik,tolerance = 1e-12,scale = 1)
  # expect_equal(fit2$loglik,fit5$loglik,tolerance = 1e-12,scale = 1)
  expect_equal(fit1$dev,fit2$dev,tolerance = 1e-6,scale = 1)
  expect_equal(fit2$dev,fit3$dev,tolerance = 1e-12,scale = 1)
  expect_equal(fit2$dev,fit4$dev,tolerance = 1e-12,scale = 1)
  # expect_equal(fit2$dev,fit5$dev,tolerance = 1e-12,scale = 1)
})

test_that(paste("When initialized \"close enough\" to a stationary point, the",
                "SCD, CCD and alt-SQP updates recover the same solution"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # First find a reasonably good initialization using the SCD updates.
  fit0 <- iterate_updates(X,F,L,40,
                          function (X,F,L) t(scd_update_factors(X,L,t(F))),
                          function (X,F,L) scd_update_loadings(X,L,t(F)))
  F0 <- fit0$F
  L0 <- fit0$L
  
  # Run 300 updates of the CCD, SCD and alt-SQP algorithms.
  numiter <- 300
  fit1 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) t(ccd_update_factors(X,L,t(F))),
                          function (X,F,L) ccd_update_loadings(X,L,t(F)))
  fit2 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) t(scd_update_factors(X,L,t(F),2)),
                          function (X,F,L) scd_update_loadings(X,L,t(F),2))
  fit3 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) altsqp_update_factors(X,F,L,4),
                          function (X,F,L) altsqp_update_loadings(X,F,L,4))

  # Check that all three algorithms recover the same solution.
  expect_equal(min(fit1$dev),min(fit2$dev),tolerance = 1e-8,scale = 1)
  expect_equal(min(fit1$dev),min(fit3$dev),tolerance = 1e-8,scale = 1)
  expect_equal(max(fit1$loglik),max(fit2$loglik),tolerance = 1e-8,scale = 1)
  expect_equal(max(fit1$loglik),max(fit3$loglik),tolerance = 1e-8,scale = 1)
})

test_that("scd updates and nnmf from NNLM package produce same result",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  k   <- 3
  out <- generate_test_data(80,100,k)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Run the SCD algorithm as it is implemented in the NNLM package.
  numiter <- 20
  fit1 <- suppressWarnings(NNLM::nnmf(X,k,init = list(W = L,H = t(F)),
                                      method = "scd",loss = "mkl",rel.tol = -1,
                                      max.iter = numiter,inner.max.iter = 4,
                                      verbose = 0,inner.rel.tol = -1))

  # Run 20 sequential coordinate descent (SCD) updates as they are
  # implemented in the fastTopics package.
  fit2 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) t(scd_update_factors(X,L,t(F),4)),
                          function (X,F,L) scd_update_loadings(X,L,t(F),4),
                          factors_first = FALSE)

  expect_equivalent(fit1$W,fit2$L,tolerance = 1e-12,scale = 1)
  expect_equivalent(t(fit1$H),fit2$F,tolerance = 1e-12,scale = 1)
})

test_that(paste("altsqp updates with dense and sparse matrices produce the",
                "same result, and monotonically increase the likelihood",
                "in all cases"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Run 20 sequential alternating SQP (alt-SQP) updates.
  numiter <- 20
  fit1 <- iterate_updates(X,F,L,numiter,
                          function (X,F,L) altsqp_update_factors(X,F,L,4),
                          function (X,F,L) altsqp_update_loadings(X,F,L,4))

  # Redo the alt-SQP updates with sparse X.
  fit2 <- iterate_updates(as(X,"dgCMatrix"),F,L,numiter,
                          function (X,F,L) altsqp_update_factors(X,F,L,4),
                          function (X,F,L) altsqp_update_loadings(X,F,L,4))

  # Redo the alt-SQP updates using multithreaded computations.
  control    <- altsqp_control_default()
  control$nc <- 4
  setThreadOptions(numThreads = control$nc)
  fit3 <- iterate_updates(X,F,L,numiter,
            function (X,F,L) altsqp_update_factors(X,F,L,4,control),
            function (X,F,L) altsqp_update_loadings(X,F,L,4,control))
  fit4 <- iterate_updates(as(X,"dgCMatrix"),F,L,numiter,
            function (X,F,L) altsqp_update_factors(X,F,L,4,control),
            function (X,F,L) altsqp_update_loadings(X,F,L,4,control))
  
  # The updated factors and loadings should be nearly the same in all
  # cases.
  expect_equivalent(fit1$F,fit2$F,tolerance = 1e-8,scale = 1)
  expect_equivalent(fit1$L,fit2$L,tolerance = 1e-8,scale = 1)
  expect_equivalent(fit1$L,fit3$L,tolerance = 1e-8,scale = 1)
  expect_equivalent(fit1$L,fit4$L,tolerance = 1e-8,scale = 1)
  
  # The alt-SQP updates should monotonically increase the likelihood
  # and decrease the deviance.
  expect_nondecreasing(fit1$loglik)
  expect_nondecreasing(fit2$loglik)
  expect_nondecreasing(fit3$loglik)
  expect_nondecreasing(fit4$loglik)
  expect_nonincreasing(fit1$dev)
  expect_nonincreasing(fit2$dev)
  expect_nonincreasing(fit3$dev)
  expect_nonincreasing(fit4$dev)
})

test_that(paste("All Poisson NMF updates recover same solution when",
                "F, L are initialized close to a stationary point"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Get a good initial estimate by running 50 iterations of the CCD
  # algorithm.
  fit0 <- iterate_updates(X,F,L,50,
                          function (X,F,L) t(ccd_update_factors(X,L,t(F))),
                          function (X,F,L) ccd_update_loadings(X,L,t(F)))
  F0   <- fit0$F
  L0   <- fit0$L

  # Now improve the fit by running 100 iterations of the CCD, SCD and
  # alt-SQP algorithms.
  numiter <- 100
  fit1 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) t(ccd_update_factors(X,L,t(F))),
                          function (X,F,L) ccd_update_loadings(X,L,t(F)))
  fit2 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) t(scd_update_factors(X,L,t(F),4)),
                          function (X,F,L) scd_update_loadings(X,L,t(F),4))
  fit3 <- iterate_updates(X,F0,L0,numiter,
                          function (X,F,L) altsqp_update_factors(X,F,L,4),
                          function (X,F,L) altsqp_update_loadings(X,F,L,4))

  # All three solution estimates should have nearly the same likelihood
  # and deviance.
  expect_equal(max(fit1$loglik),max(fit2$loglik),tolerance = 1e-4,scale = 1)
  expect_equal(max(fit1$loglik),max(fit3$loglik),tolerance = 1e-4,scale = 1)
  expect_equal(min(fit1$dev),min(fit2$dev),tolerance = 1e-4,scale = 1)
  expect_equal(min(fit1$dev),min(fit3$dev),tolerance = 1e-4,scale = 1)
  
  # All three algorithms should arrive at nearly the same solution.
  fit1 <- with(fit1,rescale.factors(F,L))
  fit2 <- with(fit2,rescale.factors(F,L))
  fit3 <- with(fit3,rescale.factors(F,L))
  expect_equivalent(fit1$F,fit2$F,tolerance = 1e-3,scale = 1)
  expect_equivalent(fit1$L,fit2$L,tolerance = 1e-3,scale = 1)
  expect_equivalent(fit1$F,fit3$F,tolerance = 1e-3,scale = 1)
  expect_equivalent(fit1$L,fit3$L,tolerance = 1e-3,scale = 1)
  expect_equivalent(fit2$F,fit3$F,tolerance = 1e-3,scale = 1)
  expect_equivalent(fit2$L,fit3$L,tolerance = 1e-3,scale = 1)
})
