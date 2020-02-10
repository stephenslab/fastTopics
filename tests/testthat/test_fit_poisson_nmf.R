context("fit_poisson_nmf")

library(Matrix)
library(RcppParallel)
library(NNLM)

test_that("Version number in fit_poisson_nmf with verbose = TRUE is correct",{
  set.seed(1)
  X   <- generate_test_data(20,40,3)$X
  out <- capture.output(fit_poisson_nmf(X,k = 3))
  x   <- unlist(strsplit(out[2],"(",fixed = TRUE))[2]
  x   <- unlist(strsplit(x,")"))[1]
  expect_equal(paste("fastTopics",packageDescription("fastTopics")$Version),x)
})


test_that(paste("multiplicative and EM updates produce same result, and",
                "monotonically increase the likelihood"),{
                    
  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out  <- generate_test_data(80,100,3)
  X    <- out$X
  fit0 <- init_poisson_nmf(X,F = out$F,L = out$L)

  # Run 20 multiplicative updates and 20 EM updates.
  numiter <- 20
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "mu"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(nc = 1)))

  # Run 20 EM updates again, this time using multithreaded computations.
  nc <- 4
  capture.output(
    fit3 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(nc = nc)))
  
  # Store the counts as a sparse matrix.
  X <- as(X,"dgCMatrix")
  
  # Run 20 EM updates a third time, this time using the sparse counts
  # matrix.
  capture.output(
    fit4 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(nc = 1)))

  # Run 20 EM updates one last time, using the sparse counts matrix,
  # and using multithreaded computations.
  capture.output(
    fit5 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(nc = 1)))

  # All the updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)
  expect_nondecreasing(fit3$progress$loglik)
  expect_nondecreasing(fit4$progress$loglik)
  expect_nondecreasing(fit5$progress$loglik)
  expect_nonincreasing(fit1$progress$dev)
  expect_nonincreasing(fit2$progress$dev)
  expect_nonincreasing(fit3$progress$dev)
  expect_nonincreasing(fit4$progress$dev)
  expect_nonincreasing(fit5$progress$dev)
  
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
  expect_equal(fit1$progress$loglik,fit2$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$loglik,fit3$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$loglik,fit4$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$loglik,fit5$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$dev,fit2$progress$dev,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$dev,fit3$progress$dev,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$dev,fit4$progress$dev,tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$dev,fit5$progress$dev,tolerance = 1e-10,scale = 1)
})

test_that(paste("ccd and scd updates produce the same result, and",
                "monotonically increase the likelihood, when initialized",
                "\"close enough\" to a stationary point"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X

  # First find a reasonably good initialization using the SCD updates.
  capture.output(
    fit0 <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F = out$F,L = out$L),
              numiter = 10,method = "scd",control = list(nc = 1)))
  
  # Run 20 cyclic co-ordinate descent (CCD) updates.
  numiter <- 20
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd"),
                            control = list(nc = 1))

  # Run 20 sequential coordinate descent (SCD) updates.
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(nc = 1)))

  # Redo the SCD updates with a sparse matrix.
  Y <- as(X,"dgCMatrix")
  capture.output(
      fit3 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "scd",
                              control = list(nc = 1)))
  
  # Run 20 sequential coordinate descent (SCD) updates using the
  # multithreaded computations.
  nc <- 4
  capture.output(
    fit4 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(nc = nc)))

  # Redo the SCD updates one more time with a sparse matrix and
  # multithreaded computations.
  capture.output(
    fit5 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(nc = nc)))

  # Redo the CCD updates with a sparse matrix, with and without
  # multithreading.
  capture.output(
    fit6 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "ccd",
                            control = list(nc = 1)))
  capture.output(
    fit7 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd",
                            control = list(nc = nc)))
  capture.output(
    fit8 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "ccd",
                            control = list(nc = nc)))
  
  # All the updates should monotonically increase the likelihood and
  # decrease the deviance.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)
  expect_nondecreasing(fit3$progress$loglik)
  expect_nondecreasing(fit4$progress$loglik)
  expect_nondecreasing(fit5$progress$loglik)
  expect_nondecreasing(fit6$progress$loglik)
  expect_nondecreasing(fit7$progress$loglik)
  expect_nondecreasing(fit8$progress$loglik)
  expect_nonincreasing(fit1$progress$dev)
  expect_nonincreasing(fit2$progress$dev)
  expect_nonincreasing(fit3$progress$dev)
  expect_nonincreasing(fit4$progress$dev)
  expect_nonincreasing(fit5$progress$dev)
  expect_nonincreasing(fit6$progress$dev)
  expect_nonincreasing(fit7$progress$dev)
  expect_nonincreasing(fit8$progress$dev)

  # The updated factors and loadings should be nearly the same.
  expect_equal(fit1$F,fit2$F,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$F,fit6$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$F,fit7$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$F,fit8$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit3$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit4$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit5$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$L,fit2$L,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$L,fit6$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$L,fit7$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$L,fit8$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit3$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit4$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit5$L,tolerance = 1e-15,scale = 1)
  
  # The likelihoods and deviances should be nearly the same.
  expect_equal(fit1$progress$loglik,fit2$progress$loglik,
               tolerance = 1e-6,scale = 1)
  expect_equal(fit1$progress$loglik,fit6$progress$loglik,
               tolerance = 1e-12,scale = 1)
  expect_equal(fit1$progress$loglik,fit7$progress$loglik,
               tolerance = 1e-12,scale = 1)
  expect_equal(fit1$progress$loglik,fit8$progress$loglik,
               tolerance = 1e-12,scale = 1)
  expect_equal(fit2$progress$loglik,fit3$progress$loglik,
               tolerance = 1e-12,scale = 1)
  expect_equal(fit2$progress$loglik,fit5$progress$loglik,
               tolerance = 1e-12,scale = 1)
  expect_equal(fit2$progress$loglik,fit4$progress$loglik,
               tolerance = 1e-12,scale = 1)
  expect_equal(fit1$progress$dev,fit2$progress$dev,tolerance = 1e-6,scale = 1)
  expect_equal(fit1$progress$dev,fit6$progress$dev,tolerance = 1e-11,scale = 1)
  expect_equal(fit1$progress$dev,fit7$progress$dev,tolerance = 1e-11,scale = 1)
  expect_equal(fit1$progress$dev,fit8$progress$dev,tolerance = 1e-11,scale = 1)
  expect_equal(fit2$progress$dev,fit3$progress$dev,tolerance = 1e-11,scale = 1)
  expect_equal(fit2$progress$dev,fit4$progress$dev,tolerance = 1e-11,scale = 1)
  expect_equal(fit2$progress$dev,fit5$progress$dev,tolerance = 1e-11,scale = 1)
})

test_that(paste("When initialized \"close enough\" to a stationary point, the",
                "SCD, CCD and alt-SQP updates recover the same solution"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X

  # First find a reasonably good initialization using the SCD updates.
  capture.output(
    fit0 <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F = out$F,L = out$L),
              numiter = 40,method = "scd"))
  
  # Run 300 updates of the CCD, SCD and alt-SQP algorithms.
  numiter <- 300
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(numiter = 2)))
  capture.output(
    fit3 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "altsqp",
                            control = list(numiter = 4)))

  # Check that all three algorithms recover the same solution.
  expect_equal(min(fit1$progress$dev),min(fit2$progress$dev),
               tolerance = 1e-8,scale = 1)
  expect_equal(min(fit1$progress$dev),min(fit3$progress$dev),
               tolerance = 1e-8,scale = 1)
  expect_equal(max(fit1$progress$loglik),max(fit2$progress$loglik),
               tolerance = 1e-8,scale = 1)
  expect_equal(max(fit1$progress$loglik),max(fit3$progress$loglik),
               tolerance = 1e-8,scale = 1)
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
  fit0    <- init_poisson_nmf(X,F = out$F,L = out$L)
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "altsqp",
                            control = list(nc = 1,numiter = 4)))

  # Redo the alt-SQP updates with sparse X.
  Y <- as(X,"dgCMatrix")
  capture.output(
    fit2 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "altsqp",
                            control = list(nc = 1,numiter = 4)))

  # Redo the alt-SQP updates using multithreaded computations.
  capture.output(
    fit3 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,
              method = "altsqp",control = list(nc = 4,numiter = 4)))
  capture.output(
    fit4 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,
                            method = "altsqp",control = list(nc = 4,numiter = 4)))
  
  # The updated factors and loadings should be nearly the same in all
  # cases.
  expect_equivalent(fit1$F,fit2$F,tolerance = 1e-8,scale = 1)
  expect_equivalent(fit1$L,fit2$L,tolerance = 1e-8,scale = 1)
  expect_equivalent(fit1$L,fit3$L,tolerance = 1e-8,scale = 1)
  expect_equivalent(fit1$L,fit4$L,tolerance = 1e-8,scale = 1)
  
  # The alt-SQP updates should monotonically increase the likelihood
  # and decrease the deviance.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)
  expect_nondecreasing(fit3$progress$loglik)
  expect_nondecreasing(fit4$progress$loglik)
  expect_nonincreasing(fit1$progress$dev)
  expect_nonincreasing(fit2$progress$dev)
  expect_nonincreasing(fit3$progress$dev)
  expect_nonincreasing(fit4$progress$dev)
})

test_that(paste("All Poisson NMF updates recover same solution when",
                "F, L are initialized close to a stationary point"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X

  # Get a good initial estimate by running 50 iterations of the CCD
  # algorithm.
  capture.output(
    fit0 <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F = out$F,L = out$L),
                            numiter = 50,method = "ccd"))

  # Now improve the fit by running 100 iterations of the CCD, SCD and
  # alt-SQP algorithms.
  numiter <- 100
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                           control = list(numiter = 4)))
  capture.output(
    fit3 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "altsqp",
                            control = list(numiter = 4)))

  # All three solution estimates should have nearly the same likelihood
  # and deviance.
  expect_equal(max(fit1$progress$loglik),max(fit2$progress$loglik),
               tolerance = 1e-4,scale = 1)
  expect_equal(max(fit1$progress$loglik),max(fit3$progress$loglik),
               tolerance = 1e-4,scale = 1)
  expect_equal(min(fit1$progress$dev),min(fit2$progress$dev),
               tolerance = 1e-4,scale = 1)
  expect_equal(min(fit1$progress$dev),min(fit3$progress$dev),
               tolerance = 1e-4,scale = 1)
  
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
