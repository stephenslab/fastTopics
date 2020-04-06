context("fit_poisson_nmf")

test_that("Version number in fit_poisson_nmf with verbose = TRUE is correct",{
  set.seed(1)
  X   <- generate_test_data(20,40,3)$X
  out <- capture.output(fit_poisson_nmf(X,k = 3))
  x   <- unlist(strsplit(out[2],"(",fixed = TRUE))[2]
  x   <- unlist(strsplit(x,")"))[1]
  expect_equal(paste("fastTopics",packageDescription("fastTopics")$Version),x)
})

test_that(paste("fit$progress$loglik and fit$progress$dev agree with",
                "loglik_poisson_nmf and deviance_poisson_nmf"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out  <- generate_test_data(80,100,3)
  X    <- out$X
  fit0 <- init_poisson_nmf(X,F = out$F,L = out$L)

  # Run 50 SCD updates.
  capture.output(
    fit <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "scd",
                           control = list(extrapolate = TRUE,nc = 1)))

  # Check that the log-likelihood and deviance calculations agree with
  # loglik_poisson_nmf and deviance_poisson_nmf.
  expect_equal(sum(loglik_poisson_nmf(X,fit)),tail(fit$progress$loglik,n = 1))
  expect_equal(sum(deviance_poisson_nmf(X,fit)),tail(fit$progress$dev,n = 1))
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
  Y <- as(X,"dgCMatrix")
  
  # Run 20 EM updates a third time, this time using the sparse counts
  # matrix.
  capture.output(
    fit4 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(nc = 1)))

  # Run 20 EM updates one last time, using the sparse counts matrix,
  # and using multithreaded computations.
  capture.output(
    fit5 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "em",
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
  expect_equal(fit1$progress$dev,fit2$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$progress$dev,fit3$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$progress$dev,fit4$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$progress$dev,fit5$progress$dev,tolerance = 1e-8,scale = 1)
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
  expect_equal(fit1$F,fit2$F,tolerance = 1e-7,scale = 1)
  expect_equal(fit1$F,fit6$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$F,fit7$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$F,fit8$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit3$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit4$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$F,fit5$F,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$L,fit2$L,tolerance = 1e-7,scale = 1)
  expect_equal(fit1$L,fit6$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$L,fit7$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit1$L,fit8$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit3$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit4$L,tolerance = 1e-15,scale = 1)
  expect_equal(fit2$L,fit5$L,tolerance = 1e-15,scale = 1)
  
  # The likelihoods and deviances should be nearly the same.
  expect_equal(fit1$progress$loglik,fit2$progress$loglik,
               tolerance = 0.01,scale = 1)
  expect_equal(fit1$progress$loglik,fit6$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$loglik,fit7$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$loglik,fit8$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit2$progress$loglik,fit3$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit2$progress$loglik,fit5$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit2$progress$loglik,fit4$progress$loglik,
               tolerance = 1e-10,scale = 1)
  expect_equal(fit1$progress$dev,fit2$progress$dev,tolerance = 0.01,scale = 1)
  expect_equal(fit1$progress$dev,fit6$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$progress$dev,fit7$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit1$progress$dev,fit8$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit2$progress$dev,fit3$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit2$progress$dev,fit4$progress$dev,tolerance = 1e-8,scale = 1)
  expect_equal(fit2$progress$dev,fit5$progress$dev,tolerance = 1e-8,scale = 1)
})

test_that(paste("When initialized \"close enough\" to a stationary point, the",
                "SCD and CCD updates recover the same solution"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X

  # First find a reasonably good initialization using the SCD updates.
  capture.output(
    fit0 <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F = out$F,L = out$L),
              numiter = 40,method = "scd"))
  
  # Run 300 updates of the CCD and SCD algorithms.
  numiter <- 300
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(numiter = 2)))

  # Check that both algorithms recover the same solution.
  fits <- list(fit1 = fit1,fit2 = fit2)
  dat  <- compare_poisson_nmf_fits(fits)
  expect_equal(dat["fit1","dev"],dat["fit2","dev"],tolerance = 0.01,scale = 1)
  expect_equal(dat["fit1","loglik"],dat["fit2","loglik"],
               tolerance = 0.01,scale = 1)
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
  skip_if_not_installed("NNLM")
  numiter <- 20
  fit1 <- suppressWarnings(NNLM::nnmf(X,k,init = list(W = L,H = t(F)),
                                      method = "scd",loss = "mkl",rel.tol = -1,
                                      max.iter = numiter,inner.max.iter = 4,
                                      verbose = 0,inner.rel.tol = -1))

  # Run 20 sequential coordinate descent (SCD) updates as they are
  # implemented in the fastTopics package.
  fit2 <- iterate_updates(X,F,L,numiter,
    function (X,F,L) t(scd_update_factors(X,L,t(F),numiter = 4,runem = FALSE)),
    function (X,F,L) scd_update_loadings(X,L,t(F),numiter = 4,runem = FALSE),
    factors_first = FALSE)

  expect_equivalent(fit1$W,fit2$L,tolerance = 1e-12,scale = 1)
  expect_equivalent(t(fit1$H),fit2$F,tolerance = 1e-12,scale = 1)
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

  # Now improve the fit by running 100 iterations of the CCD and SCD
  # algorithms.
  numiter <- 100
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                           control = list(numiter = 4)))

  # All three solution estimates should have nearly the same likelihood
  # and deviance.
  fits <- list(fit1 = fit1,fit2 = fit2)
  dat  <- compare_poisson_nmf_fits(fits)
  expect_equal(dat["fit1","loglik"],dat["fit2","loglik"],
               tolerance = 0.01,scale = 1)
  expect_equal(dat["fit1","dev"],dat["fit2","dev"],tolerance = 0.01,scale = 1)
  
  # Both algorithms should arrive at nearly the same solution.
  expect_equivalent(fit1$F,fit2$F,tolerance = 1e-3,scale = 1)
  expect_equivalent(fit1$L,fit2$L,tolerance = 1e-3,scale = 1)
})

test_that(paste("Extrapolated updates achieve solutions that are as good",
                "or better than the un-extrapolated updates, at least when",
                "initialized close to a stationary point"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out <- generate_test_data(80,100,3)
  X   <- out$X

  # Get a good initial estimate by running 100 iterations of the EM
  # algorithm.
  capture.output(
    fit0 <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F = out$F,L = out$L),
                            numiter = 100,method = "em"))

  # Run 100 non-extrapolated updates using each of the methods.
  numiter <- 100
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "mu"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(numiter = 4)))
  capture.output(
    fit3 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd"))
  capture.output(
    fit4 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(numiter = 4)))

  # Run 100 extrapolated updates using each of the methods.
  capture.output(
    fit5 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "mu",
                            control = list(extrapolate = TRUE)))
  capture.output(
    fit6 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "em",
                            control = list(numiter = 4,extrapolate = TRUE)))
  capture.output(
    fit7 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "ccd",
                            control = list(extrapolate = TRUE)))
  capture.output(
    fit8 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = "scd",
                            control = list(numiter = 4,extrapolate = TRUE)))

  # The extrapolated updates should achieve likelihoods that are
  # higher (or at least not lower) than the non-extrapolated updates.
  fits <- list(fit1 = fit1,fit2 = fit2,fit3 = fit3,fit4 = fit4,
               fit5 = fit5,fit6 = fit6,fit7 = fit7,fit8 = fit8)
  dat  <- compare_poisson_nmf_fits(fits)
  expect_lte(dat["fit1","loglik"],dat["fit5","loglik"])
  expect_lte(dat["fit2","loglik"],dat["fit6","loglik"])
  expect_lte(dat["fit3","loglik"],dat["fit7","loglik"])
  expect_lte(dat["fit4","loglik"],dat["fit8","loglik"])

  # Likewise, the deviances for the extrapolated updates should be
  # less than (or no higher than) the deviances from the non-expected
  # updates.
  expect_gte(dat["fit1","dev"],dat["fit5","dev"])
  expect_gte(dat["fit2","dev"],dat["fit6","dev"])
  expect_gte(dat["fit3","dev"],dat["fit7","dev"])
  expect_gte(dat["fit4","dev"],dat["fit8","dev"])
})

test_that("Re-fitting yields same result as one call to fit_poisson_nmf",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  k    <- 3
  out  <- generate_test_data(80,100,k)
  X    <- out$X

  # Run 20 multiplicative updates.
  fit0 <- init_poisson_nmf(X,k = k)
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 20,method = "mu"))

  # Run 10 multiplicative updates, then "re-fit" the model by running
  # another 10 multiplicative updates.
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 10,method = "mu"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit2,numiter = 10,method = "mu"))

  # The two "fits" should be exactly the same (other than the timing).
  fit1$progress$timing <- NULL
  fit2$progress$timing <- NULL
  expect_equal(fit1,fit2)

  # Repeat this test for the EM updates.
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 20,method = "em"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 10,method = "em"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit2,numiter = 10,method = "em"))
  fit1$progress$timing <- NULL
  fit2$progress$timing <- NULL
  expect_equal(fit1,fit2)

  # Repeat this test for the CCD updates.
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 20,method = "ccd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 10,method = "ccd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit2,numiter = 10,method = "ccd"))
  fit1$progress$timing <- NULL
  fit2$progress$timing <- NULL
  expect_equal(fit1,fit2)

  # Repeat this test for the SCD updates.
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 20,method = "scd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 10,method = "scd"))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit2,numiter = 10,method = "scd"))
  fit1$progress$timing <- NULL
  fit2$progress$timing <- NULL
  expect_equal(fit1,fit2)
})

test_that("Fixed factors and loadings to not change (aside from rescaling)",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  n   <- 80
  m   <- 100
  k   <- 3
  out <- generate_test_data(n,m,k)
  X   <- out$X
  Y   <- as(X,"dgCMatrix")

  # Select which factors and loadings to update.
  n1 <- 75
  m1 <- 90
  i  <- sample(n,n1)
  j  <- sample(m,m1)
  i0 <- setdiff(1:n,i)
  j0 <- setdiff(1:m,j)
  
  # Run all variants of the EM algorithm.
  numiter <- 20
  fit0    <- init_poisson_nmf(X,F = out$F,L = out$L)
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,
                            update.factors = j,update.loadings = i,
                            method = "em",control = list(nc = 1)))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,
                            update.factors = j,update.loadings = i,
                            method = "em",control = list(nc = 4)))
  capture.output(
    fit3 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "em",
                            update.factors = j,update.loadings = i,
                            control = list(nc = 1)))
  capture.output(
    fit4 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = numiter,method = "em",
                            update.factors = j,update.loadings = i,
                            control = list(nc = 4)))

  # The factors that are not selected for updating should not change,
  # aside from a rescaling.
  compare_factors_ignoring_rescaling <- function (fit1, fit2, j)
    abs(max(apply(fit2$F[j,]/fit1$F[j,],2,diff)))
  expect_equal(compare_factors_ignoring_rescaling(fit0,fit1,j0),0,
               scale = 1,tolerance = 1e-15)
  expect_equal(compare_factors_ignoring_rescaling(fit0,fit2,j0),0,
               scale = 1,tolerance = 1e-15)
  expect_equal(compare_factors_ignoring_rescaling(fit0,fit3,j0),0,
               scale = 1,tolerance = 1e-15)
  expect_equal(compare_factors_ignoring_rescaling(fit0,fit4,j0),0,
               scale = 1,tolerance = 1e-15)
  
  # The loadings that are not selected for updating should not change,
  # aside from a rescaling.
  compare_loadings_ignoring_rescaling <- function (fit1, fit2, i)
    abs(max(apply(fit2$L[i,]/fit1$L[i,],2,diff)))
  expect_equal(compare_loadings_ignoring_rescaling(fit0,fit1,i0),0,
               scale = 1,tolerance = 1e-15)
  expect_equal(compare_loadings_ignoring_rescaling(fit0,fit2,i0),0,
               scale = 1,tolerance = 1e-15)
  expect_equal(compare_loadings_ignoring_rescaling(fit0,fit3,i0),0,
               scale = 1,tolerance = 1e-15)
  expect_equal(compare_loadings_ignoring_rescaling(fit0,fit4,i0),0,
               scale = 1,tolerance = 1e-15)

})
