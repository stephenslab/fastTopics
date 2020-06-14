context("summary")

test_that("summary functions works",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out  <- generate_test_data(80,100,3)
  X    <- out$X

  # Fit a Poisson non-negative factorization.
  capture.output(fit <- fit_poisson_nmf(X,k = 3,numiter = 100))

  # Produce a summary of the model fit.
  expect_output(print(summary(fit)))
})
