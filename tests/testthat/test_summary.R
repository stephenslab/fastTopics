context("summary")

test_that("summary method and print.summary methods produce output",{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  out  <- generate_test_data(80,100,3)
  X    <- out$X

  # Fit a Poisson non-negative factorization.
  capture.output(fit <- fit_poisson_nmf(X,k = 3,numiter = 100))

  # Produce summaries of the model fit.
  expect_output(print(summary(fit)))
  expect_output(print(summary(poisson2multinom(fit))))
  expect_output(print(summary(fit),show.mixprops = TRUE,
                      show.topic.reps = TRUE))
  expect_output(print(summary(poisson2multinom(fit)),
                      show.size.factors = TRUE,
                      show.mixprops = TRUE,
                      show.topic.reps = TRUE))
})
