context("fit_poisson_nmf")

test_that("fit_topic_model successfully fits a multinomial topic model",{

  # Generate a 80 x 100 sparse count matrix to factorize.
  set.seed(1)
  out <- simulate_count_data(80,100,k = 3,sparse = TRUE)
  X   <- out$X

  # Fit a multinomial topic model to these data.
  capture.output(fit <- fit_topic_model(X,k = 3))
  expect_s3_class(fit,"multinom_topic_model_fit")
  expect_s3_class(summary(fit),"summary.multinom_topic_model_fit")
})
