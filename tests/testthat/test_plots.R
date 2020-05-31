context("plots")

test_that("Test that the plotting functions work",{
  set.seed(1)
  dat <- generate_test_data(80,100,3)
  X   <- dat$X
  capture.output(
    fit <- fit_poisson_nmf(X,k = 3,numiter = 50,method = "scd",
                           control = list(extrapolate = TRUE)))

  # Test loadings_plot.
  x <- factor(sample(1:4,80,replace = TRUE))
  expect_s3_class(loadings_plot(poisson2multinom(fit),1:3,x),"ggplot")
})
