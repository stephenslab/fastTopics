context("plots")

test_that("Test that the plotting functions work",{
  set.seed(1)
  dat  <- generate_test_data(80,100,3)
  X    <- dat$X
  capture.output(fit0 <- init_poisson_nmf(X,k = 3))
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "scd",
                            control = list(extrapolate = TRUE)))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "em",
                            control = list(extrapolate = TRUE)))

  # Test plot_progress_poisson_nmf.
  plot_progress_poisson_nmf(list(scd = fit1,em = fit2))
  
  # Test loadings_plot.
  x <- factor(sample(1:4,80,replace = TRUE))
  expect_s3_class(loadings_plot(poisson2multinom(fit1),1:3,x),"ggplot")
})
