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
  x   <- factor(sample(1:4,80,replace = TRUE))
  out <- loadings_plot(poisson2multinom(fit1),x)
  expect_s3_class(out,"ggplot")

  # Test tsne_plot.
  capture.output(out1 <- tsne_plot(fit1,color = "probability",perplexity = 10))
  capture.output(out2 <- tsne_plot(fit1,color = "loading",perplexity = 10))
  expect_s3_class(out1,"ggplot")
  expect_s3_class(out2,"ggplot")

  # Test structure_plot.
  out <- structure_plot(fit1,perplexity = 20)
})
