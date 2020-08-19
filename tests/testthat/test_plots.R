context("plots")

test_that("Test that plot_loglik_vs_rank works",{
  set.seed(1)
  dat <- generate_test_data(80,100,3)
  X   <- dat$X

  # Fit matrix factorizations with rank k = 2, 3, 5, 10.
  capture.output(fit2 <- fit_poisson_nmf(X,k = 2,numiter = 100))
  capture.output(fit3 <- fit_poisson_nmf(X,k = 3,numiter = 100))
  capture.output(fit5 <- fit_poisson_nmf(X,k = 5,numiter = 100))
  capture.output(fit10 <- fit_poisson_nmf(X,k = 10,numiter = 100))

  # Plot log-likelihood vs. rank.
  out <- plot_loglik_vs_rank(list(fit2,fit3,fit5,fit10))
  expect_s3_class(out,"ggplot")
})

test_that("Test that pca_plot works",{
  set.seed(1)
  X    <- simulate_toy_gene_data(n = 400,m = 40,k = 3,s = 1000)$X
  capture.output(fit1 <- fit_poisson_nmf(X,k = 3,numiter = 100,
                                         control = list(extrapolate = TRUE)))
  fit2 <- poisson2multinom(fit1)
  
  # Test pca_plot.
  out1 <- pca_plot(fit1)
  out2 <- pca_plot(fit2)
  expect_s3_class(out1,"ggplot")
  expect_s3_class(out2,"ggplot")
})

test_that("Test that other plotting functions work",{
  set.seed(1)
  dat <- generate_test_data(80,100,3)
  X   <- dat$X
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
  capture.output(out1 <- tsne_plot(fit1,color = "prop",perplexity = 10))
  capture.output(out2 <- tsne_plot(fit1,color = "loading",perplexity = 10))
  expect_s3_class(out1,"ggplot")
  expect_s3_class(out2,"ggplot")

  # Test structure_plot.
  g <- factor(apply(poisson2multinom(fit1)$L,1,which.max))
  capture.output(tsne <- tsne_from_topics(poisson2multinom(fit1),dims = 1,
                                          perplexity = 20))
  capture.output(out1 <- structure_plot(fit1,perplexity = 20))
  capture.output(out2 <- structure_plot(fit1,rows = order(tsne$Y)))
  out3 <- structure_plot(fit1,rows = order(tsne$Y),grouping = g,gap = 2)
  expect_s3_class(out1,"ggplot")
  expect_s3_class(out2,"ggplot")
  expect_s3_class(out3,"ggplot")

  # Test the "plot" S3 method (which creates a Structure plot).
  fit2 <- poisson2multinom(fit1)
  capture.output(out1 <- plot(fit1,perplexity = 20))
  capture.output(out2 <- plot(fit2,perplexity = 20))
  expect_s3_class(out1,"ggplot")
  expect_s3_class(out2,"ggplot")
})

