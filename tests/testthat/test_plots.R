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
  p1 <- plot_loglik_vs_rank(list(fit2,fit3,fit5,fit10))
  p2 <- plot_loglik_vs_rank(lapply(list(fit2,fit3,fit5,fit10),
                                   poisson2multinom))
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")
})

test_that("Test that pca_plot and pca_hexbin_plot work",{
  set.seed(1)
  k <- 3
  X <- simulate_toy_gene_data(n = 400,m = 40,k = k,s = 1000)$X
  capture.output(fit1 <- fit_poisson_nmf(X,k = k,numiter = 100,
                                         control = list(extrapolate = TRUE)))
  fit2 <- poisson2multinom(fit1)
  
  # Test pca_plot.
  p1 <- pca_plot(fit1)
  p2 <- pca_plot(fit2)
  p3 <- pca_hexbin_plot(fit2)
  p4 <- pca_hexbin_plot(fit2)
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")
  expect_s3_class(p3,"ggplot")
  expect_s3_class(p4,"ggplot")

  # Test the other variants of pca_plot.
  y <- factor(apply(fit2$L,1,which.max))
  levels(y) <- paste0("k",1:k)
  p5 <- pca_plot(fit1,fill = "none")
  p6 <- pca_plot(fit1,fill = fit2$L[,1])
  p7 <- pca_plot(fit1,fill = y)
  expect_s3_class(p5,"ggplot")
  expect_s3_class(p6,"ggplot")
  expect_s3_class(p7,"ggplot")
})

test_that("Test that other plotting functions work",{
  set.seed(1)
  dat <- generate_test_data(200,100,3)
  X   <- dat$X
  capture.output(fit0 <- init_poisson_nmf(X,k = 3))
  capture.output(
    fit1 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "scd",
                            control = list(extrapolate = TRUE)))
  capture.output(
    fit2 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "em",
                            control = list(extrapolate = TRUE)))

  # Test plot_progress.
  plot_progress(list(scd = fit1,em = fit2),y = "loglik")
  plot_progress(list(scd = fit1,em = fit2),y = "dev")
  plot_progress(list(scd = fit1,em = fit2),y = "res")
  plot_progress(list(scd = poisson2multinom(fit1),em = poisson2multinom(fit2)))
  
  # Test loadings_plot.
  x  <- factor(sample(1:4,200,replace = TRUE))
  p1 <- loadings_plot(fit1,x)
  p2 <- loadings_plot(poisson2multinom(fit1),x)
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")

  # Test tsne_plot and umap_plot.
  capture.output(p1 <- tsne_plot(fit1,fill = "loading"))
  capture.output(p2 <- tsne_plot(fit2,fill = "loading"))
  capture.output(p3 <- umap_plot(fit1,fill = "loading",verbose = FALSE))
  capture.output(p4 <- umap_plot(fit2,fill = "loading",verbose = FALSE))
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")
  expect_s3_class(p3,"ggplot")
  expect_s3_class(p4,"ggplot")

  # Test structure_plot.
  grouping <- factor(apply(poisson2multinom(fit1)$L,1,which.max))
  capture.output(y <- drop(tsne_from_topics(poisson2multinom(fit1),dims = 1)))
  capture.output(p1 <- structure_plot(fit1))
  capture.output(p2 <- structure_plot(fit1,grouping = grouping,gap = 5))
  capture.output(p3 <- structure_plot(poisson2multinom(fit1)$L))
  capture.output(p4 <- structure_plot(fit1$L))
  capture.output(p5 <- structure_plot(fit1,loadings_order = order(y)))
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")
  expect_s3_class(p3,"ggplot")
  expect_s3_class(p4,"ggplot")
  expect_s3_class(p5,"ggplot")

  # Test the "plot" S3 method (which creates a Structure plot).
  fit2 <- poisson2multinom(fit1)
  capture.output(p1 <- plot(fit1))
  capture.output(p2 <- plot(fit2))
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")
})
