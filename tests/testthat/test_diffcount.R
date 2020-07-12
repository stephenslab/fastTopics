context("diffcount")

test_that(paste("All variants of fit_univar_poisson_models and",
                "compute_univar_poisson_zscores produce the same, or",
                "nearly the same, result"),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 100
  m   <- 200
  k   <- 4
  s   <- 10^runif(n,-1,1)
  dat <- simulate_poisson_gene_data(n,m,k,s)
  X   <- dat$X
  Y   <- as(X,"dgCMatrix")
  L   <- dat$L

  # Fit the univariate Poisson models using optim and EM.
  capture.output(fit1 <- fit_univar_poisson_models(X,L,s,method = "optim"))
  capture.output(fit2 <- fit_univar_poisson_models(X,L,s,method = "em"))
  capture.output(fit3 <- fit_univar_poisson_models(X,L,s,method = "em-rcpp"))
  capture.output(fit4 <- fit_univar_poisson_models(Y,L,s,method = "em-rcpp"))

  # All four implementations should produce the same, or nearly the
  # same, estimates of the model parameters.
  expect_equal(fit1$F0,fit2$F0,scale = 1,tolerance = 1e-4)
  expect_equal(fit1$F1,fit2$F1,scale = 1,tolerance = 1e-4)
  expect_equal(fit2$F0,fit3$F0,scale = 1,tolerance = 1e-15)
  expect_equal(fit2$F1,fit3$F1,scale = 1,tolerance = 1e-15)
  expect_equal(fit2$F0,fit4$F0,scale = 1,tolerance = 1e-15)
  expect_equal(fit2$F1,fit4$F1,scale = 1,tolerance = 1e-15)

  # Compute the log-fold change statistics.
  F0   <- fit2$F0
  F1   <- fit2$F1
  out1 <- compute_univar_poisson_zscores(X,L,F0,F1,s)
  out2 <- compute_univar_poisson_zscores_fast(X,L,F0,F1,s)
  out3 <- compute_univar_poisson_zscores_fast(Y,L,F0,F1,s)

  # The three different implementations should produce the same, or
  # nearly the same, results.
  out1["se"] <- NULL
  out2["se"] <- NULL
  out3["se"] <- NULL
  expect_equal(out1,out2,scale = 1,tolerance = 1e-8)
  expect_equal(out1,out3,scale = 1,tolerance = 1e-8)
})
