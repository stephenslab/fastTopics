context("diffcount")

test_that("Add test(s) here",{

})

test_that(paste("All variants of fit_poisson_models should produce the",
                "same, or nearly the same, result"),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 100
  m   <- 200
  k   <- 4
  s   <- 10^runif(n,-1,1)
  dat <- simulate_poisson_gene_data(n,m,k,s)
  X   <- dat$X
  L   <- dat$L

  # Remove all-zero columns.
  X <- X[,colSums(X) > 0]

  # Compute the MLE of mu under Poisson model x ~ Pois(s*mu).
  mu <- colSums(X)/sum(s)

  # Add "pseudocounts" to the data.
  out <- add_pseudocounts(X,s*L,0.1)
  X   <- out$X
  L   <- out$L
  Y   <- as(X,"dgCMatrix")
  
  # Fit the univariate Poisson models using glm and scd.
  F1 <- fit_poisson_models(X,L,method = "glm")
  F2 <- fit_poisson_models(Y,L,method = "glm")
  F3 <- fit_poisson_models(X,L,method = "scd")
  F4 <- fit_poisson_models(Y,L,method = "scd")

  # All four implementations should produce the same, or nearly the
  # same, estimates of the model parameters.
  expect_equal(F1,F2,scale = 1,tolerance = 1e-15)
  expect_equal(F3,F4,scale = 1,tolerance = 1e-15)
  expect_equal(F1,F3,scale = 1,tolerance = 1e-5)
})

test_that(paste("diff_count_analysis with s = rowSums(X) closely recovers",
                "true probabilities (relative gene expression levels) when",
                "provided with the true topic proportions"),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 800
  m   <- 1000
  k   <- 4
  dat <- simulate_multinom_gene_data(n,m,k,sparse = TRUE)
  X   <- dat$X
  L   <- dat$L

  # Fit a Poisson model (approximating a binomial model) to each gene
  # and topic, and compute the log-fold change statistics.
  fit <- init_poisson_nmf(X,L = L,init.method = "random")
  ans <- diff_count_analysis(fit,X,verbose = FALSE)

  # The f1 estimates should be close to the multinomial probabilities
  # that were used to simulate the data.
  expect_equal(dat$F,ans$F1,scale = 1,tolerance = 1e-4)

  # Create volcano plots from the diff_count_analysis output.
  p1 <- volcano_plot(ans,y = "zscore")
  p2 <- volcano_plot(ans,y = "pvalue")
  p3 <- volcano_plotly(ans,k = 1,y = "zscore")
  p4 <- volcano_plotly(ans,k = 1,y = "pvalue")
  expect_s3_class(p1,"ggplot")
  expect_s3_class(p2,"ggplot")
  expect_s3_class(p3,"plotly")
  expect_s3_class(p4,"plotly")
})
