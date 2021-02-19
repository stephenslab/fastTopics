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

  # Fit the univariate Poisson models using glm, optim and EM.
  fit1 <- fit_univar_poisson_models(X,L,s,method = "glm",verbose = FALSE)
  fit2 <- fit_univar_poisson_models(X,L,s,method = "optim",verbose = FALSE)
  fit3 <- fit_univar_poisson_models(X,L,s,method = "em",verbose = FALSE)
  fit4 <- fit_univar_poisson_models(X,L,s,method = "em-rcpp",verbose = FALSE)
  fit5 <- fit_univar_poisson_models(Y,L,s,method = "em-rcpp",verbose = FALSE)

  # All four implementations should produce the same, or nearly the
  # same, estimates of the model parameters.
  expect_equal(fit1$F0,fit2$F0,scale = 1,tolerance = 1e-5)
  expect_equal(fit1$F1,fit2$F1,scale = 1,tolerance = 1e-5)
  expect_equal(fit1$F0,fit3$F0,scale = 1,tolerance = 1e-5)
  expect_equal(fit1$F1,fit3$F1,scale = 1,tolerance = 1e-5)
  expect_equal(fit3$F0,fit4$F0,scale = 1,tolerance = 1e-8)
  expect_equal(fit3$F1,fit4$F1,scale = 1,tolerance = 1e-8)
  expect_equal(fit3$F0,fit5$F0,scale = 1,tolerance = 1e-8)
  expect_equal(fit3$F1,fit5$F1,scale = 1,tolerance = 1e-8)

  # Additionally, the likelihood calculations should be the same, or
  # nearly the same.
  expect_equal(fit1$loglik,fit2$loglik,scale = 1,tolerance = 1e-3)
  expect_equal(fit2$loglik,fit3$loglik,scale = 1,tolerance = 1e-6)
  expect_equal(fit2$loglik,fit4$loglik,scale = 1,tolerance = 1e-6)

  # Compute the log-fold change statistics.
  F0 <- fit2$F0
  F1 <- fit2$F1
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

test_that(paste("fit_univar_poisson_models produces same result as",
                "fit_univar_poisson_models_hard when L is a matrix",
                "with hard topic assignments only (that is, entirely",
                "zeros and ones)"),{
    
  # Simulate gene expression data.
  set.seed(1)
  n   <- 100
  m   <- 200
  k   <- 4
  s   <- 10^runif(n,-1,1)
  dat <- simulate_poisson_gene_data(n,m,k,s)
  X   <- dat$X
  Y   <- as(X,"dgCMatrix")
  L   <- force_hard_topic_assignments(dat$L)
      
  # Fit the univariate Poisson models using EM.
  fit1 <- fit_univar_poisson_models(X,L,s,method = "glm",e = 1e-8,
                                    verbose = FALSE)

  # Estimate the univariate Poisson model parameters assuming all the
  # topic proportions are zeros and ones.
  fit2 <- fit_univar_poisson_models_hard(X,L,s)
  fit3 <- fit_univar_poisson_models_hard(Y,L,s)

  # Both methods should produce the same parameter estimates, and the
  # same log-likelihood values.
  expect_equal(fit1$F0,fit2$F0,scale = 1,tolerance = 1e-8)
  expect_equal(fit1$F1,fit2$F1,scale = 1,tolernace = 1e-8)
  expect_equal(fit2$F0,fit3$F0,scale = 1,tolerance = 1e-15)
  expect_equal(fit2$F1,fit3$F1,scale = 1,tolernace = 1e-15)
})

test_that(paste("When all the topic proportions are exactly zero or exactly",
                "one, the result of diff_count_analysis is nearly the same",
                "as when all the topic proportions are all almost exactly",
                "zero or one"),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 100
  m   <- 200
  k   <- 4
  s   <- 10^runif(n,-1,1)
  dat <- simulate_poisson_gene_data(n,m,k,s)
  X   <- dat$X
  L   <- force_hard_topic_assignments(dat$L)
  clusters         <- factor(apply(dat$L,1,which.max))
  levels(clusters) <- paste0("k",1:k)

  # Remove all-zero columns.
  i <- which(colSums(X > 0) >= 1)
  X <- X[,i]
  
  # Fit a Poisson model (approximating a binomial model) to each gene
  # and topic, and compute the log-fold change statistics.
  Y    <- as(X,"dgCMatrix")
  fit1 <- init_poisson_nmf(X,L = rowSums(X)*L,init.method = "random",
                           control = list(minval = 1e-8))
  out1 <- diff_count_analysis(fit1,X,verbose = FALSE)
  out2 <- diff_count_clusters(clusters,X,verbose = FALSE)
  out3 <- diff_count_clusters(clusters,Y,verbose = FALSE)

  # The outputted statistics should be the same in all calls to
  # diff_count_analysis and diff_count_clusters (ignore the standard
  # errors, as these can be unstable).
  expect_equal(out1$F0,  out2$F0,  scale = 1,tolerance = 1e-10)
  expect_equal(out1$F0,  out3$F0,  scale = 1,tolerance = 1e-10)
  expect_equal(out1$F1,  out2$F1,  scale = 1,tolerance = 1e-10)
  expect_equal(out1$F1,  out3$F1,  scale = 1,tolerance = 1e-10)
  expect_equal(out1$beta,out2$beta,scale = 1,tolerance = 1e-5)
  expect_equal(out1$beta,out3$beta,scale = 1,tolerance = 1e-5)
  expect_equal(out1$Z,   out2$Z,   scale = 1,tolerance = 1e-5)
  expect_equal(out1$Z,   out3$Z,   scale = 1,tolerance = 1e-5)
  expect_equal(out1$pval,out2$pval,scale = 1,tolerance = 1e-5)
  expect_equal(out1$pval,out3$pval,scale = 1,tolerance = 1e-5)
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
