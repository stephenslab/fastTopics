context("de_analysis")

test_that("le_diff_rcpp gives the same result as le.diff",{
  set.seed(1)
  n  <- 1000
  m  <- 20
  X  <- matrix(rnorm(n*m),n,m)
  Y1 <- t(apply(X,1,le.diff))
  Y2 <- le_diff_rcpp(X)
  expect_equal(Y1,Y2,scale = 1,tolerance = 1e-15)

  # Next check the calculations for the special case of two elements.
  m  <- 2
  X  <- matrix(rnorm(n*m),n,m)
  Y1 <- t(apply(X,1,le.diff))
  Y2 <- le_diff_rcpp(X)
  expect_equal(Y1,Y2,scale = 1,tolerance = 1e-15)
})

test_that(paste("R and C++ implementations of simulate_posterior_poisson",
                "produce the same output, and the Monte Carlo estimate of",
                "the covariance matrix is roughly similar to the covariance",
                "matrix computed using Laplace's method"),{
  
  # Simulate a Poisson data set.
  set.seed(1)
  n  <- 200
  f1 <- 0.1
  f2 <- 1
  s  <- sample(10,n,replace = TRUE)
  q  <- runif(n)
  u  <- (1-q)*f1 + q*f2
  x  <- rpois(n,s*u)
  L  <- cbind(s*(1-q),s*q)

  # Fit the Poisson glm.
  fit <- fit_poisson_glm(x,L)

  # Compute the covariance of log(f) using Laplace's method.
  f <- fit$coef
  S <- compute_poisson_covariance(x,L,f)

  # Draw samples from the posterior using the R and C++
  # implementations of the random-walk Metropolis algorithm.
  set.seed(1)
  ns   <- 1e4
  out1 <- simulate_posterior_poisson(x,L,f,ns = ns,s = 0.3)
  set.seed(1)
  i <- which(x > 0)
  D <- matrix(rnorm(2*ns),ns,2)
  U <- matrix(runif(2*ns),ns,2)
  M <- matrix(sample(2,2*ns,replace = TRUE),ns,2) - 1
  out2 <- simulate_posterior_poisson_rcpp(x,L,f,D,U,M,0.3,1e-15)
  out3 <- simulate_posterior_poisson_sparse_rcpp(x[i],L[i,],colSums(L),f,
                                                 D,U,M,0.3,1e-15)
  out2$ar <- drop(out2$ar)
  out3$ar <- drop(out3$ar)
  
  # The outputs from the R and C++ implementations should be the same.
  expect_equal(out1,out2,scale = 1,tolerance = 1e-15)
  expect_equal(out1,out3,scale = 1,tolerance = 1e-15)
  
  # The Laplace and MCMC estimates of the covariance matrix should be
  # roughly similar.
  Smc <- cov(out1$samples)
  expect_equal(S,Smc,scale = 1,tolerance = 0.01)
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

test_that(paste("de_analysis with and without multithreading, using a",
                "sparse or dense counts matrix, with shrink.method =",
                "\"none\", produces the same result"),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 800
  m   <- 1000
  k   <- 4
  dat <- simulate_multinom_gene_data(n,m,k,sparse = FALSE)
  X   <- dat$X
  Y   <- as(X,"dgCMatrix")
  L   <- dat$L

  # Run de_analysis twice, using the single-threaded computations (nc
  # = 1) and again using the multithreaded computations (nc = 2). As
  # long as the sequence of pseudorandom numbers is the same, the
  # output should be the same.
  fit <- init_poisson_nmf(X,L = L,init.method = "random")
  for (lfc.stat in c("le","vsnull",paste0("k",1:4))) {
    set.seed(1)
    capture.output(de1 <- de_analysis(fit,X,lfc.stat = lfc.stat,
                                      shrink.method = "none",
                                      control = list(nc = 1)))
    set.seed(1)
    capture.output(de2 <- de_analysis(fit,X,lfc.stat = lfc.stat,
                                      shrink.method = "none",
                                      control = list(nc = 2)))
    set.seed(1)
    capture.output(de3 <- de_analysis(fit,Y,lfc.stat = lfc.stat,
                                      shrink.method = "none",
                                      control = list(nc = 1)))
    set.seed(1)
    capture.output(de4 <- de_analysis(fit,Y,lfc.stat = lfc.stat,
                                      shrink.method = "none",
                                      control = list(nc = 2)))
    expect_equal(de1,de2,scale = 1,tolerance = 1e-5)
    expect_equal(de1,de3,scale = 1,tolerance = 1e-5)
    expect_equal(de1,de4,scale = 1,tolerance = 1e-5)
  }
})

test_that(paste("de_analysis with and without multithreading, using a",
                "sparse or dense counts matrix, with shrink.method =",
                "\"ash\", produces the same result"),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 800
  m   <- 1000
  k   <- 4
  dat <- simulate_multinom_gene_data(n,m,k,sparse = FALSE)
  X   <- dat$X
  Y   <- as(X,"dgCMatrix")
  L   <- dat$L

  # Run de_analysis twice, using the single-threaded computations (nc
  # = 1) and again using the multithreaded computations (nc = 2). As
  # long as the sequence of pseudorandom numbers is the same, the
  # output should be the same.
  fit <- init_poisson_nmf(X,L = L,init.method = "random")
  for (lfc.stat in c("le","vsnull",paste0("k",1:4))) {
    set.seed(1)
    capture.output(de1 <- de_analysis(fit,X,lfc.stat = lfc.stat,
                                      shrink.method = "ash",
                                      control = list(nc = 1)))
    set.seed(1)
    capture.output(de2 <- de_analysis(fit,X,lfc.stat = lfc.stat,
                                      shrink.method = "ash",
                                      control = list(nc = 2)))
    set.seed(1)
    capture.output(de3 <- de_analysis(fit,Y,lfc.stat = lfc.stat,
                                      shrink.method = "ash",
                                      control = list(nc = 1)))
    set.seed(1)
    capture.output(de4 <- de_analysis(fit,Y,lfc.stat = lfc.stat,
                                      shrink.method = "ash",
                                      control = list(nc = 2)))
    expect_equal(de1,de2,scale = 1,tolerance = 1e-5)
    expect_equal(de1,de3,scale = 1,tolerance = 1e-5)
    expect_equal(de1,de4,scale = 1,tolerance = 1e-5)
  }
})

test_that(paste("diff_count_analysis with s = rowSums(X) closely recovers",
                "true probabilities (relative gene expression levels) when",
                "provided with the true topic proportions"),{
  skip("test needs to be revised")                 
                    
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

test_that(paste("Pairwise and \"least extreme\" LFC statistics are correct",
                "for k = 2 topics, with shrink.method = \"none\""),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 100
  m   <- 200
  k   <- 2
  s   <- 10^runif(n,-1,1)
  dat <- simulate_poisson_gene_data(n,m,k,s)
  X   <- dat$X
  L   <- dat$L

  # Remove all-zero columns.
  X <- X[,colSums(X) > 0]
  m <- ncol(X)

  # Fit a multinomial topic model with k = 2 topics.
  capture.output(fit <- fit_topic_model(X,k = 2))
  
  # Compute "pairwise" LFC statistics.
  capture.output(de1 <- de_analysis(fit,X,lfc.stat="k1",shrink.method="none"))
  capture.output(de2 <- de_analysis(fit,X,lfc.stat="k2",shrink.method="none"))

  # By definition, LFC(k) where k is the same as lfc.stat should be
  # zero, and the z-scores and p-values should be NA.
  expect_equivalent(de1$est[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$est[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$postmean[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$postmean[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$lower[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$lower[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$upper[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$upper[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$z[,1],rep(0,m))
  expect_equivalent(de2$z[,2],rep(0,m))
  expect_equivalent(de1$lpval[,1],rep(0,m))
  expect_equivalent(de2$lpval[,2],rep(0,m))
  
  # Compute "least extreme" LFC statistics.
  capture.output(de <- de_analysis(fit,X,lfc.stat="le",shrink.method="none"))

  # By definition, LFC(1) should always be the same as -LFC(2) when
  # there are only two topics. Other quantities follow similarly.
  expect_equal(de$est[,1],-de$est[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$postmean[,1],-de$postmean[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$lower[,1],-de$upper[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$upper[,1],-de$lower[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$z[,1],-de$z[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$lpval[,1],de$lpval[,2],scale = 1,tolerance = 1e-15)
})

test_that(paste("Pairwise and \"least extreme\" LFC statistics are correct",
                "for k = 2 topics, with shrink.method = \"ash\""),{

  # Simulate gene expression data.
  set.seed(1)
  n   <- 100
  m   <- 200
  k   <- 2
  s   <- 10^runif(n,-1,1)
  dat <- simulate_poisson_gene_data(n,m,k,s)
  X   <- dat$X
  L   <- dat$L

  # Remove all-zero columns.
  X <- X[,colSums(X) > 0]
  m <- ncol(X)

  # Fit a multinomial topic model with k = 2 topics.
  capture.output(fit <- fit_topic_model(X,k = 2))
  
  # Compute "pairwise" LFC statistics.
  capture.output(de1 <- de_analysis(fit,X,lfc.stat="k1",shrink.method="ash"))
  capture.output(de2 <- de_analysis(fit,X,lfc.stat="k2",shrink.method="ash"))

  # By definition, LFC(k) where k is the same as lfc.stat should be
  # zero, and the z-scores and p-values should be NA.
  expect_equivalent(de1$est[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$est[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$postmean[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$postmean[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$lower[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$lower[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$upper[,1],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de2$upper[,2],rep(0,m),scale = 1,tolerance = 1e-5)
  expect_equivalent(de1$z[,1],rep(0,m))
  expect_equivalent(de2$z[,2],rep(0,m))
  expect_equivalent(de1$lpval[,1],rep(0,m))
  expect_equivalent(de2$lpval[,2],rep(0,m))
  
  # Compute "least extreme" LFC statistics.
  capture.output(de <- de_analysis(fit,X,lfc.stat="le",shrink.method="ash"))

  # By definition, LFC(1) should always be the same as -LFC(2) when
  # there are only two topics. Other quantities follow similarly.
  expect_equal(de$est[,1],-de$est[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$postmean[,1],-de$postmean[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$lower[,1],-de$upper[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$upper[,1],-de$lower[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$z[,1],-de$z[,2],scale = 1,tolerance = 1e-15)
  expect_equal(de$lpval[,1],de$lpval[,2],scale = 1,tolerance = 1e-15)
})
