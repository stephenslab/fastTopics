context("likelihood")

test_that(paste("R and Rcpp versions of cost function return same result",
                "for sparse and dense matrix"),{
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 1:4)  {
      
    # Generate a data set.
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    F   <- out$F
    L   <- out$L
    Y   <- as(X,"dgCMatrix")

    # Compute the loss function.
    f1 <- cost(X,L,t(F),version = "R")
    f2 <- cost(X,L,t(F),version = "Rcpp")
    f3 <- cost(Y,L,t(F),version = "R")
    f4 <- cost(Y,L,t(F),version = "Rcpp")
  
    # The cost function calculations should all give the same result.
    expect_equal(f1,f2)
    expect_equal(f1,f3)
    expect_equal(f1,f4)
  }
})

test_that(paste("loglik_poisson_nmf gives correct result for sparse and",
                "dense matrix"),{
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 1:4)  {
      
    # Generate a data set.
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    fit <- out[c("F","L")]
    class(fit) <- c("poisson_nmf_fit","list")
    
    # Compute the log-likelikhood.
    f1 <- loglik_poisson_nmf_with_dpois(X,fit)
    f2 <- loglik_poisson_nmf(X,fit,e = 0)
    f3 <- loglik_poisson_nmf(as(X,"dgCMatrix"),fit,e = 0)
    names(f1) <- rownames(X)

    # The likelihood calculations should all be the same.
    expect_equal(f1,f2)
    expect_equal(f1,f3)
  }
})

test_that(paste("loglik_multinom_topic_model gives correct result for",
                "sparse and dense matrix"),{
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 2:4) {
      
    # Generate a data set.
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    fit <- out[c("F","L")]
    class(fit) <- c("poisson_nmf_fit","list")

    # Compute the log-likelikhood.
    f1 <- loglik_multinom_topic_model_with_dmultinom(X,poisson2multinom(fit))
    f2 <- loglik_multinom_topic_model(X,fit,e = 0)
    f3 <- loglik_multinom_topic_model(as(X,"dgCMatrix"),fit,e = 0)

    # Compute the multinomial log-likelihood a different way: first
    # compute the Poisson log-likelihood, then subtract out the
    # likelihood for the "size factors".
    f4 <- loglik_poisson_nmf(X,fit,e = 0) - loglik_size_factors(X,fit$F,fit$L)
    
    # The likelihood calculations should all be the same.
    expect_equal(f1,f2)
    expect_equal(f1,f3)
    expect_equal(f1,f4)
  }
})

test_that(paste("deviance_poisson_topic_nmf gives correct result for sparse",
                "and dense matrix"),{
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 1:4) {
      
    # Generate a data set.  
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    fit <- out[c("F","L")]
    class(fit) <- c("poisson_nmf_fit","list")

    # Compute the deviances.
    d1 <- deviance_poisson_nmf_with_poisson(X,fit)
    d2 <- deviance_poisson_nmf(X,fit,e = 0)
    d3 <- deviance_poisson_nmf(as(X,"dgCMatrix"),fit,e = 0)

    # The deviance calculations should all be the same.
    expect_equal(d1,d2)
    expect_equal(d1,d3)
  }
})

test_that("poisson_nmf_kkt gives same result for sparse and dense matrices",{

  # Generate a data set.  
  set.seed(1)
  out <- simulate_count_data(10,8,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L

  # Compute the KKT residuals, and check that they are the same.
  out1 <- poisson_nmf_kkt(X,F,L)
  out2 <- poisson_nmf_kkt(as(X,"dgCMatrix"),F,L)
  expect_equal(out1,out2,tolerance = 1e-15,scale = 1)
})

test_that(paste("loglik_poisson_nmf and loglik_multinom_topic_model give",
                "the same result when provided with either a Poisson NMF",
                "or multinomial topic model fit"),{

  # Generate a data set.
  set.seed(1)
  out  <- simulate_count_data(10,8,3)
  X    <- out$X
  fit1 <- out[c("F","L")]
  class(fit1) <- c("poisson_nmf_fit","list")

  # Compute the Poisson NMF and multinomial topic model log-likelikhoods.
  fit2 <- poisson2multinom(fit1)
  f1 <- loglik_poisson_nmf(X,fit1,e = 0)
  f2 <- loglik_poisson_nmf(X,fit2,e = 0)
  f3 <- loglik_multinom_topic_model(X,fit1,e = 0)
  f4 <- loglik_multinom_topic_model(X,fit2,e = 0)

  # The Poisson NMF and multinomial topic model likelihood
  # calculations should be the same.
  expect_equal(f1,f2)
  expect_equal(f3,f4)
})

test_that(paste("loglik_poisson_nmf and loglik_multinom_topic_model plus",
                "size factor likelihoods are the same"),{

  # Generate a data set.
  set.seed(1)
  out  <- simulate_count_data(10,8,3)
  X    <- out$X

  # Fit a Poisson NMF model.
  capture.output(fit1 <- fit_poisson_nmf(X,k = 3,init.method = "random",
                                         numiter = 10))
  
  # Compute the Poisson NMF likelihoods.
  f1   <- loglik_poisson_nmf(X,fit1,e = 0)

  # Compute the multinomial topic model likelikhoods, and check that
  # they give the same result after adding the likelihoods for the
  # size factors.
  fit2 <- poisson2multinom(fit1)
  f2   <- loglik_multinom_topic_model(X,fit2,e = 0) +
          dpois(rowSums(X),fit2$s,log = TRUE)
  expect_equal(f1,f2,scale = 1,tolerance = 1e-14)
})

test_that("ldpois works for non-integer x",{
  set.seed(1)
  n <- 100
  lambda <- abs(4*rnorm(n))
  x <- rpois(n,lambda)
  y1 <- ldpois(x,lambda)
  y2 <- ldpois(x + 0.001,lambda)
  expect_equal(y1,y2,scale = 1,tolerance = 0.01)
})

test_that("ldpois gives same result as dpois for integer x",{
  set.seed(1)
  n <- 100
  lambda <- abs(4*rnorm(n))
  x <- rpois(n,lambda)
  y1 <- dpois(x,lambda,log = TRUE)
  y2 <- ldpois(x,lambda)
  expect_equal(y1,y2,scale = 1,tolerance = 1e-15)
})
