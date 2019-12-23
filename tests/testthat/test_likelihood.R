context("likelihood")

test_that(paste("R and Rcpp versions of cost function return same result",
                "for sparse and dense matrix"),{
  library(Matrix)
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
  library(Matrix)
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 1:4)  {
      
    # Generate a data set.
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    fit <- out[c("F","L")]
    Y   <- with(fit,tcrossprod(L,F))

    # Compute the log-likelikhood.
    f1 <- rowSums(dpois(X,Y,log = TRUE))
    f2 <- loglik_poisson_nmf(X,fit)
    f3 <- loglik_poisson_nmf(as(X,"dgCMatrix"),fit)
    names(f1) <- rownames(X)

    # The likelihood calculations should all be the same.
    expect_equal(f1,f2)
    expect_equal(f1,f3)
  }
})

test_that(paste("loglik_multinom_topic_model gives correct result for",
                "sparse and dense matrix"),{
  library(Matrix)
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 2:4) {
      
    # Generate a data set.
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    fit <- poisson2multinom(out[c("F","L")])
    Y   <- with(fit,tcrossprod(L,F))

    # Compute the log-likelikhood.
    f1        <- rep(0,10)
    names(f1) <- rownames(X)
    for (i in 1:10)
      f1[i] <- dmultinom(X[i,],prob = Y[i,],log = TRUE)
    f2 <- loglik_multinom_topic_model(X,fit)
    f3 <- loglik_multinom_topic_model(as(X,"dgCMatrix"),fit)

    # The likelihood calculations should all be the same.
    expect_equal(f1,f2)
    expect_equal(f1,f3)
  }
})

test_that(paste("deviance_poisson_topic_nmf gives correct result for sparse",
                "and dense matrix"),{
  library(Matrix)
  set.seed(1)

  # Repeat the tests for a few different choices of k.
  for (k in 1:4) {
      
    # Generate a data set.  
    out <- simulate_count_data(10,8,k)
    X   <- out$X
    fit <- out[c("F","L")]
    Y   <- with(fit,tcrossprod(L,F))

    # Compute the deviances.
    d1 <- rowSums(stats::poisson()$dev.resids(X,Y,1))
    d2 <- deviance_poisson_nmf(X,fit)
    d3 <- deviance_poisson_nmf(as(X,"dgCMatrix"),fit)

    # The deviance calculations should all be the same.
    expect_equal(d1,d2)
    expect_equal(d1,d3)
  }
})
