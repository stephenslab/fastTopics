context("loglik")

test_that(paste("R and Rcpp versions of cost function return same result",
                "for sparse and dense matrix"),{

  # Generate a data set.
  library(Matrix)
  set.seed(1)
  out <- simulate_count_data(10,8,3)
  X   <- out$X
  F   <- out$F
  L   <- out$L
  Y   <- as(X,"dgCMatrix")

  # Compute the loss function.
  f1 <- cost(X,L,t(F),e,version = "R")
  f2 <- cost(X,L,t(F),e,version = "Rcpp")
  f3 <- cost(Y,L,t(F),e,version = "R")
  f4 <- cost(Y,L,t(F),e,version = "Rcpp")
  
  # The cost function calculations should all give the same result.
  expect_equal(f1,f2)
  expect_equal(f1,f3)
  expect_equal(f1,f4)
})

test_that(paste("loglik_poisson_topic_model gives correct result for",
                "sparse and dense matrix"),{

  # Generate a data set.
  library(Matrix)
  set.seed(1)
  out <- simulate_count_data(10,8,3)
  X   <- out$X
  fit <- out[c("F","L")]
  Y   <- with(fit,tcrossprod(L,F))

  # Compute the log-likelikhood.
  f1 <- rowSums(dpois(X,Y,log = TRUE))
  f2 <- loglik_poisson_topic_model(X,fit)
  f3 <- loglik_poisson_topic_model(as(X,"dgCMatrix"),fit)
  names(f1) <- rownames(X)

  # The likelihood calculations should all be the same.
  expect_equal(f1,f2)
  expect_equal(f1,f3)
})

test_that(paste("loglik_multinom_topic_model gives correct result for",
                "sparse and dense matrix"),{

  # Generate a data set.
  library(Matrix)
  set.seed(1)
  out <- simulate_count_data(10,8,3)
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
})
