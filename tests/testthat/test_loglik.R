context("loglik")

test_that("loglik.poisson gives same result for sparse and dense matrix",{
  library(Matrix)
  set.seed(1)
  A   <- matrix(runif(20),4,5) > 0.5
  X   <- matrix(0:19,4,5) * A
  Y   <- as(X,"dgCMatrix")
  fit <- list(F = matrix(0:9,5,2),
              L = matrix(0:7,4,2))
  f1 <- loglik.poisson(X,fit)
  f2 <- loglik.poisson(Y,fit)
  expect_equal(f1,f2)
})

test_that("loglik.poisson gives same result for sparse and dense matrix",{
  library(Matrix)
  set.seed(1)
  A   <- matrix(runif(20),4,5) > 0.5
  X   <- matrix(0:19,4,5) * A
  Y   <- as(X,"dgCMatrix")
  fit <- list(F = matrix(0:9,5,2),
              L = matrix(0:7,4,2))
  fit <- poisson2multinom(fit)
  f1  <- loglik.multinom(X,fit)
  f2  <- loglik.multinom(Y,fit)
  expect_equal(f1,f2)
})

test_that("R and Rcpp versions of cost function return same result",{
  e  <- 1e-8
  A  <- matrix(runif(20),4,5) > 0.5
  X  <- matrix(0:19,4,5) * A
  F  <- matrix(0:9,5,2)
  L  <- matrix(0:7,4,2)
  f1 <- cost(X,L,t(F),e,"R")
  f2 <- cost(X,L,t(F),e,"Rcpp")
  expect_equal(f1,f2)

  # Next, check the calculations when X is sparse.
  X  <- as(X,"dgCMatrix")
  f1 <- cost(X,L,t(F),e,"R")
  f2 <- cost(X,L,t(F),e,"Rcpp")
  expect_equal(f1,f2)
})
