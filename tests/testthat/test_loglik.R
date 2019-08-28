context("loglik")

test_that(paste("loglik.multinom gives correct value for sparse matrix",
                "and for dense matrix"),{
  library(Matrix)
  set.seed(1)
  e   <- 1e-15
  A   <- matrix(runif(20),4,5) > 0.5
  X   <- matrix(0:19,4,5) * A
  fit <- list(F = matrix(0:9,5,2),
              L = matrix(0:7,4,2))
  fit <- poisson2multinom(fit)
  Y   <- with(fit,tcrossprod(L,F))
  f1  <- sum(X*log(Y + e))
  f2  <- loglik.multinom(X,fit,e)
  f3  <- loglik.multinom(as(X,"dgCMatrix"),fit,e)
  expect_equal(f1,f2)
  expect_equal(f1,f3)
})

test_that("loglik.poisson gives correct value",{
  e   <- 1e-15
  X   <- matrix(0:19,4,5)
  F   <- matrix(0:9,5,2)
  L   <- matrix(0:7,4,2)
  fit <- list(F = F,L = L)
  Y   <- tcrossprod(L,F)
  f1  <- sum(X*log(Y + e) - Y)
  f2  <- loglik.poisson(X,fit,e)
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
  set.seed(1)
  e  <- 1e-8
  A  <- matrix(runif(20),4,5) > 0.5
  X  <- matrix(0:19,4,5) * A
  F  <- matrix(0:9,5,2)
  L  <- matrix(0:7,4,2)
  f1 <- cost(X,L,t(F),e,version = "R")
  f2 <- cost(X,L,t(F),e,version = "Rcpp")
  expect_equal(f1,f2)

  # Next, check the calculations when X is sparse.
  X  <- as(X,"dgCMatrix")
  f1 <- cost(X,L,t(F),e,version = "R")
  f2 <- cost(X,L,t(F),e,version = "Rcpp")
  expect_equal(f1,f2)
})
