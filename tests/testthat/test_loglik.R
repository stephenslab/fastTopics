context("loglik")

test_that("loglik.poisson gives same result for sparse and dense matrix",{
  set.seed(1)
  A   <- matrix(runif(20),4,5) > 0.5
  X   <- matrix(0:19,4,5) * A
  Y   <- as(X,"sparseMatrix")
  fit <- list(L = matrix(0:7,4,2),
              F = matrix(0:9,5,2))
  f1 <- loglik.poisson(X,fit)
  f2 <- loglik.poisson(Y,fit)
  expect_equal(f1,f2)
})
