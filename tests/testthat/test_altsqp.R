context("altsqp")

test_that("altsqp gives same result for sparse and dense matrix",{
  library(Matrix)

  # Generate a 100 x 200 data matrix to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  k <- 3
  F <- matrix(0.75*runif(m*k),m,k)
  L <- matrix(0.75*runif(n*k),n,k)
  X <- matrix(rpois(n*m,L %*% t(F)),n,m)
  Y <- as(X,"dgCMatrix")
  
  # Generate random initial estimates of the factors and loadings.
  fit0 <- list(F = matrix(runif(m*k),m,k),
               L = matrix(runif(n*k),n,k))

  # Fit the non-negative matrix factorization to the sparse and dense
  # matrices, and check that the two solutions are the same.
  capture.output(fit1 <- altsqp(X,fit0,numiter = 20,version = "R"))
  capture.output(fit2 <- altsqp(X,fit0,numiter = 20,version = "Rcpp"))
  capture.output(fit3 <- altsqp(Y,fit0,numiter = 20,version = "R"))
  capture.output(fit4 <- altsqp(Y,fit0,numiter = 20,version = "Rcpp"))
  fit1$progress       <- fit1$progress[1:4]
  fit2$progress       <- fit2$progress[1:4]
  fit3$progress       <- fit3$progress[1:4]
  fit4$progress       <- fit4$progress[1:4]
  expect_equal(fit1,fit2)
  expect_equal(fit1,fit3)
  expect_equal(fit1,fit4)
})

test_that(paste("Multicore version of altsqp gives same result as",
                "single-core version"),{

  # Generate a 100 x 200 data matrix to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  k <- 3
  F <- matrix(0.75*runif(m*k),m,k)
  L <- matrix(0.75*runif(n*k),n,k)
  X <- matrix(rpois(n*m,L %*% t(F)),n,m)
  Y <- as(X,"dgCMatrix")

  # Generate random initial estimates of the factors and loadings.
  fit0 <- list(F = matrix(runif(m*k),m,k),
               L = matrix(runif(n*k),n,k))

  # Fit the non-negative matrix factorization to the sparse and dense
  # matrices, and check that the two solutions are the same.
  capture.output(fit1 <- altsqp(X,fit0,numiter = 10,version = "R"))
  capture.output(fit2 <- altsqp(X,fit0,numiter = 10,version = "R",
                                control = list(nc = 2)))
  capture.output(fit3 <- altsqp(X,fit0,numiter = 10,version = "Rcpp",
                                control = list(nc = 2)))
  fit1$progress       <- fit1$progress[1:4]
  fit2$progress       <- fit2$progress[1:4]
  fit3$progress       <- fit3$progress[1:4]
  expect_equal(fit1,fit2)
  expect_equal(fit1,fit3)
})

test_that("altsqp gives a better solution than nnmf on a sparse matrix",{
  library(Matrix)
  library(NNLM)

  # Generate a 100 x 200 data matrix to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  k <- 3
  F <- matrix(0.6*runif(m*k),m,k)
  L <- matrix(0.6*runif(n*k),n,k)
  X <- matrix(rpois(n*m,L %*% t(F)),n,m)
  X <- as(X,"dgCMatrix")

 # Generate random initial estimates of the factors and loadings.
 fit0 <- list(F = matrix(runif(m*k),m,k),
              L = matrix(runif(n*k),n,k))

 # Fit the non-negative matrix factorization using nnmf and altsqp.
 fit1 <- suppressWarnings(
   nnmf(as.matrix(X),k,init = list(W = fit0$L,H = t(fit0$F)),
        method = "scd",loss = "mkl",max.iter = 50,rel.tol = 0, 
        inner.max.iter = 4,verbose = 0))
 capture.output(fit2 <- altsqp(X,fit0,numiter = 50))

 # The altsqp solution should yield a higher likelihood than the nnmf
 # solution.
 fit1$F <- t(fit1$H)
 fit1$L <- fit1$W
 f1     <- loglik.poisson(X,fit1)
 f2     <- loglik.poisson(X,fit2)
 expect_lt(f1,f2) 
})

test_that("altsqp gives the same result as betanmf when numsqp = 0",{

  # Generate a 100 x 200 data matrix to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  k <- 3
  F <- matrix(0.75*runif(m*k),m,k)
  L <- matrix(0.75*runif(n*k),n,k)
  X <- matrix(rpois(n*m,L %*% t(F)),n,m)
  
  # Generate random initial estimates of the factors and loadings.
  fit0 <- list(F = matrix(runif(m*k),m,k),
               L = matrix(runif(n*k),n,k))

  # Fit the non-negative matrix factorization using the multiplicative
  # updates (betanmf) and using the alternating SQP updates.
  fit1 <- betanmf(X,fit0$L,t(fit0$F),numiter = 10)
  capture.output(fit2 <- altsqp(X,fit0,numiter = 10,version = "R",
                                control = list(numsqp = 0,extrapolate = Inf)))
  capture.output(fit3 <- altsqp(X,fit0,numiter = 10,version = "Rcpp",
                                control = list(numsqp = 0,extrapolate = Inf)))
  expect_equal(fit1$A,fit2$L)
  expect_equal(fit1$B,t(fit2$F))
  expect_equal(fit1$A,fit3$L)
  expect_equal(fit1$B,t(fit3$F))
})
