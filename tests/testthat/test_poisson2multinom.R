context("poisson2multinom")

test_that("poisson2multinom gives error when k = 1",{
  L   <- matrix(0:3,4,1)
  F   <- matrix(0:4,5,1)
  fit <- list(F = F,L = L)
  expect_error(poisson2multinom(fit))
})

test_that("poisson2multinom correctly scales factors and loadings",{
  L   <- matrix(0:7,4,2)
  F   <- matrix(0:9,5,2)
  rownames(L) <- paste0("i",1:4)
  rownames(F) <- paste0("j",1:5)
  colnames(L) <- paste0("k",1:2)
  colnames(F) <- paste0("k",1:2)
  fit         <- list(F = F,L = L)
  class(fit)  <- c("poisson_nmf_fit","list")
  fit         <- poisson2multinom(fit)
  expect_equivalent(colSums(fit$F),c(1,1))
  expect_equivalent(rowSums(fit$L),c(1,1,1,1))
})

test_that("multinom2poisson recovers original Poisson NMF model fit",{
  set.seed(1)
  out  <- simulate_count_data(10,20,3)
  X    <- out$X
  fit1 <- iterate_updates(X,out$F,out$L,100,
                          function (X,F,L) t(betanmf_update_factors(X,L,t(F))),
                          function (X,F,L) betanmf_update_loadings(X,L,t(F)))
  class(fit1) <- c("poisson_nmf_fit","list")
  fit2    <- poisson2multinom(fit1)
  fit2a   <- fit2
  fit2a$s <- NULL
  fit3 <- multinom2poisson(fit2)
  fit4 <- multinom2poisson(fit2a,X)
  fit5 <- multinom2poisson(fit2a,as(X,"dgCMatrix"))
  Y1   <- with(fit1,tcrossprod(L,F))
  Y3   <- with(fit3,tcrossprod(L,F))
  Y4   <- with(fit4,tcrossprod(L,F))
  Y5   <- with(fit5,tcrossprod(L,F))
  f1   <- loglik_poisson_nmf(X,fit1)
  f3   <- loglik_poisson_nmf(X,fit3)
  f4   <- loglik_poisson_nmf(X,fit4)
  f5   <- loglik_poisson_nmf(X,fit5)
  expect_equal(Y1,Y3,tolerance = 1e-15,scale = 1)
  expect_equal(Y1,Y4,tolerance = 1e-15,scale = 1)
  expect_equal(Y1,Y5,tolerance = 1e-15,scale = 1)
  expect_equal(f1,f3,tolerance = 1e-14,scale = 1)
  expect_equal(f1,f4,tolerance = 1e-14,scale = 1)
  expect_equal(f1,f5,tolerance = 1e-14,scale = 1)
})
