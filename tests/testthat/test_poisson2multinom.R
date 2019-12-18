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
  fit <- list(F = F,L = L)
  fit <- poisson2multinom(fit)
  expect_equivalent(colSums(fit$F),c(1,1))
  expect_equivalent(rowSums(fit$L),c(1,1,1,1))
})
