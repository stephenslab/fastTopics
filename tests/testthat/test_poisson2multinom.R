context("poisson2multinom")

test_that("poisson2multinom correctly scales factors and loadings",{
 L   <- matrix(0:7,4,2)
 F   <- matrix(0:9,5,2)
 fit <- list(F = F,L = L)
 out <- poisson2multinom(fit)
 expect_equal(colSums(out$F),c(1,1))
 expect_equal(rowSums(out$L),c(1,1,1,1))
})
