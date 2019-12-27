context("mixem")

test_that("poismixem and poismixem_rcpp produce same result",{

  # Generate small data set.
  set.seed(1)
  n <- 100
  m <- 8
  x <- c(1,2,0,0,0,8,0,0)
  L <- matrix(runif(n*m),n,m)
  w <- rpois(n,L %*% x)

  # Run 100 EM updates for the Poisson mixture model. The R and C++
  # implementations should give nearly the same result.
  x0 <- runif(m)
  x1 <- poismixem(L,w,x0,100)
  x2 <- drop(poismixem_rcpp(L,w,x0,100))
  expect_equal(x1,x2,tolerance = 1e-12)
})

test_that("poismixem and poismixem_rcpp work when n = 1",{

  # Generate the data set.
  set.seed(1)
  x <- c(1,2,0,0)
  L <- matrix(runif(4),1,4)
  w <- 2

  # The R and C++ implementations should give nearly the same result,
  # and should be very close to the exact solution obtained by calling
  # poismix.one.nonzero.
  # TO DO.
})
