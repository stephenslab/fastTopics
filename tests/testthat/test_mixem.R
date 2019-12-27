context("mixem")

test_that("mixem and mixem_rcpp produce same result",{

  # Generate small data set.
  set.seed(1)
  n <- 100
  m <- 8
  x <- c(1,2,0,0,0,8,0,0)
  L <- matrix(runif(n*m),n,m)
  w <- rpois(n,L %*% x)

  # Run 100 EM updates for the multinomial mixture model. The R and
  # C++ implementations should give nearly the same result.
  x0 <- runif(m)
  x1 <- mixem(L,w,x0,100)
  x2 <- drop(mixem_rcpp(L,w,x0,100))
  expect_equal(x1,x2,tolerance = 1e-12)
})

test_that("mixem and mixem_rcpp produce correct result when sum(w > 0) = 1",{

  # Generate the data set.
  set.seed(1)
  n <- 10
  m <- 4
  x <- c(1,2,0,0)
  L <- matrix(runif(n*m),n,m)
  w <- rep(0,n)
  w[8] <- 2

  # Run 100 EM updates for the multinomial mixture model.
  x0 <- runif(m)
  x1 <- mixem(L,w,x0,100)
  x2 <- drop(mixem_rcpp(L,w,x0,100))

  # The R and C implementations should give nearly the same result,
  # and should be very close to the exact solution obtained by calling
  # mixture.one.nonzero.
  x3 <- mixture.one.nonzero(L,w)
  expect_equal(x1,x2,tolerance = 1e-14)
  expect_equal(x1,x3,tolerance = 1e-14)
  expect_equal(x2,x3,tolerance = 1e-14)
})
