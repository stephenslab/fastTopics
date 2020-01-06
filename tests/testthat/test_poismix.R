context("poismix")

test_that("poismixem and poismixem_rcpp produce same result",{

  # Generate small data set.
  set.seed(1)
  out <- generate_poismix_data(100,c(1,2,0,0,0,8,0,0))
  L   <- out$L
  w   <- out$w
  
  # Run 100 EM updates for the Poisson mixture model. The R and C++
  # implementations should give nearly the same result.
  numiter <- 100
  x0 <- runif(8)
  x1 <- poismixem(L,w,x0,numiter)
  x2 <- drop(poismixem_rcpp(L,w,x0,numiter))
  expect_equal(x1,x2,tolerance = 1e-14)
})

test_that(paste("poismixem and poismixem_rcpp produce correct result",
                "when sum(w > 0) = 1"),{

  # Generate the data set.
  set.seed(1)
  n <- 10
  m <- 4
  x <- c(1,2,0,0)
  L <- matrix(runif(n*m),n,m)
  w <- rep(0,n)
  w[8] <- 2

  # Run 100 EM updates for the multinomial mixture model.
  numiter <- 100
  x0 <- runif(m)
  x1 <- poismixem(L,w,x0,numiter)
  x2 <- drop(poismixem_rcpp(L,w,x0,numiter))
  
  # The R and C++ implementations should give nearly the same result,
  # and should be very close to the exact solution obtained by calling
  # poismix.one.nonzero.
  x3 <- poismix.one.nonzero(L,w)
  expect_equal(x1,x2,tolerance = 1e-14)
  expect_equal(x1,x3,tolerance = 1e-14)
  expect_equal(x2,x3,tolerance = 1e-14)
})
