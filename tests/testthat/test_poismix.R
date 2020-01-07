context("poismix")

test_that("poismixem and poismixem_rcpp produce same result",{

  # Generate small data set.
  set.seed(1)
  out <- generate_poismix_data(100,c(1,2,0,0,0,4,0,0))
  L   <- out$L
  w   <- out$w
  
  # Run 100 EM updates for the Poisson mixture model. The R and C++
  # implementations should give nearly the same result.
  numiter <- 100
  L1 <- normalize.cols(L)
  u  <- colSums(L)
  i  <- which(w > 0)
  x0 <- runif(8)
  x1 <- poismixem(L,w,x0,numiter)
  x2 <- drop(poismixem_rcpp(L,w,x0,numiter))
  x3 <- drop(poismixem2_rcpp(L1,w,u,x0,numiter))
  x4 <- drop(poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter))
  expect_equal(x1,x2,tolerance = 1e-14)
  expect_equal(x1,x3,tolerance = 1e-14)
  expect_equal(x1,x4,tolerance = 1e-14)
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
  i <- 8
  w[i] <- 2

  # Run 100 EM updates for the multinomial mixture model.
  numiter <- 100
  x0 <- runif(m)
  x1 <- poismixem(L,w,x0,numiter)
  x2 <- drop(poismixem_rcpp(L,w,x0,numiter))

  # Run 100 EM updates another couple times, using the alternate C++
  # interfaces.
  L1 <- normalize.cols(L)
  u  <- colSums(L)
  x  <- poismixem2_rcpp(L1,w,u,x0,numiter)
  y  <- poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter)
  x3 <- drop(x)
  expect_equal(x,y,tolerance = 1e-15)
  
  # The R and C++ implementations should give nearly the same result,
  # and should be very close to the exact solution obtained by calling
  # poismix.one.nonzero.
  x4 <- poismix.one.nonzero(L,w)
  expect_equal(x1,x2,tolerance = 1e-14)
  expect_equal(x1,x4,tolerance = 1e-14)
  expect_equal(x2,x4,tolerance = 1e-14)
  expect_equal(x3,x4,tolerance = 1e-14)
})
