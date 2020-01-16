context("mixem")

test_that("mixem and mixem_rcpp produce same result",{

  # Generate small data set.
  set.seed(1)
  out <- generate_poismix_data(100,c(1,2,0,0,0,4,0,0))
  L   <- out$L
  w   <- out$w
  
  # Run 100 EM updates for the multinomial mixture model. The R and
  # C++ implementations should give nearly the same result.
  m  <- ncol(L)
  x0 <- runif(m)
  x1 <- mixem(L,w,x0,100)
  x2 <- drop(mixem_rcpp(L,w,x0,100))
  expect_equal(x1,x2,tolerance = 1e-12,scale = 1)
})

test_that("mixem and mixem_rcpp produce correct result when sum(w > 0) = 1",{

  # Generate the data set.
  set.seed(1)
  n    <- 10
  out  <- generate_poismix_data(n,c(1,2,0,0))
  L    <- out$L
  w    <- rep(0,n)
  w[8] <- 2

  # Get the exact solution.
  x0 <- mixture.one.nonzero(L,w)
  
  # Run 20 EM updates for the multinomial mixture model.
  x1 <- mixem(L,w,x0,20)
  x2 <- drop(mixem_rcpp(L,w,x0,20))

  # The solution should not change much after running the EM updates.
  expect_equal(x0,x1,tolerance = 1e-12,scale = 1)
  expect_equal(x0,x2,tolerance = 1e-12,scale = 1)
})
