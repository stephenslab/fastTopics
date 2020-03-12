context("poismix")

test_that("poismixem and poismixem_rcpp produce same result",{

  # Generate small data set.
  set.seed(1)
  out <- generate_poismix_data(100,c(1,2,0,0,0,4,0,0))
  L   <- out$L
  w   <- out$w
  
  # Run 100 EM updates for the Poisson mixture model. The R
  # implementation, and all variations of the C++ implementation,
  # should give nearly the same result.
  numiter <- 100
  m  <- ncol(L)
  L1 <- normalize.cols(L)
  u  <- colSums(L)
  i  <- which(w > 0)
  x0 <- runif(m)
  x1 <- poismixem(L,w,x0,numiter)
  x2 <- drop(poismixem_rcpp(L,w,x0,numiter))
  x3 <- drop(poismixem2_rcpp(L1,w,u,x0,numiter))
  x4 <- drop(poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter))
  expect_equal(x1,x2,tolerance = 1e-14,scale = 1)
  expect_equal(x1,x3,tolerance = 1e-14,scale = 1)
  expect_equal(x1,x4,tolerance = 1e-14,scale = 1)
})

test_that(paste("poismixem, scd_kl_update and ccd_kl_update give nearly the",
                "same solution"),{

  # Generate small data set.
  set.seed(1)
  out <- generate_poismix_data(100,c(1,2,0,0,0,4,0,0))
  L   <- out$L
  w   <- out$w

  # Run 10,000 EM updates.
  m  <- ncol(L)
  x0 <- runif(m)
  x1 <- drop(poismixem_rcpp(L,w,x0,1e4))

  # Run 100 sequential coordinate descent (SCD) updates, using both
  # C++ interfaces.
  numiter <- 100
  L1 <- normalize.cols(L)
  u  <- colSums(L)
  i  <- which(w > 0)
  x2 <- drop(scd_kl_update_rcpp(L,w,x0,numiter,1e-15))
  x3 <- drop(scd_kl_update2_rcpp(L[i,],u,w[i],x0,numiter,1e-15))

  # Run 100 cyclic coordinate descent (CCD) updates, using both C++
  # interfaces.
  x4 <- drop(ccd_kl_update_rcpp(L,w,x0,numiter,1e-15))
  x5 <- drop(ccd_kl_update2_rcpp(L[i,],u,w[i],x0,numiter,1e-15))
  
  # The coordinatewise updates should recover nearly the same solution
  # as mix-SQP, and should give the same results whether the "dense"
  # or "sparse" updates are used.
  expect_equal(x1,x2,tolerance = 1e-5,scale = 1)
  expect_equal(x1,x4,tolerance = 1e-5,scale = 1)
  expect_equal(x2,x3,tolerance = 1e-14,scale = 1)
  expect_equal(x4,x5,tolerance = 1e-14,scale = 1)
})

test_that(paste("poismixem and poismixem_rcpp produce correct result",
                "when sum(w > 0) = 1"),{

  # Generate the data set.
  set.seed(1)
  n    <- 10
  out  <- generate_poismix_data(n,c(1,2,0,0))
  L    <- out$L
  w    <- rep(0,n)
  i    <- 8
  w[i] <- 2

  # Run 100 EM updates for the multinomial mixture model.
  numiter <- 100
  m  <- ncol(L)
  x0 <- runif(m)
  x1 <- poismixem(L,w,x0,numiter)

  # Run 100 EM updates another a few times, using the different C++
  # interfaces.
  L1 <- normalize.cols(L)
  u  <- colSums(L)
  x2 <- drop(poismixem_rcpp(L,w,x0,numiter))
  x3 <- drop(poismixem2_rcpp(L1,w,u,x0,numiter))
  x4 <- drop(poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter))
  
  # The R and C++ implementations should give nearly the same result,
  # and should be very close to the exact solution obtained by calling
  # poismix.one.nonzero.
  x5 <- poismix.one.nonzero(L,w)
  expect_equal(x1,x2,tolerance = 1e-12,scale = 1)
  expect_equal(x1,x3,tolerance = 1e-12,scale = 1)
  expect_equal(x1,x4,tolerance = 1e-12,scale = 1)
  expect_equal(x1,x5,tolerance = 1e-12,scale = 1)
})
