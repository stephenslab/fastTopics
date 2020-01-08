context("mixsqp")

test_that("mixsqp_rcpp and KWDual produce same result",{

  # Generate small data set.
  set.seed(1)
  out <- generate_poismix_data(100,c(1,2,0,0,0,4,0,0))
  L   <- out$L
  w   <- out$w

  # First check that the mix-SQP updates have converged to a
  # stationary point by running a few more EM updates, and checking
  # that the solution doesn't change much.
  m   <- ncol(L)
  x0  <- runif(m)
  out <- mixsqp_rcpp(L,w,x0,10,test_mixsqp_control)
  x1  <- drop(out$x)
  x2  <- mixem(L,w,x1,100)
  expect_equal(x1,x2,tolerance = 1e-6)

  # The mix-SQP updates should monotonically increase the objective
  # function.
  expect_nonincreasing(drop(out$objective))
  
  # mix-SQP and KWDual (MOSEK) should give nearly the same solution.
  skip_if_mixkwdual_doesnt_work()
  x3 <- REBayes::KWDual(normalize.cols(L),rep(1,m),w/sum(w))$f
  expect_equal(x1,x3,tolerance = 1e-4)
})

test_that("mixsqp_rcpp produces correct result when sum(w > 0) = 1",{

  # Generate the data set.
  set.seed(1)
  n    <- 10
  out  <- generate_poismix_data(n,c(1,2,0,0))
  L    <- out$L
  w    <- rep(0,n)
  w[8] <- 2
    
  # mix-SQP should nearly recover the exact solution.
  m   <- ncol(L)
  x0  <- runif(m)
  x1  <- mixture.one.nonzero(L,w)
  out <- mixsqp_rcpp(L,w,x0,10,test_mixsqp_control)
  x2  <- drop(out$x)
  expect_equal(x1,x2,tolerance = 1e-6)
})
