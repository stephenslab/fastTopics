context("mixsqp")

test_that("Verify mixsqp and mixem on 1000 x 10 matrix",{
  
  # Load the data.
  load("mixdata.RData")
  L <- mixdata$L
  w <- mixdata$w

  # Initialize the solution estimate.
  m  <- ncol(L)
  x0 <- rep(1/m,m)

  # Fit model by iterating the SQP and EM updates.
  capture.output(fit1 <- mixsqp(L,w,x0,numiter = 14,verbose = TRUE))
  fit2 <- mixem(L,w,x0,numiter = 1000)

  # Verify the solutions. Note that the EM solution is not expected to
  # be very close to the true solution.
  expect_equal(fit1$value,mixdata$value,tolerance = 1e-8)
  expect_equal(fit2$value,mixdata$value,tolerance = 1e-4)
  expect_equal(fit1$x,mixdata$x,tolerance = 1e-8)
})

test_that("Verify mixsqp on tacks data",{

  # Load the data.
  load("tacks.RData")
  L <- tacks$L
  w <- tacks$w
  f <- mixobjective(L,w,tacks$x,0)

  # Initialize the solution estimate.
  m  <- ncol(L)
  x0 <- rep(1/m,m)

  # Fit model by iterating the SQP and EM updates.
  fit1 <- mixem(L,w,x0,numiter = 10)
  capture.output(fit2 <- mixsqp(L,w,fit1$x,numiter = 30,verbose = TRUE))

  # Verify the solution.
  expect_equal(fit2$value,f,tolerance = 1e-6)
})
