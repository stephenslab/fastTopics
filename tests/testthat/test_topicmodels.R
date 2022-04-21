context("topicmodels")

test_that("R and C++ versions of topicmodels::cost gives same result",{
  set.seed(1)
  n <- 8
  m <- 12
  k <- 4
  X <- matrix(rpois(n*m,4),n,m)
  Y <- as(X,"dgCMatrix")
  A <- matrix(abs(rnorm(n*k)),n,k)
  B <- matrix(abs(rnorm(k*m)),k,m)
  y1 <- topicmodels:::cost(X,A,B,version = "r")
  y2 <- topicmodels:::cost(Y,A,B,version = "cpp")
  expect_equal(y1,y2,scale = 1,tolerance = 1e-10)
})
