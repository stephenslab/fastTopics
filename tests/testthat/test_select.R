context("select")

test_that(paste("Select S3 method correctly subsets and re-orders the",
                "factors and loadings in a small example"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  n <- 80
  m <- 100
  k <- 3
  X <- generate_test_data(n,m,k)$X

  # Run 20 EM updates.
  capture.output(fit <- fit_poisson_nmf(X,k = k,numiter = 20,method = "em"))

  # Select and re-order factors and loadings by number.
  n0   <- 40
  m0   <- 50
  rows <- sample(n,n0)
  cols <- sample(m,m0)
  fit1 <- select(fit,factors = cols,loadings = rows)

  # Select and re-order factors and loadings by name.
  rows <- rownames(X)[rows]
  cols <- colnames(X)[cols]
  fit2 <- select(fit,factors = cols,loadings = rows)
  
  # Check the outputted Poisson NMF fits.
  expect_equal(dim(fit1$F),c(m0,k))
  expect_equal(dim(fit1$L),c(n0,k))
  expect_equal(dim(fit2$F),c(m0,k))
  expect_equal(dim(fit2$L),c(n0,k))
  expect_equal(rownames(fit2$F),cols)
  expect_equal(rownames(fit2$L),rows)

  # An error is thrown when the selected factors or loadings do not
  # exist.
  expect_error(select(fit,factors = m+1))
  expect_error(select(fit,loadings = n+1))
})
