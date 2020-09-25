context("select")

test_that(paste("Select S3 method correctly subsets and re-orders the",
                "factors and loadings in a small example; also check ",
                "merge_topics"),{

  # Generate a 80 x 100 data matrix to factorize.
  set.seed(1)
  n <- 80
  m <- 100
  k <- 3
  X <- generate_test_data(n,m,k)$X

  # Run 20 EM updates.
  capture.output(
    fit <- poisson2multinom(fit_poisson_nmf(X,k = k,numiter = 20,
                                            method = "em")))
  
  # Select and re-order factors and loadings by number (here, we use
  # the "select_loadings" function).
  n0   <- 40
  rows <- sample(n,n0)
  fit1 <- select_loadings(fit,rows)

  # Select and re-order factors and loadings by name (here, we use the
  # "select" S3 method).
  rows <- rownames(X)[rows]
  fit2 <- select(fit,rows)
  
  # Check the outputted Poisson NMF fits.
  expect_equal(dim(fit1$L),c(n0,k))
  expect_equal(dim(fit2$L),c(n0,k))
  expect_equal(rownames(fit2$L),rows)

  # An error is thrown when the selected loadings do not exist.
  expect_error(select(fit,loadings = n + 1))

  # Check that merge_topics does the right thing.
  fit3 <- merge_topics(fit2,k = 1:2)
  fit4 <- merge_topics(fit2,k = c("k1","k2"))
  expect_equal(dim(fit3$F),c(100,2))
  expect_equal(dim(fit3$L),c(40,2))
  expect_equal(colnames(fit3$F),c("k3","k1+k2"))
  expect_equal(colnames(fit3$L),c("k3","k1+k2"))
  expect_equal(fit3,fit4)
})
