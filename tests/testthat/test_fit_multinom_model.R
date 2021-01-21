context("fit_multinom_model")

test_that("fit_multinom_model gives correct factor estimates",{

  # Simulate a "toy" gene expression data set.
  set.seed(1)
  n   <- 400
  m   <- 40
  k   <- 3
  out <- simulate_toy_gene_data(n,m,k,s = 1000)
  X   <- out$X
  Y   <- as(X,"dgCMatrix")

  # Force "hard" topic assignments.
  cluster <- factor(apply(force_hard_topic_assignments(out$L),1,which.max))
  levels(cluster) <- paste0("k",1:k)

  # Fit the simple multinomial model.
  fit1 <- fit_multinom_model(cluster,X)
  fit2 <- fit_multinom_model(cluster,Y)
  
  # Both calls to fit_multinom_model should result in nearly the same
  # loadings.
  expect_equal(fit1$L,fit2$L,scale = 1,tolerance = 1e-15)

  # Check that both calls to fit_multinom_model recover the
  # maximum-likelihood estimates of the factors.
  F <- matrix(0,m,k)
  for (j in 1:k) {
    i     <- which(cluster == levels(cluster)[j])
    F[,j] <- colSums(X[i,])/sum(X[i,])
  }
  expect_equivalent(fit1$F,F,scale = 1,tolerance = 1e-15)
  expect_equivalent(fit2$F,F,scale = 1,tolerance = 1e-15)
})

