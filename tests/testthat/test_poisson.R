context("poisson")

test_that(paste("fastTopics methods recover the same model fits as glm",
                "with family = poisson(link = \"identity\")"),{

  # Simulate a Poisson data set.
  set.seed(1)
  n  <- 40
  f0 <- 0.1
  f1 <- 1
  s  <- sample(10,n,replace = TRUE)
  q  <- runif(n)
  u  <- (1-q)*f0 + q*f1
  x  <- rpois(n,s*u)

  # Fit the Poisson model using glm.
  fit <- fit_poisson_glm(x,s,q)

  # Fit the model using optim.
  out1 <- fit_poisson_optim(x,s,q)

  # Fit the model using EM.
  out2 <- fit_poisson_em(x,s,q)

  # Fit the model using the C++ EM implementation,
  out3 <- fit_poisson_em_rcpp(x,s,q,1,1,1e-15,100,1e-8)

  # Fit the model using the C++ implementation of the EM algorithm
  # that is better suited for sparse counts.
  i    <- which(x > 0)
  out4 <- fit_poisson_em_sparse_rcpp(x[i],s[i],q[i],sum(s*(1-q)),sum(s*q),
                                     1,1,1e-15,100,1e-8)

  # The glm fit and all four fastTopics methods should produce the
  # same, or nearly the same, estimates of f0 and f1.
  expect_equal(fit[c("f0","f1")],out1$par,scale = 1,tolerance = 1e-6)
  expect_equal(fit[c("f0","f1")],out2$f,scale = 1,tolerance = 1e-6)
  expect_equal(fit[c("f0","f1")],c(f0 = out3$f0,f1 = out3$f1),
               scale = 1,tolerance = 1e-6)
  expect_equal(fit[c("f0","f1")],c(f0 = out4$f0,f1 = out4$f1),
               scale = 1,tolerance = 1e-6)

  # The likelihood should be the same in the EM implementations.
  expect_equal(tail(out2$loglik,n = 1),out3$loglik,scale = 1,tolerance = 1e-8)
  expect_equal(tail(out2$loglik,n = 1),out4$loglik,scale = 1,tolerance = 1e-8)
})

test_that(paste("fit_poisson_em gives nearly the same result as fit_binom_em",
                "when the binomial probabilities are small relative to the",
                "number of binomial trials"),{

  # Simulate a binomial data set.
  set.seed(1)
  n  <- 200
  s  <- ceiling(100*runif(n))
  p0 <- 0.05
  p1 <- 0.1
  q  <- runif(n)
  p  <- q*p1 + (1-q)*p0
  x  <- rbinom(n,s,p)

  # Fit the binomial model.
  out1 <- fit_binom_em(x,s,q)

  # Fit the Poisson model
  out2 <- fit_poisson_em(x,s,q)

  # The binomial model parameter estimates should be nearly the same
  # as the Poisson model parameter estimates.
  expect_equivalent(out1$p,out2$f,scale = 1,tolerance = 1e-4)
})
