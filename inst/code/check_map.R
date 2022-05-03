library(Compositional)

# Simulate a 80 x 100 data set.
set.seed(1)
n <- 80
m <- 100
k <- 3
dat <- simulate_count_data(n,m,k)
X <- dat$X
L <- dat$L
F <- dat$F
a <- matrix(abs(rnorm(m*k)) + 1,m,k)
b <- abs(rnorm(k))

N <- 100
f0 <- rep(0,N)
f1 <- rep(0,N)
f2 <- rep(0,N)
f3 <- rep(0,N)
for (i in 1:N) {
  
  # Compute the penalized likelihood for the multinomial topic model
  # with a Dirichlet prior.
  fit <- list(L = L,F = F)
  class(fit) <- c("poisson_nmf_fit","list")
  f0[i]  <- sum(loglik_multinom_topic_model(X,poisson2multinom(fit),e = 0))
  for (j in 1:k)
    f0[i] <- f0[i] + ddiri(fit$F[,j],a[,j],logged = TRUE)
  
  # Compute the multinomial topic model likelihood with "pseudodata".
  Y      <- rbind(X,t(a - 1))
  fit2   <- fit
  u      <- colSums(a - 1)/b
  fit2$L <- rbind(fit$L,diag(1/u))
  f1[i]  <- sum(loglik_multinom_topic_model(Y,poisson2multinom(fit2),e = 0))

  # Compute the penalized Poisson NMF likelihood with a gamma prior.
  f2[i] <- sum(loglik_poisson_nmf(X,fit,e = 0))
  for (j in 1:k)
    f2[i] <- f2[i] + sum(dgamma(fit$F[,j],a[,j],b[j],log = TRUE))
  
  # Compute Poisson NMF likelihood with "pseudodata".
  f3[i] <- sum(loglik_poisson_nmf(Y,fit2,e = 0))

  # Tweak the fit.
  L <- L * matrix(exp(rnorm(n*k,sd = 0.1)),n,k)
  F <- F * matrix(exp(rnorm(m*k,sd = 0.1)),m,k)
}

# The multinomial penalized log-likelihoods and the multinomial
# log-likelihoods with pseudodata should be equal up to a constant.
plot(f0,f1,pch = 20)

# The Poisson NMF penalized log-likelihoods and the Poisson NMF
# log-likleihoods with pseudodata should be equal up to a constant.
plot(f2,f3,pch = 20)
