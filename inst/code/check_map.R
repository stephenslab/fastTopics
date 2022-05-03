library(Compositional)
set.seed(1)
n <- 80
m <- 100
k <- 3
dat <- simulate_count_data(n,m,k)
X <- dat$X
L <- dat$L
F <- dat$F
A <- matrix(abs(rnorm(m*k)) + 1,m,k)
b <- abs(rnorm(k))

N <- 100
f0 <- rep(0,N)
f1 <- rep(0,N)
f2 <- rep(0,N)
f3 <- rep(0,N)
for (i in 1:N) {
  fit <- list(L = L,F = F)
  class(fit) <- c("poisson_nmf_fit","list")
  
  # Compute multinomial topic model likelihood + Dirichlet prior.
  f0[i]  <- sum(loglik_multinom_topic_model(X,poisson2multinom(fit),e = 0))
  for (j in 1:k)
    f0[i] <- f0[i] + ddiri(fit$F[,j],A[,j],logged = TRUE)

  # Compute multinomial topic model likelihood with "pseudodata".
  Y      <- rbind(X,t(A - 1))
  fit2   <- fit
  fit2$L <- rbind(fit$L,diag(abs(rnorm(k))))
  f1[i]  <- sum(loglik_multinom_topic_model(Y,poisson2multinom(fit2),e = 0))

  # Compute Poisson NMF likelihood with "pseudodata" (here we have to
  # remove the "size factor" likelihoods).
  f2[i]  <- sum(loglik_poisson_nmf(Y,fit2,e = 0)) -
            sum(loglik_size_factors(Y,fit2$F,fit2$L))

  # Compute the Poisson NMF likelihood + gamma prior.
  f3[i]  <- sum(loglik_poisson_nmf(X,fit,e = 0)) -
            sum(loglik_size_factors(X,fit$F,fit$L))
  for (j in 1:k)
    f3[i] <- f3[i] + sum(dgamma(fit$F[,j],A[,j],b[j],log = TRUE))
  
  # Tweak the fit.
  L <- L * matrix(exp(rnorm(n*k,sd = 0.1)),n,k)
  F <- F * matrix(exp(rnorm(m*k,sd = 0.1)),m,k)
}

plot(f0,f1,pch = 20)
plot(f0,f2,pch = 20)
plot(f2,f3,pch = 20)
