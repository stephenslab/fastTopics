# A short script to check that fit_poisson_optim gives the correct
# result.

# Simulate a Poisson data set.
set.seed(1)
n  <- 200
f0 <- 1
f1 <- 4
s  <- sample(100,n,replace = TRUE)
q  <- runif(n)
u  <- (1-q)*f0 + q*f1
x  <- rpois(n,s*u)

# Fit the model parameters, f0 and f1, using optim.
out <- fit_poisson_optim(x,s,q)
