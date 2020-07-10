# A short script to check that fit_binom_optim and fit_binom_em give
# the same result, and to verify that the Poisson approximation to the
# binomial likelihood works well when the binomial success rates are
# small relative to the number of trials, and, as a result,
# fit_poisson_em should recover nearly the same parameter estimates as
# fit_binom_em.

# Simulate a binomial data set.
set.seed(1)
n  <- 200
s  <- ceiling(100*runif(n))
p0 <- 0.05
p1 <- 0.1
q  <- runif(n)
p  <- q*p1 + (1-q)*p0
x  <- rbinom(n,s,p)

# Fit the model parameters, p0 and p1, using optim.
out1 <- fit_binom_optim(x,s,q)

# Fit the model parameters using EM.
out2 <- fit_binom_em(x,s,q)

# Fit the model parameters using the Poisson approximation to the
# binomial likelihood.
out3 <- fit_poisson_em(x,s,q)

# optim and EM should give the same solution.
print(data.frame(optim = out1$par,
                 em    = out2$p,
                 pois  = out3$f))
