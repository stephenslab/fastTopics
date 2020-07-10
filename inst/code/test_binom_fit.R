# TO DO: Explain what this script does, and how to use it.

# Simulate a binomial data set.
set.seed(1)
n  <- 200
s  <- ceiling(100*runif(n))
p0 <- 0.1
p1 <- 0.4
q  <- runif(n)
p  <- q*p1 + (1-q)*p0
x  <- rbinom(n,s,p)

# Fit the model parameters, p0 and p1, using optim.
out1 <- fit_binom_optim(x,y,q)

# Fit the model parameters using EM.
out2 <- fit_binom_em(x,y,q)

# optim and EM should give the same solution.
print(out1$par)
print(out2$p)
