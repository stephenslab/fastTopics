# TO DO: Explain here what this script does, and how to use it.

# Simulate the binomial data set.
set.seed(1)
n  <- 200
s  <- ceiling(100*runif(n))
p0 <- 0.1
p1 <- 0.4
q  <- runif(n)
p  <- q*p1 + (1-q)*p0
n1 <- rbinom(n,s,p)
n0 <- s - n1    

# Fit the model parameters, p0 and p1.
out <- fit_binom_topic_model(n0,n1,q)
print(out)
