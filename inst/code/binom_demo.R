# Simulate a binomial data set.
set.seed(1)
n  <- 200
s  <- ceiling(100*runif(n))
p0 <- 0.1
p1 <- 0.4
q  <- runif(n)
p  <- q*p1 + (1-q)*p0
x  <- rbinom(n,s,p)
y  <- s - x

# Fit the model parameters, p0 and p1.
out <- fit_binom_topic_model(x,y,q)
print(out)
